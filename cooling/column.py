"""1D finite-volume heat conduction through stacked layers.

The column is a stack of homogeneous layers (bottom to top, e.g. ground and
molasses) with z pointing up and z = 0 at the bottom of the column. Cells are
cell-centred finite volumes on a graded grid: fine at the top surface and at
every interface between layers, coarse inside the layers and at the bottom
(the far field). The face between two cells conducts through both half
cells in series, which is exact for piecewise-constant conductivity.

Time stepping is implicit Euler. Each boundary exchanges heat with a
reference temperature through a coefficient, q_out = h (T_b - T_ref), where
T_b is the temperature on the boundary face; coefficients that depend on T_b
(radiation) are iterated within the step. Because the boundary fluxes are
evaluated with the new temperatures, the change of the heat content equals
the heat that crossed the boundaries, step by step, to round-off.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Protocol

import numpy as np
from scipy.linalg import solve_banded

from cooling.properties import Material


@dataclass(frozen=True)
class Layer:
    """A homogeneous layer of the column.

    Attributes:
        material: Thermal properties.
        thickness: Thickness [m].
        initial_temperature: Uniform temperature at t = 0 [K].
    """

    material: Material
    thickness: float
    initial_temperature: float

    def __post_init__(self) -> None:
        if not self.thickness > 0.0:
            raise ValueError(f"{self.material.name}: thickness must be positive (got {self.thickness})")
        if not self.initial_temperature > 0.0:
            raise ValueError(f"{self.material.name}: initial temperature must be positive (in kelvin)")


class Boundary(Protocol):
    """Heat exchange of a boundary face with its surroundings."""

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        """Return (h, T_ref) such that the outward heat flux is h (T_b - T_ref).

        Args:
            time: Time at the end of the step [s].
            boundary_temperature: Current estimate of the face temperature T_b [K].

        Returns:
            h [W/m2/K] (0 for no exchange, math.inf for a fixed temperature)
            and T_ref [K].
        """
        ...


class Adiabatic:
    """No heat flux through the boundary."""

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        return 0.0, 0.0


@dataclass(frozen=True)
class FixedTemperature:
    """Boundary held at a temperature [K]."""

    temperature: float

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        return math.inf, self.temperature


@dataclass(frozen=True)
class Convective:
    """Exchange h (T_b - T_inf) with an ambient temperature, constant or a function of time [K]."""

    h: float
    ambient: float | Callable[[float], float]

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        ambient = self.ambient(time) if callable(self.ambient) else self.ambient
        return self.h, ambient


def graded_widths(
    length: float,
    dz_min: float,
    dz_max: float,
    growth: float,
    refine_bottom: bool,
    refine_top: bool,
) -> np.ndarray:
    """Cell widths over ``length``, fine at the refined ends.

    From each refined end the widths start at about ``dz_min`` and grow by
    ``growth`` per cell up to ``dz_max``; with both ends refined the two
    halves mirror each other. The widths of each graded part are scaled so
    that they add up to the length exactly, which keeps the growth ratio but
    may make the first cell slightly smaller than ``dz_min``. Without a
    refined end the cells are uniform and at most ``dz_max`` wide.

    Args:
        length: Length to cover [m].
        dz_min: Width of the cells at a refined end [m].
        dz_max: Largest width [m].
        growth: Ratio of neighbouring widths (>= 1).
        refine_bottom: Refine at the low end (z = 0 of the layer).
        refine_top: Refine at the high end.

    Returns:
        Widths from bottom to top [m].
    """
    if not length > 0.0:
        raise ValueError(f"length must be positive (got {length})")
    if not 0.0 < dz_min <= dz_max:
        raise ValueError(f"0 < dz_min <= dz_max required (got {dz_min}, {dz_max})")
    if not growth >= 1.0:
        raise ValueError(f"growth must be >= 1 (got {growth})")

    def one_sided(span: float) -> np.ndarray:
        widths: list[float] = []
        width, total = dz_min, 0.0
        while total < span:
            widths.append(width)
            total += width
            width = min(width * growth, dz_max)
        return np.array(widths) * (span / total)

    if refine_bottom and refine_top:
        half = one_sided(0.5 * length)
        return np.concatenate([half, half[::-1]])
    if refine_bottom:
        return one_sided(length)
    if refine_top:
        return one_sided(length)[::-1]
    n = max(1, math.ceil(length / dz_max - 1e-9))
    return np.full(n, length / n)


@dataclass(frozen=True)
class ColumnResult:
    """Fields of a column run at the output times.

    Attributes:
        times: Output times, starting at 0 [s].
        centres: Cell centres, z up from the bottom of the column [m].
        widths: Cell widths [m].
        layer_index: Layer of each cell (0 = bottom layer).
        interfaces: Heights of the interfaces between layers [m].
        temperature: Cell temperatures, shape (times, cells) [K].
        surface_temperature: Temperature on the top face [K].
        interface_temperature: Temperature on each interface, shape (times, interfaces) [K].
        surface_flux: Heat flux out through the top face [W/m2].
        bottom_flux: Heat flux out through the bottom face [W/m2].
        heat_content: Sum of rho Cp dz T over the cells [J/m2].
        heat_in: Heat that entered through both boundaries since t = 0 [J/m2].
    """

    times: np.ndarray
    centres: np.ndarray
    widths: np.ndarray
    layer_index: np.ndarray
    interfaces: np.ndarray
    temperature: np.ndarray
    surface_temperature: np.ndarray
    interface_temperature: np.ndarray
    surface_flux: np.ndarray
    bottom_flux: np.ndarray
    heat_content: np.ndarray
    heat_in: np.ndarray

    def depth_in_layer(self, layer: int) -> tuple[np.ndarray, np.ndarray]:
        """Depth below the top of ``layer`` of its cell centres, and their cell indices."""
        cells = np.flatnonzero(self.layer_index == layer)
        top = self.interfaces[layer] if layer < len(self.interfaces) else self.centres[-1] + 0.5 * self.widths[-1]
        return top - self.centres[cells], cells


class Column:
    """Stack of layers with boundary conditions at the bottom and the top.

    Args:
        layers: Layers from bottom to top.
        surface: Exchange at the top face.
        bottom: Exchange at the bottom face.
        dz_min: Cell width at the top surface and at the interfaces [m].
        dz_max: Largest cell width [m].
        growth: Ratio of neighbouring cell widths.
        max_iterations: Iterations per step for boundary coefficients that depend on T_b.
        tolerance: Convergence of the boundary temperatures within a step [K].
    """

    def __init__(
        self,
        layers: Sequence[Layer],
        surface: Boundary,
        bottom: Boundary,
        dz_min: float = 5e-4,
        dz_max: float = 0.02,
        growth: float = 1.1,
        max_iterations: int = 20,
        tolerance: float = 1e-8,
    ) -> None:
        if not layers:
            raise ValueError("a column needs at least one layer")
        self.layers = tuple(layers)
        self.surface = surface
        self.bottom = bottom
        self.max_iterations = max_iterations
        self.tolerance = tolerance

        widths, index = [], []
        for i, layer in enumerate(self.layers):
            w = graded_widths(
                layer.thickness,
                dz_min,
                dz_max,
                growth,
                refine_bottom=i > 0,
                refine_top=True,
            )
            widths.append(w)
            index.append(np.full(w.size, i))
        self.widths = np.concatenate(widths)
        self.layer_index = np.concatenate(index)
        faces = np.concatenate([[0.0], np.cumsum(self.widths)])
        self.centres = 0.5 * (faces[:-1] + faces[1:])
        self.interfaces = np.cumsum([layer.thickness for layer in self.layers])[:-1]
        self._interface_cells = np.flatnonzero(np.diff(self.layer_index))  # cell below each interface

        k = np.array([self.layers[i].material.conductivity for i in self.layer_index])
        rho_cp = np.array([self.layers[i].material.volumetric_heat_capacity for i in self.layer_index])
        self.capacity = rho_cp * self.widths  # [J/m2/K]
        self._half_conductance = 2.0 * k / self.widths  # centre to face [W/m2/K]
        g = self._half_conductance
        self.face_conductance = 1.0 / (1.0 / g[:-1] + 1.0 / g[1:])  # interior faces [W/m2/K]
        self.initial_temperature = np.array([self.layers[i].initial_temperature for i in self.layer_index])

    @staticmethod
    def _coefficients(boundary: Boundary, time: float, t_face: float, g_cell: float) -> tuple[float, float, float]:
        """Coefficient h, reference temperature and the overall coefficient from the cell centre."""
        h, t_ref = boundary.coefficients(time, t_face)
        if h < 0.0:
            raise ValueError(f"negative boundary coefficient {h}")
        u = g_cell if math.isinf(h) else g_cell * h / (g_cell + h)
        return h, t_ref, u

    @staticmethod
    def _face_temperature(h: float, t_ref: float, t_cell: float, g_cell: float) -> float:
        """Temperature on a boundary face between the cell centre and the reference."""
        if math.isinf(h):
            return t_ref
        return (g_cell * t_cell + h * t_ref) / (g_cell + h)

    def run(self, t_end: float, dt: float, output_times: Sequence[float] = ()) -> ColumnResult:
        """Integrate from t = 0 to ``t_end``.

        The time step is ``dt`` except where it is shortened to land on an
        output time.

        Args:
            t_end: End time [s].
            dt: Time step [s].
            output_times: Times to store besides 0 and ``t_end`` [s].

        Returns:
            The fields at 0, the output times and ``t_end``.
        """
        if not t_end > 0.0 or not dt > 0.0:
            raise ValueError("t_end and dt must be positive")
        stops = sorted({float(t) for t in output_times if 0.0 < t < t_end} | {float(t_end)})

        n = self.widths.size
        temp = self.initial_temperature.copy()
        g = self._half_conductance
        gf = self.face_conductance
        t_top, t_bot = temp[-1], temp[0]
        time, heat_in = 0.0, 0.0

        records = [self._record(time, temp, t_top, t_bot, 0.0, 0.0, heat_in)]
        ab = np.zeros((3, n))
        ab[0, 1:] = -gf
        ab[2, :-1] = -gf
        diag_conduction = np.zeros(n)
        diag_conduction[:-1] += gf
        diag_conduction[1:] += gf

        for stop in stops:
            while time < stop - 1e-9 * max(1.0, stop):
                step = min(dt, stop - time)
                new_time = time + step if stop - time > dt else stop
                storage = self.capacity / step
                solved = None
                for _ in range(self.max_iterations):
                    h_top, ref_top, u_top = self._coefficients(self.surface, new_time, t_top, g[-1])
                    h_bot, ref_bot, u_bot = self._coefficients(self.bottom, new_time, t_bot, g[0])
                    if solved == (h_top, ref_top, h_bot, ref_bot):
                        break  # same system as the last solve
                    solved = (h_top, ref_top, h_bot, ref_bot)
                    ab[1] = storage + diag_conduction
                    ab[1, -1] += u_top
                    ab[1, 0] += u_bot
                    rhs = storage * temp
                    rhs[-1] += u_top * ref_top
                    rhs[0] += u_bot * ref_bot
                    new = solve_banded((1, 1), ab, rhs)
                    new_top = self._face_temperature(h_top, ref_top, new[-1], g[-1])
                    new_bot = self._face_temperature(h_bot, ref_bot, new[0], g[0])
                    converged = abs(new_top - t_top) < self.tolerance and abs(new_bot - t_bot) < self.tolerance
                    t_top, t_bot = new_top, new_bot
                    if converged:
                        break
                else:
                    raise RuntimeError(f"boundary iteration did not converge at t = {new_time} s")
                # Fluxes with the coefficients the system was solved with
                # (u_top, u_bot of the last iteration), so the budget closes
                q_top = u_top * (new[-1] - ref_top)
                q_bot = u_bot * (new[0] - ref_bot)
                heat_in -= step * (q_top + q_bot)
                temp, time = new, new_time
            records.append(self._record(time, temp, t_top, t_bot, q_top, q_bot, heat_in))

        fields = list(zip(*records))
        return ColumnResult(
            times=np.array(fields[0]),
            centres=self.centres,
            widths=self.widths,
            layer_index=self.layer_index,
            interfaces=self.interfaces,
            temperature=np.array(fields[1]),
            surface_temperature=np.array(fields[2]),
            interface_temperature=np.array(fields[3]).reshape(len(records), self.interfaces.size),
            surface_flux=np.array(fields[4]),
            bottom_flux=np.array(fields[5]),
            heat_content=np.array(fields[6]),
            heat_in=np.array(fields[7]),
        )

    def _record(self, time, temp, t_top, t_bot, q_top, q_bot, heat_in):
        g = self._half_conductance
        below = self._interface_cells
        t_interface = (g[below] * temp[below] + g[below + 1] * temp[below + 1]) / (g[below] + g[below + 1])
        return (time, temp.copy(), t_top, t_interface, q_top, q_bot, float(np.sum(self.capacity * temp)), heat_in)
