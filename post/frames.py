"""Case metadata and the fields of one time step, ready for plotting."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from post.foamio import (
    Grid2D,
    field_path,
    read_entry,
    read_grid,
    read_gravity,
    read_internal_field,
    read_text,
    time_dirs,
)

LIQUID = "alpha.water"


@dataclass(frozen=True)
class CaseInfo:
    """Static information about a case: grid, model switches, initial state.

    Attributes:
        path: Case directory.
        grid: 2D grid of the mesh.
        times: Time directory names, sorted.
        g: Gravitational acceleration [m/s^2].
        step: ``activeStep`` of transportProperties (None if absent).
        solve_energy: Whether the energy equation was solved.
        temperature_dependent: Whether the liquid viscosity follows mu(T).
        rho_liquid: Liquid density [kg/m^3].
        mu0_liquid: Liquid reference viscosity [Pa s].
        h0: Initial height of the liquid column [m].
        x_dam: Initial downstream end of the liquid column [m].
        t_liquid: Initial liquid temperature [K].
        t_air: Initial air temperature [K].
    """

    path: Path
    grid: Grid2D
    times: tuple[str, ...]
    g: float
    step: int | None
    solve_energy: bool
    temperature_dependent: bool
    rho_liquid: float
    mu0_liquid: float
    h0: float
    x_dam: float
    t_liquid: float
    t_air: float

    @property
    def x_wall(self) -> float:
        """Upstream wall (left end of the domain) [m]."""
        return self.grid.extent[0]

    @property
    def front_speed(self) -> float:
        """Ritter front speed 2 sqrt(g h0) [m/s]."""
        return 2.0 * float(np.sqrt(self.g * self.h0))


@dataclass(frozen=True)
class Frame:
    """Fields of one time step on the 2D grid (arrays of shape (ny, nx)).

    Attributes:
        time: Simulation time [s].
        alpha: Liquid volume fraction, clipped to [0, 1].
        ux: Velocity component along x [m/s].
        uy: Velocity component along y [m/s].
        temperature: Temperature [K], or None if not loaded.
        dy: Cell height [m].
    """

    time: float
    alpha: np.ndarray
    ux: np.ndarray
    uy: np.ndarray
    temperature: np.ndarray | None
    dy: float

    @property
    def speed(self) -> np.ndarray:
        return np.hypot(self.ux, self.uy)

    def depth(self) -> np.ndarray:
        """Liquid depth h(x) = sum of alpha dy over each column [m]."""
        return self.alpha.sum(axis=0) * self.dy

    def mean_velocity(self) -> np.ndarray:
        """Liquid depth-averaged x-velocity [m/s]; NaN in dry columns."""
        h = self.depth()
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(h > 0, (self.alpha * self.ux).sum(axis=0) * self.dy / h, np.nan)

    def froude(self, g: float, min_depth: float) -> np.ndarray:
        """Froude number u/sqrt(g h); NaN where h < ``min_depth``."""
        h = self.depth()
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(h >= min_depth, self.mean_velocity() / np.sqrt(g * h), np.nan)


def _switch(value: str) -> bool:
    return value.lower() in {"yes", "on", "true", "1"}


def load_case(path: Path) -> CaseInfo:
    """Read the grid, the model switches and the initial state of a case."""
    grid = read_grid(path)
    n = grid.n_cells
    tp = read_text(path / "constant" / "transportProperties")

    try:
        step = int(read_entry(tp, "activeStep"))
    except KeyError:
        step = None

    alpha0 = grid.to_2d(np.clip(read_internal_field(field_path(path / "0", LIQUID), n), 0.0, 1.0))
    depth0 = alpha0.sum(axis=0) * grid.dy
    h0 = float(depth0.max())
    wet = depth0 > 0.5 * h0
    x_dam = float(grid.x[wet].max() + grid.dx / 2)

    temp0 = grid.to_2d(read_internal_field(field_path(path / "0", "T"), n))
    t_liquid = float(temp0[alpha0 > 0.5].mean())
    t_air = float(temp0[alpha0 < 0.5].mean())

    return CaseInfo(
        path=path,
        grid=grid,
        times=tuple(time_dirs(path)),
        g=read_gravity(path),
        step=step,
        solve_energy=_switch(read_entry(tp, "solveEnergy", block="energy")),
        temperature_dependent=_switch(
            read_entry(tp, "temperatureDependent", block="liquidViscosity")
        ),
        rho_liquid=float(read_entry(tp, "rho", block="water")),
        mu0_liquid=float(read_entry(tp, "mu0", block="liquidViscosity")),
        h0=h0,
        x_dam=x_dam,
        t_liquid=t_liquid,
        t_air=t_air,
    )


def load_frame(info: CaseInfo, time_name: str, temperature: bool = True) -> Frame:
    """Read the fields of one time directory."""
    grid = info.grid
    n = grid.n_cells
    tdir = info.path / time_name
    u = grid.to_2d(read_internal_field(field_path(tdir, "U"), n))
    return Frame(
        time=float(time_name),
        alpha=grid.to_2d(np.clip(read_internal_field(field_path(tdir, LIQUID), n), 0.0, 1.0)),
        ux=u[..., 0],
        uy=u[..., 1],
        temperature=grid.to_2d(read_internal_field(field_path(tdir, "T"), n)) if temperature else None,
        dy=grid.dy,
    )
