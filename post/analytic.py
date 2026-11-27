"""Analytic reference solutions for the dam break."""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np


class RitterSolution(NamedTuple):
    """Ritter depth and velocity; NaN where the solution is not valid.

    Attributes:
        depth: Water depth h(x) [m].
        velocity: Depth-averaged velocity u(x) [m/s].
        x_reflected: Head of the rarefaction reflected at the wall [m];
            equal to the wall position before the reflection.
    """

    depth: np.ndarray
    velocity: np.ndarray
    x_reflected: float


def ritter(x: np.ndarray, t: float, h0: float, x_wall: float, x_dam: float, g: float) -> RitterSolution:
    """Ritter (1892) dry-bed dam-break solution of the shallow-water equations.

    The reservoir of depth ``h0`` lies between a wall at ``x_wall`` and the dam
    at ``x_dam``. Ritter's solution assumes an infinite reservoir:

        h = (2 c0 - xi/t)^2 / (9 g),  u = 2/3 (xi/t + c0),  -c0 <= xi/t <= 2 c0

    with xi = x - x_dam and c0 = sqrt(g h0). The rarefaction reaches the wall
    at t0 = L/c0 (L = x_dam - x_wall); the reflected wave travels along the
    C+ characteristics, dxi/dt = u + c = xi/(3t) + 4 c0/3, so its head is at

        x_r(t) = x_dam + 2 c0 t - 3 L (t/t0)^(1/3)

    and the solution is not valid behind it (x < x_r).

    Args:
        x: Positions [m].
        t: Time after the dam break [s].
        h0: Initial reservoir depth [m].
        x_wall: Position of the upstream wall [m].
        x_dam: Initial position of the dam (reservoir front) [m].
        g: Gravitational acceleration [m/s^2].

    Returns:
        Depth, velocity and the head of the reflected wave.
    """
    length = x_dam - x_wall
    if t <= 0:
        depth = np.where((x >= x_wall) & (x <= x_dam), h0, 0.0)
        return RitterSolution(depth, np.zeros_like(x), x_wall)

    c0 = np.sqrt(g * h0)
    s = (x - x_dam) / t
    fan = (s >= -c0) & (s <= 2 * c0)
    depth = np.where(s < -c0, h0, np.where(fan, (2 * c0 - s) ** 2 / (9 * g), 0.0))
    velocity = np.where(fan, 2.0 / 3.0 * (s + c0), 0.0)

    t0 = length / c0
    x_reflected = x_dam + 2 * c0 * t - 3 * length * (t / t0) ** (1.0 / 3.0) if t > t0 else x_wall
    valid = x >= x_reflected
    return RitterSolution(
        np.where(valid, depth, np.nan),
        np.where(valid, velocity, np.nan),
        float(x_reflected),
    )


# Huppert (1982) 2D similarity constant: the volume condition on
# phi = (3/10)^(1/3) (eta_N^2 - xi^2)^(1/3) gives
# eta_N = ((10/3)^(1/3) / I)^(3/5) with I = int_0^1 (1 - s^2)^(1/3) ds
_PROFILE_INTEGRAL = math.sqrt(math.pi) * math.gamma(4.0 / 3.0) / (2.0 * math.gamma(11.0 / 6.0))
HUPPERT_ETA_N = ((10.0 / 3.0) ** (1.0 / 3.0) / _PROFILE_INTEGRAL) ** 0.6  # 1.411


class HuppertSolution(NamedTuple):
    """Depth profile and front of a 2D viscous gravity current.

    Attributes:
        depth: Depth h(x) [m]; zero beyond the front.
        x_front: Front position x_N, measured from the wall [m].
    """

    depth: np.ndarray
    x_front: float


def huppert_front(t: np.ndarray | float, area: float, g: float, nu: float) -> np.ndarray:
    """Front x_N(t) of a 2D viscous gravity current of fixed area [m].

    Huppert (1982), J. Fluid Mech. 121, 43-58:

        x_N = eta_N (g A^3 t / (3 nu))^(1/5),  eta_N = 1.411

    for a current of area ``area`` (A, per unit width) spreading from a wall
    at x = 0 under a fluid of negligible density, at low Reynolds number
    (lubrication approximation). It forgets the release, so it applies only
    once the current is long and thin compared with the initial column.
    """
    return HUPPERT_ETA_N * (g * area**3 * np.asarray(t, dtype=float) / (3.0 * nu)) ** 0.2


def huppert(x: np.ndarray, t: float, area: float, g: float, nu: float, x_wall: float = 0.0) -> HuppertSolution:
    """Huppert (1982) similarity solution h(x, t) of a 2D viscous current.

    Solves h_t = (g / (3 nu)) (h^3 h_x)_x with no flux at the wall:

        h = (3 nu A^2 / (g t))^(1/5) (3/10)^(1/3) (eta_N^2 - xi^2)^(1/3),
        xi = (x - x_wall) / (g A^3 t / (3 nu))^(1/5) <= eta_N

    Args:
        x: Positions [m].
        t: Time [s], > 0.
        area: Area of the current per unit width A [m^2].
        g: Gravitational acceleration [m/s^2].
        nu: Kinematic viscosity [m^2/s].
        x_wall: Position of the wall [m].
    """
    scale = (g * area**3 * t / (3.0 * nu)) ** 0.2
    xi = (np.asarray(x, dtype=float) - x_wall) / scale
    phi = (0.3 ** (1.0 / 3.0)) * np.cbrt(np.clip(HUPPERT_ETA_N**2 - xi**2, 0.0, None))
    depth = np.where(xi >= 0, area / scale * phi, np.nan)
    return HuppertSolution(depth, x_wall + HUPPERT_ETA_N * scale)
