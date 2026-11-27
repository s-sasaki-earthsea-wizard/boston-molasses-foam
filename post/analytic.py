"""Analytic reference solutions for the dam break."""

from __future__ import annotations

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
