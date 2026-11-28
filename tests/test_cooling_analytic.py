"""Checks of the closed-form conduction solutions against their own equations."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cooling.analytic import contact_temperature, surface_convection, surface_step, surface_step_flux
from cooling.properties import MOLASSES

K, A = MOLASSES.conductivity, MOLASSES.diffusivity
T_I, T_INF = 283.15, 273.15


def _depth_derivative(f, d: float, step: float) -> float:
    """dT/dd at ``d`` (one-sided, second order, so it also works at d = 0)."""
    return (-3.0 * f(d) + 4.0 * f(d + step) - f(d + 2.0 * step)) / (2.0 * step)


@pytest.mark.parametrize("t", [60.0, 3600.0, 43200.0])
@pytest.mark.parametrize("h", [2.0, 20.0, 200.0])
def test_convection_meets_its_boundary_condition(t: float, h: float) -> None:
    def f(d: float) -> float:
        return float(surface_convection(d, t, T_I, T_INF, h, K, A))

    gradient = _depth_derivative(f, 0.0, 1e-7)
    assert K * gradient == pytest.approx(h * (f(0.0) - T_INF), rel=1e-5)


@pytest.mark.parametrize("d", [0.002, 0.01, 0.03])
def test_convection_solves_the_heat_equation(d: float) -> None:
    t, h, dt, dd = 3600.0, 20.0, 1.0, 1e-4

    def f(depth: float, time: float) -> float:
        return float(surface_convection(depth, time, T_I, T_INF, h, K, A))

    dtemp_dt = (f(d, t + dt) - f(d, t - dt)) / (2.0 * dt)
    d2temp_dd2 = (f(d + dd, t) - 2.0 * f(d, t) + f(d - dd, t)) / dd**2
    assert dtemp_dt == pytest.approx(A * d2temp_dd2, rel=1e-3)


def test_convection_tends_to_the_step_for_large_h() -> None:
    d = np.linspace(0.0, 0.1, 51)
    strong = surface_convection(d, 3600.0, T_I, T_INF, 1e9, K, A)
    np.testing.assert_allclose(strong, surface_step(d, 3600.0, T_I, T_INF, A), atol=1e-5)


def test_convection_starts_at_the_initial_temperature_far_down() -> None:
    assert float(surface_convection(1.0, 3600.0, T_I, T_INF, 20.0, K, A)) == pytest.approx(T_I, abs=1e-9)


def test_step_profile_and_its_flux() -> None:
    t = 3600.0
    assert float(surface_step(0.0, t, T_I, T_INF, A)) == pytest.approx(T_INF)
    assert float(surface_step(1.0, t, T_I, T_INF, A)) == pytest.approx(T_I, abs=1e-9)

    def f(d: float) -> float:
        return float(surface_step(d, t, T_I, T_INF, A))

    # Heat leaves the warm solid: q_out = k dT/dd at the surface
    assert K * _depth_derivative(f, 0.0, 1e-7) == pytest.approx(surface_step_flux(t, T_I, T_INF, K, A), rel=1e-5)


def test_contact_temperature_weights_by_effusivity() -> None:
    assert contact_temperature(10.0, 1.0, 0.0, 1.0) == pytest.approx(5.0)
    assert contact_temperature(10.0, 953.0, 0.0, 2500.0) == pytest.approx(10.0 * 953.0 / 3453.0)
    assert contact_temperature(0.0, 2.0, 10.0, 1.0) == contact_temperature(10.0, 1.0, 0.0, 2.0)


@pytest.mark.parametrize("t", [0.0, -1.0])
def test_solutions_need_positive_time(t: float) -> None:
    with pytest.raises(ValueError, match="t > 0"):
        surface_step(0.0, t, T_I, T_INF, A)
    with pytest.raises(ValueError, match="t > 0"):
        surface_step_flux(t, T_I, T_INF, K, A)
    with pytest.raises(ValueError, match="t > 0"):
        surface_convection(0.0, t, T_I, T_INF, 20.0, K, A)


def test_convection_does_not_overflow() -> None:
    # exp(H d + H^2 a t) alone overflows for these values
    h, t, d = 1e4, 43200.0, 1.0
    assert (h / K) * d + (h / K) ** 2 * A * t > math.log(np.finfo(float).max)
    assert np.isfinite(surface_convection(d, t, T_I, T_INF, h, K, A))
