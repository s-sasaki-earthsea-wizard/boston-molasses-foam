"""Tests of the 1D column solver against closed-form solutions and its heat budget."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cooling.analytic import contact_temperature, surface_convection, surface_step, surface_step_flux
from cooling.column import Adiabatic, Column, Convective, FixedTemperature, Layer, graded_widths
from cooling.properties import MOLASSES, Material

# A stone-like ground for the tests only; the ground of the runs is chosen with sources
ROCK = Material("rock", conductivity=2.5, density=2600.0, heat_capacity=850.0)
T_I, T_COLD = 283.15, 273.15
SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant [W/m2/K4]
TIMES = [600.0, 3600.0, 21600.0]
# Largest temperature errors allowed at TIMES with dt = 10 s [K]; implicit Euler
# errs most early on. Measured: 0.025 / 0.005 / 0.002 (step) and
# 0.008 / 0.0015 / 0.0012 (convection), 5e-4 at the contact.
STEP_ATOL = [0.04, 0.008, 0.004]
CONVECTION_ATOL = [0.015, 0.003, 0.003]
CONTACT_ATOL = 0.002


class Radiation:
    """Long-wave exchange e sigma (T_b^4 - T_sky^4), written as h(T_b) (T_b - T_sky)."""

    def __init__(self, emissivity: float, sky: float) -> None:
        self.emissivity, self.sky = emissivity, sky

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        t, s = boundary_temperature, self.sky
        return self.emissivity * SIGMA * (t * t + s * s) * (t + s), s


# --- grid -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("bottom", "top"), [(True, True), (False, True), (True, False), (False, False)]
)
@pytest.mark.parametrize("length", [0.02, 0.3, 1.0])
def test_graded_widths_cover_the_length(bottom: bool, top: bool, length: float) -> None:
    w = graded_widths(length, 5e-4, 0.02, 1.1, refine_bottom=bottom, refine_top=top)
    assert w.sum() == pytest.approx(length, rel=1e-12)
    assert np.all(w > 0.0) and w.max() <= 0.02 * (1 + 1e-12)
    ratio = w[1:] / w[:-1]
    assert np.all(ratio <= 1.1 * (1 + 1e-12)) and np.all(ratio >= 1 / 1.1 * (1 - 1e-12))
    if top:
        assert w[-1] <= 5e-4 * (1 + 1e-12)
    if bottom:
        assert w[0] <= 5e-4 * (1 + 1e-12)


def test_graded_widths_reject_bad_input() -> None:
    with pytest.raises(ValueError, match="length"):
        graded_widths(0.0, 5e-4, 0.02, 1.1, True, True)
    with pytest.raises(ValueError, match="dz_min"):
        graded_widths(1.0, 0.03, 0.02, 1.1, True, True)
    with pytest.raises(ValueError, match="growth"):
        graded_widths(1.0, 5e-4, 0.02, 0.9, True, True)


def test_column_refines_the_surface_and_the_interface() -> None:
    column = Column([Layer(ROCK, 1.0, T_COLD), Layer(MOLASSES, 0.3, T_I)], Adiabatic(), Adiabatic())
    assert column.widths[-1] <= 5e-4 * (1 + 1e-12)
    below = np.flatnonzero(np.diff(column.layer_index))[0]
    assert column.widths[below] <= 5e-4 * (1 + 1e-12) and column.widths[below + 1] <= 5e-4 * (1 + 1e-12)
    assert column.widths[0] == pytest.approx(0.02, rel=0.2)  # far field, coarse
    assert column.interfaces == pytest.approx([1.0])


# --- closed-form solutions ----------------------------------------------------


def _molasses_profile(result, i: int, max_depth: float):
    depth, cells = result.depth_in_layer(len(result.interfaces))
    keep = depth < max_depth
    return depth[keep], result.temperature[i, cells[keep]]


def test_surface_step_on_a_half_space() -> None:
    column = Column([Layer(MOLASSES, 0.5, T_I)], FixedTemperature(T_COLD), Adiabatic())
    result = column.run(TIMES[-1], dt=10.0, output_times=TIMES)
    for i, t in enumerate(result.times[1:], start=1):
        depth, temp = _molasses_profile(result, i, 0.3)
        expected = surface_step(depth, t, T_I, T_COLD, MOLASSES.diffusivity)
        np.testing.assert_allclose(temp, expected, atol=STEP_ATOL[i - 1])
        assert result.surface_temperature[i] == T_COLD
        q = surface_step_flux(t, T_I, T_COLD, MOLASSES.conductivity, MOLASSES.diffusivity)
        assert result.surface_flux[i] == pytest.approx(q, rel=0.01)  # measured 0.6 % at 600 s


def test_convection_on_a_half_space() -> None:
    h = 20.0
    column = Column([Layer(MOLASSES, 0.5, T_I)], Convective(h, T_COLD), Adiabatic())
    result = column.run(TIMES[-1], dt=10.0, output_times=TIMES)
    k, a = MOLASSES.conductivity, MOLASSES.diffusivity
    for i, t in enumerate(result.times[1:], start=1):
        depth, temp = _molasses_profile(result, i, 0.3)
        atol = CONVECTION_ATOL[i - 1]
        np.testing.assert_allclose(temp, surface_convection(depth, t, T_I, T_COLD, h, k, a), atol=atol)
        expected = float(surface_convection(0.0, t, T_I, T_COLD, h, k, a))
        assert result.surface_temperature[i] == pytest.approx(expected, abs=atol)
        assert result.surface_flux[i] == pytest.approx(h * (result.surface_temperature[i] - T_COLD), rel=1e-9)


def test_two_half_spaces_in_contact() -> None:
    column = Column([Layer(ROCK, 2.0, T_COLD), Layer(MOLASSES, 0.5, T_I)], Adiabatic(), Adiabatic())
    result = column.run(TIMES[-1], dt=10.0, output_times=TIMES)
    t_c = contact_temperature(T_I, MOLASSES.effusivity, T_COLD, ROCK.effusivity)
    interface = result.interfaces[0]
    for i, t in enumerate(result.times[1:], start=1):
        assert result.interface_temperature[i, 0] == pytest.approx(t_c, abs=CONTACT_ATOL)
        above = (result.layer_index == 1) & (result.centres < interface + 0.3)
        np.testing.assert_allclose(
            result.temperature[i, above],
            surface_step(result.centres[above] - interface, t, T_I, t_c, MOLASSES.diffusivity),
            atol=STEP_ATOL[i - 1],
        )
        below = (result.layer_index == 0) & (result.centres > interface - 1.0)
        np.testing.assert_allclose(
            result.temperature[i, below],
            surface_step(interface - result.centres[below], t, T_COLD, t_c, ROCK.diffusivity),
            atol=STEP_ATOL[i - 1],
        )


def test_error_falls_with_the_time_step() -> None:
    def error(dt: float) -> float:
        column = Column([Layer(MOLASSES, 0.5, T_I)], FixedTemperature(T_COLD), Adiabatic())
        result = column.run(3600.0, dt=dt)
        depth, temp = _molasses_profile(result, -1, 0.3)
        return float(np.max(np.abs(temp - surface_step(depth, 3600.0, T_I, T_COLD, MOLASSES.diffusivity))))

    # Implicit Euler is first order: four times the step, four times the error
    # (measured 0.0459 K at 120 s, 0.0117 K at 30 s)
    assert 3.5 < error(120.0) / error(30.0) < 4.5


# --- heat budget and nonlinear boundaries -----------------------------------


def test_heat_budget_closes() -> None:
    def ambient(t: float) -> float:
        return T_COLD + 5.0 * math.sin(2.0 * math.pi * t / 86400.0)

    column = Column(
        [Layer(ROCK, 1.0, 275.15), Layer(MOLASSES, 0.2, T_I)],
        surface=Convective(15.0, ambient),
        bottom=FixedTemperature(277.15),
    )
    result = column.run(43200.0, dt=60.0, output_times=[3600.0 * i for i in range(1, 12)])
    change = result.heat_content - result.heat_content[0]
    assert abs(result.heat_in[-1]) > 1e5  # the test is not trivial
    np.testing.assert_allclose(change, result.heat_in, rtol=0.0, atol=1e-6 * abs(result.heat_in[-1]))


def test_radiation_iterates_to_the_nonlinear_flux() -> None:
    sky = 250.0
    column = Column([Layer(MOLASSES, 0.3, T_I)], Radiation(0.95, sky), Adiabatic())
    result = column.run(21600.0, dt=60.0, output_times=[3600.0])
    t_s = result.surface_temperature[1:]
    np.testing.assert_allclose(result.surface_flux[1:], 0.95 * SIGMA * (t_s**4 - sky**4), rtol=1e-6)
    np.testing.assert_allclose(
        result.heat_content - result.heat_content[0], result.heat_in, rtol=0.0, atol=1e-6 * abs(result.heat_in[-1])
    )
    assert np.all(np.diff(t_s) < 0.0) and t_s[-1] > sky


def test_output_times_are_hit_exactly() -> None:
    column = Column([Layer(MOLASSES, 0.1, T_I)], Convective(10.0, T_COLD), Adiabatic())
    result = column.run(1000.0, dt=300.0, output_times=[250.0, 999.0, 0.0, 5000.0])
    assert result.times.tolist() == [0.0, 250.0, 999.0, 1000.0]
    assert result.temperature.shape == (4, column.widths.size)
    assert result.interface_temperature.shape == (4, 0)


def test_layer_and_run_validation() -> None:
    with pytest.raises(ValueError, match="thickness"):
        Layer(MOLASSES, 0.0, T_I)
    with pytest.raises(ValueError, match="kelvin"):
        Layer(MOLASSES, 0.1, -1.0)
    with pytest.raises(ValueError, match="at least one layer"):
        Column([], Adiabatic(), Adiabatic())
    column = Column([Layer(MOLASSES, 0.1, T_I)], Adiabatic(), Adiabatic())
    with pytest.raises(ValueError, match="positive"):
        column.run(100.0, dt=0.0)
