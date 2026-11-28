"""Tests of the surface heat balance: air properties, convection, sky and the boundary form."""

from __future__ import annotations

import numpy as np
import pytest

from cooling.column import Adiabatic, Column, Layer
from cooling.properties import MOLASSES
from cooling.surface import (
    SIGMA,
    BulkConvection,
    ClearConvection,
    McAdamsConvection,
    SurfaceBalance,
    air_properties,
    clear_sky_emissivity,
    natural_convection,
    saturation_vapour_pressure,
    sky_emissivity,
    wind_at_height,
)

T_AIR, T_SURFACE = 279.35, 284.35  # winter start: air 6.2 degC, molasses 11.2 degC


def test_air_properties_hit_the_table_and_refuse_outside() -> None:
    p = air_properties(280.0)
    assert (p.density, p.conductivity, p.viscosity, p.prandtl) == (1.261, 0.0249, 1.392e-5, 0.710)
    assert air_properties(285.0).conductivity == pytest.approx(0.02525)
    with pytest.raises(ValueError, match="outside"):
        air_properties(240.0)


def test_natural_convection() -> None:
    # 0.15 Ra^(1/3): the length cancels for a warm surface
    assert natural_convection(T_SURFACE, T_AIR, 1.0) == pytest.approx(natural_convection(T_SURFACE, T_AIR, 50.0))
    assert natural_convection(T_SURFACE, T_AIR, 1.0) == pytest.approx(3.22, abs=0.01)
    assert natural_convection(T_AIR, T_AIR, 1.0) == 0.0
    # a cold surface (stable) transfers less, and less over a larger surface
    cold = natural_convection(T_AIR - 5.0, T_AIR, 7.5)
    assert 0.0 < cold < natural_convection(T_SURFACE, T_AIR, 7.5)
    assert natural_convection(T_AIR - 5.0, T_AIR, 30.0) < cold


def test_clear_forced_part_averages_to_the_flat_plate_mean() -> None:
    # With natural convection damped away, the mean of 0.0296 Re_x^0.8 Pr^(1/3) k/x over
    # the fetch is 0.037 Re_L^0.8 Pr^(1/3) k/L.
    fetch, wind = 30.0, 50.0
    p = air_properties(0.5 * (T_SURFACE + T_AIR))
    re = wind * fetch / p.viscosity
    expected = 0.037 * re**0.8 * p.prandtl ** (1 / 3) * p.conductivity / fetch
    h = ClearConvection(fetch).coefficient(wind, T_SURFACE, T_AIR)
    h_n = natural_convection(T_SURFACE, T_AIR, fetch / 4.0)
    assert expected < h < expected + 0.01 * h_n
    assert h == pytest.approx(expected, rel=1e-3)


def test_clear_convection_limits() -> None:
    clear = ClearConvection(30.0)
    assert clear.coefficient(0.0, T_SURFACE, T_AIR) == pytest.approx(natural_convection(T_SURFACE, T_AIR, 7.5))
    values = [clear.coefficient(v, T_SURFACE, T_AIR) for v in (0.5, 1.0, 2.0, 4.0)]
    assert np.all(np.diff(values) > 0.0)
    # a longer fetch, a smaller mean coefficient; a rough surface, a larger one
    assert ClearConvection(100.0).coefficient(2.0, T_SURFACE, T_AIR) < clear.coefficient(2.0, T_SURFACE, T_AIR)
    assert ClearConvection(30.0, roughness=1.67).coefficient(2.0, T_SURFACE, T_AIR) > clear.coefficient(2.0, T_SURFACE, T_AIR)


def test_bulk_and_mcadams() -> None:
    p = air_properties(T_AIR)
    m10 = 3.0 * (10.0 / 3.0) ** 0.33
    assert BulkConvection().coefficient(3.0, T_SURFACE, T_AIR) == pytest.approx(p.density * p.heat_capacity * 2e-3 * m10)
    assert McAdamsConvection().coefficient(3.0, T_SURFACE, T_AIR) == pytest.approx(17.1)
    assert wind_at_height(4.8, 57.0, 3.0) == pytest.approx(1.82, abs=0.01)


def test_sky_emissivities() -> None:
    assert clear_sky_emissivity("swinbank", 280.0) * SIGMA * 280.0**4 == pytest.approx(5.31e-13 * 280.0**6)
    assert clear_sky_emissivity("idso_jackson", 273.16) == pytest.approx(0.739)
    e = 0.7 * saturation_vapour_pressure(280.0)
    # the three clear-sky forms agree within a few hundredths at 7 degC and 70 %
    values = [clear_sky_emissivity("swinbank", 280.0), clear_sky_emissivity("idso_jackson", 280.0),
              clear_sky_emissivity("prata", 280.0, e)]
    assert max(values) - min(values) < 0.03
    with pytest.raises(ValueError, match="vapour"):
        clear_sky_emissivity("prata", 280.0)
    assert sky_emissivity(0.7, 0.0) == 0.7
    assert sky_emissivity(0.7, 1.0) == pytest.approx(0.16 * 0.7 + 0.84)


def test_saturation_vapour_pressure_fao56_example() -> None:
    # FAO-56, chapter 3, example: e(24.5 degC) = 3.075 kPa, e(15 degC) = 1.705 kPa
    assert saturation_vapour_pressure(273.15 + 24.5) == pytest.approx(3.075, abs=5e-4)
    assert saturation_vapour_pressure(273.15 + 15.0) == pytest.approx(1.705, abs=5e-4)


def _balance(convection, **kwargs) -> SurfaceBalance:
    return SurfaceBalance(
        air=lambda t: T_AIR - 1e-4 * t,
        vapour_pressure=lambda t: 0.7,
        solar=lambda t: 400.0,
        wind=2.0,
        convection=convection,
        **kwargs,
    )


@pytest.mark.parametrize("convection", [ClearConvection(30.0), BulkConvection(), McAdamsConvection()])
@pytest.mark.parametrize("t_s", [270.0, 284.35, 300.0])
def test_boundary_form_reproduces_the_balance(convection, t_s: float) -> None:
    balance = _balance(convection, cloud=0.4, sky_model="prata")
    h, t_ref = balance.coefficients(600.0, t_s)
    assert h * (t_s - t_ref) == pytest.approx(balance.terms(600.0, t_s).net_loss, rel=1e-9, abs=1e-9)


def test_mcadams_adds_only_the_sky_deficit() -> None:
    balance = _balance(McAdamsConvection())
    terms = balance.terms(0.0, T_SURFACE)
    assert terms.longwave == pytest.approx(0.95 * (SIGMA * T_AIR**4 - terms.sky))
    assert terms.convection == pytest.approx(13.3 * (T_SURFACE - T_AIR))  # 5.7 + 3.8 x 2 m/s


def test_column_with_the_surface_balance_closes_its_budget() -> None:
    column = Column([Layer(MOLASSES, 0.1, T_SURFACE)], _balance(ClearConvection(30.0)), Adiabatic())
    result = column.run(7200.0, dt=30.0, output_times=[1800.0, 3600.0])
    change = result.heat_content - result.heat_content[0]
    np.testing.assert_allclose(change, result.heat_in, rtol=0.0, atol=1e-6 * abs(result.heat_in[-1]))
    # with 400 W/m2 of sun the surface warms although the air is colder
    assert result.surface_temperature[-1] > T_SURFACE
