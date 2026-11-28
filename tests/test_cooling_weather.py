"""Tests of the weather of the two days of the column runs."""

from __future__ import annotations

import pytest

from cooling.weather import SUMMER, WINTER_1919, fahrenheit_to_kelvin


def test_molasses_starts_5_k_above_the_air_at_12_40() -> None:
    # Winter: 43 F at noon, 44 F at 3 p.m. -> 43.2 F = 6.2 degC; summer: 79.6 F = 26.4 degC
    assert WINTER_1919.air(0.0) - 273.15 == pytest.approx(6.2, abs=0.05)
    assert SUMMER.air(0.0) - 273.15 == pytest.approx(26.4, abs=0.1)
    for day in (WINTER_1919, SUMMER):
        assert day.molasses_initial - day.air(0.0) == pytest.approx(5.0, abs=0.1)


def test_air_follows_the_readings_and_holds_after_midnight() -> None:
    hours = 3600.0
    # 18:00 is 5 h 20 min after 12:40
    assert WINTER_1919.air(5 * hours + 1200.0) == pytest.approx(float(fahrenheit_to_kelvin(39.0)))
    assert WINTER_1919.air(12 * hours) == pytest.approx(float(fahrenheit_to_kelvin(31.0)))
    assert SUMMER.air(12 * hours) == pytest.approx(float(fahrenheit_to_kelvin(70.6 + 0.0)), abs=0.4)


def test_sun_sets_within_the_runs() -> None:
    hours = 3600.0
    assert WINTER_1919.solar(0.0) > 300.0
    assert WINTER_1919.solar(5 * hours) == 0.0  # 17:40, after sunset
    assert SUMMER.solar(6 * hours) > 0.0 and SUMMER.solar(8 * hours) == 0.0  # 18:40 and 20:40


def test_humidity() -> None:
    assert 0.0 < WINTER_1919.vapour_pressure(0.0) < 1.0
    assert SUMMER.vapour_pressure(0.0) == pytest.approx(SUMMER.vapour_pressure(9 * 3600.0))
    assert SUMMER.vapour_pressure(0.0) == pytest.approx(1.85, abs=0.02)  # dew point 61.4 F = 16.3 degC
