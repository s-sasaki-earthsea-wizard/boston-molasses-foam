"""Tests of the solar position and the clear-sky irradiance."""

from __future__ import annotations

import math
from datetime import date

import pytest

from cooling import sun

WINTER, SUMMER = date(1919, 1, 15), date(1919, 7, 15)


def _elevation(day: date, hour: float) -> float:
    return 90.0 - sun.zenith_angle(day, hour)


@pytest.mark.parametrize(
    ("day", "noon", "highest", "elevations"),
    [
        # Values on the wiki page Cooling-column-inputs, "Short-wave radiation"
        (WINTER, 11 + 53 / 60, 26.4, (25.4, 16.8, 9.1)),
        (SUMMER, 11 + 50 / 60, 69.3, (66.8, 50.8, 40.0)),
    ],
)
def test_sun_of_the_two_days(day: date, noon: float, highest: float, elevations: tuple[float, ...]) -> None:
    assert sun.solar_noon(day) == pytest.approx(noon, abs=0.5 / 60)
    assert _elevation(day, sun.solar_noon(day)) == pytest.approx(highest, abs=0.05)
    for hour, expected in zip((12 + 40 / 60, 14.5, 15.5), elevations, strict=True):
        assert _elevation(day, hour) == pytest.approx(expected, abs=0.05)


def test_sunset_of_15_january_1919() -> None:
    # 16:35 EST computed for the wiki; the Boston American printed "sets, 4:37"
    assert sun.sunset(WINTER) == pytest.approx(16 + 35 / 60, abs=1.0 / 60)
    assert sun.zenith_angle(WINTER, sun.sunset(WINTER)) == pytest.approx(sun.SUNSET_ZENITH, abs=0.1)


def test_noon_is_the_highest_sun_and_hours_wrap_to_the_next_day() -> None:
    noon = sun.solar_noon(WINTER)
    assert sun.zenith_angle(WINTER, noon) < min(sun.zenith_angle(WINTER, noon - 0.2), sun.zenith_angle(WINTER, noon + 0.2))
    nxt = date(1919, 1, 16)
    assert sun.zenith_angle(WINTER, 24.0 + 13.0) == pytest.approx(sun.zenith_angle(nxt, 13.0), abs=1e-9)


def test_haurwitz_and_the_cloud_factor() -> None:
    assert sun.haurwitz_ghi(0.0) == pytest.approx(1098.0 * math.exp(-0.057))
    assert sun.haurwitz_ghi(90.0) == 0.0 and sun.haurwitz_ghi(120.0) == 0.0
    assert sun.cloud_factor(0.0) == 1.0
    assert sun.cloud_factor(1.0) == pytest.approx(0.25)
    with pytest.raises(ValueError, match="cloud"):
        sun.cloud_factor(1.5)
