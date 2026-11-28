"""Air, humidity and sun of the two days of the column runs, from the collapse at 12:40 EST.

Values and sources: wiki page Cooling-column-inputs.

- 15 January 1919: air from the readings at Thompson's Spa printed in the
  Boston Daily Globe of 16 January 1919, p. 1 (3-hourly, deg F); relative
  humidity from the monthly mean of January 1919 (Monthly Weather Review
  47(1), Table I).
- 15 July: hourly normals 1991-2020 of Boston Logan (USW00014739), in local
  standard time; dew point from the normals at 12 and 13 LST.
- Sun: clear-sky irradiance from the solar position of the day (cooling.sun).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

import numpy as np

from cooling import sun
from cooling.surface import saturation_vapour_pressure

COLLAPSE_HOUR = 12.0 + 40.0 / 60.0  # 12:40 EST, in decimal hours of the clock


def fahrenheit_to_kelvin(value: float | np.ndarray) -> float | np.ndarray:
    return (np.asarray(value, dtype=float) - 32.0) / 1.8 + 273.15


@dataclass(frozen=True)
class Day:
    """Weather of one day of the runs.

    Attributes:
        name: Label.
        day: Date (the clock is Eastern Standard Time).
        hours: Clock hours of the air readings (may pass 24) [h].
        air_f: Air temperature at ``hours`` [deg F]; held constant beyond the ends.
        molasses_initial: Molasses at the collapse [K]: air at 12:40 + 5 K, as decided in #9.
        ground_initial: Uniform ground at the collapse [K], base value.
        ground_range: Range of the ground temperature for the sensitivity runs [K].
        relative_humidity: Relative humidity (0-1), if the humidity is given this way.
        dew_point_f: Dew point [deg F], if the humidity is given this way.
    """

    name: str
    day: date
    hours: tuple[float, ...]
    air_f: tuple[float, ...]
    molasses_initial: float
    ground_initial: float
    ground_range: tuple[float, float]
    relative_humidity: float | None = None
    dew_point_f: float | None = None
    site: sun.Site = field(default_factory=sun.Site)

    def clock(self, time: float) -> float:
        """Clock hour at ``time`` seconds after the collapse."""
        return COLLAPSE_HOUR + time / 3600.0

    def air(self, time: float) -> float:
        """Air temperature [K] at ``time`` seconds after the collapse (linear between readings)."""
        return float(fahrenheit_to_kelvin(np.interp(self.clock(time), self.hours, self.air_f)))

    def vapour_pressure(self, time: float) -> float:
        """Vapour pressure [kPa]."""
        if self.dew_point_f is not None:
            return saturation_vapour_pressure(float(fahrenheit_to_kelvin(self.dew_point_f)))
        if self.relative_humidity is None:
            raise ValueError(f"{self.name}: no humidity given")
        return self.relative_humidity * saturation_vapour_pressure(self.air(time))

    def zenith(self, time: float) -> float:
        """Solar zenith angle [deg]."""
        return sun.zenith_angle(self.day, self.clock(time), self.site)

    def solar(self, time: float) -> float:
        """Clear-sky global horizontal irradiance [W/m2]."""
        return sun.haurwitz_ghi(self.zenith(time))


# Thompson's Spa, Boston Daily Globe, 16 Jan 1919, p. 1: 3 a.m. to midnight
WINTER_1919 = Day(
    name="15 January 1919",
    day=date(1919, 1, 15),
    hours=(3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0),
    air_f=(42.0, 41.0, 41.0, 43.0, 44.0, 39.0, 34.0, 31.0),
    molasses_initial=273.15 + 11.2,
    ground_initial=273.15 + 1.0,
    ground_range=(273.15 - 0.5, 273.15 + 5.0),
    relative_humidity=0.70,
)

# Hourly normals of 15 July, 00-23 LST; the value of 00 LST closes the day at 24 h
_JULY_15_F = (70.6, 70.0, 69.5, 69.0, 68.5, 68.4, 69.8, 71.9, 74.1, 76.0, 77.5, 78.6,
              79.3, 79.8, 79.8, 79.5, 79.1, 78.2, 77.0, 75.5, 74.1, 73.1, 72.2, 71.4)
SUMMER = Day(
    name="15 July",
    day=date(1919, 7, 15),
    hours=tuple(float(h) for h in range(25)),
    air_f=_JULY_15_F + (_JULY_15_F[0],),
    molasses_initial=273.15 + 31.4,
    ground_initial=273.15 + 24.0,
    ground_range=(273.15 + 21.0, 273.15 + 26.0),
    dew_point_f=61.4,
)

DAYS = {"winter": WINTER_1919, "summer": SUMMER}
