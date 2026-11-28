"""Position of the sun and clear-sky global irradiance on a horizontal surface.

The solar position follows NOAA Global Monitoring Laboratory, "General Solar
Position Calculations" (fractional year, equation of time, declination, hour
angle, zenith). The clear-sky irradiance is Haurwitz (1945) as given by
Reno, Hansen & Stein (2012, eq. 18), and the cloud factor is Kasten &
Czeplak (1980) as given by Flerchinger et al. (2009, eq. 5). Sources and
values: wiki page Cooling-column-inputs, section "Short-wave radiation".
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import date

# North End waterfront, Boston (approximate site of the tank)
BOSTON_LATITUDE = 42.3683  # [deg N]
BOSTON_LONGITUDE = -71.0560  # [deg, east positive]
EASTERN_STANDARD_TIME = -5.0  # [h from UTC]
SUNSET_ZENITH = 90.833  # [deg], refraction and solar disk at sunrise and sunset


@dataclass(frozen=True)
class Site:
    """Where and on which clock the sun is computed.

    Attributes:
        latitude: Latitude [deg, north positive].
        longitude: Longitude [deg, east positive].
        timezone: Offset of the clock from UTC [h].
    """

    latitude: float = BOSTON_LATITUDE
    longitude: float = BOSTON_LONGITUDE
    timezone: float = EASTERN_STANDARD_TIME


def _fractional_year(day: date, hour: float) -> float:
    """Fractional year gamma [rad] at ``hour`` (local clock) of ``day``."""
    days_in_year = 366 if (day.year % 4 == 0 and day.year % 100 != 0) or day.year % 400 == 0 else 365
    return 2.0 * math.pi / days_in_year * (day.timetuple().tm_yday - 1 + (hour - 12.0) / 24.0)


def equation_of_time(gamma: float) -> float:
    """Equation of time [min] at the fractional year ``gamma`` [rad]."""
    return 229.18 * (
        0.000075
        + 0.001868 * math.cos(gamma)
        - 0.032077 * math.sin(gamma)
        - 0.014615 * math.cos(2.0 * gamma)
        - 0.040849 * math.sin(2.0 * gamma)
    )


def declination(gamma: float) -> float:
    """Solar declination [rad] at the fractional year ``gamma`` [rad]."""
    return (
        0.006918
        - 0.399912 * math.cos(gamma)
        + 0.070257 * math.sin(gamma)
        - 0.006758 * math.cos(2.0 * gamma)
        + 0.000907 * math.sin(2.0 * gamma)
        - 0.002697 * math.cos(3.0 * gamma)
        + 0.00148 * math.sin(3.0 * gamma)
    )


def zenith_angle(day: date, hour: float, site: Site = Site()) -> float:
    """Solar zenith angle [deg] at ``hour`` (decimal hours on the clock of ``site``) of ``day``.

    ``hour`` may run past 24 into the next day.
    """
    day_offset, clock = divmod(hour, 24.0)
    day = date.fromordinal(day.toordinal() + int(day_offset))
    gamma = _fractional_year(day, clock)
    offset = equation_of_time(gamma) + 4.0 * site.longitude - 60.0 * site.timezone  # [min]
    hour_angle = math.radians((clock * 60.0 + offset) / 4.0 - 180.0)
    lat, decl = math.radians(site.latitude), declination(gamma)
    cos_z = math.sin(lat) * math.sin(decl) + math.cos(lat) * math.cos(decl) * math.cos(hour_angle)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_z))))


def solar_noon(day: date, site: Site = Site()) -> float:
    """Solar noon [decimal hours on the clock of ``site``]."""
    gamma = _fractional_year(day, 12.0)
    return (720.0 - 4.0 * site.longitude - equation_of_time(gamma)) / 60.0 + site.timezone


def sunset(day: date, site: Site = Site()) -> float:
    """Sunset [decimal hours on the clock of ``site``], zenith 90.833 deg."""
    gamma = _fractional_year(day, 12.0)
    lat, decl = math.radians(site.latitude), declination(gamma)
    cos_ha = math.cos(math.radians(SUNSET_ZENITH)) / (math.cos(lat) * math.cos(decl)) - math.tan(lat) * math.tan(decl)
    hour_angle = -math.degrees(math.acos(cos_ha))  # negative: sunset
    return (720.0 - 4.0 * (site.longitude + hour_angle) - equation_of_time(gamma)) / 60.0 + site.timezone


def haurwitz_ghi(zenith: float) -> float:
    """Clear-sky global horizontal irradiance [W/m2], Haurwitz: 1098 cos z exp(-0.057 / cos z).

    Zero with the sun at or below the horizon.
    """
    cos_z = math.cos(math.radians(zenith))
    if cos_z <= 0.0:
        return 0.0
    return 1098.0 * cos_z * math.exp(-0.057 / cos_z)


def cloud_factor(cloud: float, b1: float = 0.75, b2: float = 3.4) -> float:
    """Kasten & Czeplak transmission under a cloud fraction ``cloud`` (0-1): 1 - b1 c^b2."""
    if not 0.0 <= cloud <= 1.0:
        raise ValueError(f"cloud fraction must be in [0, 1] (got {cloud})")
    return 1.0 - b1 * cloud**b2
