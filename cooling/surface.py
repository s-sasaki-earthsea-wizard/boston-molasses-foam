"""Heat exchange of the free surface of the molasses with the air, the sky and the sun.

The surface loses

    q_out = h_c (T_s - T_a) + eps sigma T_s^4 - eps L_down - (1 - albedo) S_down

with the convective coefficient h_c, the downward long-wave radiation of the
sky L_down and the global short-wave irradiance S_down. ``SurfaceBalance``
writes this as h (T_s - T_ref) for ``cooling.column.Column``; h depends on
T_s through the radiation and the natural convection, which the column
iterates within each step.

Choices and sources: wiki pages Cooling-column-inputs and
Why-convection-coefficients-differ.

- Convection, base: Clear, Gartland & Winkelmann (2001, LBNL-47275), fitted
  on horizontal roofs of 2,400-2,900 m2 outdoors: flat-plate natural
  convection, damped by eta when the wind dominates, plus turbulent
  flat-plate forced convection from the upwind edge (their Table 8),
  averaged over the fetch.
- Convection, alternatives: the bulk transfer form of meteorology (Stull
  2015, eq. 3.35) and McAdams (5.7 + 3.8 v, as cast by Palyvos 2008), whose
  coefficient already contains the radiation to surroundings at air
  temperature and is used only as an upper bound.
- Sky: Swinbank (1963), Idso & Jackson (1969), Prata (1996), with the cloud
  correction of Unsworth & Monteith (1975), as given by Flerchinger et al.
  (2009).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol

import numpy as np

SIGMA = 5.670374419e-8  # Stefan-Boltzmann constant [W/m2/K4]
GRAVITY = 9.80665  # [m/s2]
TOWN_WIND_EXPONENT = 0.33  # power-law exponent of the wind profile, towns and cities (EnergyPlus Table 3.4)

# Air at atmospheric pressure, Lienhard & Lienhard (2024), Table A.6:
# T [K], rho [kg/m3], cp [J/kg/K], nu [m2/s], k [W/m/K], alpha [m2/s], Pr
_AIR = np.array(
    [
        [250, 1.413, 1006, 1.135e-5, 0.0226, 1.59e-5, 0.715],
        [260, 1.359, 1006, 1.218e-5, 0.0233, 1.71e-5, 0.713],
        [270, 1.308, 1006, 1.304e-5, 0.0241, 1.83e-5, 0.711],
        [280, 1.261, 1006, 1.392e-5, 0.0249, 1.96e-5, 0.710],
        [290, 1.218, 1006, 1.482e-5, 0.0256, 2.09e-5, 0.708],
        [300, 1.177, 1006, 1.575e-5, 0.0264, 2.23e-5, 0.707],
        [310, 1.139, 1007, 1.670e-5, 0.0271, 2.37e-5, 0.706],
        [320, 1.103, 1007, 1.766e-5, 0.0279, 2.51e-5, 0.705],
        [330, 1.070, 1008, 1.865e-5, 0.0286, 2.65e-5, 0.704],
        [340, 1.038, 1009, 1.966e-5, 0.0293, 2.80e-5, 0.703],
        [350, 1.009, 1009, 2.069e-5, 0.0300, 2.95e-5, 0.702],
    ]
)


@dataclass(frozen=True)
class AirProperties:
    """Properties of air at one temperature (linear in the table)."""

    density: float
    heat_capacity: float
    viscosity: float  # kinematic [m2/s]
    conductivity: float
    diffusivity: float
    prandtl: float


def air_properties(temperature: float) -> AirProperties:
    """Air at ``temperature`` [K], interpolated in Lienhard's Table A.6 (250-350 K)."""
    t = _AIR[:, 0]
    if not t[0] <= temperature <= t[-1]:
        raise ValueError(f"air temperature {temperature:.1f} K outside the table ({t[0]:.0f}-{t[-1]:.0f} K)")
    values = [float(np.interp(temperature, t, _AIR[:, j])) for j in range(1, 7)]
    return AirProperties(*values)


def wind_at_height(speed: float, height: float, to_height: float, exponent: float = TOWN_WIND_EXPONENT) -> float:
    """Wind speed at ``to_height`` from one at ``height`` on a power-law profile in the same terrain."""
    return speed * (to_height / height) ** exponent


# --- convection -------------------------------------------------------------


class Convection(Protocol):
    """A convective heat transfer coefficient."""

    includes_radiation: bool

    def coefficient(self, wind: float, surface: float, air: float) -> float:
        """h_c [W/m2/K] for the wind [m/s] and the surface and air temperatures [K]."""
        ...


def natural_convection(surface: float, air: float, length: float) -> float:
    """Natural convection from the top of a horizontal surface [W/m2/K].

    Surface warmer than the air: Nu = 0.15 Ra^(1/3) (turbulent, the length
    cancels); colder: Nu = 0.27 Ra_L^(1/4) with L = area / perimeter
    (Clear et al. 2001, Table 5 and eqs. 5a-5b). Gas properties at
    T_w - 0.83 (T_w - T_inf) for the heated case (Lienhard & Lienhard 2024,
    eq. 8.36), at the mean of the two otherwise; beta = 1 / T_inf.
    """
    dt = surface - air
    if dt == 0.0:
        return 0.0
    t_ref = surface - 0.83 * dt if dt > 0.0 else 0.5 * (surface + air)
    p = air_properties(t_ref)
    gb = GRAVITY * abs(dt) / air / (p.viscosity * p.diffusivity)  # Ra / L^3
    if dt > 0.0:
        return 0.15 * p.conductivity * gb ** (1.0 / 3.0)
    return 0.27 * p.conductivity * (gb * length**3) ** 0.25 / length


def _eta(gr_over_re2: float) -> float:
    """Damping of natural convection by the wind, ln(1 + G) / (1 + ln(1 + G)) (Clear et al. 2001, eq. 4b)."""
    if math.isinf(gr_over_re2):
        return 1.0
    log = math.log1p(gr_over_re2)
    return log / (1.0 + log)


# Gauss-Legendre nodes on [0, 1] for the averages over the fetch
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(48)
_NODES, _WEIGHTS = 0.5 * (_NODES + 1.0), 0.5 * _WEIGHTS


@dataclass(frozen=True)
class ClearConvection:
    """Natural plus forced convection over a horizontal surface, after Clear et al. (2001), Table 8.

    The local coefficient at a distance x from the upwind edge is

        h_x = eta(x) h_n + R_f (k / x) Nu_f,x

    with Nu_f,x = 0.0296 Re_x^(4/5) Pr^(1/3) where the flow is turbulent and
    0.332 Re_x^(1/2) Pr^(1/3) where it is laminar. For a surface warmer than
    the air there is no laminar part (Re_x,c ~ 0, Clear et al. 2001, p. 19);
    for a colder surface laminar flow holds up to Re_x,c = 5e5. The
    coefficient returned is the average of h_x over the fetch: the forced
    part in closed form on each side of the transition, so that it varies
    smoothly with the transition point, and the natural part by quadrature.

    Args:
        fetch: Length of the surface along the wind [m].
        roughness: Multiplier R_f of the forced part (1.0 for glass and smooth surfaces).
        critical_reynolds: Re_x,c for a surface colder than the air.
    """

    fetch: float
    roughness: float = 1.0
    critical_reynolds: float = 5e5
    includes_radiation: bool = False

    def coefficient(self, wind: float, surface: float, air: float) -> float:
        if not self.fetch > 0.0:
            raise ValueError(f"fetch must be positive (got {self.fetch})")
        dt = surface - air
        length_n = self.fetch / 4.0  # area / perimeter of a square of side ``fetch``
        h_n = natural_convection(surface, air, length_n)
        if wind <= 0.0:
            return h_n
        p = air_properties(0.5 * (surface + air))
        beta_g = GRAVITY * abs(dt) / air
        x_c = 0.0 if dt >= 0.0 else self.critical_reynolds * p.viscosity / wind
        # Forced part: the fetch means of 0.332 Re_x^(1/2) up to x_c and of
        # 0.0296 Re_x^(4/5) beyond, (k/L) [0.664 Re_c^(1/2) + 0.037 (Re_L^(4/5) - Re_c^(4/5))] Pr^(1/3)
        re_l = wind * self.fetch / p.viscosity
        re_c = min(wind * x_c / p.viscosity, re_l)
        nu_f = (0.664 * np.sqrt(re_c) + 0.037 * (re_l**0.8 - re_c**0.8)) * p.prandtl ** (1.0 / 3.0)
        h_f = self.roughness * p.conductivity * nu_f / self.fetch
        # Natural part, weighted by eta(x); x = fetch * s^5 clusters the nodes at the edge
        s = _NODES
        x = self.fetch * s**5
        jac = 5.0 * s**4  # dx / (fetch ds)
        eta = np.array([_eta(g) for g in beta_g * x / wind**2])
        return float(np.sum(_WEIGHTS * jac * eta) * h_n + h_f)


@dataclass(frozen=True)
class BulkConvection:
    """Bulk transfer of meteorology, h_c = rho cp C_H M (Stull 2015, eq. 3.35).

    M is the wind at 10 m, obtained from the wind at ``wind_height`` on the
    town profile. Natural convection is not added: C_H is a fit to the total
    turbulent exchange.

    Args:
        transfer_coefficient: C_H (2e-3 over smooth surfaces).
        wind_height: Height of the wind given to ``coefficient`` [m].
    """

    transfer_coefficient: float = 2e-3
    wind_height: float = 3.0
    includes_radiation: bool = False

    def coefficient(self, wind: float, surface: float, air: float) -> float:
        p = air_properties(air)
        m10 = wind_at_height(wind, self.wind_height, 10.0)
        return p.density * p.heat_capacity * self.transfer_coefficient * m10


@dataclass(frozen=True)
class McAdamsConvection:
    """McAdams, h = 5.7 + 3.8 v (Palyvos 2008, via the EnergyPlus Engineering Reference, eq. 3.104).

    Fitted by Juerges (1924) on a heated plate of 0.5 m in a wind tunnel; the
    coefficient contains the radiation to the tunnel walls. ``SurfaceBalance``
    therefore adds only the part of the long-wave loss beyond the exchange with
    surroundings at air temperature. Used as an upper bound only.
    """

    includes_radiation: bool = True

    def coefficient(self, wind: float, surface: float, air: float) -> float:
        return 5.7 + 3.8 * wind


# --- sky -------------------------------------------------------------------


def saturation_vapour_pressure(temperature: float) -> float:
    """Saturation vapour pressure over water [kPa] at ``temperature`` [K] (FAO-56, eq. 11)."""
    t = temperature - 273.15
    return 0.6108 * math.exp(17.27 * t / (t + 237.3))


def clear_sky_emissivity(model: str, air: float, vapour_pressure: float | None = None) -> float:
    """Clear-sky emissivity at screen-level air temperature [K] and vapour pressure [kPa].

    Formulas as listed by Flerchinger et al. (2009, Table 1).
    """
    if model == "swinbank":
        return 5.31e-13 * air**6 / (SIGMA * air**4)
    if model == "idso_jackson":
        return 1.0 - 0.261 * math.exp(-0.00077 * (air - 273.16) ** 2)
    if model == "prata":
        if vapour_pressure is None:
            raise ValueError("Prata needs the vapour pressure")
        # Precipitable water in cm. Flerchinger et al. print w = 4650 e/T with
        # e in kPa (Table 1, notes a and c), which gives w in mm: at 280 K and
        # 0.69 kPa it gives 11.5 and an emissivity of 0.97, against 0.73-0.75
        # from Swinbank and Idso & Jackson. 465 gives w in cm.
        w = 465.0 * vapour_pressure / air
        return 1.0 - (1.0 + w) * math.exp(-math.sqrt(1.2 + 3.0 * w))
    raise ValueError(f"unknown sky model {model!r}")


def sky_emissivity(clear: float, cloud: float) -> float:
    """Cloud correction of Unsworth & Monteith: (1 - 0.84 c) eps_clear + 0.84 c."""
    if not 0.0 <= cloud <= 1.0:
        raise ValueError(f"cloud fraction must be in [0, 1] (got {cloud})")
    return (1.0 - 0.84 * cloud) * clear + 0.84 * cloud


# --- surface balance -------------------------------------------------------


@dataclass(frozen=True)
class SurfaceTerms:
    """The terms of the surface balance at one time and surface temperature [W/m2 and W/m2/K]."""

    convection: float  # h_c (T_s - T_a), or for McAdams h (T_s - T_a)
    longwave: float  # eps sigma T_s^4 - eps L_down, or for McAdams eps (sigma T_a^4 - L_down)
    shortwave: float  # (1 - albedo) S_down, absorbed
    h_convection: float
    air: float
    sky: float  # downward long-wave L_down

    @property
    def net_loss(self) -> float:
        return self.convection + self.longwave - self.shortwave


@dataclass(frozen=True)
class SurfaceBalance:
    """Boundary condition of the free surface (``cooling.column.Boundary``).

    Args:
        air: Air temperature [K] against time since the start [s].
        vapour_pressure: Vapour pressure [kPa] against time [s].
        solar: Global horizontal irradiance under a clear sky [W/m2] against time [s].
        wind: Wind speed [m/s] at the height the convection model expects.
        convection: Convective coefficient.
        sky_model: "swinbank", "idso_jackson" or "prata".
        cloud: Cloud fraction (0-1), for the sky radiation and the sun.
        emissivity: Long-wave emissivity of the surface.
        albedo: Short-wave reflectance of the surface.
        cloud_transmission: Factor on the clear-sky irradiance for ``cloud``.
    """

    air: Callable[[float], float]
    vapour_pressure: Callable[[float], float]
    solar: Callable[[float], float]
    wind: float
    convection: Convection
    sky_model: str = "swinbank"
    cloud: float = 0.0
    emissivity: float = 0.95
    albedo: float = 0.07
    cloud_transmission: float = 1.0

    def terms(self, time: float, surface: float) -> SurfaceTerms:
        t_a = self.air(time)
        e = self.vapour_pressure(time) if self.sky_model == "prata" else None
        l_down = sky_emissivity(clear_sky_emissivity(self.sky_model, t_a, e), self.cloud) * SIGMA * t_a**4
        s_abs = (1.0 - self.albedo) * self.cloud_transmission * self.solar(time)
        h_c = self.convection.coefficient(self.wind, surface, t_a)
        if self.convection.includes_radiation:
            longwave = self.emissivity * (SIGMA * t_a**4 - l_down)
        else:
            longwave = self.emissivity * (SIGMA * surface**4 - l_down)
        return SurfaceTerms(h_c * (surface - t_a), longwave, s_abs, h_c, t_a, l_down)

    def coefficients(self, time: float, boundary_temperature: float) -> tuple[float, float]:
        """(h, T_ref) with q_out = h (T_s - T_ref), from the terms at the current T_s."""
        t_s = boundary_temperature
        terms = self.terms(time, t_s)
        h_c, t_a = terms.h_convection, terms.air
        if self.convection.includes_radiation:
            # q = h_c (T_s - T_a) + LW_deficit - S
            return h_c, t_a - (terms.longwave - terms.shortwave) / h_c
        # eps sigma (T_s^4 - T_e^4) = h_r (T_s - T_e) with the sky temperature T_e
        t_e = (terms.sky / SIGMA) ** 0.25
        h_r = self.emissivity * SIGMA * (t_s * t_s + t_e * t_e) * (t_s + t_e)
        h = h_c + h_r
        return h, (h_c * t_a + h_r * t_e + terms.shortwave) / h
