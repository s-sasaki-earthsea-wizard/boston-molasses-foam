"""Closed-form solutions of 1D conduction in a semi-infinite solid.

They are the references of the column tests. In all of them the solid
occupies depth d >= 0 below its surface, starts at a uniform temperature, and
its surface condition switches on at t = 0. These are the classical results
for a semi-infinite solid (e.g. Carslaw & Jaeger, Conduction of Heat in
Solids); the tests check each one against its differential equation and
boundary condition, so they do not rest on the citation alone.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.special import erf, erfc, erfcx


def _similarity(depth: np.ndarray | float, t: float, diffusivity: float) -> np.ndarray:
    if not t > 0.0:
        raise ValueError(f"t > 0 required (got {t})")
    return np.asarray(depth, dtype=float) / (2.0 * math.sqrt(diffusivity * t))


def surface_step(
    depth: np.ndarray | float, t: float, initial: float, surface: float, diffusivity: float
) -> np.ndarray:
    """Temperature after a step of the surface temperature.

        T(d, t) = T_s + (T_i - T_s) erf(d / (2 sqrt(a t)))

    Two semi-infinite solids brought into contact follow the same profile on
    each side, with T_s the contact temperature (see ``contact_temperature``).

    Args:
        depth: Depth below the surface d [m].
        t: Time since the step [s].
        initial: Initial temperature T_i [K].
        surface: Surface temperature T_s for t > 0 [K].
        diffusivity: Thermal diffusivity a [m2/s].

    Returns:
        Temperature [K], same shape as ``depth``.
    """
    return surface + (initial - surface) * erf(_similarity(depth, t, diffusivity))


def surface_step_flux(t: float, initial: float, surface: float, conductivity: float, diffusivity: float) -> float:
    """Heat flux out of the solid through its surface after a temperature step.

        q = k (T_i - T_s) / sqrt(pi a t)

    Positive when the solid loses heat (T_i > T_s).
    """
    if not t > 0.0:
        raise ValueError(f"t > 0 required (got {t})")
    return conductivity * (initial - surface) / math.sqrt(math.pi * diffusivity * t)


def contact_temperature(t1: float, e1: float, t2: float, e2: float) -> float:
    """Constant contact temperature of two semi-infinite solids.

        T_c = (e1 T1 + e2 T2) / (e1 + e2),  e = sqrt(k rho Cp)

    Args:
        t1, t2: Initial temperatures of the two solids [K].
        e1, e2: Their effusivities [J/m2/K/s^(1/2)].
    """
    return (e1 * t1 + e2 * t2) / (e1 + e2)


def surface_convection(
    depth: np.ndarray | float,
    t: float,
    initial: float,
    ambient: float,
    h: float,
    conductivity: float,
    diffusivity: float,
) -> np.ndarray:
    """Temperature with convective exchange at the surface, k dT/dd = h (T - T_inf) at d = 0.

        (T - T_i) / (T_inf - T_i) = erfc(eta) - exp(H d + H^2 a t) erfc(eta + H sqrt(a t))

    with eta = d / (2 sqrt(a t)) and H = h / k. The second term is evaluated
    as exp(-eta^2) erfcx(eta + H sqrt(a t)), which is the same expression
    without the overflow of exp(H d + H^2 a t).

    Args:
        depth: Depth below the surface d [m].
        t: Time since the exchange started [s].
        initial: Initial temperature T_i [K].
        ambient: Ambient temperature T_inf [K].
        h: Heat-transfer coefficient [W/m2/K].
        conductivity: Thermal conductivity k [W/m/K].
        diffusivity: Thermal diffusivity a [m2/s].

    Returns:
        Temperature [K], same shape as ``depth``.
    """
    eta = _similarity(depth, t, diffusivity)
    beta = h / conductivity * math.sqrt(diffusivity * t)
    fraction = erfc(eta) - np.exp(-(eta**2)) * erfcx(eta + beta)
    return initial + (ambient - initial) * fraction
