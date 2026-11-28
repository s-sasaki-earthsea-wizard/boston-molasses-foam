"""Thermal properties and the Andrade viscosity of the column materials."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Material:
    """Thermal properties of a homogeneous material.

    Attributes:
        name: Label used in outputs.
        conductivity: Thermal conductivity k [W/m/K].
        density: Density rho [kg/m3].
        heat_capacity: Specific heat capacity Cp [J/kg/K].
    """

    name: str
    conductivity: float
    density: float
    heat_capacity: float

    def __post_init__(self) -> None:
        for field in ("conductivity", "density", "heat_capacity"):
            value = getattr(self, field)
            if not value > 0.0:
                raise ValueError(f"{self.name}: {field} must be positive (got {value})")

    @property
    def volumetric_heat_capacity(self) -> float:
        """rho Cp [J/m3/K]."""
        return self.density * self.heat_capacity

    @property
    def diffusivity(self) -> float:
        """Thermal diffusivity a = k / (rho Cp) [m2/s]."""
        return self.conductivity / self.volumetric_heat_capacity

    @property
    def effusivity(self) -> float:
        """Thermal effusivity sqrt(k rho Cp) [J/m2/K/s^(1/2)]."""
        return math.sqrt(self.conductivity * self.volumetric_heat_capacity)


# Ground under the molasses (wiki Cooling-column-inputs, "Ground"). Only k and
# rho Cp enter the conduction.
# Granite: diffusivity 1.6e-6 m2/s and effusivity 2,177 J/m2/K/s^0.5 (USGS
# OFR 88-441, Table 13), so k = 2.754 W/m/K and rho Cp = 1.721 MJ/m3/K; the
# split of rho Cp into rho and Cp is nominal.
GRANITE = Material("granite", conductivity=2.754, density=2650.0, heat_capacity=649.5)
# Concrete with stone aggregate: k = 0.75 Btu/(ft h F) = 1.30 W/m/K and
# 140 lb/ft3 = 2,240 kg/m3 (TM 5-852-6, Table 2-2), Cp = 0.20 Btu/(lb F) =
# 837 J/kg/K (Table 2-1, "Concrete (avg. stone)").
CONCRETE = Material("concrete", conductivity=1.30, density=2240.0, heat_capacity=837.0)

# Cane final molasses of Steps 4-5 (scripts/set_step.sh, MOLASSES_*);
# tests/test_cooling_properties.py keeps the two in step.
MOLASSES = Material("molasses", conductivity=0.31, density=1430.0, heat_capacity=2050.0)
MOLASSES_B = 8.5e3  # [K]
MOLASSES_T0 = 283.15  # [K], reference temperature of mu0 and initial temperature
MOLASSES_MU0 = {"low": 20.0, "central": 55.0, "high": 200.0}  # [Pa s] at T0


def andrade_viscosity(temperature: np.ndarray | float, mu0: float, b: float, t0: float) -> np.ndarray:
    """Andrade (Arrhenius) viscosity mu = mu0 exp(B (1/T - 1/T0)).

    Unlike the VOF solver (TMin, muMax in constant/transportProperties), the
    temperature is not clipped and the viscosity is not capped: callers flag
    temperatures outside the range of the fit instead.

    Args:
        temperature: Temperature T [K].
        mu0: Viscosity at the reference temperature [Pa s].
        b: Activation temperature B [K].
        t0: Reference temperature T0 [K].

    Returns:
        Dynamic viscosity [Pa s], same shape as ``temperature``.
    """
    temperature = np.asarray(temperature, dtype=float)
    if np.any(temperature <= 0.0):
        raise ValueError("temperature must be positive (in kelvin)")
    return mu0 * np.exp(b * (1.0 / temperature - 1.0 / t0))
