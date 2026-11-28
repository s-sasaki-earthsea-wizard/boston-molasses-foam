"""Tests of the column materials and the Andrade viscosity."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from cooling.properties import (
    MOLASSES,
    MOLASSES_B,
    MOLASSES_MU0,
    MOLASSES_T0,
    Material,
    andrade_viscosity,
)

SET_STEP = Path(__file__).resolve().parents[1] / "scripts" / "set_step.sh"


def _set_step_value(name: str) -> float:
    match = re.search(rf"^{name}=([0-9.eE+-]+)", SET_STEP.read_text(), flags=re.MULTILINE)
    assert match, f"{name} not found in {SET_STEP}"
    return float(match.group(1))


def test_molasses_matches_the_vof_case() -> None:
    assert MOLASSES.density == _set_step_value("MOLASSES_RHO")
    assert MOLASSES.heat_capacity == _set_step_value("MOLASSES_CP")
    assert MOLASSES.conductivity == _set_step_value("MOLASSES_KAPPA")
    assert MOLASSES_B == _set_step_value("MOLASSES_B")
    assert MOLASSES_T0 == _set_step_value("MOLASSES_T0")
    for variant, mu0 in MOLASSES_MU0.items():
        assert mu0 == _set_step_value(f"MOLASSES_MU0_{variant.upper()}")


def test_derived_properties() -> None:
    assert MOLASSES.diffusivity == pytest.approx(1.0575e-7, rel=1e-4)
    assert MOLASSES.effusivity == pytest.approx(953.29, rel=1e-4)


def test_andrade_reference_and_the_factor_three() -> None:
    mu0 = MOLASSES_MU0["central"]
    assert float(andrade_viscosity(MOLASSES_T0, mu0, MOLASSES_B, MOLASSES_T0)) == pytest.approx(mu0)
    # Sharp, Kennedy & Rubinstein (2016): x3 from 10 to 0 degC
    ratio = andrade_viscosity(273.15, mu0, MOLASSES_B, MOLASSES_T0) / mu0
    assert float(ratio) == pytest.approx(3.0, rel=0.01)
    temps = np.array([263.15, 273.15, 283.15, 298.15])
    assert np.all(np.diff(andrade_viscosity(temps, mu0, MOLASSES_B, MOLASSES_T0)) < 0.0)


def test_andrade_needs_kelvin() -> None:
    with pytest.raises(ValueError, match="kelvin"):
        andrade_viscosity(np.array([283.15, -5.0]), 55.0, MOLASSES_B, MOLASSES_T0)


@pytest.mark.parametrize("field", ["conductivity", "density", "heat_capacity"])
def test_material_needs_positive_properties(field: str) -> None:
    values = {"conductivity": 1.0, "density": 1.0, "heat_capacity": 1.0, field: 0.0}
    with pytest.raises(ValueError, match=field):
        Material("x", **values)
