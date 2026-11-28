"""Tests of the scenario runs of the cooling column: set-up, viscosity measures, summary and figures."""

from __future__ import annotations

import json

import numpy as np
import pytest

from cooling.analytic import contact_temperature
from cooling.column import Adiabatic
from cooling.figures import _spread, profile_figure
from cooling.properties import GRANITE, MOLASSES
from cooling.scenario import (
    THICKNESSES,
    Case,
    all_cases,
    build_column,
    doubling_hours,
    dump_summary,
    flow_weighted,
    run_case,
    sensitivity_cases,
    summarise,
    viscosity_ratio,
)
from cooling.surface import SurfaceBalance
from cooling.weather import DAYS


@pytest.fixture(scope="module")
def winter_5cm():
    return run_case(Case("winter", 0.05))


def test_case_rejects_unknown_inputs() -> None:
    with pytest.raises(ValueError, match="day"):
        Case("spring", 0.1)
    with pytest.raises(ValueError, match="path"):
        Case("winter", 0.1, path="air")


def test_columns_of_the_three_paths() -> None:
    both = build_column(Case("winter", 0.1))
    assert [layer.material.name for layer in both.layers] == ["concrete", "granite", "molasses"]
    assert sum(layer.thickness for layer in both.layers[:-1]) == pytest.approx(2.0)
    assert isinstance(both.surface, SurfaceBalance)
    assert isinstance(both.bottom, Adiabatic)
    assert both.layers[0].initial_temperature == DAYS["winter"].ground_initial
    assert both.layers[-1].initial_temperature == DAYS["winter"].molasses_initial

    surface = build_column(Case("winter", 0.1, path="surface"))
    assert [layer.material.name for layer in surface.layers] == ["molasses"]
    assert isinstance(build_column(Case("winter", 0.1, path="ground")).surface, Adiabatic)

    warm = build_column(Case("summer", 0.1, ground_initial=299.15, ground="granite"))
    assert [(layer.material.name, layer.initial_temperature) for layer in warm.layers[:-1]] == [("granite", 299.15)]


def test_viscosity_ratio() -> None:
    assert viscosity_ratio(np.array([284.35]), 284.35)[0] == pytest.approx(1.0)
    # B was fitted to a factor of 3 between 10 and 0 degC
    assert viscosity_ratio(np.array([273.15]), 283.15)[0] == pytest.approx(3.0, rel=1e-3)


def test_flow_weighted_viscosity() -> None:
    h, n = 0.1, 1000
    widths = np.full(n, h / n)
    height = (np.arange(n) + 0.5) * h / n
    assert flow_weighted(np.full((1, n), 3.0), height, widths, h)[0] == pytest.approx(3.0)
    # (H - z)^2 puts 7/8 of the weight on the bottom half
    stiff_bottom = np.where(height < h / 2, 2.0, 1.0)[None, :]
    assert flow_weighted(stiff_bottom, height, widths, h)[0] == pytest.approx(1.0 / (7 / 16 + 1 / 8), rel=1e-4)
    stiff_top = np.where(height < h / 2, 1.0, 2.0)[None, :]
    assert flow_weighted(stiff_top, height, widths, h)[0] == pytest.approx(1.0 / (7 / 8 + 1 / 16), rel=1e-4)


def test_doubling_hours() -> None:
    hours = np.array([0.0, 1.0, 2.0, 3.0])
    assert doubling_hours(hours, np.array([1.0, 1.5, 2.5, 3.0])) == pytest.approx(1.5)
    assert doubling_hours(hours, np.array([2.2, 2.0, 2.5, 3.0])) == 0.0
    assert doubling_hours(hours, np.array([1.0, 1.9, 1.5, 1.2])) is None


def test_run_starts_from_rest_and_closes_its_budget(winter_5cm) -> None:
    r = winter_5cm
    day = DAYS["winter"]
    assert r.budget_error < 1e-9
    assert r.hours[0] == 0.0
    assert r.hours[-1] == pytest.approx(12.0)
    assert r.ratio_flow[0] == pytest.approx(1.0)
    assert r.ratio_surface[0] == pytest.approx(1.0)
    # On contact the molasses at the ground takes the temperature of two semi-infinite bodies ...
    expected = contact_temperature(day.molasses_initial, MOLASSES.effusivity, day.ground_initial, GRANITE.effusivity)
    assert r.contact[0] == pytest.approx(expected)
    # ... and keeps it while both still look semi-infinite (10 min)
    assert r.contact[1] == pytest.approx(expected, abs=0.1)
    assert r.profile_ratio.shape == (4, r.profile_depth.size)
    assert set(r.surface_terms) == {"convection", "longwave", "shortwave"}


def test_paths() -> None:
    ground = run_case(Case("summer", 0.02, path="ground"))
    assert ground.budget_error < 1e-9
    assert not ground.surface_terms
    # The flow weighting looks at the bottom, which a 1 m layer keeps at its start for 12 h
    surface = run_case(Case("winter", 1.0, path="surface"))
    assert surface.ratio_flow[-1] == pytest.approx(1.0, abs=0.01)
    assert surface.ratio_surface[-1] > 2.0


def test_sensitivity_cases() -> None:
    cases = sensitivity_cases("winter", 0.1)
    assert cases[0] == Case("winter", 0.1)
    labels = [c.label for c in cases]
    assert len(set(labels)) == len(labels) == 18
    assert all(c.path == "both" and c.day == "winter" and c.thickness == 0.1 for c in cases)
    assert Case("winter", 0.1, ground_initial=272.65, label="ground cold") in cases
    assert len(all_cases()) == len(DAYS) * len(THICKNESSES) * (18 + 2)


def test_summary_lists_one_run_per_line_with_its_changes(winter_5cm) -> None:
    other = run_case(Case("winter", 0.05, wind=4.0, label="wind 4 m/s"))
    text = dump_summary(summarise([winter_5cm, other], "abc1234"))
    loaded = json.loads(text)
    assert loaded["revision"] == "abc1234"
    assert loaded["start_viscosity_Pa_s"]["winter"]["central"] == pytest.approx(48.46, abs=0.01)
    assert [run["case"] for run in loaded["runs"]] == [
        {"day": "winter", "thickness": 0.05, "path": "both", "label": "base"},
        {"day": "winter", "thickness": 0.05, "path": "both", "wind": 4.0, "label": "wind 4 m/s"},
    ]
    assert "series" in loaded["runs"][0]
    assert "series" not in loaded["runs"][1]
    assert loaded["runs"][0]["doubling_hours"]["contact"] == 0.0
    assert text.count('\n  {"case"') == 2


def test_spread_keeps_order_gap_and_bounds() -> None:
    values = np.random.default_rng(1).uniform(0.0, 1.0, 6)
    pos = _spread(values, 0.1, 0.0, 1.0)
    np.testing.assert_array_equal(np.argsort(pos), np.argsort(values))
    assert np.min(np.diff(np.sort(pos))) >= 0.1 - 1e-12
    assert -1e-12 <= pos.min() and pos.max() <= 1.0 + 1e-12
    far = np.array([0.1, 0.5, 0.9])
    np.testing.assert_allclose(_spread(far, 0.1, 0.0, 1.0), far)


def test_figures_are_byte_reproducible(tmp_path) -> None:
    results = [run_case(Case(day, 0.2)) for day in DAYS]
    first = profile_figure(results, tmp_path / "a.svg", "abc1234").read_bytes()
    second = profile_figure(results, tmp_path / "b.svg", "abc1234").read_bytes()
    assert first == second
    assert b"<dc:date>" not in first
    assert b"abc1234" in first  # the revision of the footer
