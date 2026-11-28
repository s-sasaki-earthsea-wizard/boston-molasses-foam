"""Tests of post.metrics on small synthetic cases (tests/foamcase.py)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from post import metrics
from post.metrics import evaluate
from tests.foamcase import column_alpha, write_case

X_FACES = np.linspace(0.0, 40.0, 81)  # 0.5 m cells
Y_FACES = np.linspace(0.0, 4.0, 9)  # 0.5 m cells
COLUMN = 2.0  # [m] initial liquid depth over x < 10 m


def _x(x_faces: np.ndarray = X_FACES) -> np.ndarray:
    return 0.5 * (x_faces[:-1] + x_faces[1:])


def _front_case(path: Path, fronts: dict[str, float], x_faces: np.ndarray = X_FACES,
                y_faces: np.ndarray = Y_FACES, depth: float = 1.0, **kwargs) -> Path:
    """A layer of ``depth`` over x < front at each time ("0": the initial column)."""
    x = _x(x_faces)
    fields = {}
    for name, x_front in fronts.items():
        alpha = column_alpha(y_faces, np.where(x < x_front, COLUMN if name == "0" else depth, 0.0))
        fields[name] = (alpha, np.full_like(alpha, 0.1))
    return write_case(path, x_faces, y_faces, fields, **kwargs)


def test_front_speed_and_arrival(tmp_path: Path) -> None:
    case = _front_case(tmp_path / "case", {"0": 10.0, "3": 20.0, "8": 30.0, "20": 35.0}, step=4)
    m = evaluate(case)
    assert m.front_x == [9.75, 19.75, 29.75, 34.75]  # last cell centre with h >= 0.25 m
    assert m.front_speed == pytest.approx(2.0)
    assert m.arrival == {"50": None, "100": None, "150": None}
    assert m.volume0 == pytest.approx(10.0 * COLUMN)


def test_volume_from_unclipped_alpha(tmp_path: Path) -> None:
    """The volume keeps the solver's bounding errors; depth and front use clipped alpha."""
    alpha0 = column_alpha(Y_FACES, np.where(_x() < 10.0, COLUMN, 0.0))
    alpha = alpha0.copy()
    alpha[0, 0] = 1.0 + 1e-3  # overshoot in the column
    alpha[-1, -1] = -5e-4  # undershoot in the air
    case = write_case(tmp_path / "case", X_FACES, Y_FACES,
                      {"0": (alpha0, np.zeros_like(alpha0)), "10": (alpha, np.zeros_like(alpha))}, step=4)
    m = evaluate(case)
    assert m.volume_dev == pytest.approx((1e-3 - 5e-4) * 0.5 * 0.5)
    assert m.front_x == [9.75, 9.75]


def test_outlet_wetting_time(tmp_path: Path) -> None:
    """The first written time with liquid in the outlet column bounds the volume check."""
    case = _front_case(tmp_path / "case", {"0": 10.0, "5": 30.0, "10": 40.0, "15": 40.0}, step=4)
    m = evaluate(case)
    assert m.outlet_wetting_time == 10.0
    # 30 m² at 5 s against 20 m² at 0 s; the 40 m² at 10 and 15 s are left out
    assert m.volume_dev == pytest.approx(10.0)
    assert "| outlet column first holds alpha > 0.001 (written times) [s] | 10 |" in metrics.table([m])


def test_outlet_wetted_from_the_first_time(tmp_path: Path) -> None:
    """No written time precedes the wetting of the outlet column: no volume deviation."""
    case = _front_case(tmp_path / "case", {"0": 41.0, "5": 41.0}, step=4)
    m = evaluate(case)
    assert m.outlet_wetting_time == 0.0
    assert m.volume_dev is None
    assert "| largest deviation of V from V₀ before that [m²] | – |" in metrics.table([m])


def _strict_json(text: str):
    def reject(name: str):
        raise ValueError(f"invalid JSON constant {name}")

    return json.loads(text, parse_constant=reject)


def test_json_writes_null_for_undefined_values(tmp_path: Path) -> None:
    """No front at 5 s and no front speed (both NaN in the metrics) become null."""
    case = _front_case(tmp_path / "case", {"0": 10.0, "5": 0.0}, step=4, finished=False)
    out = tmp_path / "metrics.json"
    assert metrics.main([str(case), "--json", str(out)]) == 0
    (m,) = _strict_json(out.read_text())
    assert m["front_speed"] is None
    assert m["front_x"] == [9.75, None]
    assert m["outlet_wetting_time"] is None
