"""Tests of post.metrics on small synthetic cases (tests/foamcase.py)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

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
