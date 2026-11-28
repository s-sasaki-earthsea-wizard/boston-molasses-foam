"""Tests of post.metrics on small synthetic cases (tests/foamcase.py)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from post.metrics import evaluate
from tests.foamcase import column_alpha, write_case

X_FACES = np.linspace(0.0, 40.0, 81)  # 0.5 m cells
Y_FACES = np.concatenate([np.linspace(0.0, 2.0, 41), [2.5, 3.0, 4.0]])  # graded above 2 m
COLUMN = 2.0  # [m] initial liquid depth over x < 10 m


def _x() -> np.ndarray:
    return 0.5 * (X_FACES[:-1] + X_FACES[1:])


def _initial() -> tuple[np.ndarray, np.ndarray]:
    alpha = column_alpha(Y_FACES, np.where(_x() < 10.0, COLUMN, 0.0))
    return alpha, np.zeros_like(alpha)


def test_evaluate_initial_fields_only(tmp_path: Path) -> None:
    """A case with only 0/ (just initialised) is evaluated without Huppert's values."""
    case = write_case(tmp_path / "case", X_FACES, Y_FACES, {"0": _initial()}, finished=False)
    m = evaluate(case)
    assert m.times == (0.0, 0.0)
    assert math.isnan(m.huppert_wall_depth_end)
    assert m.volume0 == pytest.approx(10.0 * COLUMN)
