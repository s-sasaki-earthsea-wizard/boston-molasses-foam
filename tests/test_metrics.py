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


def _film_case(tmp_path: Path, alpha_body: float, first_alpha: float | None = None) -> Path:
    """A current 0.8-0.4 m deep over 100 m, half-parabolic, at t = 60 s.

    With ``alpha_body`` < 1 the liquid is spread as that fraction over the
    cells below 0.8/alpha_body of the depth, so that no cell is liquid
    (alpha >= 0.5) when ``alpha_body`` < 0.5. ``first_alpha`` replaces the
    alpha of the first cell under the current (air at the bottom).
    """
    x_faces = np.linspace(0.0, 120.0, 121)
    x = 0.5 * (x_faces[:-1] + x_faces[1:])
    depth = np.where(x < 100.0, 0.8 - 0.4 * x / 100.0, 0.0)
    y = 0.5 * (Y_FACES[:-1] + Y_FACES[1:])
    alpha = column_alpha(Y_FACES, depth / alpha_body) * alpha_body
    if first_alpha is not None:
        alpha[0, depth > 0] = first_alpha
    spread = np.maximum(depth / alpha_body, 1e-9)
    ux = np.where(y[:, None] < spread[None, :],
                  1e-2 * (2.0 * y[:, None] / spread[None, :] - (y[:, None] / spread[None, :]) ** 2), 0.0)
    alpha0 = column_alpha(Y_FACES, np.where(x < 10.0, COLUMN, 0.0))
    return write_case(tmp_path / "case", x_faces, Y_FACES,
                      {"0": (alpha0, np.zeros_like(alpha0)), "60": (alpha, ux)})


def test_film_check_wetted_body(tmp_path: Path) -> None:
    m = evaluate(_film_case(tmp_path, 1.0), front_at=(60.0,))
    film = m.film["60"]
    assert film.columns > 0
    assert film.bed_alpha == pytest.approx(1.0)
    assert film.bed_air == pytest.approx(0.0, abs=1e-12)
    assert film.surface_velocity == pytest.approx(1.0, abs=0.05)
    # Half-parabola sampled at dy1/2: 1 - dy1/(4h) of the wall value
    assert film.wall_shear_rate == pytest.approx(1.0, abs=0.03)
    assert film.wall_stress == pytest.approx(film.wall_shear_rate, rel=1e-9)


def test_film_check_stress_with_air_in_the_first_cell(tmp_path: Path) -> None:
    """Half air in the first cell halves its viscosity, so the stress against the shear rate."""
    film = evaluate(_film_case(tmp_path, 1.0, first_alpha=0.5), front_at=(60.0,)).film["60"]
    assert film.bed_alpha == pytest.approx(0.5)
    assert film.wall_stress == pytest.approx(0.5 * film.wall_shear_rate, rel=1e-3)


def test_film_check_without_liquid_cells(tmp_path: Path) -> None:
    """Diffuse liquid (alpha 0.4 everywhere) has no surface: no air, no surface speed."""
    m = evaluate(_film_case(tmp_path, 0.4), front_at=(60.0,))
    film = m.film["60"]
    assert film.columns > 0
    assert film.bed_alpha == pytest.approx(0.4)
    assert math.isnan(film.surface_velocity)
    assert math.isnan(film.bed_air)
