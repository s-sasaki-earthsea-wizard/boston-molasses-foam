"""Tests of the grid recovered from cell centres and mesh points."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from post.foamio import Grid2D, read_grid
from tests.foamcase import write_case


def _wall_graded_faces(first: float, ratio: float, top: float, largest: float) -> np.ndarray:
    """Faces from 0 to ``top``: cells growing by ``ratio`` from ``first`` up to ``largest``."""
    faces, dy = [0.0], first
    while faces[-1] < top - 1e-12:
        faces.append(min(faces[-1] + dy, top))
        dy = min(dy * ratio, largest)
    return np.array(faces)


def _centres(x_faces: np.ndarray, y_faces: np.ndarray) -> np.ndarray:
    cx, cy = np.meshgrid(0.5 * (x_faces[:-1] + x_faces[1:]), 0.5 * (y_faces[:-1] + y_faces[1:]))
    return np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, 0.5)])


# A 400 m long domain with 0.3 mm cells at the bottom (a resolved thermal
# boundary layer): the cell centres near the wall are closer than 1e-6 of
# the domain length
X_FACES = np.linspace(0.0, 400.0, 801)
Y_FACES = _wall_graded_faces(3e-4, 1.2, 20.0, 0.5)


def test_thin_cells_of_a_long_domain_from_centres() -> None:
    grid = Grid2D.from_cell_centres(_centres(X_FACES, Y_FACES), X_FACES, Y_FACES)
    assert (grid.nx, grid.ny) == (X_FACES.size - 1, Y_FACES.size - 1)
    np.testing.assert_allclose(grid.dy_cells[:3], [3e-4, 3.6e-4, 4.32e-4])


def test_thin_cells_of_a_long_domain_from_points(tmp_path: Path) -> None:
    x_faces = np.linspace(0.0, 400.0, 81)
    ny = Y_FACES.size - 1
    case = write_case(tmp_path / "case", x_faces, Y_FACES,
                      {"0": (np.zeros((ny, x_faces.size - 1)), np.zeros((ny, x_faces.size - 1)))})
    grid = read_grid(case)
    assert grid.ny == ny
    assert grid.dy_cells[0] == pytest.approx(3e-4)
