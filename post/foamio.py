"""Minimal readers for OpenFOAM ASCII case data.

Only what the animations and the metrics need is supported: the
internalField of volScalarField / volVectorField files (uniform or nonuniform,
plain or gzip-compressed), flat dictionary entries, time directory names, the
mesh points and the single-layer 2D grid (uniform or graded) recovered from the
cell centres written by ``postProcess -func writeCellCentres`` and the points.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_TIME_NAME = re.compile(r"^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$")
_UNIFORM = re.compile(r"internalField\s+uniform\s+(\([^)]*\)|[^;\s]+)\s*;")
_NONUNIFORM = re.compile(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*(\d+)\s*\(")
_N_COMPONENTS = {"scalar": 1, "vector": 3}


def read_text(path: Path) -> str:
    """Read a plain or gzip-compressed text file."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            return f.read()
    return path.read_text()


def field_path(time_dir: Path, name: str) -> Path:
    """Return the file of a field in a time directory (``name.gz`` or ``name``).

    Raises:
        FileNotFoundError: If neither file exists.
    """
    for candidate in (time_dir / f"{name}.gz", time_dir / name):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"field '{name}' not found in {time_dir}")


def read_internal_field(path: Path, n_cells: int | None = None) -> np.ndarray:
    """Read the internalField of a volScalarField or volVectorField.

    Args:
        path: Field file, plain or gzip-compressed.
        n_cells: Number of cells; required to expand a uniform field.

    Returns:
        Array of shape (n_cells,) for scalars or (n_cells, 3) for vectors.

    Raises:
        ValueError: If no supported internalField is found, or a uniform
            field is read without ``n_cells``.
    """
    text = read_text(path)

    m = _NONUNIFORM.search(text)
    if m:
        ncomp = _N_COMPONENTS[m.group(1)]
        n = int(m.group(2))
        body = text[m.end():]
        if ncomp == 1:
            body = body[: body.index(")")]
        else:
            body = body.replace("(", " ").replace(")", " ")
        values = np.array(body.split(None, n * ncomp)[: n * ncomp], dtype=float)
        return values.reshape(n, ncomp) if ncomp > 1 else values

    m = _UNIFORM.search(text)
    if m:
        if n_cells is None:
            raise ValueError(f"{path}: a uniform internalField needs n_cells")
        value = np.array(m.group(1).strip("()").split(), dtype=float)
        if value.size == 1:
            return np.full(n_cells, value[0])
        return np.tile(value, (n_cells, 1))

    raise ValueError(f"{path}: no supported internalField found")


def strip_comments(text: str) -> str:
    """Remove C/C++ style comments from dictionary text."""
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def read_entry(text: str, key: str, block: str | None = None) -> str:
    """Return the raw value of a ``key value;`` dictionary entry.

    Args:
        text: Dictionary file contents.
        key: Entry keyword.
        block: Name of a flat (non-nested) sub-dictionary to search in. Without
            it the first matching entry anywhere in the file is returned.

    Raises:
        KeyError: If the sub-dictionary or the entry is not found.
    """
    text = strip_comments(text)
    if block is not None:
        m = re.search(rf"^\s*{re.escape(block)}\s*\{{([^{{}}]*)\}}", text, re.M)
        if not m:
            raise KeyError(f"sub-dictionary '{block}' not found")
        text = m.group(1)
    m = re.search(rf"^\s*{re.escape(key)}\s+([^;]+);", text, re.M)
    if not m:
        where = f"'{block}'" if block else "the dictionary"
        raise KeyError(f"entry '{key}' not found in {where}")
    return m.group(1).strip()


def time_dirs(case: Path) -> list[str]:
    """Return the time directory names of a case (including 0), sorted by time."""
    names = [p.name for p in case.iterdir() if p.is_dir() and _TIME_NAME.match(p.name)]
    return sorted(names, key=float)


def read_gravity(case: Path) -> float:
    """Return the magnitude of the gravitational acceleration in constant/g."""
    value = read_entry(read_text(case / "constant" / "g"), "value")
    return float(np.linalg.norm(np.array(value.strip("()").split(), dtype=float)))


@dataclass(frozen=True)
class Grid2D:
    """Single-layer (x-y) structured grid, uniform or graded.

    Attributes:
        x: Cell-centre coordinates along x, shape (nx,).
        y: Cell-centre coordinates along y, shape (ny,).
        order: Cell indices in row-major (y, x) order, shape (nx*ny,).
        x_faces: Cell-face coordinates along x, shape (nx + 1,).
        y_faces: Cell-face coordinates along y, shape (ny + 1,).
    """

    x: np.ndarray
    y: np.ndarray
    order: np.ndarray
    x_faces: np.ndarray
    y_faces: np.ndarray

    @property
    def nx(self) -> int:
        return self.x.size

    @property
    def ny(self) -> int:
        return self.y.size

    @property
    def n_cells(self) -> int:
        return self.order.size

    @property
    def dx_cells(self) -> np.ndarray:
        """Cell widths along x, shape (nx,)."""
        return np.diff(self.x_faces)

    @property
    def dy_cells(self) -> np.ndarray:
        """Cell heights along y, shape (ny,)."""
        return np.diff(self.y_faces)

    @property
    def uniform(self) -> bool:
        """Whether the cells have one width and one height."""
        return _is_uniform(self.dx_cells) and _is_uniform(self.dy_cells)

    @property
    def dx(self) -> float:
        """Cell width along x.

        Raises:
            ValueError: If the widths vary (graded along x).
        """
        return _uniform_spacing(self.dx_cells, "x")

    @property
    def dy(self) -> float:
        """Cell height along y.

        Raises:
            ValueError: If the heights vary (graded along y).
        """
        return _uniform_spacing(self.dy_cells, "y")

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """Domain bounds (x0, x1, y0, y1), as used by imshow."""
        return (
            float(self.x_faces[0]),
            float(self.x_faces[-1]),
            float(self.y_faces[0]),
            float(self.y_faces[-1]),
        )

    def cell_areas(self) -> np.ndarray:
        """Cell areas in the x-y plane, shape (ny, nx)."""
        return self.dy_cells[:, None] * self.dx_cells[None, :]

    def column_sum(self, values: np.ndarray) -> np.ndarray:
        """Sum of ``values`` times the cell height over each column, shape (nx,).

        Args:
            values: Per-cell values of shape (ny, nx).
        """
        dy = self.dy_cells
        if _is_uniform(dy):
            # Scalar product keeps uniform-grid results bit-identical to sum * dy
            return values.sum(axis=0) * float(dy[0])
        return (values * dy[:, None]).sum(axis=0)

    def integrate(self, values: np.ndarray) -> float:
        """Area integral of per-cell values of shape (ny, nx)."""
        if self.uniform:
            return float(values.sum() * (self.dx * self.dy))
        return float((values * self.cell_areas()).sum())

    def spacing_label(self) -> str:
        """Cell size for captions: "Δ = 0.5 m", or per axis with ranges if graded."""
        def size(widths: np.ndarray) -> str:
            lo, hi = float(widths.min()), float(widths.max())
            return f"{lo:.3g} m" if _is_uniform(widths) else f"{lo:.3g}–{hi:.3g} m"

        dx, dy = size(self.dx_cells), size(self.dy_cells)
        return f"Δ = {dx}" if dx == dy and self.uniform else f"Δx = {dx}, Δy = {dy}"

    def to_2d(self, values: np.ndarray) -> np.ndarray:
        """Reorder per-cell values into shape (ny, nx) or (ny, nx, ncomp)."""
        return values[self.order].reshape(self.ny, self.nx, *values.shape[1:])

    @classmethod
    def from_cell_centres(
        cls,
        centres: np.ndarray,
        x_faces: np.ndarray | None = None,
        y_faces: np.ndarray | None = None,
    ) -> Grid2D:
        """Recover the grid from cell centres of shape (n_cells, 3).

        Args:
            centres: Cell centres.
            x_faces: Face coordinates along x; derived from uniformly spaced
                centres when omitted.
            y_faces: Face coordinates along y; as ``x_faces``.

        Raises:
            ValueError: If the cells do not form a single-layer x-y grid, or
                the faces do not enclose the centres.
            NotImplementedError: If the spacing is not uniform and no faces
                are given (a graded mesh needs the mesh points).
        """
        axes = []
        for name, c, faces in (("x", centres[:, 0], x_faces), ("y", centres[:, 1], y_faces)):
            # Per axis: a long domain must not merge thin cells along the other
            tol = _merge_tolerance(c)
            unique = _unique_sorted(c, tol)
            if unique.size < 2:
                raise ValueError(f"the mesh has a single cell along {name}")
            if faces is None:
                spacing = np.diff(unique)
                if not np.allclose(spacing, spacing[0], rtol=1e-6, atol=tol):
                    raise NotImplementedError(
                        f"non-uniform cell spacing along {name}: the face coordinates are needed"
                    )
                faces = np.concatenate(([unique[0] - spacing[0] / 2], unique + spacing[0] / 2))
            elif faces.size != unique.size + 1 or not np.allclose(
                0.5 * (faces[:-1] + faces[1:]), unique, rtol=0, atol=10 * tol
            ):
                raise ValueError(f"the faces along {name} do not enclose the cell centres")
            axes.append((unique, np.asarray(faces, dtype=float)))
        (xs, xf), (ys, yf) = axes
        if xs.size * ys.size != len(centres):
            raise ValueError(
                f"{len(centres)} cells do not form a single-layer {xs.size} x {ys.size} grid"
            )

        ix = np.searchsorted(xf, centres[:, 0]) - 1
        iy = np.searchsorted(yf, centres[:, 1]) - 1
        linear = iy * xs.size + ix
        if np.unique(linear).size != linear.size:
            raise ValueError("cell centres do not map one-to-one onto the grid")
        return cls(xs, ys, np.argsort(linear), xf, yf)


def _merge_tolerance(values: np.ndarray) -> float:
    """Coordinates closer than this along one axis are the same [m]."""
    return 1e-6 * float(np.ptp(values))


def _unique_sorted(values: np.ndarray, tol: float) -> np.ndarray:
    """Sorted distinct values, merging those closer than ``tol``."""
    s = np.sort(values)
    return s[np.concatenate(([True], np.diff(s) > tol))]


def _is_uniform(widths: np.ndarray) -> bool:
    return bool(np.allclose(widths, widths[0], rtol=1e-6, atol=0))


def _uniform_spacing(widths: np.ndarray, name: str) -> float:
    if not _is_uniform(widths):
        raise ValueError(f"the cells are graded along {name}; use the per-cell sizes")
    return float(widths[0])


def read_points(case: Path) -> np.ndarray:
    """Read the mesh points ``constant/polyMesh/points``, shape (n_points, 3)."""
    path = None
    for candidate in ("points.gz", "points"):
        if (case / "constant" / "polyMesh" / candidate).is_file():
            path = case / "constant" / "polyMesh" / candidate
            break
    if path is None:
        raise FileNotFoundError(f"no mesh points in {case / 'constant' / 'polyMesh'}")
    text = strip_comments(read_text(path))
    body = text[text.index("}") + 1:]  # skip the FoamFile header
    m = re.search(r"(\d+)\s*\(", body)
    if not m:
        raise ValueError(f"{path}: no point list found")
    n = int(m.group(1))
    values = body[m.end():].replace("(", " ").replace(")", " ").split(None, 3 * n)[: 3 * n]
    return np.array(values, dtype=float).reshape(n, 3)


def read_grid(case: Path) -> Grid2D:
    """Read the cell centres from ``0/C`` and the faces from the mesh points.

    Raises:
        FileNotFoundError: If ``0/C`` is missing; it is written by
            ``postProcess -func writeCellCentres -time 0`` (``make anim`` and
            ``make metrics`` run it).
    """
    try:
        path = field_path(case / "0", "C")
    except FileNotFoundError as err:
        raise FileNotFoundError(
            f"{err}; run 'postProcess -func writeCellCentres -time 0' in the case first"
        ) from None
    centres = read_internal_field(path)
    try:
        points = read_points(case)
    except FileNotFoundError:
        return Grid2D.from_cell_centres(centres)
    x, y = points[:, 0], points[:, 1]
    return Grid2D.from_cell_centres(
        centres, _unique_sorted(x, _merge_tolerance(x)), _unique_sorted(y, _merge_tolerance(y))
    )
