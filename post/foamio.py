"""Minimal readers for OpenFOAM ASCII case data.

Only what the animations need is supported: the internalField of
volScalarField / volVectorField files (uniform or nonuniform, plain or
gzip-compressed), flat dictionary entries, time directory names and the
uniform 2D grid recovered from the cell centres written by
``postProcess -func writeCellCentres``.
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
    """Uniform single-layer (x-y) structured grid.

    Attributes:
        x: Cell-centre coordinates along x, shape (nx,).
        y: Cell-centre coordinates along y, shape (ny,).
        order: Cell indices in row-major (y, x) order, shape (nx*ny,).
    """

    x: np.ndarray
    y: np.ndarray
    order: np.ndarray

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
    def dx(self) -> float:
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        return float(self.y[1] - self.y[0])

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """Cell-face bounds (x0, x1, y0, y1), as used by imshow."""
        return (
            float(self.x[0] - self.dx / 2),
            float(self.x[-1] + self.dx / 2),
            float(self.y[0] - self.dy / 2),
            float(self.y[-1] + self.dy / 2),
        )

    def to_2d(self, values: np.ndarray) -> np.ndarray:
        """Reorder per-cell values into shape (ny, nx) or (ny, nx, ncomp)."""
        return values[self.order].reshape(self.ny, self.nx, *values.shape[1:])

    @classmethod
    def from_cell_centres(cls, centres: np.ndarray) -> Grid2D:
        """Recover the grid from cell centres of shape (n_cells, 3).

        Raises:
            ValueError: If the cells do not form a single-layer x-y grid.
            NotImplementedError: If the spacing is not uniform (graded mesh).
        """
        tol = 1e-6 * float(np.ptp(centres[:, :2], axis=0).max())
        axes = []
        for name, c in (("x", centres[:, 0]), ("y", centres[:, 1])):
            s = np.sort(c)
            unique = s[np.concatenate(([True], np.diff(s) > tol))]
            if unique.size < 2:
                raise ValueError(f"the mesh has a single cell along {name}")
            spacing = np.diff(unique)
            if not np.allclose(spacing, spacing[0], rtol=1e-6, atol=tol):
                raise NotImplementedError(
                    f"non-uniform cell spacing along {name}: graded meshes are not supported yet"
                )
            axes.append(unique)
        xs, ys = axes
        if xs.size * ys.size != len(centres):
            raise ValueError(
                f"{len(centres)} cells do not form a single-layer {xs.size} x {ys.size} grid"
            )

        ix = np.rint((centres[:, 0] - xs[0]) / (xs[1] - xs[0])).astype(int)
        iy = np.rint((centres[:, 1] - ys[0]) / (ys[1] - ys[0])).astype(int)
        linear = iy * xs.size + ix
        if np.unique(linear).size != linear.size:
            raise ValueError("cell centres do not map one-to-one onto the grid")
        return cls(xs, ys, np.argsort(linear))


def read_grid(case: Path) -> Grid2D:
    """Read the cell centres from ``0/C`` and recover the 2D grid.

    Raises:
        FileNotFoundError: If ``0/C`` is missing; it is written by
            ``postProcess -func writeCellCentres -time 0`` (``make anim``
            runs it).
    """
    try:
        path = field_path(case / "0", "C")
    except FileNotFoundError as err:
        raise FileNotFoundError(
            f"{err}; run 'postProcess -func writeCellCentres -time 0' in the case first"
        ) from None
    return Grid2D.from_cell_centres(read_internal_field(path))
