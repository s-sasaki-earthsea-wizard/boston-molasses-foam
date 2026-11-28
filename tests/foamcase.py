"""Write small synthetic OpenFOAM cases for the tests of post/.

A case has the inputs of case/ switched to a step by scripts/set_step.sh, the
mesh points and cell centres of a single-layer grid, and the fields
alpha.water, U and T of the given time directories, in ASCII.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
_HEADER = "FoamFile\n{{\n    format      ascii;\n    class       {cls};\n    object      {obj};\n}}\n\n"


def _scalars(name: str, values: np.ndarray) -> str:
    body = "\n".join(f"{v:.10g}" for v in values)
    return (_HEADER.format(cls="volScalarField", obj=name)
            + f"internalField   nonuniform List<scalar>\n{values.size}\n(\n{body}\n)\n;\n")


def _vectors(name: str, cls: str, values: np.ndarray, field: bool = True) -> str:
    body = "\n".join(f"({x:.10g} {y:.10g} {z:.10g})" for x, y, z in values)
    if not field:
        return _HEADER.format(cls=cls, obj=name) + f"{len(values)}\n(\n{body}\n)\n"
    return (_HEADER.format(cls=cls, obj=name)
            + f"internalField   nonuniform List<vector>\n{len(values)}\n(\n{body}\n)\n;\n")


def write_case(
    path: Path,
    x_faces: np.ndarray,
    y_faces: np.ndarray,
    fields: dict[str, tuple[np.ndarray, np.ndarray]],
    step: int = 5,
    variant: str = "central",
    finished: bool = True,
    temperature: float = 283.15,
) -> Path:
    """Write a case with the given grid and fields.

    Args:
        path: Case directory to create.
        x_faces: Cell-face coordinates along x [m].
        y_faces: Cell-face coordinates along y [m].
        fields: Time name -> (alpha, ux), each of shape (ny, nx); "0" is
            required (it sets the initial column).
        step: Step passed to scripts/set_step.sh.
        variant: Viscosity variant passed to scripts/set_step.sh.
        finished: Write a solver log that ends with "End".
        temperature: Uniform temperature of every time [K].
    """
    for name in ("system", "0.orig"):
        shutil.copytree(ROOT / "case" / name, path / name)
    (path / "constant" / "polyMesh").mkdir(parents=True)
    for name in ("g", "transportProperties"):
        shutil.copy(ROOT / "case" / "constant" / name, path / "constant" / name)
    subprocess.run(["bash", str(ROOT / "scripts/set_step.sh"), str(step), str(path), variant],
                   check=True, capture_output=True, env=dict(os.environ))

    xc = 0.5 * (x_faces[:-1] + x_faces[1:])
    yc = 0.5 * (y_faces[:-1] + y_faces[1:])
    px, py = np.meshgrid(x_faces, y_faces)
    points = np.concatenate([np.column_stack([px.ravel(), py.ravel(), np.full(px.size, z)])
                             for z in (0.0, 1.0)])
    (path / "constant" / "polyMesh" / "points").write_text(_vectors("points", "vectorField", points, False))
    cx, cy = np.meshgrid(xc, yc)  # row-major (y, x), as OpenFOAM numbers a blockMesh grid
    centres = np.column_stack([cx.ravel(), cy.ravel(), np.full(cx.size, 0.5)])

    for name, (alpha, ux) in fields.items():
        tdir = path / name
        tdir.mkdir()
        u = np.column_stack([np.asarray(ux, dtype=float).ravel(), np.zeros(alpha.size), np.zeros(alpha.size)])
        (tdir / "alpha.water").write_text(_scalars("alpha.water", np.asarray(alpha, dtype=float).ravel()))
        (tdir / "U").write_text(_vectors("U", "volVectorField", u))
        (tdir / "T").write_text(_HEADER.format(cls="volScalarField", obj="T")
                                + f"internalField   uniform {temperature};\n")
    (path / "0" / "C").write_text(_vectors("C", "volVectorField", centres))

    lines = [f"Time = {name}" for name in fields if name != "0"]
    if finished:
        lines.append("End")
    (path / f"log.interThermalFoam.step{step}").write_text("\n".join(lines) + "\n")
    return path


def column_alpha(y_faces: np.ndarray, depth: np.ndarray) -> np.ndarray:
    """Alpha of shape (ny, nx) holding liquid of ``depth`` [m] in each column."""
    lo, hi = y_faces[:-1, None], y_faces[1:, None]
    return np.clip((np.asarray(depth)[None, :] - lo) / (hi - lo), 0.0, 1.0)
