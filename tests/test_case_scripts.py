"""Regression tests of the case scripts: when make rebuilds the mesh and 0/.

The stamps of make/case.mk stand for a mesh or initial fields made from the
inputs older than them. These tests fake the generated files (no OpenFOAM)
and ask make whether they are up to date.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MESH_STAMP = Path("constant/polyMesh/.stamp")
INIT_STAMP = Path("0/.stamp")


def _copy_inputs(dst: Path) -> None:
    """Copy the pristine inputs of case/ (no mesh, no 0/) to ``dst``."""
    for name in ("0.orig", "system"):
        shutil.copytree(ROOT / "case" / name, dst / name)
    (dst / "constant").mkdir(parents=True)
    for f in (ROOT / "case" / "constant").iterdir():
        if f.name != "polyMesh":
            shutil.copy2(f, dst / "constant" / f.name)


def _set_times(case: Path, mtime: float) -> None:
    """Set the modification time of every input file of ``case``."""
    for name in ("0.orig", "system", "constant"):
        for f in (case / name).rglob("*"):
            if f.is_file():
                os.utime(f, (mtime, mtime))


def _set_step(case: Path, step: int, **env: str) -> None:
    subprocess.run(["bash", str(ROOT / "scripts/set_step.sh"), str(step), str(case)],
                   check=True, capture_output=True, env={**os.environ, **env})


def _fake_build(case: Path, mtime: float) -> None:
    """Stand in for blockMesh and setFields: write both stamps at ``mtime``."""
    for stamp in (MESH_STAMP, INIT_STAMP):
        path = case / stamp
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        os.utime(path, (mtime, mtime))


def _up_to_date(case: Path, stamp: Path) -> bool:
    result = subprocess.run(["make", "-q", f"CASE_DIR={case}", str(case / stamp)], cwd=ROOT,
                            capture_output=True, text=True)
    assert result.returncode in (0, 1), result.stderr
    return result.returncode == 0


def test_anim_case_rebuilds_after_another_step(tmp_path: Path) -> None:
    """Step 5 then Step 1 in the same ANIM_DIR must not reuse the Step 5 mesh."""
    now = time.time()
    src, dst = tmp_path / "src", tmp_path / "anim"
    _copy_inputs(src)
    _set_step(src, 1)
    _set_times(src, now - 1000)  # pristine inputs, older than any build

    subprocess.run(["bash", "scripts/prepare_anim_case.sh", str(src), str(dst)], cwd=ROOT,
                   check=True, capture_output=True)
    _set_step(dst, 5)
    _set_times(dst, now - 200)
    _fake_build(dst, now - 100)
    assert _up_to_date(dst, MESH_STAMP) and _up_to_date(dst, INIT_STAMP)

    # Step 1 inputs again: unchanged by set_step.sh, older than the Step 5 build
    subprocess.run(["bash", "scripts/prepare_anim_case.sh", str(src), str(dst)], cwd=ROOT,
                   check=True, capture_output=True)
    _set_step(dst, 1)
    assert "lengthX         180;" in (dst / "system/blockMeshDict").read_text()
    assert not _up_to_date(dst, MESH_STAMP)
    assert not _up_to_date(dst, INIT_STAMP)


def test_refinement_change_rebuilds_the_mesh(tmp_path: Path) -> None:
    """Another MESH_REFINE of the same step changes blockMeshDict, so the mesh."""
    now = time.time()
    case = tmp_path / "case"
    _copy_inputs(case)
    _set_step(case, 5)
    _set_times(case, now - 200)
    _fake_build(case, now - 100)
    assert _up_to_date(case, MESH_STAMP)

    _set_step(case, 5, MESH_REFINE="0.5")
    assert not _up_to_date(case, MESH_STAMP)
    assert not _up_to_date(case, INIT_STAMP)
