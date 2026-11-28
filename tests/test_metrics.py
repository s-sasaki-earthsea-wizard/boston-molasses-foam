"""Tests of post.metrics on small synthetic cases (tests/foamcase.py)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from post import metrics
from post.frames import load_case
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


def _run(tmp_path: Path, name: str, times: tuple[str, ...], finished: bool) -> Path:
    fronts = {t: 10.0 + float(t) for t in times}
    return _front_case(tmp_path / name, fronts, step=4, finished=finished)


@pytest.mark.parametrize(("kind", "message"), [
    ("initialised", "initialised: no solver log (log.interThermalFoam.step*)"),
    ("interrupted", "interrupted: log.interThermalFoam.step4 does not end with End"),
    ("short", "short: last written time 10 s, endTime 20 s"),
])
def test_require_finished_fails_before_any_output(tmp_path: Path, capsys: pytest.CaptureFixture,
                                                  kind: str, message: str) -> None:
    """Reference figures need finished runs; anything else fails before output."""
    finished = _run(tmp_path, "finished", ("0", "10", "20"), finished=True)
    case = {
        "initialised": _run(tmp_path, "initialised", ("0",), finished=False),
        "interrupted": _run(tmp_path, "interrupted", ("0", "10"), finished=False),
        "short": _run(tmp_path, "short", ("0", "10"), finished=True),  # End, but endTime is 20 s
    }[kind]
    if kind == "initialised":
        next(case.glob("log.*")).unlink()  # make init writes no solver log
    out, data = tmp_path / "fronts.svg", tmp_path / "metrics.json"
    args = [str(finished), str(case), "--plot", str(out), "--json", str(data)]
    assert metrics.main(args) == 0  # the generic metrics take partial runs
    capsys.readouterr()
    out.unlink()
    data.unlink()
    assert metrics.main([*args, "--require-finished"]) == 1
    err = capsys.readouterr()
    assert message in err.err
    assert err.out == ""
    assert not out.exists() and not data.exists()


def test_require_finished_passes_finished_runs(tmp_path: Path) -> None:
    case = _run(tmp_path, "finished", ("0", "10", "20"), finished=True)
    out = tmp_path / "fronts.svg"
    assert metrics.main([str(case), "--plot", str(out), "--require-finished"]) == 0
    assert out.exists()


@pytest.mark.parametrize(("times", "speed"), [
    (("0", "2.5", "3.5", "7.5", "8.5", "20"), 2.0),  # 3 and 8 s fall between written times
    (("0", "3.0001", "7.9999", "20"), 2.0),  # written a little off 3 and 8 s
    (("0", "2.5", "3.5", "6"), None),  # ends before 8 s: no extrapolation
])
def test_front_speed_between_written_times(tmp_path: Path, times: tuple[str, ...], speed: float | None) -> None:
    """The front moves at 2 m/s; its speed over 3-8 s needs no output at exactly 3 and 8 s."""
    case = _front_case(tmp_path / "case", {t: 10.0 + 2.0 * float(t) for t in times}, step=4)
    m = evaluate(case)
    if speed is None:
        assert np.isnan(m.front_speed)
    else:
        assert m.front_speed == pytest.approx(speed, rel=1e-3)


def test_table_tells_meshes_apart_by_both_spacings(tmp_path: Path) -> None:
    """Same Δx, other Δy: the meshes differ, and neither is labelled Δx = Δy."""
    square = _front_case(tmp_path / "square", {"0": 10.0, "20": 30.0}, step=4)
    flat = _front_case(tmp_path / "flat", {"0": 10.0, "20": 30.0}, y_faces=np.linspace(0.0, 4.0, 41), step=4)
    text = metrics.table([evaluate(square), evaluate(flat)])
    head, _, mesh_row = text.splitlines()[:3]
    assert "Step 4 central (55 Pa s), Δ = 0.5 m |" in head
    assert "Step 4 central (55 Pa s), Δx = 0.5 m, Δy = 0.1 m |" in head
    assert mesh_row == "| mesh (cells) | Δ = 0.5 m (80 × 8) | Δx = 0.5 m, Δy = 0.1 m (80 × 40) |"


@pytest.mark.parametrize(("alpha_max", "bounds"), [
    ("1.000043", "[-3.9e-05, 1 + 4.3e-05]"),
    ("0.99998", "[-3.9e-05, 1 - 2.0e-05]"),  # no overshoot
])
def test_alpha_bounds_with_sign(tmp_path: Path, alpha_max: str, bounds: str) -> None:
    case = _front_case(tmp_path / "case", {"0": 10.0, "20": 30.0}, step=4)
    log = case / "log.interThermalFoam.step4"
    log.write_text(f"Time = 20\nMin(alpha.water) = -3.9e-05  Max(alpha.water) = {alpha_max}\nEnd\n")
    assert f"| alpha bounds (log) | {bounds} |" in metrics.table([evaluate(case)])


def test_plot_reads_each_case_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cases = [_front_case(tmp_path / name, {"0": 10.0, "20": 30.0}, step=4) for name in ("a", "b")]
    read: list[Path] = []
    load = metrics.load_case
    monkeypatch.setattr(metrics, "load_case", lambda path: read.append(path) or load(path))
    assert metrics.main([*map(str, cases), "--plot", str(tmp_path / "fronts.svg")]) == 0
    assert read == cases


@pytest.fixture
def figures(monkeypatch: pytest.MonkeyPatch) -> list:
    """Figures plot_fronts would save."""
    from matplotlib.figure import Figure

    saved: list = []
    monkeypatch.setattr(Figure, "savefig", lambda fig, *args, **kwargs: saved.append(fig))
    return saved


def _line(ax, label: str):
    return next(line for line in ax.get_lines() if line.get_label() == label)


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_fronts_stop_at_their_own_outlet(tmp_path: Path, figures: list, order: tuple[int, int]) -> None:
    """A short domain must not cut the front of a longer one, in either order."""
    short = _front_case(tmp_path / "short", {"0": 10.0, "10": 30.0, "20": 40.0}, variant="low", step=4)
    long = _front_case(tmp_path / "long", {"0": 10.0, "10": 60.0, "20": 100.0}, x_faces=np.linspace(0.0, 120.0, 241),
                       variant="high", step=4)
    cases = [(short, long)[i] for i in order]
    results = [evaluate(c) for c in cases]
    metrics.plot_fronts(results, [load_case(c) for c in cases], tmp_path / "fronts.svg")
    (ax,) = figures[0].axes
    label = {m.variant: metrics.case_label(m, False) for m in results}
    assert max(_line(ax, label["high"]).get_ydata()) == 99.75
    assert max(_line(ax, label["low"]).get_ydata()) == 39.75  # stops at its outlet
    assert ax.get_ylim()[1] >= 120.0
    outlets = sorted(line.get_ydata()[0] for line in ax.get_lines() if line.get_label().startswith("_"))
    assert outlets == [40.0, 120.0]


def test_front_plot_needs_one_initial_column(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """Ritter's front is drawn once, so cases from other columns fail before any output."""
    a = _front_case(tmp_path / "a", {"0": 10.0, "20": 30.0}, step=4)
    b = _front_case(tmp_path / "b", {"0": 15.0, "20": 30.0}, step=4)
    out = tmp_path / "fronts.svg"
    assert metrics.main([str(a), str(b), "--plot", str(out)]) == 1
    err = capsys.readouterr()
    assert "b starts from another column than" in err.err
    assert "x_dam = 15 against 10 m" in err.err
    assert err.out == "" and not out.exists()
    assert metrics.main([str(a), str(b)]) == 0  # the table alone does not need one
    with pytest.raises(ValueError, match="another column"):
        metrics.plot_fronts([evaluate(a), evaluate(b)], [load_case(a), load_case(b)], out)
