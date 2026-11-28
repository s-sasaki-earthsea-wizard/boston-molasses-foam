"""Tests of post.metrics on small synthetic cases (tests/foamcase.py)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from post import metrics
from post.analytic import huppert
from post.frames import load_case, load_frame
from post.metrics import evaluate
from tests.foamcase import column_alpha, write_case

X_FACES = np.linspace(0.0, 40.0, 81)  # 0.5 m cells
Y_FACES = np.concatenate([np.linspace(0.0, 2.0, 41), [2.5, 3.0, 4.0]])  # graded above 2 m
SQUARE_Y = np.linspace(0.0, 4.0, 9)  # 0.5 m cells, square with X_FACES
COLUMN = 2.0  # [m] initial liquid depth over x < 10 m


def _x(x_faces: np.ndarray = X_FACES) -> np.ndarray:
    return 0.5 * (x_faces[:-1] + x_faces[1:])


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


def _front_case(path: Path, length: float, fronts: dict[str, float], step: int, depth: float = 1.0,
                y_faces: np.ndarray = Y_FACES, **kwargs) -> Path:
    """A current of ``depth`` over x < front at each time, on a domain of ``length``.

    "0" is the initial column; ``kwargs`` go to ``write_case``.
    """
    x_faces = np.linspace(0.0, length, int(length / 0.5) + 1)
    x = _x(x_faces)
    fields = {}
    for name, x_front in fronts.items():
        alpha = column_alpha(y_faces, np.where(x < x_front, COLUMN if name == "0" else depth, 0.0))
        fields[name] = (alpha, np.full_like(alpha, 0.1))
    return write_case(path, x_faces, y_faces, fields, step=step, **kwargs)


@pytest.fixture
def figures(monkeypatch: pytest.MonkeyPatch) -> list:
    """Figures the plot functions would save."""
    saved: list = []
    monkeypatch.setattr(metrics, "_save", lambda fig, results, out: saved.append(fig))
    return saved


def _lines(ax, label: str):
    return next(line for line in ax.get_lines() if line.get_label() == label)


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_fronts_stop_at_their_own_outlet(tmp_path: Path, figures: list, order: tuple[int, int]) -> None:
    """A 180 m-style short domain must not cut the front of a longer domain."""
    short = _front_case(tmp_path / "short", 40.0, {"0": 10.0, "10": 30.0, "20": 40.0}, step=4)
    long = _front_case(tmp_path / "long", 120.0, {"0": 10.0, "10": 60.0, "20": 100.0}, step=5)
    cases = [(short, long)[i] for i in order]
    results = [evaluate(c) for c in cases]
    metrics.plot_fronts(results, [load_case(c) for c in cases], tmp_path / "fronts.svg")
    (ax,) = figures[0].axes
    by_step = {m.step: metrics.case_label(m, False) for m in results}
    long_line = _lines(ax, by_step[5])
    assert max(long_line.get_xdata()) == 20.0
    assert max(long_line.get_ydata()) == 99.75
    assert max(_lines(ax, by_step[4]).get_ydata()) == 39.75  # stops at its outlet
    assert ax.get_ylim()[1] >= 120.0
    outlets = sorted(line.get_ydata()[0] for line in ax.get_lines() if line.get_label().startswith("_"))
    assert outlets == [40.0, 120.0]


def test_huppert_front_starts_at_huppert_from(tmp_path: Path, figures: list) -> None:
    case = _front_case(tmp_path / "case", 120.0, {"0": 10.0, "300": 60.0, "600": 100.0}, step=5)
    metrics.plot_fronts([evaluate(case)], [load_case(case)], tmp_path / "fronts.svg", with_huppert=True)
    (ax,) = figures[0].axes
    huppert_lines = [line for line in ax.get_lines()
                     if line.get_label().startswith("_") and np.size(line.get_xdata()) > 2]
    assert len(huppert_lines) == 1
    t = huppert_lines[0].get_xdata()
    assert t[0] == pytest.approx(metrics.HUPPERT_FROM)
    assert t[-1] == pytest.approx(600.0)


def test_depth_axis_covers_huppert_profiles(tmp_path: Path, figures: list) -> None:
    """At 30 s Huppert's current is deeper at the wall than the simulated layer."""
    case = _front_case(tmp_path / "case", 120.0, {"0": 10.0, "30": 60.0}, step=5, depth=0.3)
    m = evaluate(case, profile_times=(30.0,))
    metrics.plot_profiles([m], [load_case(case)], tmp_path / "profiles.svg", with_huppert=True)
    (ax,) = figures[0].axes
    wall_depth = huppert(np.array([0.0]), 30.0, m.volume0, 9.81, m.mu0 / m.rho).depth[0]
    assert wall_depth > max(m.profiles["30"])  # the case this test is about
    assert ax.get_ylim()[1] >= wall_depth


def _legend_texts(fig) -> list[str]:
    return [t.get_text() for ax in fig.axes if ax.get_legend() for t in ax.get_legend().get_texts()]


@pytest.mark.parametrize(("t_end", "expected"), [(2.0, False), (5.0, False), (600.0, True)])
def test_huppert_front_legend_only_with_a_curve(tmp_path: Path, figures: list, t_end: float,
                                                expected: bool) -> None:
    """A run that ends by HUPPERT_FROM draws no Huppert front, so the legend has no entry for it."""
    case = _front_case(tmp_path / "case", 120.0, {"0": 10.0, f"{t_end:g}": 60.0}, step=5)
    metrics.plot_fronts([evaluate(case)], [load_case(case)], tmp_path / "fronts.svg", with_huppert=True)
    (ax,) = figures[0].axes
    assert any("Huppert" in t for t in _legend_texts(figures[0])) == expected
    assert any(np.size(line.get_xdata()) > 2 for line in ax.get_lines()
               if line.get_label().startswith("_")) == expected


@pytest.mark.parametrize(("step", "times", "expected"), [(1, (30.0,), False), (5, (0.0,), False),
                                                         (5, (30.0,), True)])
def test_huppert_legend_only_with_a_curve(tmp_path: Path, figures: list, step: int,
                                          times: tuple[float, ...], expected: bool) -> None:
    """Water, or only t = 0, draws no Huppert profile, so the legend has no entry for it."""
    case = _front_case(tmp_path / "case", 120.0, {"0": 10.0, "30": 60.0}, step=step, depth=0.3)
    m = evaluate(case, profile_times=times)
    metrics.plot_profiles([m], [load_case(case)], tmp_path / "profiles.svg", with_huppert=True)
    assert any("Huppert" in t for t in _legend_texts(figures[0])) == expected


def _two_cases(tmp_path: Path, finished: bool = True) -> tuple[Path, Path]:
    """Two Step 5 cases; the second stopped before 300 s."""
    full = _front_case(tmp_path / "full", 120.0, {"0": 10.0, "30": 40.0, "300": 60.0}, step=5, depth=0.3)
    short = write_case(tmp_path / "short", np.linspace(0.0, 120.0, 241), Y_FACES,
                       {name: fields for name, fields in _fields_of(full).items() if name != "300"},
                       step=5, variant="high", finished=finished)
    return full, short


def _fields_of(case: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    info = load_case(case)
    frames = {name: load_frame(info, name, temperature=False) for name in info.times}
    return {name: (f.alpha, f.ux) for name, f in frames.items()}


@pytest.mark.parametrize("order", [(0, 1), (1, 0)])
def test_missing_profile_time_fails_before_output(tmp_path: Path, capsys: pytest.CaptureFixture,
                                                  order: tuple[int, int]) -> None:
    pair = _two_cases(tmp_path)
    cases = [pair[i] for i in order]
    out = tmp_path / "profiles.svg"
    status = metrics.main([*map(str, cases), "--profiles", str(out), "--profile-times", "30", "300"])
    err = capsys.readouterr()
    assert status == 1
    assert not out.exists()
    assert "short: no output at t = 300 s" in err.err
    assert err.out == ""


def test_allow_missing_names_the_case_in_the_panel(tmp_path: Path, figures: list,
                                                   capsys: pytest.CaptureFixture) -> None:
    short, full = _two_cases(tmp_path)[::-1]
    status = metrics.main([str(short), str(full), "--profiles", str(tmp_path / "p.svg"),
                           "--profile-times", "30", "300", "--allow-missing"])
    assert status == 0
    assert "Warning:" in capsys.readouterr().err
    (fig,) = figures
    assert len(fig.axes) == 2  # both requested times, though the first case lacks 300 s
    texts = [t.get_text() for t in fig.axes[1].texts]
    assert any(t.startswith("no output: Step 5 high") for t in texts)
    assert sum("Step 5" in t for t in _legend_texts(fig)) == 2


def test_require_finished(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    full, short = _two_cases(tmp_path, finished=False)
    assert metrics.main([str(full), str(short), "--front-at", "30"]) == 0
    capsys.readouterr()
    assert metrics.main([str(full), str(short), "--front-at", "30", "--require-finished"]) == 1
    assert "short: log.interThermalFoam.step5 does not end with End" in capsys.readouterr().err


def test_huppert_values_only_for_viscous_cases(tmp_path: Path) -> None:
    water = _front_case(tmp_path / "water", 120.0, {"0": 10.0, "300": 60.0, "600": 80.0}, step=1, depth=0.3)
    m = evaluate(water, front_at=(300.0, 600.0))
    assert math.isnan(m.huppert_at["600"])
    assert math.isnan(m.huppert_exponent)
    assert math.isnan(m.huppert_wall_depth_end)
    assert "79.8, –" in metrics.table([m])  # front at 600 s; no Huppert front


def test_huppert_exponent_ignores_the_outlet(tmp_path: Path) -> None:
    """Huppert's front past the end of a short domain still has a growth exponent."""
    case = _front_case(tmp_path / "short", 40.0, {"0": 10.0, "300": 39.0, "600": 39.5}, step=5, depth=0.3)
    m = evaluate(case, front_at=(300.0, 600.0))
    assert m.huppert_at["300"] > 40.0  # beyond the outlet
    assert math.isnan(m.growth_exponent)  # the simulated front reaches the outlet
    assert 0.0 < m.huppert_exponent < 0.2


# Wave stage (PR #4): volume, outlet, input checks, JSON, table and front plot


def test_front_speed_and_arrival(tmp_path: Path) -> None:
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "3": 20.0, "8": 30.0, "20": 35.0}, step=4)
    m = evaluate(case)
    assert m.front_x == [9.75, 19.75, 29.75, 34.75]  # last cell centre with h >= 0.25 m
    assert m.front_speed == pytest.approx(2.0)
    assert m.arrival == {"50": None, "100": None, "150": None}
    assert m.volume0 == pytest.approx(10.0 * COLUMN)


def test_volume_from_unclipped_alpha(tmp_path: Path) -> None:
    """The volume keeps the solver's bounding errors; depth and front use clipped alpha."""
    alpha0 = column_alpha(SQUARE_Y, np.where(_x() < 10.0, COLUMN, 0.0))
    alpha = alpha0.copy()
    alpha[0, 0] = 1.0 + 1e-3  # overshoot in the column
    alpha[-1, -1] = -5e-4  # undershoot in the air
    case = write_case(tmp_path / "case", X_FACES, SQUARE_Y,
                      {"0": (alpha0, np.zeros_like(alpha0)), "10": (alpha, np.zeros_like(alpha))}, step=4)
    m = evaluate(case)
    assert m.volume_dev == pytest.approx((1e-3 - 5e-4) * 0.5 * 0.5)
    assert m.front_x == [9.75, 9.75]


def test_outlet_wetting_time(tmp_path: Path) -> None:
    """The first written time with liquid in the outlet column bounds the volume check."""
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "5": 30.0, "10": 40.0, "15": 40.0}, step=4)
    m = evaluate(case)
    assert m.outlet_wetting_time == 10.0
    # 30 m² at 5 s against 20 m² at 0 s; the 40 m² at 10 and 15 s are left out
    assert m.volume_dev == pytest.approx(10.0)
    assert "| outlet column first holds alpha > 0.001 (written times) [s] | 10 |" in metrics.table([m])


def test_outlet_wetted_from_the_first_time(tmp_path: Path) -> None:
    """No written time precedes the wetting of the outlet column: no volume deviation."""
    case = _front_case(tmp_path / "case", 40.0, {"0": 41.0, "5": 41.0}, step=4)
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
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "5": 0.0}, step=4, finished=False)
    out = tmp_path / "metrics.json"
    assert metrics.main([str(case), "--json", str(out)]) == 0
    (m,) = _strict_json(out.read_text())
    assert m["front_speed"] is None
    assert m["front_x"] == [9.75, None]
    assert m["outlet_wetting_time"] is None


def _run(tmp_path: Path, name: str, times: tuple[str, ...], finished: bool) -> Path:
    fronts = {t: 10.0 + float(t) for t in times}
    return _front_case(tmp_path / name, 40.0, fronts, step=4, finished=finished)


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
    case = _front_case(tmp_path / "case", 40.0, {t: 10.0 + 2.0 * float(t) for t in times}, step=4)
    m = evaluate(case)
    if speed is None:
        assert np.isnan(m.front_speed)
    else:
        assert m.front_speed == pytest.approx(speed, rel=1e-3)


def test_table_tells_meshes_apart_by_both_spacings(tmp_path: Path) -> None:
    """Same Δx, other Δy: the meshes differ, and neither is labelled Δx = Δy."""
    square = _front_case(tmp_path / "square", 40.0, {"0": 10.0, "20": 30.0}, y_faces=SQUARE_Y, step=4)
    flat = _front_case(tmp_path / "flat", 40.0, {"0": 10.0, "20": 30.0}, y_faces=np.linspace(0.0, 4.0, 41), step=4)
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
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "20": 30.0}, step=4)
    log = case / "log.interThermalFoam.step4"
    log.write_text(f"Time = 20\nMin(alpha.water) = -3.9e-05  Max(alpha.water) = {alpha_max}\nEnd\n")
    assert f"| alpha bounds (log) | {bounds} |" in metrics.table([evaluate(case)])


def test_plot_reads_each_case_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cases = [_front_case(tmp_path / name, 40.0, {"0": 10.0, "20": 30.0}, step=4) for name in ("a", "b")]
    read: list[Path] = []
    load = metrics.load_case
    monkeypatch.setattr(metrics, "load_case", lambda path: read.append(path) or load(path))
    assert metrics.main([*map(str, cases), "--plot", str(tmp_path / "fronts.svg")]) == 0
    assert read == cases


def test_front_plot_needs_one_initial_column(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    """Ritter's front is drawn once, so cases from other columns fail before any output."""
    a = _front_case(tmp_path / "a", 40.0, {"0": 10.0, "20": 30.0}, step=4)
    b = _front_case(tmp_path / "b", 40.0, {"0": 15.0, "20": 30.0}, step=4)
    out = tmp_path / "fronts.svg"
    assert metrics.main([str(a), str(b), "--plot", str(out)]) == 1
    err = capsys.readouterr()
    assert "b starts from another column than" in err.err
    assert "x_dam = 15 against 10 m" in err.err
    assert err.out == "" and not out.exists()
    assert metrics.main([str(a), str(b)]) == 0  # the table alone does not need one
    with pytest.raises(ValueError, match="another column"):
        metrics.plot_fronts([evaluate(a), evaluate(b)], [load_case(a), load_case(b)], out)


def test_line_styles_tell_liquid_and_mesh_apart(tmp_path: Path, figures: list) -> None:
    """Three meshes and water on two of them: no two plotted fronts look alike."""
    meshes = {"a": SQUARE_Y, "b": np.linspace(0.0, 4.0, 41), "c": np.linspace(0.0, 4.0, 17)}
    cases = [_front_case(tmp_path / f"{name}{step}", 40.0, {"0": 10.0, "20": 30.0}, y_faces=y, step=step)
             for name, y in meshes.items() for step in (4, 1) if (name, step) != ("c", 1)]
    results = [evaluate(c) for c in cases]
    metrics.plot_fronts(results, [load_case(c) for c in cases], tmp_path / "fronts.svg")
    (ax,) = figures[0].axes
    # Line2D keeps the dash pattern it was given only privately
    lines = [_lines(ax, metrics.case_label(m, True)) for m in results]
    looks = {(line.get_color(), str(line._unscaled_dash_pattern), line.get_linewidth()) for line in lines}
    assert len(looks) == len(results)
    assert lines[0]._unscaled_dash_pattern == (0, None)  # molasses on the first mesh: solid
    assert lines[1]._unscaled_dash_pattern == metrics._WATER_STYLE["ls"]  # water on it: dotted


def test_plotted_front_is_the_metrics_front(tmp_path: Path, figures: list) -> None:
    """The plot draws the front series as computed, t = 0 included (a cell centre)."""
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "10": 20.0, "20": 30.0}, step=4)
    m = evaluate(case)
    metrics.plot_fronts([m], [load_case(case)], tmp_path / "fronts.svg")
    line = _lines(figures[0].axes[0], metrics.case_label(m, False))
    assert list(line.get_xdata()) == m.front_time
    assert list(line.get_ydata()) == m.front_x == [9.75, 19.75, 29.75]


def test_outlet_label_above_the_outlet(tmp_path: Path, figures: list) -> None:
    """The fronts stay below the outlet line, so its label goes above it, inside the axes."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "10": 39.5, "20": 40.0}, step=4)
    metrics.plot_fronts([evaluate(case)], [load_case(case)], tmp_path / "fronts.svg")
    fig = figures[0]
    renderer = FigureCanvasAgg(fig).get_renderer()
    (ax,) = fig.axes
    (label,) = [t for t in ax.texts if t.get_text().startswith("outlet")]
    box = label.get_window_extent(renderer)
    assert box.y0 > ax.transData.transform((0.0, 40.0))[1]
    assert box.y1 <= ax.get_window_extent(renderer).y1


def test_front_plot_states_where_the_lines_end(tmp_path: Path, figures: list) -> None:
    """A front short of the last cell centre (39.75 m) ends at its first point within 1 m of the outlet."""
    case = _front_case(tmp_path / "case", 40.0, {"0": 10.0, "10": 39.0, "15": 39.5, "20": 39.5}, step=4)
    m = evaluate(case)
    metrics.plot_fronts([m], [load_case(case)], tmp_path / "fronts.svg")
    (ax,) = figures[0].axes
    line = _lines(ax, metrics.case_label(m, False))
    assert list(line.get_xdata()) == [0.0, 10.0, 15.0]
    assert list(line.get_ydata()) == [9.75, 38.75, 39.25]
    (label,) = [t.get_text() for t in ax.texts if t.get_text().startswith("outlet")]
    assert label == "outlet (x = 40 m); lines stop at the first front within 1 m"
