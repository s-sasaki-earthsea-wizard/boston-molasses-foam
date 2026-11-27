"""Reference-run metrics of dam-break cases: front, volume, bounds, cooling.

Usage (from the project root; ``make metrics`` wraps this):

    uv run python -m post.metrics CASE_DIR [CASE_DIR ...] [--json FILE] [--plot FILE]

Each case needs the cell centres in ``0/C`` (``postProcess -func
writeCellCentres -time 0``) and its solver log ``log.interThermalFoam.stepN``.

The front is the largest cell-centre x where the liquid depth
h(x) = sum(alpha dy) is at least ``FRONT_DEPTH`` = 0.25 m. Thinner liquid
ahead of it is left out on purpose: liquid fractions of a few per cent in the
bottom cells are dragged ahead of the front by the air, because a cell has
one velocity (see the README, "Reference runs").
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from post.foamio import read_entry, read_text
from post.frames import CaseInfo, Frame, load_case, load_frame

FRONT_DEPTH = 0.25  # [m] front: last column with h >= FRONT_DEPTH
SPEED_WINDOW = (3.0, 8.0)  # [s] front speed (x(t1) - x(t0)) / (t1 - t0)
ARRIVAL_X = (50.0, 100.0, 150.0)  # [m] arrival times of the front
OUTLET_ALPHA = 1e-3  # outflow starts when the last column holds more liquid
CORE_ALPHA = 0.99  # cells counted as liquid for the cooling
OUTLET_MARGIN = 1.0  # [m] plotted fronts end this close to the outlet
READ_THREADS = 16  # concurrent reads (network shares are latency-bound)

_ALPHA = re.compile(r"Min\(alpha\.water\) = (\S+)\s+Max\(alpha\.water\) = (\S+)")
_T = re.compile(r"T min/max = (\S+), (\S+)")
_CLOCK = re.compile(r"ExecutionTime = (\S+) s\s+ClockTime = (\S+) s")


@dataclass
class LogStats:
    """Summary of a solver log (bounds over every time step)."""

    path: str
    finished: bool
    time_steps: int
    execution_s: float | None
    clock_s: float | None
    alpha_min: float | None
    alpha_max: float | None
    t_min: float | None
    t_max: float | None


@dataclass
class CaseMetrics:
    """Metrics of one case; NaN or None where not defined."""

    case: str
    step: int | None
    variant: str | None
    rho: float
    mu0: float
    dx: float
    dy: float
    times: tuple[float, float]
    reynolds: float
    log: LogStats | None
    front_time: list[float] = field(repr=False)
    front_x: list[float] = field(repr=False)
    front_speed: float
    arrival: dict[str, float | None]
    outflow_start: float | None
    volume0: float
    volume_dev: float
    volume_end: float
    cooling: float | None
    mu_change: float | None
    wall_shear: float


def read_log(case: Path) -> LogStats | None:
    """Parse the newest ``log.interThermalFoam.step*`` of a case."""
    logs = sorted(case.glob("log.interThermalFoam.step*"), key=lambda p: p.stat().st_mtime)
    if not logs:
        return None
    a_min = a_max = t_min = t_max = clock = None
    steps, finished = 0, False
    for line in logs[-1].read_text().splitlines():
        if line.startswith("Time = "):
            steps += 1
        elif m := _ALPHA.search(line):
            lo, hi = float(m[1]), float(m[2])
            a_min = lo if a_min is None else min(a_min, lo)
            a_max = hi if a_max is None else max(a_max, hi)
        elif m := _T.search(line):
            lo, hi = float(m[1]), float(m[2])
            t_min = lo if t_min is None else min(t_min, lo)
            t_max = hi if t_max is None else max(t_max, hi)
        elif m := _CLOCK.search(line):
            clock = (float(m[1]), float(m[2]))
        elif line.strip() == "End":
            finished = True
    return LogStats(
        path=str(logs[-1]),
        finished=finished,
        time_steps=steps,
        execution_s=clock[0] if clock else None,
        clock_s=clock[1] if clock else None,
        alpha_min=a_min,
        alpha_max=a_max,
        t_min=t_min,
        t_max=t_max,
    )


def front_position(frame: Frame, x: np.ndarray) -> float:
    """Front x [m]: last cell centre with h >= FRONT_DEPTH (NaN if none)."""
    wet = frame.depth() >= FRONT_DEPTH
    return float(x[wet].max()) if wet.any() else float("nan")


def _at(times: np.ndarray, values: np.ndarray, t: float) -> float:
    hit = np.isclose(times, t)
    return float(values[hit][0]) if hit.any() else float("nan")


def evaluate(path: Path) -> CaseMetrics:
    """Compute the metrics of one case from its written time steps."""
    info: CaseInfo = load_case(path)
    grid = info.grid
    tp = read_text(path / "constant" / "transportProperties")
    try:
        variant = read_entry(tp, "viscosityVariant")
    except KeyError:
        variant = None
    b_visc = float(read_entry(tp, "B", block="liquidViscosity"))
    t0_visc = float(read_entry(tp, "T0", block="liquidViscosity"))
    mu_air = float(read_entry(tp, "rho", block="air")) * float(read_entry(tp, "nu", block="air"))

    names = list(info.times)
    need_t = info.solve_energy or info.temperature_dependent
    with ThreadPoolExecutor(max_workers=READ_THREADS) as pool:
        frames = list(pool.map(lambda n: load_frame(info, n, temperature=need_t), names))
    times = np.array([f.time for f in frames])
    area = grid.dx * grid.dy

    front = np.array([front_position(f, grid.x) for f in frames])
    volume = np.array([f.alpha.sum() * area for f in frames])
    outlet = np.array([f.alpha[:, -1].max() > OUTLET_ALPHA for f in frames])
    outflow_start = float(times[outlet][0]) if outlet.any() else None
    before = times < outflow_start if outflow_start is not None else np.ones_like(times, dtype=bool)

    arrival = {}
    for x_a in ARRIVAL_X:
        reached = front >= x_a
        arrival[f"{x_a:g}"] = float(times[reached][0]) if reached.any() else None

    def wall_shear(f: Frame) -> float:
        """Largest first-cell wall shear where h >= FRONT_DEPTH.

        Mixture viscosity as in updateViscosity.H, without its T clipping and
        muMax cap (neither is reached in these runs).
        """
        mu_liquid = np.full(grid.nx, info.mu0_liquid)
        if info.temperature_dependent:
            mu_liquid = info.mu0_liquid * np.exp(b_visc * (1.0 / f.temperature[0] - 1.0 / t0_visc))
        alpha = f.alpha[0]
        tau = (alpha * mu_liquid + (1.0 - alpha) * mu_air) * np.abs(f.ux[0]) / (grid.dy / 2.0)
        wet = f.depth() >= FRONT_DEPTH
        return float(tau[wet].max()) if wet.any() else 0.0

    t0, t1 = SPEED_WINDOW
    speed = (_at(times, front, t1) - _at(times, front, t0)) / (t1 - t0)

    cooling = mu_change = None
    if info.solve_energy:
        drops = [info.t_liquid - float(f.temperature[f.alpha > CORE_ALPHA].min())
                 for f in frames if (f.alpha > CORE_ALPHA).any()]
        cooling = max(drops) if drops else 0.0
        if info.temperature_dependent:
            t_i = info.t_liquid
            mu_change = float(np.exp(b_visc * (1.0 / (t_i - cooling) - 1.0 / t_i)) - 1.0)

    nu0 = info.mu0_liquid / info.rho_liquid
    return CaseMetrics(
        case=str(path),
        step=info.step,
        variant=variant,
        rho=info.rho_liquid,
        mu0=info.mu0_liquid,
        dx=grid.dx,
        dy=grid.dy,
        times=(float(times[0]), float(times[-1])),
        reynolds=float(np.sqrt(info.g * info.h0) * info.h0 / nu0),
        log=read_log(path),
        front_time=times.tolist(),
        front_x=front.tolist(),
        front_speed=speed,
        arrival=arrival,
        outflow_start=outflow_start,
        volume0=float(volume[0]),
        volume_dev=float(np.abs(volume[before] - volume[0]).max()),
        volume_end=float(volume[-1]),
        cooling=cooling,
        mu_change=mu_change,
        wall_shear=max(wall_shear(f) for f in frames),
    )


def case_label(m: CaseMetrics, with_mesh: bool) -> str:
    """Short label, e.g. "Step 4 central (55 Pa s)"."""
    if m.step is None:
        label = Path(m.case).name
    elif m.step == 4:
        label = f"Step 4 {m.variant or ''} ({m.mu0:g} Pa s)".replace("  ", " ")
    else:
        label = f"Step {m.step} (water)"
    return f"{label}, Δ = {m.dx:g} m" if with_mesh else label


def _fmt(value: float | None, spec: str) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "–"
    return format(value, spec)


def table(results: list[CaseMetrics]) -> str:
    """Markdown table, one column per case."""
    with_mesh = len({m.dx for m in results}) > 1
    rows: list[tuple[str, list[str]]] = [
        ("Δx = Δy [m]", [_fmt(m.dx, "g") for m in results]),
        ("written times [s]", [f"{m.times[0]:g}–{m.times[1]:g}" for m in results]),
        ("Re = √(gH) H / ν₀", [_fmt(m.reynolds, ".2g") for m in results]),
        ("finished / time steps", [f"{'yes' if m.log and m.log.finished else 'no'} / "
                                   f"{m.log.time_steps if m.log else '–'}" for m in results]),
        ("execution / clock time [s]", [f"{_fmt(m.log.execution_s, 'g')} / {_fmt(m.log.clock_s, 'g')}"
                                        if m.log else "–" for m in results]),
        ("alpha bounds (log)", [f"[{_fmt(m.log.alpha_min, '.1e')}, 1 + {_fmt(m.log.alpha_max - 1, '.1e')}]"
                                if m.log and m.log.alpha_max is not None else "–" for m in results]),
        ("T bounds (log) [K]", [f"[{_fmt(m.log.t_min, '.2f')}, {_fmt(m.log.t_max, '.2f')}]"
                                if m.log and m.log.t_min is not None else "–" for m in results]),
        (f"front speed {SPEED_WINDOW[0]:g}–{SPEED_WINDOW[1]:g} s [m/s]",
         [_fmt(m.front_speed, ".1f") for m in results]),
    ]
    for key in results[0].arrival:
        rows.append((f"front at x = {key} m [s]", [_fmt(m.arrival[key], "g") for m in results]))
    rows += [
        ("outflow starts [s]", [_fmt(m.outflow_start, "g") for m in results]),
        ("largest deviation of V from V₀ before outflow [m²]", [_fmt(m.volume_dev, ".1e") for m in results]),
        ("V₀, V(end) [m²]", [f"{m.volume0:.3f}, {m.volume_end:.1f}" for m in results]),
        ("max cooling, alpha > 0.99 [K]", [_fmt(m.cooling, ".1e") for m in results]),
        ("μ change for that cooling [%]",
         [_fmt(None if m.mu_change is None else m.mu_change * 100, ".2f") for m in results]),
        ("max wall shear, first cell [Pa]", [_fmt(m.wall_shear, ".2g") for m in results]),
    ]
    head = "| | " + " | ".join(case_label(m, with_mesh) for m in results) + " |"
    sep = "|---|" + "---|" * len(results)
    return "\n".join([head, sep] + [f"| {name} | " + " | ".join(vals) + " |" for name, vals in rows])


def plot_fronts(results: list[CaseMetrics], infos: list[CaseInfo], out: Path) -> None:
    """Front position over time against Ritter's front, saved as SVG or PNG."""
    import matplotlib as mpl

    mpl.use("Agg")
    from matplotlib.figure import Figure

    from post import plot as pp

    mpl.rcParams.update(pp._RC)
    mpl.rcParams["svg.hashsalt"] = "boston-molasses-fronts"

    # Ordinal one-hue ramp for the molasses viscosities (validated light->dark)
    ramp = {1: ["#2370cd"], 2: ["#6aabf1", "#0e3d7e"], 3: ["#6aabf1", "#2370cd", "#0e3d7e"]}
    mu_levels = sorted({m.mu0 for m in results if m.step == 4})
    colours = ramp.get(len(mu_levels)) or [
        mpl.colors.to_hex(c) for c in mpl.colors.LinearSegmentedColormap.from_list(
            "ramp", ["#6aabf1", "#0e3d7e"])(np.linspace(0, 1, len(mu_levels)))
    ]
    with_mesh = len({m.dx for m in results}) > 1
    finest = min(m.dx for m in results)

    fig = Figure(figsize=(8.0, 4.4), facecolor=pp.SURFACE)
    ax = fig.add_subplot()
    info = infos[0]
    t_end = max(m.times[1] for m in results)
    t_ref = np.linspace(0.0, t_end, 201)
    x_end = info.grid.extent[1]
    ax.plot(t_ref, info.x_dam + info.front_speed * t_ref, color=pp.REFERENCE, lw=1.0, ls=pp.DASH,
            label=f"Ritter front {info.x_dam:g} m + 2√(gH) t (frictionless)")
    for m, inf in zip(results, infos):
        style: dict = {"lw": 2.0}
        if m.step == 4:
            style["color"] = colours[mu_levels.index(m.mu0)]
        else:
            style.update(color=pp.INK_MUTED, lw=1.4, ls=(0, (1.5, 2)))
        if with_mesh and m.dx == finest:
            style.update(ls=(0, (4, 2)), lw=1.6)
        t = np.array([0.0, *m.front_time[1:]])
        x = np.array([inf.x_dam, *m.front_x[1:]])
        # Stop at the outlet: the last columns drain over the open boundary,
        # so the front definition no longer tracks a front there
        hit = np.nonzero(x >= x_end - OUTLET_MARGIN)[0]
        stop = hit[0] + 1 if hit.size else x.size
        ax.plot(t[:stop], x[:stop], label=case_label(m, with_mesh), **style)
    ax.axhline(x_end, color=pp.AXIS, lw=0.8)
    ax.text(t_end * 0.99, x_end - 3, f"outlet (x = {x_end:g} m)", ha="right", va="top",
            color=pp.INK_SECONDARY, fontsize=8)
    ax.set(xlim=(0, t_end), ylim=(0, x_end * 1.05), xlabel="time [s]",
           ylabel=f"front position, depth ≥ {FRONT_DEPTH:g} m [m]")
    ax.grid(True, color=pp.GRIDLINE, lw=0.6)
    ax.legend(loc="lower right", fontsize=7.5)
    runs = ", ".join(Path(m.case).name for m in results)
    fig.text(0.01, 0.01, f"OpenFOAM v2312, interThermalFoam; runs: {runs}", color=pp.INK_MUTED,
             fontsize=6.5, va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"Date": None} if out.suffix == ".svg" else None
    fig.savefig(out, dpi=150, metadata=metadata)


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="python -m post.metrics", description=__doc__.split("\n")[0])
    p.add_argument("cases", type=Path, nargs="+", help="case directories with results")
    p.add_argument("--json", type=Path, help="also write the metrics (with the front series) as JSON")
    p.add_argument("--plot", type=Path, help="plot the front positions to this file (.svg or .png)")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    results = [evaluate(c) for c in args.cases]
    print(f"Front: last column with depth >= {FRONT_DEPTH:g} m\n")
    print(table(results))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps([asdict(m) for m in results], indent=1))
        print(f"\nwrote {args.json}")
    if args.plot:
        plot_fronts(results, [load_case(c) for c in args.cases], args.plot)
        print(f"wrote {args.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
