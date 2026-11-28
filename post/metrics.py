"""Reference-run metrics of dam-break cases: front, volume, bounds, cooling.

Usage (from the project root; ``make metrics`` wraps this):

    uv run python -m post.metrics CASE_DIR [CASE_DIR ...] [--json FILE] [--plot FILE]
        [--arrival X ...] [--front-at T ...] [--huppert]
        [--profiles FILE --profile-times T ...]

Each case needs the cell centres in ``0/C`` (``postProcess -func
writeCellCentres -time 0``), its mesh points and its solver log
``log.interThermalFoam.stepN``. Uniform and graded meshes are supported.

The front is the largest cell-centre x where the liquid depth
h(x) = sum(alpha dy) is at least ``FRONT_DEPTH`` = 0.25 m. Thinner liquid
ahead of it is left out on purpose: liquid fractions of a few per cent in the
bottom cells are dragged ahead of the front by the air, because a cell has
one velocity (see the README, "Reference runs").

The spreading stage is compared with Huppert's (1982) similarity solution of a
2D viscous current of the initial area (``post.analytic.huppert``), with its
front taken at the same depth threshold as the simulated front.
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

from post.analytic import huppert, huppert_front
from post.foamio import read_entry, read_text
from post.frames import CaseInfo, Frame, load_case, load_frame

FRONT_DEPTH = 0.25  # [m] front: last column with h >= FRONT_DEPTH
SPEED_WINDOW = (3.0, 8.0)  # [s] front speed (x(t1) - x(t0)) / (t1 - t0)
ARRIVAL_X = (50.0, 100.0, 150.0)  # [m] arrival times of the front (default)
HUPPERT_FROM = 5.0  # [s] Huppert curves are drawn from this time on
OUTLET_ALPHA = 1e-3  # outflow starts when the last column holds more liquid
CORE_ALPHA = 0.99  # cells counted as liquid for the cooling
OUTLET_MARGIN = 1.0  # [m] plotted fronts end this close to the outlet
READ_THREADS = 16  # concurrent reads (network shares are latency-bound)
FILM_MARGIN = 20.0  # [m] film check: columns this far from the wall and the front
FILM_SLOPE_HALF = 10.0  # [m] film check: half window of the surface-slope fit
FILM_MIN_SPEED = 1e-4  # [m/s] film check: columns at rest are left out

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
class FilmCheck:
    """Body of the current against a viscous film of the same depth and flux.

    A lubrication flow has the half-parabolic profile
    u(y) = 3/2 ubar (2 y/h - y^2/h^2), so h and the mean velocity ubar fix the
    wall shear rate 3 ubar / h, the wall stress 3 mu0 ubar / h and the surface
    velocity 3/2 ubar; its flux is -g h^3 (dh/dx) / (3 nu). Each ratio (1 for
    such a film) is the median over the body columns: h >= FRONT_DEPTH, at
    least FILM_MARGIN from the wall and the front, and |ubar| > FILM_MIN_SPEED.
    NaN where there are no such columns.

    Shear rate and stress differ where the first cell holds air: its
    alpha-weighted viscosity is lower than mu0, so a weak bottom layer shows as
    a shear rate above 1 and a stress below 1 (an apparent slip). Both are
    first-cell estimates at the cell centre (dy1/2 from the wall).

    Attributes:
        columns: Number of body columns.
        depth: Median depth h [m].
        cells: Median number of cells over the depth.
        reynolds: Median ubar h / nu.
        wall_shear_rate: First-cell shear rate u1 / (dy1/2) against 3 ubar / h.
        wall_stress: First-cell stress mu u1 / (dy1/2), with the mixture
            viscosity mu of the solver (``MixtureViscosity``), against
            3 mu0 ubar / h.
        surface_velocity: Velocity of the top liquid cell (alpha >= 0.5)
            against 3/2 ubar; columns without such a cell (diffuse liquid)
            are left out.
        flux: Flux sum(alpha u dy) against -g h^3 (dh/dx) / (3 nu), with dh/dx
            fitted over +-FILM_SLOPE_HALF; negative where the flux runs up
            the fitted slope.
        bed_alpha: Median alpha of the first cell; 1 where the liquid wets the
            bottom.
        bed_air: Median air thickness sum((1 - alpha) dy) [m] below the top
            liquid cell, i.e. air trapped under or inside the current;
            columns without a liquid cell are left out.
    """

    columns: int
    depth: float
    cells: float
    reynolds: float
    wall_shear_rate: float
    wall_stress: float
    surface_velocity: float
    flux: float
    bed_alpha: float
    bed_air: float


    @classmethod
    def empty(cls) -> FilmCheck:
        """No body columns: every median NaN."""
        nan = float("nan")
        return cls(0, nan, nan, nan, nan, nan, nan, nan, nan, nan)


@dataclass(frozen=True)
class MixtureViscosity:
    """Cell viscosity of the solver, as in updateViscosity.H.

    mu = alpha mu_liquid(T) + (1 - alpha) mu_air with the Andrade
    mu_liquid(T) = mu0 exp(B (1/T - 1/T0)) when it is temperature dependent;
    without the T clipping and the muMax cap of the solver, which the
    reference runs do not reach. alpha is clipped to [0, 1] by load_frame.
    """

    mu0: float
    b: float
    t0: float
    mu_air: float
    temperature_dependent: bool

    @classmethod
    def of_case(cls, info: CaseInfo, transport_properties: str) -> MixtureViscosity:
        tp = transport_properties
        return cls(
            mu0=info.mu0_liquid,
            b=float(read_entry(tp, "B", block="liquidViscosity")),
            t0=float(read_entry(tp, "T0", block="liquidViscosity")),
            mu_air=float(read_entry(tp, "rho", block="air")) * float(read_entry(tp, "nu", block="air")),
            temperature_dependent=info.temperature_dependent,
        )

    def __call__(self, alpha: np.ndarray, temperature: np.ndarray | None) -> np.ndarray:
        """Viscosity [Pa s] of cells with liquid fraction ``alpha`` at ``temperature`` [K]."""
        if self.temperature_dependent:
            mu_liquid = self.mu0 * np.exp(self.b * (1.0 / temperature - 1.0 / self.t0))
        else:
            mu_liquid = np.full(np.shape(alpha), self.mu0)
        return alpha * mu_liquid + (1.0 - alpha) * self.mu_air


@dataclass
class CaseMetrics:
    """Metrics of one case; NaN or None where not defined."""

    case: str
    step: int | None
    variant: str | None
    rho: float
    mu0: float
    mesh: str
    cells: tuple[int, int]
    dy_wall: float
    times: tuple[float, float]
    reynolds: float
    log: LogStats | None
    front_time: list[float] = field(repr=False)
    front_x: list[float] = field(repr=False)
    front_speed: float
    arrival: dict[str, float | None]
    front_at: dict[str, float]
    huppert_at: dict[str, float]
    growth_exponent: float
    huppert_exponent: float
    wall_depth_end: float
    huppert_wall_depth_end: float
    film: dict[str, FilmCheck]
    profiles: dict[str, list[float]] = field(repr=False)
    profile_x: list[float] = field(repr=False)
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


def _index(times: np.ndarray, t: float) -> int | None:
    hit = np.flatnonzero(np.isclose(times, t))
    return int(hit[0]) if hit.size else None


def _at(times: np.ndarray, values: np.ndarray, t: float) -> float:
    i = _index(times, t)
    return float(values[i]) if i is not None else float("nan")


def _median(values: np.ndarray) -> float:
    """Median of the finite values; NaN if there are none."""
    finite = values[np.isfinite(values)]
    return float(np.median(finite)) if finite.size else float("nan")


def film_check(frame: Frame, info: CaseInfo, viscosity: MixtureViscosity) -> FilmCheck:
    """Compare the body of the current with a viscous film (see FilmCheck).

    Args:
        frame: Fields of one time step (with the temperature if the
            viscosity depends on it).
        info: Case information.
        viscosity: Mixture viscosity of the case, for the wall stress.
    """
    grid = info.grid
    nu = info.mu0_liquid / info.rho_liquid
    h = frame.depth()
    q = grid.column_sum(frame.alpha * frame.ux)
    wet = np.flatnonzero(h >= FRONT_DEPTH)
    if wet.size == 0:
        return FilmCheck.empty()
    x_front = grid.x[wet.max()]
    with np.errstate(invalid="ignore", divide="ignore"):
        moving = np.abs(np.where(h > 0, q / h, 0.0)) > FILM_MIN_SPEED
    body = ((grid.x >= info.x_wall + FILM_MARGIN) & (grid.x <= x_front - FILM_MARGIN)
            & (h >= FRONT_DEPTH) & moving)
    cols = np.flatnonzero(body)
    if cols.size == 0:
        return FilmCheck.empty()
    hb = h[cols]
    ubar = q[cols] / hb
    liquid = frame.alpha[:, cols] >= 0.5
    # Diffuse columns without a liquid cell have no surface: argmax would
    # give the top of the domain
    has_liquid = liquid.any(axis=0)
    top = grid.ny - 1 - np.argmax(liquid[::-1], axis=0)
    u_surface = np.where(has_liquid, frame.ux[top, cols], np.nan)
    below_top = np.arange(grid.ny)[:, None] < top[None, :]
    air = np.where(has_liquid,
                   ((1.0 - frame.alpha[:, cols]) * below_top * grid.dy_cells[:, None]).sum(axis=0), np.nan)
    slope = np.array([np.polyfit(grid.x[w], h[w], 1)[0] for w in
                      (np.abs(grid.x - grid.x[i]) <= FILM_SLOPE_HALF for i in cols)])
    rate = frame.ux[0, cols] / (grid.dy_cells[0] / 2.0)
    film_rate = 3.0 * ubar / hb
    temperature = frame.temperature[0, cols] if frame.temperature is not None else None
    mu_first = viscosity(frame.alpha[0, cols], temperature)
    with np.errstate(invalid="ignore", divide="ignore"):
        return FilmCheck(
            columns=int(cols.size),
            depth=float(np.median(hb)),
            cells=float(np.median(np.searchsorted(grid.y_faces, hb))),
            reynolds=float(np.median(np.abs(ubar) * hb / nu)),
            wall_shear_rate=float(np.median(rate / film_rate)),
            wall_stress=float(np.median(mu_first * rate / (info.mu0_liquid * film_rate))),
            surface_velocity=_median(u_surface / (1.5 * ubar)),
            flux=float(np.median(q[cols] / (-info.g * hb**3 * slope / (3.0 * nu)))),
            bed_alpha=float(np.median(frame.alpha[0, cols])),
            bed_air=_median(air),
        )


def growth_exponent(times: np.ndarray, front: np.ndarray, x_wall: float, x_end: float) -> float:
    """Slope of ln(x_f - x_wall) against ln t over the second half of the run.

    NaN if the front comes within ``OUTLET_MARGIN`` of the outlet by the end,
    since the outflow then sets the front. Huppert's nose x_N grows with 1/5;
    his front at a depth threshold grows more slowly as the current thins.
    """
    if not np.isfinite(front).all() or (front >= x_end - OUTLET_MARGIN).any():
        return float("nan")
    late = times >= times[-1] / 2.0
    if late.sum() < 2 or times[late][0] <= 0:
        return float("nan")
    return float(np.polyfit(np.log(times[late]), np.log(front[late] - x_wall), 1)[0])


def evaluate(
    path: Path,
    arrival_x: tuple[float, ...] = ARRIVAL_X,
    front_at: tuple[float, ...] = (),
    profile_times: tuple[float, ...] = (),
) -> CaseMetrics:
    """Compute the metrics of one case from its written time steps.

    Args:
        path: Case directory.
        arrival_x: Positions [m] whose front arrival times are reported.
        front_at: Times [s] at which the front and Huppert's front (at the same
            depth threshold) are reported.
        profile_times: Times [s] whose depth profiles are kept.
    """
    info: CaseInfo = load_case(path)
    grid = info.grid
    tp = read_text(path / "constant" / "transportProperties")
    try:
        variant = read_entry(tp, "viscosityVariant")
    except KeyError:
        variant = None
    viscosity = MixtureViscosity.of_case(info, tp)

    names = list(info.times)
    need_t = info.solve_energy or info.temperature_dependent
    with ThreadPoolExecutor(max_workers=READ_THREADS) as pool:
        frames = list(pool.map(lambda n: load_frame(info, n, temperature=need_t), names))
    times = np.array([f.time for f in frames])

    front = np.array([front_position(f, grid.x) for f in frames])
    volume = np.array([grid.integrate(f.alpha) for f in frames])
    outlet = np.array([f.alpha[:, -1].max() > OUTLET_ALPHA for f in frames])
    outflow_start = float(times[outlet][0]) if outlet.any() else None
    before = times < outflow_start if outflow_start is not None else np.ones_like(times, dtype=bool)

    arrival = {}
    for x_a in arrival_x:
        reached = front >= x_a
        arrival[f"{x_a:g}"] = float(times[reached][0]) if reached.any() else None

    def wall_shear(f: Frame) -> float:
        """Largest first-cell wall shear stress where h >= FRONT_DEPTH."""
        temperature = f.temperature[0] if f.temperature is not None else None
        tau = viscosity(f.alpha[0], temperature) * np.abs(f.ux[0]) / (dy_wall / 2.0)
        wet = f.depth() >= FRONT_DEPTH
        return float(tau[wet].max()) if wet.any() else 0.0

    dy_wall = float(grid.dy_cells[0])
    t0, t1 = SPEED_WINDOW
    speed = (_at(times, front, t1) - _at(times, front, t0)) / (t1 - t0)

    cooling = mu_change = None
    if info.solve_energy:
        drops = [info.t_liquid - float(f.temperature[f.alpha > CORE_ALPHA].min())
                 for f in frames if (f.alpha > CORE_ALPHA).any()]
        cooling = max(drops) if drops else 0.0
        if info.temperature_dependent:
            t_i = info.t_liquid
            mu_change = float(np.exp(viscosity.b * (1.0 / (t_i - cooling) - 1.0 / t_i)) - 1.0)

    nu0 = info.mu0_liquid / info.rho_liquid
    area0 = float(volume[0])
    x_end = grid.extent[1]
    wall_depth = frames[-1].depth()[0]
    t_last = float(times[-1])
    huppert_wall = (float(huppert(np.array([info.x_wall]), t_last, area0, info.g, nu0, info.x_wall).depth[0])
                    if t_last > 0 else float("nan"))
    profiles = {f"{t:g}": frames[i].depth().tolist()
                for t in profile_times if (i := _index(times, t)) is not None}
    return CaseMetrics(
        case=str(path),
        step=info.step,
        variant=variant,
        rho=info.rho_liquid,
        mu0=info.mu0_liquid,
        mesh=grid.spacing_label(),
        cells=(grid.nx, grid.ny),
        dy_wall=dy_wall,
        times=(float(times[0]), float(times[-1])),
        reynolds=float(np.sqrt(info.g * info.h0) * info.h0 / nu0),
        log=read_log(path),
        front_time=times.tolist(),
        front_x=front.tolist(),
        front_speed=speed,
        arrival=arrival,
        front_at={f"{t:g}": _at(times, front, t) for t in front_at},
        huppert_at={f"{t:g}": info.x_wall + float(huppert_front(t, area0, info.g, nu0, FRONT_DEPTH))
                    for t in front_at},
        growth_exponent=growth_exponent(times, front, info.x_wall, x_end),
        huppert_exponent=growth_exponent(
            times, info.x_wall + huppert_front(times, area0, info.g, nu0, FRONT_DEPTH), info.x_wall, x_end),
        wall_depth_end=float(wall_depth),
        huppert_wall_depth_end=huppert_wall,
        film={f"{t:g}": film_check(frames[i], info, viscosity)
              for t in front_at if (i := _index(times, t)) is not None},
        profiles=profiles,
        profile_x=grid.x.tolist(),
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
    elif m.step >= 4:
        label = f"Step {m.step} {m.variant or ''} ({m.mu0:g} Pa s)".replace("  ", " ")
    else:
        label = f"Step {m.step} (water)"
    return f"{label}, {m.mesh}" if with_mesh else label


def _film(c: FilmCheck | None) -> str:
    if c is None or c.columns == 0:
        return "–"
    return (f"{c.reynolds:.2g}; {c.flux:.2f}; {c.wall_shear_rate:.2f}; {c.wall_stress:.2f}; "
            f"{c.surface_velocity:.2f}")


def _bed(c: FilmCheck | None) -> str:
    if c is None or c.columns == 0:
        return "–"
    return f"{c.bed_alpha:.2f}, {100.0 * c.bed_air:.1f}"


def _fmt(value: float | None, spec: str) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "–"
    return format(value, spec)


def table(results: list[CaseMetrics]) -> str:
    """Markdown table, one column per case.

    The spreading-stage rows (front against Huppert, growth exponent, depth at
    the wall) appear when the front was asked for at given times.
    """
    with_mesh = len({m.mesh for m in results}) > 1
    rows: list[tuple[str, list[str]]] = [
        ("mesh (cells)", [f"{m.mesh} ({m.cells[0]} × {m.cells[1]})" for m in results]),
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
        rows.append((f"front at x = {key} m [s]", [_fmt(m.arrival.get(key), "g") for m in results]))
    for key in results[0].front_at:
        rows.append((f"front at t = {key} s, Huppert at h = {FRONT_DEPTH:g} m [m]",
                     [f"{_fmt(m.front_at.get(key), '.1f')}, {_fmt(m.huppert_at.get(key), '.1f')}"
                      for m in results]))
    if results[0].front_at:
        rows += [
            (f"growth exponent d ln x / d ln t, 2nd half; Huppert at h = {FRONT_DEPTH:g} m (x_N: 0.2)",
             [f"{_fmt(m.growth_exponent, '.3f')}; {_fmt(m.huppert_exponent, '.3f')}" for m in results]),
            ("depth at the wall at the end, Huppert [m]",
             [f"{_fmt(m.wall_depth_end, '.2f')}, {_fmt(m.huppert_wall_depth_end, '.2f')}" for m in results]),
        ]
        for key in results[0].front_at:
            rows.append((f"film at t = {key} s: Re; q/q_lub; γ̇₁/(3ū/h); τ₁/(3μ₀ū/h); u_s/(3ū/2) "
                         "(ratios 1 for a viscous film)",
                         [_film(m.film.get(key)) for m in results]))
        for key in results[0].front_at:
            rows.append((f"bed at t = {key} s: first-cell α, air below the surface [cm] (1, 0 when wetted)",
                         [_bed(m.film.get(key)) for m in results]))
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


# Ordinal one-hue ramp for the molasses viscosities (validated light->dark)
_RAMP = {1: ["#2370cd"], 2: ["#6aabf1", "#0e3d7e"], 3: ["#6aabf1", "#2370cd", "#0e3d7e"]}
_WATER_STYLE = {"lw": 1.4, "ls": (0, (1.5, 2))}
_FINER_MESH_STYLE = {"ls": (0, (4, 2)), "lw": 1.6}
_HUPPERT_LS = (0, (1, 1.5))


def _styles(results: list[CaseMetrics]) -> list[dict]:
    """Line style per case: viscous liquids on the ramp, water muted; other meshes dashed."""
    import matplotlib as mpl

    from post import plot as pp

    mu_levels = sorted({m.mu0 for m in results if m.step is not None and m.step >= 4})
    colours = _RAMP.get(len(mu_levels)) or [
        mpl.colors.to_hex(c) for c in mpl.colors.LinearSegmentedColormap.from_list(
            "ramp", ["#6aabf1", "#0e3d7e"])(np.linspace(0, 1, len(mu_levels)))
    ]
    with_mesh = len({m.mesh for m in results}) > 1
    styles = []
    for m in results:
        style: dict = {"lw": 2.0}
        if m.mu0 in mu_levels and m.step is not None and m.step >= 4:
            style["color"] = colours[mu_levels.index(m.mu0)]
        else:
            style.update(color=pp.INK_MUTED, **_WATER_STYLE)
        if with_mesh and m.mesh != results[0].mesh:
            style.update(_FINER_MESH_STYLE)
        styles.append(style)
    return styles


def _figure(salt: str, size: tuple[float, float], rows: int = 1):
    import matplotlib as mpl

    mpl.use("Agg")
    from matplotlib.figure import Figure

    from post import plot as pp

    mpl.rcParams.update(pp._RC)
    mpl.rcParams["svg.hashsalt"] = salt
    fig = Figure(figsize=size, facecolor=pp.SURFACE)
    axes = fig.subplots(rows, 1, sharex=True, squeeze=False)[:, 0]
    return fig, axes


def _save(fig, results: list[CaseMetrics], out: Path) -> None:
    from post import plot as pp

    runs = ", ".join(Path(m.case).name for m in results)
    fig.text(0.01, 0.01, f"OpenFOAM v2312, interThermalFoam; runs: {runs}", color=pp.INK_MUTED,
             fontsize=6.5, va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    metadata = {"Date": None} if out.suffix == ".svg" else None
    fig.savefig(out, dpi=150, metadata=metadata)


def _huppert_handle(lw: float, what: str = "viscous current"):
    from matplotlib.lines import Line2D

    from post import plot as pp

    return Line2D([], [], color=pp.INK_SECONDARY, lw=lw, ls=_HUPPERT_LS,
                  label=f"Huppert (1982) {what}, same colour (dotted)")


def plot_fronts(
    results: list[CaseMetrics], infos: list[CaseInfo], out: Path, with_huppert: bool = False
) -> None:
    """Front position over time against Ritter's front, saved as SVG or PNG.

    Each front stops at the outlet of its own domain; every outlet is drawn.

    Args:
        results: Metrics of the cases.
        infos: Case information, in the same order.
        out: Output file (.svg or .png).
        with_huppert: Also draw Huppert's viscous front of each viscous case,
            taken at the same depth threshold (``FRONT_DEPTH``),
            from ``HUPPERT_FROM`` on.
    """
    from post import plot as pp

    with_mesh = len({m.mesh for m in results}) > 1
    styles = _styles(results)
    fig, (ax,) = _figure("boston-molasses-fronts", (8.0, 4.4))
    info = infos[0]
    t_end = max(m.times[1] for m in results)
    t_ref = np.linspace(0.0, t_end, 201)
    # Denser early on, where Huppert's front rises fastest
    t_h = np.geomspace(HUPPERT_FROM, t_end, 201) if t_end > HUPPERT_FROM else np.empty(0)
    outlets = sorted({inf.grid.extent[1] for inf in infos})
    ax.plot(t_ref, info.x_dam + info.front_speed * t_ref, color=pp.REFERENCE, lw=1.0, ls=pp.DASH,
            label=f"Ritter front {info.x_dam:g} m + 2√(gH) t (frictionless)")
    drawn: set[float] = set()
    for m, inf, style in zip(results, infos, styles):
        t = np.array([0.0, *m.front_time[1:]])
        x = np.array([inf.x_dam, *m.front_x[1:]])
        # Stop at the outlet: the last columns drain over the open boundary,
        # so the front definition no longer tracks a front there
        hit = np.nonzero(x >= inf.grid.extent[1] - OUTLET_MARGIN)[0]
        stop = hit[0] + 1 if hit.size else x.size
        ax.plot(t[:stop], x[:stop], label=case_label(m, with_mesh), **style)
        if with_huppert and m.step is not None and m.step >= 4 and m.mu0 not in drawn:
            drawn.add(m.mu0)
            ax.plot(t_h, inf.x_wall + huppert_front(t_h, m.volume0, inf.g, m.mu0 / m.rho, FRONT_DEPTH),
                    color=style["color"], lw=1.2, ls=_HUPPERT_LS)
    for x_end in outlets:
        ax.axhline(x_end, color=pp.AXIS, lw=0.8)
        ax.text(t_end * 0.99, x_end - 3, f"outlet (x = {x_end:g} m)", ha="right", va="top",
                color=pp.INK_SECONDARY, fontsize=8)
    ax.set(xlim=(0, t_end), ylim=(0, outlets[-1] * 1.05), xlabel="time [s]",
           ylabel=f"front position, depth ≥ {FRONT_DEPTH:g} m [m]")
    ax.grid(True, color=pp.GRIDLINE, lw=0.6)
    handles, _ = ax.get_legend_handles_labels()
    if drawn:
        handles.append(_huppert_handle(1.2, f"viscous current, depth ≥ {FRONT_DEPTH:g} m"))
    ax.legend(handles=handles, loc="lower right", fontsize=7.5)
    _save(fig, results, out)


def plot_profiles(
    results: list[CaseMetrics], infos: list[CaseInfo], out: Path, with_huppert: bool = False
) -> None:
    """Depth profiles h(x) at the kept times, one panel per time.

    The x axis spans the longest domain.

    Args:
        results: Metrics of the cases, with the same profile times.
        infos: Case information, in the same order.
        out: Output file (.svg or .png).
        with_huppert: Also draw Huppert's profile of each viscous case.
    """
    from post import plot as pp

    keys = list(results[0].profiles)
    if not keys:
        raise ValueError("no depth profiles kept: give --profile-times that were written")
    with_mesh = len({m.mesh for m in results}) > 1
    styles = _styles(results)
    fig, axes = _figure("boston-molasses-profiles", (8.0, 1.0 + 1.6 * len(keys)), rows=len(keys))
    x_wall = min(inf.x_wall for inf in infos)
    x_end = max(inf.grid.extent[1] for inf in infos)
    h_max = max(max(m.profiles[k]) for m in results for k in keys if k in m.profiles)
    for ax, key in zip(axes, keys):
        t = float(key)
        drawn: set[float] = set()
        for m, inf, style in zip(results, infos, styles):
            if key not in m.profiles:
                continue
            ax.plot(m.profile_x, m.profiles[key], label=case_label(m, with_mesh), **style)
            if with_huppert and t > 0 and m.step is not None and m.step >= 4 and m.mu0 not in drawn:
                drawn.add(m.mu0)
                x = np.linspace(inf.x_wall, x_end, 2001)
                ax.plot(x, huppert(x, t, m.volume0, inf.g, m.mu0 / m.rho, inf.x_wall).depth,
                        color=style["color"], lw=1.2, ls=_HUPPERT_LS)
        ax.text(0.99, 0.92, f"t = {key} s", transform=ax.transAxes, ha="right", va="top",
                color=pp.INK_SECONDARY, fontsize=8)
        ax.set(xlim=(x_wall, x_end), ylim=(0, 1.1 * h_max), ylabel="depth h [m]")
        ax.grid(True, color=pp.GRIDLINE, lw=0.6)
    axes[-1].set_xlabel("x [m]")
    handles, _ = axes[0].get_legend_handles_labels()
    if with_huppert:
        handles.append(_huppert_handle(1.2))
    # Above the first panel: inside, it covered the deep early profiles
    axes[0].legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02), fontsize=7.5,
                   ncols=2, frameon=False)
    _save(fig, results, out)


def _floats(text: str) -> tuple[float, ...]:
    return tuple(float(v) for v in text.replace(",", " ").split())


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="python -m post.metrics", description=__doc__.split("\n")[0])
    p.add_argument("cases", type=Path, nargs="+", help="case directories with results")
    p.add_argument("--json", type=Path, help="also write the metrics (with the front series) as JSON")
    p.add_argument("--plot", type=Path, help="plot the front positions to this file (.svg or .png)")
    p.add_argument("--arrival", type=float, nargs="+", default=list(ARRIVAL_X), metavar="X",
                   help=f"positions [m] of the arrival times (default {' '.join(f'{x:g}' for x in ARRIVAL_X)})")
    p.add_argument("--front-at", type=float, nargs="+", default=[], metavar="T",
                   help="times [s] of the front against Huppert's at the same depth (adds the spreading rows)")
    p.add_argument("--huppert", action="store_true",
                   help="draw Huppert's viscous solution in the plots of viscous cases")
    p.add_argument("--profiles", type=Path, help="plot depth profiles at --profile-times to this file")
    p.add_argument("--profile-times", type=float, nargs="+", default=[], metavar="T",
                   help="times [s] of the depth profiles")
    args = p.parse_args(argv)
    if args.profiles and not args.profile_times:
        p.error("--profiles needs --profile-times")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    results = [evaluate(c, tuple(args.arrival), tuple(args.front_at), tuple(args.profile_times))
               for c in args.cases]
    print(f"Front: last column with depth >= {FRONT_DEPTH:g} m\n")
    print(table(results))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps([asdict(m) for m in results], indent=1))
        print(f"\nwrote {args.json}")
    if args.plot or args.profiles:
        infos = [load_case(c) for c in args.cases]
        if args.plot:
            plot_fronts(results, infos, args.plot, args.huppert)
            print(f"wrote {args.plot}")
        if args.profiles:
            plot_profiles(results, infos, args.profiles, args.huppert)
            print(f"wrote {args.profiles}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
