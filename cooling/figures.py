"""SVG figures of the cooling-column runs (#9).

The figures are written byte for byte the same on every run of the same
code: a fixed ``svg.hashsalt`` per figure and no date in the metadata. The
chart chrome and the blue ramp are those of ``post.plot``; the three hues of
the surface terms are the first three categorical slots of the dataviz
reference palette (validated on the surface #fcfcfb: all pairs pass the CVD
and normal-vision checks; orange and aqua are below 3:1 contrast, so their
lines carry direct labels).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from cooling import sun
from cooling.scenario import THICKNESSES, CaseResult
from cooling.weather import COLLAPSE_HOUR, DAYS

PROFILE_THICKNESS = 0.2  # [m]
SENSITIVITY_THICKNESS = 0.1  # [m]
SENSITIVITY_HOUR = 6.0  # [h] after the collapse: 18:40, two hours after the winter sunset
TERMS_THICKNESS = 0.1  # [m]
TERM_COLOURS = {"convection": "#257de0", "longwave": "#f47230", "shortwave": "#00b894"}
BASE_NOTE = (
    "Base case: Clear et al. (2001) convection over a 30 m fetch, wind 2 m/s at 3 m, clear sky (Swinbank), "
    "emissivity 0.95, albedo 0.07, 0.15 m granite on concrete; Andrade viscosity, B = 8.5e3 K"
)


def _pp():
    import matplotlib as mpl

    mpl.use("Agg")
    from post import plot as pp

    return mpl, pp


def _figure(salt: str, size: tuple[float, float], rows: int, cols: int, **kwargs):
    mpl, pp = _pp()
    from matplotlib.figure import Figure

    mpl.rcParams.update(pp._RC)
    mpl.rcParams["svg.hashsalt"] = salt
    fig = Figure(figsize=size, facecolor=pp.SURFACE)
    axes = fig.subplots(rows, cols, squeeze=False, **kwargs)
    for ax in axes.flat:
        ax.grid(True, color=pp.GRIDLINE, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    return fig, axes


def _save(fig, out: Path, revision: str, title: str) -> Path:
    _, pp = _pp()
    fig.suptitle(title, x=0.01, ha="left", color=pp.INK, fontsize=10)
    source = f"1D cooling column, cooling/scenario.py at {revision or 'working tree'}. {BASE_NOTE}."
    fig.text(0.01, 0.01, source, color=pp.INK_MUTED, fontsize=6.5, va="bottom", wrap=True)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, metadata={"Date": None})
    return out


def _ramp(n: int, start: float = 0.3) -> list:
    _, pp = _pp()
    return [pp.SPEED_CMAP(x) for x in np.linspace(start, 1.0, n)]


def _clock_axis(ax, label: bool = True) -> None:
    ticks = np.arange(14, 25, 2)
    ax.set_xticks(ticks - COLLAPSE_HOUR, [f"{h % 24:02d}:00" for h in ticks])
    ax.set_xlim(0.0, 12.0)
    if label:
        ax.set_xlabel("Clock time, EST (collapse at 12:40)")


def _sunset(ax, day: str) -> None:
    _, pp = _pp()
    d = DAYS[day]
    x = sun.sunset(d.day, d.site) - COLLAPSE_HOUR
    ax.axvline(x, color=pp.INK_MUTED, lw=0.9, ls=(0, (1, 2)), zorder=1)
    ax.annotate("sunset", (x, 1.0), xycoords=("data", "axes fraction"), xytext=(3, -3),
                textcoords="offset points", ha="left", va="top", fontsize=7.5, color=pp.INK_MUTED)


def _spread(values: np.ndarray, gap: float, low: float, high: float) -> np.ndarray:
    """Positions near ``values``, at least ``gap`` apart and in the same order, inside [low, high].

    Labels that would overlap are merged into a cluster centred on the mean of
    their targets, until no two clusters overlap.
    """
    order = np.argsort(values)
    target = np.clip(np.asarray(values, dtype=float)[order], low, high)
    clusters = [[i, i + 1] for i in range(target.size)]  # [first, end) in sorted order

    def place(c: list[int]) -> np.ndarray:
        n = c[1] - c[0]
        return target[c[0]:c[1]].mean() + gap * (np.arange(n) - 0.5 * (n - 1))

    merged = True
    while merged:
        merged = False
        for k in range(len(clusters) - 1):
            if place(clusters[k + 1])[0] - place(clusters[k])[-1] < gap:
                clusters[k] = [clusters[k][0], clusters[k + 1][1]]
                del clusters[k + 1]
                merged = True
                break
    pos = np.concatenate([place(c) for c in clusters])
    pos += max(0.0, low - pos[0]) - max(0.0, pos[-1] - high)
    out = np.empty_like(pos)
    out[order] = pos
    return out


def _end_labels(ax, ends: np.ndarray, labels: list[str], log: bool, x: float = 12.0) -> None:
    """Direct labels right of the line ends, spread so that they do not overlap."""
    _, pp = _pp()
    lo, hi = ax.get_ylim()
    if log:
        ends, lo, hi = np.log10(ends), np.log10(lo), np.log10(hi)
    pos = _spread(ends, 0.05 * (hi - lo), lo, hi)
    for y, text in zip(pos, labels, strict=True):
        ax.annotate(text, (x, 10**y if log else y), xytext=(4, 0), textcoords="offset points", va="center",
                    ha="left", fontsize=7.5, color=pp.INK_SECONDARY, annotation_clip=False)


def _ratio_axis(ax) -> None:
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator([0.125, 0.25, 0.5, 1, 2, 4, 8]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"×{v:g}"))


def _by(results: list[CaseResult], **match) -> CaseResult:
    found = [r for r in results if all(getattr(r.case, k) == v for k, v in match.items())]
    if len(found) != 1:
        raise LookupError(f"{len(found)} runs match {match}")
    return found[0]


def _cm(h: float) -> str:
    return f"{h * 100:g} cm" if h < 1.0 else f"{h:g} m"


def viscosity_figure(results: list[CaseResult], out: Path, revision: str) -> Path:
    """Flow-weighted and surface viscosity against time, for every thickness, both days."""
    _, pp = _pp()
    fig, axes = _figure("cooling-viscosity", (pp.FIG_WIDTH, 6.6), 2, 2, sharex=True, sharey="row")
    colours = _ramp(len(THICKNESSES))
    rows = (("ratio_flow", "Flow-weighted viscosity\nrelative to the start"),
            ("ratio_surface", "Viscosity at the surface\nrelative to the start"))
    for i, (attr, ylabel) in enumerate(rows):
        series = {(day, h): getattr(_by(results, day=day, thickness=h, label="base"), attr)
                  for day in DAYS for h in THICKNESSES}
        values = np.concatenate(list(series.values()))
        for j, day in enumerate(DAYS):
            ax = axes[i, j]
            hours = _by(results, day=day, thickness=THICKNESSES[0], label="base").hours
            for h, c in zip(THICKNESSES, colours, strict=True):
                ax.plot(hours, series[(day, h)], color=c, lw=2.0, label=_cm(h), solid_capstyle="round")
            ax.axhline(2.0, color=pp.INK_SECONDARY, lw=0.9, ls=pp.DASH, zorder=1)
            _ratio_axis(ax)
            ax.set_ylim(0.8 * values.min(), 1.25 * values.max())
            _sunset(ax, day)
            _clock_axis(ax, label=i == 1)
            if i == 0:
                ax.set_title(DAYS[day].name if day == "winter" else f"{DAYS[day].name} (normals of 1991–2020)",
                             loc="left", fontsize=9, color=pp.INK)
            if j == 0:
                ax.set_ylabel(ylabel)
                ax.annotate("doubled", (0.0, 2.0), xytext=(3, 2), textcoords="offset points", fontsize=7.5,
                            color=pp.INK_SECONDARY, va="bottom")
            if attr == "ratio_flow":  # at the surface the lines run together; the legend names them
                _end_labels(ax, np.array([series[(day, h)][-1] for h in THICKNESSES]),
                            [_cm(h) for h in THICKNESSES], log=True)
    axes[0, 1].legend(title="Layer", title_fontsize=7.5, loc="upper left", bbox_to_anchor=(1.12, 1.0))
    return _save(fig, out, revision, "A layer of molasses at rest: viscosity over the first 12 h")


def profile_figure(results: list[CaseResult], out: Path, revision: str) -> Path:
    """Viscosity against height in a 0.2 m layer at the profile times, both days."""
    from cooling.scenario import PROFILE_HOURS

    _, pp = _pp()
    fig, axes = _figure("cooling-profiles", (pp.FIG_WIDTH, 3.9), 1, 2, sharey=True)
    colours = _ramp(len(PROFILE_HOURS), start=0.35)
    for j, day in enumerate(DAYS):
        ax = axes[0, j]
        r = _by(results, day=day, thickness=PROFILE_THICKNESS, label="base")
        height = (PROFILE_THICKNESS - r.profile_depth) * 100.0
        ax.axvline(1.0, color=pp.INK_MUTED, lw=0.9, ls=(0, (1, 2)), zorder=1)
        ax.axvline(2.0, color=pp.INK_SECONDARY, lw=0.9, ls=pp.DASH, zorder=1)
        for h, ratio, c in zip(PROFILE_HOURS, r.profile_ratio, colours, strict=True):
            clock = (COLLAPSE_HOUR + h) % 24
            text = f"{h:g} h ({int(clock):02d}:{round(clock % 1 * 60):02d})"
            ax.plot(ratio, height, color=c, lw=2.0, label=text, solid_capstyle="round")
        ax.set_xscale("log")
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

        ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2, 4, 8]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"×{v:g}"))
        ax.set_ylim(0.0, PROFILE_THICKNESS * 100.0)
        ax.set_xlabel("Viscosity relative to the start")
        ax.set_title(DAYS[day].name, loc="left", fontsize=9, color=pp.INK, pad=16)
        ax.legend(title="After the collapse", title_fontsize=7.5, loc="lower left")
        # Direct labels at the free surface, spread along the axis
        lo, hi = (np.log10(v) for v in ax.get_xlim())
        tops = np.log10(r.profile_ratio[:, -1])  # cells run up from the ground
        for x, h in zip(_spread(tops, 0.07 * (hi - lo), lo, hi), PROFILE_HOURS, strict=True):
            ax.annotate(f"{h:g} h", (10**x, PROFILE_THICKNESS * 100.0), xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=7.5, color=pp.INK_SECONDARY, annotation_clip=False)
    axes[0, 0].set_ylabel("Height above the street [cm]")
    return _save(fig, out, revision, f"Viscosity through a {PROFILE_THICKNESS * 100:g} cm layer: "
                                     "the street cools the bottom at once, the sky the top after sunset")


def sensitivity_figure(results: list[CaseResult], out: Path, revision: str) -> Path:
    """Flow-weighted and surface viscosity at SENSITIVITY_HOUR for one change at a time."""
    _, pp = _pp()
    labels = [r.case.label for r in results
              if r.case.day == "winter" and r.case.thickness == SENSITIVITY_THICKNESS and r.case.path == "both"
              and r.case.label != "base"]
    fig, axes = _figure("cooling-sensitivity", (pp.FIG_WIDTH, 5.4), 1, 2, sharex=True, sharey=True)
    metrics = (("ratio_flow", "flow-weighted", TERM_COLOURS["convection"], "o"),
               ("ratio_surface", "at the surface", TERM_COLOURS["longwave"], "D"))
    y = np.arange(len(labels))[::-1]
    for j, day in enumerate(DAYS):
        ax = axes[0, j]
        base = _by(results, day=day, thickness=SENSITIVITY_THICKNESS, label="base")
        for attr, name, colour, marker in metrics:
            b = float(np.interp(SENSITIVITY_HOUR, base.hours, getattr(base, attr)))
            ax.axvline(b, color=colour, lw=1.2, ls=pp.DASH, zorder=1)
            vals = [float(np.interp(SENSITIVITY_HOUR, r.hours, getattr(r, attr)))
                    for r in (_by(results, day=day, thickness=SENSITIVITY_THICKNESS, label=lab) for lab in labels)]
            ax.hlines(y, b, vals, color=colour, lw=1.0, alpha=0.55, zorder=2)
            ax.plot(vals, y, ls="none", marker=marker, ms=5.5, color=colour, mec=pp.SURFACE, mew=1.0, zorder=3,
                    label=f"{name} (dashed: base case, ×{b:.2f})")
        ax.set_xscale("log")
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

        ax.xaxis.set_major_locator(FixedLocator([0.5, 0.75, 1, 1.5, 2, 3]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"×{v:g}"))
        ax.set_yticks(y, labels)
        ax.tick_params(axis="y", length=0)
        ax.grid(False, axis="y")
        ax.set_ylim(-0.7, len(labels) - 0.3)
        clock = COLLAPSE_HOUR + SENSITIVITY_HOUR
        ax.set_xlabel(f"Viscosity at {int(clock):02d}:{round(clock % 1 * 60):02d} relative to the start")
        ax.set_title(DAYS[day].name, loc="left", fontsize=9, color=pp.INK)
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=1, fontsize=7.5)
    return _save(fig, out, revision, f"One input at a time, {SENSITIVITY_THICKNESS * 100:g} cm layer, "
                                     f"{SENSITIVITY_HOUR:g} h after the collapse")


def terms_figure(results: list[CaseResult], out: Path, revision: str) -> Path:
    """Terms of the surface balance against time for one thickness, both days."""
    _, pp = _pp()
    fig, axes = _figure("cooling-surface-terms", (pp.FIG_WIDTH, 3.9), 1, 2, sharey=True)
    names = {"convection": "convection", "longwave": "net long-wave", "shortwave": "absorbed sun (gain)"}
    for j, day in enumerate(DAYS):
        ax = axes[0, j]
        r = _by(results, day=day, thickness=TERMS_THICKNESS, label="base")
        t = r.surface_terms
        flux = {"convection": t["convection"], "longwave": t["longwave"], "shortwave": -t["shortwave"]}
        net = t["convection"] + t["longwave"] - t["shortwave"]
        ax.axhline(0.0, color=pp.AXIS, lw=0.9, zorder=1)
        for key, values in flux.items():
            ax.plot(r.hours, values, color=TERM_COLOURS[key], lw=2.0, label=names[key], solid_capstyle="round")
        ax.plot(r.hours, net, color=pp.INK, lw=2.0, label="net loss", solid_capstyle="round")
        _sunset(ax, day)
        _clock_axis(ax)
        ax.set_title(DAYS[day].name, loc="left", fontsize=9, color=pp.INK)
        ends = np.array([flux["convection"][-1], flux["longwave"][-1], net[-1]])
        _end_labels(ax, ends, ["convection", "long-wave", "net"], log=False)
        k = int(np.argmin(flux["shortwave"]))
        ax.annotate("sun", (r.hours[k], flux["shortwave"][k]), xytext=(6, -2), textcoords="offset points",
                    va="top", fontsize=7.5, color=pp.INK_SECONDARY)
    axes[0, 0].set_ylabel("Heat flux out of the surface [W/m²]\n(negative: into the molasses)")
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.12, 1.0))
    return _save(fig, out, revision, f"Surface heat balance of a {TERMS_THICKNESS * 100:g} cm layer (base case)")


def write_all(results: list[CaseResult], out_dir: Path, revision: str = "") -> list[Path]:
    """Write the four figures into ``out_dir``."""
    return [
        viscosity_figure(results, out_dir / "cooling-viscosity.svg", revision),
        profile_figure(results, out_dir / "cooling-profiles.svg", revision),
        sensitivity_figure(results, out_dir / "cooling-sensitivity.svg", revision),
        terms_figure(results, out_dir / "cooling-surface-terms.svg", revision),
    ]
