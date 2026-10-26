"""Panels and figure layout of the dam-break animation frames.

A frame is a vertical stack of panels sharing the x axis:

- "strip" panels show a field over the x-y plane at the true aspect ratio
  (no vertical exaggeration); the y range is cropped to the initial column
  height.
- "line" panels show quantities along x (depth, Froude number) with
  analytic references.

Panels are registered in ``PANELS`` and selected by name, so a layout is
just a list of names (e.g. the README GIF uses flow, temperature, depth).
Sizes are given in inches and fonts in points: rendering the same figure at
a different dpi scales it uniformly.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import matplotlib as mpl
import numpy as np
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Colormap, LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from post.analytic import ritter
from post.frames import CaseInfo, Frame

KELVIN = 273.15

# Chart chrome (light surface)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"
NEUTRAL = "#f0efec"  # air, and regions without a valid reference solution
SERIES = "#2a78d6"  # simulation
REFERENCE = INK_SECONDARY  # analytic solutions, dashed
DASH = (0, (4, 3))

# Sequential one-hue ramps, light -> dark
SPEED_CMAP = LinearSegmentedColormap.from_list(
    "speed", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]
)
TEMPERATURE_CMAP = LinearSegmentedColormap.from_list(
    "temperature", ["#fbe6dc", "#f4b393", "#eb6834", "#c24e1e", "#8a3510"]
)

# Layout [in]
FIG_WIDTH = 9.6
MARGIN_LEFT = 0.62
MARGIN_RIGHT = 0.88
HEADER = 0.56
STRIP_TITLE = 0.28
PANEL_GAP = 0.14
LINE_HEIGHT = 1.25
FOOTER = 0.5
CBAR_GAP = 0.12
CBAR_WIDTH = 0.1

FROUDE_MAX = 4.0
FROUDE_MIN_DEPTH = 0.1  # [m] below this the Froude number is not shown

_RC = {
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.labelcolor": INK_SECONDARY,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.facecolor": SURFACE,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.color": INK_MUTED,
    "ytick.color": INK_MUTED,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "legend.labelcolor": INK_SECONDARY,
}


@dataclass(frozen=True)
class ColourScale:
    """Colour mapping of a strip panel, drawn as its colour bar."""

    cmap: Colormap
    norm: Normalize
    label: str


@dataclass(frozen=True)
class Panel:
    """A registered panel.

    Attributes:
        kind: "strip" (x-y field at true aspect ratio) or "line" (along x).
        draw: Draws the panel; strip panels return their colour scale.
        title: Caption above a strip panel.
        needs_temperature: Whether the panel reads the temperature field.
    """

    kind: Literal["strip", "line"]
    draw: Callable[[Axes, Frame, CaseInfo], ColourScale | None]
    title: str = ""
    needs_temperature: bool = False


def strip_top(info: CaseInfo) -> float:
    """Upper y limit of the strip panels: the initial column height, rounded up."""
    return min(float(np.ceil(info.h0)), info.grid.extent[3])


def _blend_liquid(alpha: np.ndarray, rgba: np.ndarray) -> np.ndarray:
    """Blend liquid colours with the neutral air colour by the volume fraction."""
    a = alpha[..., None]
    return a * rgba[..., :3] + (1.0 - a) * np.array(to_rgb(NEUTRAL))


def _draw_strip(ax: Axes, rgb: np.ndarray, frame: Frame, info: CaseInfo) -> None:
    grid = info.grid
    ax.imshow(rgb, origin="lower", extent=grid.extent, interpolation="nearest", aspect="equal")
    if frame.alpha.min() < 0.5 < frame.alpha.max():
        ax.contour(grid.x, grid.y, frame.alpha, levels=[0.5], colors=INK, linewidths=0.8)
    ax.set_xlim(grid.extent[0], grid.extent[1])
    ax.set_ylim(grid.extent[2], strip_top(info))
    ax.set_ylabel("y [m]")


def draw_flow(ax: Axes, frame: Frame, info: CaseInfo) -> ColourScale:
    """Liquid coloured by the speed |U|; air is not shown."""
    norm = Normalize(0.0, 2.0 * np.ceil(info.front_speed / 2.0))
    _draw_strip(ax, _blend_liquid(frame.alpha, SPEED_CMAP(norm(frame.speed))), frame, info)
    return ColourScale(SPEED_CMAP, norm, "|U| [m/s]")


def draw_temperature(ax: Axes, frame: Frame, info: CaseInfo) -> ColourScale:
    """Liquid coloured by temperature; air is not shown.

    With a single mixture temperature, air in cells that held traces of liquid
    takes the liquid temperature while holding a negligible amount of heat, so
    showing T in the air would highlight a model artefact.
    """
    lo, hi = sorted((info.t_air - KELVIN, info.t_liquid - KELVIN))
    norm = Normalize(lo, hi)
    rgba = TEMPERATURE_CMAP(norm(frame.temperature - KELVIN))
    _draw_strip(ax, _blend_liquid(frame.alpha, rgba), frame, info)
    return ColourScale(TEMPERATURE_CMAP, norm, "T [°C]")


def _style_line_panel(ax: Axes, info: CaseInfo) -> None:
    ax.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(info.grid.extent[0], info.grid.extent[1])


def _ritter(frame: Frame, info: CaseInfo, x: np.ndarray):
    return ritter(x, frame.time, info.h0, info.x_wall, info.x_dam, info.g)


def _reflection_band(ax: Axes, info: CaseInfo, x_reflected: float) -> Patch:
    if x_reflected > info.x_wall:
        ax.axvspan(info.x_wall, x_reflected, color=NEUTRAL, lw=0, zorder=0)
    return Patch(facecolor=NEUTRAL, edgecolor=AXIS, lw=0.6,
                 label="Ritter not valid (wave reflected at the wall)")


def _line_legend(ax: Axes, band: Patch) -> None:
    handles = [
        Line2D([], [], color=SERIES, lw=1.8, label="simulation"),
        Line2D([], [], color=REFERENCE, lw=1.2, ls=DASH, label="Ritter (shallow water, analytic)"),
        band,
    ]
    ax.legend(handles=handles, loc="upper right")


def draw_depth(ax: Axes, frame: Frame, info: CaseInfo) -> None:
    """Liquid depth along x against Ritter's solution."""
    _style_line_panel(ax, info)
    x = np.linspace(info.grid.extent[0], info.grid.extent[1], 4 * info.grid.nx + 1)
    sol = _ritter(frame, info, x)
    band = _reflection_band(ax, info, sol.x_reflected)
    ax.plot(x, sol.depth, color=REFERENCE, lw=1.2, ls=DASH)
    ax.plot(info.grid.x, frame.depth(), color=SERIES, lw=1.8)
    ax.set_ylim(0.0, 1.05 * info.h0)
    ax.set_ylabel("depth h [m]")
    _line_legend(ax, band)


def draw_froude(ax: Axes, frame: Frame, info: CaseInfo) -> None:
    """Froude number of the depth-averaged flow against Ritter's solution."""
    _style_line_panel(ax, info)
    x = np.linspace(info.grid.extent[0], info.grid.extent[1], 4 * info.grid.nx + 1)
    sol = _ritter(frame, info, x)
    band = _reflection_band(ax, info, sol.x_reflected)
    with np.errstate(invalid="ignore", divide="ignore"):
        fr_ref = np.where(sol.depth > 1e-3, sol.velocity / np.sqrt(info.g * sol.depth), np.nan)
    ax.axhline(1.0, color=INK_MUTED, lw=0.8)
    ax.axvline(info.x_dam, color=INK_MUTED, lw=0.8)
    ax.text(info.x_dam + 1.5, 0.04 * FROUDE_MAX, "dam site", color=INK_SECONDARY, va="bottom")
    ax.plot(x, fr_ref, color=REFERENCE, lw=1.2, ls=DASH)
    ax.plot(info.grid.x, frame.froude(info.g, FROUDE_MIN_DEPTH), color=SERIES, lw=1.8)
    ax.set_ylim(0.0, FROUDE_MAX)
    ax.set_ylabel("Froude number")
    _line_legend(ax, band)


PANELS: dict[str, Panel] = {
    "flow": Panel("strip", draw_flow, "Liquid coloured by speed |U| (air not shown)"),
    "temperature": Panel(
        "strip",
        draw_temperature,
        "Liquid coloured by temperature T (air not shown)",
        needs_temperature=True,
    ),
    "depth": Panel("line", draw_depth),
    "froude": Panel("line", draw_froude),
}


def describe_case(info: CaseInfo) -> tuple[str, str]:
    """Default title and subtitle of the frames."""
    step = f"Step {info.step}" if info.step is not None else info.path.name
    mu0 = (
        f"{info.mu0_liquid * 1e3:.4g} mPa·s" if info.mu0_liquid < 0.1 else f"{info.mu0_liquid:.4g} Pa·s"
    )
    title = f"2D dam break — {step}"
    subtitle = (
        f"liquid ρ = {info.rho_liquid:g} kg/m³, μ₀ = {mu0}"
        f" · energy equation {'on' if info.solve_energy else 'off'}"
        f" · viscosity {'μ(T)' if info.temperature_dependent else 'constant'}"
        f" · 2D VOF, Δ = {info.grid.dx:g} m"
    )
    return title, subtitle


def render_frame(
    info: CaseInfo,
    frame: Frame,
    panels: Sequence[str],
    path: Path,
    dpi: float,
    title: str,
    subtitle: str,
) -> None:
    """Render one frame with the named panels and save it as PNG."""
    x0, x1, y0, _ = info.grid.extent
    plot_w = FIG_WIDTH - MARGIN_LEFT - MARGIN_RIGHT
    strip_h = plot_w * (strip_top(info) - y0) / (x1 - x0)

    # (space above, height) of each panel, from the top
    slots = []
    for name in panels:
        if PANELS[name].kind == "strip":
            slots.append((STRIP_TITLE, strip_h))
        else:
            slots.append((PANEL_GAP, LINE_HEIGHT))
    # 0.01 in steps keep the pixel height even at 200 dpi (yuv420p needs it)
    fig_h = round(HEADER + sum(a + h for a, h in slots) + FOOTER, 2)

    # Figure without pyplot: no GUI backend, safe in worker processes
    with mpl.rc_context(_RC):
        fig = Figure(figsize=(FIG_WIDTH, fig_h), dpi=dpi, facecolor=SURFACE)
        fig.text(MARGIN_LEFT / FIG_WIDTH, 1 - 0.2 / fig_h, title,
                 color=INK, fontsize=11, fontweight="bold", va="center")
        fig.text(MARGIN_LEFT / FIG_WIDTH, 1 - 0.42 / fig_h, subtitle,
                 color=INK_SECONDARY, fontsize=8, va="center")
        fig.text((MARGIN_LEFT + plot_w) / FIG_WIDTH, 1 - 0.2 / fig_h, f"t = {frame.time:5.2f} s",
                 color=INK, fontsize=11, family="monospace", ha="right", va="center")

        top = HEADER
        axes = []
        for name, (above, height) in zip(panels, slots):
            top += above
            bottom = (fig_h - top - height) / fig_h
            ax = fig.add_axes((MARGIN_LEFT / FIG_WIDTH, bottom, plot_w / FIG_WIDTH, height / fig_h))
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
            ax.tick_params(length=3, width=0.6)

            panel = PANELS[name]
            scale = panel.draw(ax, frame, info)
            if panel.title:
                ax.set_title(panel.title, loc="left", color=INK_SECONDARY, fontsize=8.5, pad=3)
            if scale is not None:
                cax = fig.add_axes((
                    (MARGIN_LEFT + plot_w + CBAR_GAP) / FIG_WIDTH,
                    bottom,
                    CBAR_WIDTH / FIG_WIDTH,
                    height / fig_h,
                ))
                cbar = fig.colorbar(ScalarMappable(scale.norm, scale.cmap), cax=cax)
                # Strips are short: label both ends and the middle only
                cbar.set_ticks(np.linspace(scale.norm.vmin, scale.norm.vmax, 3))
                cbar.ax.yaxis.set_major_formatter("{x:.3g}")
                cbar.set_label(scale.label, color=INK_SECONDARY, fontsize=8)
                cbar.ax.tick_params(labelsize=7, length=2, width=0.6)
                cbar.outline.set_visible(False)
            axes.append(ax)
            top += height

        for ax in axes[:-1]:
            ax.tick_params(labelbottom=False)
        axes[-1].set_xlabel("x [m]")

        fig.savefig(path, dpi=dpi, facecolor=SURFACE)
