"""Render a dam-break case into PNG frames, an mp4 video and a GIF.

Usage (from the project root; ``make anim`` wraps this):

    uv run python -m post.animate CASE_DIR [--panels flow,temperature,depth] [...]

The case needs the cell centres in ``0/C`` (``postProcess -func
writeCellCentres -time 0``). Frames are rendered once at ``--dpi`` (1920 px
wide at the default 200 dpi) into ``OUT/frames``. The mp4 (H.264, e.g. for
slides) plays them at ``--speed`` times real time; the GIF (for the README) is
derived from the same frames, downscaled to ``--gif-width`` and decimated to
``--gif-fps`` at the same speed. Both hold the last frame for ``--hold`` s.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from post.frames import CaseInfo, load_case, load_frame
from post.plot import PANELS, describe_case, render_frame

DEFAULT_PANELS = "flow,temperature,depth"
GIF_SIZE_BUDGET = 10 * 1024**2  # [bytes] keep README GIFs light
FRAME_PATTERN = "frame_%05d.png"

# Per-process render settings, set once by the pool initializer
_settings: dict = {}


def _init_worker(info: CaseInfo, panels: list[str], frames_dir: Path, dpi: float,
                 title: str, subtitle: str) -> None:
    _settings.update(
        info=info,
        panels=panels,
        frames_dir=frames_dir,
        dpi=dpi,
        title=title,
        subtitle=subtitle,
        temperature=any(PANELS[p].needs_temperature for p in panels),
    )


def _render(job: tuple[int, str]) -> None:
    index, time_name = job
    s = _settings
    frame = load_frame(s["info"], time_name, temperature=s["temperature"])
    render_frame(s["info"], frame, s["panels"], s["frames_dir"] / (FRAME_PATTERN % index),
                 s["dpi"], s["title"], s["subtitle"])


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(prog="python -m post.animate", description=__doc__.split("\n")[0])
    p.add_argument("case", type=Path, help="case directory with results")
    p.add_argument("--panels", default=DEFAULT_PANELS,
                   help=f"comma-separated panels, top to bottom: {', '.join(PANELS)} "
                        f"(default: {DEFAULT_PANELS})")
    p.add_argument("--out", type=Path, help="output directory (default: CASE/postProcessing/anim)")
    p.add_argument("--name", default="animation", help="base name of the mp4/GIF (default: animation)")
    p.add_argument("--t-max", type=float, help="last time to include [s] (default: all)")
    p.add_argument("--speed", type=float, default=1.0,
                   help="playback speed relative to simulated time (default: 1)")
    p.add_argument("--dpi", type=float, default=200, help="frame resolution (default: 200 -> 1920 px)")
    p.add_argument("--gif-fps", type=float, default=12.5, help="GIF frame rate (default: 12.5)")
    p.add_argument("--gif-width", type=int, default=960, help="GIF width [px] (default: 960)")
    p.add_argument("--hold", type=float, default=1.5, help="hold the last frame [s] (default: 1.5)")
    p.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1),
                   help="parallel render processes")
    p.add_argument("--title", help="title (default: derived from the case)")
    p.add_argument("--subtitle", help="subtitle (default: derived from the case)")
    return p.parse_args(argv)


def _ffmpeg(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def encode_mp4(ffmpeg: str, frames_dir: Path, fps: float, hold: float, out: Path) -> None:
    """Encode the frames as H.264 mp4 (yuv420p, playable in browsers)."""
    _ffmpeg([
        ffmpeg, "-y", "-loglevel", "error",
        "-framerate", f"{fps:g}", "-i", str(frames_dir / FRAME_PATTERN),
        "-vf", f"tpad=stop_mode=clone:stop_duration={hold:g},format=yuv420p",
        "-c:v", "libx264", "-preset", "slow", "-crf", "18", "-movflags", "+faststart",
        str(out),
    ])


def encode_gif(ffmpeg: str, frames_dir: Path, fps: float, gif_fps: float, width: int,
               hold: float, out: Path) -> None:
    """Encode the frames as a looping GIF with a palette computed from all frames."""
    graph = (
        f"fps={gif_fps:g},scale={width}:-1:flags=lanczos,"
        f"tpad=stop_mode=clone:stop_duration={hold:g},split[a][b];"
        "[a]palettegen=stats_mode=full[p];[b][p]paletteuse=dither=sierra2_4a"
    )
    _ffmpeg([
        ffmpeg, "-y", "-loglevel", "error",
        "-framerate", f"{fps:g}", "-i", str(frames_dir / FRAME_PATTERN),
        "-filter_complex", graph, "-loop", "0", str(out),
    ])


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        print("Error: ffmpeg not found on PATH", file=sys.stderr)
        return 1

    panels = [p.strip() for p in args.panels.split(",") if p.strip()]
    unknown = [p for p in panels if p not in PANELS]
    if unknown:
        print(f"Error: unknown panel(s) {', '.join(unknown)}; choose from {', '.join(PANELS)}",
              file=sys.stderr)
        return 1

    try:
        info = load_case(args.case)
    except (OSError, KeyError, ValueError, NotImplementedError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 1
    if not info.solve_energy:
        dropped = [p for p in panels if PANELS[p].needs_temperature]
        if dropped:
            print(f"Warning: the energy equation is off, T is the frozen initial field; "
                  f"dropping {', '.join(dropped)}", file=sys.stderr)
            panels = [p for p in panels if p not in dropped]
    if not panels:
        print("Error: no panels left to draw", file=sys.stderr)
        return 1

    times = [t for t in info.times if args.t_max is None or float(t) <= args.t_max + 1e-9]
    if len(times) < 2:
        print(f"Error: {args.case} has fewer than two time directories to animate", file=sys.stderr)
        return 1
    steps = np.diff([float(t) for t in times])
    interval = float(np.median(steps))
    if not np.allclose(steps, interval, rtol=1e-3):
        print("Warning: the write interval is not uniform; playback speed will vary",
              file=sys.stderr)
    fps = args.speed / interval
    gif_fps = min(args.gif_fps, fps)

    out = args.out or args.case / "postProcessing" / "anim"
    frames_dir = out / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    for old in frames_dir.glob("frame_*.png"):
        old.unlink()

    title, subtitle = describe_case(info)
    title = args.title or title
    subtitle = args.subtitle or subtitle

    print(f"Rendering {len(times)} frames ({times[0]}-{times[-1]} s, panels: {', '.join(panels)}) "
          f"with {args.jobs} processes -> {frames_dir}", flush=True)
    with ProcessPoolExecutor(
        max_workers=args.jobs,
        initializer=_init_worker,
        initargs=(info, panels, frames_dir, args.dpi, title, subtitle),
    ) as pool:
        for done, _ in enumerate(pool.map(_render, enumerate(times), chunksize=4), start=1):
            if done % 50 == 0 or done == len(times):
                print(f"  {done}/{len(times)}", flush=True)

    mp4 = out / f"{args.name}.mp4"
    gif = out / f"{args.name}.gif"
    print(f"Encoding {mp4} ({fps:g} fps, {args.speed:g}x real time)", flush=True)
    encode_mp4(ffmpeg, frames_dir, fps, args.hold, mp4)
    print(f"Encoding {gif} ({gif_fps:g} fps, {args.gif_width} px wide)", flush=True)
    encode_gif(ffmpeg, frames_dir, fps, gif_fps, args.gif_width, args.hold, gif)

    for path in (mp4, gif):
        print(f"  {path}: {path.stat().st_size / 1024**2:.1f} MB")
    if gif.stat().st_size > GIF_SIZE_BUDGET:
        print(f"Warning: the GIF exceeds {GIF_SIZE_BUDGET / 1024**2:.0f} MB; consider a lower "
              "--gif-fps, --gif-width or --t-max", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
