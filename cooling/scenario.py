"""Scenario runs of the cooling column (#9).

A layer of molasses at rest on the street, from the collapse at 12:40 EST for
12 h, on 15 January 1919 and on 15 July, for layer thicknesses from 2 cm to
1 m. Each base run is repeated with only the surface exchange (insulated
bottom) and with only the ground exchange (insulated surface), and with
one uncertain input at a time moved to the end of its range. The run
conditions are listed on the wiki page Cooling-column-inputs, "Values
adopted for the column".

The viscosity is reported as a ratio to the viscosity at the start,
mu(T) / mu(T_start) = exp(B (1/T - 1/T_start)) with the Andrade constant B of
Steps 4-5, so it does not depend on mu0; absolute values scale with mu0.
Besides the viscosity at the surface and at the ground, two depth averages:

- mean: the arithmetic mean of mu over the depth;
- flow: the viscosity of a uniform layer that would carry the same flux in
  a slow gravity-driven flow, where the shear stress grows linearly with
  depth: H^3 / (3 * integral (H - z)^2 / mu dz), z up from the ground. It
  weights the bottom, where the shear is largest.

At t = 0 the contact is given the temperature of two semi-infinite bodies
brought together, (e_m T_m + e_g T_g) / (e_m + e_g) with the effusivities
e, which is what the molasses at the ground takes at once; a doubling time
of 0 means the viscosity there doubles on contact.

B was fitted between 0 and 10 degC (wiki, Molasses-properties): below 0 degC
the ratios are lower bounds (the thickening "becomes more extreme"), and the
July runs extrapolate it to 24-31 degC.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np

from cooling.analytic import contact_temperature
from cooling.column import Adiabatic, Column, ColumnResult, Layer
from cooling.properties import CONCRETE, GRANITE, MOLASSES, MOLASSES_B, MOLASSES_MU0, MOLASSES_T0, andrade_viscosity
from cooling.sun import cloud_factor
from cooling.surface import BulkConvection, ClearConvection, McAdamsConvection, SurfaceBalance
from cooling.weather import DAYS

HOUR = 3600.0
DURATION = 12.0 * HOUR
TIME_STEP = 30.0  # [s]
OUTPUT_EVERY = 600.0  # [s]
PROFILE_HOURS = (1.0, 3.0, 6.0, 12.0)
REPORT_HOURS = (1.0, 3.0, 6.0, 12.0)
THICKNESSES = (0.02, 0.05, 0.1, 0.2, 0.5, 1.0)  # [m]
GROUND_DEPTH = 2.0  # [m], insulated below
GRANITE_DEPTH = 0.15  # [m]
PATHS = ("both", "surface", "ground")


@dataclass(frozen=True)
class Case:
    """One run of the column. Temperatures in kelvin; None takes the value of the day."""

    day: str
    thickness: float
    path: str = "both"
    convection: str = "clear"
    fetch: float = 30.0
    wind: float = 2.0
    sky: str = "swinbank"
    cloud: float = 0.0
    emissivity: float = 0.95
    albedo: float = 0.07
    ground: str = "granite_concrete"
    ground_initial: float | None = None
    label: str = "base"

    def __post_init__(self) -> None:
        if self.day not in DAYS:
            raise ValueError(f"unknown day {self.day!r}")
        if self.path not in PATHS:
            raise ValueError(f"unknown path {self.path!r}")


def _convection(case: Case):
    if case.convection == "clear":
        return ClearConvection(case.fetch)
    if case.convection == "bulk":
        return BulkConvection()
    if case.convection == "mcadams":
        return McAdamsConvection()
    raise ValueError(f"unknown convection {case.convection!r}")


def _ground_layers(case: Case, temperature: float) -> list[Layer]:
    if case.ground == "granite_concrete":
        return [Layer(CONCRETE, GROUND_DEPTH - GRANITE_DEPTH, temperature), Layer(GRANITE, GRANITE_DEPTH, temperature)]
    if case.ground == "granite":
        return [Layer(GRANITE, GROUND_DEPTH, temperature)]
    if case.ground == "concrete":
        return [Layer(CONCRETE, GROUND_DEPTH, temperature)]
    raise ValueError(f"unknown ground {case.ground!r}")


def build_column(case: Case) -> Column:
    """The column of ``case``: ground layers (unless the path is surface only) under the molasses."""
    day = DAYS[case.day]
    molasses = Layer(MOLASSES, case.thickness, day.molasses_initial)
    if case.path == "ground":
        surface = Adiabatic()
    else:
        surface = SurfaceBalance(
            air=day.air,
            vapour_pressure=day.vapour_pressure,
            solar=day.solar,
            wind=case.wind,
            convection=_convection(case),
            sky_model=case.sky,
            cloud=case.cloud,
            emissivity=case.emissivity,
            albedo=case.albedo,
            cloud_transmission=cloud_factor(case.cloud),
        )
    if case.path == "surface":
        return Column([molasses], surface, Adiabatic())
    ground_t = day.ground_initial if case.ground_initial is None else case.ground_initial
    return Column([*_ground_layers(case, ground_t), molasses], surface, Adiabatic())


@dataclass(frozen=True)
class CaseResult:
    """Time series of a run (at ``hours``) and profiles of the molasses at ``PROFILE_HOURS``."""

    case: Case
    hours: np.ndarray
    air: np.ndarray  # [K]
    surface: np.ndarray  # temperature of the free surface [K]
    contact: np.ndarray  # temperature at the bottom of the molasses [K]
    mean_temperature: np.ndarray  # depth mean of the molasses [K]
    ratio_surface: np.ndarray
    ratio_contact: np.ndarray
    ratio_mean: np.ndarray
    ratio_flow: np.ndarray
    profile_depth: np.ndarray  # below the free surface [m]
    profile_ratio: np.ndarray  # shape (PROFILE_HOURS, cells)
    profile_temperature: np.ndarray
    min_temperature: float  # lowest molasses temperature in the run [K]
    max_temperature: float  # highest [K]
    budget_error: float  # |change of heat content - heat in| / |heat lost by the molasses|, at the end
    surface_terms: dict[str, np.ndarray]  # W/m2 at ``hours`` (empty for the ground path)

    def doubling_hours(self, ratio: np.ndarray) -> float | None:
        """First time [h] at which ``ratio`` (a series of this run) reaches 2."""
        return doubling_hours(self.hours, ratio)


def doubling_hours(hours: np.ndarray, ratio: np.ndarray) -> float | None:
    """First time at which ``ratio`` reaches 2, interpolated between samples; None if it never does."""
    above = np.flatnonzero(ratio >= 2.0)
    if above.size == 0:
        return None
    i = above[0]
    if i == 0:
        return float(hours[0])
    r0, r1 = ratio[i - 1], ratio[i]
    return float(hours[i - 1] + (2.0 - r0) / (r1 - r0) * (hours[i] - hours[i - 1]))


def viscosity_ratio(temperature: np.ndarray, start: float) -> np.ndarray:
    """mu(T) / mu(T_start) for the Andrade viscosity of the molasses."""
    mu0 = MOLASSES_MU0["central"]  # cancels
    return andrade_viscosity(temperature, mu0, MOLASSES_B, MOLASSES_T0) / float(
        andrade_viscosity(start, mu0, MOLASSES_B, MOLASSES_T0)
    )


def flow_weighted(ratio: np.ndarray, height: np.ndarray, widths: np.ndarray, thickness: float) -> np.ndarray:
    """H^3 / (3 integral (H - z)^2 / r dz), per row of ``ratio`` (cells along the last axis).

    H^3 / 3 is taken as the same cell sum with r = 1, so that a uniform layer
    returns its own ratio exactly on any grid.
    """
    weights = (thickness - height) ** 2 * widths
    return np.sum(weights) / np.sum(weights / ratio, axis=-1)


def run_case(case: Case) -> CaseResult:
    """Run ``case`` for 12 h from the collapse and reduce it to the reported quantities."""
    day = DAYS[case.day]
    column = build_column(case)
    outputs = np.arange(OUTPUT_EVERY, DURATION, OUTPUT_EVERY)
    result: ColumnResult = column.run(DURATION, TIME_STEP, outputs)
    molasses_layer = len(column.layers) - 1
    depth, cells = result.depth_in_layer(molasses_layer)
    height = case.thickness - depth  # above the ground
    widths = result.widths[cells]
    temp = result.temperature[:, cells]
    start = day.molasses_initial

    if case.path == "surface":
        contact = temp[:, 0]  # the bottom is insulated: face and cell agree
    else:
        contact = result.interface_temperature[:, -1].copy()
        ground = column.layers[-2]
        contact[0] = contact_temperature(
            start, MOLASSES.effusivity, ground.initial_temperature, ground.material.effusivity
        )
    ratio_cells = viscosity_ratio(temp, start)
    hours = result.times / HOUR
    keep = [int(np.argmin(np.abs(hours - h))) for h in PROFILE_HOURS]

    terms: dict[str, np.ndarray] = {}
    if case.path != "ground":
        boundary = column.surface
        rows = [boundary.terms(t, ts) for t, ts in zip(result.times, result.surface_temperature, strict=True)]
        terms = {
            "convection": np.array([r.convection for r in rows]),
            "longwave": np.array([r.longwave for r in rows]),
            "shortwave": np.array([r.shortwave for r in rows]),
        }

    change = result.heat_content[-1] - result.heat_content[0]
    lost = np.sum(column.capacity[cells] * (temp[0] - temp[-1]))
    budget = abs(change - result.heat_in[-1]) / abs(lost)
    return CaseResult(
        case=case,
        hours=hours,
        air=np.array([day.air(t) for t in result.times]),
        surface=result.surface_temperature,
        contact=contact,
        mean_temperature=np.sum(temp * widths, axis=1) / case.thickness,
        ratio_surface=viscosity_ratio(result.surface_temperature, start),
        ratio_contact=viscosity_ratio(contact, start),
        ratio_mean=np.sum(ratio_cells * widths, axis=1) / case.thickness,
        ratio_flow=flow_weighted(ratio_cells, height, widths, case.thickness),
        profile_depth=depth,
        profile_ratio=ratio_cells[keep],
        profile_temperature=temp[keep],
        min_temperature=float(min(temp.min(), result.surface_temperature.min(), contact.min())),
        max_temperature=float(max(temp.max(), result.surface_temperature.max(), contact.max())),
        budget_error=float(budget),
        surface_terms=terms,
    )


def sensitivity_cases(day: str, thickness: float) -> list[Case]:
    """The base case and one-at-a-time changes of the uncertain inputs."""
    base = Case(day, thickness)
    low_g, high_g = DAYS[day].ground_range
    changes = [
        ("wind 1 m/s", {"wind": 1.0}),
        ("wind 4 m/s", {"wind": 4.0}),
        ("fetch 10 m", {"fetch": 10.0}),
        ("fetch 100 m", {"fetch": 100.0}),
        ("bulk convection", {"convection": "bulk"}),
        ("McAdams (upper bound)", {"convection": "mcadams"}),
        ("sky Idso & Jackson", {"sky": "idso_jackson"}),
        ("sky Prata", {"sky": "prata"}),
        ("cloud 0.65", {"cloud": 0.65}),
        ("emissivity 0.93", {"emissivity": 0.93}),
        ("emissivity 0.98", {"emissivity": 0.98}),
        ("albedo 0.04", {"albedo": 0.04}),
        ("albedo 0.10", {"albedo": 0.10}),
        ("ground cold", {"ground_initial": low_g}),
        ("ground warm", {"ground_initial": high_g}),
        ("ground all granite", {"ground": "granite"}),
        ("ground all concrete", {"ground": "concrete"}),
    ]
    return [base] + [replace(base, label=label, **change) for label, change in changes]


def all_cases() -> list[Case]:
    cases: list[Case] = []
    for day in DAYS:
        for h in THICKNESSES:
            cases += sensitivity_cases(day, h)
            cases += [Case(day, h, path=p, label=f"{p} only") for p in PATHS[1:]]
    return cases


def _at_hours(result: CaseResult, values: np.ndarray, digits: int = 3) -> dict[str, float]:
    return {f"{h:g} h": round(float(np.interp(h, result.hours, values)), digits) for h in REPORT_HOURS}


def _case_entry(case: Case) -> dict:
    """Day, thickness, path, label and the inputs that differ from the base case."""
    base = asdict(Case(case.day, case.thickness))
    keep = ("day", "thickness", "path", "label")
    return {k: v for k, v in asdict(case).items() if k in keep or v != base[k]}


def start_viscosity(day: str) -> dict[str, float]:
    """Viscosity of the molasses at the start of ``day`` for the mu0 levels of Steps 4-5 [Pa s]."""
    t = DAYS[day].molasses_initial
    return {k: round(float(andrade_viscosity(t, mu0, MOLASSES_B, MOLASSES_T0)), 2) for k, mu0 in MOLASSES_MU0.items()}


def summarise(results: list[CaseResult], revision: str = "") -> dict:
    """Summary of all runs for docs/cooling/summary.json (time series only for the base runs)."""
    runs = []
    for r in results:
        entry = {
            "case": _case_entry(r.case),
            "molasses_start_C": round(DAYS[r.case.day].molasses_initial - 273.15, 2),
            "min_temperature_C": round(r.min_temperature - 273.15, 3),
            "max_temperature_C": round(r.max_temperature - 273.15, 3),
            "below_fit_range": bool(r.min_temperature < 273.15),
            "budget_error": float(f"{r.budget_error:.2e}"),
            "doubling_hours": {
                k: (None if (v := r.doubling_hours(getattr(r, f"ratio_{k}"))) is None else round(v, 3))
                for k in ("surface", "contact", "mean", "flow")
            },
            "ratio": {k: _at_hours(r, getattr(r, f"ratio_{k}")) for k in ("surface", "contact", "mean", "flow")},
            "temperature_C": {k: _at_hours(r, getattr(r, k) - 273.15, 2) for k in ("surface", "contact")},
        }
        if r.case.label == "base":
            entry["series"] = {
                "hours": [round(float(x), 4) for x in r.hours],
                **{k: [round(float(x), 3) for x in getattr(r, f"ratio_{k}")] for k in ("surface", "contact", "flow")},
                "surface_C": [round(float(x) - 273.15, 2) for x in r.surface],
                "contact_C": [round(float(x) - 273.15, 2) for x in r.contact],
            }
        runs.append(entry)
    return {
        "description": "Cooling column of #9: 12 h from 12:40 EST; viscosity as a ratio to the start (Andrade, B = 8.5e3 K)",
        "revision": revision,
        "duration_h": DURATION / HOUR,
        "time_step_s": TIME_STEP,
        "output_every_s": OUTPUT_EVERY,
        "thicknesses_m": list(THICKNESSES),
        "base_case": {k: v for k, v in asdict(Case("winter", THICKNESSES[0])).items()
                      if k not in ("day", "thickness", "label")},
        "ground_layers": f"granite {GRANITE_DEPTH} m on concrete to {GROUND_DEPTH} m, insulated below",
        "start_viscosity_Pa_s": {day: start_viscosity(day) for day in DAYS},
        "runs": runs,
    }


def dump_summary(summary: dict) -> str:
    """JSON text of ``summary`` with one line per run, so that a change of the runs reads as a line diff."""
    head = {k: v for k, v in summary.items() if k != "runs"}
    lines = json.dumps(head, indent=1)[:-2]
    runs = ",\n  ".join(json.dumps(run, separators=(",", ":")) for run in summary["runs"])
    return f'{lines},\n "runs": [\n  {runs}\n ]\n}}\n'


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", type=Path, help="write the summary of all runs here")
    parser.add_argument("--figures", type=Path, help="write the SVG figures into this directory")
    parser.add_argument("--revision", default="", help="label of the code revision, for the summary and the figures")
    parser.add_argument("--jobs", type=int, default=1, help="runs in parallel")
    args = parser.parse_args(argv)

    cases = all_cases()
    if args.jobs > 1:
        with ProcessPoolExecutor(args.jobs) as pool:
            results = list(pool.map(run_case, cases, chunksize=4))
    else:
        results = [run_case(c) for c in cases]
    worst = max(r.budget_error for r in results)
    if not worst < 1e-6:
        raise RuntimeError(f"heat budget does not close (relative error {worst:.1e})")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(dump_summary(summarise(results, args.revision)))
    if args.figures:
        from cooling import figures

        figures.write_all(results, args.figures, args.revision)
    base = [r for r in results if r.case.label == "base"]
    for r in base:
        d = r.doubling_hours(r.ratio_flow)
        print(f"{r.case.day:6s} H = {r.case.thickness:4.2f} m: flow-weighted mu x{r.ratio_flow[-1]:.2f} at 12 h, "
              f"doubles after {'-' if d is None else f'{d:.1f} h'}; min {r.min_temperature - 273.15:.1f} degC")
    print(f"{len(results)} runs, largest heat budget error {worst:.1e}")


if __name__ == "__main__":
    main()
