#!/usr/bin/env python3
"""
lrc_n14_n18_tournament_feedback_s490.py

codex-2026-06-01 S490

Alternate between the n=14 and n=18 Lonely Runner frontiers through the
repo's Tournament Analysis lens.  When one side gets stuck, inject a "noise
card" from adjacent graph/geometry literature and test whether it suggests a
new computable feature on the other side.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S452 = SourceFileLoader(
    "lrc_runner_distance_tournament_s452",
    str(ROOT / "04-computation" / "lrc_runner_distance_tournament_s452.py"),
).load_module()
S470 = SourceFileLoader(
    "lrc_pairwise_tournament_s470",
    str(ROOT / "04-computation" / "lrc_pairwise_tournament_s470.py"),
).load_module()


@dataclass(frozen=True)
class LadderProbe:
    n: int
    scale: int
    skip: int
    speeds: tuple[int, ...]
    gap_ratio: Fraction
    forbidden_gap: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    pressure_density: Fraction | None
    classification: str


@dataclass(frozen=True)
class TournamentProbe:
    label: str
    tag: str
    time: Fraction
    origin_pair: tuple[Fraction, Fraction]
    best_speed: int
    best_pair: tuple[Fraction, Fraction]
    pressure_arcs: int
    pressure_ties: int
    pressure_triangles: int
    pressure_largest_scc: int
    pressure_sources: int
    pressure_sinks: int
    half_ties: int
    half_score_width: int
    half_triangles: int


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def factorization(n: int) -> str:
    x = n
    out: list[str] = []
    p = 2
    while p * p <= x:
        if x % p == 0:
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            out.append(f"{p}^{e}" if e > 1 else str(p))
        p += 1 if p == 2 else 2
    if x > 1:
        out.append(str(x))
    return "*".join(out)


def divisors(n: int) -> tuple[int, ...]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return tuple(sorted(out))


def proper_scales(n: int) -> tuple[int, ...]:
    return tuple(d for d in divisors(n) if 1 < d <= n)


def phi(n: int) -> int:
    return sum(1 for a in range(1, n) if gcd(a, n) == 1)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def scale_speeds(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        raise ValueError((n, scale, skip, speeds))
    return speeds


def probe_ladder(n: int, scale: int, skip: int) -> LadderProbe:
    speeds = scale_speeds(n, scale, skip)
    report = S356.report(f"n={n},d={scale},skip={skip}", list(speeds))
    summary = S360.summarize(list(speeds))
    gap_ratio = report.max_gap / report.threshold
    pressure_density = None
    if gap_ratio:
        pressure_density = Fraction(summary.unprotected_count, 1) / gap_ratio / (n * n)
    return LadderProbe(
        n=n,
        scale=scale,
        skip=skip,
        speeds=speeds,
        gap_ratio=gap_ratio,
        forbidden_gap=1 - report.forbidden_length,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        pressure_density=pressure_density,
        classification=summary.classification,
    )


def best_ladder(n: int, scale: int) -> LadderProbe:
    scored: list[tuple[Fraction, int, Fraction, int]] = []
    for skip in range(1, n):
        speeds = scale_speeds(n, scale, skip)
        report = S356.report(f"n={n},d={scale},skip={skip}", list(speeds))
        scored.append(
            (
                report.max_gap / report.threshold,
                report.boundary_witness_count,
                1 - report.forbidden_length,
                skip,
            )
        )
    best_skip = min(scored)[3]
    return probe_ladder(n, scale, best_skip)


def selected_tournament_probes(label: str, speeds: tuple[int, ...]) -> tuple[TournamentProbe, ...]:
    if len(speeds) + 1 <= 14:
        rows, _candidate_count = S470.selected_rows(speeds)
    else:
        rows = bounded_selected_rows(speeds)
    threshold = Fraction(1, len(speeds) + 1)
    out: list[TournamentProbe] = []
    for tag, row in rows:
        out.append(
            TournamentProbe(
                label=label,
                tag=tag,
                time=row.t,
                origin_pair=(row.origin_d1 / threshold, row.origin_d2 / threshold),
                best_speed=row.best_speed,
                best_pair=(row.best_d1 / threshold, row.best_d2 / threshold),
                pressure_arcs=row.pressure.strict_arcs,
                pressure_ties=row.pressure.ties,
                pressure_triangles=row.pressure.strict_triangles,
                pressure_largest_scc=row.pressure.largest_scc,
                pressure_sources=row.pressure.source_count,
                pressure_sinks=row.pressure.sink_count,
                half_ties=row.half.collision_ties + row.half.antipodal_ties,
                half_score_width=row.half.score_width,
                half_triangles=row.half.strict_triangles,
            )
        )
    return tuple(out)


def bounded_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    n = len(speeds) + 1
    report = S356.report("bounded", list(speeds))
    raw = S356.forbidden_intervals(speeds)
    components = S356.merge_intervals(raw)
    gaps = S356.circular_gaps(components)
    endpoints = sorted({S356.circle_point(point) for interval in raw for point in interval})
    stride = max(1, len(endpoints) // 36)
    sampled_endpoints = endpoints[::stride][:40]
    candidates = {
        Fraction(0),
        Fraction(1, n),
        Fraction(1, 2 * n),
        Fraction(1, 2),
    }
    if report.witness is not None:
        candidates.add(report.witness)
    if report.boundary_witness is not None:
        candidates.add(report.boundary_witness)
    for gap in gaps:
        candidates.add(S356.midpoint_on_circle(gap))
        candidates.add(S356.circle_point(gap[0]))
        candidates.add(S356.circle_point(gap[1]))
    candidates.update(sampled_endpoints)
    return tuple(sorted(candidates))


def unique_tagged(rows: list[tuple[str, object]]) -> list[tuple[str, object]]:
    seen: set[Fraction] = set()
    out: list[tuple[str, object]] = []
    for tag, row in rows:
        t = row.t
        if t in seen:
            continue
        seen.add(t)
        out.append((tag, row))
    return out


def bounded_selected_rows(speeds: tuple[int, ...]) -> list[tuple[str, object]]:
    analyses = [S470.analyze_time(speeds, t) for t in bounded_times(speeds)]
    report = S356.report("bounded", list(speeds))
    rows: list[tuple[str, object]] = []
    if report.witness is not None:
        rows.append(("origin gap mid", S470.analyze_time(speeds, report.witness)))
    if report.boundary_witness is not None:
        rows.append(("origin boundary", S470.analyze_time(speeds, report.boundary_witness)))
    rows.extend(
        [
            ("best origin d1", max(analyses, key=lambda row: (row.origin_d1, row.origin_d2, -row.t))),
            ("best mobile d1", max(analyses, key=lambda row: (row.best_d1, row.best_d2, -row.t))),
            (
                "max pressure core",
                max(
                    analyses,
                    key=lambda row: (
                        row.pressure.largest_scc,
                        row.pressure.strict_triangles,
                        row.pressure.strict_arcs,
                        -row.t,
                    ),
                ),
            ),
            (
                "max pressure cycles",
                max(
                    analyses,
                    key=lambda row: (
                        row.pressure.strict_triangles,
                        row.pressure.largest_scc,
                        row.pressure.strict_arcs,
                        -row.t,
                    ),
                ),
            ),
            (
                "max semicircle cycles",
                max(
                    analyses,
                    key=lambda row: (
                        row.half.strict_triangles,
                        -row.half.score_width,
                        row.half.strict_arcs,
                        -row.t,
                    ),
                ),
            ),
        ]
    )
    return unique_tagged(rows)


def safe_gap_report(label: str, speeds: tuple[int, ...]) -> tuple[str, str, int, int, int]:
    report = S356.report(label, list(speeds))
    t = report.witness or report.boundary_witness or Fraction(0)
    snap = S452.snapshot(label, speeds, t)
    safe_count = snap.safe_gap_mask.count("1")
    adjacent_safe = len(snap.lonely_vertices)
    longest_unsafe = max((len(part) for part in snap.safe_gap_mask.split("1")), default=0)
    return (label, S356.fmt_frac(t), safe_count, adjacent_safe, longest_unsafe)


def independent_cycle_count(n: int) -> int:
    return S452.cycle_independence(n, 1)


def hamiltonian_lower_noise(n: int) -> int:
    # Busch's strong-tournament lower bound is at least 5^((n-1)/3).
    # Use a conservative integer floor for feature-card reporting only.
    return int(5 ** Fraction(n - 1, 3))


def print_noise_cards() -> None:
    print("NOISE CARDS USED BETWEEN STUCK PHASES")
    print("=" * 108)
    print(
        "1. circular-arc independent sets: no-lonely static masks are cycle "
        "independent sets; use edge-cover language for unsafe gaps."
    )
    print(
        "2. strong tournament Hamiltonian paths: if a pressure SCC appears, "
        "it should carry a large Hamiltonian-path obligation, not just one cycle."
    )
    print(
        "3. zonotope/covering radius algorithms: treat n=18 as a product-tree "
        "fundamental-domain problem, not a scalar gap search."
    )
    print()


def print_ladder_table(probes: tuple[LadderProbe, ...]) -> None:
    print("SCALE-LADDER ENDPOINT PRESSURE")
    print("=" * 108)
    print(
        f"{'n':>3} {'fac':>7} {'scale':>5} {'skip':>4} {'class':>12} "
        f"{'gap/th':>10} {'1-meas':>10} {'bdy':>5} {'unprot':>7} "
        f"{'p/n^2':>10} {'first exposed'}"
    )
    print("-" * 108)
    for row in probes:
        print(
            f"{row.n:>3} {factorization(row.n):>7} {row.scale:>5} {row.skip:>4} "
            f"{row.classification:>12} {fmt(row.gap_ratio):>10} "
            f"{fmt(row.forbidden_gap):>10} {row.boundary_witnesses:>5} "
            f"{row.unprotected:>7} {fmt_float(row.pressure_density):>10} "
            f"{fmt(row.first_unprotected)}"
        )
    print()


def print_tournament_table(rows: tuple[TournamentProbe, ...]) -> None:
    print("TOURNAMENT PRESSURE SNAPSHOTS")
    print("=" * 108)
    print(
        f"{'case':<18} {'tag':<20} {'t':>12} {'origin':>13} {'best':>13} "
        f"{'p a/t':>8} {'p tri':>5} {'p scc':>5} {'src/sink':>8} "
        f"{'half ties':>9} {'hwidth':>6} {'htri':>5}"
    )
    print("-" * 108)
    for row in rows:
        origin = f"{fmt(row.origin_pair[0])},{fmt(row.origin_pair[1])}"
        best = f"{row.best_speed}:{fmt(row.best_pair[0])},{fmt(row.best_pair[1])}"
        print(
            f"{row.label:<18} {row.tag:<20} {fmt(row.time):>12} "
            f"{origin:>13} {best:>13} "
            f"{row.pressure_arcs:>3}/{row.pressure_ties:<3} "
            f"{row.pressure_triangles:>5} {row.pressure_largest_scc:>5} "
            f"{row.pressure_sources:>3}/{row.pressure_sinks:<3} "
            f"{row.half_ties:>9} {row.half_score_width:>6} {row.half_triangles:>5}"
        )
    print()


def print_safe_mask_table(rows: tuple[tuple[str, str, int, int, int], ...]) -> None:
    print("SAFE-GAP MASKS AT WITNESS TIMES")
    print("=" * 108)
    print(
        f"{'case':<22} {'t':>12} {'safe gaps':>9} "
        f"{'adjacent safe pairs':>19} {'longest unsafe run':>20}"
    )
    print("-" * 108)
    for label, t, safe_count, adjacent_safe, longest_unsafe in rows:
        print(
            f"{label:<22} {t:>12} {safe_count:>9} "
            f"{adjacent_safe:>19} {longest_unsafe:>20}"
        )
    print()


def print_alternation(probes: tuple[LadderProbe, ...], trows: tuple[TournamentProbe, ...]) -> None:
    n14 = [p for p in probes if p.n == 14]
    n18 = [p for p in probes if p.n == 18]
    n14_scc = max(row.pressure_largest_scc for row in trows if row.label.startswith("n14"))
    n18_scc = max(row.pressure_largest_scc for row in trows if row.label.startswith("n18"))
    best18 = min(n18, key=lambda row: row.gap_ratio)
    print("ALTERNATING WORK LEDGER")
    print("=" * 108)
    ledger = [
        (
            "n=14",
            "gate fan proof route",
            f"stuck at all-times certificate; sampled pressure SCC max={n14_scc}",
            "strong-tournament noise says any real SCC core should be expensive",
            "try pressure-leaf + endpoint-private-row induction first",
        ),
        (
            "n=18",
            "mixed torsion disproof pressure",
            f"best scale d={best18.scale} has gap/th={fmt(best18.gap_ratio)} and debt={best18.unprotected}",
            "circular-arc noise says unsafe gaps are an edge cover",
            "track safe-mask edge-cover transitions over the d=9 and d=18 ladders",
        ),
        (
            "n=14",
            "metric-switch tournament route",
            f"scale debts {', '.join(f'd={p.scale}:u={p.unprotected}' for p in n14)}",
            "zonotope noise says use fundamental-domain walls, not scalar extrema",
            "separate fixed gate walls from mobile chord-switch walls",
        ),
        (
            "n=18",
            "proof/disproof fork",
            f"scales tested={','.join(str(p.scale) for p in n18)}; phi/(n-1)={fmt(Fraction(phi(18),17))}",
            "n=14 leafiness suggests search for first non-leafy pressure SCC",
            "if no SCC appears, n=18 becomes a proof-lab rather than disproof-lab",
        ),
    ]
    print(f"{'side':<6} {'route':<28} {'stuck point':<44} next idea")
    print("-" * 108)
    for side, route, stuck, noise, action in ledger:
        print(f"{side:<6} {route:<28} {stuck:<44} {noise}")
        print(f"{'':<6} {'':<28} {'':<44} -> {action}")
    print()
    print("Current read:")
    print(
        "  n=14 remains the better proof target because its known hard rows look "
        "pressure-peelable.  n=18 is the better forced-random target because "
        "its 2*3^2 quotient structure creates more wall-crossing noise without "
        "yet showing a disproof-like pressure SCC."
    )


def main() -> None:
    print("LRC n=14 / n=18 Tournament Analysis feedback loop (codex-2026-06-01 S490)")
    print()
    print_noise_cards()

    n14_scales = (2, 7, 14)
    n18_scales = (2, 3, 6, 9, 18)
    probes = tuple(best_ladder(n, d) for n, scales in ((14, n14_scales), (18, n18_scales)) for d in scales)
    print_ladder_table(probes)

    named_cases = [
        ("n14 d=7", next(p for p in probes if p.n == 14 and p.scale == 7).speeds),
        ("n14 d=14", next(p for p in probes if p.n == 14 and p.scale == 14).speeds),
        ("n18 d=3", next(p for p in probes if p.n == 18 and p.scale == 3).speeds),
        ("n18 d=9", next(p for p in probes if p.n == 18 and p.scale == 9).speeds),
        ("n18 d=18", next(p for p in probes if p.n == 18 and p.scale == 18).speeds),
    ]
    trows = tuple(row for label, speeds in named_cases for row in selected_tournament_probes(label, speeds))
    print_tournament_table(trows)

    safe_rows = tuple(safe_gap_report(label, speeds) for label, speeds in named_cases)
    print_safe_mask_table(safe_rows)

    print("STATIC MASK SPACE")
    print("=" * 108)
    for n in (14, 18):
        print(
            f"n={n}: I(C_n,1)={independent_cycle_count(n)} "
            f"feasible no-lonely masks={independent_cycle_count(n)-1} "
            f"strong-tournament HP noise floor~{hamiltonian_lower_noise(n)}"
        )
    print()

    print_alternation(probes, trows)


if __name__ == "__main__":
    main()
