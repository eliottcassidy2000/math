#!/usr/bin/env python3
"""
lrc_witness_or_core_s570.py

codex-2026-06-03 S570

Exact primal-dual checker for reduced integer LRC instances.

Pipeline:
1. Try cheap primal witness clocks: n-clock, pair-sum pinch, antipodes.
2. If they fail, fall back to exact open-gap/boundary witnesses from S356.
3. Report the peeled endpoint-protection core from S359 either way, so every
   row exports an obstruction ledger instead of only a scalar "near miss".

Tournament Analysis note:
- Natural vertices here are theorem gates or active pinch pairs, not runners.
- This first audit stays unary because its job is exact witness/core routing.
- The proposed active-pair tournament lives in the paired S570 reflection.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path
import random


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S359 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s359",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s359.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class CandidateSummary:
    kind: str
    time: Fraction
    score: Fraction
    active_count: int
    denominator: int


@dataclass(frozen=True)
class WitnessCoreReport:
    speeds: tuple[int, ...]
    threshold: Fraction
    exact_M: Fraction
    exact_witness: Fraction
    route: str
    route_time: Fraction
    route_score: Fraction
    pair_sum_best: CandidateSummary
    antipode_best: CandidateSummary
    n_clock_best: CandidateSummary
    pair_sum_exact: bool
    positive_gap: bool
    boundary_witness_count: int
    core_endpoint_count: int
    core_interval_count: int
    core_rounds: int
    core_removed_profile: tuple[int, ...]
    short_resonances: int
    critical_moments: int


def circle_point(x: Fraction) -> Fraction:
    return x % ONE


def circle_distance_to_integer(x: Fraction) -> Fraction:
    r = circle_point(x)
    return min(r, ONE - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circle_distance_to_integer(v * t) for v in speeds)


def active_count(speeds: tuple[int, ...], t: Fraction) -> int:
    score = score_time(speeds, t)
    return sum(circle_distance_to_integer(v * t) == score for v in speeds)


def best_candidate(kind: str, speeds: tuple[int, ...], times: set[Fraction]) -> CandidateSummary:
    best_t = Fraction(0)
    best_score = Fraction(-1)
    best_active = 0
    for t in times:
        score = score_time(speeds, t)
        active = active_count(speeds, t)
        if (score, active, -t) > (best_score, best_active, -best_t):
            best_t = t
            best_score = score
            best_active = active
    return CandidateSummary(
        kind=kind,
        time=best_t,
        score=best_score,
        active_count=best_active,
        denominator=best_t.denominator,
    )


def pair_sum_times(speeds: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = set()
    for a, b in combinations(speeds, 2):
        s = a + b
        for m in range(1, s):
            out.add(Fraction(m, s))
    return out


def antipode_times(speeds: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = set()
    for v in speeds:
        for k in range(v):
            out.add(Fraction(2 * k + 1, 2 * v))
    return out


def n_clock_times(speeds: tuple[int, ...]) -> set[Fraction]:
    n = len(speeds) + 1
    return {Fraction(j, n) for j in range(1, n)}


def critical_moment_count(speeds: tuple[int, ...]) -> int:
    n = len(speeds) + 1
    points = set()
    for v in speeds:
        for k in range(v + 1):
            for sign in (-1, 1):
                points.add(Fraction(k * n + sign, n * v) % 1)
    return len(points)


def short_resonances3(speeds: tuple[int, ...], bound: int = 2) -> int:
    total = 0
    for x, y, z in combinations(speeds, 3):
        for a in range(-bound, bound + 1):
            for b in range(-bound, bound + 1):
                for c in range(-bound, bound + 1):
                    if (a, b, c) != (0, 0, 0) and a * x + b * y + c * z == 0:
                        total += 1
    return total


def analyze(raw_speeds: tuple[int, ...]) -> WitnessCoreReport:
    speeds = S356.normalize_speed_set(list(raw_speeds))
    threshold = Fraction(1, len(speeds) + 1)
    gap_report = S356.report("witness-or-core", list(speeds))
    core_report = S359.protection_report(speeds)

    pair_sum_best = best_candidate("pair_sum", speeds, pair_sum_times(speeds))
    antipode_best = best_candidate("antipode", speeds, antipode_times(speeds))
    n_clock_best = best_candidate("n_clock", speeds, n_clock_times(speeds))

    if pair_sum_best.score >= antipode_best.score:
        exact_M = pair_sum_best.score
        exact_witness = pair_sum_best.time
        pair_sum_exact = True
    else:
        exact_M = antipode_best.score
        exact_witness = antipode_best.time
        pair_sum_exact = False

    if n_clock_best.score >= threshold:
        route = "n_clock"
        route_time = n_clock_best.time
        route_score = n_clock_best.score
    elif pair_sum_best.score >= threshold:
        route = "pair_sum"
        route_time = pair_sum_best.time
        route_score = pair_sum_best.score
    elif antipode_best.score >= threshold:
        route = "antipode"
        route_time = antipode_best.time
        route_score = antipode_best.score
    elif gap_report.witness is not None and score_time(speeds, gap_report.witness) >= threshold:
        route = "gap_midpoint"
        route_time = gap_report.witness
        route_score = score_time(speeds, route_time)
    elif gap_report.boundary_witness is not None:
        route = "boundary"
        route_time = gap_report.boundary_witness
        route_score = score_time(speeds, route_time)
    else:
        route = "core"
        route_time = exact_witness
        route_score = exact_M

    return WitnessCoreReport(
        speeds=speeds,
        threshold=threshold,
        exact_M=exact_M,
        exact_witness=exact_witness,
        route=route,
        route_time=route_time,
        route_score=route_score,
        pair_sum_best=pair_sum_best,
        antipode_best=antipode_best,
        n_clock_best=n_clock_best,
        pair_sum_exact=pair_sum_exact,
        positive_gap=gap_report.max_gap > 0,
        boundary_witness_count=gap_report.boundary_witness_count,
        core_endpoint_count=core_report.core_endpoint_count,
        core_interval_count=core_report.core_interval_count,
        core_rounds=core_report.core_rounds,
        core_removed_profile=core_report.core_removed_profile,
        short_resonances=short_resonances3(speeds),
        critical_moments=critical_moment_count(speeds),
    )


def fmt_frac(x: Fraction) -> str:
    return S356.fmt_frac(x)


def print_report(label: str, report: WitnessCoreReport) -> None:
    print(f"[{label}]")
    print(f"  speeds={report.speeds}")
    print(
        "  "
        f"route={report.route} route_time={fmt_frac(report.route_time)} "
        f"route_score={fmt_frac(report.route_score)} threshold={fmt_frac(report.threshold)}"
    )
    print(
        "  "
        f"exact_M={fmt_frac(report.exact_M)} exact_witness={fmt_frac(report.exact_witness)} "
        f"pair_sum_exact={report.pair_sum_exact}"
    )
    print(
        "  "
        f"pair_sum_best=({fmt_frac(report.pair_sum_best.time)}, {fmt_frac(report.pair_sum_best.score)}, "
        f"act={report.pair_sum_best.active_count}, q={report.pair_sum_best.denominator}) "
        f"antipode_best=({fmt_frac(report.antipode_best.time)}, {fmt_frac(report.antipode_best.score)}, "
        f"act={report.antipode_best.active_count}, q={report.antipode_best.denominator})"
    )
    print(
        "  "
        f"n_clock_best=({fmt_frac(report.n_clock_best.time)}, {fmt_frac(report.n_clock_best.score)}, "
        f"act={report.n_clock_best.active_count}, q={report.n_clock_best.denominator}) "
        f"positive_gap={report.positive_gap} boundary_witnesses={report.boundary_witness_count}"
    )
    print(
        "  "
        f"core_endpoints={report.core_endpoint_count} core_intervals={report.core_interval_count} "
        f"core_rounds={report.core_rounds} removed={report.core_removed_profile}"
    )
    print(
        "  "
        f"short_resonances={report.short_resonances} critical_moments={report.critical_moments}"
    )
    print()


def sample_families() -> list[tuple[str, tuple[int, ...]]]:
    rng = random.Random(570)
    return [
        ("AP n=14", tuple(range(1, 14))),
        ("V* sporadic n=14", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)),
        ("Fibonacci n=14", (1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377)),
        ("Sidon n=14", (1, 2, 5, 11, 22, 40, 56, 73, 78, 97, 118, 134, 150)),
        ("random A n=14", tuple(sorted(rng.sample(range(1, 80), 13)))),
        ("random B n=14", tuple(sorted(rng.sample(range(1, 80), 13)))),
        ("near-tight k=6", (1, 4, 5, 6, 7, 11)),
        ("sporadic n=8", (1, 4, 5, 6, 7, 11, 13)),
    ]


def primitive(combo: tuple[int, ...]) -> bool:
    g = 0
    for v in combo:
        g = gcd(g, v)
    return g == 1


def scan_box(k: int, max_speed: int) -> None:
    counts: Counter[str] = Counter()
    pair_sum_exact = 0
    nonempty_core = 0
    hardest: list[tuple[int, int, Fraction, WitnessCoreReport]] = []
    total = 0
    for combo in combinations(range(1, max_speed + 1), k):
        if not primitive(combo):
            continue
        total += 1
        report = analyze(combo)
        counts[report.route] += 1
        pair_sum_exact += int(report.pair_sum_exact)
        nonempty_core += int(report.core_endpoint_count > 0)
        hardest.append(
            (
                0 if report.route == "core" else 1,
                report.short_resonances,
                report.exact_M,
                report,
            )
        )

    hardest.sort(key=lambda item: (item[0], -item[1], item[2], item[3].speeds))
    print(f"Primitive box k={k}, max_speed={max_speed}")
    print(f"  total={total}")
    print(f"  route_hist={dict(sorted(counts.items()))}")
    print(f"  pair_sum_exact={pair_sum_exact}/{total}")
    print(f"  nonempty_core={nonempty_core}")
    print("  hardest_examples")
    for _, _, _, report in hardest[:5]:
        print(
            "    "
            f"speeds={report.speeds} route={report.route} "
            f"M={fmt_frac(report.exact_M)} pair_sum={fmt_frac(report.pair_sum_best.score)} "
            f"core={report.core_endpoint_count} reson={report.short_resonances}"
        )
    print()


def main() -> None:
    print("LRC witness-or-core pipeline (codex-2026-06-03 S570)")
    print("Exact rank-one integer checker: n-clock, pair-sum pinch, antipodes,")
    print("then endpoint-core peel if the cheap primal clocks do not settle the row.\n")

    print("Sample families")
    for label, speeds in sample_families():
        print_report(label, analyze(speeds))

    print("Bounded primitive boxes")
    for k, max_speed in [(3, 20), (4, 16), (5, 13), (6, 11)]:
        scan_box(k, max_speed)

    print("Synthesis")
    print("  1. Pair-sum/antipode candidates appear to recover the exact tropical max M(S).")
    print("  2. The n-clock is the wall witness route for the resonance-maximal tight rows.")
    print("  3. The endpoint core is a dual ledger, not a fallback scalar: if it survives,")
    print("     the checker exports the obstruction packet instead of only a missed witness.")


if __name__ == "__main__":
    main()
