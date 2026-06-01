#!/usr/bin/env python3
"""
lrc_n14_n18_tournament_pingpong_s492.py

codex-2026-06-01 S492

Alternate between the n=14 and n=18 Lonely Runner frontiers and test whether
the current tournament-analysis lifts see a genuine cyclic pressure core.

This script deliberately keeps the "ping-pong" structure visible:

* compare the initial, quotient-ladder, gate-ladder, and single-gate-repair
  rows for n=14 and n=18;
* at exact endpoint/gap times, compute several pairwise tournament shadows;
* look for a strict SCC or directed 3-cycle in any shadow before proposing a
  proof route.

The forced-random inspiration for this pass was the k-nearest-neighbor graph:
the old S470 pressure graph used nearest-neighbor deletion relief.  Here we
also use two-neighbor deletion relief and threshold-deficit relief.  If those
still peel, the pressure-peeling hypothesis becomes more credible.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S470 = SourceFileLoader(
    "lrc_pairwise_tournament_s470",
    str(ROOT / "04-computation" / "lrc_pairwise_tournament_s470.py"),
).load_module()


ONE = Fraction(1, 1)
SCAN_LIMIT = 640


@dataclass(frozen=True)
class Example:
    label: str
    n: int
    speeds: tuple[int, ...]
    role: str


@dataclass(frozen=True)
class GraphMetrics:
    strict_arcs: int
    ties: int
    strict_triangles: int
    largest_scc: int
    source_count: int
    sink_count: int
    outdegree_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class SafeGapMetrics:
    safe_gaps: int
    lonely_vertices: int
    max_safe_run: int
    origin_left: Fraction
    origin_right: Fraction
    min_gap: Fraction
    max_gap: Fraction


@dataclass(frozen=True)
class TimeDigest:
    tag: str
    t: Fraction
    origin_d1: Fraction
    origin_d2: Fraction
    best_speed: int
    best_d1: Fraction
    best_d2: Fraction
    half_triangles: int
    half_ties: int
    k1: GraphMetrics
    k2: GraphMetrics
    deficit: GraphMetrics
    safe: SafeGapMetrics


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def ratio(value: Fraction, unit: Fraction) -> str:
    return fmt(value / unit)


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or not primitive(speeds):
        raise ValueError(f"bad ladder n={n} scale={scale} skip={skip}")
    return speeds


def largest_proper_divisor(n: int) -> int:
    for d in range(n // 2, 1, -1):
        if n % d == 0:
            return d
    raise ValueError(f"{n} has no proper divisor")


@lru_cache(maxsize=None)
def report_cached(speeds: tuple[int, ...]):
    return S356.report("cached", list(speeds))


@lru_cache(maxsize=None)
def summary_cached(speeds: tuple[int, ...]):
    return S360.summarize(list(speeds))


def best_ladder(n: int, scale: int) -> tuple[int, tuple[int, ...]]:
    candidates = []
    threshold = Fraction(1, n)
    for skip in range(1, n):
        speeds = ladder(n, scale, skip)
        report = report_cached(speeds)
        candidates.append(
            (
                report.max_gap / threshold,
                report.boundary_witness_count,
                skip,
                speeds,
            )
        )
    _, _, skip, speeds = min(candidates)
    return skip, speeds


def graph_metrics(vertex_count: int, arcs: set[tuple[int, int]], ties: int) -> GraphMetrics:
    out = [0] * vertex_count
    indeg = [0] * vertex_count
    for a, b in arcs:
        out[a] += 1
        indeg[b] += 1
    sizes = S470.scc_sizes(vertex_count, arcs)
    return GraphMetrics(
        strict_arcs=len(arcs),
        ties=ties,
        strict_triangles=S470.triangle_count(vertex_count, arcs),
        largest_scc=sizes[0] if sizes else 0,
        source_count=sum(1 for v in range(vertex_count) if indeg[v] == 0 and out[v] > 0),
        sink_count=sum(1 for v in range(vertex_count) if out[v] == 0 and indeg[v] > 0),
        outdegree_hist=tuple(sorted(Counter(out).items())),
    )


def nearest_sum(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    k: int,
    removed: int | None = None,
) -> Fraction:
    values = sorted(
        dist[i][j] for j in range(len(dist)) if j != i and j != removed
    )
    return sum(values[:k], Fraction(0))


def k_relief_metrics(dist: tuple[tuple[Fraction, ...], ...], k: int) -> GraphMetrics:
    vertex_count = len(dist)
    base = [nearest_sum(dist, i, k) for i in range(vertex_count)]
    arcs: set[tuple[int, int]] = set()
    ties = 0
    for i, j in combinations(range(vertex_count), 2):
        relief_i_from_j = nearest_sum(dist, i, k, removed=j) - base[i]
        relief_j_from_i = nearest_sum(dist, j, k, removed=i) - base[j]
        if relief_i_from_j > relief_j_from_i:
            arcs.add((j, i))
        elif relief_j_from_i > relief_i_from_j:
            arcs.add((i, j))
        else:
            ties += 1
    return graph_metrics(vertex_count, arcs, ties)


def first_distances(
    dist: tuple[tuple[Fraction, ...], ...],
    i: int,
    k: int,
    removed: int | None = None,
) -> tuple[Fraction, ...]:
    values = sorted(
        dist[i][j] for j in range(len(dist)) if j != i and j != removed
    )
    return tuple(values[:k])


def deficit_score(values: tuple[Fraction, ...], threshold: Fraction) -> Fraction:
    return sum((threshold - value for value in values if value < threshold), Fraction(0))


def deficit_relief_metrics(
    dist: tuple[tuple[Fraction, ...], ...], threshold: Fraction
) -> GraphMetrics:
    vertex_count = len(dist)
    base = [
        deficit_score(first_distances(dist, i, 2), threshold)
        for i in range(vertex_count)
    ]
    arcs: set[tuple[int, int]] = set()
    ties = 0
    for i, j in combinations(range(vertex_count), 2):
        after_i = deficit_score(first_distances(dist, i, 2, removed=j), threshold)
        after_j = deficit_score(first_distances(dist, j, 2, removed=i), threshold)
        relief_i_from_j = base[i] - after_i
        relief_j_from_i = base[j] - after_j
        if relief_i_from_j > relief_j_from_i:
            arcs.add((j, i))
        elif relief_j_from_i > relief_i_from_j:
            arcs.add((i, j))
        else:
            ties += 1
    return graph_metrics(vertex_count, arcs, ties)


def longest_cyclic_true_run(mask: tuple[bool, ...]) -> int:
    if not mask:
        return 0
    if all(mask):
        return len(mask)
    doubled = mask + mask
    best = 0
    cur = 0
    for value in doubled:
        if value:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return min(best, len(mask))


def safe_gap_metrics(
    runner_speeds: tuple[int, ...], t: Fraction, threshold: Fraction
) -> SafeGapMetrics:
    pos = S470.positions(runner_speeds, t)
    ordered = sorted((position, idx) for idx, position in enumerate(pos))
    n = len(ordered)
    gaps = []
    origin_pos_index = 0
    for idx, (_, runner_index) in enumerate(ordered):
        if runner_index == 0:
            origin_pos_index = idx
        nxt = ordered[(idx + 1) % n][0]
        cur = ordered[idx][0]
        gaps.append((nxt - cur) % ONE)
    safe = tuple(gap >= threshold for gap in gaps)
    lonely_vertices = sum(
        1 for idx in range(n) if safe[idx - 1] and safe[idx]
    )
    return SafeGapMetrics(
        safe_gaps=sum(1 for value in safe if value),
        lonely_vertices=lonely_vertices,
        max_safe_run=longest_cyclic_true_run(safe),
        origin_left=gaps[origin_pos_index - 1],
        origin_right=gaps[origin_pos_index],
        min_gap=min(gaps),
        max_gap=max(gaps),
    )


def digest_time(tag: str, speeds: tuple[int, ...], t: Fraction) -> TimeDigest:
    runner_speeds = (0,) + speeds
    threshold = Fraction(1, len(runner_speeds))
    pos = S470.positions(runner_speeds, t)
    dist = S470.distance_matrix(pos)
    profs = S470.profiles(runner_speeds, t)
    best = max(profs, key=lambda row: (row.d1, row.d2, -row.speed))
    half = S470.half_circle_metrics(pos)
    return TimeDigest(
        tag=tag,
        t=t,
        origin_d1=profs[0].d1,
        origin_d2=profs[0].d2,
        best_speed=best.speed,
        best_d1=best.d1,
        best_d2=best.d2,
        half_triangles=half.strict_triangles,
        half_ties=half.collision_ties + half.antipodal_ties,
        k1=k_relief_metrics(dist, 1),
        k2=k_relief_metrics(dist, 2),
        deficit=deficit_relief_metrics(dist, threshold),
        safe=safe_gap_metrics(runner_speeds, t, threshold),
    )


def candidate_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    return S470.candidate_times(speeds)


def scanned_candidate_times(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, ...], int]:
    candidates = candidate_times(speeds)
    if len(candidates) <= SCAN_LIMIT:
        return candidates, len(candidates)

    n = len(speeds) + 1
    report = report_cached(speeds)
    essential = {
        Fraction(0),
        Fraction(1, n),
        Fraction(1, 2 * n),
        Fraction(1, 2),
    }
    if report.witness is not None:
        essential.add(report.witness)
    if report.boundary_witness is not None:
        essential.add(report.boundary_witness)

    budget = max(1, SCAN_LIMIT - len(essential))
    stride = max(1, len(candidates) // budget)
    chosen = set(essential)
    chosen.update(candidates[::stride])
    chosen.add(candidates[-1])

    if len(chosen) > SCAN_LIMIT:
        nonessential = sorted(chosen - essential)
        keep_count = max(0, SCAN_LIMIT - len(essential))
        if keep_count:
            stride2 = max(1, len(nonessential) // keep_count)
            chosen = set(essential)
            chosen.update(nonessential[::stride2][:keep_count])
        else:
            chosen = set(sorted(essential)[:SCAN_LIMIT])
    return tuple(sorted(chosen)), len(candidates)


def unique_by_time(rows: list[TimeDigest]) -> tuple[TimeDigest, ...]:
    seen: set[Fraction] = set()
    out: list[TimeDigest] = []
    for row in rows:
        if row.t in seen:
            continue
        seen.add(row.t)
        out.append(row)
    return tuple(out)


def retag(row: TimeDigest, tag: str) -> TimeDigest:
    return TimeDigest(
        tag=tag,
        t=row.t,
        origin_d1=row.origin_d1,
        origin_d2=row.origin_d2,
        best_speed=row.best_speed,
        best_d1=row.best_d1,
        best_d2=row.best_d2,
        half_triangles=row.half_triangles,
        half_ties=row.half_ties,
        k1=row.k1,
        k2=row.k2,
        deficit=row.deficit,
        safe=row.safe,
    )


def selected_digests_clean(example: Example) -> tuple[tuple[TimeDigest, ...], int, int]:
    speeds = example.speeds
    report = report_cached(speeds)
    candidates, total_candidates = scanned_candidate_times(speeds)
    all_rows = [digest_time("scan", speeds, t) for t in candidates]
    selected: list[TimeDigest] = []
    if report.witness is not None:
        selected.append(digest_time("origin gap mid", speeds, report.witness))
    if report.boundary_witness is not None:
        selected.append(digest_time("origin boundary", speeds, report.boundary_witness))
    selected.extend(
        [
            retag(
                max(
                    all_rows,
                    key=lambda row: (
                        min(row.safe.origin_left, row.safe.origin_right),
                        row.origin_d1,
                        -row.t,
                    ),
                ),
                "best origin bracket",
            ),
            retag(
                max(
                    all_rows,
                    key=lambda row: (
                        row.k1.largest_scc,
                        row.k1.strict_triangles,
                        row.k1.strict_arcs,
                        -row.t,
                    ),
                ),
                "max k1 pressure",
            ),
            retag(
                max(
                    all_rows,
                    key=lambda row: (
                        row.k2.largest_scc,
                        row.k2.strict_triangles,
                        row.k2.strict_arcs,
                        -row.t,
                    ),
                ),
                "max k2 pressure",
            ),
            retag(
                max(
                    all_rows,
                    key=lambda row: (
                        row.deficit.largest_scc,
                        row.deficit.strict_triangles,
                        row.deficit.strict_arcs,
                        -row.t,
                    ),
                ),
                "max deficit pressure",
            ),
            retag(
                max(
                    all_rows,
                    key=lambda row: (
                        row.safe.lonely_vertices,
                        row.safe.safe_gaps,
                        row.safe.max_safe_run,
                        -row.t,
                    ),
                ),
                "max mobile lonely",
            ),
        ]
    )
    return unique_by_time(selected), len(candidates), total_candidates


def make_examples() -> tuple[Example, ...]:
    n14_lpd_skip, n14_lpd = best_ladder(14, largest_proper_divisor(14))
    n18_lpd_skip, n18_lpd = best_ladder(18, largest_proper_divisor(18))
    n14_gate_skip, n14_gate = best_ladder(14, 14)
    n18_gate_skip, n18_gate = best_ladder(18, 18)
    n14_repair = tuple(sorted((set(initial(14)) - {6}) | {14 * 16}))
    n18_repair = tuple(sorted((set(initial(18)) - {8}) | {18 * 18}))
    return (
        Example("initial n=14", 14, initial(14), "boundary unit skeleton"),
        Example("initial n=18", 18, initial(18), "boundary unit skeleton"),
        Example(
            "n=14 lpd ladder",
            14,
            n14_lpd,
            f"scale 7 skip {n14_lpd_skip}; classic seven ladder",
        ),
        Example(
            "n=18 lpd ladder",
            18,
            n18_lpd,
            f"scale 9 skip {n18_lpd_skip}; square-odd-core ladder",
        ),
        Example(
            "n=14 gate ladder",
            14,
            n14_gate,
            f"scale 14 skip {n14_gate_skip}; S380-style gate depth",
        ),
        Example(
            "n=18 gate ladder",
            18,
            n18_gate,
            f"scale 18 skip {n18_gate_skip}; doubled 9-core gate depth",
        ),
        Example(
            "n=14 single-gate repair",
            14,
            n14_repair,
            "HYP-1828 low-exposure repair: replace 6 by 14*16",
        ),
        Example(
            "n=18 single-gate repair",
            18,
            n18_repair,
            "n=18 analogue: replace lpd skip 8 by 18*18",
        ),
    )


def print_header() -> None:
    print("codex-2026-06-01-S492")
    print("n=14 / n=18 LRC tournament-analysis ping-pong")
    print("=" * 110)
    print()
    print("Forced-random import used in this pass:")
    print("  k-nearest-neighbor graphs suggest replacing one-neighbor deletion relief")
    print("  by k=2 relief and by threshold-deficit relief.  Those are still finite")
    print("  incomplete tournaments on the runners at exact LRC endpoint/gap times.")
    print()


def print_example_summary(example: Example, scanned_count: int, total_count: int) -> None:
    report = report_cached(example.speeds)
    summary = summary_cached(example.speeds)
    threshold = Fraction(1, example.n)
    print(f"[{example.label}]")
    print(f"  role={example.role}")
    print(f"  speeds={example.speeds}")
    print(
        "  "
        f"class={summary.classification} gap/th={ratio(report.max_gap, threshold)} "
        f"forbidden={fmt(report.forbidden_length)} endpoints={summary.unique_endpoint_count} "
        f"unprotected={summary.unprotected_count} first={fmt(summary.first_unprotected)} "
        f"sampled_times={scanned_count}/{total_count}"
    )


def graph_text(row: GraphMetrics) -> str:
    return (
        f"{row.strict_arcs}/{row.ties} "
        f"tri={row.strict_triangles} scc={row.largest_scc} "
        f"src/sink={row.source_count}/{row.sink_count}"
    )


def print_time_digest(example: Example, row: TimeDigest) -> None:
    threshold = Fraction(1, example.n)
    print(
        f"  {row.tag:<21} t={fmt(row.t):>14} "
        f"origin={ratio(row.origin_d1, threshold):>7},{ratio(row.origin_d2, threshold):<7} "
        f"best={row.best_speed:>4}:{ratio(row.best_d1, threshold):>7},{ratio(row.best_d2, threshold):<7} "
        f"safe={row.safe.safe_gaps:>2} lonely={row.safe.lonely_vertices:>2} "
        f"br={ratio(min(row.safe.origin_left, row.safe.origin_right), threshold):>7} "
        f"half=tie{row.half_ties:>2}/tri{row.half_triangles:<4}"
    )
    print(
        f"    k1       {graph_text(row.k1)}"
    )
    print(
        f"    k2       {graph_text(row.k2)}"
    )
    print(
        f"    deficit  {graph_text(row.deficit)}"
    )


def print_synthesis(examples: tuple[Example, ...], rows_by_label: dict[str, tuple[TimeDigest, ...]]) -> None:
    print("CROSS-CASE SYNTHESIS")
    print("=" * 110)
    totals: Counter[tuple[int, str]] = Counter()
    cyclic_rows: list[tuple[str, str, str, GraphMetrics]] = []
    for example in examples:
        for row in rows_by_label[example.label]:
            for name, metrics in (
                ("k1", row.k1),
                ("k2", row.k2),
                ("deficit", row.deficit),
            ):
                totals[(example.n, name)] += 1
                if metrics.largest_scc > 1 or metrics.strict_triangles > 0:
                    cyclic_rows.append((example.label, row.tag, name, metrics))

    for n in (14, 18):
        print(f"n={n} selected tournament rows:")
        for name in ("k1", "k2", "deficit"):
            count = totals[(n, name)]
            bad = sum(
                1
                for label, _, metric_name, metrics in cyclic_rows
                if metric_name == name
                and next(example.n for example in examples if example.label == label) == n
                and (metrics.largest_scc > 1 or metrics.strict_triangles > 0)
            )
            print(f"  {name:<7} cyclic-or-SCC rows {bad}/{count}")
    print()

    if cyclic_rows:
        print("Rows that resisted immediate pressure peeling:")
        for label, tag, name, metrics in cyclic_rows[:12]:
            print(
                f"  {label} / {tag} / {name}: "
                f"arcs={metrics.strict_arcs} ties={metrics.ties} "
                f"tri={metrics.strict_triangles} scc={metrics.largest_scc}"
            )
    else:
        print("No selected row produced a strict pressure SCC or directed 3-cycle.")
    print()

    print("Working conclusion")
    print("-" * 110)
    print("1. The n=14 seven/gate ladders remain pressure-peelable even after the")
    print("   k=2 and threshold-deficit tournament lifts.")
    print("2. The n=18 lpd and gate ladders behave the same way in these exact samples:")
    print("   they increase endpoint debt and pairwise crowding, but not a mobile")
    print("   pressure core.")
    print("3. The best next proof move is therefore not to hunt arbitrary speed-level")
    print("   cycles.  It is to combine pressure leaves with labelled endpoint")
    print("   handoff leaves.  A counterexample must survive both peelings.")
    print("4. If future search finds a strict k=2 or deficit SCC for n=18 before n=14,")
    print("   that SCC is the first place to try THM-365-style arithmetic labels.")


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    print_header()
    examples = make_examples()
    rows_by_label: dict[str, tuple[TimeDigest, ...]] = {}
    for example in examples:
        rows, scanned_count, total_count = selected_digests_clean(example)
        rows_by_label[example.label] = rows
        print_example_summary(example, scanned_count, total_count)
        for row in rows:
            print_time_digest(example, row)
        print()
    print_synthesis(examples, rows_by_label)


if __name__ == "__main__":
    main()
