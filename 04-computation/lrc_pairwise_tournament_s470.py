#!/usr/bin/env python3
"""
lrc_pairwise_tournament_s470.py

codex-2026-05-31 S470

Explore how tournament structure appears when the Lonely Runner problem is
viewed through pairwise runner distances rather than only origin distance.

Two time-slice objects are built for the runners with speeds {0} union V:

1. semicircle orientation:
   i -> j when runner j lies in the clockwise open half-circle from i.
   Collisions and antipodal pairs are missing arcs, so this is an incomplete
   circular tournament at degenerate times.

2. blocker-pressure graph:
   For each unordered pair {i,j}, compare how much runner i's nearest-neighbor
   moat improves if j is deleted versus how much j's moat improves if i is
   deleted.  The stronger blocker points to the more blocked runner.  This
   uses nearest-neighbor distances with multiplicity; ties are left as
   missing arcs.

The point is proof technology: LRC endpoint protection is already a labelled
tournament cut-protection system.  This script asks whether the mobile
pairwise geometry has a complementary tournament-like protection core.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
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


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class RunnerProfile:
    index: int
    speed: int
    position: Fraction
    d1: Fraction
    d2: Fraction
    nearest: tuple[int, ...]
    second: tuple[int, ...]


@dataclass(frozen=True)
class HalfCircleMetrics:
    strict_arcs: int
    collision_ties: int
    antipodal_ties: int
    score_width: int
    score_hist: tuple[tuple[int, int], ...]
    strict_triangles: int


@dataclass(frozen=True)
class PressureMetrics:
    strict_arcs: int
    ties: int
    strict_triangles: int
    largest_scc: int
    source_count: int
    sink_count: int
    outdegree_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class TimeAnalysis:
    t: Fraction
    origin_d1: Fraction
    origin_d2: Fraction
    best_speed: int
    best_d1: Fraction
    best_d2: Fraction
    unique_nearest_vertices: int
    mutual_nearest_pairs: int
    half: HalfCircleMetrics
    pressure: PressureMetrics


@dataclass(frozen=True)
class Example:
    label: str
    speeds: tuple[int, ...]
    note: str


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def ratio(value: Fraction, unit: Fraction) -> str:
    return fmt(value / unit)


def circle(value: Fraction) -> Fraction:
    return value % ONE


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = circle(a - b)
    return min(delta, ONE - delta)


def clockwise_delta(a: Fraction, b: Fraction) -> Fraction:
    return circle(b - a)


def positions(runner_speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in runner_speeds)


def distance_matrix(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], ...]:
    n = len(pos)
    return tuple(
        tuple(Fraction(0) if i == j else circular_distance(pos[i], pos[j]) for j in range(n))
        for i in range(n)
    )


def profiles(
    runner_speeds: tuple[int, ...], t: Fraction
) -> tuple[RunnerProfile, ...]:
    pos = positions(runner_speeds, t)
    dist = distance_matrix(pos)
    out: list[RunnerProfile] = []
    for i, speed in enumerate(runner_speeds):
        grouped: dict[Fraction, list[int]] = defaultdict(list)
        ordered_neighbors: list[tuple[Fraction, int]] = []
        for j in range(len(runner_speeds)):
            if i != j:
                grouped[dist[i][j]].append(j)
                ordered_neighbors.append((dist[i][j], j))
        ordered_neighbors.sort()
        d1 = ordered_neighbors[0][0]
        d2 = ordered_neighbors[1][0]
        out.append(
            RunnerProfile(
                index=i,
                speed=speed,
                position=pos[i],
                d1=d1,
                d2=d2,
                nearest=tuple(sorted(grouped[d1])),
                second=tuple(sorted(j for distance, j in ordered_neighbors if distance == d2)),
            )
        )
    return tuple(out)


def triangle_count(vertex_count: int, arcs: set[tuple[int, int]]) -> int:
    count = 0
    for a, b, c in combinations(range(vertex_count), 3):
        present = {
            (a, b): (a, b) in arcs,
            (b, a): (b, a) in arcs,
            (a, c): (a, c) in arcs,
            (c, a): (c, a) in arcs,
            (b, c): (b, c) in arcs,
            (c, b): (c, b) in arcs,
        }
        if (
            (present[(a, b)] and present[(b, c)] and present[(c, a)])
            or (present[(a, c)] and present[(c, b)] and present[(b, a)])
        ):
            count += 1
    return count


def half_circle_metrics(pos: tuple[Fraction, ...]) -> HalfCircleMetrics:
    n = len(pos)
    arcs: set[tuple[int, int]] = set()
    collision_ties = 0
    antipodal_ties = 0
    score = [0] * n
    for i, j in combinations(range(n), 2):
        delta = clockwise_delta(pos[i], pos[j])
        if delta == 0:
            collision_ties += 1
        elif delta == HALF:
            antipodal_ties += 1
        elif delta < HALF:
            arcs.add((i, j))
            score[i] += 1
        else:
            arcs.add((j, i))
            score[j] += 1
    return HalfCircleMetrics(
        strict_arcs=len(arcs),
        collision_ties=collision_ties,
        antipodal_ties=antipodal_ties,
        score_width=max(score) - min(score),
        score_hist=tuple(sorted(Counter(score).items())),
        strict_triangles=triangle_count(n, arcs),
    )


def nearest_without(
    dist: tuple[tuple[Fraction, ...], ...], i: int, removed: int
) -> Fraction:
    return min(dist[i][j] for j in range(len(dist)) if j not in (i, removed))


def scc_sizes(vertex_count: int, arcs: set[tuple[int, int]]) -> tuple[int, ...]:
    adj = {v: [] for v in range(vertex_count)}
    radj = {v: [] for v in range(vertex_count)}
    for a, b in arcs:
        adj[a].append(b)
        radj[b].append(a)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(vertex_count):
        if v not in seen:
            dfs(v)

    sizes: list[int] = []
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for w in radj[v]:
                if w not in seen:
                    seen.add(w)
                    todo.append(w)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def pressure_metrics(
    profs: tuple[RunnerProfile, ...], dist: tuple[tuple[Fraction, ...], ...]
) -> PressureMetrics:
    n = len(profs)
    arcs: set[tuple[int, int]] = set()
    ties = 0
    out = [0] * n
    indeg = [0] * n
    for i, j in combinations(range(n), 2):
        relief_i_from_j = nearest_without(dist, i, j) - profs[i].d1
        relief_j_from_i = nearest_without(dist, j, i) - profs[j].d1
        if relief_i_from_j > relief_j_from_i:
            arc = (j, i)
        elif relief_j_from_i > relief_i_from_j:
            arc = (i, j)
        else:
            ties += 1
            continue
        arcs.add(arc)
        out[arc[0]] += 1
        indeg[arc[1]] += 1
    sizes = scc_sizes(n, arcs)
    return PressureMetrics(
        strict_arcs=len(arcs),
        ties=ties,
        strict_triangles=triangle_count(n, arcs),
        largest_scc=sizes[0] if sizes else 0,
        source_count=sum(1 for v in range(n) if indeg[v] == 0 and (indeg[v] + out[v]) > 0),
        sink_count=sum(1 for v in range(n) if out[v] == 0 and (indeg[v] + out[v]) > 0),
        outdegree_hist=tuple(sorted(Counter(out).items())),
    )


def mutual_nearest_count(profs: tuple[RunnerProfile, ...]) -> int:
    count = 0
    for i, j in combinations(range(len(profs)), 2):
        if j in profs[i].nearest and i in profs[j].nearest:
            count += 1
    return count


def analyze_time(speeds: tuple[int, ...], t: Fraction) -> TimeAnalysis:
    runner_speeds = (0,) + speeds
    pos = positions(runner_speeds, t)
    dist = distance_matrix(pos)
    profs = profiles(runner_speeds, t)
    best = max(profs, key=lambda row: (row.d1, row.d2, -row.speed))
    return TimeAnalysis(
        t=t,
        origin_d1=profs[0].d1,
        origin_d2=profs[0].d2,
        best_speed=best.speed,
        best_d1=best.d1,
        best_d2=best.d2,
        unique_nearest_vertices=sum(1 for row in profs if len(row.nearest) == 1),
        mutual_nearest_pairs=mutual_nearest_count(profs),
        half=half_circle_metrics(pos),
        pressure=pressure_metrics(profs, dist),
    )


def candidate_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    report = S356.report("candidate", list(speeds))
    components = S356.merge_intervals(S356.forbidden_intervals(speeds))
    gaps = S356.circular_gaps(components)
    n = len(speeds) + 1
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
    for endpoint in S360.endpoints(speeds):
        candidates.add(endpoint.value)
    for gap in gaps:
        candidates.add(S356.midpoint_on_circle(gap))
        candidates.add(circle(gap[0]))
        candidates.add(circle(gap[1]))
    return tuple(sorted(candidates))


def unique_rows(rows: list[tuple[str, TimeAnalysis]]) -> list[tuple[str, TimeAnalysis]]:
    seen: set[Fraction] = set()
    out: list[tuple[str, TimeAnalysis]] = []
    for tag, row in rows:
        if row.t in seen:
            continue
        seen.add(row.t)
        out.append((tag, row))
    return out


def selected_rows(speeds: tuple[int, ...]) -> tuple[list[tuple[str, TimeAnalysis]], int]:
    candidates = candidate_times(speeds)
    analyses = [analyze_time(speeds, t) for t in candidates]
    threshold = Fraction(1, len(speeds) + 1)
    report = S356.report("selected", list(speeds))
    rows: list[tuple[str, TimeAnalysis]] = []
    if report.witness is not None:
        rows.append(("origin gap mid", analyze_time(speeds, report.witness)))
    if report.boundary_witness is not None:
        rows.append(("origin boundary", analyze_time(speeds, report.boundary_witness)))
    rows.extend(
        [
            (
                "best origin d1",
                max(analyses, key=lambda row: (row.origin_d1, row.origin_d2, -row.t)),
            ),
            (
                "best mobile d1",
                max(analyses, key=lambda row: (row.best_d1, row.best_d2, -row.t)),
            ),
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
    return unique_rows(rows), len(candidates)


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def examples() -> tuple[Example, ...]:
    return (
        Example(
            "initial n=14",
            initial(14),
            "tight boundary row; origin sees the perfectly spaced unit skeleton",
        ),
        Example(
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
            "classic near-disproof ladder with gap/th=5/924",
        ),
        Example(
            "n14 S380 gate ladder",
            (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
            "gate-heavy composite-denominator anomaly with doubled endpoint debt",
        ),
        Example(
            "initial n=15",
            initial(15),
            "odd-denominator baseline beside the n=14 composite anomaly",
        ),
        Example(
            "initial n=16",
            initial(16),
            "pure dyadic baseline with a unique local gate-fan cover in older sessions",
        ),
    )


def print_methodology() -> None:
    print("PAIRWISE TOURNAMENT REFRAMES")
    print("=" * 96)
    print("1. Semicircle orientation: positions at time t define an incomplete circular")
    print("   tournament.  Missing arcs are exactly collisions or antipodal ties.")
    print("2. Two-neighbor pressure: i -> j means i is the more irreplaceable blocker")
    print("   of j's nearest-neighbor moat.  This turns first/second-nearest data into")
    print("   a labelled incomplete tournament, not just a scalar min-distance ledger.")
    print("3. A pressure SCC is the mobile analogue of a protected endpoint core: no")
    print("   runner in the component can be peeled without immediately changing")
    print("   someone else's local moat.")
    print()


def print_example(example: Example) -> None:
    speeds = S356.normalize_speed_set(list(example.speeds))
    report = S356.report(example.label, list(speeds))
    threshold = Fraction(1, len(speeds) + 1)
    rows, candidate_count = selected_rows(speeds)
    print(f"[{example.label}]")
    print(f"  note={example.note}")
    print(
        "  "
        f"n={len(speeds)+1} threshold={fmt(threshold)} "
        f"classification={S360.summarize(list(speeds)).classification} "
        f"gap/th={ratio(report.max_gap, threshold)} "
        f"boundary_witnesses={report.boundary_witness_count} "
        f"sampled_times={candidate_count}"
    )
    print(
        "  tag                  t              origin d1,d2   best speed d1,d2"
        "   nn uniq/mut   pressure arcs/ties tri scc src/sink   half ties width tri"
    )
    for tag, row in rows:
        print(
            f"  {tag:<20} {fmt(row.t):>14} "
            f"{ratio(row.origin_d1, threshold):>5},{ratio(row.origin_d2, threshold):<5} "
            f"{row.best_speed:>5} {ratio(row.best_d1, threshold):>5},{ratio(row.best_d2, threshold):<5} "
            f"{row.unique_nearest_vertices:>3}/{row.mutual_nearest_pairs:<3} "
            f"{row.pressure.strict_arcs:>4}/{row.pressure.ties:<4} "
            f"{row.pressure.strict_triangles:>3} {row.pressure.largest_scc:>3} "
            f"{row.pressure.source_count:>3}/{row.pressure.sink_count:<3} "
            f"{row.half.collision_ties:>2}+{row.half.antipodal_ties:<2} "
            f"{row.half.score_width:>3} {row.half.strict_triangles:>4}"
        )
    print()


def print_cross_example_findings() -> None:
    print("SYNTHESIS")
    print("=" * 96)
    print("The origin-only LRC ledger forgets two pieces of structure:")
    print("* whether the tight runner has one blocker or a two-sided pair of blockers;")
    print("* whether blocker dependence peels like a DAG or closes into an SCC.")
    print()
    print("For a proof search, the useful invariant is not just max_i d1(i).")
    print("Keep the ordered pair (d1(i), d2(i)) and the deletion relief")
    print("  relief_i(j) = nearest_i(after deleting j) - nearest_i.")
    print("Those reliefs form an incomplete tournament whose cyclic core can be")
    print("compared directly with the old endpoint-protection core.")
    print()
    print("A concrete n=14 target suggested by this pass:")
    print("  Any attempted all-protected endpoint row should induce either")
    print("  (a) a pressure leaf, giving a private runner/endpoint to peel, or")
    print("  (b) a strict pressure SCC whose labels must satisfy the same kind of")
    print("      arithmetic cycle constraint as THM-365's endpoint cycles.")
    print()
    print("The second-nearest runner distance is the missing datum.  Without it,")
    print("deleting a nearest neighbor has no measurable relief, so nearest-neighbor")
    print("graphs alone cannot distinguish a rigid two-sided moat from a fragile")
    print("one-blocker moat.")


def main() -> None:
    print("codex-2026-05-31-S470")
    print("LRC pairwise distance tournaments and two-neighbor blocker pressure")
    print()
    print_methodology()
    for example in examples():
        print_example(example)
    print_cross_example_findings()


if __name__ == "__main__":
    main()
