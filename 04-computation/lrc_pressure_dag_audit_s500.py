#!/usr/bin/env python3
"""
lrc_pressure_dag_audit_s500.py

codex-2026-06-01 S500

Investigate the concept that LRC pressure searches are returning DAGs.

The pressure relation is the two-neighbor deletion-relief graph from S470:
orient a -> b when deleting a improves b's nearest-neighbor moat more than
deleting b improves a's moat.  If the strict pressure graph is a DAG, the
search did not merely fail to find a cycle; it returned a topological
certificate candidate, with sources as structural blockers and sinks as
fragile blocked runners.
"""

from __future__ import annotations

from collections import Counter, deque
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
S470 = SourceFileLoader(
    "lrc_pairwise_tournament_s470",
    str(ROOT / "04-computation" / "lrc_pairwise_tournament_s470.py"),
).load_module()
S490 = SourceFileLoader(
    "lrc_n14_n18_tournament_feedback_s490",
    str(ROOT / "04-computation" / "lrc_n14_n18_tournament_feedback_s490.py"),
).load_module()


@dataclass(frozen=True)
class DagMetrics:
    t: Fraction
    n: int
    arcs: int
    ties: int
    is_dag: bool
    largest_scc: int
    source_count: int
    sink_count: int
    height: int
    layer_widths: tuple[int, ...]
    longest_path_edges: int
    d1_descent: Fraction
    relief_gap_average: Fraction
    source_speeds: tuple[int, ...]
    sink_speeds: tuple[int, ...]


@dataclass(frozen=True)
class CaseSummary:
    label: str
    n: int
    sample_count: int
    dag_count: int
    cyclic_count: int
    max_largest_scc: int
    max_height: int
    max_longest_path_edges: int
    max_arcs: int
    source_range: tuple[int, int]
    sink_range: tuple[int, int]
    mean_d1_descent: Fraction
    mean_relief_gap: Fraction
    representative: DagMetrics


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_pct(value: Fraction) -> str:
    return f"{float(100 * value):5.1f}%"


def pressure_arcs(
    speeds: tuple[int, ...], t: Fraction
) -> tuple[set[tuple[int, int]], int, tuple[object, ...], tuple[tuple[Fraction, ...], ...], list[Fraction]]:
    runner_speeds = (0,) + speeds
    pos = S470.positions(runner_speeds, t)
    dist = S470.distance_matrix(pos)
    profs = S470.profiles(runner_speeds, t)
    arcs: set[tuple[int, int]] = set()
    relief_gaps: list[Fraction] = []
    ties = 0
    for i, j in combinations(range(len(runner_speeds)), 2):
        relief_i_from_j = S470.nearest_without(dist, i, j) - profs[i].d1
        relief_j_from_i = S470.nearest_without(dist, j, i) - profs[j].d1
        if relief_i_from_j > relief_j_from_i:
            arcs.add((j, i))
            relief_gaps.append(relief_i_from_j - relief_j_from_i)
        elif relief_j_from_i > relief_i_from_j:
            arcs.add((i, j))
            relief_gaps.append(relief_j_from_i - relief_i_from_j)
        else:
            ties += 1
    return arcs, ties, profs, dist, relief_gaps


def topological_layers(n: int, arcs: set[tuple[int, int]]) -> tuple[bool, tuple[tuple[int, ...], ...], int]:
    adj = [[] for _ in range(n)]
    indeg = [0] * n
    for a, b in arcs:
        adj[a].append(b)
        indeg[b] += 1

    q = deque([i for i, d in enumerate(indeg) if d == 0])
    layers: list[tuple[int, ...]] = []
    seen = 0
    longest = [0] * n
    while q:
        layer = tuple(sorted(q))
        layers.append(layer)
        next_q: deque[int] = deque()
        while q:
            v = q.popleft()
            seen += 1
            for w in adj[v]:
                if longest[w] < longest[v] + 1:
                    longest[w] = longest[v] + 1
                indeg[w] -= 1
                if indeg[w] == 0:
                    next_q.append(w)
        q = next_q
    return seen == n, tuple(layers), max(longest) if longest else 0


def source_sink_speeds(
    runner_speeds: tuple[int, ...], n: int, arcs: set[tuple[int, int]]
) -> tuple[tuple[int, ...], tuple[int, ...], int, int]:
    out = [0] * n
    indeg = [0] * n
    for a, b in arcs:
        out[a] += 1
        indeg[b] += 1
    active = [i for i in range(n) if out[i] + indeg[i] > 0]
    sources = tuple(runner_speeds[i] for i in active if indeg[i] == 0)
    sinks = tuple(runner_speeds[i] for i in active if out[i] == 0)
    return sources, sinks, len(sources), len(sinks)


def dag_metrics(speeds: tuple[int, ...], t: Fraction) -> DagMetrics:
    runner_speeds = (0,) + speeds
    n = len(runner_speeds)
    arcs, ties, profs, _dist, relief_gaps = pressure_arcs(speeds, t)
    is_dag, layers, longest = topological_layers(n, arcs)
    scc = S470.scc_sizes(n, arcs)
    sources, sinks, source_count, sink_count = source_sink_speeds(runner_speeds, n, arcs)
    d1_descent_count = sum(1 for a, b in arcs if profs[a].d1 >= profs[b].d1)
    d1_descent = Fraction(d1_descent_count, len(arcs)) if arcs else Fraction(1)
    relief_gap_average = sum(relief_gaps, Fraction(0)) / len(relief_gaps) if relief_gaps else Fraction(0)
    return DagMetrics(
        t=t,
        n=n,
        arcs=len(arcs),
        ties=ties,
        is_dag=is_dag,
        largest_scc=scc[0] if scc else 0,
        source_count=source_count,
        sink_count=sink_count,
        height=len(layers),
        layer_widths=tuple(len(layer) for layer in layers),
        longest_path_edges=longest,
        d1_descent=d1_descent,
        relief_gap_average=relief_gap_average,
        source_speeds=sources,
        sink_speeds=sinks,
    )


def sample_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    n = len(speeds) + 1
    if n <= 14:
        return S470.candidate_times(speeds)
    return S490.bounded_times(speeds)


def summarize_case(label: str, speeds: tuple[int, ...]) -> CaseSummary:
    rows = [dag_metrics(speeds, t) for t in sample_times(speeds)]
    representative = max(
        rows,
        key=lambda row: (
            row.largest_scc,
            row.longest_path_edges,
            row.height,
            row.arcs,
            -row.t,
        ),
    )
    return CaseSummary(
        label=label,
        n=len(speeds) + 1,
        sample_count=len(rows),
        dag_count=sum(1 for row in rows if row.is_dag),
        cyclic_count=sum(1 for row in rows if not row.is_dag),
        max_largest_scc=max(row.largest_scc for row in rows),
        max_height=max(row.height for row in rows),
        max_longest_path_edges=max(row.longest_path_edges for row in rows),
        max_arcs=max(row.arcs for row in rows),
        source_range=(min(row.source_count for row in rows), max(row.source_count for row in rows)),
        sink_range=(min(row.sink_count for row in rows), max(row.sink_count for row in rows)),
        mean_d1_descent=sum((row.d1_descent for row in rows), Fraction(0)) / len(rows),
        mean_relief_gap=sum((row.relief_gap_average for row in rows), Fraction(0)) / len(rows),
        representative=representative,
    )


def print_method() -> None:
    print("PRESSURE-DAG METHOD")
    print("=" * 108)
    print("Pressure edges compare deletion reliefs: a -> b means a is the more")
    print("irreplaceable blocker of b's nearest-neighbor moat.")
    print()
    print("If the strict pressure graph is a DAG, the search has returned:")
    print("  source layer: blockers that can be charged first")
    print("  sink layer: most blocked runners")
    print("  topological height: length of dependency cascade")
    print("  no SCC core: no mutual pressure knot to support a disproof candidate")
    print()


def print_summary(rows: tuple[CaseSummary, ...]) -> None:
    print("DAG AUDIT ACROSS PRESSURE SEARCH SAMPLES")
    print("=" * 108)
    print(
        f"{'case':<14} {'n':>3} {'samples':>7} {'DAGs':>9} {'cyclic':>7} "
        f"{'maxSCC':>6} {'height':>6} {'longest':>7} {'max arcs':>8} "
        f"{'src range':>11} {'sink range':>11} {'d1-desc'}"
    )
    print("-" * 108)
    for row in rows:
        print(
            f"{row.label:<14} {row.n:>3} {row.sample_count:>7} "
            f"{row.dag_count:>4}/{row.sample_count:<4} {row.cyclic_count:>7} "
            f"{row.max_largest_scc:>6} {row.max_height:>6} "
            f"{row.max_longest_path_edges:>7} {row.max_arcs:>8} "
            f"{str(row.source_range):>11} {str(row.sink_range):>11} "
            f"{fmt_pct(row.mean_d1_descent):>7}"
        )
    print()


def print_representatives(rows: tuple[CaseSummary, ...]) -> None:
    print("MAXIMAL DEPENDENCY DAGS")
    print("=" * 108)
    print(
        f"{'case':<14} {'t':>12} {'arcs/ties':>10} {'height':>6} "
        f"{'longest':>7} {'layers':<24} {'sources':<18} {'sinks'}"
    )
    print("-" * 108)
    for summary in rows:
        row = summary.representative
        layers = ",".join(str(x) for x in row.layer_widths[:8])
        if len(row.layer_widths) > 8:
            layers += ",..."
        sources = ",".join(str(x) for x in row.source_speeds[:6])
        sinks = ",".join(str(x) for x in row.sink_speeds[:6])
        print(
            f"{summary.label:<14} {fmt(row.t):>12} "
            f"{row.arcs:>4}/{row.ties:<5} {row.height:>6} "
            f"{row.longest_path_edges:>7} {layers:<24} "
            f"{sources:<18} {sinks}"
        )
    print()


def print_interpretation(rows: tuple[CaseSummary, ...]) -> None:
    cyclic_total = sum(row.cyclic_count for row in rows)
    sample_total = sum(row.sample_count for row in rows)
    print("INTERPRETATION")
    print("=" * 108)
    print(f"Across these cases, cyclic pressure samples = {cyclic_total}/{sample_total}.")
    print(
        "A pressure search returning DAGs should be read as a failed counterexample "
        "search but a successful certificate search: the obstruction has a "
        "topological order."
    )
    print()
    print("Proof-useful reframe:")
    print("  1. pressure sources are chargeable blockers;")
    print("  2. pressure sinks are exposed blocked runners;")
    print("  3. topological layers give an induction order;")
    print("  4. the first non-DAG sample would be the first serious disproof-like core.")
    print()
    print(
        "The d1-descent column is not perfect, so the DAG is not merely sorting by "
        "nearest distance.  It is a real relation-level residue of the two-neighbor "
        "geometry."
    )


def main() -> None:
    print("LRC pressure DAG audit (codex-2026-06-01 S500)")
    print()
    print_method()
    cases = (
        ("n14-d7", S490.scale_speeds(14, 7, 6)),
        ("n14-d14", S490.scale_speeds(14, 14, 6)),
        ("n18-d3", S490.scale_speeds(18, 3, 8)),
        ("n18-d9", S490.scale_speeds(18, 9, 8)),
        ("n18-d18", S490.scale_speeds(18, 18, 8)),
    )
    summaries = tuple(summarize_case(label, speeds) for label, speeds in cases)
    print_summary(summaries)
    print_representatives(summaries)
    print_interpretation(summaries)


if __name__ == "__main__":
    main()
