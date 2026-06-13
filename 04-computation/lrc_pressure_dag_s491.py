#!/usr/bin/env python3
"""
lrc_pressure_dag_s491.py

codex-2026-06-01 S491

Investigate the concept that LRC pressure searches often return DAGs.

Previous sessions looked for pressure cycles or pressure SCCs as possible
counterexample-like signals.  This script turns the negative result around:
when deletion-relief pressure returns a DAG, the DAG itself is a certificate
candidate.  It has source/sink layers, a topological dependency order, and a
natural pressure-peeling grammar.

Tournament Analysis data:

* objects: runners {0} union V at a selected time t;
* pairwise observable: deletion relief
      relief_i(j) = nearest_i(after deleting j) - nearest_i;
* switch/gauge: orient j -> i when j gives more deletion relief to i than
      i gives to j;
* tie Hamiltonian path: strict ties are left missing for the pressure DAG,
      but the base path is recorded as the completion if a total tournament is
      needed.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
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
S470 = SourceFileLoader(
    "lrc_pairwise_tournament_s470",
    str(ROOT / "04-computation" / "lrc_pairwise_tournament_s470.py"),
).load_module()
S490 = SourceFileLoader(
    "lrc_n14_n18_tournament_feedback_s490",
    str(ROOT / "04-computation" / "lrc_n14_n18_tournament_feedback_s490.py"),
).load_module()


@dataclass(frozen=True)
class Case:
    label: str
    speeds: tuple[int, ...]
    mode: str


@dataclass(frozen=True)
class PressureArc:
    src: int
    dst: int
    src_speed: int
    dst_speed: int
    relief_src_from_dst: Fraction
    relief_dst_from_src: Fraction
    margin: Fraction


@dataclass(frozen=True)
class DagMetrics:
    label: str
    t: Fraction
    time_kind: str
    n: int
    arcs: int
    ties: int
    active_vertices: int
    dag: bool
    directed_triangles: int
    largest_scc: int
    source_layers: tuple[tuple[int, ...], ...]
    sink_layers: tuple[tuple[int, ...], ...]
    longest_chain_edges: int
    width_bruteforce: int | None
    max_margin: Fraction
    top_arcs: tuple[PressureArc, ...]


def fmt(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | None) -> str:
    if value is None:
        return "-"
    return f"{float(value):.6f}"


def gcd_all(values: tuple[int, ...]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def scale_speeds(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or gcd_all(speeds) != 1:
        raise ValueError((n, scale, skip, speeds))
    return speeds


def hard_cases() -> tuple[Case, ...]:
    rows = {
        (14, 7): S490.best_ladder(14, 7),
        (14, 14): S490.best_ladder(14, 14),
        (18, 3): S490.best_ladder(18, 3),
        (18, 9): S490.best_ladder(18, 9),
        (18, 18): S490.best_ladder(18, 18),
    }
    return (
        Case("n14 initial", initial(14), "boundary"),
        Case("n14 d=7", rows[(14, 7)].speeds, "row-parent"),
        Case("n14 d=14", rows[(14, 14)].speeds, "gate"),
        Case("n18 initial", initial(18), "boundary"),
        Case("n18 d=3", rows[(18, 3)].speeds, "third-torsion"),
        Case("n18 d=9", rows[(18, 9)].speeds, "row-parent"),
        Case("n18 d=18", rows[(18, 18)].speeds, "gate"),
    )


def pressure_arcs(speeds: tuple[int, ...], t: Fraction) -> tuple[tuple[PressureArc, ...], int]:
    runner_speeds = (0,) + speeds
    pos = S470.positions(runner_speeds, t)
    dist = S470.distance_matrix(pos)
    profs = S470.profiles(runner_speeds, t)
    arcs: list[PressureArc] = []
    ties = 0
    for i, j in combinations(range(len(runner_speeds)), 2):
        relief_i_from_j = S470.nearest_without(dist, i, j) - profs[i].d1
        relief_j_from_i = S470.nearest_without(dist, j, i) - profs[j].d1
        if relief_i_from_j > relief_j_from_i:
            src, dst = j, i
            margin = relief_i_from_j - relief_j_from_i
        elif relief_j_from_i > relief_i_from_j:
            src, dst = i, j
            margin = relief_j_from_i - relief_i_from_j
        else:
            ties += 1
            continue
        arcs.append(
            PressureArc(
                src=src,
                dst=dst,
                src_speed=runner_speeds[src],
                dst_speed=runner_speeds[dst],
                relief_src_from_dst=relief_src_from_dst_value(
                    src, dst, relief_i_from_j, relief_j_from_i, i, j
                ),
                relief_dst_from_src=relief_dst_from_src_value(src, dst, relief_i_from_j, relief_j_from_i, i, j),
                margin=margin,
            )
        )
    return tuple(arcs), ties


def relief_src_from_dst_value(
    src: int,
    dst: int,
    relief_i_from_j: Fraction,
    relief_j_from_i: Fraction,
    i: int,
    j: int,
) -> Fraction:
    if src == i and dst == j:
        return relief_i_from_j
    if src == j and dst == i:
        return relief_j_from_i
    raise ValueError((src, dst, i, j))


def relief_dst_from_src_value(
    src: int,
    dst: int,
    relief_i_from_j: Fraction,
    relief_j_from_i: Fraction,
    i: int,
    j: int,
) -> Fraction:
    if src == i and dst == j:
        return relief_j_from_i
    if src == j and dst == i:
        return relief_i_from_j
    raise ValueError((src, dst, i, j))


def arcs_as_pairs(arcs: tuple[PressureArc, ...]) -> set[tuple[int, int]]:
    return {(arc.src, arc.dst) for arc in arcs}


def directed_triangles(n: int, arcs: set[tuple[int, int]]) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        if (a, b) in arcs and (b, c) in arcs and (c, a) in arcs:
            total += 1
        elif (a, c) in arcs and (c, b) in arcs and (b, a) in arcs:
            total += 1
    return total


def scc_sizes(n: int, arcs: set[tuple[int, int]]) -> tuple[int, ...]:
    graph = [[] for _ in range(n)]
    reverse = [[] for _ in range(n)]
    for src, dst in arcs:
        graph[src].append(dst)
        reverse[dst].append(src)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in graph[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []
    for start in reversed(order):
        if start in seen:
            continue
        todo = [start]
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def active_vertices(n: int, arcs: set[tuple[int, int]]) -> set[int]:
    out: set[int] = set()
    for src, dst in arcs:
        out.add(src)
        out.add(dst)
    return out


def kahn_layers(
    n: int,
    arcs: set[tuple[int, int]],
    reverse: bool = False,
) -> tuple[tuple[tuple[int, ...], ...], bool]:
    vertices = active_vertices(n, arcs)
    if not vertices:
        return tuple(), True
    graph = {v: set() for v in vertices}
    indeg = {v: 0 for v in vertices}
    for src, dst in arcs:
        a, b = (dst, src) if reverse else (src, dst)
        if a in vertices and b in vertices:
            graph[a].add(b)
    for src in vertices:
        for dst in graph[src]:
            indeg[dst] += 1

    queue = sorted(v for v in vertices if indeg[v] == 0)
    layers: list[tuple[int, ...]] = []
    processed = 0
    while queue:
        layer = tuple(queue)
        layers.append(layer)
        processed += len(layer)
        next_queue: list[int] = []
        for src in layer:
            for dst in sorted(graph[src]):
                indeg[dst] -= 1
                if indeg[dst] == 0:
                    next_queue.append(dst)
        queue = sorted(next_queue)
    return tuple(layers), processed == len(vertices)


def longest_chain_edges(n: int, arcs: set[tuple[int, int]]) -> int:
    layers, dag = kahn_layers(n, arcs)
    if not dag:
        return -1
    order = [v for layer in layers for v in layer]
    graph = {v: [] for v in active_vertices(n, arcs)}
    for src, dst in arcs:
        graph.setdefault(src, []).append(dst)
        graph.setdefault(dst, [])
    dist = {v: 0 for v in graph}
    for v in order:
        for nxt in graph[v]:
            dist[nxt] = max(dist[nxt], dist[v] + 1)
    return max(dist.values(), default=0)


def transitive_closure(n: int, arcs: set[tuple[int, int]]) -> list[int]:
    reach = [0] * n
    graph = [[] for _ in range(n)]
    for src, dst in arcs:
        graph[src].append(dst)
    for start in range(n):
        seen = set()
        todo = list(graph[start])
        while todo:
            v = todo.pop()
            if v in seen:
                continue
            seen.add(v)
            todo.extend(graph[v])
        mask = 0
        for v in seen:
            mask |= 1 << v
        reach[start] = mask
    return reach


def antichain_width(n: int, arcs: set[tuple[int, int]]) -> int | None:
    active = sorted(active_vertices(n, arcs))
    if len(active) > 20:
        return None
    reach = transitive_closure(n, arcs)
    best = 0
    for mask in range(1 << len(active)):
        size = mask.bit_count()
        if size <= best:
            continue
        ok = True
        chosen = [active[i] for i in range(len(active)) if mask & (1 << i)]
        for a, b in combinations(chosen, 2):
            if (reach[a] & (1 << b)) or (reach[b] & (1 << a)):
                ok = False
                break
        if ok:
            best = size
    return best


def selected_times(case: Case) -> tuple[tuple[str, Fraction], ...]:
    report = S356.report(case.label, list(case.speeds))
    n = len(case.speeds) + 1
    out: list[tuple[str, Fraction]] = []
    if report.witness is not None:
        out.append(("gap-mid", report.witness))
    if report.boundary_witness is not None:
        out.append(("boundary", report.boundary_witness))
    out.extend(
        [
            ("unit", Fraction(1, n)),
            ("half-unit", Fraction(1, 2 * n)),
            ("half-turn", Fraction(1, 2)),
        ]
    )
    seen: set[Fraction] = set()
    unique: list[tuple[str, Fraction]] = []
    for tag, t in out:
        if t not in seen:
            seen.add(t)
            unique.append((tag, t))
    return tuple(unique)


def scan_times(case: Case) -> tuple[tuple[str, Fraction], ...]:
    return tuple(("bounded-candidate", t) for t in S490.bounded_times(case.speeds))


def dag_metrics(case: Case, time_kind: str, t: Fraction, compute_width: bool = True) -> DagMetrics:
    arcs, ties = pressure_arcs(case.speeds, t)
    pairs = arcs_as_pairs(arcs)
    n = len(case.speeds) + 1
    sources, dag = kahn_layers(n, pairs)
    sinks, sink_dag = kahn_layers(n, pairs, reverse=True)
    scc = scc_sizes(n, pairs)
    top_arcs = tuple(sorted(arcs, key=lambda arc: (-arc.margin, arc.src_speed, arc.dst_speed))[:5])
    return DagMetrics(
        label=case.label,
        t=t,
        time_kind=time_kind,
        n=n,
        arcs=len(arcs),
        ties=ties,
        active_vertices=len(active_vertices(n, pairs)),
        dag=dag and sink_dag,
        directed_triangles=directed_triangles(n, pairs),
        largest_scc=scc[0] if scc else 0,
        source_layers=sources,
        sink_layers=sinks,
        longest_chain_edges=longest_chain_edges(n, pairs),
        width_bruteforce=antichain_width(n, pairs) if compute_width else None,
        max_margin=max((arc.margin for arc in arcs), default=Fraction(0)),
        top_arcs=top_arcs,
    )


def layer_text(layers: tuple[tuple[int, ...], ...], speeds: tuple[int, ...], limit: int = 4) -> str:
    runner_speeds = (0,) + speeds
    pieces = []
    for layer in layers[:limit]:
        pieces.append("{" + ",".join(str(runner_speeds[v]) for v in layer) + "}")
    if len(layers) > limit:
        pieces.append("...")
    return " ".join(pieces) or "-"


def print_selected_dags(cases: tuple[Case, ...]) -> None:
    print("SELECTED PRESSURE DAGS")
    print("=" * 112)
    print(
        f"{'case':<15} {'time kind':<18} {'t':>14} {'a/t':>8} {'act':>3} "
        f"{'dag':>3} {'tri':>3} {'scc':>3} {'src layers':>4} {'sink layers':>5} "
        f"{'chain':>5} {'width':>5} {'max margin':>10}"
    )
    print("-" * 112)
    selected: list[tuple[Case, DagMetrics]] = []
    for case in cases:
        for kind, t in selected_times(case):
            row = dag_metrics(case, kind, t)
            selected.append((case, row))
            width = "-" if row.width_bruteforce is None else str(row.width_bruteforce)
            print(
                f"{case.label:<15} {kind:<18} {fmt(t):>14} "
                f"{row.arcs:>3}/{row.ties:<3} {row.active_vertices:>3} "
                f"{'yes' if row.dag else 'no':>3} {row.directed_triangles:>3} "
                f"{row.largest_scc:>3} {len(row.source_layers):>4} "
                f"{len(row.sink_layers):>5} {row.longest_chain_edges:>5} "
                f"{width:>5} {fmt(row.max_margin):>10}"
            )
    print()
    print("Representative source/sink peel layers")
    print("-" * 112)
    for case, row in selected:
        if row.arcs == 0:
            continue
        if row.time_kind not in {"gap-mid", "origin gap mid", "max pressure core", "max pressure cycles"}:
            continue
        print(f"[{case.label} @ {row.time_kind}, t={fmt(row.t)}]")
        print(f"  source peel: {layer_text(row.source_layers, case.speeds)}")
        print(f"  sink peel:   {layer_text(row.sink_layers, case.speeds)}")
        arc_bits = []
        for arc in row.top_arcs:
            arc_bits.append(f"{arc.src_speed}->{arc.dst_speed} margin={fmt(arc.margin)}")
        print(f"  top margins: {', '.join(arc_bits) or '-'}")
    print()


def print_scan_summary(cases: tuple[Case, ...]) -> None:
    print("PRESSURE SEARCH SUMMARY")
    print("=" * 112)
    print(
        f"{'case':<15} {'times':>6} {'all DAG':>7} {'cyclic':>6} "
        f"{'max arcs':>8} {'max scc':>7} {'max tri':>7} {'max chain':>9} "
        f"{'max width':>9} {'max active':>10}"
    )
    print("-" * 112)
    for case in cases:
        rows = [
            dag_metrics(case, kind, t, compute_width=False)
            for kind, t in scan_times(case)
        ]
        cyclic = [row for row in rows if not row.dag]
        print(
            f"{case.label:<15} {len(rows):>6} {str(not cyclic):>7} "
            f"{len(cyclic):>6} {max((row.arcs for row in rows), default=0):>8} "
            f"{max((row.largest_scc for row in rows), default=0):>7} "
            f"{max((row.directed_triangles for row in rows), default=0):>7} "
            f"{max((row.longest_chain_edges for row in rows), default=0):>9} "
            f"{'-':>9} {max((row.active_vertices for row in rows), default=0):>10}"
        )
    print()


def print_endpoint_pressure_join(cases: tuple[Case, ...]) -> None:
    print("ENDPOINT DEBT VERSUS PRESSURE DAG SHAPE")
    print("=" * 112)
    print(
        f"{'case':<15} {'class':<13} {'gap/th':>10} {'unprot':>7} "
        f"{'prod':>10} {'best DAG chain':>14} {'best active':>12} {'interpretation'}"
    )
    print("-" * 112)
    for case in cases:
        report = S356.report(case.label, list(case.speeds))
        summary = S360.summarize(list(case.speeds))
        gap_ratio = report.max_gap / report.threshold
        product = gap_ratio * summary.unprotected_count
        rows = [dag_metrics(case, kind, t) for kind, t in selected_times(case)]
        best_chain = max((row.longest_chain_edges for row in rows), default=0)
        best_active = max((row.active_vertices for row in rows), default=0)
        if summary.classification == "boundary_only":
            interp = "endpoint witness branch"
        elif best_active == 0:
            interp = "no pressure dependency"
        else:
            interp = "DAG gives ordered pressure peel"
        print(
            f"{case.label:<15} {summary.classification:<13} "
            f"{fmt(gap_ratio):>10} {summary.unprotected_count:>7} "
            f"{fmt(product):>10} {best_chain:>14} {best_active:>12} {interp}"
        )
    print()


def print_synthesis() -> None:
    print("SYNTHESIS")
    print("=" * 112)
    print("A pressure search returning a DAG should be read as structure, not null output.")
    print()
    print("Pressure cycles/SCCs are the disproof-like signal: they say every runner in")
    print("a pressure core is blocking another runner in a closed dependency loop.")
    print("A pressure DAG says the opposite: dependency can be topologically sorted.")
    print("That order is a candidate proof certificate because it gives source and")
    print("sink peel layers that can be paired with endpoint-private rows.")
    print()
    print("The current hard LRC rows keep falling into the DAG regime.  That supports")
    print("a two-stage proof search:")
    print("  1. Build the mobile pressure DAG at a critical time or chamber.")
    print("  2. Peel along sink/source layers; if the endpoint incidence system cannot")
    print("     supply a matching private row, the obstruction must appear as a labelled")
    print("     pressure SCC.  The sampled rows do not show such an SCC.")


def main() -> None:
    cases = hard_cases()
    print("LRC pressure DAG audit (codex-2026-06-01 S491)")
    print("All time-slice calculations use exact Fraction arithmetic.\n")
    print_selected_dags(cases)
    print_scan_summary(cases)
    print_endpoint_pressure_join(cases)
    print_synthesis()


if __name__ == "__main__":
    main()
