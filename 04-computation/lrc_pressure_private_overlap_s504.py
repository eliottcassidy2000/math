#!/usr/bin/env python3
"""
lrc_pressure_private_overlap_s504.py

codex-2026-06-01 S504

Connect the S500 pressure-DAG observation to the S362 endpoint-peeling
certificate.

Tournament Analysis declaration:

* pairwise observable: for runners a,b at time t, compare their deletion
  reliefs on nearest-neighbor distance;
* switch/gauge: orient a -> b when deleting a helps b more than deleting b
  helps a;
* tie Hamiltonian path: ties are not forced here, because the certificate
  uses the strict pressure DAG and its topological layers.  For displayed
  layer signatures, speeds are listed in the fixed numerical Hamiltonian path.

The script extracts transitive reductions of strict pressure DAGs and compares
their source/sink layers with the first endpoint-peeling owner layer.  A high
overlap would support an induction proof route; a low overlap means the current
pressure gauge is not yet the endpoint-core realization promised by THM-380.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()
S490 = SourceFileLoader(
    "lrc_n14_n18_tournament_feedback_s490",
    str(ROOT / "04-computation" / "lrc_n14_n18_tournament_feedback_s490.py"),
).load_module()
S500 = SourceFileLoader(
    "lrc_pressure_dag_audit_s500",
    str(ROOT / "04-computation" / "lrc_pressure_dag_audit_s500.py"),
).load_module()


@dataclass(frozen=True)
class EndpointLayer:
    step: int
    endpoint_count: int
    interval_count: int
    endpoint_owner_speeds: tuple[int, ...]
    interval_owner_speeds: tuple[int, ...]
    private_owner_speeds: tuple[int, ...]
    shared_owner_speeds: tuple[int, ...]


@dataclass(frozen=True)
class PressureLayerRow:
    t: Fraction
    arcs: int
    transitive_reduction_arcs: int
    is_dag: bool
    height: int
    longest_path: int
    source_speeds: tuple[int, ...]
    sink_speeds: tuple[int, ...]
    first_pressure_layer: tuple[int, ...]
    last_pressure_layer: tuple[int, ...]
    source_interval_hit: bool
    source_private_hit: bool
    sink_interval_hit: bool
    sink_private_hit: bool


@dataclass(frozen=True)
class CaseReport:
    label: str
    n: int
    sample_count: int
    dag_count: int
    mean_tr_ratio: Fraction
    source_interval_hits: int
    source_private_hits: int
    sink_interval_hits: int
    sink_private_hits: int
    first_endpoint_layer: EndpointLayer
    source_speed_hist: tuple[tuple[int, int], ...]
    sink_speed_hist: tuple[tuple[int, int], ...]
    representative: PressureLayerRow


def fmt(value: Fraction | None) -> str:
    return S500.fmt(value)


def fmt_pct(value: Fraction) -> str:
    return f"{float(100 * value):5.1f}%"


def speed_tuple(values) -> tuple[int, ...]:
    return tuple(sorted(set(int(value) for value in values)))


def endpoint_peel_layers(speeds: tuple[int, ...]) -> tuple[EndpointLayer, ...]:
    endpoints, intervals, owners, protectors, boundary = S362.build_endpoint_system(speeds)
    remaining_endpoints = set(endpoints)
    remaining_intervals = set(intervals)
    layers: list[EndpointLayer] = []
    step = 0

    while True:
        dead_endpoints = {
            endpoint
            for endpoint in remaining_endpoints
            if not (protectors[endpoint] & remaining_intervals)
        }
        if not dead_endpoints:
            break

        dead_intervals = {
            interval
            for interval in remaining_intervals
            if boundary[interval] & dead_endpoints
        }

        endpoint_owner_speeds = []
        private_owner_speeds = []
        shared_owner_speeds = []
        for endpoint in dead_endpoints:
            active_owners = owners[endpoint] & remaining_intervals
            endpoint_owner_speeds.extend(interval.speed for interval in active_owners)
            if len(active_owners) == 1:
                private_owner_speeds.extend(interval.speed for interval in active_owners)
            elif len(active_owners) > 1:
                shared_owner_speeds.extend(interval.speed for interval in active_owners)

        layers.append(
            EndpointLayer(
                step=step,
                endpoint_count=len(dead_endpoints),
                interval_count=len(dead_intervals),
                endpoint_owner_speeds=speed_tuple(endpoint_owner_speeds),
                interval_owner_speeds=speed_tuple(interval.speed for interval in dead_intervals),
                private_owner_speeds=speed_tuple(private_owner_speeds),
                shared_owner_speeds=speed_tuple(shared_owner_speeds),
            )
        )

        remaining_endpoints -= dead_endpoints
        remaining_intervals -= dead_intervals
        step += 1

    if not layers:
        return (
            EndpointLayer(
                step=0,
                endpoint_count=0,
                interval_count=0,
                endpoint_owner_speeds=tuple(),
                interval_owner_speeds=tuple(),
                private_owner_speeds=tuple(),
                shared_owner_speeds=tuple(),
            ),
        )
    return tuple(layers)


def adjacency(vertex_count: int, arcs: set[tuple[int, int]]) -> list[list[int]]:
    adj = [[] for _ in range(vertex_count)]
    for a, b in arcs:
        adj[a].append(b)
    return adj


def reaches_without(
    vertex_count: int,
    arcs: set[tuple[int, int]],
    start: int,
    target: int,
    forbidden_edge: tuple[int, int],
) -> bool:
    adj = adjacency(vertex_count, arcs - {forbidden_edge})
    todo = deque([start])
    seen = {start}
    while todo:
        v = todo.popleft()
        for w in adj[v]:
            if w == target:
                return True
            if w not in seen:
                seen.add(w)
                todo.append(w)
    return False


def transitive_reduction(
    vertex_count: int, arcs: set[tuple[int, int]]
) -> set[tuple[int, int]]:
    return {
        edge
        for edge in arcs
        if not reaches_without(vertex_count, arcs, edge[0], edge[1], edge)
    }


def layer_speeds(layer: tuple[int, ...], runner_speeds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(runner_speeds[index] for index in layer if runner_speeds[index] != 0)


def pressure_row(
    speeds: tuple[int, ...], t: Fraction, endpoint_layer: EndpointLayer
) -> PressureLayerRow:
    runner_speeds = (0,) + speeds
    arcs, _ties, _profs, _dist, _relief_gaps = S500.pressure_arcs(speeds, t)
    is_dag, layers, longest = S500.topological_layers(len(runner_speeds), arcs)
    tr = transitive_reduction(len(runner_speeds), arcs) if is_dag else set(arcs)
    sources, sinks, _source_count, _sink_count = S500.source_sink_speeds(
        runner_speeds, len(runner_speeds), arcs
    )
    source_speeds = tuple(speed for speed in sources if speed != 0)
    sink_speeds = tuple(speed for speed in sinks if speed != 0)
    first_layer = layer_speeds(layers[0], runner_speeds) if layers else tuple()
    last_layer = layer_speeds(layers[-1], runner_speeds) if layers else tuple()

    source_set = set(source_speeds)
    sink_set = set(sink_speeds)
    interval_set = set(endpoint_layer.interval_owner_speeds)
    private_set = set(endpoint_layer.private_owner_speeds)

    return PressureLayerRow(
        t=t,
        arcs=len(arcs),
        transitive_reduction_arcs=len(tr),
        is_dag=is_dag,
        height=len(layers),
        longest_path=longest,
        source_speeds=source_speeds,
        sink_speeds=sink_speeds,
        first_pressure_layer=first_layer,
        last_pressure_layer=last_layer,
        source_interval_hit=bool(source_set & interval_set),
        source_private_hit=bool(source_set & private_set),
        sink_interval_hit=bool(sink_set & interval_set),
        sink_private_hit=bool(sink_set & private_set),
    )


def summarize_case(label: str, speeds: tuple[int, ...]) -> CaseReport:
    endpoint_layer = endpoint_peel_layers(speeds)[0]
    rows = [pressure_row(speeds, t, endpoint_layer) for t in S500.sample_times(speeds)]
    source_hist: Counter[int] = Counter()
    sink_hist: Counter[int] = Counter()
    for row in rows:
        source_hist.update(row.source_speeds)
        sink_hist.update(row.sink_speeds)
    tr_ratio_sum = sum(
        (
            Fraction(row.transitive_reduction_arcs, row.arcs)
            if row.arcs
            else Fraction(0)
            for row in rows
        ),
        Fraction(0),
    )
    representative = max(
        rows,
        key=lambda row: (
            row.longest_path,
            row.height,
            row.transitive_reduction_arcs,
            row.arcs,
            -row.t,
        ),
    )
    return CaseReport(
        label=label,
        n=len(speeds) + 1,
        sample_count=len(rows),
        dag_count=sum(1 for row in rows if row.is_dag),
        mean_tr_ratio=tr_ratio_sum / len(rows),
        source_interval_hits=sum(1 for row in rows if row.source_interval_hit),
        source_private_hits=sum(1 for row in rows if row.source_private_hit),
        sink_interval_hits=sum(1 for row in rows if row.sink_interval_hit),
        sink_private_hits=sum(1 for row in rows if row.sink_private_hit),
        first_endpoint_layer=endpoint_layer,
        source_speed_hist=tuple(sorted(source_hist.items())),
        sink_speed_hist=tuple(sorted(sink_hist.items())),
        representative=representative,
    )


def hist_top(hist: tuple[tuple[int, int], ...], limit: int = 5) -> str:
    top = sorted(hist, key=lambda item: (-item[1], item[0]))[:limit]
    return ",".join(f"{speed}:{count}" for speed, count in top)


def print_method() -> None:
    print("LRC pressure/private-overlap audit (codex-2026-06-01 S504)")
    print("=" * 112)
    print("Question: do strict pressure DAG sources already point at endpoint-peel")
    print("owners, or is the current pressure gauge missing endpoint labels?")
    print()
    print("For each hard row:")
    print("  1. compute the first exact endpoint-peeling layer;")
    print("  2. sample the strict deletion-relief pressure DAG;")
    print("  3. compute the DAG transitive reduction;")
    print("  4. compare source/sink speeds with interval/private owner speeds.")
    print()


def print_endpoint_layers(rows: tuple[CaseReport, ...]) -> None:
    print("FIRST ENDPOINT-PEEL OWNER LAYERS")
    print("=" * 112)
    print(
        f"{'case':<12} {'n':>3} {'deadE':>6} {'deadI':>6} "
        f"{'interval owners':<30} {'private owners'}"
    )
    print("-" * 112)
    for row in rows:
        layer = row.first_endpoint_layer
        interval = ",".join(str(v) for v in layer.interval_owner_speeds[:14])
        private = ",".join(str(v) for v in layer.private_owner_speeds[:14])
        if len(layer.interval_owner_speeds) > 14:
            interval += ",..."
        if len(layer.private_owner_speeds) > 14:
            private += ",..."
        print(
            f"{row.label:<12} {row.n:>3} {layer.endpoint_count:>6} "
            f"{layer.interval_count:>6} {interval:<30} {private}"
        )
    print()


def print_overlap_summary(rows: tuple[CaseReport, ...]) -> None:
    print("PRESSURE DAG OVERLAP SUMMARY")
    print("=" * 112)
    print(
        f"{'case':<12} {'samples':>7} {'DAGs':>9} {'mean TR':>8} "
        f"{'src-I':>8} {'src-P':>8} {'sink-I':>8} {'sink-P':>8} "
        f"{'top sources':<22} {'top sinks'}"
    )
    print("-" * 112)
    for row in rows:
        total = row.sample_count
        print(
            f"{row.label:<12} {total:>7} {row.dag_count:>4}/{total:<4} "
            f"{fmt_pct(row.mean_tr_ratio):>8} "
            f"{row.source_interval_hits:>4}/{total:<3} "
            f"{row.source_private_hits:>4}/{total:<3} "
            f"{row.sink_interval_hits:>4}/{total:<3} "
            f"{row.sink_private_hits:>4}/{total:<3} "
            f"{hist_top(row.source_speed_hist):<22} {hist_top(row.sink_speed_hist)}"
        )
    print()


def print_representatives(rows: tuple[CaseReport, ...]) -> None:
    print("REPRESENTATIVE TRANSITIVE-REDUCED PRESSURE DAGS")
    print("=" * 112)
    print(
        f"{'case':<12} {'t':>12} {'arcs/TR':>10} {'height':>6} "
        f"{'longest':>7} {'sources':<22} {'sinks':<22} {'first->last'}"
    )
    print("-" * 112)
    for summary in rows:
        row = summary.representative
        sources = ",".join(str(v) for v in row.source_speeds[:8])
        sinks = ",".join(str(v) for v in row.sink_speeds[:8])
        first = ",".join(str(v) for v in row.first_pressure_layer[:8])
        last = ",".join(str(v) for v in row.last_pressure_layer[:8])
        print(
            f"{summary.label:<12} {fmt(row.t):>12} "
            f"{row.arcs:>4}/{row.transitive_reduction_arcs:<5} "
            f"{row.height:>6} {row.longest_path:>7} "
            f"{sources:<22} {sinks:<22} {first}->{last}"
        )
    print()


def print_interpretation(rows: tuple[CaseReport, ...]) -> None:
    total = sum(row.sample_count for row in rows)
    source_private = sum(row.source_private_hits for row in rows)
    source_interval = sum(row.source_interval_hits for row in rows)
    sink_interval = sum(row.sink_interval_hits for row in rows)
    sink_private = sum(row.sink_private_hits for row in rows)
    print("INTERPRETATION")
    print("=" * 112)
    print(
        f"Source/private overlap across all samples is {source_private}/{total}; "
        f"source/interval overlap is {source_interval}/{total}."
    )
    print(
        f"Sink/interval overlap is {sink_interval}/{total}; sink/private overlap "
        f"is {sink_private}/{total}.  Thus the current strict deletion-relief "
        "gauge is a strong acyclic pressure certificate, but it is not yet a "
        "direct endpoint-owner realization certificate."
    )
    print()
    print("Proof consequence:")
    print("  pressure DAG acyclicity is real, but THM-380 still needs a labelled")
    print("  lift from endpoint protection incidences to pressure arrows.  The")
    print("  transitive reduction should be used as the minimal relation to label,")
    print("  not a scalar collapse of the pressure shadow.")


def main() -> None:
    print_method()
    cases = (
        ("n14-d7", S490.scale_speeds(14, 7, 6)),
        ("n14-d14", S490.scale_speeds(14, 14, 6)),
        ("n18-d3", S490.scale_speeds(18, 3, 8)),
        ("n18-d9", S490.scale_speeds(18, 9, 8)),
        ("n18-d18", S490.scale_speeds(18, 18, 8)),
    )
    reports = tuple(summarize_case(label, speeds) for label, speeds in cases)
    print_endpoint_layers(reports)
    print_overlap_summary(reports)
    print_representatives(reports)
    print_interpretation(reports)


if __name__ == "__main__":
    main()
