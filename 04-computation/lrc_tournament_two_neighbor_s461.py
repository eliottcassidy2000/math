#!/usr/bin/env python3
"""
lrc_tournament_two_neighbor_s461.py

codex-2026-05-31 S461

Explore where tournament structure can enter the Lonely Runner endpoint
program if we keep more than the scalar nearest-runner distance.

The lift used here is deliberately finite:

* At every forbidden endpoint, record the two-sided nearest neighbors of the
  stationary runner.
* Build a speed-level handoff relation: owner speed a points to protector
  speed b when b strictly protects an endpoint owned by a.
* Compare each unordered pair of speeds by the number of endpoint labels each
  protects for the other.  Decisive comparisons form a weighted tournament
  shadow; ties mark missing arithmetic data.

This does not try to prove LRC.  It tests whether the endpoint-protection
object has the same kind of incidence/cut/handoff grammar as tournament
good-cut protection.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class DistanceRecord:
    speed: int
    position: Fraction
    distance: Fraction
    side: str


@dataclass(frozen=True)
class EndpointLift:
    value: Fraction
    owners: tuple[int, ...]
    inside: tuple[int, ...]
    on_boundary: tuple[int, ...]
    nearest_distance: Fraction
    second_distance: Fraction | None
    left_speed: int | None
    left_distance: Fraction | None
    right_speed: int | None
    right_distance: Fraction | None
    state: str


@dataclass(frozen=True)
class HandoffSummary:
    label: str
    n: int
    speeds: tuple[int, ...]
    classification: str
    threshold: Fraction
    unique_endpoints: int
    unprotected: int
    endpoint_state_hist: tuple[tuple[str, int], ...]
    owner_out_hist: tuple[tuple[int, int], ...]
    protector_in_hist: tuple[tuple[int, int], ...]
    decisive_pairs: int
    tied_pairs: int
    majority_cycles: int
    scc_sizes: tuple[int, ...]
    first_witness: Fraction | None


def fmt(x: Fraction | None) -> str:
    if x is None:
        return "-"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def circle(x: Fraction) -> Fraction:
    return x % ONE


def dist_to_zero(x: Fraction) -> Fraction:
    y = circle(x)
    return min(y, ONE - y)


def side_of_position(x: Fraction) -> str:
    y = circle(x)
    if y == 0:
        return "zero"
    if y <= Fraction(1, 2):
        return "right"
    return "left"


def distance_records(speeds: tuple[int, ...], t: Fraction) -> tuple[DistanceRecord, ...]:
    out: list[DistanceRecord] = []
    for speed in speeds:
        pos = circle(speed * t)
        out.append(
            DistanceRecord(
                speed=speed,
                position=pos,
                distance=dist_to_zero(pos),
                side=side_of_position(pos),
            )
        )
    return tuple(sorted(out, key=lambda row: (row.distance, row.speed)))


def distinct_second_distance(records: tuple[DistanceRecord, ...]) -> Fraction | None:
    if not records:
        return None
    first = records[0].distance
    for row in records:
        if row.distance > first:
            return row.distance
    return None


def side_nearest(
    records: tuple[DistanceRecord, ...], side: str
) -> tuple[int | None, Fraction | None]:
    candidates = [row for row in records if row.side == side]
    if not candidates:
        return None, None
    best = min(candidates, key=lambda row: (row.distance, row.speed))
    return best.speed, best.distance


def endpoint_lifts(speeds: tuple[int, ...]) -> tuple[EndpointLift, ...]:
    threshold = Fraction(1, len(speeds) + 1)
    by_value: dict[Fraction, list] = defaultdict(list)
    for endpoint in S360.endpoints(speeds):
        by_value[endpoint.value].append(endpoint)

    rows: list[EndpointLift] = []
    for value, labels in sorted(by_value.items()):
        records = distance_records(speeds, value)
        inside = tuple(row.speed for row in records if row.distance < threshold)
        on_boundary = tuple(row.speed for row in records if row.distance == threshold)
        owners = tuple(sorted({label.speed for label in labels}))
        left_speed, left_distance = side_nearest(records, "left")
        right_speed, right_distance = side_nearest(records, "right")
        if not inside:
            state = "witness"
        elif len(inside) == 1:
            state = "single_handoff"
        else:
            state = "redundant_handoff"
        rows.append(
            EndpointLift(
                value=value,
                owners=owners,
                inside=inside,
                on_boundary=on_boundary,
                nearest_distance=records[0].distance,
                second_distance=distinct_second_distance(records),
                left_speed=left_speed,
                left_distance=left_distance,
                right_speed=right_speed,
                right_distance=right_distance,
                state=state,
            )
        )
    return tuple(rows)


def protection_counts(speeds: tuple[int, ...]) -> dict[tuple[int, int], int]:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for endpoint in S360.endpoints(speeds):
        owner = endpoint.speed
        for protector in speeds:
            if protector != owner and S360.direct_protects(speeds, protector, endpoint.value):
                counts[(owner, protector)] += 1
    return counts


def majority_edges(
    speeds: tuple[int, ...], counts: dict[tuple[int, int], int]
) -> tuple[set[tuple[int, int]], int]:
    edges: set[tuple[int, int]] = set()
    ties = 0
    for a, b in combinations(speeds, 2):
        ab = counts.get((a, b), 0)
        ba = counts.get((b, a), 0)
        if ab > ba:
            edges.add((a, b))
        elif ba > ab:
            edges.add((b, a))
        else:
            ties += 1
    return edges, ties


def cycle_count_3(speeds: tuple[int, ...], edges: set[tuple[int, int]]) -> int:
    total = 0
    edge_set = set(edges)
    for a, b, c in combinations(speeds, 3):
        tri = {(a, b), (b, c), (c, a)}
        rev = {(a, c), (c, b), (b, a)}
        if tri <= edge_set or rev <= edge_set:
            total += 1
    return total


def scc_sizes(speeds: tuple[int, ...], counts: dict[tuple[int, int], int]) -> tuple[int, ...]:
    graph = {speed: set() for speed in speeds}
    rgraph = {speed: set() for speed in speeds}
    for (owner, protector), count in counts.items():
        if count <= 0:
            continue
        graph[owner].add(protector)
        rgraph[protector].add(owner)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for speed in speeds:
        if speed not in seen:
            dfs(speed)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: int) -> int:
        seen.add(v)
        total = 1
        for w in rgraph[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for speed in reversed(order):
        if speed not in seen:
            sizes.append(rdfs(speed))
    return tuple(sorted(sizes, reverse=True))


def summarize(label: str, speeds: tuple[int, ...]) -> HandoffSummary:
    base = S360.summarize(list(speeds))
    lifts = endpoint_lifts(base.speeds)
    counts = protection_counts(base.speeds)
    edges, ties = majority_edges(base.speeds, counts)

    owner_out = Counter()
    protector_in = Counter()
    for (owner, protector), count in counts.items():
        if count:
            owner_out[owner] += 1
            protector_in[protector] += 1

    return HandoffSummary(
        label=label,
        n=len(base.speeds) + 1,
        speeds=base.speeds,
        classification=base.classification,
        threshold=base.threshold,
        unique_endpoints=base.unique_endpoint_count,
        unprotected=base.unprotected_count,
        endpoint_state_hist=tuple(sorted(Counter(row.state for row in lifts).items())),
        owner_out_hist=tuple(sorted(Counter(owner_out.values()).items())),
        protector_in_hist=tuple(sorted(Counter(protector_in.values()).items())),
        decisive_pairs=len(edges),
        tied_pairs=ties,
        majority_cycles=cycle_count_3(base.speeds, edges),
        scc_sizes=scc_sizes(base.speeds, counts),
        first_witness=base.first_unprotected,
    )


def two_sided_gap_text(row: EndpointLift) -> str:
    left = "none" if row.left_speed is None else f"{row.left_speed}@{fmt(row.left_distance)}"
    right = "none" if row.right_speed is None else f"{row.right_speed}@{fmt(row.right_distance)}"
    second = fmt(row.second_distance)
    return (
        f"t={fmt(row.value)} state={row.state} owners={row.owners} "
        f"inside={row.inside} on={row.on_boundary} nearest={fmt(row.nearest_distance)} "
        f"second={second} left={left} right={right}"
    )


def print_summary(row: HandoffSummary) -> None:
    print(f"[{row.label}]")
    print(f"  n={row.n} speeds={row.speeds}")
    print(
        f"  classification={row.classification} threshold={fmt(row.threshold)} "
        f"unique_endpoints={row.unique_endpoints} unprotected={row.unprotected} "
        f"first_witness={fmt(row.first_witness)}"
    )
    print(f"  endpoint_state_hist={dict(row.endpoint_state_hist)}")
    print(
        f"  handoff_scc_sizes={row.scc_sizes} "
        f"owner_out_hist={dict(row.owner_out_hist)} "
        f"protector_in_hist={dict(row.protector_in_hist)}"
    )
    print(
        f"  majority_handoff_pairs decisive={row.decisive_pairs} "
        f"ties={row.tied_pairs} directed_3cycles={row.majority_cycles}"
    )
    print()


def print_endpoint_examples(label: str, speeds: tuple[int, ...], limit: int = 4) -> None:
    lifts = endpoint_lifts(tuple(sorted(speeds)))
    print(f"{label}: two-neighbor endpoint samples")
    print("-" * 96)
    witnesses = [row for row in lifts if row.state == "witness"]
    singles = [row for row in lifts if row.state == "single_handoff"]
    redundant = [row for row in lifts if row.state == "redundant_handoff"]
    for family_name, family in (
        ("witness", witnesses[:limit]),
        ("single_handoff", singles[:limit]),
        ("redundant_handoff", redundant[:limit]),
    ):
        if not family:
            continue
        print(f"  {family_name}")
        for row in family:
            print(f"    {two_sided_gap_text(row)}")
    print()


def semicircle_scores(speeds: tuple[int, ...], t: Fraction) -> tuple[tuple[int, int], ...]:
    labels = (0,) + tuple(sorted(speeds))
    positions = {label: circle(label * t) for label in labels}
    scores: dict[int, int] = {label: 0 for label in labels}
    for a, b in combinations(labels, 2):
        delta = circle(positions[b] - positions[a])
        if delta == Fraction(1, 2):
            continue
        if 0 < delta < Fraction(1, 2):
            scores[a] += 1
        else:
            scores[b] += 1
    return tuple(sorted(scores.items(), key=lambda item: (-item[1], item[0])))


def print_round_tournament_probe(label: str, speeds: tuple[int, ...]) -> None:
    summary = S360.summarize(list(speeds))
    if summary.first_unprotected is None:
        return
    t = summary.first_unprotected
    records = distance_records(tuple(sorted(speeds)), t)
    left_speed, left_distance = side_nearest(records, "left")
    right_speed, right_distance = side_nearest(records, "right")
    print(f"{label}: round-tournament probe at first witness t={fmt(t)}")
    print("-" * 96)
    print(
        f"  stationary two-sided neighbors: "
        f"left={left_speed}@{fmt(left_distance)} right={right_speed}@{fmt(right_distance)}"
    )
    print(f"  top semicircle scores including stationary runner: {semicircle_scores(speeds, t)[:8]}")
    print()


def print_synthesis() -> None:
    print("S461 synthesis")
    print("=" * 96)
    print(
        "Nearest-runner distance is the scalar shadow.  The endpoint proof object "
        "starts when we also remember which runner is second at a boundary and "
        "which side of the stationary runner each blocker occupies."
    )
    print(
        "At an LRC endpoint, a full cover needs a strict inside runner while an "
        "owner sits exactly on the threshold.  That is a handoff, not a pointwise "
        "minimum.  Handoff data naturally forms a labelled directed graph and a "
        "pairwise protection-dominance tournament shadow."
    )
    print(
        "The useful tournament import is therefore not arbitrary tournaments of "
        "speed sets.  It is good-cut/SCC technology for handoff cores, plus "
        "Omega-style conflict graphs for compatible repair packets."
    )
    print(
        "Zeckendorf enters at the next layer: if the handoff core quotients to a "
        "path-like carry graph, adjacent repairs are illegal and the residual "
        "endpoint debt has a no-adjacent normal form."
    )
    print()


def main() -> None:
    examples = [
        ("initial n=14", tuple(range(1, 14))),
        (
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
        ),
        (
            "n14 S380 gate ladder",
            (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
        ),
        ("initial n=16", tuple(range(1, 16))),
    ]

    print("=" * 96)
    print("LRC TOURNAMENT AND TWO-NEIGHBOR LIFT")
    print("=" * 96)
    print("codex-2026-05-31-S461")
    print()

    for label, speeds in examples:
        print_summary(summarize(label, speeds))

    print_endpoint_examples("initial n=14", tuple(range(1, 14)))
    print_endpoint_examples(
        "n14 S380 gate ladder",
        (1, 14, 28, 42, 56, 70, 98, 112, 126, 140, 154, 168, 182),
    )
    print_round_tournament_probe("initial n=14", tuple(range(1, 14)))
    print_round_tournament_probe(
        "n14 seven-ladder",
        (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
    )
    print_synthesis()


if __name__ == "__main__":
    main()
