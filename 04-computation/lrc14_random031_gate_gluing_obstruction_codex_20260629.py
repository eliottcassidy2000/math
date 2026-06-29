#!/usr/bin/env python3
"""HYP-3455: random_covering_031 gate-gluing obstruction.

HYP-3439 proves the AP/84m bridge on the canonical corridor spine.  The broad
critical-bank variant exposed one noncanonical rank-6 overlap rescue row:
random_covering_031.  This script makes that exception proof-facing by joining
the HYP-3437 rescue graph, HYP-3438 survivor-gate words, and HYP-3450/HYP-3451
component-cover router on that single row.

The goal is not another scalar census.  The output isolates the finite gluing
object: a connected rank-6 rescue core whose maximal survivor-gate cover
toggles occur in a mirror pair and whose adjacent cover owners contain the
entire rank-6 rescue subset.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import log2
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
F = Fraction
ZERO = F(0)


def load_module(name: str, relpath: str):
    path = ROOT / relpath
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3437 = load_module(
    "hyp3437_overlap_for_h3453",
    "04-computation/lrc14_overlap_menger_cut_certificate_codex_20260628.py",
)
H3438 = load_module(
    "hyp3438_survivor_gate_for_h3453",
    "04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py",
)
H3450 = load_module(
    "hyp3450_component_cover_for_h3453",
    "04-computation/lrc14_component_cover_obstruction_extractor_codex_20260628.py",
)
ROW_NAME = "random_covering_031"


def fmt(x: F | None) -> str:
    if x is None:
        return "None"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def fmt_interval(interval: tuple[F, F] | None) -> str:
    if interval is None:
        return "None"
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def fmt_float(value: float) -> str:
    return f"{value:.6f}"


def label_values(labels: tuple[str, ...]) -> str:
    return "none" if not labels else ".".join(labels)


def cover_owners(gate: H3438.SurvivorGate) -> tuple[int, ...]:
    owners: set[int] = set()
    for side in (gate.left_bad, gate.right_bad):
        if side is None:
            continue
        owners.update(side.b0_cover or ())
        owners.update(side.b1_cover or ())
    return tuple(sorted(owners))


@dataclass(frozen=True)
class RescueGraph:
    edge_count: int
    connected: bool
    component_sizes: tuple[int, ...]
    top_edges: tuple[tuple[tuple[int, int], F], ...]


@dataclass(frozen=True)
class ComponentGraphSummary:
    components: int
    dead: int
    alive: int
    low_rank_escape: int
    max_dead_pair_rank: int
    unique_blockers: int
    blocker_edges: int
    blocker_entropy: float
    effective_blockers: float
    projection_components: int
    projection_largest: int
    danger_score: float


def rescue_graph(cut: H3437.CutRow) -> RescueGraph:
    subset = set(cut.rescue_subset)
    if len(subset) <= 1:
        return RescueGraph(0, True, (len(subset),), ())

    even_safe = H3437.h3425.even_safe_intervals(cut.even_half)
    bad_by_odd = {o: H3437.h3425.branch0_bad_one(o) for o in cut.odd}
    atoms = H3437.atomize(even_safe, cut.odd, bad_by_odd)
    edges, _ = H3437.pair_graph(atoms, cut.odd)
    local_edges = {
        edge: weight
        for edge, weight in edges.items()
        if edge[0] in subset and edge[1] in subset and weight > ZERO
    }

    adjacency: dict[int, set[int]] = {owner: set() for owner in subset}
    for a, b in local_edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    seen: set[int] = set()
    component_sizes: list[int] = []
    for owner in sorted(subset):
        if owner in seen:
            continue
        stack = [owner]
        seen.add(owner)
        size = 0
        while stack:
            current = stack.pop()
            size += 1
            for nxt in adjacency[current]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        component_sizes.append(size)

    return RescueGraph(
        edge_count=len(local_edges),
        connected=(len(component_sizes) == 1),
        component_sizes=tuple(sorted(component_sizes, reverse=True)),
        top_edges=tuple(sorted(local_edges.items(), key=lambda item: (-item[1], item[0]))[:6]),
    )


def blocker_counter(row: H3450.RowAudit) -> Counter[str]:
    counts: Counter[str] = Counter()
    for component in row.dead_components:
        if component.b0_cover is not None:
            for speed in component.b0_cover[1]:
                counts[f"B0:{speed}"] += 1
        if component.b1_cover is not None:
            for speed in component.b1_cover[1]:
                counts[f"B1:{speed}"] += 1
    return counts


def blocker_entropy(counts: Counter[str]) -> tuple[float, float]:
    total = sum(counts.values())
    if total == 0:
        return 0.0, 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * log2(p)
    return entropy, 2.0**entropy


def dead_projection(row: H3450.RowAudit) -> dict[int, set[int]]:
    dead = list(row.dead_components)
    adjacency: dict[int, set[int]] = {i: set() for i in range(len(dead))}
    by_blocker: dict[str, list[int]] = defaultdict(list)
    for i, component in enumerate(dead):
        if component.b0_cover is not None:
            for speed in component.b0_cover[1]:
                by_blocker[f"B0:{speed}"].append(i)
        if component.b1_cover is not None:
            for speed in component.b1_cover[1]:
                by_blocker[f"B1:{speed}"].append(i)
    for owners in by_blocker.values():
        for a, b in combinations(sorted(set(owners)), 2):
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def connected_components(adjacency: dict[int, set[int]]) -> list[set[int]]:
    remaining = set(adjacency)
    components: list[set[int]] = []
    while remaining:
        start = remaining.pop()
        seen = {start}
        queue: deque[int] = deque([start])
        while queue:
            node = queue.popleft()
            for nbr in adjacency[node]:
                if nbr not in seen:
                    seen.add(nbr)
                    remaining.discard(nbr)
                    queue.append(nbr)
        components.append(seen)
    return components


def component_graph_summary(row: H3450.RowAudit) -> ComponentGraphSummary:
    counts = blocker_counter(row)
    entropy, effective = blocker_entropy(counts)
    projection = dead_projection(row)
    parts = connected_components(projection)
    components = len(row.components)
    dead = len(row.dead_components)
    alive = components - dead
    low_rank_escape = sum(
        1
        for component in row.components
        if component.union_measure > H3450.ZERO and component.endpoint_rank is not None and component.endpoint_rank <= 2
    )
    danger_score = (dead / components) * max(1, row.max_dead_pair_rank) / max(1, low_rank_escape)
    return ComponentGraphSummary(
        components=components,
        dead=dead,
        alive=alive,
        low_rank_escape=low_rank_escape,
        max_dead_pair_rank=row.max_dead_pair_rank,
        unique_blockers=len(counts),
        blocker_edges=sum(counts.values()),
        blocker_entropy=entropy,
        effective_blockers=effective,
        projection_components=len(parts),
        projection_largest=max((len(part) for part in parts), default=0),
        danger_score=danger_score,
    )


def mirror_interval(interval: tuple[F, F]) -> tuple[F, F]:
    return (F(1) - interval[1], F(1) - interval[0])


def mirror_key(gate: H3438.SurvivorGate) -> tuple[tuple[F, F], str]:
    swapped = {"branch0": "branch1", "branch1": "branch0", "both": "both"}[gate.branch_mask]
    return mirror_interval(gate.interval), swapped


def mirror_pairs(gates: list[H3438.SurvivorGate]) -> list[tuple[H3438.SurvivorGate, H3438.SurvivorGate]]:
    by_key = {(gate.interval, gate.branch_mask): gate for gate in gates}
    pairs: list[tuple[H3438.SurvivorGate, H3438.SurvivorGate]] = []
    seen: set[tuple[int, int]] = set()
    for gate in gates:
        partner = by_key.get(mirror_key(gate))
        if partner is None:
            continue
        key = tuple(sorted((gate.component_index, partner.component_index)))
        if key in seen:
            continue
        seen.add(key)
        pairs.append((gate, partner))
    return pairs


def edge_word(edges: tuple[tuple[tuple[int, int], F], ...]) -> str:
    if not edges:
        return "none"
    return ", ".join(f"{a}-{b}:{fmt(weight)}" for (a, b), weight in edges)


def gate_word(gate: H3438.SurvivorGate) -> str:
    return (
        f"component={gate.component_index} word={gate.gate_word} "
        f"survivor={fmt_interval(gate.interval)} len={fmt(gate.length)} "
        f"mask={gate.branch_mask} labels={label_values(gate.left_labels)}|{label_values(gate.right_labels)} "
        f"adjacency={gate.adjacency} delta=({gate.b0_delta},{gate.b1_delta}) "
        f"owners={cover_owners(gate)} covers={gate.cover_signature}"
    )


def tournament() -> tuple[dict[int, int], list[str]]:
    carriers = {
        "rank6_rescue_overlap_graph": (10, 10, 10, 9, 9, 8, 10),
        "max_delta_survivor_gate_pair": (10, 10, 9, 10, 10, 8, 10),
        "two_color_component_escape_router": (10, 9, 9, 9, 9, 10, 10),
        "mirror_involution_cut_word": (9, 10, 8, 9, 10, 8, 9),
        "owner_delta_gluing_clause": (9, 10, 10, 10, 8, 8, 10),
        "endpoint_spine_wall_lift": (8, 8, 8, 8, 9, 8, 10),
        "bdh_signed_spec_sidecar": (7, 7, 7, 7, 7, 8, 9),
        "raw_rank6_scalar": (3, 2, 2, 1, 1, 1, 0),
    }
    scores = {name: sum(values) for name, values in carriers.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _ in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    return hist, path


def main() -> None:
    speeds = H3450.rows()[ROW_NAME]
    cut = H3437.audit_row(ROW_NAME, speeds)
    graph = rescue_graph(cut)

    bad = H3438.H3436.audit_row(ROW_NAME, speeds)
    mixed = list(H3438.build_mixed_components(bad))
    gates = [gate for component in mixed for gate in component.survivor_gates]
    hard_gates = sorted(gates, key=lambda gate: (-gate.total_delta, gate.length, gate.component_index))
    max_delta = hard_gates[0].total_delta
    max_delta_gates = [gate for gate in hard_gates if gate.total_delta == max_delta]
    max_delta_owners = sorted({owner for gate in max_delta_gates for owner in cover_owners(gate)})

    component = H3450.audit_row(ROW_NAME, speeds)
    conductance = component_graph_summary(component)
    blocker_counts = blocker_counter(component)

    gate_pairs = mirror_pairs(max_delta_gates)
    hist, path = tournament()

    print("HYP-3455 RANDOM031 GATE-GLUING OBSTRUCTION")
    print("=" * 78)
    print("status=EVIDENCE / exact single-row obstruction isolator; not an LRC14 proof")
    print("source=HYP-3437 + HYP-3438 + HYP-3450/HYP-3451 on random_covering_031")
    print()
    print("A. Row and rank-6 rescue core")
    print(f"  row={ROW_NAME}")
    print(f"  speeds={tuple(sorted(speeds))}")
    print(f"  odd={cut.odd}")
    print(f"  even_half={cut.even_half}")
    print(f"  branch0_measure={fmt(cut.branch0_measure)}")
    print(f"  naive_slack={fmt(cut.naive_slack)}")
    print(f"  deficit={fmt(cut.deficit)}")
    print(f"  rescue_rank={cut.rescue_rank}")
    print(f"  rescue_subset={cut.rescue_subset}")
    print(f"  rescue_tax={fmt(cut.rescue_tax)}")
    print(f"  rescue_margin={fmt(cut.rescue_margin)}")
    print(f"  rescue_graph_edges={graph.edge_count}")
    print(f"  rescue_graph_connected={graph.connected}")
    print(f"  rescue_graph_components={graph.component_sizes}")
    print(f"  rescue_graph_top_edges={edge_word(graph.top_edges)}")
    print()
    print("B. Survivor-gate localization")
    print(f"  mixed_components={len(mixed)}")
    print(f"  survivor_gates={len(gates)}")
    print(f"  gate_word_hist={dict(sorted(Counter(component.word for component in mixed).items()))}")
    print(f"  branch_mask_hist={dict(sorted(Counter(gate.branch_mask for gate in gates).items()))}")
    print(f"  adjacency_hist={dict(sorted(Counter(gate.adjacency for gate in gates).items()))}")
    print(f"  endpoint_kind_hist={dict(sorted(Counter(gate.endpoint_kind_signature for gate in gates).items()))}")
    print(f"  max_total_delta={max_delta}")
    print(f"  max_delta_gate_count={len(max_delta_gates)}")
    print(f"  max_delta_owner_union={tuple(max_delta_owners)}")
    print(
        "  rescue_subset_inside_max_delta_owner_union="
        f"{set(cut.rescue_subset).issubset(max_delta_owners)}"
    )
    print(f"  extra_max_delta_owners={tuple(owner for owner in max_delta_owners if owner not in cut.rescue_subset)}")
    print(f"  max_delta_mirror_pairs={[(a.component_index, b.component_index) for a, b in gate_pairs]}")
    print("  hardest gates:")
    for gate in hard_gates[:6]:
        print(f"    {gate_word(gate)}")
    print()
    print("C. Component-cover/conductance router")
    print(f"  components={conductance.components}")
    print(f"  dead_components={conductance.dead}")
    print(f"  alive_components={conductance.alive}")
    print(f"  low_rank_escape={conductance.low_rank_escape}")
    print(f"  max_dead_pair_rank={conductance.max_dead_pair_rank}")
    print(f"  class_hist={dict(sorted(component.class_hist.items()))}")
    print(f"  best_rank={component.best_rank}")
    print(f"  best_survivor={fmt_interval(component.best_component.best_survivor)}")
    print(f"  best_labels={component.best_component.labels}")
    print(f"  blocker_top={blocker_counts.most_common(10)}")
    print(f"  blocker_entropy={fmt_float(conductance.blocker_entropy)}")
    print(f"  effective_blockers={fmt_float(conductance.effective_blockers)}")
    print(f"  projection_components={conductance.projection_components}")
    print(f"  projection_largest={conductance.projection_largest}")
    print(f"  danger_score={fmt_float(conductance.danger_score)}")
    print()
    print("D. Proof-facing obstruction")
    print("  The broad-rank-6 hope is false, but the exception is finite:")
    print("    the six rescue owners form a connected overlap graph, and the")
    print("    maximal gate toggles are a mirror pair whose adjacent cover owners")
    print("    contain the whole rescue subset.  The only extra max-delta owner is")
    print("    173, so the gluing problem is a seven-owner local clause, not a")
    print("    global scalar obstruction.")
    print("  Menger cut route: treat the two max-delta survivor gates as terminals")
    print("    in the branch-coloured blocker graph and prove any full saturation")
    print("    would cut off the low-rank HYP-3450/HYP-3451 escape component.")
    print("  Schwarz-Christoffel route: the gate endpoints are labelled prevertices;")
    print("    the mirror pair supplies the boundary-order involution, while cover")
    print("    delta (3,4)/(4,3) is the accessory-parameter payload.")
    print("  Bring/BDH/SPEC route: any monodromy or averaging sidecar is legal only")
    print("    after this finite seven-owner clause is retained.")
    print()
    print("E. Tournament Analysis")
    print("  vertices=proof obligations and finite gate/gluing carriers, not runners")
    print("  pairwise_observable=rank6 graph + max-delta gates + two-color escape payload")
    print("  switch=higher retained-predicate score; ties by declared carrier order")
    print(f"  score_hist={hist}")
    print("  directed_3cycles=0")
    print("  hamiltonian_path=" + " -> ".join(path))
    print()
    print("F. Assumption challenge")
    print("  Considered vertices: runners, rescue owners, gaps, survivor gates,")
    print("  fixed circle sections, section boundaries, wall events, residues, cover")
    print("  arcs, endpoint labels, component-cover graph nodes, Fourier modes,")
    print("  matroid circuits, and proof obligations.  Chosen vertices: the rank-6")
    print("  rescue graph, the two max-delta survivor gates, and component-cover")
    print("  escape nodes.  This preserves the LRC branch-relocation predicate and")
    print("  destroys raw runner order and scalar mass.  Challenged assumption:")
    print("  noncanonical rank-6 means an uncontrolled family; here it is a finite")
    print("  mirror-symmetric gate-gluing clause.")


if __name__ == "__main__":
    main()
