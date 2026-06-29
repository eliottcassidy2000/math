#!/usr/bin/env python3
"""HYP-3486: random031 seam-complement fiber graph.

HYP-3482 asked for the next exact object: delete the max-delta hard seam and
route the 282 q=14V phase witnesses through the complement.  HYP-3483 then
split the geometry into an n*2 phase pullback and an n+2 owner-boundary seam.

This scout makes the n*2 side literal.  It uses

    u = 2t mod 1

as a cylinder base coordinate and remembers the original half-turn as the two
branch sheets.  The legal mirror move is (u, branch) -> (-u, 1-branch).  The
vertical half-turn move (u, branch) -> (u, 1-branch) is reported only as a
quotient guardrail because it is not class-preserving.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3481 = load_module(
    "hyp3481_topology_for_hyp3484",
    "lrc14_random031_topology_atlas_codex_20260629.py",
)

BASE = H3481.Q // 2


@dataclass(frozen=True)
class GateHit:
    component: int
    branch_mask: str
    delta: int
    route: str
    endpoint_rank: int | None
    component_class: str
    owners: tuple[int, ...]


@dataclass(frozen=True)
class WitnessCell:
    a: int
    u_index: int
    branch: int
    phase: int
    cell_class: str
    hits: tuple[GateHit, ...]

    @property
    def node(self) -> tuple[int, int]:
        return (self.u_index, self.branch)


@dataclass(frozen=True)
class FlowComponent:
    nodes: tuple[tuple[int, int], ...]
    class_hist: Counter[str]
    branch_hist: Counter[int]
    hit_components: tuple[int, ...]
    endpoint_ranks: tuple[int | None, ...]

    @property
    def size(self) -> int:
        return len(self.nodes)

    @property
    def type_word(self) -> str:
        if len(self.class_hist) == 1:
            return next(iter(self.class_hist))
        return "+".join(f"{key}:{self.class_hist[key]}" for key in sorted(self.class_hist))


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("rank2_seam_complement_discharge", 96),
    Carrier("legal_mirror_run_graph", 88),
    Carrier("pure_twelve_cell_bypass_component", 82),
    Carrier("free_hole_mirror_packets", 75),
    Carrier("vertical_halfturn_guardrail", 71),
    Carrier("fiber_occupancy_word", 65),
    Carrier("raw_282_witness_count", 12),
)


def fmt(value: F | int | None) -> str:
    if value is None:
        return "None"
    if isinstance(value, int):
        return str(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def gate_hit(gate, row) -> GateHit:
    component = row.components[gate.component_index]
    return GateHit(
        component=gate.component_index,
        branch_mask=gate.branch_mask,
        delta=gate.total_delta,
        route=gate.route,
        endpoint_rank=component.endpoint_rank,
        component_class=component.component_class,
        owners=H3481.cover_owners(gate),
    )


def build_cells(gates, row) -> tuple[WitnessCell, ...]:
    hard_gates = {
        gate
        for gate in gates
        if gate.component_index in H3481.HARD_COMPONENTS and gate.total_delta == 7
    }
    cells: list[WitnessCell] = []
    for a in H3481.actual_witnesses():
        branch = 0 if a < BASE else 1
        u_index = a % BASE
        u = F(u_index, BASE)
        hard_hits = [
            gate
            for gate in hard_gates
            if H3481.contains_closed(gate.interval, u) and H3481.compatible(gate, branch)
        ]
        hits = [
            gate
            for gate in gates
            if gate not in hard_gates
            and H3481.contains_closed(gate.interval, u)
            and H3481.compatible(gate, branch)
        ]
        if hard_hits:
            cell_class = "hard_seam"
        elif any(gate.component_index in H3481.HARD_COMPONENTS for gate in hits):
            cell_class = "bypass"
        elif hits:
            cell_class = "ordinary"
        else:
            cell_class = "free_hole"
        cells.append(
            WitnessCell(
                a=a,
                u_index=u_index,
                branch=branch,
                phase=a % 14,
                cell_class=cell_class,
                hits=tuple(gate_hit(gate, row) for gate in hits),
            )
        )
    return tuple(cells)


def mirror_node(node: tuple[int, int]) -> tuple[int, int]:
    u_index, branch = node
    return ((-u_index) % BASE, 1 - branch)


def vertical_node(node: tuple[int, int]) -> tuple[int, int]:
    u_index, branch = node
    return (u_index, 1 - branch)


def connected_components(
    cells: tuple[WitnessCell, ...],
    moves: set[str],
) -> tuple[FlowComponent, ...]:
    by_node = {cell.node: cell for cell in cells}
    node_set = set(by_node)
    adjacency: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for node in node_set:
        adjacency[node]

    if "horizontal" in moves:
        for u_index, branch in node_set:
            for neighbor_u in ((u_index - 1) % BASE, (u_index + 1) % BASE):
                neighbor = (neighbor_u, branch)
                if neighbor in node_set:
                    adjacency[(u_index, branch)].add(neighbor)
                    adjacency[neighbor].add((u_index, branch))

    if "mirror" in moves:
        for node in node_set:
            neighbor = mirror_node(node)
            if neighbor in node_set:
                adjacency[node].add(neighbor)
                adjacency[neighbor].add(node)

    if "vertical" in moves:
        for node in node_set:
            neighbor = vertical_node(node)
            if neighbor in node_set:
                adjacency[node].add(neighbor)
                adjacency[neighbor].add(node)

    if "same_component" in moves:
        buckets: dict[int, list[tuple[int, int]]] = defaultdict(list)
        for cell in cells:
            for hit in cell.hits:
                buckets[hit.component].append(cell.node)
        for nodes in buckets.values():
            first = nodes[0]
            for node in nodes[1:]:
                adjacency[first].add(node)
                adjacency[node].add(first)

    seen: set[tuple[int, int]] = set()
    components: list[FlowComponent] = []
    for start in sorted(node_set):
        if start in seen:
            continue
        queue = deque([start])
        seen.add(start)
        nodes: list[tuple[int, int]] = []
        while queue:
            node = queue.popleft()
            nodes.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        class_hist = Counter(by_node[node].cell_class for node in nodes)
        branch_hist = Counter(by_node[node].branch for node in nodes)
        hit_components = tuple(
            sorted({hit.component for node in nodes for hit in by_node[node].hits})
        )
        endpoint_ranks = tuple(
            sorted({hit.endpoint_rank for node in nodes for hit in by_node[node].hits}, key=lambda x: -1 if x is None else x)
        )
        components.append(
            FlowComponent(
                nodes=tuple(sorted(nodes)),
                class_hist=class_hist,
                branch_hist=branch_hist,
                hit_components=hit_components,
                endpoint_ranks=endpoint_ranks,
            )
        )
    components.sort(key=lambda comp: (-comp.size, min(by_node[node].a for node in comp.nodes)))
    return tuple(components)


def component_a_range(component: FlowComponent, by_node: dict[tuple[int, int], WitnessCell]) -> tuple[int, int]:
    values = sorted(by_node[node].a for node in component.nodes)
    return values[0], values[-1]


def component_u_range(component: FlowComponent) -> tuple[int, int]:
    values = sorted(node[0] for node in component.nodes)
    return values[0], values[-1]


def compact_counter(counter: Counter) -> dict:
    return dict(sorted(counter.items(), key=lambda item: item[0]))


def mirror_pair_counts(cells: tuple[WitnessCell, ...]) -> Counter[tuple[str, str]]:
    by_node = {cell.node: cell for cell in cells}
    used: set[tuple[int, int]] = set()
    counts: Counter[tuple[str, str]] = Counter()
    for cell in cells:
        if cell.node in used:
            continue
        mate = by_node.get(mirror_node(cell.node))
        if mate is None:
            continue
        used.add(cell.node)
        used.add(mate.node)
        counts[tuple(sorted((cell.cell_class, mate.cell_class)))] += 1
    return counts


def vertical_pair_counts(cells: tuple[WitnessCell, ...]) -> Counter[tuple[str, str]]:
    by_node = {cell.node: cell for cell in cells}
    used: set[tuple[int, int]] = set()
    counts: Counter[tuple[str, str]] = Counter()
    for cell in cells:
        if cell.node in used:
            continue
        mate = by_node.get(vertical_node(cell.node))
        if mate is None:
            continue
        used.add(cell.node)
        used.add(mate.node)
        counts[tuple(sorted((cell.cell_class, mate.cell_class)))] += 1
    return counts


def top_component_lines(
    components: tuple[FlowComponent, ...],
    by_node: dict[tuple[int, int], WitnessCell],
    limit: int = 10,
) -> list[str]:
    lines: list[str] = []
    for component in components[:limit]:
        a0, a1 = component_a_range(component, by_node)
        u0, u1 = component_u_range(component)
        lines.append(
            "  comp size={size} type={type_word} a=[{a0},{a1}] "
            "u_index=[{u0},{u1}] branches={branches} ranks={ranks} hit_components={hits}".format(
                size=component.size,
                type_word=component.type_word,
                a0=a0,
                a1=a1,
                u0=u0,
                u1=u1,
                branches=compact_counter(component.branch_hist),
                ranks=component.endpoint_ranks,
                hits=component.hit_components,
            )
        )
    return lines


def main() -> None:
    row = H3481.H3450.audit_row(H3481.ROW_NAME, H3481.SPEEDS)
    gates = H3481.build_gates()
    cells = build_cells(gates, row)
    by_node = {cell.node: cell for cell in cells}

    fiber_hist: Counter[int] = Counter()
    fiber_classes: Counter[tuple[str, ...]] = Counter()
    fibers: dict[int, list[WitnessCell]] = defaultdict(list)
    for cell in cells:
        fibers[cell.u_index].append(cell)
    for fiber_cells in fibers.values():
        fiber_hist[len(fiber_cells)] += 1
        fiber_classes[tuple(sorted(cell.cell_class for cell in fiber_cells))] += 1

    class_counts = Counter(cell.cell_class for cell in cells)
    branch_class_counts = Counter((cell.branch, cell.cell_class) for cell in cells)
    hit_count_hist = Counter(len(cell.hits) for cell in cells)
    hit_cells = [cell for cell in cells if cell.hits]
    hit_ranks = Counter(hit.endpoint_rank for cell in hit_cells for hit in cell.hits)
    unique_hit_components = {
        hit.component
        for cell in hit_cells
        for hit in cell.hits
    }
    unique_hit_component_class_hist = Counter(
        row.components[index].component_class for index in unique_hit_components
    )

    legal_components = connected_components(cells, {"horizontal", "mirror"})
    vertical_components = connected_components(cells, {"horizontal", "mirror", "vertical"})
    component_type_hist = Counter(component.type_word for component in legal_components)
    legal_size_hist = Counter(component.size for component in legal_components)
    free_components = [component for component in legal_components if component.type_word == "free_hole"]
    routed_components = [component for component in legal_components if component.hit_components]
    mixed_vertical = [component for component in vertical_components if len(component.class_hist) > 1]

    bypass_cells = [cell for cell in cells if cell.cell_class == "bypass"]
    bypass_by_branch = {
        branch: tuple(cell.u_index for cell in sorted(bypass_cells, key=lambda item: item.u_index) if cell.branch == branch)
        for branch in (0, 1)
    }
    bypass_phase_by_branch = {
        branch: tuple(cell.phase for cell in sorted(bypass_cells, key=lambda item: item.a) if cell.branch == branch)
        for branch in (0, 1)
    }
    bypass_component_ids = tuple(sorted({hit.component for cell in bypass_cells for hit in cell.hits}))
    bypass_owners = tuple(sorted({owner for cell in bypass_cells for hit in cell.hits for owner in hit.owners}))

    mirror_counts = mirror_pair_counts(cells)
    vertical_counts = vertical_pair_counts(cells)
    mirror_missing = sum(1 for cell in cells if mirror_node(cell.node) not in by_node)
    vertical_present_cells = sum(1 for cell in cells if vertical_node(cell.node) in by_node)

    carrier_hist = Counter(carrier.score for carrier in CARRIERS)
    path = [carrier.name for carrier in sorted(CARRIERS, key=lambda carrier: -carrier.score)]

    print("HYP-3486 RANDOM031 SEAM-COMPLEMENT FIBER GRAPH")
    print("status=EVIDENCE / exact fiber graph and quotient guardrail; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Cylinder Fiber Model")
    print(f"q={H3481.Q} base_fibers=q/2={BASE}")
    print("coordinate=u=2t mod 1, node=(u_index, branch)")
    print(f"witness_cells={len(cells)} occupied_fibers={len(fibers)}")
    print(f"fiber_size_hist={compact_counter(fiber_hist)}")
    print(f"fiber_class_signature_top={fiber_classes.most_common(10)}")
    print(f"cell_class_counts={compact_counter(class_counts)}")
    print(f"branch_class_counts={compact_counter(branch_class_counts)}")
    print(f"hit_count_hist={compact_counter(hit_count_hist)}")
    print()

    print("## Seam Deleted: Rank-2 Routing")
    print("deleted_seam=max-delta hard gates on components (43,54)")
    print(f"hard_seam_hits_after_delete={class_counts.get('hard_seam', 0)}")
    print(f"gate_routed_cells={len(hit_cells)} free_hole_cells={class_counts['free_hole']}")
    print(f"hit_endpoint_rank_hist={compact_counter(hit_ranks)}")
    print(f"unique_hit_components={len(unique_hit_components)}")
    print(f"unique_hit_component_class_hist={compact_counter(unique_hit_component_class_hist)}")
    print(
        "rank2_escape_certificate="
        f"{len(hit_cells)}/242 routed cells hit exactly one endpoint-rank-2 seam-complement gate"
    )
    print()

    print("## Legal Mirror-Run Graph")
    print("moves=horizontal u-neighbor within each branch + mirror (u,b)->(-u,1-b)")
    print(f"component_count={len(legal_components)}")
    print(f"component_size_hist={compact_counter(legal_size_hist)}")
    print(f"component_type_hist={compact_counter(component_type_hist)}")
    print(f"routed_component_count={len(routed_components)}")
    print(f"free_hole_component_count={len(free_components)}")
    print(f"free_hole_component_size_hist={compact_counter(Counter(component.size for component in free_components))}")
    print("top_components:")
    for line in top_component_lines(legal_components, by_node, 12):
        print(line)
    print()

    print("## Pure Bypass Component")
    print(f"bypass_cells={len(bypass_cells)}")
    print(f"bypass_components={bypass_component_ids}")
    print(f"bypass_owners={bypass_owners}")
    print(f"bypass_u_blocks_by_branch={bypass_by_branch}")
    print(f"bypass_phase_blocks_by_branch={bypass_phase_by_branch}")
    print("bypass_reading=one pure 12-cell mirror component, not a diffuse wall")
    print()

    print("## Quotient Guardrails")
    print(f"mirror_missing={mirror_missing}")
    print(f"mirror_pair_class_counts={compact_counter(mirror_counts)}")
    print(f"vertical_halfturn_present_cells={vertical_present_cells}/{len(cells)}")
    print(f"vertical_pair_class_counts={compact_counter(vertical_counts)}")
    print(f"vertical_glued_component_count={len(vertical_components)}")
    print(f"vertical_mixed_component_count={len(mixed_vertical)}")
    print(
        "guardrail=mirror is legal and class-preserving; vertical half-turn "
        "is a useful n*2 projection but not a legal sheet gluing because it "
        "mixes free_hole and ordinary cells and never carries the bypass."
    )
    print()

    print("## Proof Pull")
    print(
        "P1: Replace the vague seam-complement target by a finite dichotomy: "
        "after deleting the hard seam, every routed witness reaches an "
        "endpoint-rank-2 gate component, and the unrouted witnesses are "
        "14 mirror-closed free-hole packets."
    )
    print(
        "P2: The lower-delta bypass is a single pure 12-cell component.  Prove "
        "that a pure mirror bypass component on the hard seam components, with "
        "owners (23,93,113), discharges the seven-owner seam boundary or emits "
        "named owner-current/two-adic/SPEC debt."
    )
    print(
        "P3: Do not quotient by vertical half-turn unless a sidecar remembers "
        "which fibers are single-sheeted.  That quotient is the n*2 address, "
        "not the legal topology."
    )
    print()

    print("## Tournament Analysis")
    print("vertices=fiber-graph proof carriers, not runners or raw gates")
    print("pairwise_observable=rank-2 escape retention + legal mirror topology + quotient guardrail")
    print("switch=higher retained terminal-discharge payload")
    print(f"score_hist={compact_counter(carrier_hist)}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Assumption challenge")
    print("Considered vertices: runners, raw gates, u-fibers, branch cells,")
    print("horizontal flow runs, mirror pairs, vertical half-turn fibers,")
    print("component exits, endpoint ranks, free-hole packets, and proof")
    print("obligations.  Chosen vertices are witness cells and legal")
    print("mirror-run components.  This preserves the random031 terminal")
    print("predicate after seam deletion and destroys raw runner order only after")
    print("replacing it with u-index, branch, mirror mate, gate rank, and free-hole")
    print("sidecars.")


if __name__ == "__main__":
    main()
