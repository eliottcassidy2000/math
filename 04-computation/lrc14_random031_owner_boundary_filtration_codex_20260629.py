#!/usr/bin/env python3
"""HYP-3522: random031 owner-boundary bracket filtration.

This script refines HYP-3520's four-owner owner-current boundary word on the
owner side of the HYP-3494/HYP-3510/HYP-3511 random031 route.  The question is
not whether every seam owner is transported through the lower-delta bypass.
The useful finite question is sharper:

    which owners are transport charge,
    which owners are branch-boundary bracket charge,
    and which owners remain true seam-boundary debt?

The output is evidence and experiment design, not an LRC14 proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
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


H3486 = load_module(
    "hyp3486_for_hyp3512",
    "lrc14_random031_seam_complement_fiber_graph_codex_20260629.py",
)
H3510 = load_module(
    "hyp3510_for_hyp3512",
    "lrc14_random031_seam_complement_flow_codex_20260629.py",
)
H3511 = load_module(
    "hyp3511_for_hyp3512",
    "lrc14_random031_free_hole_bracket_atlas_codex_20260629.py",
)


SEAM_OWNERS = H3486.H3481.SEVEN_OWNER_CORE
RESCUE_CORE = H3486.H3481.RESCUE_CORE
HARD_COMPONENTS = H3486.H3481.HARD_COMPONENTS


@dataclass(frozen=True)
class BranchBoundary:
    branch: int
    bypass_u: tuple[int, ...]
    bypass_phases: tuple[int, ...]
    left_cell: object
    right_cell: object


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("transport_owner_word_constant", 98),
    Carrier("mirror_pair_owner_word_persistence", 94),
    Carrier("branch_boundary_bracket_lift", 87),
    Carrier("residual_puncture_apex_debt", 82),
    Carrier("free_hole_bracket_separation", 76),
    Carrier("coarse_connected_phase_carrier", 70),
    Carrier("raw_seven_owner_shadow", 18),
)


def fmt(value) -> str:
    if value is None:
        return "None"
    if isinstance(value, int):
        return str(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_text(interval) -> str:
    return f"[{fmt(interval[0])},{fmt(interval[1])}]"


def owner_union_from_hits(hits) -> tuple[int, ...]:
    return tuple(sorted({owner for hit in hits for owner in hit.owners}))


def owner_union_from_cells(cells) -> tuple[int, ...]:
    return tuple(sorted({owner for cell in cells for hit in cell.hits for owner in hit.owners}))


def gate_owner_union(gates) -> tuple[int, ...]:
    return tuple(sorted({owner for gate in gates for owner in H3486.H3481.cover_owners(gate)}))


def cover_speeds(cover) -> tuple[int, ...]:
    return tuple() if cover is None else tuple(cover[1])


def component_owner_pair(component) -> tuple[int, ...]:
    return tuple(sorted(cover_speeds(component.b0_cover) + cover_speeds(component.b1_cover)))


def residue_word(owners: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(owner % 14 for owner in owners)


def compact_counter(counter: Counter) -> dict:
    return dict(sorted(counter.items(), key=lambda item: item[0]))


def build_context():
    row = H3486.H3481.H3450.audit_row(H3486.H3481.ROW_NAME, H3486.H3481.SPEEDS)
    gates = H3486.H3481.build_gates()
    cells = H3486.build_cells(gates, row)
    by_node = {cell.node: cell for cell in cells}
    legal = H3486.connected_components(cells, {"horizontal", "mirror"})
    bypass_component = [component for component in legal if component.type_word == "bypass"]
    if len(bypass_component) != 1:
        raise RuntimeError(f"expected one pure bypass component, got {len(bypass_component)}")
    bypass_component = bypass_component[0]
    bypass_cells = tuple(by_node[node] for node in bypass_component.nodes)
    hard_seam_gates = [
        gate
        for gate in gates
        if gate.component_index in HARD_COMPONENTS and gate.total_delta == 7
    ]
    lower_bypass_gates = [
        gate
        for gate in gates
        if gate.component_index in HARD_COMPONENTS and gate.total_delta < 7
    ]
    return row, gates, cells, by_node, legal, bypass_component, bypass_cells, hard_seam_gates, lower_bypass_gates


def branch_boundaries(cells, by_node, bypass_component) -> tuple[BranchBoundary, ...]:
    branch_seq = {
        branch: tuple(sorted((cell for cell in cells if cell.branch == branch), key=lambda item: item.u_index))
        for branch in (0, 1)
    }
    branch_pos = {
        (cell.branch, cell.u_index): index
        for branch in (0, 1)
        for index, cell in enumerate(branch_seq[branch])
    }
    boundaries: list[BranchBoundary] = []
    for branch in (0, 1):
        bypass_u = tuple(sorted(u for (u, b) in bypass_component.nodes if b == branch))
        positions = sorted(branch_pos[(branch, u)] for u in bypass_u)
        left = branch_seq[branch][positions[0] - 1]
        right = branch_seq[branch][(positions[-1] + 1) % len(branch_seq[branch])]
        boundaries.append(
            BranchBoundary(
                branch=branch,
                bypass_u=bypass_u,
                bypass_phases=tuple(by_node[(u, branch)].phase for u in bypass_u),
                left_cell=left,
                right_cell=right,
            )
        )
    return tuple(boundaries)


def mirror_pairs(by_node, bypass_component):
    pairs = []
    node_set = set(bypass_component.nodes)
    seen: set[tuple[int, int]] = set()
    for node in sorted(node_set):
        mate = H3486.mirror_node(node)
        if node in seen or mate not in node_set:
            continue
        seen.add(node)
        seen.add(mate)
        pairs.append((by_node[node], by_node[mate]))
    return tuple(pairs)


def free_hole_bracket_stats():
    packets, node_to_packet, branch_seq = H3511.packet_summaries()
    clusters = H3511.half_open_clusters(packets, node_to_packet, branch_seq)
    return {
        "packet_count": len(packets),
        "packet_size_hist": Counter(packet.size for packet in packets),
        "bracketed_count": sum(1 for packet in packets if not packet.half_open),
        "half_open_count": sum(1 for packet in packets if packet.half_open),
        "half_open_cluster_count": len(clusters),
        "half_open_cluster_packet_size_hist": Counter(len(cluster) for cluster in clusters),
        "half_open_cluster_cell_size_hist": Counter(sum(packets[idx].size for idx in cluster) for cluster in clusters),
        "clusters": tuple(tuple(cluster) for cluster in clusters),
    }


def seam_connectivity_stats():
    data = H3510.build_graph()
    pure_branch, incidence, branch, mirror = H3510.graph_stats(data)
    return {
        "pure_branch_cycle_sizes": tuple(sorted(len(comp) for comp in pure_branch)),
        "incidence_components": len(incidence),
        "branch_order_components": len(branch),
        "mirror_completed_components": len(mirror),
        "branch_summaries": tuple(H3510.component_summary(comp) for comp in branch),
        "mirror_summaries": tuple(H3510.component_summary(comp) for comp in mirror),
    }


def owner_layer(owner, transport, branch_boundary, dead_owners, outer_rescue, apex):
    in_transport = owner in transport
    in_boundary = owner in branch_boundary
    if in_transport and in_boundary:
        return "transport+branch_boundary"
    if in_transport:
        return "transport_only"
    if in_boundary:
        return "branch_boundary_only"
    if owner in apex:
        return "residual_apex_boundary"
    if owner in dead_owners:
        return "residual_dead_island_boundary"
    if owner in outer_rescue:
        return "residual_outer_rescue_boundary"
    return "residual_boundary"


def main() -> None:
    (
        row,
        gates,
        cells,
        by_node,
        legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        lower_bypass_gates,
    ) = build_context()
    boundaries = branch_boundaries(cells, by_node, bypass_component)
    mirror = mirror_pairs(by_node, bypass_component)
    free_stats = free_hole_bracket_stats()
    conn_stats = seam_connectivity_stats()

    seam_owners = gate_owner_union(hard_seam_gates)
    transport_owners = owner_union_from_cells(bypass_cells)
    lower_bypass_owners = gate_owner_union(lower_bypass_gates)
    boundary_cells = tuple(cell for boundary in boundaries for cell in (boundary.left_cell, boundary.right_cell))
    branch_boundary_owners = owner_union_from_cells(boundary_cells)
    transport_plus_boundary = tuple(sorted(set(transport_owners) | set(branch_boundary_owners)))
    residual_after_transport = tuple(owner for owner in seam_owners if owner not in transport_owners)
    residual_after_branch_boundary = tuple(owner for owner in seam_owners if owner not in transport_plus_boundary)
    bracket_lift = tuple(owner for owner in branch_boundary_owners if owner not in transport_owners)
    transport_only = tuple(owner for owner in transport_owners if owner not in branch_boundary_owners)
    dead_owners = tuple(sorted({owner for component in row.dead_components for owner in component_owner_pair(component)}))
    outer_rescue = tuple(owner for owner in RESCUE_CORE if owner not in dead_owners)
    apex = tuple(owner for owner in seam_owners if owner not in RESCUE_CORE)

    bypass_owner_words = Counter(owner_union_from_hits(cell.hits) for cell in bypass_cells)
    mirror_owner_persistent = all(
        owner_union_from_hits(left.hits) == owner_union_from_hits(right.hits)
        for left, right in mirror
    )
    hard_gate_hits = sum(1 for cell in cells if cell.cell_class == "hard_seam")

    print("HYP-3522 RANDOM031 OWNER-BOUNDARY BRACKET FILTRATION AUDIT")
    print("status=EVIDENCE / exact owner-filtration certificate; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Owner Filtration")
    print(f"hard_seam_owners={seam_owners}")
    print(f"dead_island_owners={dead_owners}")
    print(f"outer_rescue_owners={outer_rescue}")
    print(f"apex_owner_layer={apex}")
    print(f"pure_bypass_transport_owners={transport_owners}")
    print(f"branch_boundary_bracket_owners={branch_boundary_owners}")
    print(f"transport_plus_branch_boundary={transport_plus_boundary}")
    print(f"residual_after_transport={residual_after_transport}")
    print(f"bracket_lift_owners={bracket_lift}")
    print(f"transport_only_owners={transport_only}")
    print(f"residual_after_branch_boundary={residual_after_branch_boundary}")
    print(
        "reading=the seven-owner seam filters as transport (23,93,113), "
        "branch-boundary lift (147,169), and residual puncture/apex debt (45,173)."
    )
    print()

    print("## Gate-Level Persistence Tests")
    print(f"hard_seam_gate_count={len(hard_seam_gates)}")
    print(f"hard_gate_phase_hits={hard_gate_hits}")
    for gate in hard_seam_gates:
        print(
            "  seam_gate component={component} branch={branch} delta={delta} "
            "interval={interval} owners={owners}".format(
                component=gate.component_index,
                branch=gate.branch_mask,
                delta=gate.total_delta,
                interval=interval_text(gate.interval),
                owners=H3486.H3481.cover_owners(gate),
            )
        )
    print(f"lower_bypass_gate_count={len(lower_bypass_gates)}")
    print(f"lower_bypass_gate_owner_union={lower_bypass_owners}")
    for gate in lower_bypass_gates:
        print(
            "  bypass_gate component={component} branch={branch} delta={delta} "
            "interval={interval} owners={owners}".format(
                component=gate.component_index,
                branch=gate.branch_mask,
                delta=gate.total_delta,
                interval=interval_text(gate.interval),
                owners=H3486.H3481.cover_owners(gate),
            )
        )
    print(f"lower_bypass_owners_equal_transport={lower_bypass_owners == transport_owners}")
    print(f"seam_owners_avoid_transport_debt={all(owner not in transport_owners for owner in residual_after_transport)}")
    print()

    print("## Pure Bypass Stalk")
    print(f"component_size={bypass_component.size}")
    print(f"branch_hist={compact_counter(bypass_component.branch_hist)}")
    print(f"hit_components={bypass_component.hit_components}")
    print(f"endpoint_ranks={bypass_component.endpoint_ranks}")
    print(f"bypass_owner_word_hist={compact_counter(bypass_owner_words)}")
    print(f"mirror_pair_count={len(mirror)}")
    print(f"mirror_owner_word_persistent={mirror_owner_persistent}")
    for left, right in mirror:
        print(
            "  mirror_pair left=(a={la},u={lu},b={lb},p={lp}) "
            "right=(a={ra},u={ru},b={rb},p={rp}) owners={owners}".format(
                la=left.a,
                lu=left.u_index,
                lb=left.branch,
                lp=left.phase,
                ra=right.a,
                ru=right.u_index,
                rb=right.branch,
                rp=right.phase,
                owners=owner_union_from_hits(left.hits),
            )
        )
    print()

    print("## Branch-Order Boundary Brackets")
    for boundary in boundaries:
        left_owners = owner_union_from_hits(boundary.left_cell.hits)
        right_owners = owner_union_from_hits(boundary.right_cell.hits)
        print(
            "branch={branch} bypass_u={u} phases={phases} "
            "left=(a={la},u={lu},p={lp},class={lc},owners={lo}) "
            "right=(a={ra},u={ru},p={rp},class={rc},owners={ro})".format(
                branch=boundary.branch,
                u=boundary.bypass_u,
                phases=boundary.bypass_phases,
                la=boundary.left_cell.a,
                lu=boundary.left_cell.u_index,
                lp=boundary.left_cell.phase,
                lc=boundary.left_cell.cell_class,
                lo=left_owners,
                ra=boundary.right_cell.a,
                ru=boundary.right_cell.u_index,
                rp=boundary.right_cell.phase,
                rc=boundary.right_cell.cell_class,
                ro=right_owners,
            )
        )
    boundary_intersections = {}
    for boundary in boundaries:
        left_owners = set(owner_union_from_hits(boundary.left_cell.hits))
        right_owners = set(owner_union_from_hits(boundary.right_cell.hits))
        boundary_intersections[boundary.branch] = tuple(sorted(left_owners & right_owners))
    print(f"branch_boundary_owner_intersection_by_branch={boundary_intersections}")
    for branch, intersection in boundary_intersections.items():
        print(f"  branch={branch} intersection={intersection}")
    print()

    print("## Owner Table")
    print("owner residue layer dead_island outer_rescue apex transport branch_boundary residual_after_boundary")
    layer_hist = Counter()
    for owner in seam_owners:
        layer = owner_layer(owner, transport_owners, branch_boundary_owners, dead_owners, outer_rescue, apex)
        layer_hist[layer] += 1
        print(
            f"{owner:>5} {owner % 14:>7} {layer:<31} "
            f"{owner in dead_owners!s:<11} {owner in outer_rescue!s:<12} "
            f"{owner in apex!s:<5} {owner in transport_owners!s:<9} "
            f"{owner in branch_boundary_owners!s:<15} {owner in residual_after_branch_boundary}"
        )
    print(f"owner_layer_hist={compact_counter(layer_hist)}")
    print()

    print("## HYP-3510 / HYP-3511 Separation")
    print(f"coarse_branch_cycle_sizes={conn_stats['pure_branch_cycle_sizes']}")
    print(f"coarse_branch_order_components={conn_stats['branch_order_components']}")
    print(f"coarse_mirror_completed_components={conn_stats['mirror_completed_components']}")
    print(f"coarse_branch_summaries={conn_stats['branch_summaries']}")
    print(f"free_hole_packets={free_stats['packet_count']}")
    print(f"free_hole_packet_size_hist={compact_counter(free_stats['packet_size_hist'])}")
    print(f"free_hole_bracketed_count={free_stats['bracketed_count']}")
    print(f"free_hole_half_open_count={free_stats['half_open_count']}")
    print(f"free_hole_half_open_clusters={free_stats['clusters']}")
    print(
        "reading=HYP-3510 supplies connected phase transport and HYP-3511 "
        "turns no-gate cells into a local bracket lemma; neither forces the "
        "residual owners (45,173) into the pure bypass transport word."
    )
    print()

    print("## Route / Compression Guardrails")
    print(
        "firewall_route_guardrail=HYP-3513: private-firewall status is pure "
        "for existing colored axes C/F/N/T, but the full HYP-3490 route needs "
        "sidecar R unless route reconstruction is proved."
    )
    print(
        "seam_sheaf_safe_quotients="
        "('flow_class','allowed_exit','owner_union','sheet_pgf_bucket')"
    )
    print(
        "seam_sheaf_unsafe_quotients="
        "('owner_union_size','endpoint_ranks','branch_hist','size','mirror_closed')"
    )
    print(
        "reading=HYP-3522 may compress only after retaining transport, "
        "branch-boundary, residual, and route sidecars."
    )
    print()

    print("## Proof Pull")
    print(
        "P1: Prove a transport-word constancy lemma: every pure-bypass cell and "
        "mirror mate carries exactly owners (23,93,113)."
    )
    print(
        "P2: Prove a branch-boundary bracket lemma: the ordinary cells adjacent "
        "to the bypass add exactly owners (147,169) to the local seam owner "
        "coverage while leaving residual debt (45,173)."
    )
    print(
        "P3: Combine P1/P2 with HYP-3510/HYP-3511: rank-2 routes and free-hole "
        "brackets discharge the phase carrier; the only remaining owner-boundary "
        "debt is the puncture/apex pair (45,173), not the full seven-owner seam."
    )
    print(
        "P4: A terminal random031 proof can now target a two-owner residual "
        "boundary lemma with HYP-3513 route sidecar R, instead of an "
        "undifferentiated seven-owner gluing clause."
    )
    print()

    print("## Tournament Analysis")
    print("vertices=owner-boundary proof carriers, not runners or raw gate arcs")
    print(
        "pairwise_observable=owner charge retained + transport/boundary split + "
        "mirror persistence + quotient safety"
    )
    print("switch=higher retained terminal proof payload; ties use filtration order")
    print(f"score_hist={compact_counter(Counter(carrier.score for carrier in CARRIERS))}")
    print("directed_3cycles=0")
    print(
        "hamiltonian_path="
        + " -> ".join(carrier.name for carrier in sorted(CARRIERS, key=lambda item: -item.score))
    )
    print()

    print("## Assumption Challenge")
    print(
        "Considered vertices: runners, hard gates, lower-delta bypass gates, owner "
        "labels, dead islands, branch-boundary ordinary cells, mirror pairs, "
        "free-hole packets, coarse phase components, and proof obligations."
    )
    print(
        "Chosen vertices are owner-boundary carriers.  This preserves the "
        "random031 terminal predicate while deliberately destroying raw runner "
        "order only after replacing it with transport, bracket, and residual "
        "owner sidecars."
    )


if __name__ == "__main__":
    main()
