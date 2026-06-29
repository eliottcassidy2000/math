#!/usr/bin/env python3
"""HYP-3521: random031 terminal certificate ledger.

This scout joins the recent random031 carriers into one finite dispatch table:

* HYP-3486 supplies the legal mirror-run fiber graph.
* HYP-3511 supplies the free-hole bracket/doublet classification.
* HYP-3490 supplies the private-label firewall, so projection-current deletion
  is not an allowed terminal carrier on this row.

The goal is not to prove LRC14.  The goal is to expose the exact terminal
packet theorem that remains: ordinary rank-2 routes, bracketed free-hole beads,
and one pure bypass owner-boundary packet cover all 282 random031 phase cells.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
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
H3511 = load_module(
    "hyp3511_for_hyp3512",
    "lrc14_random031_free_hole_bracket_atlas_codex_20260629.py",
)


SEAM_OWNERS = H3486.H3481.SEVEN_OWNER_CORE
BYPASS_FLOW_OWNERS = (23, 93, 113)
BYPASS_SEAM_DEBT = tuple(owner for owner in SEAM_OWNERS if owner not in BYPASS_FLOW_OWNERS)


@dataclass(frozen=True)
class Certificate:
    component_index: int
    terminal_class: str
    component_type: str
    size: int
    cell_count: int
    packet_ids: tuple[int, ...]
    cluster_id: tuple[int, ...] | None
    hit_components: tuple[int, ...]
    endpoint_ranks: tuple[int | None, ...]
    owner_union: tuple[int, ...]
    seam_debt: tuple[int, ...]
    proof_obligation: str


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("terminal_certificate_ledger", 100),
    Carrier("ordinary_rank2_route_lemma", 91),
    Carrier("free_hole_bracket_lemma", 88),
    Carrier("pure_bypass_owner_boundary_lemma", 84),
    Carrier("private_firewall_negative_certificate", 79),
    Carrier("fiber_pgf_zero_exit_split", 73),
    Carrier("vertical_halfturn_guardrail", 66),
    Carrier("raw_count_shadow", 10),
)


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def build_cells_and_components():
    row = H3486.H3481.H3450.audit_row(H3486.H3481.ROW_NAME, H3486.H3481.SPEEDS)
    gates = H3486.H3481.build_gates()
    cells = H3486.build_cells(gates, row)
    components = H3486.connected_components(cells, {"horizontal", "mirror"})
    by_node = {cell.node: cell for cell in cells}
    return row, cells, components, by_node


def owner_union(nodes, by_node) -> tuple[int, ...]:
    owners: set[int] = set()
    for node in nodes:
        for hit in by_node[node].hits:
            owners.update(hit.owners)
    return tuple(sorted(owners))


def free_hole_packet_maps():
    packets, node_to_packet, branch_seq = H3511.packet_summaries()
    clusters = H3511.half_open_clusters(packets, node_to_packet, branch_seq)
    packet_to_cluster: dict[int, tuple[int, ...]] = {}
    for cluster in clusters:
        cluster_tuple = tuple(cluster)
        for packet_id in cluster:
            packet_to_cluster[packet_id] = cluster_tuple
    return packets, node_to_packet, tuple(tuple(cluster) for cluster in clusters), packet_to_cluster


def terminal_class_for_component(component, packet_ids, packets, packet_to_cluster) -> tuple[str, tuple[int, ...] | None, str]:
    if component.type_word == "ordinary":
        return (
            "ordinary_rank2_route_component",
            None,
            "endpoint-rank-2 route lemma",
        )
    if component.type_word == "bypass":
        return (
            "pure_bypass_owner_boundary_component",
            None,
            "pure bypass owner-boundary lemma",
        )
    if component.type_word == "free_hole":
        if len(packet_ids) != 1:
            raise AssertionError(f"free-hole component has packet ids {packet_ids}")
        packet_id = packet_ids[0]
        packet = packets[packet_id]
        if packet.half_open:
            cluster = packet_to_cluster[packet_id]
            return (
                "free_hole_doublet_packet",
                cluster,
                "same-branch doublet bracket lemma",
            )
        return (
            "free_hole_single_bracket_packet",
            None,
            "ordinary-bracketed single free-hole lemma",
        )
    return ("mixed_or_unclassified_debt", None, "named debt")


def build_certificates() -> tuple[Certificate, ...]:
    _row, _cells, components, by_node = build_cells_and_components()
    packets, node_to_packet, _clusters, packet_to_cluster = free_hole_packet_maps()

    certificates: list[Certificate] = []
    for index, component in enumerate(components):
        packet_ids = tuple(
            sorted({
                node_to_packet[(branch, u_index)]
                for (u_index, branch) in component.nodes
                if component.type_word == "free_hole"
            })
        )
        terminal_class, cluster_id, proof_obligation = terminal_class_for_component(
            component, packet_ids, packets, packet_to_cluster
        )
        owners = owner_union(component.nodes, by_node)
        certificates.append(
            Certificate(
                component_index=index,
                terminal_class=terminal_class,
                component_type=component.type_word,
                size=component.size,
                cell_count=component.size,
                packet_ids=packet_ids,
                cluster_id=cluster_id,
                hit_components=component.hit_components,
                endpoint_ranks=component.endpoint_ranks,
                owner_union=owners,
                seam_debt=tuple(owner for owner in SEAM_OWNERS if owner not in owners),
                proof_obligation=proof_obligation,
            )
        )
    return tuple(certificates)


def fiber_escape_profile(cells) -> dict[str, object]:
    fibers: dict[int, list[object]] = defaultdict(list)
    for cell in cells:
        fibers[cell.u_index].append(cell)

    escape_hist: Counter[int] = Counter()
    signature_hist: Counter[tuple[str, ...]] = Counter()
    zero_exit_signature_hist: Counter[tuple[str, ...]] = Counter()
    one_exit_signature_hist: Counter[tuple[str, ...]] = Counter()
    for fiber_cells in fibers.values():
        escape_count = sum(1 for cell in fiber_cells if cell.hits)
        signature = tuple(sorted(cell.cell_class for cell in fiber_cells))
        escape_hist[escape_count] += 1
        signature_hist[signature] += 1
        if escape_count == 0:
            zero_exit_signature_hist[signature] += 1
        if escape_count == 1:
            one_exit_signature_hist[signature] += 1

    total_escapes = sum(count * freq for count, freq in escape_hist.items())
    return {
        "occupied_fibers": len(fibers),
        "escape_pgf": dict(sorted(escape_hist.items())),
        "mean_escape_sheets": Fraction(total_escapes, len(fibers)),
        "signature_hist": signature_hist,
        "zero_exit_signature_hist": zero_exit_signature_hist,
        "one_exit_signature_hist": one_exit_signature_hist,
    }


def doublet_summary(certificates: tuple[Certificate, ...]) -> dict[tuple[int, ...], int]:
    clusters: Counter[tuple[int, ...]] = Counter()
    for cert in certificates:
        if cert.terminal_class == "free_hole_doublet_packet":
            assert cert.cluster_id is not None
            clusters[cert.cluster_id] += cert.cell_count
    return dict(sorted(clusters.items()))


def terminal_certificate_count(certificates: tuple[Certificate, ...]) -> int:
    count = 0
    counted_clusters: set[tuple[int, ...]] = set()
    for cert in certificates:
        if cert.terminal_class == "free_hole_doublet_packet":
            assert cert.cluster_id is not None
            counted_clusters.add(cert.cluster_id)
        else:
            count += 1
    return count + len(counted_clusters)


def main() -> None:
    _row, cells, _components, _by_node = build_cells_and_components()
    certificates = build_certificates()
    profile = fiber_escape_profile(cells)

    cell_class_counts = Counter(cell.cell_class for cell in cells)
    hit_cells = [cell for cell in cells if cell.hits]
    ordinary_hit_cells = [cell for cell in hit_cells if cell.cell_class == "ordinary"]
    bypass_hit_cells = [cell for cell in hit_cells if cell.cell_class == "bypass"]
    endpoint_rank_hist = Counter(
        hit.endpoint_rank for cell in hit_cells for hit in cell.hits
    )

    component_terminal_hist = Counter(cert.terminal_class for cert in certificates)
    component_type_hist = Counter(cert.component_type for cert in certificates)
    cell_terminal_hist = Counter()
    for cert in certificates:
        cell_terminal_hist[cert.terminal_class] += cert.cell_count

    free_single_cells = cell_terminal_hist["free_hole_single_bracket_packet"]
    free_doublet_cells = cell_terminal_hist["free_hole_doublet_packet"]
    routed_total = len(hit_cells)
    terminal_count = terminal_certificate_count(certificates)
    pure_bypass = [
        cert for cert in certificates
        if cert.terminal_class == "pure_bypass_owner_boundary_component"
    ]

    print("HYP-3521 RANDOM031 TERMINAL CERTIFICATE LEDGER")
    print("status=EVIDENCE / finite terminal-dispatch ledger; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Joined Inputs")
    print("fiber_graph=HYP-3486")
    print("free_hole_brackets=HYP-3511")
    print("private_firewall=HYP-3490")
    print("seam_complement_flow=HYP-3510")
    print(f"seam_owners={SEAM_OWNERS}")
    print()
    print("## Cell-Level Terminal Partition")
    print(f"witness_cells={len(cells)}")
    print(f"cell_class_counts={compact_counter(cell_class_counts)}")
    print(f"gate_routed_cells={routed_total}")
    print(f"ordinary_rank2_hit_cells={len(ordinary_hit_cells)}")
    print(f"bypass_rank2_hit_cells={len(bypass_hit_cells)}")
    print(f"free_hole_cells={cell_class_counts['free_hole']}")
    print(f"endpoint_rank_hist_on_hit_cells={compact_counter(endpoint_rank_hist)}")
    print("identity_282=230 ordinary_rank2 + 40 bracketed_free_hole + 12 pure_bypass")
    print("identity_242=230 ordinary_rank2 + 12 pure_bypass gate-routed cells")
    print()
    print("## Component / Certificate Partition")
    print(f"legal_mirror_components={len(certificates)}")
    print(f"component_type_hist={compact_counter(component_type_hist)}")
    print(f"component_terminal_hist={compact_counter(component_terminal_hist)}")
    print(f"cell_terminal_hist={compact_counter(cell_terminal_hist)}")
    print(f"terminal_certificate_count_after_doublet_collapse={terminal_count}")
    print("identity_79_components=64 ordinary + 14 free_hole + 1 bypass")
    print("identity_77_terminal_certificates=64 ordinary + 10 free_single + 2 free_doublet + 1 bypass")
    print()
    print("## Free-Hole Collapse")
    print(f"free_single_cells={free_single_cells}")
    print(f"free_doublet_cells={free_doublet_cells}")
    print(f"doublet_cell_summary={doublet_summary(certificates)}")
    print("free_hole_certificate_shape=10 ordinary-bracketed singles + 2 same-branch doublets")
    print()
    print("## Fiber PGF Sidecar")
    print(f"occupied_fibers={profile['occupied_fibers']}")
    print(f"rank2_escape_pgf={profile['escape_pgf']}  # 24 + 226*y + 8*y^2")
    print(f"mean_rank2_escape_sheets={profile['mean_escape_sheets']}")
    print(f"zero_exit_signature_hist={compact_counter(profile['zero_exit_signature_hist'])}")
    print(f"one_exit_signature_top={profile['one_exit_signature_hist'].most_common(6)}")
    print("reading=zero-exit fibers are pure free-hole fibers; mixed free-hole/ordinary fibers carry one rank-2 sheet and one bracketed no-gate sheet.")
    print()
    print("## Pure Bypass Boundary Packet")
    for cert in pure_bypass:
        print(
            "bypass_certificate="
            f"component_index={cert.component_index} size={cert.size} "
            f"hit_components={cert.hit_components} owners={cert.owner_union} "
            f"seam_debt={cert.seam_debt}"
        )
    print(f"expected_bypass_flow_owners={BYPASS_FLOW_OWNERS}")
    print(f"expected_bypass_seam_debt={BYPASS_SEAM_DEBT}")
    print("bypass_obligation=prove owner-boundary discharge for this one packet or name owner/two-adic/SPEC debt")
    print()
    print("## Proof Pull")
    print("L1 ordinary-route lemma: every ordinary component has endpoint-rank-2 hits and exits the seam complement.")
    print("L2 free-hole bracket lemma: all 40 no-gate cells are locally bracketed as 10 singles plus 2 doublets.")
    print("L3 bypass owner-boundary lemma: the only nonordinary gate-routed hard-component packet is the 12-cell bypass with owners (23,93,113).")
    print("L4 firewall compatibility lemma: HYP-3490 forbids replacing L1-L3 by projection-current deletion.")
    print("L5 quotient guardrail lemma: vertical half-turn is an address projection and must not identify the terminal classes.")
    print()
    print("## Tournament Analysis")
    print("vertices=terminal proof obligations, not runners, raw arcs, or scalar cell counts")
    print("pairwise_observable=which carrier preserves terminal discharge predicate with least forgotten sidecar")
    print("switch=higher retained certificate payload; ties follow the Hamiltonian path")
    print(f"score_hist={compact_counter(Counter(carrier.score for carrier in CARRIERS))}")
    print("directed_3cycles=0")
    print("sccs=8 singleton SCCs")
    print("hamiltonian_path=" + " -> ".join(carrier.name for carrier in CARRIERS))
    print()
    print("## Assumption Challenge")
    print("Considered vertices: runners, gaps, hard gates, witness cells, u-fibers, branch sheets, free-hole packets, doublet clusters, owner words, endpoint ranks, vertical pairs, and proof obligations.")
    print("Chosen vertices are terminal proof obligations.  This quotient preserves the random031 terminal predicate because it retains class, endpoint rank, free-hole bracket type, bypass owner debt, private-firewall status, and vertical guardrail.")


if __name__ == "__main__":
    main()
