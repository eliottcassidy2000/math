#!/usr/bin/env python3
"""HYP-3528: executable random031 proof-contract router.

The recent random031 stack has converged to a contract problem.  HYP-3521
gives terminal certificates, HYP-3523 streams them with bounded carry, HYP-3524
and HYP-3522 split the owner tail, HYP-3525 supplies the guarded-emission rule,
and HYP-3526 says terminal dispatch still needs route sidecar R.

This script turns those facts into a finite proof contract.  Rows are theorem
clauses and required sidecars, not runners, arcs, or scalar counts.  The goal
is to expose exactly what remains before a Lean-facing random031 terminal
packet can be stated cleanly.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
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


H3523 = load_module(
    "hyp3523_for_hyp3527",
    "lrc14_random031_certificate_spigot_codex_20260629.py",
)
H3526 = load_module(
    "hyp3526_for_hyp3527",
    "lrc14_spigot_route_sidecar_dispatch_codex_20260629.py",
)


@dataclass(frozen=True)
class Clause:
    name: str
    status: str
    consumed_certificates: int
    covered_components: int
    covered_cells: int
    required_sidecars: tuple[str, ...]
    emits: tuple[str, ...]
    remaining_tail: tuple[str, ...]
    forbidden_shortcuts: tuple[str, ...]
    theorem_shape: str
    proof_debt: str
    readiness: int
    hidden_tail_risk: int
    scalar_forgetting_risk: int

    @property
    def score(self) -> int:
        return (
            5 * self.readiness
            - 3 * self.hidden_tail_risk
            - 2 * self.scalar_forgetting_risk
            - 2 * len(self.remaining_tail)
            - len(self.required_sidecars)
        )

    @property
    def open_debt(self) -> bool:
        return "OPEN" in self.status or "conditional" in self.status.lower()


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def route_sidecar_facts() -> dict[str, object]:
    records = H3526.H3513.build_records()
    p_i_private = H3526.purity(
        records,
        ("I",),
        "private_firewall_bit",
        lambda record: record.private_firewall,
    )
    p_q_private = H3526.purity(
        records,
        ("Q",),
        "private_firewall_bit",
        lambda record: record.private_firewall,
    )
    p_i_status = H3526.purity(
        records,
        ("I",),
        "private_dead_status",
        lambda record: record.private_dead_status,
    )
    p_q_status = H3526.purity(
        records,
        ("Q",),
        "private_dead_status",
        lambda record: record.private_dead_status,
    )
    p_iq_route = H3526.purity(
        records,
        ("I", "Q"),
        "h3490_route",
        lambda record: record.route,
    )
    p_all_iq_route = H3526.purity(
        records,
        H3526.EXISTING_AXES + ("I", "Q"),
        "h3490_route",
        lambda record: record.route,
    )
    p_r_route = H3526.purity(
        records,
        ("R",),
        "h3490_route",
        lambda record: record.route,
    )
    existing_route_reconstructs, existing_route_axes = H3526.any_route_reconstruction(
        records,
        H3526.EXISTING_AXES,
    )
    return {
        "private_status_by_I": p_i_private.pure and p_i_status.pure,
        "private_status_by_Q": p_q_private.pure and p_q_status.pure,
        "route_by_IQ": p_iq_route.pure,
        "route_by_all_existing_plus_IQ": p_all_iq_route.pure,
        "route_by_R": p_r_route.pure,
        "existing_route_reconstructs": existing_route_reconstructs,
        "existing_route_axes": existing_route_axes,
        "all_existing_plus_IQ_mixed": (p_all_iq_route.mixed_fibers, p_all_iq_route.mixed_rows),
        "all_existing_plus_IQ_fibers": p_all_iq_route.fibers,
        "all_existing_plus_IQ_max_fiber": p_all_iq_route.max_fiber,
    }


def vertical_guardrail_facts() -> dict[str, object]:
    h3486 = H3523.H3521.H3486
    row = h3486.H3481.H3450.audit_row(h3486.H3481.ROW_NAME, h3486.H3481.SPEEDS)
    gates = h3486.H3481.build_gates()
    cells = h3486.build_cells(gates, row)
    by_node = {cell.node: cell for cell in cells}
    vertical_components = h3486.connected_components(cells, {"horizontal", "mirror", "vertical"})
    mixed_vertical = [component for component in vertical_components if len(component.class_hist) > 1]
    vertical_present_cells = sum(1 for cell in cells if h3486.vertical_node(cell.node) in by_node)
    return {
        "vertical_halfturn_present_cells": vertical_present_cells,
        "witness_cells": len(cells),
        "vertical_pair_class_counts": h3486.vertical_pair_counts(cells),
        "vertical_glued_component_count": len(vertical_components),
        "vertical_mixed_component_count": len(mixed_vertical),
    }


def build_contract() -> tuple[tuple[Clause, ...], dict[str, object]]:
    windows, components, certificates, _by_node = H3523.component_windows()
    boundary = H3523.boundary_info(components)
    stream = H3523.simulate_stream(windows, certificates, boundary)
    route = route_sidecar_facts()
    vertical = vertical_guardrail_facts()

    component_terminal_hist = Counter(window.terminal_class for window in windows)
    cell_terminal_hist: Counter[str] = Counter()
    for window in windows:
        cell_terminal_hist[window.terminal_class] += window.cell_count
    action_hist = Counter(record.action_kind for record in stream.records)
    pending_width_hist = Counter(record.pending_doublet_components for record in stream.records)
    owner_width_hist = Counter(len(record.owner_carry) for record in stream.records)

    boundary_ranks = tuple(
        window.rank
        for window in windows
        if window.component_index in boundary.boundary_component_indices
    )
    bypass_rank = next(
        window.rank
        for window in windows
        if window.terminal_class == "pure_bypass_owner_boundary_component"
    )
    doublet_rank_gaps = tuple(
        end_rank - start_rank
        for _cluster, start_rank, end_rank, _members in stream.doublet_events
    )
    owner_rank_gap = (
        stream.owner_events[-1][0] - stream.owner_events[0][0]
        if len(stream.owner_events) >= 2
        else 0
    )

    clauses = (
        Clause(
            name="ordinary_route_emit",
            status="FINITE-CLOSED / formal ordinary-route lemma remains",
            consumed_certificates=component_terminal_hist["ordinary_rank2_route_component"],
            covered_components=component_terminal_hist["ordinary_rank2_route_component"],
            covered_cells=cell_terminal_hist["ordinary_rank2_route_component"],
            required_sidecars=("endpoint_rank2", "seam_complement_route", "R"),
            emits=("ordinary_terminal_certificate",),
            remaining_tail=(),
            forbidden_shortcuts=("raw_230_count", "endpoint_rank_scalar_only"),
            theorem_shape="ordinary_rank2_route_component -> terminal_discharge R",
            proof_debt="formalize the HYP-3486/HYP-3510 rank-2 route discharge",
            readiness=88,
            hidden_tail_risk=5,
            scalar_forgetting_risk=11,
        ),
        Clause(
            name="free_hole_single_emit",
            status="FINITE-CLOSED / formal bracket lemma remains",
            consumed_certificates=component_terminal_hist["free_hole_single_bracket_packet"],
            covered_components=component_terminal_hist["free_hole_single_bracket_packet"],
            covered_cells=cell_terminal_hist["free_hole_single_bracket_packet"],
            required_sidecars=("HYP3511_bracket_packet", "ordinary_boundary_neighbor", "R"),
            emits=("free_hole_single_certificate",),
            remaining_tail=(),
            forbidden_shortcuts=("raw_40_free_hole_count", "no_gate_scalar"),
            theorem_shape="ordinary_bracketed_free_hole_single -> terminal_discharge R",
            proof_debt="formalize the HYP-3511 exposed-boundary bracket packet",
            readiness=84,
            hidden_tail_risk=7,
            scalar_forgetting_risk=13,
        ),
        Clause(
            name="free_hole_doublet_buffer_emit",
            status="FINITE-CLOSED / one-step predigit buffer",
            consumed_certificates=len(stream.doublet_events),
            covered_components=component_terminal_hist["free_hole_doublet_packet"],
            covered_cells=cell_terminal_hist["free_hole_doublet_packet"],
            required_sidecars=("doublet_cluster_id", "one_component_predigit_buffer", "R"),
            emits=("free_hole_doublet_certificate",),
            remaining_tail=(),
            forbidden_shortcuts=("terminal_class_without_spigot_state", "raw_doublet_count"),
            theorem_shape="half_open_same_branch_doublet closes after one component event",
            proof_debt="state the two HYP-3511 doublets as a bounded-buffer lemma",
            readiness=82,
            hidden_tail_risk=8,
            scalar_forgetting_risk=12,
        ),
        Clause(
            name="bypass_transport_emit",
            status="FINITE-CLOSED / transport word constancy lemma remains",
            consumed_certificates=0,
            covered_components=component_terminal_hist["pure_bypass_owner_boundary_component"],
            covered_cells=cell_terminal_hist["pure_bypass_owner_boundary_component"],
            required_sidecars=("seam_owner_word", "bypass_owner_word", "R"),
            emits=tuple(f"owner_{owner}" for owner in boundary.transport_owners),
            remaining_tail=tuple(f"owner_{owner}" for owner in boundary.residual_after_transport),
            forbidden_shortcuts=("bypass_owner_only", "raw_12_bypass_count"),
            theorem_shape="pure_bypass_owner_word emits transport (23,93,113)",
            proof_debt="formalize HYP-3522 transport-word constancy on the pure bypass",
            readiness=76,
            hidden_tail_risk=16,
            scalar_forgetting_risk=14,
        ),
        Clause(
            name="bypass_bracket_lift_emit",
            status="FINITE-CLOSED / adjacent branch-boundary lift lemma remains",
            consumed_certificates=1,
            covered_components=len(boundary.boundary_component_indices),
            covered_cells=0,
            required_sidecars=("branch_boundary_ordinary_neighbors", "bracket_lift_word", "R"),
            emits=tuple(f"owner_{owner}" for owner in boundary.bracket_lift_owners),
            remaining_tail=tuple(f"owner_{owner}" for owner in boundary.residual_after_branch_boundary),
            forbidden_shortcuts=("branch_hist_scalar", "owner_count_shadow"),
            theorem_shape="adjacent ordinary branch-boundary packets emit lift (147,169)",
            proof_debt="formalize the rank-46 branch-boundary bracket lift",
            readiness=73,
            hidden_tail_risk=18,
            scalar_forgetting_risk=15,
        ),
        Clause(
            name="private_firewall_route_sidecar",
            status="GUARDRAIL-CLOSED / route reconstruction remains false",
            consumed_certificates=0,
            covered_components=0,
            covered_cells=0,
            required_sidecars=("I", "Q", "R"),
            emits=("private_firewall_status", "h3490_route_when_R_is_retained"),
            remaining_tail=(),
            forbidden_shortcuts=("I_Q_as_route", "KNTSFCMAIQ_as_route"),
            theorem_shape="I/Q prove private cut; R supplies terminal route",
            proof_debt="prove the row-free multiplicity-one incidence-cut lemma for I/Q",
            readiness=70,
            hidden_tail_risk=14,
            scalar_forgetting_risk=18,
        ),
        Clause(
            name="vertical_halfturn_guard",
            status="GUARDRAIL-CLOSED / forbid n*2 sheet-gluing shortcut",
            consumed_certificates=0,
            covered_components=vertical["vertical_glued_component_count"],
            covered_cells=vertical["witness_cells"],
            required_sidecars=("u_index", "branch_sheet", "mirror_legality"),
            emits=("vertical_quotient_forbidden_cut",),
            remaining_tail=(),
            forbidden_shortcuts=("vertical_halfturn_gluing", "n_times_2_address_projection"),
            theorem_shape="vertical halfturn may address fibers but cannot identify terminal classes",
            proof_debt="carry the HYP-3486 vertical guard into every terminal theorem statement",
            readiness=66,
            hidden_tail_risk=12,
            scalar_forgetting_risk=21,
        ),
        Clause(
            name="residual_pair_close_tail",
            status="OPEN-PROOF-DEBT / final non-streamed random031 tail",
            consumed_certificates=0,
            covered_components=1,
            covered_cells=cell_terminal_hist["pure_bypass_owner_boundary_component"],
            required_sidecars=("residual_pair_(45,173)", "R", "no_hidden_tail_guard"),
            emits=("close_owner_45", "close_owner_173"),
            remaining_tail=("owner_45", "owner_173"),
            forbidden_shortcuts=("residue_mod14", "centered_residue", "sliced_box_volume", "raw_owner_count"),
            theorem_shape="transport + bracket lift + R -> residual (45,173) cannot hide downstream",
            proof_debt="prove the two-owner residual boundary lemma",
            readiness=36,
            hidden_tail_risk=40,
            scalar_forgetting_risk=30,
        ),
        Clause(
            name="raw_count_shadow",
            status="CHECKSUM-ONLY / never theorem-facing",
            consumed_certificates=stream.emitted_certificates,
            covered_components=len(windows),
            covered_cells=sum(window.cell_count for window in windows),
            required_sidecars=(),
            emits=("282", "242", "79", "77", "40", "12"),
            remaining_tail=("all_route_owner_and_vertical_sidecars",),
            forbidden_shortcuts=("count_as_proof_token",),
            theorem_shape="counts verify completed clauses but cannot discharge one",
            proof_debt="none as theorem input; use only after sidecar clauses close",
            readiness=4,
            hidden_tail_risk=45,
            scalar_forgetting_risk=50,
        ),
    )

    facts = {
        "component_tokens": len(windows),
        "terminal_certificate_count": H3523.H3521.terminal_certificate_count(certificates),
        "emitted_certificate_count": stream.emitted_certificates,
        "component_terminal_hist": component_terminal_hist,
        "cell_terminal_hist": cell_terminal_hist,
        "action_hist": action_hist,
        "pending_width_hist": pending_width_hist,
        "owner_width_hist": owner_width_hist,
        "max_pending_doublet_components": stream.max_pending_doublet_components,
        "max_owner_carry_width": stream.max_owner_carry_width,
        "final_owner_carry": stream.final_owner_carry,
        "doublet_rank_gaps": doublet_rank_gaps,
        "owner_rank_gap": owner_rank_gap,
        "boundary_component_ranks": boundary_ranks,
        "bypass_rank": bypass_rank,
        "seam_owners": boundary.seam_owners,
        "transport_owners": boundary.transport_owners,
        "residual_after_transport": boundary.residual_after_transport,
        "bracket_lift_owners": boundary.bracket_lift_owners,
        "residual_after_branch_boundary": boundary.residual_after_branch_boundary,
        "route": route,
        "vertical": vertical,
    }
    return clauses, facts


def beats(left: Clause, right: Clause) -> bool:
    left_key = (
        left.score,
        left.readiness,
        -left.hidden_tail_risk,
        -left.scalar_forgetting_risk,
        left.name,
    )
    right_key = (
        right.score,
        right.readiness,
        -right.hidden_tail_risk,
        -right.scalar_forgetting_risk,
        right.name,
    )
    return left_key > right_key


def directed_three_cycles(clauses: tuple[Clause, ...]) -> int:
    total = 0
    for a, b, c in combinations(clauses, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(b, a) and beats(c, b) and beats(a, c)
        ):
            total += 1
    return total


def main() -> None:
    clauses, facts = build_contract()
    ranked = tuple(sorted(clauses, key=lambda clause: (clause.score, clause.name), reverse=True))
    open_clauses = tuple(clause.name for clause in clauses if clause.open_debt)
    sidecar_hist = Counter(sidecar for clause in clauses for sidecar in clause.required_sidecars)
    forbidden_hist = Counter(shortcut for clause in clauses for shortcut in clause.forbidden_shortcuts)

    print("HYP-3528 RANDOM031 PROOF-CONTRACT ONE-RED-CLAUSE AUDIT")
    print("status=EVIDENCE / executable finite theorem-interface router; not an LRC14 proof")
    print("row=random_covering_031")
    print("source=HYP-3521 + HYP-3523 + HYP-3524 + HYP-3525 + HYP-3526")
    print()

    print("## Contract Thesis")
    print("router=terminal clauses consume sidecars, emit proof tokens, and leave named tail debt")
    print("geometry=random031 is a mirror-punctured cylinder whose forbidden seam is a carry boundary")
    print("phase_flow=the seam complement streams 77 terminal certificates from 79 component events")
    print("only_open_tail=residual_pair_(45,173) plus route sidecar R and no-hidden-tail guard")
    print("namespace_note=two HYP-3525 entries exist upstream; HYP-3528 executes their shared router idea")
    print()

    print("## Imported Finite Stream Facts")
    print(f"component_tokens={facts['component_tokens']}")
    print(f"terminal_certificate_count={facts['terminal_certificate_count']}")
    print(f"emitted_certificate_count={facts['emitted_certificate_count']}")
    print(f"component_terminal_hist={compact_counter(facts['component_terminal_hist'])}")
    print(f"cell_terminal_hist={compact_counter(facts['cell_terminal_hist'])}")
    print(f"action_hist={compact_counter(facts['action_hist'])}")
    print(f"pending_doublet_width_hist={compact_counter(facts['pending_width_hist'])}")
    print(f"owner_carry_width_hist={compact_counter(facts['owner_width_hist'])}")
    print(f"max_pending_doublet_components={facts['max_pending_doublet_components']}")
    print(f"max_owner_carry_width={facts['max_owner_carry_width']}")
    print(f"final_owner_carry={facts['final_owner_carry']}")
    print(f"doublet_rank_gaps={facts['doublet_rank_gaps']}")
    print(f"owner_rank_gap={facts['owner_rank_gap']}")
    print()

    print("## Owner Carry Geometry")
    print(f"seam_owners={facts['seam_owners']}")
    print(f"transport_owners={facts['transport_owners']}")
    print(f"residual_after_transport={facts['residual_after_transport']}")
    print(f"bracket_lift_owners={facts['bracket_lift_owners']}")
    print(f"residual_after_branch_boundary={facts['residual_after_branch_boundary']}")
    print(f"bypass_rank={facts['bypass_rank']}")
    print(f"boundary_component_ranks={facts['boundary_component_ranks']}")
    print("reading=transport opens at rank 45, branch-boundary lift fires at rank 46, residual is (45,173)")
    print()

    print("## Route and Vertical Guardrails")
    route = facts["route"]
    print(f"private_status_by_I={route['private_status_by_I']}")
    print(f"private_status_by_Q={route['private_status_by_Q']}")
    print(f"route_by_IQ={route['route_by_IQ']}")
    print(f"route_by_all_existing_plus_IQ={route['route_by_all_existing_plus_IQ']}")
    print(f"route_by_R={route['route_by_R']}")
    print(f"existing_route_reconstructs={route['existing_route_reconstructs']}")
    print(f"all_existing_plus_IQ_route_mixed={route['all_existing_plus_IQ_mixed']}")
    print(f"all_existing_plus_IQ_fibers={route['all_existing_plus_IQ_fibers']}")
    print(f"all_existing_plus_IQ_max_fiber={route['all_existing_plus_IQ_max_fiber']}")
    vertical = facts["vertical"]
    print(f"vertical_halfturn_present_cells={vertical['vertical_halfturn_present_cells']}/{vertical['witness_cells']}")
    print(f"vertical_pair_class_counts={compact_counter(vertical['vertical_pair_class_counts'])}")
    print(f"vertical_glued_component_count={vertical['vertical_glued_component_count']}")
    print(f"vertical_mixed_component_count={vertical['vertical_mixed_component_count']}")
    print()

    print("## Contract Clauses")
    for clause in ranked:
        print(
            f"clause={clause.name} score={clause.score} status={clause.status} "
            f"certificates={clause.consumed_certificates} components={clause.covered_components} "
            f"cells={clause.covered_cells}"
        )
        print(f"  required_sidecars={clause.required_sidecars}")
        print(f"  emits={clause.emits}")
        print(f"  remaining_tail={clause.remaining_tail}")
        print(f"  forbidden_shortcuts={clause.forbidden_shortcuts}")
        print(f"  theorem_shape={clause.theorem_shape}")
        print(f"  proof_debt={clause.proof_debt}")
    print()

    print("## Router Summary")
    print(f"contract_clause_count={len(clauses)}")
    print(f"open_clause_names={open_clauses}")
    print(f"sidecar_hist={compact_counter(sidecar_hist)}")
    print(f"forbidden_shortcut_count={sum(forbidden_hist.values())}")
    print(f"top_forbidden_shortcuts={forbidden_hist.most_common(10)}")
    print("formalization_order=ordinary_route_emit -> free_hole_single_emit -> free_hole_doublet_buffer_emit -> bypass_transport_emit -> bypass_bracket_lift_emit -> private_firewall_route_sidecar -> vertical_halfturn_guard -> residual_pair_close_tail")
    print("creative_reframe=the proof ABI is not rows-to-routes, but sidecars-to-safe-emission; residual (45,173) is the last untyped-looking tail after all typed carry is routed")
    print()

    print("## Tournament Analysis")
    print("vertices=proof-contract clauses and sidecar obligations, not runners or arcs")
    print("pairwise_observable=readiness - hidden-tail-risk - scalar-forgetting-risk")
    print("switch=higher contract score; ties use readiness, lower risk, and clause name")
    print(f"score_hist={compact_counter(Counter(clause.score for clause in clauses))}")
    print(f"directed_3cycles={directed_three_cycles(clauses)}")
    print("sccs=9 singleton SCCs")
    print("hamiltonian_path=" + " -> ".join(clause.name for clause in ranked))
    print()

    print("## Lean-Facing Skeleton")
    print("def Random031TerminalContract := ordinary_route_emit + free_hole_single_emit + free_hole_doublet_buffer_emit + bypass_transport_emit + bypass_bracket_lift_emit + private_firewall_route_sidecar + vertical_halfturn_guard + residual_pair_close_tail")
    print("lemma random031_contract_streams_all_certificates: component_events=79 -> terminal_certificates=77")
    print("lemma random031_contract_requires_R: I/Q proves private cut but terminal route uses R")
    print("lemma random031_contract_final_tail: after transport and bracket lift, remaining owner tail is exactly (45,173)")
    print("remaining_goal=random031_residual_pair_close_tail")
    print()

    print("## Assumption Challenge")
    print("candidate_vertices=runners,arcs,gaps,fixed_sections,section_boundaries,wall_crossings,residue_buckets,cover_arcs,Fourier_modes,matroid_circuits,stream_states,proof_contracts")
    print("chosen_vertices=proof_contracts_and_sidecar_obligations")
    print("preserves=random031 terminal discharge, owner-tail emission, route R, vertical guard, and forbidden quotient cuts")
    print("destroys=raw row identity, component index after emission, scalar counts, residue chambers, and hydrotope volumes")
    print("challenged_assumption=a proof contract should route rows directly; the better interface routes typed sidecars and only then emits terminal tokens")


if __name__ == "__main__":
    main()
