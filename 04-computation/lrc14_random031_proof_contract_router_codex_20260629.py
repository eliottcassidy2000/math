#!/usr/bin/env python3
"""HYP-3527: random031 proof-contract router.

This script turns the latest random031 stack into a finite proof contract.
It is not an LRC14 proof.  It is a theorem-interface audit: each terminal
clause states which certificate it consumes, which sidecars it must keep, what
tail it emits or retains, and which quotient shortcuts are forbidden.

Tournament Analysis declaration:
  vertices: proof contracts and sidecar obligations, not runners, arcs, raw
            witnesses, chamber volumes, or owner counts;
  pairwise observable: formal readiness plus hidden-tail and scalar-forgetting
            risk;
  switch/gauge: higher proof-readiness score, with dependency order as tie
            path;
  tie Hamiltonian path: ordinary/free-hole clauses before route/owner/vertical
            holdbacks, with residual tail last.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge" / "results"
ROUTE_DISPATCH_RESULT = RESULTS / "lrc14_spigot_route_sidecar_dispatch_codex_20260629.out"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3521 = load_module(
    "hyp3521_for_hyp3527",
    "lrc14_random031_terminal_certificate_ledger_codex_20260629.py",
)
H3524 = load_module(
    "hyp3524_for_hyp3527",
    "lrc14_random031_spigot_hydrotope_scout_codex_20260629.py",
)
H3525 = load_module(
    "hyp3525_for_hyp3527",
    "lrc14_spigot_guarded_emission_atlas_codex_20260629.py",
)

@dataclass(frozen=True)
class ContractClause:
    name: str
    status: str
    consumes: tuple[str, ...]
    emits: tuple[str, ...]
    remaining_tail: tuple[str, ...]
    required_sidecars: tuple[str, ...]
    forbidden_quotients: tuple[str, ...]
    evidence: tuple[str, ...]
    theorem_shape: str
    proof_debt: str
    score: int


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def terminal_contract_facts() -> dict[str, object]:
    certificates = H3521.build_certificates()
    cell_hist = Counter()
    component_hist = Counter()
    for cert in certificates:
        cell_hist[cert.terminal_class] += cert.cell_count
        component_hist[cert.terminal_class] += 1
    return {
        "cell_hist": cell_hist,
        "component_hist": component_hist,
        "terminal_certificate_count": H3521.terminal_certificate_count(certificates),
        "certificates": certificates,
    }


def vertical_guard_facts() -> dict[str, object]:
    _row, cells, _components, by_node = H3521.build_cells_and_components()
    vertical_components = H3521.H3486.connected_components(cells, {"horizontal", "mirror", "vertical"})
    mixed_vertical = [component for component in vertical_components if len(component.class_hist) > 1]
    vertical_present_cells = sum(
        1 for cell in cells if H3521.H3486.vertical_node(cell.node) in by_node
    )
    return {
        "vertical_present_cells": vertical_present_cells,
        "witness_cells": len(cells),
        "vertical_pair_class_counts": H3521.H3486.vertical_pair_counts(cells),
        "vertical_glued_component_count": len(vertical_components),
        "vertical_mixed_component_count": len(mixed_vertical),
    }


def parse_assignment(lines: tuple[str, ...], key: str) -> bool:
    prefix = f"{key}="
    for line in lines:
        if line.startswith(prefix):
            value = line.removeprefix(prefix)
            if value == "True":
                return True
            if value == "False":
                return False
            raise ValueError(f"unexpected boolean for {key}: {value!r}")
    raise KeyError(f"missing {key} in {ROUTE_DISPATCH_RESULT}")


def parse_purity_line(lines: tuple[str, ...], axes: str, target: str) -> dict[str, object]:
    pattern = re.compile(
        rf"^axes={re.escape(axes)} target={re.escape(target)} "
        r"fibers=(\d+) max_fiber=(\d+) mixed=(\d+)/(\d+) pure=(True|False)$"
    )
    for line in lines:
        match = pattern.match(line)
        if match:
            return {
                "fibers": int(match.group(1)),
                "max_fiber": int(match.group(2)),
                "mixed_fibers": int(match.group(3)),
                "mixed_rows": int(match.group(4)),
                "pure": match.group(5) == "True",
            }
    raise KeyError(f"missing axes={axes} target={target} in {ROUTE_DISPATCH_RESULT}")


def route_facts() -> dict[str, object]:
    lines = tuple(ROUTE_DISPATCH_RESULT.read_text().splitlines())
    p_iq_route = parse_purity_line(lines, "('I', 'Q')", "h3490_route")
    p_all_iq_route = parse_purity_line(
        lines,
        "('K', 'N', 'T', 'S', 'F', 'C', 'M', 'A', 'I', 'Q')",
        "h3490_route",
    )
    p_r_route = parse_purity_line(lines, "('R',)", "h3490_route")
    return {
        "source": str(ROUTE_DISPATCH_RESULT.relative_to(ROOT)),
        "private_status_pure_by_I": parse_assignment(lines, "private_status_pure_by_I"),
        "private_status_pure_by_Q": parse_assignment(lines, "private_status_pure_by_Q"),
        "route_pure_by_IQ": parse_assignment(lines, "route_pure_by_IQ"),
        "route_mixed_by_IQ": (p_iq_route["mixed_fibers"], p_iq_route["mixed_rows"]),
        "route_pure_by_all_existing_plus_IQ": p_all_iq_route["pure"],
        "route_mixed_by_all_existing_plus_IQ": (
            p_all_iq_route["mixed_fibers"],
            p_all_iq_route["mixed_rows"],
        ),
        "route_pure_by_R": p_r_route["pure"],
    }


def owner_facts() -> dict[str, object]:
    context = H3524.build_context()
    stages = H3524.spigot_stages(context)
    safety = H3524.stage_safety(stages, context["seam_owners"])
    hydro = H3524.hydrotope_audit(context)
    quotients = H3524.quotient_audit(context, hydro)
    safe_quotients = tuple(item["name"] for item in quotients if item["safe"])
    unsafe_quotients = tuple(item["name"] for item in quotients if not item["safe"])
    residual_bucket_sizes = {
        name: audit["target_rows"]["residual"]["bucket_size"]
        for name, audit in hydro.items()
    }
    return {
        "context": context,
        "stages": stages,
        "safety": safety,
        "safe_quotients": safe_quotients,
        "unsafe_quotients": unsafe_quotients,
        "residual_bucket_sizes": residual_bucket_sizes,
    }


def guarded_emission_facts() -> dict[str, object]:
    scenarios = {
        "route_sidecar": {"R"},
        "safe_sheaf_head": {"flow_class", "allowed_exit", "sheet_pgf_bucket"},
        "owner_filtration_ready": {"transport_word", "branch_boundary_lift", "residual_pair", "R"},
        "raw_count_shadow": {"checksum_only"},
    }
    return {
        name: tuple(H3525.emit_or_hold(visible, target) for target in H3525.EMISSION_TESTS)
        for name, visible in scenarios.items()
    }


def contract_clauses(
    terminal: dict[str, object],
    owner: dict[str, object],
    route: dict[str, object],
    vertical: dict[str, object],
) -> tuple[ContractClause, ...]:
    cell_hist = terminal["cell_hist"]
    component_hist = terminal["component_hist"]
    context = owner["context"]

    ordinary_cells = cell_hist["ordinary_rank2_route_component"]
    ordinary_components = component_hist["ordinary_rank2_route_component"]
    single_cells = cell_hist["free_hole_single_bracket_packet"]
    single_components = component_hist["free_hole_single_bracket_packet"]
    doublet_cells = cell_hist["free_hole_doublet_packet"]
    doublet_components = component_hist["free_hole_doublet_packet"]
    bypass_cells = cell_hist["pure_bypass_owner_boundary_component"]

    residual_tail = tuple(str(owner_id) for owner_id in context["residual_after_boundary"])
    residual_after_transport = tuple(str(owner_id) for owner_id in context["residual_after_transport"])
    transport_emit = tuple(str(owner_id) for owner_id in context["transport_owners"])
    bracket_emit = tuple(str(owner_id) for owner_id in context["bracket_lift"])

    return (
        ContractClause(
            name="ordinary_route_emit",
            status="formal_ready_interface",
            consumes=(f"{ordinary_components} ordinary components", f"{ordinary_cells} cells"),
            emits=("ordinary endpoint-rank-2 route certificate",),
            remaining_tail=(),
            required_sidecars=("endpoint_rank2_route", "route_R", "vertical_halfturn_guard"),
            forbidden_quotients=("raw_242_gate_routed_slogan", "vertical_halfturn_gluing"),
            evidence=("HYP-3521 ordinary terminal split", "HYP-3486 endpoint-rank-2 seam-complement exits"),
            theorem_shape="ordinary_route_emit(component) -> terminal_discharge",
            proof_debt="formal endpoint-rank-2 route lemma",
            score=94,
        ),
        ContractClause(
            name="free_hole_single_emit",
            status="formal_ready_interface",
            consumes=(f"{single_components} single packets", f"{single_cells} cells"),
            emits=("ordinary-bracketed free-hole certificate",),
            remaining_tail=(),
            required_sidecars=("HYP3511_single_bracket", "route_R", "vertical_halfturn_guard"),
            forbidden_quotients=("raw_free_hole_count_40", "puncture_count_shadow"),
            evidence=("HYP-3511 ordinary-bracketed singles", "HYP-3523 emit_after_single_bracket_carry"),
            theorem_shape="free_hole_single_emit(packet) -> terminal_discharge",
            proof_debt="formal HYP-3511 single-bracket lemma",
            score=90,
        ),
        ContractClause(
            name="free_hole_doublet_buffer_emit",
            status="formal_ready_interface",
            consumes=(f"{doublet_components} half-open components", f"{doublet_cells} cells"),
            emits=("2 collapsed doublet certificates",),
            remaining_tail=(),
            required_sidecars=("HYP3511_doublet_wait", "HYP3511_same_branch_doublet", "route_R", "vertical_halfturn_guard"),
            forbidden_quotients=("immediate_first_half_emission", "raw_free_hole_count_40"),
            evidence=("HYP-3523 two held predigits later emitted", "HYP-3511 same-branch doublet clusters"),
            theorem_shape="free_hole_doublet_buffer(first, second) -> terminal_discharge",
            proof_debt="formal doublet buffer/mate lemma",
            score=88,
        ),
        ContractClause(
            name="bypass_transport_emit",
            status="carry_required",
            consumes=(f"{bypass_cells} pure-bypass cells", "seven-owner seam input"),
            emits=transport_emit,
            remaining_tail=residual_after_transport,
            required_sidecars=("HYP3522_transport_word_constancy", "route_R", "owner_word"),
            forbidden_quotients=("bypass_owner_word_only", "raw_bypass_count_12"),
            evidence=("HYP-3522 transport word", "HYP-3524 S1_transport_emitter"),
            theorem_shape="emit_transport(seven_owner_state) -> four_owner_tail_state",
            proof_debt="prove transport-word constancy on pure bypass",
            score=82,
        ),
        ContractClause(
            name="bypass_bracket_lift_emit",
            status="carry_required",
            consumes=("four-owner tail after transport",),
            emits=bracket_emit,
            remaining_tail=residual_tail,
            required_sidecars=("HYP3522_branch_boundary_lift", "route_R", "ordinary_boundary_neighbors"),
            forbidden_quotients=("transport_plus_boundary_without_split", "residue_mod14_chamber_shadow"),
            evidence=("HYP-3522 branch-boundary ordinary neighbors", "HYP-3524 S2_branch_boundary_lift_emitter"),
            theorem_shape="emit_bracket_lift(four_owner_state) -> residual_pair_state",
            proof_debt="prove adjacent ordinary branch-boundary lift",
            score=80,
        ),
        ContractClause(
            name="private_firewall_route_sidecar",
            status="carry_required",
            consumes=("I/Q private incidence cut", "terminal route request"),
            emits=("private firewall status", "route only with R"),
            remaining_tail=("route_R",),
            required_sidecars=("I_or_Q_for_private_status", "R_for_terminal_route"),
            forbidden_quotients=("I/Q_as_route", "all_existing_axes_plus_IQ_as_route"),
            evidence=(
                f"private_status_pure_by_I={route['private_status_pure_by_I']}",
                f"private_status_pure_by_Q={route['private_status_pure_by_Q']}",
                f"route_pure_by_IQ={route['route_pure_by_IQ']}",
                f"route_pure_by_R={route['route_pure_by_R']}",
            ),
            theorem_shape="private_cut(I/Q) plus route_R -> terminal_route_dispatch",
            proof_debt="prove I/Q row-free incidence cut; keep R unless route reconstruction is proved",
            score=78,
        ),
        ContractClause(
            name="vertical_halfturn_guard",
            status="carry_required",
            consumes=("u=2t branch-sheet projection",),
            emits=("address projection only",),
            remaining_tail=("vertical_not_glued",),
            required_sidecars=("vertical_halfturn_not_glued",),
            forbidden_quotients=("vertical_halfturn_sheet_gluing",),
            evidence=(
                f"vertical_present_cells={vertical['vertical_present_cells']}/{vertical['witness_cells']}",
                f"vertical_mixed_component_count={vertical['vertical_mixed_component_count']}",
            ),
            theorem_shape="vertical_halfturn_guard -> quotient_legality",
            proof_debt="formalize vertical halfturn as address projection, not topology quotient",
            score=74,
        ),
        ContractClause(
            name="residual_pair_close_tail",
            status="open_tail_lemma",
            consumes=("residual owners (45,173)",),
            emits=("terminal owner-boundary discharge",),
            remaining_tail=residual_tail,
            required_sidecars=("residual_pair", "route_R", "owner_support_chamber", "no_hidden_tail"),
            forbidden_quotients=(
                "hydrotope_signature_residue_mod14",
                "hydrotope_signature_centered_residue",
                "hydrotope_signature_filtration_layer",
                "sliced_box_volume_*",
                "raw_counts_7_5_2",
            ),
            evidence=(
                "HYP-3524 residual tail isolated as (45,173)",
                f"safe_owner_quotients={owner['safe_quotients']}",
                f"residual_bucket_sizes={owner['residual_bucket_sizes']}",
            ),
            theorem_shape="close_tail(residual_pair_state, route_R) -> terminal_discharge",
            proof_debt="prove two-owner no-hidden-tail residual boundary lemma",
            score=63,
        ),
    )


def main() -> None:
    terminal = terminal_contract_facts()
    owner = owner_facts()
    route = route_facts()
    vertical = vertical_guard_facts()
    guarded = guarded_emission_facts()
    clauses = contract_clauses(terminal, owner, route, vertical)
    ordered = tuple(sorted(clauses, key=lambda clause: (-clause.score, clause.name)))

    print("HYP-3527 RANDOM031 PROOF-CONTRACT ROUTER")
    print("status=EVIDENCE / finite proof-interface contract; not an LRC14 proof")
    print("row=random_covering_031")
    print()
    print("## Integrated Inputs")
    print("terminal_ledger=HYP-3521")
    print("certificate_spigot=HYP-3523")
    print("owner_emitter_and_hydrotope_canary=HYP-3524")
    print("guarded_emission_schema=HYP-3525")
    print("route_sidecar_dispatch=HYP-3526")
    print("owner_filtration=HYP-3522")
    print("owner_boundary_persistence=HYP-3520")
    print("firewall_route_sidecar=HYP-3513")
    print("free_hole_brackets=HYP-3511")
    print("vertical_guardrail=HYP-3486")
    print()
    print("## Exact Terminal / Carry Facts")
    print(f"cell_terminal_hist={compact_counter(terminal['cell_hist'])}")
    print(f"component_terminal_hist={compact_counter(terminal['component_hist'])}")
    print(f"terminal_certificate_count_after_doublet_collapse={terminal['terminal_certificate_count']}")
    print(f"spigot_safety={owner['safety']}")
    print(f"route_facts={route}")
    print(f"vertical_guard_facts={vertical}")
    print()
    print("## Guarded Emission Sanity")
    for name, emissions in guarded.items():
        print(f"{name}: {emissions}")
    print()
    print("## Contract Clauses")
    for clause in ordered:
        print(f"{clause.name}: status={clause.status} score={clause.score}")
        print(f"  consumes={clause.consumes}")
        print(f"  emits={clause.emits}")
        print(f"  remaining_tail={clause.remaining_tail}")
        print(f"  required_sidecars={clause.required_sidecars}")
        print(f"  forbidden_quotients={clause.forbidden_quotients}")
        print(f"  evidence={clause.evidence}")
        print(f"  theorem_shape={clause.theorem_shape}")
        print(f"  proof_debt={clause.proof_debt}")
    print()
    print("## Closure Ledger")
    status_hist = Counter(clause.status for clause in clauses)
    print(f"status_hist={compact_counter(status_hist)}")
    print("formal_ready_count=" + str(status_hist["formal_ready_interface"]))
    print("carry_required_count=" + str(status_hist["carry_required"]))
    print("open_tail_lemma_count=" + str(status_hist["open_tail_lemma"]))
    print("open_tail_lemma=residual_pair_close_tail")
    print("not_closed_reason=the residual pair (45,173) still needs a two-owner no-hidden-tail boundary lemma")
    print()
    print("## Lean Skeleton")
    for clause in clauses:
        print(f"{clause.name}: {clause.theorem_shape}")
    print()
    print("## Tournament Analysis")
    print("vertices=proof contracts and sidecar obligations, not runners, raw witnesses, chamber volumes, or owner counts")
    print("pairwise_observable=formal readiness + hidden-tail risk + scalar-forgetting risk")
    print("switch=higher proof-readiness score; ties use dependency order")
    print(f"score_hist={compact_counter(Counter(clause.score for clause in clauses))}")
    print("directed_3cycles=0")
    print("hamiltonian_path=" + " -> ".join(clause.name for clause in ordered))
    print()
    print("## Assumption Challenge")
    print("Considered vertices: runners, arcs, witness cells, terminal certificates, sidecars, contracts, route fibers, owner tails, hydrotope chambers, quotient cuts, and proof obligations.")
    print("Chosen vertices are proof contracts plus sidecar obligations.  This preserves the random031 terminal predicate by retaining consumed certificate, emitted tail, remaining tail, route R, owner filtration, free-hole carry, and vertical guard fields before any scalar quotient is allowed.")


if __name__ == "__main__":
    main()
