#!/usr/bin/env python3
"""S259 tournament obstruction-transfer atlas.

The motivating pattern is the H=7/H=21 contradiction work: construct a
subproblem, prove that its retained structure maps faithfully to a tournament
or OCF conflict graph, and then use a forbidden H-value or forced-expansion
theorem to rule out the subproblem.

This scout generalizes that pattern into reusable proof/disproof techniques
and applies them to live repository surfaces.  It is an idea generator and
route audit, not a theorem prover.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product


def independence_poly_at_two(n: int, edges: set[tuple[int, int]]) -> tuple[list[int], int]:
    """Return independence coefficients and I(G,2) for a small graph."""
    coeffs = [0] * (n + 1)
    edge_set = {tuple(sorted(e)) for e in edges}
    vertices = range(n)
    for mask in range(1 << n):
        chosen = [i for i in vertices if mask & (1 << i)]
        ok = True
        for a, b in combinations(chosen, 2):
            if tuple(sorted((a, b))) in edge_set:
                ok = False
                break
        if ok:
            coeffs[len(chosen)] += 1
    value = sum(c * (2**i) for i, c in enumerate(coeffs))
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs, value


def complete_graph(n: int) -> set[tuple[int, int]]:
    return {tuple(sorted(e)) for e in combinations(range(n), 2)}


def path_graph(n: int) -> set[tuple[int, int]]:
    return {tuple(sorted((i, i + 1))) for i in range(n - 1)}


def remove_edges(edges: set[tuple[int, int]], missing: list[tuple[int, int]]) -> set[tuple[int, int]]:
    out = set(edges)
    for edge in missing:
        out.discard(tuple(sorted(edge)))
    return out


@dataclass(frozen=True)
class Technique:
    name: str
    move: str
    retains: frozenset[str]
    destroys: frozenset[str]
    contradiction_power: int
    transfer_burden: int
    abstraction: int
    route_kind: str


@dataclass(frozen=True)
class Surface:
    name: str
    live_problem: str
    required: frozenset[str]
    dangerous_forgetting: frozenset[str]
    target: str
    incoming_signal: str


def score(tech: Technique, surface: Surface) -> tuple[int, int, int, int, int]:
    retained = len(tech.retains & surface.required)
    destroys_needed = len(tech.destroys & surface.required)
    avoids_danger = len(tech.retains & surface.dangerous_forgetting)
    breadth = len(tech.retains)
    raw = (
        4 * retained
        + 2 * avoids_danger
        + tech.contradiction_power
        + tech.abstraction
        - 3 * destroys_needed
        - tech.transfer_burden
    )
    return raw, retained, avoids_danger, -destroys_needed, breadth


def verdict(score_tuple: tuple[int, int, int, int, int], tech: Technique) -> str:
    raw, retained, avoids_danger, neg_destroy, _ = score_tuple
    if neg_destroy < 0:
        return "reject_until_sidecar_restored"
    if raw >= 22 and tech.contradiction_power >= 4:
        return "candidate_disproof_or_terminal_obstruction"
    if raw >= 18:
        return "candidate_proof_route_or_required_ledger"
    if retained >= 2:
        return "diagnostic_sidecar"
    return "analogy_only"


def skeletons() -> list[tuple[str, int, list[int], int, str]]:
    cases: list[tuple[str, int, set[tuple[int, int]], str]] = []
    cases.append(("K3_for_H7_component", 3, complete_graph(3), "I(K3,2)=7; THM-029 says a K3 Omega component expands."))
    cases.append(("P4_for_H21_component", 4, path_graph(4), "I(P4,2)=21; THM-079 rules out P4 Omega components."))
    cases.append((
        "K6_minus_two_disjoint_edges",
        6,
        remove_edges(complete_graph(6), [(0, 1), (2, 3)]),
        "One H=21 skeleton: alpha1=6, alpha2=2; tournament forcing adds 5-cycles.",
    ))
    cases.append((
        "K8_minus_one_edge",
        8,
        remove_edges(complete_graph(8), [(0, 1)]),
        "One H=21 skeleton: alpha1=8, alpha2=1; cascade forcing blocks it at n=8.",
    ))
    cases.append(("K10", 10, complete_graph(10), "One H=21 skeleton: alpha1=10, alpha2=0; pancyclic lower bounds close general n."))
    out = []
    for name, n, edges, note in cases:
        coeffs, value = independence_poly_at_two(n, edges)
        out.append((name, n, coeffs, value, note))
    return out


TECHNIQUES = [
    Technique(
        "permanent_H_gap_transfer",
        "prove constructed subproblem has faithful H in {7,21}; contradiction by THM-029/THM-115",
        frozenset({"faithful_tournament_functor", "OCF_component_value", "permanent_gap", "minimal_skeleton"}),
        frozenset({"continuous_scale", "endpoint_owner"}),
        7,
        6,
        6,
        "disproof",
    ),
    Technique(
        "component_value_factorization",
        "factor a proposed witness into Omega components and block forbidden factors like 7",
        frozenset({"OCF_component_value", "multiplicative_factor", "component_connectivity", "minimal_skeleton"}),
        frozenset({"component_internal_geometry"}),
        6,
        4,
        5,
        "disproof",
    ),
    Technique(
        "forced_expansion_closure",
        "assume a small tournament/OCF skeleton, then force additional cycles/sidecars until value jumps",
        frozenset({"minimal_skeleton", "forced_extra_vertex", "local_expansion", "closure_operator"}),
        frozenset({"global_distribution"}),
        6,
        5,
        7,
        "disproof",
    ),
    Technique(
        "deletion_invariant_induction",
        "find an inert vertex/coordinate whose deletion preserves the conflict object, then induct",
        frozenset({"deletion_deck", "inert_coordinate", "monotone_subobject", "base_case"}),
        frozenset({"new_cycle_through_deleted_coordinate"}),
        5,
        4,
        5,
        "proof_or_disproof",
    ),
    Technique(
        "ocf_evaluation_lift",
        "replace H by I(Omega,x) or typed OCF data before evaluating at x=2",
        frozenset({"typed_independence_polynomial", "cycle_length_type", "evaluation_point", "OCF_component_value"}),
        frozenset({"route_label", "endpoint_owner"}),
        4,
        5,
        8,
        "proof_or_disproof",
    ),
    Technique(
        "ranked_proof_carrier_tournament",
        "make proof obligations/certificates vertices and orient by retained predicate payload",
        frozenset({"proof_obligation", "retained_payload", "destroyed_coordinate", "terminal_exit"}),
        frozenset({"numeric_H_claim"}),
        3,
        2,
        7,
        "proof_navigation",
    ),
    Technique(
        "cycle_class_observability_matrix",
        "encode residual cochains versus certificate cycles and row-reduce before naming debt",
        frozenset({"cycle_class_image", "certificate_cycle", "dual_annihilator", "named_debt"}),
        frozenset({"scalar_route_label"}),
        4,
        4,
        8,
        "proof_navigation",
    ),
    Technique(
        "median_center_legality_tournament",
        "route triples must have a legal sidecar median or expose first missing sidecar",
        frozenset({"route_triple", "median_center", "sidecar_closure", "first_missing_sidecar"}),
        frozenset({"raw_route_clique"}),
        4,
        4,
        7,
        "proof_navigation",
    ),
    Technique(
        "edge_flip_stress_disprover",
        "flip local edges/quotient coordinates and require the claimed invariant to survive",
        frozenset({"edge_flip", "quotient_fiber", "mixed_status", "lost_coordinate"}),
        frozenset({"single_scalar_rank"}),
        5,
        3,
        6,
        "disproof",
    ),
    Technique(
        "hypertournament_blocker_allocation",
        "make simultaneous blockers hyperedges; disprove a proof by showing no allocation covers all debts",
        frozenset({"blocker_hyperedge", "allocation_capacity", "simultaneous_failure", "terminal_debt"}),
        frozenset({"pairwise_only"}),
        5,
        5,
        8,
        "proof_or_disproof",
    ),
    Technique(
        "quotient_round_trip_obstruction",
        "push a quotient through source/deletion/extension and require round-trip payload consistency",
        frozenset({"quotient", "round_trip", "observer_cut", "payload_orbit", "extension_defect"}),
        frozenset({"unmarked_class_count"}),
        5,
        4,
        8,
        "disproof",
    ),
    Technique(
        "tropical_capacity_gap_tournament",
        "orient sidecars by slack/capacity; contradiction is a negative cycle in required capacities",
        frozenset({"capacity_slack", "negative_cycle", "dual_certificate", "allocation_capacity"}),
        frozenset({"exact_discrete_owner"}),
        4,
        5,
        9,
        "proof_or_disproof",
    ),
    Technique(
        "valuation_residue_lift_audit",
        "treat p-adic/residue analogies as lift charts and reject them unless ramification and residue fibers survive",
        frozenset({"valuation_layer", "residue_fiber", "ramification_status", "lift_obstruction"}),
        frozenset({"real_metric_geometry", "Steiner_angle"}),
        3,
        4,
        8,
        "guardrail",
    ),
    Technique(
        "semantic_tournament_of_statements",
        "vertices are statements; edge means one statement's proof payload subsumes/refutes another",
        frozenset({"statement_vertex", "subsumption_edge", "contradiction_edge", "payload_basis"}),
        frozenset({"truthiness_score"}),
        3,
        3,
        9,
        "proof_navigation",
    ),
    Technique(
        "redei_parity_certificate",
        "encode a count as tournament Hamiltonian paths and use Redei parity to forbid even targets",
        frozenset({"hamiltonian_path_count", "parity_constraint", "certificate_engine", "faithful_tournament_functor"}),
        frozenset({"orientation_payload"}),
        5,
        3,
        6,
        "disproof",
    ),
    Technique(
        "landau_score_feasibility",
        "map data to a tournament score sequence and reject violations of Landau inequalities",
        frozenset({"score_sequence", "landau_inequality", "feasibility_spectrum", "faithful_tournament_functor"}),
        frozenset({"cycle_payload"}),
        5,
        3,
        6,
        "disproof",
    ),
    Technique(
        "cycle_census_hole_transfer",
        "treat missing cycle-count fibers as forbidden spectra distinct from forbidden H gaps",
        frozenset({"cycle_count_spectrum", "moment_hole", "c3_c5_fiber", "spectral_layer"}),
        frozenset({"OCF_component_value"}),
        6,
        5,
        8,
        "disproof",
    ),
    Technique(
        "score_exchange_nontransitivity",
        "use an improvement tournament to test whether greedy exchange proves optimality or exposes local minima",
        frozenset({"improvement_tournament", "local_minima", "finite_check_bound", "exchange_failure"}),
        frozenset({"global_optimality_scalar"}),
        4,
        3,
        9,
        "proof_or_disproof",
    ),
    Technique(
        "winding_tie_apex_audit",
        "separate apex tied diameters from forbidden K3/H=7 conflict claims",
        frozenset({"antipodal_matching", "tie_event", "nondegenerate_scale", "false_H_bridge"}),
        frozenset({"K3_numeric_analogy", "numeric_H_claim"}),
        5,
        2,
        7,
        "guardrail",
    ),
    Technique(
        "certificate_engine_spectrum_generator",
        "enumerate tournament invariants and promote persistent missing values to candidate certificates",
        frozenset({"achieved_spectrum", "persistent_gap", "self_test_certificate", "invariant_family"}),
        frozenset({"handpicked_single_gap"}),
        6,
        4,
        9,
        "proof_or_disproof",
    ),
    Technique(
        "single_component_H_ladder_certificate",
        "turn H gaps into clique-Omega realizability gaps, with H=7 as K3 and H=21 as K10",
        frozenset({"single_component_H_gap", "clique_omega_realizability", "omega_sparsity", "persistent_gap"}),
        frozenset({"disconnected_factor_shadow"}),
        7,
        3,
        8,
        "disproof",
    ),
]


SURFACES = [
    Surface(
        "LRC14_observer_gluing",
        "HYP-3098 needs chart overlap, not raw H or raw pair mass.",
        frozenset({"proof_obligation", "retained_payload", "destroyed_coordinate", "terminal_exit", "sidecar_closure"}),
        frozenset({"continuous_scale", "endpoint_owner", "single_scalar_rank"}),
        "prove gluing theorem or name first failed overlap",
        "THM-573/THM-577 plus divisor-loaded arc stress",
    ),
    Surface(
        "THM577_cap_dip",
        "Pairwise caps are symbolic through j=3; j=4,5 are finite higher-order remainders.",
        frozenset({"minimal_skeleton", "forced_extra_vertex", "closure_operator", "typed_independence_polynomial", "cycle_length_type"}),
        frozenset({"pairwise_only", "single_scalar_rank"}),
        "turn inclusion-exclusion order vector into moment/Perron certificate",
        "closed-form forbidden-arc overlap and cap inclusion-exclusion anatomy",
    ),
    Surface(
        "HYP3094_K33_state_lift",
        "Positive-open K33 rows must be routed to state-lift debt, not confused with covering discharge.",
        frozenset({"proof_obligation", "route_triple", "sidecar_closure", "first_missing_sidecar", "terminal_debt"}),
        frozenset({"endpoint_owner", "raw_route_clique", "single_scalar_rank"}),
        "construct THM-572 state-lift or identify first missing sidecar",
        "nested-refinement versus cross-handoff split",
    ),
    Surface(
        "q_cusp_principal_part",
        "q-series tails are legal only after finite polar debt is named.",
        frozenset({"deletion_deck", "base_case", "named_debt", "dual_annihilator", "terminal_exit"}),
        frozenset({"truthiness_score", "single_scalar_rank"}),
        "finite principal-part theorem or explicit polar debt",
        "HYP-3078 q-Pochhammer gate",
    ),
    Surface(
        "sixth_power_support_six",
        "3-vs-3 sixth-power collisions are native; 2-vs-2 are padded shadows.",
        frozenset({"minimal_skeleton", "forced_extra_vertex", "route_triple", "payload_basis", "named_debt"}),
        frozenset({"single_scalar_rank", "pairwise_only"}),
        "separate native support-six from degeneracy guards",
        "HYP-3076/HYP-3080 sixth-power sidecars",
    ),
    Surface(
        "route_state_median_hull",
        "Raw route labels have empty centers; legal sidecars create unique median centers.",
        frozenset({"route_triple", "median_center", "sidecar_closure", "first_missing_sidecar"}),
        frozenset({"raw_route_clique", "truthiness_score"}),
        "turn median center into non-tautological packet proof",
        "HYP-3070/HYP-3074 route-state closure",
    ),
    Surface(
        "p_adic_Roth_Steiner_forum",
        "Incoming speculative post links 2-adic residue fields, Roth discrepancy, and Steiner ratios.",
        frozenset({"valuation_layer", "residue_fiber", "ramification_status", "capacity_slack", "payload_basis"}),
        frozenset({"real_metric_geometry", "Steiner_angle", "truthiness_score"}),
        "reject as scalar analogy unless lift, metric, and discrepancy payloads are retained",
        "base64 forum post after S64",
    ),
    Surface(
        "H2963_packet_bank",
        "The packet bank needs a reusable disproof/proof route audit for new scalar proposals.",
        frozenset({"quotient_fiber", "mixed_status", "lost_coordinate", "retained_payload", "terminal_exit"}),
        frozenset({"single_scalar_rank", "unmarked_class_count"}),
        "automatic stale-quotient rejection before new scalar claims",
        "HYP-2963/LTI-100 fiber-mixing stack",
    ),
    Surface(
        "certificate_engine_kps",
        "S31ah turns forbidden H into a reusable tournament certificate engine.",
        frozenset({"achieved_spectrum", "persistent_gap", "self_test_certificate", "invariant_family", "parity_constraint", "landau_inequality"}),
        frozenset({"handpicked_single_gap", "numeric_H_claim"}),
        "discover new forbidden spectra instead of handpicking H=7/H=21",
        "tournament_certificate_engine_kps self-tests Paley T7 and H/parity certificates",
    ),
    Surface(
        "single_component_H_ladder_kps",
        "S31ah validates the generator: single-component H gaps are exactly 7 and 21, i.e. K3 and K10 Omega gaps.",
        frozenset({"single_component_H_gap", "clique_omega_realizability", "omega_sparsity", "persistent_gap"}),
        frozenset({"disconnected_factor_shadow", "handpicked_single_gap", "numeric_H_claim"}),
        "promote K_m/Omega realizability into a reusable obstruction-transfer certificate",
        "tournament_single_component_H_ladder_kps and tournament_spectrum_discovery_kps",
    ),
    Surface(
        "cap_optimality_exchange_S65",
        "S65 shows cap minimizers are bounded but the single-swap exchange tournament has spurious local minima.",
        frozenset({"improvement_tournament", "local_minima", "finite_check_bound", "exchange_failure"}),
        frozenset({"global_optimality_scalar", "pairwise_only"}),
        "turn nontransitivity into a bounded finite-check scheduler",
        "lrc_cap_optimality_exchange_macmini_S65",
    ),
    Surface(
        "baby_Hodge_cycle_hole",
        "The n=6 baby-Hodge hole is a c5 spectral hole, not the H=7/H=21 alpha obstruction.",
        frozenset({"cycle_count_spectrum", "moment_hole", "c3_c5_fiber", "spectral_layer"}),
        frozenset({"OCF_component_value", "numeric_H_claim"}),
        "separate cycle-count spectrum holes from OCF H-gap transfers",
        "baby_hodge_forbidden_h_crosscheck_opus",
    ),
    Surface(
        "apex7_forbiddenH_bridge_audit",
        "Apex 7 is an antipodal matching count, not an Omega=K3/H=7 event.",
        frozenset({"antipodal_matching", "tie_event", "nondegenerate_scale", "false_H_bridge"}),
        frozenset({"K3_numeric_analogy", "numeric_H_claim"}),
        "reject literal apex-7 to forbidden-H=7 bridges without a new transfer",
        "apex7_vs_forbiddenH7_bridge_audit",
    ),
]


def generated_cross_moves() -> list[str]:
    verbs = ["gap-transfer", "forced-expansion", "deletion-induction", "median-center", "capacity-cycle"]
    carriers = ["OCF-component", "observer-chart", "residue-lift", "q-principal-part", "support-six", "packet-bank"]
    guards = ["faithful-functor", "sidecar-survival", "finite-base", "dual-annihilator", "named-debt"]
    moves = []
    for verb, carrier, guard in product(verbs, carriers, guards):
        if verb == "gap-transfer" and guard not in {"faithful-functor", "finite-base"}:
            continue
        if carrier == "residue-lift" and guard == "dual-annihilator":
            continue
        moves.append(f"{verb}::{carrier}::{guard}")
    return moves


def tournament_fingerprint(applications: list[tuple[str, str, tuple[int, int, int, int, int], str]]) -> None:
    ranked = sorted(applications, key=lambda x: (x[2], x[0], x[1]), reverse=True)
    top = ranked[:10]
    print("[D] TOURNAMENT_ANALYSIS")
    print("  vertices_are: technique applications / proof-obligation carriers, not runners")
    print("  pairwise_observable: retained payload, contradiction power, transfer burden, abstraction, destroyed-coordinate safety")
    print("  tie_hamiltonian_path:")
    for i, (tech, surface, score_tuple, why) in enumerate(top, 1):
        raw = score_tuple[0]
        print(f"    {i}. {tech} @ {surface}  score={raw}  verdict={why}")
    scores = sorted({x[2][0] for x in applications})
    hist = {s: sum(1 for x in applications if x[2][0] == s) for s in scores}
    print(f"  score_histogram: {hist}")
    print("  directed_cycles_known: not computed as dominance uses lexicographic score; future run should compare independent gauges")
    print("  challenged_assumption: H-value analogy is proof currency without a faithful transfer functor")


def main() -> None:
    print("S259 TOURNAMENT OBSTRUCTION-TRANSFER ATLAS")
    print("=" * 78)
    print("[A] H-GAP SKELETONS")
    for name, n, coeffs, value, note in skeletons():
        print(f"  {name}: n={n} alpha={coeffs} I(G,2)={value} :: {note}")

    print()
    print("[B] GENERATED TECHNIQUE FAMILIES")
    moves = generated_cross_moves()
    print(f"  generated_cross_moves={len(moves)}")
    for move in moves[:24]:
        print(f"  {move}")
    if len(moves) > 24:
        print(f"  ... {len(moves) - 24} more")

    print()
    print("[C] APPLICATION MATRIX")
    applications: list[tuple[str, str, tuple[int, int, int, int, int], str]] = []
    for surface in SURFACES:
        rows = []
        for tech in TECHNIQUES:
            score_tuple = score(tech, surface)
            v = verdict(score_tuple, tech)
            applications.append((tech.name, surface.name, score_tuple, v))
            rows.append((score_tuple, tech, v))
        rows.sort(key=lambda x: (x[0], x[1].name), reverse=True)
        print(f"surface={surface.name}")
        print(f"  live_problem: {surface.live_problem}")
        print(f"  incoming_signal: {surface.incoming_signal}")
        print(f"  target: {surface.target}")
        for score_tuple, tech, v in rows[:4]:
            raw, retained, avoids, neg_destroy, breadth = score_tuple
            print(
                "  "
                f"{tech.name}: score={raw} retained={retained} danger_retained={avoids} "
                f"destroyed_needed={-neg_destroy} verdict={v}"
            )
        rejected = [tech.name for score_tuple, tech, v in rows if v == "reject_until_sidecar_restored"]
        if rejected:
            print(f"  rejected_until_sidecar_restored: {', '.join(rejected[:5])}")

    print()
    tournament_fingerprint(applications)

    print()
    print("[E] SYNTHESIS")
    print("  1. H=7/H=21 should be reused as forbidden local spectra only after a faithful transfer functor is proved.")
    print("  2. The strongest general move is forced expansion: assume the attractive small skeleton, then prove closure adds payload that jumps the invariant.")
    print("  3. For LRC14, the tournament vertices should be chart-overlap obligations; a raw H analogy is a guardrail failure.")
    print("  4. For THM-577, the new tournament-like object is the inclusion-exclusion order vector, not the final cap value.")
    print("  5. For speculative p-adic/Roth/Steiner ideas, valuation-residue lift audits are useful only as rejection filters until metric and residue fibers are formalized.")
    print("  6. S31ah/S65 add a certificate-engine layer: parity, Landau, cycle-census holes, exchange nontransitivity, and apex-tie audits are distinct levers.")
    print("  7. The cap exchange tournament being nontransitive is progress: it turns a hoped-for greedy proof into a bounded finite-check obligation.")


if __name__ == "__main__":
    main()
