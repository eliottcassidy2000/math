#!/usr/bin/env python3
"""HYP-3141 scout: tournament edges as recursive tip/tail witnesses.

This is a synthesis scout, not a proof.  It treats a directed tournament edge
as an information-bearing witness packet:

  tail deletion payload + tip extension payload + observer-cut orbit
  + commutator defect + terminal proof exit.

Tournament Analysis uses edge-witness carrier types as vertices, not runners,
raw arcs, or scalar tournament values.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


FEATURE_WEIGHTS = {
    "lrc_predicate": 5,
    "tail_recursion": 4,
    "tip_recursion": 4,
    "cut_payload_orbit": 5,
    "endpoint_role": 4,
    "cross_orientation": 4,
    "coordinate_resurrection": 4,
    "proof_circuit_input": 3,
    "finite_address_exit": 5,
    "observer_gluing_exit": 5,
    "exact_gap": 4,
    "decorrelation_floor": 4,
    "minkowski_relation": 4,
    "ising_wall": 4,
    "phi4_quartic": 4,
    "de_moivre_fold": 3,
    "global_consistency_quotient": 2,
    "resolvent_middle_payload": 2,
    "recursion_commutator": 4,
    "information_gain": 3,
    "raw_arc_penalty": -8,
}


@dataclass(frozen=True)
class EdgeCarrier:
    name: str
    sources: str
    definition: str
    tail_payload: str
    tip_payload: str
    witness_condition: str
    destroyed_if_scalarized: str
    features: dict[str, int]

    def score(self) -> int:
        return sum(FEATURE_WEIGHTS[key] * value for key, value in self.features.items())


CARRIERS = [
    EdgeCarrier(
        "recursive_tip_tail_edge_witness",
        "HYP-3049, HYP-3054, HYP-3056, HYP-3118",
        "A directed edge e=(u->v) is a witness when tail deletion and tip extension commute up to a named repair.",
        "delete/condition on the tail u; record what endpoint-owner, finite-address, or route data becomes lost.",
        "extend/observe at the tip v; record incident word, sector deck, ear payload, or next proof demand.",
        "commutator(delete_tail, extend_tip) is zero, reconstructed, dual-annihilated, descended, boundary-stopped, or named debt.",
        "old/new endpoint role, cross-sector orientation, and the changed proof predicate.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "endpoint_role": 1,
            "cross_orientation": 1,
            "coordinate_resurrection": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
            "finite_address_exit": 1,
            "observer_gluing_exit": 1,
        },
    ),
    EdgeCarrier(
        "observer_cut_payload_orbit",
        "HYP-3054, HYP-3056, HYP-3120",
        "The edge is the next observer/cut orbit inside a coarse quotient fiber.",
        "boundary slice and visible-fiber automorphism orbit before the cut.",
        "extended object shadow after the observer is inserted.",
        "the orbit separates every route/status-changing pair or has an explicit discharge mode.",
        "the first forgotten coordinate that changes boundary/open status or route.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "coordinate_resurrection": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
            "observer_gluing_exit": 1,
        },
    ),
    EdgeCarrier(
        "directed_edge_perspective",
        "HYP-3049, HYP-3106, HYP-3133, HYP-3134",
        "A000568 first-failure edge: old root to new observer, then quotient to directed-edge perspective.",
        "old-root role and deletion-parent cache.",
        "new-observer incident word, ordered-pair sector deck, A000568 extension shadow, envelope position, global-consistency class, and cross-sector orientation.",
        "sector cross/full deck separates the first converse pair and HYP-3134 names when the A000568 middle quotient may forget paired child payload.",
        "old/new endpoint role, A000568 extension shadow, edge-envelope position, child-gluing status, and cross-sector orientation word.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "endpoint_role": 1,
            "cross_orientation": 1,
            "global_consistency_quotient": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "source_sink_apex_arc",
        "HYP-2008",
        "The source-sink arc is the apex tile, block-extreme switch, and LRC observer-gap witness.",
        "source-side transitive/semicircle branch and largest-gap deletion.",
        "sink-side wrap/back-arc branch and AP-hard boundary cycle closure.",
        "apex orientation selects the correct loneliness regime and names the hard wrap branch.",
        "whether the arc is an observer gap, a back-arc, or only a raw tournament edge.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "endpoint_role": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "finite_address_phi_edge",
        "HYP-3083, HYP-3116, HYP-3117, HYP-3120",
        "An edge witnesses proof progress when it transports a finite-address Phi/P packet to a terminal receiver.",
        "finite address and endpoint activation before edge contraction.",
        "observer-gluing, endpoint Phi/P, or Lean frontier receiver after edge transport.",
        "transport preserves exact gap or emits missing-input vector and terminal debt.",
        "Phi/P activation vector, finite address, and observer-gluing receiver.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "proof_circuit_input": 1,
            "finite_address_exit": 1,
            "observer_gluing_exit": 1,
            "exact_gap": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "asano_tip_contraction_edge",
        "HYP-3140, HYP-3137, HYP-3136, HYP-3135, HYP-3134, HYP-3133, HYP-3132, HYP-3131, HYP-3130, HYP-3129, HYP-3128, HYP-3127, HYP-3126, HYP-3125, HYP-3121, HYP-3124",
        "A multi-far tip edge transports one single-far Lee-Yang factor through Asano contraction, wide decoupling, or SPEC-certified obstruction splitting.",
        "bounded-core floor, contraction context, and already-contracted tip factors before deleting or averaging the tail.",
        "single-far factor, zero-free region certificate, SPEC/equidistribution rate, minorant tail, far-zero-push status, obstruction status, A000568 global-consistency class, resolvent middle payload, GF payload layer, fiber-PGF conditional first moment, and next contraction order at the tip.",
        "contract-tip-first and decouple-wide-tip-first agree within the HYP-3126/HYP-3129/HYP-3130 budget, exploit HYP-3131 zero push, reduce to the HYP-3132 bounded-core resolvent, use the HYP-3133/HYP-3134 extension-envelope shadow, preserve the HYP-3135 middle resolvent packet, route through the HYP-3136 integrated floor factorization, HYP-3137 GF payload atlas, and HYP-3140 fiber-PGF certificate, or discharge to HYP-3128 obstruction debt.",
        "the Asano contraction order, single-far factor identity, zero-free certificate, SPEC/minorant certificate, far-zero-push status, bounded-core resolvent status, A000568 envelope/gluing status, resolvent pair/triple payload, GF coefficient/root-locus/log-derivative payload, fiber-PGF conditional moment, obstruction status, and wide-decoupling error budget.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "finite_address_exit": 1,
            "observer_gluing_exit": 1,
            "exact_gap": 1,
            "decorrelation_floor": 1,
            "phi4_quartic": 1,
            "global_consistency_quotient": 1,
            "resolvent_middle_payload": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "coordinate_resurrection_edge",
        "HYP-3118",
        "The edge is a quotient/section pair that resurrects exactly the coordinate demanded by the next operation.",
        "destroyed coordinate vector at the tail quotient.",
        "minimal legal repair cover or adjoint section at the tip.",
        "repair cover rank strictly decreases missing-coordinate debt or reaches a proof minterm.",
        "base stalk, live section, repair cover rank, and adjoint status.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "coordinate_resurrection": 1,
            "proof_circuit_input": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
            "observer_gluing_exit": 1,
        },
    ),
    EdgeCarrier(
        "lee_yang_ear_motion_edge",
        "HYP-3137, HYP-3109, HYP-3112, HYP-3113",
        "A one-runner ear edge witnesses root motion only with the hidden A_t payload retained.",
        "base PGF root curve and root-stratum state before adding the runner.",
        "tip ear payload A_t(E,a), full PGF coefficient/root-locus payload, root-motion reconstruction, and root-wall contact status.",
        "G_{E+a}(z)=G_E(z)+(z^-1-1)A(z) reconstructs the motion or names root-collision debt.",
        "ear payload vector, coefficient layer, root locus, parity/nesting status, and root trajectory.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "coordinate_resurrection": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "ising_domain_wall_edge",
        "HYP-3115, HYP-3109, HYP-3112",
        "A one-swap edge across the #real-root spin boundary is a domain-wall witness, not just a graph edge.",
        "tail root-stratum spin, base PGF roots, and relation/ear payload before the swap.",
        "tip root-stratum spin, partition-zero packet, root-collision side, and legal ear or observer-gluing exit after the swap.",
        "the wall is classified as legal partition-zero/root-collision ear, observer-gluing exit, finite-address exit, or forbidden debt.",
        "which of the 10084 HYP-3115 wall edges is a legal proof transition versus a forbidden phase boundary, plus the finite-volume partition-zero locus that made the spin analogy legal.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "coordinate_resurrection": 1,
            "ising_wall": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "minkowski_relation_wall_edge",
        "HYP-3062, HYP-3111, HYP-3115",
        "A relation-rich edge is useful only when it retains the lattice body, low-height wall, and successive-minima proxy.",
        "relation-lattice body, short-wall owner, and covolume/successive-minima data before tail deletion.",
        "new short-relation class, low-height wall deletion status, covolume threshold, and convex-body forcing sidecar at the tip.",
        "Minkowski pressure becomes proof currency only after the wall is deleted, reconstructed, or named as residual debt.",
        "the convex body, relation owner, low-height wall class, covolume threshold, and why raw volume pressure was legal.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "minkowski_relation": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "de_moivre_stationary_quintic_edge",
        "HYP-3139, HYP-3138, HYP-3135, HYP-3132, HYP-3110, HYP-3111, HYP-3115",
        "A root-collision edge carries the translated G'(z) quintic residual and branch orbit as an algebraic stress packet.",
        "stationary quintic coefficients, auxiliary quadratic choice, branch choice, and de Moivre residual before crossing the wall.",
        "new branch orbit, auxiliary-quadratic status, biquadratic resolvent status, reflection-core block status, k=8 reflection-fold adjoint status, odd-coordinate resurrection table, pair/triple elementary-symmetric payload, forbidden y^2 term, linear-gap defect, and finite-address/observer-gluing handoff.",
        "the fold supplies a solvable local model, the HYP-3132 bounded-core biquadratic discharge, the HYP-3135 middle-layer packet, the HYP-3138 reflection-fold adjoint, the HYP-3139 inner-shell block bound, or a named branch-coordinate debt; it is not an extremality scalar.",
        "fifth-root branch data, auxiliary-quadratic branch, biquadratic resolvent status, reflection-core block data, k=8 fold adjoint, odd-coordinate resurrection data, pair/triple payload, and stationary-collision geometry behind a single residual number.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "coordinate_resurrection": 1,
            "de_moivre_fold": 1,
            "resolvent_middle_payload": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "phi4_quartic_stabilizer_edge",
        "HYP-3139, HYP-3138, HYP-3135, HYP-3132, HYP-3122 phi4, HYP-3103, HYP-3113",
        "A binding-row edge carries the cap dip as a quartic-cumulant stabilizer, not as a scalar correction.",
        "quadratic pair-Pascal cap shadow and pre-dip miss-count cumulants at the tail.",
        "kappa4 sign, effective lambda, S4 dip payload, biquadratic resolvent status, reflection-core block, center/boundary leakage, reflection-fold adjoint status, middle elementary-symmetric layer, and odd/even ear-cumulant handoff at the tip.",
        "a uniform quartic-cumulant bound, HYP-3132 biquadratic resolvent, HYP-3135 packet theorem, HYP-3138 fold-adjoint certificate, or HYP-3139 block-page leakage ceiling discharges the dip, or the edge names phi4/Lee-Yang residual debt.",
        "the non-pairwise cap dip, kappa4 sign, reflection-block leakage, fold-adjoint legality, resolvent middle payload, and whether lambda>0 is a legal stabilizer.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "coordinate_resurrection": 1,
            "phi4_quartic": 1,
            "resolvent_middle_payload": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "endpoint_owner_transfer_edge",
        "HYP-3045, HYP-3034, HYP-3037",
        "The edge carries endpoint-owner current through a route/status-changing cut.",
        "coarse endpoint word and owner-deletion persistence at the tail.",
        "external endpoint-owner strip, owner-transfer delta, or residual capacitor cut at the tip.",
        "owner transfer splits a mixed route/status fiber or names K33/F7/THM-572 debt.",
        "external owner names collapsed by coarse endpoint counts.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "endpoint_role": 1,
            "cut_payload_orbit": 1,
            "observer_gluing_exit": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "rectangle_hourglass_flow_edge",
        "HYP-3053",
        "The edge is a diagonal-layer flow line whose rectangle/hourglass residues certify compression.",
        "layer-potential tail data and spanning-tree line basis.",
        "tip bridge line plus rectangle/hourglass cycle defect.",
        "all defects vanish, are reconstructed from potentials, or become named hidden-coordinate witnesses.",
        "cycle defects, owner/barcode support, and nonlinear LRC payload on the line.",
        {
            "tail_recursion": 1,
            "tip_recursion": 1,
            "cut_payload_orbit": 1,
            "cross_orientation": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "proof_circuit_missing_input_edge",
        "HYP-3116, HYP-3117",
        "An edge is a proof-circuit transition that must reduce the missing-input vector.",
        "current essential-input set, missing inputs, and size/depth/fanin/uniformity budget before applying the edge.",
        "certificate minterm, repaired sidecar, next missing input, or uniformity-guard discharge after the edge.",
        "missing vector strictly decreases under a named size/depth/fanin/uniformity guard or the edge is rejected as scalar telemetry.",
        "which input was actually repaired by the edge, and whether the circuit family remains uniform enough to be proof currency.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "proof_circuit_input": 1,
            "coordinate_resurrection": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "node3_decorrelation_edge",
        "HYP-3121, THM-571, HYP-2968",
        "The edge witnesses quasi-independence between R-safe and Q-lonely lift events.",
        "small-part safe event at the tail, with resonance frequencies retained.",
        "apex/lift lonely event at the tip, with decorrelation floor R' retained.",
        "intersection is positive by a uniform R' floor, not by Bonferroni.",
        "low-frequency resonance correction and the fact that Bonferroni failed.",
        {
            "lrc_predicate": 1,
            "tail_recursion": 1,
            "tip_recursion": 1,
            "decorrelation_floor": 1,
            "recursion_commutator": 1,
            "information_gain": 1,
        },
    ),
    EdgeCarrier(
        "raw_tournament_arc",
        "negative control",
        "A bare orientation u->v without tail/tip payload, observer orbit, or proof predicate.",
        "none.",
        "none.",
        "not a proof witness until sidecars are attached.",
        "everything except the binary orientation.",
        {
            "raw_arc_penalty": 1,
        },
    ),
]


TIE_PATH = [
    "recursive_tip_tail_edge_witness",
    "asano_tip_contraction_edge",
    "finite_address_phi_edge",
    "coordinate_resurrection_edge",
    "observer_cut_payload_orbit",
    "ising_domain_wall_edge",
    "directed_edge_perspective",
    "endpoint_owner_transfer_edge",
    "lee_yang_ear_motion_edge",
    "phi4_quartic_stabilizer_edge",
    "proof_circuit_missing_input_edge",
    "de_moivre_stationary_quintic_edge",
    "minkowski_relation_wall_edge",
    "node3_decorrelation_edge",
    "source_sink_apex_arc",
    "rectangle_hourglass_flow_edge",
    "raw_tournament_arc",
]


NEW_SIGNALS = [
    "edge_witness_packet_id",
    "edge_tail_payload_word",
    "edge_tip_payload_word",
    "tail_delete_recursion_depth",
    "tip_extend_recursion_depth",
    "tip_tail_commutator_defect",
    "edge_cut_payload_orbit_id",
    "old_new_endpoint_role",
    "cross_sector_orientation_word",
    "edge_gf_carrier_type",
    "edge_coefficient_payload_layer",
    "edge_pgf_root_locus_status",
    "edge_log_derivative_cumulant_status",
    "edge_fiber_pgf_word",
    "edge_q_masked_fiber_pgf_word",
    "edge_conditional_first_moment_floor_status",
    "edge_spec_resonance_lattice_status",
    "edge_a000568_extension_shadow",
    "edge_a000568_envelope_position",
    "edge_global_consistency_class",
    "edge_child_gluing_status",
    "edge_information_gain_rank",
    "edge_predicate_delta",
    "edge_coordinate_resurrection_cover",
    "edge_missing_input_delta",
    "edge_proof_circuit_size_depth_fanin",
    "edge_circuit_uniformity_guard",
    "edge_circuit_uniformity_guard_status",
    "edge_phi_p_activation_delta",
    "edge_minkowski_relation_wall_class",
    "edge_minkowski_covolume_threshold_status",
    "edge_successive_minima_proxy",
    "edge_ising_domain_wall_id",
    "edge_ising_partition_zero_locus_status",
    "edge_domain_wall_legal_exit",
    "edge_ear_payload_vector",
    "edge_phi4_quartic_cumulant_delta",
    "edge_phi4_lambda_sign",
    "edge_de_moivre_quintic_residual_delta",
    "edge_de_moivre_auxiliary_quadratic_status",
    "edge_de_moivre_biquadratic_resolvent_status",
    "edge_k8_reflection_fold_adjoint_status",
    "edge_odd_coordinate_resurrection_status",
    "edge_reflection_core_block_status",
    "edge_inner_shell_bound_status",
    "edge_center_boundary_leakage_status",
    "edge_resolvent_middle_payload_status",
    "edge_resolvent_pair_triple_layer",
    "edge_signed_spec_resolvent_packet_status",
    "edge_de_moivre_branch_orbit_word",
    "edge_rectangle_hourglass_residue",
    "edge_decorrelation_floor_status",
    "edge_asano_contraction_order_word",
    "edge_single_far_factor_id",
    "edge_zero_free_region_status",
    "edge_wide_decoupling_rate_bound",
    "edge_spec_certificate_status",
    "edge_lee_yang_obstruction_status",
    "edge_minorant_tail_certificate_status",
    "edge_far_zero_push_status",
    "edge_bounded_core_floor_exit",
    "edge_terminal_exit_or_debt",
]


def orient(a: EdgeCarrier, b: EdgeCarrier) -> tuple[EdgeCarrier, EdgeCarrier]:
    if a.score() > b.score():
        return a, b
    if b.score() > a.score():
        return b, a
    rank = {name: i for i, name in enumerate(TIE_PATH)}
    return (a, b) if rank[a.name] < rank[b.name] else (b, a)


def build_edges() -> dict[str, set[str]]:
    edges = {carrier.name: set() for carrier in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        winner, loser = orient(a, b)
        edges[winner.name].add(loser.name)
    return edges


def directed_3cycles(edges: dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(edges, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            count += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            count += 1
    return count


def sccs(edges: dict[str, set[str]]) -> list[list[str]]:
    names = list(edges)
    reverse = {name: set() for name in names}
    for src, dsts in edges.items():
        for dst in dsts:
            reverse[dst].add(src)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in edges[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, comp)

    for name in reversed(order):
        if name not in seen:
            comp: list[str] = []
            rdfs(name, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(edges: dict[str, set[str]]) -> tuple[int, list[str]]:
    names = list(edges)
    count = 0
    first_path: list[str] = []

    def backtrack(path: list[str], used: set[str]) -> None:
        nonlocal count, first_path
        if len(path) == len(names):
            count += 1
            if not first_path:
                first_path = path[:]
            return
        for candidate in names:
            if candidate in used:
                continue
            if path and candidate not in edges[path[-1]]:
                continue
            used.add(candidate)
            path.append(candidate)
            backtrack(path, used)
            path.pop()
            used.remove(candidate)

    for start in names:
        backtrack([start], {start})
    return count, first_path


def print_definition() -> None:
    print("MAP 1: recursive edge-witness definition")
    print("-" * 72)
    print("For a directed edge e=(tail -> tip) in a proof quotient q:")
    print("  Tail(e) = payload needed to delete/condition on the tail.")
    print("  Tip(e)  = payload needed to extend/observe at the tip.")
    print("  Orbit(e)= observer-cut payload orbit modulo visible automorphisms.")
    print("  Defect(e)= extend_tip(delete_tail(x)) - delete_tail(extend_tip(x)).")
    print("The edge is a witness only if Defect(e) is zero, reconstructed,")
    print("dual-annihilated, descended, boundary-stopped, or named residual debt.")
    print()
    print("Information view: a scalar quotient is legal only when the edge channel")
    print("retains the information needed by the next LRC predicate.  The raw arc")
    print("orientation is telemetry; the witness is the transported sidecar.")
    print()


def print_carriers() -> None:
    print("MAP 2: edge-witness carrier atlas")
    print("-" * 72)
    for carrier in sorted(CARRIERS, key=lambda c: (-c.score(), c.name)):
        print(f"{carrier.name} [{carrier.sources}] score={carrier.score()}")
        print(f"  definition: {carrier.definition}")
        print(f"  tail: {carrier.tail_payload}")
        print(f"  tip: {carrier.tip_payload}")
        print(f"  witness: {carrier.witness_condition}")
        print(f"  destroyed_if_scalarized: {carrier.destroyed_if_scalarized}")
    print()


def print_tournament() -> None:
    edges = build_edges()
    score_hist = Counter(carrier.score() for carrier in CARRIERS)
    path_count, first_path = hamiltonian_path_count(edges)
    print("TOURNAMENT ANALYSIS")
    print("-" * 72)
    print("vertices=edge-witness carrier types, not runners or raw arcs")
    print(
        "pairwise observable=which carrier preserves more LRC predicate "
        "information across both tail deletion and tip extension, while "
        "naming the observer-cut orbit and commutator defect"
    )
    print("switch/gauge=weighted proof-readiness score; ties follow the recursive receiver path")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycles(edges)}")
    print(f"scc_sizes={sorted(len(comp) for comp in sccs(edges))}")
    print(f"hamiltonian_path_count={path_count}")
    print("priority_path=" + " -> ".join(first_path))
    print()


def print_signals() -> None:
    print("NEW SIGNALS TO MEASURE")
    print("-" * 72)
    for signal in NEW_SIGNALS:
        print(signal)
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("-" * 72)
    considered = [
        "runners",
        "raw tournament vertices",
        "raw arcs",
        "directed edges",
        "ordered pairs",
        "A000568 extension shadows",
        "A000568 edge-envelope global-consistency quotients",
        "tip/tail endpoint roles",
        "observer-cut payload orbits",
        "sector-deck matrices",
        "source-sink apex arcs",
        "Lee-Yang ear edges",
        "Asano contraction tip edges",
        "phi4 quartic-stabilizer edges",
        "Minkowski relation-wall edges",
        "Ising domain-wall edges",
        "stationary de Moivre quintic branches",
        "bounded-core biquadratic resolvents",
        "k8 reflection-fold adjoints",
        "reflection-block proof pages",
        "resolvent middle-layer packets",
        "proof-circuit gates",
        "generating-function payload carriers",
        "fiber-PGF coefficient carriers",
        "finite-address proof edges",
        "coordinate-resurrection sections",
        "Node-3 decorrelation event pairs",
        "proof obligations",
    ]
    print("candidate vertices considered=" + ", ".join(considered))
    print("chosen vertices=edge-witness carrier types")
    print(
        "preserved predicate=progress toward LRC14 through exact Phi/P, "
        "finite address, observer gluing, coordinate resurrection, "
        "quasi-independence, or named residual debt"
    )
    print(
        "destroyed data=raw labels, old/new endpoint role, incident words, "
        "cross-sector orientation, root-motion payload, line-cycle defects, "
        "relation-wall owners, Ising wall legality, phi4 cumulant signs, "
        "de Moivre branch data, circuit uniformity, and low-frequency resonance "
        "corrections, Asano contraction order, single-far factor identity, "
        "zero-free certificates, SPEC/minorant certificates, far-zero-push "
        "status, A000568 extension shadows, bounded-core resolvent status, "
        "edge-envelope gluing status, resolvent middle-layer payload, "
        "GF coefficient/root-locus/log-derivative payload, fiber-PGF "
        "conditional moments, k=8 reflection-block leakage, fold-adjoint and "
        "odd-coordinate resurrection data, and wide-decoupling budgets when "
        "only the arc sign remains"
    )
    print()


def main() -> None:
    print("HYP-3141 tournament edge-witness tip/tail scout")
    print("=" * 72)
    print_definition()
    print_carriers()
    print_tournament()
    print_signals()
    print_assumption_challenge()
    print("SUMMARY")
    print("-" * 72)
    print(
        "Treating tournament edges as witnesses means keeping the edge as an "
        "information channel.  The witness is not u->v; it is the tail payload, "
        "tip payload, observer-cut orbit, and commutator discharge that let "
        "the LRC predicate survive recursive deletion and extension."
    )


if __name__ == "__main__":
    main()
