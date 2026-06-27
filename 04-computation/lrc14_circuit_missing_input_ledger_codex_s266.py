#!/usr/bin/env python3
"""HYP-3116 circuit lower-bound / missing-input synthesis ledger.

This is a static, auditable route-miner rather than a numeric LRC solver.  It
turns prior repo artifacts into proof-gate records, then scores the records by
how much LRC14 proof payload they retain.  The intent is to make "circuit
complexity" proof-facing: a shortcut is weak exactly when its missing-input
vector contains an undisclosed endpoint, packet, or observer-gluing coordinate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


METRIC_ORDER = (
    "exact_gap",
    "kernel_exclusion",
    "frontier_adjacency",
    "sidecar_retention",
    "quotient_leakage_control",
    "formalizable",
    "global_lift",
    "uniformity",
    "finite_checkability",
)

WEIGHTS = {
    "exact_gap": 5,
    "kernel_exclusion": 5,
    "frontier_adjacency": 4,
    "sidecar_retention": 4,
    "quotient_leakage_control": 4,
    "formalizable": 3,
    "global_lift": 3,
    "uniformity": 3,
    "finite_checkability": 2,
}


@dataclass(frozen=True)
class Gate:
    name: str
    label: str
    kind: str
    anchors: tuple[str, ...]
    input_basis: tuple[str, ...]
    essential_inputs: tuple[str, ...]
    certificate_minterms: tuple[tuple[str, ...], ...]
    missing_inputs: tuple[str, ...]
    reconstructible: tuple[str, ...]
    required_sidecars: tuple[str, ...]
    destroyed_coordinate: str
    terminal_exit: str
    metrics: dict[str, int]

    @property
    def score(self) -> int:
        return sum(self.metrics[k] * WEIGHTS[k] for k in METRIC_ORDER)

    @property
    def vector(self) -> tuple[int, ...]:
        return tuple(self.metrics[k] for k in METRIC_ORDER)


GATES = (
    Gate(
        name="endpoint_phi_sum_gap",
        label="SUM activation Phi(C) = exact gap",
        kind="ReLU/ramp endpoint-cover proof circuit",
        anchors=("HYP-2112", "HYP-2108", "S576"),
        input_basis=(
            "endpoint_cover_activation_vector",
            "component_midpoint_phase",
            "component_length",
            "multiple_of_n_configuration",
            "endpoint_owner_word",
        ),
        essential_inputs=(
            "endpoint_cover_activation_vector",
            "component_midpoint_phase",
            "component_length",
            "multiple_of_n_configuration",
        ),
        certificate_minterms=(
            ("Phi_sum_positive",),
            ("Phi_kernel_empty_on_residual_packet",),
        ),
        missing_inputs=(
            "finite_address_packet",
            "observer_gluing_certificate",
            "residual_kernel_exclusion_certificate",
            "endpoint_period_numerator_sidecar",
        ),
        reconstructible=("P_max_activation_from_phi_vector",),
        required_sidecars=("finite_address_packet", "observer_gluing_certificate"),
        destroyed_coordinate="raw runner labels and raw time after endpoint activations",
        terminal_exit="Phi(C)>0 gives the exact loose gap; kernel debt routes to named residual",
        metrics={
            "exact_gap": 5,
            "kernel_exclusion": 5,
            "frontier_adjacency": 4,
            "sidecar_retention": 4,
            "quotient_leakage_control": 5,
            "formalizable": 4,
            "global_lift": 4,
            "uniformity": 5,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="endpoint_P_max_activation",
        label="MAX activation P(S)>0 endpoint cover gate",
        kind="max-of-ramps cover positivity test",
        anchors=("HYP-2108", "S575"),
        input_basis=(
            "endpoint_cover_activation_vector",
            "midpoint_phase",
            "component_length",
            "endpoint_owner_word",
        ),
        essential_inputs=("midpoint_phase", "component_length", "v_equals_nw_clock"),
        certificate_minterms=(("one_endpoint_activation_positive",),),
        missing_inputs=(
            "residual_kernel_exclusion_certificate",
            "finite_address_packet",
            "observer_gluing_certificate",
            "all_v_clock_family_status",
        ),
        reconstructible=("Phi_sum_lower_bound_when_all_activations_known",),
        required_sidecars=("endpoint_owner_word", "v_equals_nw_clock"),
        destroyed_coordinate="which activation stayed zero across the closed cover",
        terminal_exit="P(S)>0 iff loose for the tested cover circuit",
        metrics={
            "exact_gap": 4,
            "kernel_exclusion": 4,
            "frontier_adjacency": 3,
            "sidecar_retention": 4,
            "quotient_leakage_control": 4,
            "formalizable": 4,
            "global_lift": 3,
            "uniformity": 4,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="lean_proof_frontier_dag",
        label="Lean proof-frontier certificate DAG",
        kind="monotone proof-state circuit",
        anchors=("HYP-3107", "TournamentH7.LRCProofFrontier"),
        input_basis=(
            "q_witness_gate",
            "level7_sieve",
            "dyadic_lift",
            "coverage_extremality",
            "finite_address_packet",
            "observer_gluing_certificate",
            "fine_scale_winding_transfer",
        ),
        essential_inputs=(
            "coverage_extremality",
            "finite_address_packet",
            "observer_gluing_certificate",
            "fine_scale_winding_transfer",
        ),
        certificate_minterms=(
            ("q_witness_gate",),
            ("level7_sieve", "dyadic_lift"),
            ("finite_address_packet", "observer_gluing_certificate"),
        ),
        missing_inputs=(
            "coverage_extremality_producer",
            "node3_effective_peel",
            "finite_ruler_glue",
            "fine_modp_winding_transfer",
        ),
        reconstructible=("terminal_Mreach_from_frontier_coverage",),
        required_sidecars=("finite_address_packet", "observer_gluing_certificate"),
        destroyed_coordinate="none intentionally; this is the current formal interface",
        terminal_exit="lrc14_from_bleeding_edge_frontier style theorem",
        metrics={
            "exact_gap": 2,
            "kernel_exclusion": 4,
            "frontier_adjacency": 5,
            "sidecar_retention": 5,
            "quotient_leakage_control": 5,
            "formalizable": 5,
            "global_lift": 4,
            "uniformity": 5,
            "finite_checkability": 3,
        },
    ),
    Gate(
        name="finite_address_observer_gluing",
        label="finite-address plus observer-gluing exit",
        kind="terminal proof interface",
        anchors=("HYP-3098", "HYP-3107", "HYP-3111"),
        input_basis=(
            "finite_address_packet",
            "observer_gluing_certificate",
            "chart_overlap_certificate",
            "route_legality_status",
        ),
        essential_inputs=(
            "finite_address_packet",
            "observer_gluing_certificate",
            "route_legality_status",
        ),
        certificate_minterms=(("finite_address_packet", "observer_gluing_certificate"),),
        missing_inputs=(
            "endpoint_cover_activation_vector",
            "carrier_certificate",
            "chart_overlap_payload_difference",
        ),
        reconstructible=("packet_sheaf_legal_exit",),
        required_sidecars=("chart_overlap_certificate", "route_legality_status"),
        destroyed_coordinate="raw quotient chart after legal gluing data is retained",
        terminal_exit="observer-gluing theorem or finite-address branch theorem",
        metrics={
            "exact_gap": 1,
            "kernel_exclusion": 4,
            "frontier_adjacency": 5,
            "sidecar_retention": 5,
            "quotient_leakage_control": 5,
            "formalizable": 5,
            "global_lift": 4,
            "uniformity": 5,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="route_state_median_closure",
        label="route-state median/Horn closure circuit",
        kind="Boolean sidecar closure circuit",
        anchors=("HYP-3074", "HYP-3070", "HYP-3069"),
        input_basis=(
            "packet_state_bits",
            "route_state_bits",
            "certificate_bits",
            "sidecar_bits",
            "discharge_bits",
        ),
        essential_inputs=("legal_sidecar_tree", "median_center_signature", "discharge_atom"),
        certificate_minterms=(("unique_legal_median_center", "specific_discharge_atom"),),
        missing_inputs=(
            "non_tautological_packet_witness_floor",
            "endpoint_cover_activation_vector",
            "residual_debt_basis_atom",
        ),
        reconstructible=("raw_route_center_from_sidecar_tree",),
        required_sidecars=("median_center_signature", "specific_discharge_atom"),
        destroyed_coordinate="raw route label clique center",
        terminal_exit="center-control coverage implies LRC14 only after real packet witnesses",
        metrics={
            "exact_gap": 1,
            "kernel_exclusion": 3,
            "frontier_adjacency": 5,
            "sidecar_retention": 5,
            "quotient_leakage_control": 5,
            "formalizable": 4,
            "global_lift": 3,
            "uniformity": 5,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="boolean_mobius_low_depth_cut",
        label="low-depth signed Boolean-Mobius cut",
        kind="threshold circuit over missed-sector atom types",
        anchors=("HYP-2791", "HYP-2744"),
        input_basis=("T1_atom", "T2sep_atom", "T2adj_atom", "AP_orbit_status"),
        essential_inputs=("T1_atom", "T2sep_atom", "T2adj_atom"),
        certificate_minterms=(("21*T1+57*T2sep+2*T2adj_positive",),),
        missing_inputs=(
            "endpoint_period_numerator_sidecar",
            "speed_owner_sidecar",
            "global_lift_from_k8_bank",
            "height_phase_sidecar",
        ),
        reconstructible=("bounded_bank_cut_margin",),
        required_sidecars=("AP_orbit_status", "endpoint_period_numerator_sidecar"),
        destroyed_coordinate="height, speed ownership, and far-element phase",
        terminal_exit="sharp k=8 bank cut; global lift still open",
        metrics={
            "exact_gap": 3,
            "kernel_exclusion": 3,
            "frontier_adjacency": 3,
            "sidecar_retention": 2,
            "quotient_leakage_control": 2,
            "formalizable": 5,
            "global_lift": 2,
            "uniformity": 4,
            "finite_checkability": 5,
        },
    ),
    Gate(
        name="full_boolean_mobius_hierarchy",
        label="full Boolean-Mobius/Delsarte atom law",
        kind="exact Boolean atom basis",
        anchors=("HYP-2744",),
        input_basis=("all_missed_sector_atoms", "containment_law", "mobius_inversion"),
        essential_inputs=("all_missed_sector_atoms", "mobius_inversion"),
        certificate_minterms=(("small_signed_functional_found",),),
        missing_inputs=(
            "small_signed_functional",
            "endpoint_phase_connection",
            "global_packet_lift",
        ),
        reconstructible=("exact_atom_law_from_containment_law",),
        required_sidecars=("dihedral_type", "packet_key"),
        destroyed_coordinate="speed ownership and endpoint phase if only atom totals survive",
        terminal_exit="basis for finite signed cuts, not a terminal scalar",
        metrics={
            "exact_gap": 3,
            "kernel_exclusion": 2,
            "frontier_adjacency": 2,
            "sidecar_retention": 3,
            "quotient_leakage_control": 3,
            "formalizable": 4,
            "global_lift": 2,
            "uniformity": 3,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="fourier_toeplitz_psd_dual",
        label="Fourier/Toeplitz PSD dual certificate",
        kind="linear-algebraic dual circuit",
        anchors=("HYP-2974", "HYP-2981", "LTI-056"),
        input_basis=("Toeplitz_matrix", "Fejer_vector", "packet_key", "divisor_curried_coefficients"),
        essential_inputs=("Toeplitz_matrix", "Fejer_vector", "packet_key"),
        certificate_minterms=(("negative_quadratic_form",),),
        missing_inputs=(
            "interval_arithmetic_certificate",
            "endpoint_owner_decoding",
            "finite_address_packet",
        ),
        reconstructible=("dual_trig_square_multiplier",),
        required_sidecars=("packet_key", "endpoint_owner_word"),
        destroyed_coordinate="which endpoint owner explains the negative mode",
        terminal_exit="PSD failure forces positive safe interval unless boundary-owned",
        metrics={
            "exact_gap": 3,
            "kernel_exclusion": 3,
            "frontier_adjacency": 3,
            "sidecar_retention": 4,
            "quotient_leakage_control": 4,
            "formalizable": 3,
            "global_lift": 3,
            "uniformity": 4,
            "finite_checkability": 3,
        },
    ),
    Gate(
        name="haar_zipper_cocycle",
        label="Haar zipper fixed-margin cocycle",
        kind="local switch / cocycle circuit",
        anchors=("HYP-2991", "HYP-2989", "HYP-2997"),
        input_basis=("row_margin", "column_margin", "haar_zeta", "owner_strip_label"),
        essential_inputs=("haar_zeta", "owner_strip_label"),
        certificate_minterms=(("fixed_margin_plus_zeta_separates_fiber",),),
        missing_inputs=(
            "endpoint_cover_activation_vector",
            "packet_route_key",
            "dual_annihilator_status",
        ),
        reconstructible=("local_fixed_margin_fiber_after_zeta",),
        required_sidecars=("haar_zeta", "owner_strip_label"),
        destroyed_coordinate="mixed Haar sign if only margins are kept",
        terminal_exit="quotient legality or first cocycle debt",
        metrics={
            "exact_gap": 2,
            "kernel_exclusion": 3,
            "frontier_adjacency": 3,
            "sidecar_retention": 4,
            "quotient_leakage_control": 5,
            "formalizable": 4,
            "global_lift": 3,
            "uniformity": 4,
            "finite_checkability": 5,
        },
    ),
    Gate(
        name="lee_yang_root_ear_payload",
        label="Lee-Yang root curve plus ear payload",
        kind="PGF/root-motion proof sidecar",
        anchors=("HYP-3112", "HYP-3109", "HYP-3108"),
        input_basis=(
            "miss_count_pgf_coefficients",
            "root_multiset",
            "ear_payload_A_vector",
            "root_motion_reconstruction_status",
        ),
        essential_inputs=("root_multiset", "ear_payload_A_vector"),
        certificate_minterms=(("root_wall_classified", "ear_payload_reconstructs_extension"),),
        missing_inputs=(
            "endpoint_owner_word",
            "finite_address_packet",
            "observer_gluing_certificate",
            "root_collision_wall_label",
        ),
        reconstructible=("q_full_from_q_base_and_A_vector",),
        required_sidecars=("nested_ear_status", "root_motion_reconstruction_status"),
        destroyed_coordinate="sector incidence and next-extension root motion if A is dropped",
        terminal_exit="root wall, legal ear discharge, or named root-collision debt",
        metrics={
            "exact_gap": 2,
            "kernel_exclusion": 2,
            "frontier_adjacency": 3,
            "sidecar_retention": 4,
            "quotient_leakage_control": 4,
            "formalizable": 3,
            "global_lift": 3,
            "uniformity": 3,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="minkowski_circuit_ising_bridge",
        label="Minkowski/circuit/Ising/De Moivre cut-payload bridge",
        kind="source-carrier sidecar bundle",
        anchors=("HYP-3111", "HYP-3115"),
        input_basis=(
            "q_body_inequality_word",
            "proof_circuit_missing_input_vector",
            "ising_zero_arc_signature",
            "demoivre_branch_orbit_word",
        ),
        essential_inputs=("proof_circuit_missing_input_vector", "q_body_inequality_word"),
        certificate_minterms=(("finite_address_packet", "observer_gluing_certificate", "carrier_sidecar"),),
        missing_inputs=(
            "endpoint_cover_activation_vector",
            "proof_uniformity_schema",
            "root_ear_lattice_legal_exit",
        ),
        reconstructible=("duodecimal_preserve_destroy_handoff_audit",),
        required_sidecars=(
            "q_body_inequality_word",
            "proof_circuit_missing_input_vector",
            "ising_zero_arc_signature",
            "demoivre_branch_orbit_word",
        ),
        destroyed_coordinate="outside-theorem route legality if cut payload is omitted",
        terminal_exit="pressure gauge; proof exit only through finite address and observer gluing",
        metrics={
            "exact_gap": 1,
            "kernel_exclusion": 2,
            "frontier_adjacency": 3,
            "sidecar_retention": 3,
            "quotient_leakage_control": 3,
            "formalizable": 2,
            "global_lift": 2,
            "uniformity": 2,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="endpoint_period_warning",
        label="endpoint-period numerator warning",
        kind="guardrail against scalar Boolean transfer",
        anchors=("HYP-2790",),
        input_basis=("endpoint_period_numerator", "q0_cover_slack", "type_slack"),
        essential_inputs=("endpoint_period_numerator",),
        certificate_minterms=(("endpoint_period_numerator_retained",),),
        missing_inputs=(
            "endpoint_period_numerator_sidecar",
            "endpoint_cover_activation_vector",
            "speed_owner_sidecar",
        ),
        reconstructible=("refutation_of_scalar_transfer",),
        required_sidecars=("endpoint_period_numerator_sidecar",),
        destroyed_coordinate="phase numerator when Boolean type is isolated",
        terminal_exit="guardrail: Boolean/type scalars do not transfer alone",
        metrics={
            "exact_gap": 1,
            "kernel_exclusion": 2,
            "frontier_adjacency": 2,
            "sidecar_retention": 3,
            "quotient_leakage_control": 4,
            "formalizable": 3,
            "global_lift": 2,
            "uniformity": 3,
            "finite_checkability": 4,
        },
    ),
    Gate(
        name="finite_bank_apex7_literal",
        label="finite-bank apex7_error <= 5 literal",
        kind="nonuniform fitted classifier",
        anchors=("HYP-3115",),
        input_basis=("apex7_error",),
        essential_inputs=("apex7_error",),
        certificate_minterms=(("apex7_error_le_5_on_anchored_bank",),),
        missing_inputs=(
            "proof_uniformity_schema",
            "endpoint_cover_activation_vector",
            "finite_address_packet",
            "observer_gluing_certificate",
            "train_test_family_split",
        ),
        reconstructible=("bounded_bank_max_p0_signal",),
        required_sidecars=("proof_uniformity_schema",),
        destroyed_coordinate="everything outside the fitted anchored bank",
        terminal_exit="warning signal only; not a proof circuit",
        metrics={
            "exact_gap": 1,
            "kernel_exclusion": 1,
            "frontier_adjacency": 1,
            "sidecar_retention": 1,
            "quotient_leakage_control": 0,
            "formalizable": 1,
            "global_lift": 0,
            "uniformity": 0,
            "finite_checkability": 5,
        },
    ),
)


TIE_PATH = tuple(g.name for g in GATES)
TIE_RANK = {name: i for i, name in enumerate(TIE_PATH)}


def beats(a: Gate, b: Gate) -> bool:
    if a.score != b.score:
        return a.score > b.score
    if a.vector != b.vector:
        return a.vector > b.vector
    return TIE_RANK[a.name] < TIE_RANK[b.name]


def adjacency(gates: Iterable[Gate]) -> dict[str, set[str]]:
    nodes = list(gates)
    adj = {g.name: set() for g in nodes}
    for a, b in combinations(nodes, 2):
        if beats(a, b):
            adj[a.name].add(b.name)
        else:
            adj[b.name].add(a.name)
    return adj


def score_hist(adj: dict[str, set[str]]) -> dict[int, int]:
    return dict(sorted(Counter(len(v) for v in adj.values()).items()))


def directed_3cycles(adj: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles = []
    for a, b, c in combinations(adj, 3):
        edges = ((b in adj[a]), (c in adj[b]), (a in adj[c]))
        if all(edges) or not any(edges):
            cycles.append((a, b, c))
    return cycles


def strongly_connected_components(adj: dict[str, set[str]]) -> list[list[str]]:
    nodes = list(adj)
    rev = {n: set() for n in nodes}
    for u, outs in adj.items():
        for v in outs:
            rev[v].add(u)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(u: str) -> None:
        seen.add(u)
        for v in adj[u]:
            if v not in seen:
                dfs(v)
        order.append(u)

    for n in nodes:
        if n not in seen:
            dfs(n)

    comps: list[list[str]] = []
    seen.clear()

    def rdfs(u: str, comp: list[str]) -> None:
        seen.add(u)
        comp.append(u)
        for v in rev[u]:
            if v not in seen:
                rdfs(v, comp)

    for n in reversed(order):
        if n not in seen:
            comp: list[str] = []
            rdfs(n, comp)
            comps.append(sorted(comp, key=lambda x: TIE_RANK[x]))
    return comps


def hamiltonian_path_count(adj: dict[str, set[str]]) -> int:
    nodes = list(adj)
    n = len(nodes)
    index = {name: i for i, name in enumerate(nodes)}
    dp: dict[tuple[int, int], int] = {}
    for name in nodes:
        dp[(1 << index[name], index[name])] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if count == 0:
                continue
            u = nodes[last]
            for v in adj[u]:
                j = index[v]
                if not (mask >> j) & 1:
                    dp[(mask | (1 << j), j)] = dp.get((mask | (1 << j), j), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def priority_path(adj: dict[str, set[str]]) -> list[str]:
    return sorted(adj, key=lambda name: (-len(adj[name]), -by_name[name].score, TIE_RANK[name]))


def edge_flips_against_novelty_gauge(gates: Iterable[Gate], adj: dict[str, set[str]]) -> int:
    # Novelty gauge intentionally overvalues fashionable outside analogies and
    # finite-bank checkability, to expose how many edges the proof gauge flips.
    novelty_order = sorted(
        gates,
        key=lambda g: (
            -("HYP-3115" in g.anchors or "HYP-3111" in g.anchors),
            -g.metrics["finite_checkability"],
            g.metrics["frontier_adjacency"],
            TIE_RANK[g.name],
        ),
    )
    novelty_rank = {g.name: i for i, g in enumerate(novelty_order)}
    flips = 0
    for a, b in combinations([g.name for g in gates], 2):
        proof_a_beats = b in adj[a]
        novelty_a_beats = novelty_rank[a] < novelty_rank[b]
        if proof_a_beats != novelty_a_beats:
            flips += 1
    return flips


def format_tuple(values: Iterable[str]) -> str:
    vals = list(values)
    if not vals:
        return "none"
    return ", ".join(vals)


by_name = {g.name: g for g in GATES}


def main() -> None:
    adj = adjacency(GATES)
    cycles = directed_3cycles(adj)
    comps = strongly_connected_components(adj)
    path = priority_path(adj)

    missing_counter = Counter()
    essential_counter = Counter()
    anchor_counter = Counter()
    for gate in GATES:
        missing_counter.update(gate.missing_inputs)
        essential_counter.update(gate.essential_inputs)
        anchor_counter.update(gate.anchors)

    missing_by_gate = defaultdict(list)
    for gate in GATES:
        for item in gate.missing_inputs:
            missing_by_gate[item].append(gate.name)

    print("HYP-3116 circuit lower-bound / missing-input ledger")
    print("=====================================================")
    print()
    print(f"proof gates mined: {len(GATES)}")
    print(f"metric keys: {', '.join(METRIC_ORDER)}")
    print(
        "pairwise observable: weighted proof-payload retention, then exact metric vector, "
        "then declared tie Hamiltonian path"
    )
    print(
        "assumption challenge: vertices are proof carriers and sidecar ledgers, "
        "not runners, arcs, Boolean gates, roots, or raw scalar values"
    )
    print()

    print("Core synthesis")
    print("--------------")
    print(
        "The strongest old circuit bridge is HYP-2108/HYP-2112: P(S) is a MAX "
        "activation gate and Phi(C) is the SUM activation gate equal to the exact gap."
    )
    print(
        "Thus the LRC14 circuit task is a kernel-exclusion problem for an "
        "endpoint-cover activation circuit, not a generic P/poly lower-bound slogan."
    )
    print(
        "HYP-2791 supplies a small threshold-style Boolean cut, but HYP-2790 says "
        "the endpoint-period numerator and ownership sidecars must travel with it."
    )
    print()

    print("Gate ledger")
    print("-----------")
    for gate in sorted(GATES, key=lambda g: (-g.score, TIE_RANK[g.name])):
        print(f"{gate.name}: {gate.label}")
        print(f"  anchors: {format_tuple(gate.anchors)}")
        print(f"  kind: {gate.kind}")
        print(f"  score: {gate.score}; vector={gate.vector}")
        print(f"  essential: {format_tuple(gate.essential_inputs)}")
        print(f"  missing: {format_tuple(gate.missing_inputs)}")
        print(
            "  minterms: "
            + "; ".join(" & ".join(term) for term in gate.certificate_minterms)
        )
        print(f"  terminal: {gate.terminal_exit}")
    print()

    print("Top missing inputs")
    print("------------------")
    for item, count in missing_counter.most_common():
        gate_names = ", ".join(missing_by_gate[item])
        print(f"{item}: {count} gates -> {gate_names}")
    print()

    print("Top essential inputs")
    print("--------------------")
    for item, count in essential_counter.most_common(18):
        print(f"{item}: {count}")
    print()

    print("Most reused anchors")
    print("-------------------")
    for item, count in anchor_counter.most_common(18):
        print(f"{item}: {count}")
    print()

    print("Tournament Analysis")
    print("-------------------")
    print(f"score_hist = {score_hist(adj)}")
    print(f"directed_3cycles = {len(cycles)}")
    if cycles:
        for cyc in cycles[:8]:
            print("  cycle:", " -> ".join(cyc))
    print(f"scc_sizes = {[len(c) for c in comps]}")
    print(f"sccs = {comps}")
    print(f"hamiltonian_path_count = {hamiltonian_path_count(adj)}")
    print(f"edge_flips_against_novelty_gauge = {edge_flips_against_novelty_gauge(GATES, adj)}")
    print("priority_path =")
    for i, name in enumerate(path, 1):
        print(f"  {i}. {name}")
    print()

    print("Circuit lower-bound readout")
    print("---------------------------")
    print(
        "Any low-depth LRC14 shortcut omitting endpoint_cover_activation_vector has "
        "not reached the HYP-2108/HYP-2112 gap circuit."
    )
    print(
        "Any shortcut omitting both finite_address_packet and observer_gluing_certificate "
        "does not feed the HYP-3107 proof frontier."
    )
    print(
        "Any Boolean/type shortcut omitting endpoint_period_numerator_sidecar repeats "
        "the HYP-2790 scalar-transfer failure mode."
    )
    print(
        "The proof-facing lower bound is therefore: every terminal circuit must carry "
        "Phi/Phi-kernel data plus packet/gluing sidecars, or explicitly route the "
        "first missing input to named residual debt."
    )
    print()

    print("Proposed packet fields")
    print("----------------------")
    fields = (
        "endpoint_cover_activation_vector",
        "phi_gap_sum",
        "phi_kernel_status",
        "P_max_activation",
        "simultaneous_resonance_winding_word",
        "boolean_mobius_low_depth_cut",
        "endpoint_period_numerator_sidecar",
        "proof_circuit_missing_input_vector",
        "proof_uniformity_schema",
        "finite_address_packet",
        "observer_gluing_certificate",
        "walsh_degree_support_profile",
    )
    for field in fields:
        print(f"- {field}")


if __name__ == "__main__":
    main()
