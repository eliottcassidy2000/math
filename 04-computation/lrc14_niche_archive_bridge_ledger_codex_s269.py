#!/usr/bin/env python3
"""Niche archive bridge scout for the LRC14 proof frontier.

This is intentionally a synthesis scout: it does not recompute HYP-2963.
It turns older, already-computed archive threads into a ranked list of packet
columns to add to the current HYP-3114/HYP-3115 proof obligations.
"""

from collections import Counter
from itertools import combinations, permutations


AXES = [
    ("preserves_lrc_predicate", 4),
    ("repairs_known_failure", 4),
    ("finite_checkability", 3),
    ("packet_compression", 3),
    ("integrates_hyp3114_hyp3115", 3),
    ("closure_potential", 4),
    ("destroyed_coordinate_control", 3),
    ("implementation_next", 2),
]


CARRIERS = [
    {
        "id": "endpoint_phi_p_activation_circuit",
        "label": "HYP-2108/HYP-2112 endpoint Phi/P activation circuit, sharpened by HYP-3116",
        "anchors": [
            "HYP-2108",
            "HYP-2112",
            "HYP-2790",
            "HYP-2791",
            "HYP-3107",
            "HYP-3116",
            "HYP-3117",
        ],
        "axes": {
            "preserves_lrc_predicate": 4,
            "repairs_known_failure": 4,
            "finite_checkability": 3,
            "packet_compression": 3,
            "integrates_hyp3114_hyp3115": 3,
            "closure_potential": 4,
            "destroyed_coordinate_control": 3,
            "implementation_next": 3,
        },
        "fields": [
            "endpoint_cover_activation_vector",
            "phi_gap_sum",
            "phi_kernel_status",
            "P_max_activation",
            "Phi_gap",
            "P_sign",
            "LMR_terminal_state",
            "magnitude_cocycle",
            "horn_sidecar_closure_status",
            "protected_branch_graph_status",
            "no_naked_bridge_certificate",
            "endpoint_period_numerator_sidecar",
            "finite_address_packet",
            "observer_gluing_certificate",
            "residual_kernel_exclusion_certificate",
            "proof_uniformity_schema",
            "uniform_family_parameter",
            "essential_input_set",
            "minimal_certificate_minterm",
            "legal_sidecar_repair_debt",
            "destroyed_coordinate_discharge",
        ],
        "next": "Treat every archive shortcut as legal only when it emits, lower-bounds, or repairs the exact Phi/P endpoint activation circuit, then passes magnitude, Horn-closure, and protected-branch gates or names the first missing input as residual debt.",
    },
    {
        "id": "normalized_interval_denominator_center",
        "label": "THM-565 normalized interval plus HYP-2866 denominator-center budget",
        "anchors": [
            "HYP-2866",
            "HYP-3088",
            "HYP-3089",
            "HYP-3098",
            "HYP-3114",
            "THM-565",
            "THM-575",
        ],
        "axes": {
            "preserves_lrc_predicate": 3,
            "repairs_known_failure": 3,
            "finite_checkability": 3,
            "packet_compression": 3,
            "integrates_hyp3114_hyp3115": 3,
            "closure_potential": 3,
            "destroyed_coordinate_control": 3,
            "implementation_next": 3,
        },
        "fields": [
            "normalized_interval_floor_status",
            "slow_ruler_component_word",
            "denominator_center_prefix_profile",
            "largest_component_farey_center",
            "all_denominator_grid_bound",
        ],
        "next": "Run the HYP-3114 interval-margin scout in THM-565 slow/ruler coordinates and attach the HYP-2866 denominator-center prefix profile.",
    },
    {
        "id": "et_hensel_fiber_zipper",
        "label": "HYP-3020/HYP-3024 ET+unit discrepancy and Hensel fiber zipper",
        "anchors": [
            "HYP-3020",
            "HYP-3024",
            "HYP-3030",
            "HYP-3033",
            "HYP-3035",
            "HYP-3071",
            "HYP-3098",
        ],
        "axes": {
            "preserves_lrc_predicate": 3,
            "repairs_known_failure": 3,
            "finite_checkability": 3,
            "packet_compression": 3,
            "integrates_hyp3114_hyp3115": 2,
            "closure_potential": 3,
            "destroyed_coordinate_control": 3,
            "implementation_next": 3,
        },
        "fields": [
            "ET_unit_fiber_key",
            "hensel_unit_root_status",
            "zero_root_scale_debt",
            "residue_discrepancy_bin",
            "coarse_to_exact_fiber_zipper",
        ],
        "next": "Join ET+unit and Hensel fields to the observer-gluing ledger before reading interval or root signals as route-separating.",
    },
    {
        "id": "crt_level7_gear",
        "label": "old CRT gear descent for the 14 -> 7 -> 2 composite modulus",
        "anchors": [
            "HYP-2072",
            "HYP-2073",
            "HYP-2081",
            "THM-573",
            "HYP-3094",
            "HYP-3098",
            "HYP-3114",
        ],
        "axes": {
            "preserves_lrc_predicate": 3,
            "repairs_known_failure": 3,
            "finite_checkability": 3,
            "packet_compression": 2,
            "integrates_hyp3114_hyp3115": 2,
            "closure_potential": 3,
            "destroyed_coordinate_control": 3,
            "implementation_next": 2,
        },
        "fields": [
            "crt_gear_state_2x7",
            "crt_c7_lift_status",
            "crt_c2_dyadic_lift_status",
            "apex_stuck_flag",
            "product_preserved_dyadic_lift",
        ],
        "next": "Promote the paper-style c=7,c=2 lift check into a packet field instead of treating it as a standalone computational fallback.",
    },
    {
        "id": "finite_l7_resonance_odometer",
        "label": "HYP-2730 finite L7 resonance atlas plus odometer rowdef",
        "anchors": [
            "HYP-2676",
            "HYP-2730",
            "HYP-2737",
            "HYP-3114",
            "THM-573",
        ],
        "axes": {
            "preserves_lrc_predicate": 2,
            "repairs_known_failure": 2,
            "finite_checkability": 3,
            "packet_compression": 3,
            "integrates_hyp3114_hyp3115": 2,
            "closure_potential": 2,
            "destroyed_coordinate_control": 3,
            "implementation_next": 2,
        },
        "fields": [
            "l7_ratio_resonance_id",
            "exceptional_low_denominator_resonance",
            "odometer_rowdef_word",
            "signed_ET_resonance_tail",
            "low_growth_ruszsa_model_id",
        ],
        "next": "Reuse the finite low-denominator resonance atlas as the exception list for approximation and Bravais peak sidecars.",
    },
    {
        "id": "anti_bohr_boundary_cocycle",
        "label": "anti-Bohr boundary theorem plus lost-coordinate cocycle",
        "anchors": [
            "T345",
            "T490",
            "HYP-2995",
            "HYP-3054",
            "HYP-3102",
            "HYP-3114",
            "HYP-3118",
        ],
        "axes": {
            "preserves_lrc_predicate": 3,
            "repairs_known_failure": 2,
            "finite_checkability": 2,
            "packet_compression": 2,
            "integrates_hyp3114_hyp3115": 2,
            "closure_potential": 2,
            "destroyed_coordinate_control": 3,
            "implementation_next": 2,
        },
        "fields": [
            "destroyed_coordinate_vector",
            "coordinate_resurrection_cover",
            "adjoint_section_status",
            "repair_cover_rank",
            "concept_lattice_intent_id",
            "core_stalk_presence",
            "live_section_type",
            "anti_bohr_boundary_core_id",
            "endpoint_owner_cocycle_id",
            "lost_coordinate_cocycle_id",
            "coordinate_resurrection_repair_cover",
            "observer_extension_cut_payload",
            "boundary_h1_owner_support",
        ],
        "next": "Use AP13's no-positive-component result as an endpoint-owner/cocycle demand, then apply HYP-3118 destroyed-coordinate vectors, cover ranks, adjoint sections, and concept intents before treating it as approximation failure alone.",
    },
    {
        "id": "relation_lattice_ising_wall",
        "label": "HYP-3115 relation lattice plus Ising root-stratum domain walls",
        "anchors": [
            "HYP-3062",
            "HYP-3109",
            "HYP-3116",
            "HYP-3117",
            "HYP-3111",
            "HYP-3112",
            "HYP-3113",
            "HYP-3115",
        ],
        "axes": {
            "preserves_lrc_predicate": 2,
            "repairs_known_failure": 2,
            "finite_checkability": 3,
            "packet_compression": 2,
            "integrates_hyp3114_hyp3115": 3,
            "closure_potential": 2,
            "destroyed_coordinate_control": 2,
            "implementation_next": 2,
        },
        "fields": [
            "short_relation_wall_class",
            "minkowski_successive_minima_proxy",
            "proof_circuit_missing_input_vector",
            "proof_circuit_packet_id",
            "proof_carrier_gate_stack_id",
            "ising_domain_wall_id",
            "root_collision_exit_type",
            "stationary_quintic_residual_class",
        ],
        "next": "Classify HYP-3115 one-swap domain-wall edges by relation wall, HYP-3116 missing-input vector, HYP-3116 proof-carrier gate stack, HYP-3117 proof-circuit packet, root collision, observer-gluing exit, or forbidden debt.",
    },
    {
        "id": "ostrowski_automatic_shadow",
        "label": "Ostrowski/Zeckendorf/automatic-gap shadow with no-adjacent-carry guard",
        "anchors": [
            "HYP-3009",
            "HYP-3010",
            "T494",
            "T495",
            "HYP-3114",
        ],
        "axes": {
            "preserves_lrc_predicate": 2,
            "repairs_known_failure": 1,
            "finite_checkability": 3,
            "packet_compression": 2,
            "integrates_hyp3114_hyp3115": 2,
            "closure_potential": 1,
            "destroyed_coordinate_control": 2,
            "implementation_next": 2,
        },
        "fields": [
            "ostrowski_endpoint_debt_word",
            "zeckendorf_no_adjacent_carry_flag",
            "automatic_gap_language_id",
            "power_lift_guard_status",
            "same_shadow_decoy_flag",
        ],
        "next": "Use automatic/Ostrowski words only after endpoint-debt quotients become path-like; otherwise mark them as decoy shadows.",
    },
    {
        "id": "raw_direct_time_named_constants",
        "label": "raw direct-time intervals and named irrational/transcendental constants",
        "anchors": [
            "HYP-3114",
            "HYP-3088",
            "HYP-3089",
        ],
        "axes": {
            "preserves_lrc_predicate": 1,
            "repairs_known_failure": 0,
            "finite_checkability": 2,
            "packet_compression": 0,
            "integrates_hyp3114_hyp3115": 1,
            "closure_potential": 0,
            "destroyed_coordinate_control": 0,
            "implementation_next": 1,
        },
        "fields": [
            "raw_time_component_length",
            "named_constant_first_hit",
            "continued_fraction_denominator",
        ],
        "next": "Keep as diagnostic output only; it does not close a proof unless upgraded to a normalized interval-margin packet.",
    },
]


TIE_PATH = [carrier["id"] for carrier in CARRIERS]


def weighted_score(carrier):
    axis_values = carrier["axes"]
    return sum(axis_values[name] * weight for name, weight in AXES)


def compare(left, right):
    left_score = weighted_score(left)
    right_score = weighted_score(right)
    if left_score != right_score:
        return left if left_score > right_score else right
    return left if TIE_PATH.index(left["id"]) < TIE_PATH.index(right["id"]) else right


def tournament_edges(carriers):
    edges = {carrier["id"]: set() for carrier in carriers}
    by_id = {carrier["id"]: carrier for carrier in carriers}
    for left_id, right_id in combinations(by_id, 2):
        winner = compare(by_id[left_id], by_id[right_id])
        loser_id = right_id if winner["id"] == left_id else left_id
        edges[winner["id"]].add(loser_id)
    return edges


def directed_three_cycles(edges):
    count = 0
    vertices = list(edges)
    for a, b, c in combinations(vertices, 3):
        triples = [(a, b, c), (a, c, b)]
        for x, y, z in triples:
            if y in edges[x] and z in edges[y] and x in edges[z]:
                count += 1
    return count


def strongly_connected_components(edges):
    vertices = list(edges)
    reverse = {v: set() for v in vertices}
    for source, targets in edges.items():
        for target in targets:
            reverse[target].add(source)

    seen = set()
    order = []

    def dfs_forward(vertex):
        seen.add(vertex)
        for target in edges[vertex]:
            if target not in seen:
                dfs_forward(target)
        order.append(vertex)

    for vertex in vertices:
        if vertex not in seen:
            dfs_forward(vertex)

    seen.clear()
    components = []

    def dfs_reverse(vertex, component):
        seen.add(vertex)
        component.append(vertex)
        for target in reverse[vertex]:
            if target not in seen:
                dfs_reverse(target, component)

    for vertex in reversed(order):
        if vertex not in seen:
            component = []
            dfs_reverse(vertex, component)
            components.append(component)
    return components


def hamiltonian_path_count(edges):
    vertices = list(edges)
    count = 0
    sample = None
    for path in permutations(vertices):
        if all(path[i + 1] in edges[path[i]] for i in range(len(path) - 1)):
            count += 1
            if sample is None:
                sample = path
    return count, sample


def main():
    edges = tournament_edges(CARRIERS)
    outdegrees = {carrier_id: len(targets) for carrier_id, targets in edges.items()}
    score_hist = dict(sorted(Counter(outdegrees.values()).items()))
    sccs = strongly_connected_components(edges)
    h_count, h_sample = hamiltonian_path_count(edges)
    priority = sorted(CARRIERS, key=lambda carrier: (-outdegrees[carrier["id"]], TIE_PATH.index(carrier["id"])))

    print("HYP-3119 NICHE ARCHIVE BRIDGE LEDGER")
    print("claim=older failed scalar routes become useful when read as packet columns feeding the exact endpoint Phi/P activation circuit")
    print("pairwise_observable=which archive carrier better preserves the LRC predicate, repairs a known quotient failure, feeds or lower-bounds Phi/P activation, and can be joined to HYP-3114/HYP-3115 without deleting endpoint or route data")
    print()

    print("TOURNAMENT_ANALYSIS")
    print(f"vertices={len(CARRIERS)}")
    print(f"weighted_scores={{{', '.join(f'{c['id']}:{weighted_score(c)}' for c in CARRIERS)}}}")
    print(f"score_hist={score_hist}")
    print(f"directed_3cycles={directed_three_cycles(edges)}")
    print(f"scc_sizes={[len(component) for component in sccs]}")
    print(f"hamiltonian_path_count={h_count}")
    print("priority_path=" + " -> ".join(carrier["id"] for carrier in priority))
    if h_sample:
        print("tie_hamiltonian_path=" + " -> ".join(h_sample))
    print()

    print("ARCHIVE_BRIDGES")
    for carrier in priority:
        print(f"- {carrier['id']}")
        print(f"  label={carrier['label']}")
        print(f"  anchors={','.join(carrier['anchors'])}")
        print(f"  fields={','.join(carrier['fields'])}")
        print(f"  next={carrier['next']}")
    print()

    print("TWO_BACK_AND_FORTH_TASKS")
    print("1. endpoint_interval_back_and_forth: attach endpoint Phi/P activation fields to HYP-3114 rows, upgrade direct-time intervals to THM-565 slow/ruler coordinates, then test HYP-2866 denominator-center prefix majorization as a lower bound on Phi.")
    print("2. et_hensel_crt_resonance_wall_merge: add ET+unit/Hensel fiber keys, CRT 2x7 gear state, finite L7 resonance ids, HYP-3116 missing-input vectors and proof-carrier gates, HYP-3117 proof-circuit packet ids, and HYP-3118 coordinate-resurrection cover ranks/adjoint sections/concept intents to HYP-3098/HYP-3112 packets, then compare those fields against HYP-3115 Ising domain-wall edges and endpoint Phi/P kernels.")
    print()

    print("ASSUMPTION_CHALLENGE")
    print("vertices_considered=runners,gaps,fixed_sections,section_boundaries,wall_crossings,residue_packets,cover_arcs,fourier_modes,matroid_circuits,proof_obligations")
    print("chosen_vertices=archive proof carriers and packet columns")
    print("preserved_predicate=legal LRC14 witness, observer-gluing, or finite-address route survives the quotient")
    print("destroyed_if_unretained=raw time, endpoint owner, endpoint activation vector, Phi kernel status, LMR wall state, magnitude cocycle, Horn sidecar closure, protected branch status, normalized coordinate, Hensel unit root, CRT gear state, low-denominator resonance id, relation wall class, proof-circuit missing-input vector, proof-circuit packet id, coordinate-resurrection cover rank, adjoint section, concept intent, and exceptional approximant list")


if __name__ == "__main__":
    main()
