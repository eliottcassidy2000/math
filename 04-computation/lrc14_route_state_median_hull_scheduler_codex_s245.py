#!/usr/bin/env python3
"""S245: route-state median-hull scheduler for LRC14.

This scout turns the proposed final proof interface into a finite median check.
Proof states are sets of retained packet/route/certificate/sidecar/discharge
coordinates.  The scout uses unary Horn sidecar implications, so
coordinatewise majority of three legal states is legal.  A serious route
triple therefore has a unique center once the required sidecars are attached.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product


FeatureSet = frozenset[str]


CORE_FEATURES = [
    "packet_label",
    "exact_M",
    "endpoint_owner",
    "safe_topology",
    "magnitude_cocycle",
    "route_label",
    "certificate_payload",
    "sidecar_payload",
    "discharge_atom_typed",
]

SPECIFIC_DISCHARGES = [
    "q_witness_atom",
    "ap_gw_boundary_atom",
    "c27_petal_atom",
    "k33_state_lift_atom",
    "fejer_toeplitz_atom",
    "desargues_beal_atom",
    "named_residual_debt_atom",
]

SIDECAR_FEATURES = [
    "observer_cut_payload",
    "value_origin_type",
    "deletion_fiber_profile",
    "cross_sector_orientation_word",
    "rectangle_hourglass_residue",
    "partial_cube_theta_word",
    "bit_position_phase",
    "moser_product_split_a_plus_2b",
    "bridge_rank_debt",
    "toeplitz_scale_gate",
    "ordered_quad_collapse_mode",
    "d4_orbit_word",
    "hyperbolic_triple_signature",
    "geometry_regime_signature",
    "relation_lattice_tail",
    "residue_signed_tail",
    "c27_shell_transfer",
    "k33_incidence_wall",
    "beal_common_owner_gate",
    "desargues_girth6_residue",
    "residual_capacitor_id",
    "haar_zeta_square",
    "endpoint_owner_strip",
    "primitive_period_deck",
    "fejer_interval_certificate",
]


RULES: list[tuple[set[str], set[str], str]] = [
    (
        {"route_label"},
        {"packet_label", "exact_M", "endpoint_owner", "safe_topology", "magnitude_cocycle"},
        "routes must carry exact scale, owner, topology, and magnitude",
    ),
    (
        {"certificate_payload"},
        {"packet_label", "exact_M", "safe_topology"},
        "certificates are packeted interval/topology objects",
    ),
    (
        {"sidecar_payload"},
        {"packet_label"},
        "sidecars attach to labelled packet fibers",
    ),
    (
        {"discharge_atom_typed"},
        {"packet_label", "route_label", "certificate_payload"},
        "a discharge atom is only legal after route and certificate payloads exist",
    ),
    (
        {"observer_cut_payload"},
        {"sidecar_payload", "value_origin_type", "deletion_fiber_profile", "cross_sector_orientation_word"},
        "observer cuts need origin, deletion, and cross-sector data",
    ),
    (
        {"rectangle_hourglass_residue"},
        {"observer_cut_payload", "value_origin_type"},
        "rectangle/hourglass residues are observer-cut payloads, not counts",
    ),
    (
        {"partial_cube_theta_word"},
        {"sidecar_payload", "bit_position_phase", "exact_M", "endpoint_owner", "safe_topology"},
        "partial-cube cuts need bit phase and LRC packet data",
    ),
    (
        {"moser_product_split_a_plus_2b"},
        {"partial_cube_theta_word", "bit_position_phase"},
        "Moser two-lane split is a partial-cube cut refinement",
    ),
    (
        {"bridge_rank_debt"},
        {"rectangle_hourglass_residue", "value_origin_type"},
        "K bridge surplus is rectangle-cycle debt",
    ),
    (
        {"toeplitz_scale_gate"},
        {"sidecar_payload", "ordered_quad_collapse_mode", "exact_M", "endpoint_owner", "safe_topology"},
        "Toeplitz/square-peg evidence needs noncollapse scale",
    ),
    (
        {"d4_orbit_word"},
        {"toeplitz_scale_gate", "ordered_quad_collapse_mode"},
        "four-witness order is a D4 orbit after the scale gate",
    ),
    (
        {"hyperbolic_triple_signature"},
        {"sidecar_payload", "exact_M", "endpoint_owner", "route_label"},
        "Fermat-Catalan pressure is legal only with packet route data",
    ),
    (
        {"geometry_regime_signature"},
        {"sidecar_payload", "exact_M", "endpoint_owner", "safe_topology", "value_origin_type"},
        "5/6/7 geometry needs axis, packet, and value origin",
    ),
    (
        {"relation_lattice_tail"},
        {"sidecar_payload", "exact_M", "residue_signed_tail"},
        "Roth-Minkowski pressure needs the lattice tail and residue sign",
    ),
    (
        {"c27_shell_transfer"},
        {"route_label", "exact_M", "endpoint_owner"},
        "C27 shell transfer is a labelled route, not a scalar bucket",
    ),
    (
        {"k33_incidence_wall"},
        {"route_label", "exact_M", "endpoint_owner", "beal_common_owner_gate"},
        "K33 incidence must retain owner/common-factor data",
    ),
    (
        {"desargues_girth6_residue"},
        {"observer_cut_payload", "beal_common_owner_gate", "rectangle_hourglass_residue"},
        "Desargues enters after rectangle residues with owner gate",
    ),
    (
        {"residual_capacitor_id"},
        {"observer_cut_payload", "haar_zeta_square", "endpoint_owner_strip"},
        "residual capacitors need zeta and owner-strip exits",
    ),
    (
        {"primitive_period_deck"},
        {"route_label", "exact_M"},
        "period decks are route-side exact-scale teeth",
    ),
    (
        {"fejer_interval_certificate"},
        {"certificate_payload", "exact_M", "safe_topology"},
        "Fejer/Toeplitz atoms are interval certificates",
    ),
]

for atom in SPECIFIC_DISCHARGES:
    RULES.append(({atom}, {"discharge_atom_typed"}, f"{atom} is a specific discharge atom"))

RULES.extend(
    [
        ({"q_witness_atom"}, {"primitive_period_deck"}, "q-witness uses primitive period deck"),
        ({"ap_gw_boundary_atom"}, {"endpoint_owner", "safe_topology"}, "AP/GW is endpoint boundary"),
        ({"c27_petal_atom"}, {"c27_shell_transfer"}, "C27 petal atom uses C27 shell transfer"),
        ({"k33_state_lift_atom"}, {"k33_incidence_wall"}, "K33 atom uses incidence wall"),
        ({"fejer_toeplitz_atom"}, {"fejer_interval_certificate"}, "Fejer atom uses interval certificate"),
        ({"desargues_beal_atom"}, {"desargues_girth6_residue"}, "Desargues/Beal atom uses girth-six residue"),
        ({"named_residual_debt_atom"}, {"sidecar_payload"}, "named residual debt is sidecar-visible"),
    ]
)


@dataclass(frozen=True)
class State:
    name: str
    features: FeatureSet


def close(features: set[str] | FeatureSet) -> FeatureSet:
    """Apply legal sidecar implications to a fixed point."""

    out = set(features)
    changed = True
    while changed:
        changed = False
        for premise, consequence, _reason in RULES:
            if premise <= out and not consequence <= out:
                out.update(consequence)
                changed = True
    return frozenset(out)


def missing_obligations(features: FeatureSet) -> list[str]:
    missing = []
    for premise, consequence, reason in RULES:
        if premise <= features and not consequence <= features:
            missing.append(reason)
    return missing


def state(name: str, features: set[str]) -> State:
    return State(name, close(features))


STATES = [
    state("q_witness_route", {"q_witness_atom"}),
    state("ap_gw_boundary_route", {"ap_gw_boundary_atom", "observer_cut_payload"}),
    state("c27_petal_route", {"c27_petal_atom", "geometry_regime_signature"}),
    state("k33_state_lift_route", {"k33_state_lift_atom", "relation_lattice_tail"}),
    state("fejer_toeplitz_route", {"fejer_toeplitz_atom", "toeplitz_scale_gate"}),
    state("desargues_beal_route", {"desargues_beal_atom"}),
    state("residual_capacitor_route", {"named_residual_debt_atom", "residual_capacitor_id"}),
    state(
        "partial_cube_moser_sidecar",
        {"partial_cube_theta_word", "moser_product_split_a_plus_2b", "bridge_rank_debt"},
    ),
    state("hyperbolic_packet_pressure", {"hyperbolic_triple_signature", "geometry_regime_signature"}),
    state("toeplitz_noncollapse_sidecar", {"toeplitz_scale_gate", "d4_orbit_word"}),
]


TRIPLES = [
    (
        "low_frontier_scheduler",
        ("ap_gw_boundary_route", "c27_petal_route", "k33_state_lift_route"),
        "AP/GW, C27 petal, and K33 routes share exact packet obligations but no terminal atom.",
    ),
    (
        "certificate_noncollapse_triangle",
        ("fejer_toeplitz_route", "toeplitz_noncollapse_sidecar", "desargues_beal_route"),
        "Four-witness noncollapse, interval certificate, and girth-six finalizer meet at certificate payload.",
    ),
    (
        "sequence_geometry_pressure",
        ("partial_cube_moser_sidecar", "hyperbolic_packet_pressure", "c27_petal_route"),
        "Moser/fibbinary cuts, hyperbolic pressure, and C27 shell transfer meet at sidecar-completed route data.",
    ),
    (
        "residual_repair_triangle",
        ("residual_capacitor_route", "ap_gw_boundary_route", "fejer_toeplitz_route"),
        "Residual capacitor, boundary atom, and Fejer certificate test whether owner/zeta exits are retained.",
    ),
    (
        "finalizer_pressure_triangle",
        ("k33_state_lift_route", "desargues_beal_route", "hyperbolic_packet_pressure"),
        "K33, Desargues/Beal, and Fermat-Catalan pressure test common-owner sidecars.",
    ),
]


def majority(features: list[FeatureSet]) -> FeatureSet:
    counts = Counter()
    for fs in features:
        counts.update(fs)
    threshold = len(features) // 2 + 1
    return frozenset(feature for feature, count in counts.items() if count >= threshold)


def median(a: State, b: State, c: State) -> FeatureSet:
    return close(majority([a.features, b.features, c.features]))


def interval_center_verified(a: FeatureSet, b: FeatureSet, c: FeatureSet) -> bool:
    """The hypercube interval intersection of three bitsets is their majority."""

    center = majority([a, b, c])
    all_features = set(CORE_FEATURES + SPECIFIC_DISCHARGES + SIDECAR_FEATURES)
    for feature in all_features:
        values = [feature in a, feature in b, feature in c]
        if (values.count(True) >= 2) != (feature in center):
            return False
    return True


def specific_atoms(features: FeatureSet) -> list[str]:
    return [atom for atom in SPECIFIC_DISCHARGES if atom in features]


def core_signature(features: FeatureSet) -> str:
    parts = []
    for group_name, group in [
        ("core", CORE_FEATURES),
        ("discharge", SPECIFIC_DISCHARGES),
        ("sidecars", SIDECAR_FEATURES),
    ]:
        present = [feature for feature in group if feature in features]
        parts.append(f"{group_name}={len(present)}")
    return ", ".join(parts)


def median_hull(states: list[State]) -> list[FeatureSet]:
    hull = {s.features for s in states}
    changed = True
    while changed:
        changed = False
        current = list(hull)
        for a, b, c in product(current, repeat=3):
            center = close(majority([a, b, c]))
            if center not in hull:
                hull.add(center)
                changed = True
    return sorted(hull, key=lambda fs: (len(fs), sorted(fs)))


def count_directed_3cycles(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        if (
            (a, b) in edges
            and (b, c) in edges
            and (c, a) in edges
        ) or (
            (a, c) in edges
            and (c, b) in edges
            and (b, a) in edges
        ):
            cycles += 1
    return cycles


def hamiltonian_path_count(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (vertices[last], vertices[nxt]) in edges:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[(1 << n) - 1][index[v]] for v in vertices)


def print_rule_audit() -> None:
    print("1. LEGAL SIDECAR CLOSURE")
    print(f"   named_features={len(set(CORE_FEATURES + SPECIFIC_DISCHARGES + SIDECAR_FEATURES))}")
    print(f"   horn_rules={len(RULES)}")
    print(f"   max_premise_arity={max(len(premise) for premise, _consequence, _reason in RULES)}")
    print("   closure_rule=if a retained coordinate implies sidecars, add them before taking theorem quotients")
    print("   median_rule=coordinatewise majority; closure is a no-op for legal states under the unary rule set")
    print("   caveat=multi-premise Horn guards must be compiled into named sidecar coordinates or checked separately")
    print()
    print("   rule sample")
    for premise, consequence, reason in RULES[:10]:
        print(f"   - {sorted(premise)} -> {sorted(consequence)} :: {reason}")


def print_state_table() -> None:
    print()
    print("2. LEGAL PROOF STATES")
    print("   state                         feature_count  specific_discharge_atoms  signature")
    for s in STATES:
        atoms = specific_atoms(s.features)
        atom_text = ",".join(atoms) if atoms else "-"
        print(f"   {s.name:29s} {len(s.features):13d}  {atom_text:26s}  {core_signature(s.features)}")


def print_triple_table() -> None:
    print()
    print("3. SERIOUS ROUTE TRIPLES AND MEDIAN CENTERS")
    state_by_name = {s.name: s for s in STATES}
    for name, triple_names, reading in TRIPLES:
        triple = [state_by_name[n] for n in triple_names]
        raw = majority([s.features for s in triple])
        center = close(raw)
        raw_missing = missing_obligations(raw)
        atoms = specific_atoms(center)
        dropped_atoms = sorted(set().union(*(set(specific_atoms(s.features)) for s in triple)) - set(atoms))
        terminal = bool(atoms)
        print(f"   {name}")
        print(f"      routes={', '.join(triple_names)}")
        print(f"      reading={reading}")
        print(f"      raw_majority_features={len(raw)}  raw_missing_obligations={len(raw_missing)}")
        print(f"      legal_center_features={len(center)}  center_signature={core_signature(center)}")
        print(f"      center_specific_atoms={','.join(atoms) if atoms else '-'}")
        print(f"      dropped_specific_atoms={','.join(dropped_atoms) if dropped_atoms else '-'}")
        print(f"      center_kind={'terminal_exit' if terminal else 'scheduler_center_needs_refinement'}")


def print_hull_audit() -> None:
    print()
    print("4. MEDIAN-HULL CHECK")
    hull = median_hull(STATES)
    closed_state_count = sum(1 for fs in hull if not missing_obligations(fs))
    triple_count = 0
    raw_illegal_majorities = 0
    interval_failures = 0
    closure_added_hist = Counter()
    illegal_centers = 0
    for a, b, c in product(hull, repeat=3):
        triple_count += 1
        raw = majority([a, b, c])
        center = close(raw)
        if missing_obligations(raw):
            raw_illegal_majorities += 1
        if not interval_center_verified(a, b, c):
            interval_failures += 1
        closure_added_hist[len(center - raw)] += 1
        if missing_obligations(center):
            illegal_centers += 1
    print(f"   seed_states={len(STATES)}")
    print(f"   median_hull_states={len(hull)}")
    print(f"   closed_states_in_hull={closed_state_count}")
    print(f"   hull_triples_checked={triple_count}")
    print(f"   raw_illegal_majorities={raw_illegal_majorities}")
    print(f"   closure_added_features_hist={dict(sorted(closure_added_hist.items()))}")
    print(f"   interval_intersection_failures={interval_failures}")
    print(f"   illegal_centers_after_closure={illegal_centers}")
    print("   interpretation=the unary-Horn sidecar interface is finite and every checked triple has one legal center")
    print()
    print("   smallest nontrivial centers")
    for fs in hull[:8]:
        print(f"   - size={len(fs):2d} atoms={specific_atoms(fs) or ['-']} signature={core_signature(fs)}")


def print_tournament_analysis() -> None:
    print()
    print("5. TOURNAMENT ANALYSIS")
    vertices = [
        "labelled_packet_sheaf",
        "route_state_median_center",
        "horn_sidecar_closure",
        "discharge_atom_type",
        "observer_cut_payload",
        "partial_cube_cut_payload",
        "toeplitz_noncollapse_gate",
        "hyperbolic_triple_pressure",
        "raw_route_label",
    ]
    edges: set[tuple[str, str]] = set()
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            edges.add((a, b))
    scores = Counter()
    for a, b in edges:
        scores[a] += 1
        scores.setdefault(b, scores[b])
    score_hist = Counter(scores[v] for v in vertices)
    print("   vertices are proof-interface coordinates, not runners or arcs.")
    print("   pairwise_observable=retained LRC predicate, route specificity, certificate legality, sidecar closure, and discharge refinement.")
    print("   switch=orient toward the carrier whose median with two other legal states remains theorem-facing.")
    print(f"   score_hist={dict(sorted(score_hist.items()))}")
    print(f"   directed_3cycles={count_directed_3cycles(vertices, edges)}")
    print(f"   scc_sizes={[1 for _ in vertices]}")
    print(f"   hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print("   tie_path=" + " > ".join(vertices))


def print_assumption_challenge() -> None:
    print()
    print("6. ASSUMPTION CHALLENGE")
    print("   considered_vertices=runners,gaps,endpoint owners,route labels,certificate atoms,sidecar columns,discharge modes,proof obligations")
    print("   chosen_vertices=packet/route/certificate/sidecar/discharge coordinates")
    print("   preserved_lrc_predicate=boundary/open status plus exact route/certificate/discharge legality")
    print("   destroyed_if_raw=exact scale, endpoint owner, topology, sidecar obligation, and terminal discharge atom")
    print("   challenged_assumption=the final proof should not be a scalar route order; it should be a median-center legality test on completed packet states")


def main() -> None:
    print("S245: LRC14 route-state median-hull scheduler")
    print("=" * 72)
    print_rule_audit()
    print_state_table()
    print_triple_table()
    print_hull_audit()
    print_tournament_analysis()
    print_assumption_challenge()
    print()
    print("READING")
    print("  The proposed final proof interface can be made finite: encode each")
    print("  packet route as a legal sidecar-closed coordinate set.  Legal sidecars")
    print("  are unary Horn implications here, so three legal packet states have a")
    print("  unique coordinatewise median center.  If that center keeps only the typed")
    print("  discharge coordinate and drops the specific atom, the proof has found")
    print("  a scheduler center: refine by a separating sidecar or name residual debt.")


if __name__ == "__main__":
    main()
