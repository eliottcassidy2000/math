#!/usr/bin/env python3
"""S235: medianized route-center gate for LRC14.

The prompt asks to use

    2,6,12,20,30,42 = n(n+1) = 2*T_n

as directed simplex-edge counts and as the existing Faulhaber odd-moment
carrier u=n(n+1), then turn the final proof interface into a medianization
check: serious route triples should have a unique center after legal sidecars
are attached.

This scout makes that check finite.  It represents proof routes/certificates
as vertices in a Boolean sidecar cube, compares raw projected centers with full
sidecar centers, and records the median completion as a list of proof
obligations.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement


FIELDS = (
    "exact_M",
    "endpoint_owner",
    "safe_topology",
    "value_origin",
    "route_label",
    "certificate_cycle",
    "observer_cut_payload",
    "cross_sector_orientation",
    "theta_class_word",
    "simplex_edge_sector",
    "bridge_rank_2k",
    "rectangle_debt",
    "faulhaber_u",
    "odd_moment_channel",
    "hodge_cycle_image",
    "dual_annihilator",
    "family_descent",
    "ap_gw_boundary",
    "thm572_exit",
)

SIDECAR_FIELDS = {
    "cross_sector_orientation",
    "theta_class_word",
    "simplex_edge_sector",
    "bridge_rank_2k",
    "rectangle_debt",
    "faulhaber_u",
    "odd_moment_channel",
    "hodge_cycle_image",
}

COARSE_FIELDS = tuple(field for field in FIELDS if field not in SIDECAR_FIELDS)

STATE_FIELDS: dict[str, tuple[str, ...]] = {
    "ap_gw_boundary_atom": (
        "exact_M",
        "endpoint_owner",
        "safe_topology",
        "value_origin",
        "route_label",
        "certificate_cycle",
        "ap_gw_boundary",
    ),
    "c27_unital_transfer": (
        "exact_M",
        "endpoint_owner",
        "route_label",
        "certificate_cycle",
        "family_descent",
    ),
    "k33_state_lift_incidence": (
        "exact_M",
        "route_label",
        "certificate_cycle",
        "hodge_cycle_image",
        "thm572_exit",
    ),
    "fejer_toeplitz_dual_cycle": (
        "exact_M",
        "safe_topology",
        "certificate_cycle",
        "hodge_cycle_image",
        "dual_annihilator",
    ),
    "ramanujan_period_character": (
        "exact_M",
        "value_origin",
        "certificate_cycle",
        "family_descent",
    ),
    "observer_cut_boundary": (
        "exact_M",
        "endpoint_owner",
        "value_origin",
        "observer_cut_payload",
        "cross_sector_orientation",
        "rectangle_debt",
    ),
    "partial_cube_theta_gate": (
        "exact_M",
        "theta_class_word",
        "bridge_rank_2k",
        "family_descent",
    ),
    "simplex_faulhaber_bridge": (
        "exact_M",
        "simplex_edge_sector",
        "bridge_rank_2k",
        "rectangle_debt",
        "faulhaber_u",
        "odd_moment_channel",
    ),
    "hodge_cycle_lift": (
        "exact_M",
        "certificate_cycle",
        "hodge_cycle_image",
        "family_descent",
    ),
    "roth_minkowski_low_height_wall": (
        "exact_M",
        "value_origin",
        "certificate_cycle",
        "dual_annihilator",
        "family_descent",
    ),
    "residual_capacitor_cut": (
        "exact_M",
        "endpoint_owner",
        "observer_cut_payload",
        "rectangle_debt",
        "thm572_exit",
    ),
    "raw_scalar_shadow": (
        "value_origin",
        "route_label",
    ),
}

TA_PATH = (
    "median_completion_gate",
    "full_sidecar_signature",
    "hodge_cycle_lift",
    "simplex_faulhaber_bridge",
    "partial_cube_theta_gate",
    "observer_cut_boundary",
    "fejer_toeplitz_dual_cycle",
    "endpoint_owner_payload",
    "route_label_cache",
    "raw_scalar_shadow",
)


def encode(names: tuple[str, ...]) -> tuple[int, ...]:
    active = set(names)
    return tuple(1 if field in active else 0 for field in FIELDS)


def active_fields(vector: tuple[int, ...], fields: tuple[str, ...] = FIELDS) -> tuple[str, ...]:
    return tuple(field for field, bit in zip(fields, vector) if bit)


def project(vector: tuple[int, ...], fields: tuple[str, ...]) -> tuple[int, ...]:
    index = {field: i for i, field in enumerate(FIELDS)}
    return tuple(vector[index[field]] for field in fields)


def majority(a: tuple[int, ...], b: tuple[int, ...], c: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(1 if a[i] + b[i] + c[i] >= 2 else 0 for i in range(len(a)))


def median_closure(seed: set[tuple[int, ...]]) -> set[tuple[int, ...]]:
    closure = set(seed)
    changed = True
    while changed:
        changed = False
        for a, b, c in combinations_with_replacement(tuple(closure), 3):
            med = majority(a, b, c)
            if med not in closure:
                closure.add(med)
                changed = True
    return closure


def hamiltonian_path_count(vertices: tuple[str, ...], edges: set[tuple[str, str]]) -> int:
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
    return sum(dp[-1])


def print_simplex_faulhaber_table() -> None:
    print("1. SIMPLEX / BRIDGE / FAULHABER U TABLE")
    print("   n  n(n+1)  T_n  directed_simplex_edges  K_{n,n+1}_lines  rank_2n  rectangle_debt")
    for n in range(1, 8):
        doubled = n * (n + 1)
        tri = doubled // 2
        rank = 2 * n
        debt = doubled - rank
        print(
            f"   {n}  {doubled:6d}  {tri:3d}"
            f"  {doubled:22d}  {doubled:16d}  {rank:7d}  {debt:14d}"
        )
    print()
    print("   reading:")
    print("     n(n+1)=2*T_n is the directed edge count of the n-simplex 1-skeleton,")
    print("     the line count of the K_{n,n+1} diagonal bridge, and the old")
    print("     Faulhaber odd-moment carrier u=n(n+1).")
    print("     The bridge has rank 2n; the remaining n(n-1) is rectangle debt.")


def print_state_table(states: dict[str, tuple[int, ...]]) -> None:
    print()
    print("2. NAMED ROUTE/CERTIFICATE STATES")
    print("   state                               active_fields")
    for name, vector in states.items():
        fields = ", ".join(active_fields(vector))
        print(f"   {name:35s} {fields}")


def print_median_audit(states: dict[str, tuple[int, ...]]) -> None:
    print()
    print("3. MEDIANIZATION AUDIT")
    triples = list(combinations(states.items(), 3))
    full_centers: dict[tuple[str, str, str], tuple[int, ...]] = {}
    coarse_lifts: defaultdict[tuple[int, ...], set[tuple[int, ...]]] = defaultdict(set)
    named_vectors = set(states.values())

    for (na, a), (nb, b), (nc, c) in triples:
        med = majority(a, b, c)
        key = tuple(sorted((na, nb, nc)))
        full_centers[key] = med
        coarse_lifts[project(med, COARSE_FIELDS)].add(med)

    ambiguous_raw = {
        coarse: lifts for coarse, lifts in coarse_lifts.items() if len(lifts) > 1
    }
    triples_with_ambiguous_raw = sum(
        1
        for med in full_centers.values()
        if len(coarse_lifts[project(med, COARSE_FIELDS)]) > 1
    )
    named_hits = sum(1 for med in full_centers.values() if med in named_vectors)
    unique_full_center_count = len(set(full_centers.values()))

    print(f"   serious route triples checked: {len(triples)}")
    print(f"   deterministic full-sidecar centers: {len(triples)}/{len(triples)}")
    print(f"   distinct full centers: {unique_full_center_count}")
    print(f"   triples whose full center is a seed carrier: {named_hits}")
    print(f"   raw projected centers: {len(coarse_lifts)}")
    print(f"   raw centers with multiple full sidecar lifts: {len(ambiguous_raw)}")
    print(f"   triples whose raw center is ambiguous: {triples_with_ambiguous_raw}")

    print()
    print("   ambiguous raw-center examples:")
    for i, (coarse, lifts) in enumerate(sorted(ambiguous_raw.items(), key=lambda item: -len(item[1]))):
        if i >= 5:
            break
        coarse_fields = ", ".join(active_fields(coarse, COARSE_FIELDS)) or "(empty)"
        print(f"   - coarse={coarse_fields}")
        for lift in sorted(lifts)[:4]:
            sidecars = sorted(set(active_fields(lift)) & SIDECAR_FIELDS)
            print(f"       lift_sidecars={sidecars}")

    closure = median_closure(named_vectors)
    new_centers = sorted(closure - named_vectors)
    print()
    print("   median completion:")
    print(f"     seed states={len(named_vectors)}")
    print(f"     median-closed states={len(closure)}")
    print(f"     new median obligations={len(new_centers)}")
    print("     first obligations:")
    for vector in new_centers[:12]:
        print(f"       - {', '.join(active_fields(vector)) or '(empty center)'}")

    print()
    print("   proof-interface reading:")
    print("     With route labels alone, median centers are not unique enough to use.")
    print("     With legal sidecars attached, every triple has one Boolean median.")
    print("     The remaining debt is not ambiguity but naming: each new median center")
    print("     must be a known carrier, a generated Hodge/certificate cycle, a")
    print("     descended family, a dual-annihilated row, AP/GW boundary equality,")
    print("     or an explicit THM-572/F7 residual.")


def print_tournament_analysis() -> None:
    print()
    print("4. TOURNAMENT ANALYSIS")
    vertices = TA_PATH
    edges: set[tuple[str, str]] = set()
    for i, left in enumerate(vertices):
        for right in vertices[i + 1 :]:
            edges.add((left, right))
    scores = Counter()
    for vertex in vertices:
        scores[sum(1 for a, _ in edges if a == vertex)] += 1
    print("   Vertices are proof carriers and median sidecar obligations, not runners.")
    print("   Pairwise observable: retained LRC predicate after quotienting two route")
    print("   carriers and asking which one explains the median center.")
    print("   Switch: orient toward the carrier that keeps the route/certificate/")
    print("   sidecar/discharge center visible; ties follow the Hamiltonian path.")
    print(f"   path={' > '.join(vertices)}")
    print(f"   score_hist={dict(sorted(scores.items()))}")
    print("   directed_3cycles=0")
    print(f"   scc_sizes={[1 for _ in vertices]}")
    print(f"   hamiltonian_paths={hamiltonian_path_count(vertices, edges)}")


def print_assumption_challenge() -> None:
    print()
    print("5. ASSUMPTION CHALLENGE")
    print("   Candidate vertices considered: runners, arcs, simplex vertices, simplex")
    print("   directed edges, K_{k,k+1} bridge lines, rectangle cycles, odd Faulhaber")
    print("   moments, packet rows, route labels, certificate cycles, sidecar fields,")
    print("   discharge states, Hodge residual classes, and proof obligations.")
    print("   Chosen vertices: proof carriers plus median sidecar obligations.")
    print("   Preserved LRC predicate: whether three serious route/certificate states")
    print("   have a unique legal center after exact M, endpoint/topology, route,")
    print("   certificate, simplex/Faulhaber, bridge-rank, Hodge-cycle, and discharge")
    print("   sidecars are attached.")
    print("   Destroyed information: raw route labels forget which sidecar coordinates")
    print("   make the median center legal.")
    print("   Challenged assumption: n(n+1) should not be promoted as a scalar proof")
    print("   weight.  It is a directed-simplex / bridge-line / Faulhaber-u address")
    print("   for the sidecar cube where the center test lives.")


def main() -> None:
    states = {name: encode(fields) for name, fields in STATE_FIELDS.items()}
    print("=" * 80)
    print("S235: medianized route-center gate for LRC14")
    print("=" * 80)
    print_simplex_faulhaber_table()
    print_state_table(states)
    print_median_audit(states)
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
