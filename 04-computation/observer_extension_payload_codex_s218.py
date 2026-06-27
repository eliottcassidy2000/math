#!/usr/bin/env python3
"""S218: observer-extension/cut payload ledger.

This script synthesizes the first A000568 perspective failure with the
controlled-forgetting carriers from S211/S213/S216/S217.

The point is deliberately small and exact: at the first shifted failure,

    R(5)=48, U(6)=56, defect=8.

The number 12 is not the additive complement of 48 in 56.  Instead it recurs
as a parent/fold/fixed-locus coordinate: U(5)=12, R(4)=12, the 5->6 source and
sink deletion slices each have size 12, and the number of self-converse
6-classes is 12.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, gcd


U = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
}

# R(n) is the number of rooted/node-perspective n-tournaments.
R = {
    1: 1,
    2: 2,
    3: 4,
    4: 12,
    5: 48,
    6: 296,
    7: 3040,
    8: 54256,
}

# Self-converse unlabeled tournaments from the existing S216/S217 small audit.
SELF_CONVERSE = {
    3: 2,
    4: 2,
    5: 8,
    6: 12,
}

VERTEX_ORBIT_PROFILE = {
    3: {1: 1, 3: 1},
    4: {2: 2, 4: 2},
    5: {1: 1, 3: 4, 5: 7},
    6: {2: 5, 4: 10, 6: 41},
}

# S211 Burnside terms for U(6).  The [3,3] term has no fixed vertex and is the
# rootless symmetry warning behind the first shifted failure.
BURNSIDE_U6_TERMS = [
    ("1^6", 1, 32768, 6),
    ("3,1,1,1", 40, 128, 3),
    ("3,3", 40, 32, 0),
    ("5,1", 144, 8, 1),
]

S213 = {
    "rooted5": 48,
    "u6": 56,
    "gap": 8,
    "extension_states": 1408,
    "ordered_pair_perspectives": 1408,
    "directed_edge_perspectives": 704,
    "unordered_pair_perspectives": 704,
    "sector_deck_unique": {
        "size": (55, 56),
        "internal": (55, 56),
        "cross": (56, 56),
        "full": (56, 56),
    },
    "only_size_internal_collision": (344, 345),
}

S216 = {
    "parent_classes_5": 12,
    "raw_extensions_5_to_6": 384,
    "layer_orbits_5_to_6": 296,
    "rooted_6": 296,
    "child_classes_6": 56,
    "source_children": 12,
    "sink_children": 12,
    "parent_aut_hist_5": {1: 7, 3: 4, 5: 1},
    "word_orbit_size_hist_5": {1: 258, 3: 32, 5: 6},
    "rooted_orbits_per_6_class": {2: 5, 4: 10, 6: 41},
    "unique_deletion_parent_hist_6": {1: 3, 2: 6, 3: 13, 4: 13, 5: 17, 6: 4},
    "max_repeated_deletion_parent_hist_6": {1: 4, 2: 27, 3: 17, 4: 3, 5: 2, 6: 3},
}


def frac(num: int, den: int) -> str:
    value = Fraction(num, den)
    return f"{value.numerator}/{value.denominator}"


def weighted_sum(hist: dict[int, int]) -> int:
    return sum(value * count for value, count in hist.items())


def local_flow(k: int) -> dict[str, int]:
    return {
        "k": k,
        "lines": k * (k + 1),
        "rank": 2 * k,
        "rectangle_redundancy": k * (k - 1),
        "line_profile_count": (k + 1) * (k + 2),
        "aligned_profile_count": 2 * comb(k + 3, 3),
    }


def global_flow(n: int) -> dict[str, int]:
    full_tiles = comb(n, 2)
    full_lines = 2 * comb(n, 3)
    full_rank = full_tiles - 1 if n >= 3 else 0
    full_redundancy = full_lines - full_rank
    fixed_tiles = comb(n - 1, 2)
    fixed_lines = 2 * comb(n - 1, 3)
    fixed_rank = fixed_tiles - 1 if n >= 4 else max(0, fixed_tiles - 1)
    fixed_redundancy = fixed_lines - fixed_rank
    return {
        "n": n,
        "full_lines": full_lines,
        "full_rank": full_rank,
        "full_redundancy": full_redundancy,
        "fixed_lines": fixed_lines,
        "fixed_rank": fixed_rank,
        "fixed_redundancy": fixed_redundancy,
        "hourglass_cycles": comb(n - 2, 2),
        "rectangle_cycles": 2 * comb(n - 1, 3),
    }


def directed_three_cycles(adj: list[list[int]]) -> int:
    total = 0
    n = len(adj)
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                total += int(adj[a][b] and adj[b][c] and adj[c][a])
                total += int(adj[a][c] and adj[c][b] and adj[b][a])
    return total


def scc_sizes(adj: list[list[int]]) -> list[int]:
    n = len(adj)

    def reach(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u, bit in enumerate(graph[v]):
                if bit and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    sizes = []
    while remaining:
        start = min(remaining)
        comp = reach(start, adj) & reach(start, rev)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for seen in range(1 << n):
        for v in range(n):
            if not dp[seen][v]:
                continue
            for u in range(n):
                if seen & (1 << u):
                    continue
                if adj[v][u]:
                    dp[seen | (1 << u)][u] += dp[seen][v]
    return sum(dp[-1])


def carrier_tournament() -> tuple[list[tuple[str, dict[str, int]]], list[list[int]]]:
    carriers = [
        ("raw_node_perspective", dict(observer=1, incident=0, edge=0, cross=0, deletion=0, cycle=0, owner=0, proof=0, cost=1)),
        ("incident_word_extension", dict(observer=1, incident=3, edge=0, cross=0, deletion=0, cycle=0, owner=0, proof=0, cost=2)),
        ("ordered_pair_edge_sector", dict(observer=2, incident=3, edge=2, cross=1, deletion=0, cycle=0, owner=0, proof=0, cost=3)),
        ("cross_sector_orientation", dict(observer=2, incident=3, edge=2, cross=3, deletion=0, cycle=0, owner=0, proof=0, cost=4)),
        ("deletion_parent_fiber", dict(observer=1, incident=3, edge=1, cross=2, deletion=3, cycle=0, owner=0, proof=0, cost=4)),
        ("rectangle_hourglass_residue", dict(observer=1, incident=2, edge=1, cross=2, deletion=1, cycle=3, owner=0, proof=0, cost=4)),
        ("endpoint_owner_payload", dict(observer=2, incident=3, edge=2, cross=3, deletion=2, cycle=2, owner=3, proof=0, cost=5)),
        ("proof_obligation_sidecar", dict(observer=2, incident=3, edge=2, cross=3, deletion=3, cycle=3, owner=3, proof=3, cost=6)),
    ]

    def score(data: dict[str, int]) -> tuple[int, int]:
        retained = (
            2 * data["observer"]
            + 3 * data["incident"]
            + 2 * data["edge"]
            + 3 * data["cross"]
            + 3 * data["deletion"]
            + 3 * data["cycle"]
            + 4 * data["owner"]
            + 4 * data["proof"]
        )
        return retained, -data["cost"]

    n = len(carriers)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            si = score(carriers[i][1])
            sj = score(carriers[j][1])
            if si > sj or (si == sj and carriers[i][0] < carriers[j][0]):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return carriers, adj


def one_hamiltonian_path(adj: list[list[int]], names: list[str]) -> list[str]:
    n = len(adj)
    path: list[int] = []
    used = [False] * n

    def dfs(v: int) -> bool:
        path.append(v)
        used[v] = True
        if len(path) == n:
            return True
        for u in range(n):
            if not used[u] and adj[v][u] and dfs(u):
                return True
        used[v] = False
        path.pop()
        return False

    for start in range(n):
        if dfs(start):
            return [names[i] for i in path]
    return []


def print_shift_arithmetic() -> None:
    print("1. SHIFTED PERSPECTIVE ARITHMETIC")
    print("   Correct first-failure arithmetic: 48 + 8 = 56; 48 + 12 = 60.")
    print("   The value 12 is a recurring fixed-locus/parent count, not the additive defect.")
    print()
    print("   m  U(m)  R(m)  U(m+1)  defect  R/U_next  defect/U_next  U(m)/U_next  SC(m+1)/U_next")
    for m in range(3, 7):
        defect = U[m + 1] - R[m]
        sc_next = SELF_CONVERSE.get(m + 1, -1)
        sc_frac = "-" if sc_next < 0 else frac(sc_next, U[m + 1])
        print(
            f"   {m}  {U[m]:4d}  {R[m]:4d}  {U[m+1]:7d}  {defect:6d}  "
            f"{frac(R[m], U[m+1]):>8s}  {frac(defect, U[m+1]):>13s}  "
            f"{frac(U[m], U[m+1]):>11s}  {sc_frac:>15s}"
        )
    print()
    print("   At m=5: R(5)/U(6)=6/7, defect/U(6)=1/7, and U(5)/U(6)=SC(6)/U(6)=3/14.")
    print("   R(5) itself decomposes by vertex-orbit profile as 7*5 + 4*3 + 1*1 = 48.")


def print_twelve_sources() -> None:
    print()
    print("2. WHERE THE 12S COME FROM")
    print("   role                                      value  interpretation")
    rows = [
        ("R(4)=P(4)", R[4], "node perspectives on 4-classes, equal to U(5)"),
        ("U(5)", U[5], "parent isomorphism classes for the 5->6 extension"),
        ("source slice 5->6", S216["source_children"], "all-zero incident word deletion slice"),
        ("sink slice 5->6", S216["sink_children"], "all-one incident word deletion slice"),
        ("SC(6)", SELF_CONVERSE[6], "self-converse/converse-fold branch locus"),
    ]
    for name, value, note in rows:
        print(f"   {name:38s} {value:5d}  {note}")
    print()
    print(f"   U(5) parent automorphism histogram: {S216['parent_aut_hist_5']}")
    print("   This is the same profile as vertex orbit counts reversed: {5:7, 3:4, 1:1}.")


def print_burnside_terms() -> None:
    print()
    print("3. BURNSIDE ROOTLESS WARNING FOR U(6)")
    total = sum(class_size * fixed for _, class_size, fixed, _ in BURNSIDE_U6_TERMS)
    print("   cycle_type  class_size  fixed_tournaments  fixed_vertices  contribution_to_U6")
    for cycle_type, class_size, fixed, fixed_vertices in BURNSIDE_U6_TERMS:
        numerator = class_size * fixed
        print(
            f"   {cycle_type:10s} {class_size:10d} {fixed:18d} "
            f"{fixed_vertices:15d} {str(Fraction(numerator, 720)):>18s}"
        )
    print(f"   Burnside numerator={total}, divided by 6! gives {total // 720}.")
    print("   The [3,3] symmetry has fixed_vertices=0: it is exactly the kind of observerless cut payload a node-root quotient cannot name.")


def print_extension_and_sector_lifts() -> None:
    print()
    print("4. INCIDENT WORDS, EDGE SECTORS, AND DELETION FIBERS")
    raw = S216["raw_extensions_5_to_6"]
    rooted = S216["rooted_6"]
    child = S216["child_classes_6"]
    print(
        f"   5->6 raw incident words: U(5)*2^5 = {raw}; "
        f"parent-aut word orbits = rooted R(6) = {rooted}; child sinks U(6) = {child}."
    )
    print(
        f"   Fractions: raw/rooted={frac(raw, rooted)}, rooted/U(6)={frac(rooted, child)}, raw/U(6)={frac(raw, child)}."
    )
    print(f"   word orbit size histogram under parent automorphisms: {S216['word_orbit_size_hist_5']}")
    print(f"   rooted extension orbits per 6-class: {S216['rooted_orbits_per_6_class']}")
    print(f"   unique deletion-parent classes per 6-class: {S216['unique_deletion_parent_hist_6']}")
    print(f"   max repeated deletion-parent multiplicity per 6-class: {S216['max_repeated_deletion_parent_hist_6']}")
    print()
    print(
        f"   S213 lift: rooted-5 + incident word states = {S213['extension_states']} = "
        f"exact ordered-pair perspectives; directed-edge perspectives = {S213['directed_edge_perspectives']}."
    )
    print("   Sector decks: size/internal separate 55/56 classes; cross/full separate 56/56.")
    print(
        f"   The only size/internal collision is masks {S213['only_size_internal_collision']}; "
        "cross-sector orientation is the chirality repair."
    )


def print_s217_residue_laws() -> None:
    print()
    print("5. S217 RECTANGLE/HOURGLASS RESIDUES")
    print("   k  K_edges  xor_rank  rectangle_red  line_profiles  aligned_profiles")
    for k in range(1, 7):
        row = local_flow(k)
        print(
            f"   {k}  {row['lines']:7d}  {row['rank']:8d}  {row['rectangle_redundancy']:13d} "
            f"{row['line_profile_count']:14d} {row['aligned_profile_count']:17d}"
        )
    row6 = global_flow(6)
    print()
    print(
        "   n=6 full adjacent-layer flow: "
        f"lines={row6['full_lines']}, rank={row6['full_rank']}, redundancy={row6['full_redundancy']} "
        f"= rectangles {row6['rectangle_cycles']} + hourglass {row6['hourglass_cycles']}."
    )
    print(
        "   n=6 fixed-path flow: "
        f"lines={row6['fixed_lines']}, rank={row6['fixed_rank']}, redundancy={row6['fixed_redundancy']}."
    )
    print("   Interpretation: a vanished rectangle/hourglass residue means the quotient descends to potentials; a nonzero residue names hidden payload.")


def print_controlled_forgetting_synthesis() -> None:
    print()
    print("6. ABSTRACT OBSERVER-EXTENSION/CUT PAYLOAD RULE")
    print("   Payload = the minimal coordinate that turns an unsafe observer quotient into a recoverable one.")
    print("   In tournaments it appears as:")
    print("     node root -> incident word -> ordered pair / edge sector -> cross-sector orientation -> deletion fiber -> rectangle/hourglass residue.")
    print("   In the LRC stack the same rule reads:")
    applications = [
        ("pair-good decoys", "blocker generator + barcode bar + active owner support before counting decoys"),
        ("endpoint-owner transfer", "old/new endpoint role is repaired by external owner strips and owner-transfer deltas"),
        ("fiber zipper / Hensel", "automatic/residue fibers need ET bins, unit-root clocks, and magnitude cocycles"),
        ("residual capacitor / Haar square", "rectangle zeta and first nonzero cochain are cut residues"),
        ("status-topology gate", "closed-arc H1 owner support is the observerless cut payload for AP/GW"),
        ("matrix observability", "sidecar columns must separate, reconstruct, annihilate, descend, or name debt"),
    ]
    for name, rule in applications:
        print(f"     {name:30s} -> {rule}")


def print_tournament_analysis() -> None:
    print()
    print("7. TOURNAMENT ANALYSIS OVER PAYLOAD CARRIERS")
    carriers, adj = carrier_tournament()
    names = [name for name, _ in carriers]
    scores = [sum(row) for row in adj]
    print("   vertices are carriers/proof obligations, not runners")
    print(f"   vertices={names}")
    print(f"   score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"   directed_3_cycles={directed_three_cycles(adj)}")
    print(f"   scc_sizes={scc_sizes(adj)}")
    print(f"   hamiltonian_paths={hamiltonian_path_count(adj)}")
    print("   one_hamiltonian_path=" + " -> ".join(one_hamiltonian_path(adj, names)))


def main() -> None:
    print("=" * 80)
    print("S218: observer-extension/cut payload ledger")
    print("=" * 80)
    print_shift_arithmetic()
    print_twelve_sources()
    print_burnside_terms()
    print_extension_and_sector_lifts()
    print_s217_residue_laws()
    print_controlled_forgetting_synthesis()
    print_tournament_analysis()
    print()
    print("READING")
    print("  The first failure is not 48+12=56; it is 48+8=56.")
    print("  The value 12 is a stable fold/parent count: it names U(5), R(4), source/sink deletion slices, and SC(6).")
    print("  The defect 8 is the first observer-extension/cut payload: node perspectives cannot see a rootless coupling, so the proof carrier must add incident words, edge-sector cross orientation, deletion fibers, and finally rectangle/hourglass residues when line-flow payload is attached.")
    print("  Controlled forgetting should therefore be stated as a payload law: a quotient is legal only when the forgotten coordinate is constant, reconstructed, dual-annihilated, descended to potentials, boundary-stopped, or routed to a named residual sidecar.")


if __name__ == "__main__":
    main()
