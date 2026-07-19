#!/usr/bin/env python3
"""Exact replay for THM-1161: the local three-tooth factor is sharply one.

All arithmetic is ``fractions.Fraction`` arithmetic.  The program checks the
core-safe certificate, the legal infinite family, the complete local wall
chronology, both equal-split normalizations, and a small tournament audit.
It is a proof replay of displayed identities, not a search for examples.
"""

from fractions import Fraction as Q
from itertools import permutations


P = (1, 2, 3, 5, 6, 7, 9, 10)
CORE = (Q(29, 126), Q(27, 98))


def tooth(k: int, q: int) -> tuple[Q, Q]:
    """Closed threshold-1/14 tooth of comb k centered at q/k."""
    return (Q(q, k) - Q(1, 14 * k), Q(q, k) + Q(1, 14 * k))


def intersects(x: tuple[Q, Q], y: tuple[Q, Q]) -> bool:
    return max(x[0], y[0]) <= min(x[1], y[1])


def core_certificate() -> tuple[tuple[int, int, Q, Q], ...]:
    """Certify p*CORE lies in one safe unit sector [n+1/14,n+13/14]."""
    rows = []
    for p in P:
        left, right = p * CORE[0], p * CORE[1]
        candidates = [
            n
            for n in range(0, 11)
            if Q(n) + Q(1, 14) <= left
            and right <= Q(n) + Q(13, 14)
        ]
        assert len(candidates) == 1
        rows.append((p, candidates[0], left, right))
    return tuple(rows)


def family(K: int) -> dict[str, object]:
    """Return and exactly verify the K-th member of the counterfamily."""
    assert K % 4 == 0 and K >= 132
    assert K > 13 * max(P)  # the legal-killer threshold used by THM-1140
    assert all(K + d > 13 * max(P) for d in range(4))

    m = K // 4
    A = Q(1, 4) + Q(1, 14 * K)
    B = Q(1, 4) + Q(13, 14 * K)
    gap = (A, B)
    assert B - A == Q(6, 7 * K)
    assert CORE[0] < A < B < CORE[1]

    # The sole tooth of comb K+d meeting the gap has index m+1.  Check this
    # directly against every potentially nearby integer, independently of the
    # closed formulas used below.
    unique_indices = {}
    teeth = {}
    for d in (1, 2, 3):
        k = K + d
        lower = k * A - Q(1, 14)
        upper = k * B + Q(1, 14)
        assert Q(m) < lower < Q(m + 1) < upper < Q(m + 2)
        candidates = tuple(
            q
            for q in range(lower.numerator // lower.denominator - 1,
                           upper.numerator // upper.denominator + 2)
            if intersects(gap, tooth(k, q))
        )
        assert candidates == (m + 1,)
        unique_indices[d] = candidates[0]
        teeth[d] = tooth(k, m + 1)

        expected_left = Q(1, 4) + Q(26 - 7 * d, 28 * (K + d))
        expected_right = Q(1, 4) + Q(30 - 7 * d, 28 * (K + d))
        assert teeth[d] == (expected_left, expected_right)

    # The complete local event carrier.  It alternates surviving pieces and
    # teeth; hence it also proves containment, disjointness, and chronology.
    walls = (
        A,
        teeth[3][0],
        teeth[3][1],
        teeth[2][0],
        teeth[2][1],
        teeth[1][0],
        teeth[1][1],
        B,
    )
    assert all(x < y for x, y in zip(walls, walls[1:]))

    pieces = (
        walls[1] - walls[0],
        walls[3] - walls[2],
        walls[5] - walls[4],
        walls[7] - walls[6],
    )
    formulas = (
        Q(3 * (K - 2), 28 * K * (K + 3)),
        Q(3 * (K + 6), 28 * (K + 2) * (K + 3)),
        Q(3 * K + 22, 28 * (K + 1) * (K + 2)),
        Q(3 * K + 26, 28 * K * (K + 1)),
    )
    assert pieces == formulas
    assert Q(0) < pieces[0] < pieces[1] < pieces[2] < pieces[3]

    exact_mean = sum(pieces, Q(0)) / 4
    exact_mean_formula = Q(
        3 * K**3 + 24 * K**2 + 55 * K + 36,
        28 * K * (K + 1) * (K + 2) * (K + 3),
    )
    crude = (Q(6, 7 * K) - Q(3, 7 * (K + 1))) / 4
    assert exact_mean == exact_mean_formula
    assert crude == Q(3 * (K + 2), 28 * K * (K + 1))
    assert crude <= exact_mean <= pieces[3]

    crude_ratio = pieces[3] / crude
    mean_ratio = pieces[3] / exact_mean
    assert crude_ratio == Q(3 * K + 26, 3 * (K + 2))
    assert crude_ratio - 1 == Q(20, 3 * K + 6)
    assert mean_ratio - 1 == Q(
        17 * K**2 + 93 * K + 120,
        3 * K**3 + 24 * K**2 + 55 * K + 36,
    )

    return {
        "K": K,
        "gap": gap,
        "indices": unique_indices,
        "teeth": teeth,
        "walls": walls,
        "pieces": pieces,
        "crude": crude,
        "exact_mean": exact_mean,
        "crude_ratio": crude_ratio,
        "mean_ratio": mean_ratio,
        "threshold_ratio": pieces[3] * 7 * (K + 3),
    }


def tournament_fingerprint(pieces: tuple[Q, ...]) -> dict[str, object]:
    """Tournament on surviving pieces, oriented by increasing exact length."""
    n = len(pieces)
    edge = {(i, j): pieces[i] < pieces[j] for i in range(n) for j in range(n) if i != j}
    scores = tuple(sum(edge[i, j] for j in range(n) if i != j) for i in range(n))
    triangles = sum(
        (edge[i, j] and edge[j, k] and edge[k, i])
        or (edge[j, i] and edge[i, k] and edge[k, j])
        for i in range(n)
        for j in range(i + 1, n)
        for k in range(j + 1, n)
    )
    paths = tuple(
        perm
        for perm in permutations(range(n))
        if all(edge[perm[i], perm[i + 1]] for i in range(n - 1))
    )
    reach = [[i == j or (i != j and edge[i, j]) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc_sizes_list = []
    while unseen:
        i = min(unseen)
        component = {j for j in unseen if reach[i][j] and reach[j][i]}
        scc_sizes_list.append(len(component))
        unseen -= component
    scc_sizes = tuple(sorted(scc_sizes_list, reverse=True))
    assert scores == (3, 2, 1, 0)
    assert triangles == 0
    assert paths == ((0, 1, 2, 3),)
    return {
        "scores": scores,
        "score_histogram": tuple(sorted(scores)),
        "directed_triangles": triangles,
        "scc_sizes": scc_sizes,
        "hamiltonian_paths": paths,
    }


def main() -> None:
    print("THM-1161 exact replay: the local three-tooth factor is sharply one")
    print(f"core P={P}, max(P)={max(P)}, legal threshold k>={13 * max(P) + 1}")
    print(f"core-safe component J=[{CORE[0]},{CORE[1]}]")
    print("safe-sector certificate p : n ; pJ subset [n+1/14,n+13/14]")
    for p, n, left, right in core_certificate():
        print(f"  {p:2d} : {n} ; [{left},{right}]")

    rows = [family(K) for K in (132, 136, 400, 1000, 10000)]
    print("\nexact family (K=0 mod 4, K>=132; killers K,K+1,K+2,K+3)")
    for row in rows:
        K = row["K"]
        print(f"K={K}")
        print(f"  K-gap=[{row['gap'][0]},{row['gap'][1]}]")
        print(
            "  unique foreign tooth indices (d:index)="
            + ",".join(f"{d}:{q}" for d, q in row["indices"].items())
        )
        print("  pieces=" + ",".join(str(x) for x in row["pieces"]))
        print(
            f"  crude={row['crude']} exact_mean={row['exact_mean']} "
            f"max/crude={row['crude_ratio']} ({float(row['crude_ratio']):.12f}) "
            f"max/mean={row['mean_ratio']} ({float(row['mean_ratio']):.12f})"
        )
        print(
            f"  7*(K+3)*max={row['threshold_ratio']} "
            f"({float(row['threshold_ratio']):.12f})"
        )

    first = rows[0]
    assert first["crude_ratio"] < Q(1295, 1000) < Q(4, 3)
    assert first["mean_ratio"] < Q(1295, 1000)
    print("\nsharpness identities")
    print("  max/crude = 1 + 20/(3K+6) -> 1")
    print(
        "  max/exact_mean = 1 + (17K^2+93K+120)/"
        "(3K^3+24K^2+55K+36) -> 1"
    )
    print(
        "  K=132 already gives max/crude="
        f"{first['crude_ratio']} < 1.295 < 4/3"
    )
    print("  pigeonhole gives factor >=1; the displayed family makes 1 sharp")

    # A larger exact sweep is only a replay guard: the proof is the formulas.
    sweep = [family(K) for K in range(132, 4097, 4)]
    assert all(row["crude_ratio"] < Q(1295, 1000) for row in sweep)
    assert all(
        sweep[i + 1]["crude_ratio"] < sweep[i]["crude_ratio"]
        for i in range(len(sweep) - 1)
    )
    first_edges = tuple(
        first["pieces"][i] < first["pieces"][j]
        for i in range(4) for j in range(i + 1, 4)
    )
    assert all(
        tuple(row["pieces"][i] < row["pieces"][j]
              for i in range(4) for j in range(i + 1, 4)) == first_edges
        for row in sweep
    )
    print(f"\nexact sweep: {len(sweep)} legal family members; all formula checks passed")

    fp = tournament_fingerprint(first["pieces"])
    print("\ntournament / alternate-carrier audit")
    print("  vertices: four surviving pieces p0,p1,p2,p3")
    print("  observable: Delta(i,j)=length(pj)-length(pi)")
    print("  switch/gauge: orient i->j iff Delta(i,j)>0; sign reversal reverses all edges")
    print("  tie Hamiltonian path: p0,p1,p2,p3 (there are no ties in the family)")
    print(
        f"  scores={fp['scores']} histogram={fp['score_histogram']} "
        f"triangles={fp['directed_triangles']} SCCs={fp['scc_sizes']} "
        f"Hamiltonian_paths={len(fp['hamiltonian_paths'])} edge_flips_across_family=0"
    )
    print("  piece quotient preserves the local max/mean factor but loses tooth ownership")
    print("  tooth carrier preserves comb chronology but loses the two boundary pieces")
    print("  eight wall events A,L3,R3,L2,R2,L1,R1,B preserve the full local predicate")
    print("  runner vertices alone lose the chosen gap phase and cannot decide the factor")
    print("DONE")


if __name__ == "__main__":
    main()
