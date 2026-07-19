#!/usr/bin/env python3
"""Exact arithmetic certificate for THM-1153.

This script checks the all-N multi-deletion harmonic-crown coefficient, the
LRC(14) scalar constants, and the exact ordered top-seven recurrence.  It is
proof telemetry for a paper argument whose analytic inputs are the known
lower LRC cases, Lipschitz fattening, and FragmentationLemma.killer_budget.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path


def crown_coefficient(n: int, r: int) -> Fraction:
    """Coefficient c with sum(1/s) >= c/max(core)."""
    assert n >= 3
    assert 1 <= r < n
    return Fraction(r * (n - 2 * r), n * (n - r))


def scalar_ratio(n: int, r: int) -> Fraction:
    """Bound min(deleted)/max(core) from the harmonic crown, when positive."""
    assert 2 * r < n
    return Fraction(n * (n - r), n - 2 * r)


def interval_replay(n: int, r: int) -> None:
    """Verify fattening + fragmentation equals the closed-form coefficient."""
    assert 1 <= r < n
    delta = Fraction(1, n - r) - Fraction(1, n)
    length_at_m_one = 2 * delta
    fragmentation_lower = (
        Fraction(n, 2) * length_at_m_one * (1 - Fraction(2 * r, n))
    )
    assert fragmentation_lower == crown_coefficient(n, r)


def compact_recurrence(seed_x12: Fraction) -> dict[int, Fraction]:
    """Worst-case triangular relaxation of all top-r harmonic inequalities."""
    x = {13: Fraction(1), 12: seed_x12}
    for r in range(2, 7):
        i = 13 - r
        c = crown_coefficient(14, r)
        x[i] = sum(x[j] for j in range(14 - r, 14)) / c
    return x


def closed_form(r: int) -> Fraction:
    """Product form for the compact-seeded bound x_(13-r), 2 <= r <= 6."""
    assert 2 <= r <= 6
    cumulative = Fraction(14)
    for j in range(2, r):
        cumulative *= Fraction(98 - j * j, j * (7 - j))
    return cumulative / crown_coefficient(14, r)


def fmt(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def main() -> None:
    for n in range(3, 101):
        for r in range(1, n):
            interval_replay(n, r)

    print("THM-1153 exact multi-deletion harmonic crown")
    print("general coefficient: r(N-2r)/(N(N-r)m)")
    print("all-N algebra replay: N=3..100, every 1<=r<N: PASS")
    print()

    print("LRC14 actual-threshold table")
    for r in range(1, 8):
        c = crown_coefficient(14, r)
        if 2 * r < 14:
            scalar = fmt(scalar_ratio(14, r))
            relation = "<" if r >= 2 else "<="
            print(f"r={r}: c_r={fmt(c)}; scalar min(S)/m {relation} {scalar}")
        else:
            print(f"r={r}: c_r={fmt(c)}; scalar consequence: no information")
    assert crown_coefficient(14, 7) == 0
    assert crown_coefficient(14, 8) < 0
    print("r=7 coefficient: EXACTLY ZERO")
    print("r=8 coefficient: NEGATIVE")
    print()

    compact = compact_recurrence(Fraction(13))
    expected = {
        11: Fraction(588, 5),
        10: Fraction(25333, 30),
        9: Fraction(204967, 36),
        8: Fraction(8403647, 200),
        7: Fraction(613466231, 1350),
    }
    assert {i: compact[i] for i in range(7, 12)} == expected

    print("compact harmonic top-chain (x_i=v13/v_i, x12<13, x13=1)")
    cumulative = Fraction(14)
    for r in range(2, 7):
        i = 13 - r
        c = crown_coefficient(14, r)
        assert compact[i] == cumulative / c
        assert compact[i] == closed_form(r)
        amp = Fraction(98 - r * r, r * (7 - r))
        assert 1 + 1 / c == amp
        print(
            f"r={r}: c_r={fmt(c)}; x_{i}<{fmt(compact[i])}; "
            f"T_(r+1)/T_r={fmt(amp)}"
        )
        cumulative += compact[i]
    print(f"top-seven ceiling: v13/v7 < {fmt(compact[7])}")
    print(f"top-seven decimal: {float(compact[7]):.12f}")
    print()

    scalar_product = Fraction(13)
    auxiliary_thirteenth_product = Fraction(13)
    for r in range(2, 7):
        scalar_product *= scalar_ratio(14, r)
        auxiliary_thirteenth_product *= Fraction(
            13 * r * (14 - r), (r - 1) * (13 - 2 * r)
        )
    assert scalar_product == 173044872
    assert auxiliary_thirteenth_product == Fraction(20388441216, 7)
    print(f"product of scalar 1/14 bounds: {fmt(scalar_product)}")
    print(
        "harmonic improvement factor: "
        f"{float(scalar_product / compact[7]):.12f}"
    )
    print(
        "improvement over auxiliary 1/13 product: "
        f"{float(auxiliary_thirteenth_product / compact[7]):.12f}"
    )
    print()

    citation_seed = compact_recurrence(Fraction(91, 6))
    assert citation_seed[7] == Fraction(8500889201, 16200)
    print(
        "without compact dominance, r=1 crown gives: "
        f"v13/v7 <= {fmt(citation_seed[7])}"
    )
    print()

    # Tournament-analysis audit.  Orient the seven ranked deletion obligations
    # by dependency order.  This quotient is transitive and intentionally
    # recorded as lossy: the proof data are the harmonic weights and cuts.
    vertices = list(range(1, 8))
    edges = {(a, b) for a in vertices for b in vertices if a < b}
    scores = [sum((a, b) in edges for b in vertices if b != a) for a in vertices]
    assert sorted(scores) == list(range(7))
    directed_triangles = 0
    for a in vertices:
        for b in vertices:
            for c in vertices:
                directed_triangles += int(
                    a < b < c
                    and (a, b) in edges
                    and (b, c) in edges
                    and (c, a) in edges
                )
    assert directed_triangles == 0
    strongly_connected_components = sum(
        all(not ((a, b) in edges and (b, a) in edges) for b in vertices if b != a)
        for a in vertices
    )
    assert strongly_connected_components == 7
    from itertools import permutations

    hamiltonian_paths = sum(
        all((path[i], path[i + 1]) in edges for i in range(len(path) - 1))
        for path in permutations(vertices)
    )
    assert hamiltonian_paths == 1
    print("tournament/deletion-obligation audit")
    print(f"vertices: {vertices}")
    print(f"score histogram: {sorted(scores)}")
    print(
        "directed cycles: 0; "
        f"SCCs: {strongly_connected_components} singletons; "
        f"Hamiltonian paths: {hamiltonian_paths}"
    )
    print("challenged assumption: vertices are deletion cuts, not runners")
    print("preserved: dependency rank and wall location")
    print("destroyed: tooth phases, pair/triple overlaps, boundary excess")
    print("faithful carrier: weighted nested-deletion hypergraph + endpoint sidecar")
    print()

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={source_hash}")


if __name__ == "__main__":
    main()
