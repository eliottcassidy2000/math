#!/usr/bin/env python3
"""Exact replay for THM-1156, the radius-1/14 tooth-seam / chi7 law.

All endpoint arithmetic uses ``fractions.Fraction``.  The finite boxes are
not evidence for the theorem: they replay identities proved symbolically in
the accompanying Lean module and expose the induced Paley/Fano fingerprint.

Tournament-analysis declaration
--------------------------------
Vertices are the seven *section residues* in F_7, not runners or teeth.
The pair observable is chi7(y-x); the switch x -> y is chi7(y-x)=+1.
There are no ties, and (0,1,2,3,4,5,6) is the canonical Hamiltonian path.
This quotient preserves the seam's quadratic sign flip and Fano incidence.
It destroys gcd scale, tooth indices, endpoint adjacency, containment, and
the third-support predicate, so it cannot by itself prove a cover theorem.
The arbitrary-label Fano same-colour count is audited separately over all
``2^7`` colourings; it is not inferred from the Paley quotient.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import gcd


LIMIT = 84
INDEX_BOX = 9
QR7 = frozenset({1, 2, 4})


def chi7(x: int) -> int:
    r = x % 7
    if r == 0:
        return 0
    return 1 if r in QR7 else -1


def vp(n: int, p: int) -> int:
    assert n > 0
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e


def xgcd(a: int, b: int) -> tuple[int, int, int]:
    """Return g,x,y with ax+by=g>0 for positive a,b."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t


def left(a: int, m: int) -> Fraction:
    return Fraction(14 * m - 1, 14 * a)


def right(a: int, m: int) -> Fraction:
    return Fraction(14 * m + 1, 14 * a)


def seam_num(a: int, b: int, m: int, n: int) -> int:
    return 14 * (a * n - b * m) - (a + b)


def compatible(a: int, b: int) -> bool:
    return (a + b) % (14 * gcd(a, b)) == 0


def exact_seam(a: int, b: int) -> tuple[int, int]:
    """Construct m,n with R_(a,m)=L_(b,n), assuming compatibility."""
    g = gcd(a, b)
    assert compatible(a, b)
    A, B = a // g, b // g
    gg, x, y = xgcd(A, B)
    assert gg == 1 and A * x + B * y == 1
    c = (a + b) // (14 * g)
    n = x * c
    m = -y * c
    assert seam_num(a, b, m, n) == 0
    return m, n


def determinant_lift(a: int, b: int, z: int) -> tuple[int, int]:
    """Construct m,n with a*n-b*m=gcd(a,b)*z."""
    g = gcd(a, b)
    A, B = a // g, b // g
    gg, x, y = xgcd(A, B)
    assert gg == 1 and A * x + B * y == 1
    n = x * z
    m = -y * z
    assert a * n - b * m == g * z
    return m, n


def r_plus(a: int, b: int) -> int:
    """Least positive numerator on the strict-penetration side."""
    g = gcd(a, b)
    M = 14 * g
    return a + b - M * ((a + b - 1) // M)


def interval_overlap(a: int, m: int, b: int, n: int) -> Fraction:
    lo = max(left(a, m), left(b, n))
    hi = min(right(a, m), right(b, n))
    return max(Fraction(0), hi - lo)


def circle_distance(x: Fraction) -> Fraction:
    r = x % 1
    return min(r, 1 - r)


def paley_edges() -> set[tuple[int, int]]:
    return {(x, y) for x in range(7) for y in range(7) if x != y and chi7(y - x) == 1}


def strongly_connected_components(vertices: tuple[int, ...], edges: set[tuple[int, int]]) -> list[tuple[int, ...]]:
    # Tiny exact Kosaraju replay.
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in vertices:
            if (v, w) in edges and w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    out: list[tuple[int, ...]] = []

    def rdfs(v: int, block: list[int]) -> None:
        seen.add(v)
        block.append(v)
        for w in vertices:
            if (w, v) in edges and w not in seen:
                rdfs(w, block)

    for v in reversed(order):
        if v not in seen:
            block: list[int] = []
            rdfs(v, block)
            out.append(tuple(sorted(block)))
    return sorted(out)


def main() -> None:
    identity_checks = 0
    congruence_checks = 0
    quantum_checks = 0
    directed_quantum_checks = 0
    crossing_checks = 0
    containment_checks = 0
    compatible_pairs: list[tuple[int, int]] = []
    primitive_pairs: list[tuple[int, int]] = []
    seam_support_counts: Counter[int] = Counter()

    for a in range(1, LIMIT + 1):
        for b in range(1, LIMIT + 1):
            g = gcd(a, b)
            M = 14 * g
            s = a + b
            q = (s - 1) // M
            r = r_plus(a, b)
            assert 1 <= r <= M and r == s - M * q
            qm, qn = determinant_lift(a, b, q)
            assert seam_num(a, b, qm, qn) == -r
            assert right(a, qm) - left(b, qn) == Fraction(r, 14 * a * b)
            # q is the largest determinant lattice coordinate with strict
            # penetration.  This checks both minimality and attainability of
            # (6)--(7), for every pair in the box, not only compatible pairs.
            for z in range(q - 3, q + 4):
                N = M * z - s
                if N < 0:
                    assert -N >= r
                elif z > q:
                    assert N >= 0
                directed_quantum_checks += 1

            # Check the exact crossing/containment caveat for every possible
            # nonnegative determinant coordinate on the overlap side.
            for z in range(q + 1):
                m0, n0 = determinant_lift(a, b, z)
                d = g * z
                assert 14 * d < s
                overlap = interval_overlap(a, m0, b, n0)
                penetration = Fraction(s - 14 * d, 14 * a * b)
                if abs(a - b) < 14 * d:
                    assert overlap == penetration
                    crossing_checks += 1
                else:
                    assert 14 * d <= abs(a - b)
                    assert overlap == Fraction(1, 7 * max(a, b))
                    containment_checks += 1

            for m in range(-INDEX_BOX, INDEX_BOX + 1):
                for n in range(-INDEX_BOX, INDEX_BOX + 1):
                    N = seam_num(a, b, m, n)
                    assert left(b, n) - right(a, m) == Fraction(N, 14 * a * b)
                    identity_checks += 1
                    assert (N + a + b) % (14 * g) == 0
                    congruence_checks += 1
                    if compatible(a, b) and N > 0:
                        assert N % (14 * g) == 0
                        assert N >= 14 * g
                        assert left(b, n) - right(a, m) >= Fraction(g, a * b)
                        quantum_checks += 1

            if not compatible(a, b):
                continue
            compatible_pairs.append((a, b))
            m, n = exact_seam(a, b)
            p = right(a, m)
            assert p == left(b, n)
            # The seam is a boundary point of both open danger combs.
            assert circle_distance(a * p) == Fraction(1, 14)
            assert circle_distance(b * p) == Fraction(1, 14)

            A, B = a // g, b // g
            assert gcd(A, B) == 1 and (A + B) % 14 == 0
            assert A % 7 and B % 7 and (A + B) % 7 == 0
            assert chi7(B) == -chi7(A)
            if g == 1:
                primitive_pairs.append((a, b))

            # Literal v7 stripping gives the same sign flip, even when the
            # common gcd has a non-seven part.
            ea, eb = vp(a, 7), vp(b, 7)
            assert ea == eb == vp(g, 7)
            aa, bb = a // (7**ea), b // (7**eb)
            assert aa % 7 and bb % 7 and (aa + bb) % 7 == 0
            assert chi7(bb) == -chi7(aa)

            # Bounded diagnostic only: count third speeds whose open comb
            # actually contains this seam.  The theorem merely says a cover
            # must select one; it does not promise one in this finite bank.
            supporters = [
                c
                for c in range(1, LIMIT + 1)
                if c not in (a, b) and circle_distance(c * p) < Fraction(1, 14)
            ]
            seam_support_counts[len(supporters)] += 1

    # Exact criterion replay: compatibility is equivalent to the constructive
    # Bezout seam throughout the box; incompatibility forbids zero because all
    # seam numerators occupy one nonzero residue class modulo 14g.
    for a in range(1, LIMIT + 1):
        for b in range(1, LIMIT + 1):
            has_constructed = False
            if compatible(a, b):
                m, n = exact_seam(a, b)
                has_constructed = seam_num(a, b, m, n) == 0
            assert has_constructed == compatible(a, b)

    vertices = tuple(range(7))
    edges = paley_edges()
    assert len(edges) == 21
    score_hist = Counter(sum((v, w) in edges for w in vertices if w != v) for v in vertices)
    cycles = {
        frozenset((a, b, c))
        for a in vertices
        for b in vertices
        for c in vertices
        if len({a, b, c}) == 3 and (a, b) in edges and (b, c) in edges and (c, a) in edges
    }
    fano_lines = {frozenset(((t + q) % 7 for q in QR7)) for t in vertices}
    fano_cycles = cycles & fano_lines
    hamiltonian_paths = sum(
        all((path[i], path[i + 1]) in edges for i in range(6))
        for path in permutations(vertices)
    )
    canonical_path = tuple(range(7))
    assert all((canonical_path[i], canonical_path[i + 1]) in edges for i in range(6))
    sccs = strongly_connected_components(vertices, edges)
    negation_flips = sum(
        ((-x) % 7, (-y) % 7) not in edges
        for x, y in edges
    )
    assert score_hist == Counter({3: 7})
    assert len(cycles) == 14 and len(fano_lines) == 7 and len(fano_cycles) == 7
    assert sccs == [vertices]
    assert negation_flips == 21

    # Independent arbitrary-label Fano audit.  This is the elementary
    # two-colour/pigeonhole law used by THM-1156, not a Paley consequence.
    fano_same_pair_hist: Counter[int] = Counter()
    fano_blocked_line_hist: Counter[int] = Counter()
    for mask in range(1 << 7):
        colour = tuple((mask >> v) & 1 for v in vertices)
        p = sum(colour)
        same_pairs = p * (p - 1) // 2 + (7 - p) * (6 - p) // 2
        blocked_lines = sum(
            any(colour[x] == colour[y] for x in line for y in line if x < y)
            for line in fano_lines
        )
        assert same_pairs >= 9
        assert blocked_lines == 7
        fano_same_pair_hist[same_pairs] += 1
        fano_blocked_line_hist[blocked_lines] += 1

    seam_pairs_mod7 = {tuple(sorted((r, (-r) % 7))) for r in range(1, 7)}
    assert seam_pairs_mod7 == {(1, 6), (2, 5), (3, 4)}
    assert all(chi7(x) == -chi7(y) for x, y in seam_pairs_mod7)

    # The actual compatibility graph on the finite runner bank is bipartite;
    # unlike the Paley sidecar it is not a tournament.
    speed_colour = {v: chi7(v // (7 ** vp(v, 7))) for v in range(1, LIMIT + 1)}
    seam_graph_edges = {
        (a, b)
        for a in range(1, LIMIT + 1)
        for b in range(a + 1, LIMIT + 1)
        if compatible(a, b)
    }
    assert all(speed_colour[a] == -speed_colour[b] for a, b in seam_graph_edges)
    seam_triangles = sum(
        (a, b) in seam_graph_edges and (a, c) in seam_graph_edges and (b, c) in seam_graph_edges
        for a in range(1, LIMIT + 1)
        for b in range(a + 1, LIMIT + 1)
        for c in range(b + 1, LIMIT + 1)
    )
    assert seam_triangles == 0

    print("THM-1156 exact tooth-seam / chi7 replay")
    print(f"box: 1 <= a,b <= {LIMIT}; tooth indices |m|,|n| <= {INDEX_BOX}")
    print(f"endpoint identities checked: {identity_checks}")
    print(f"gcd congruences checked: {congruence_checks}")
    print(f"positive compatible seam quanta checked: {quantum_checks}")
    print(f"arbitrary-pair directed r_plus checks: {directed_quantum_checks}")
    print(f"crossing-regime identities checked: {crossing_checks}")
    print(f"containment-regime identities checked: {containment_checks}")
    print(f"ordered compatible pairs: {len(compatible_pairs)}")
    print(f"ordered primitive compatible pairs: {len(primitive_pairs)}")
    print("criterion: exists exact seam iff 14*gcd(a,b) divides a+b [constructive Bezout replay PASS]")
    print("r_plus law: least strict penetration is attained and minimal for every pair [PASS]")
    print("caveat: penetration equals intersection only in crossing; containment plateau checked separately [PASS]")
    print("v7 law: compatible pairs have equal v7 and opposite chi7 after stripping [PASS]")
    print("open seam law: both owner combs attain distance exactly 1/14, so neither covers the point [PASS]")
    print(f"bounded third-support count histogram (diagnostic only): {sorted(seam_support_counts.items())}")
    print()
    print("Exact seam graph / arbitrary-label Fano audit")
    print(f"finite speed-bank seam edges: {len(seam_graph_edges)}; chi7-crossing edges: {len(seam_graph_edges)}")
    print(f"finite speed-bank seam triangles: {seam_triangles}")
    print(f"all 128 Fano colourings same-colour-pair histogram: {sorted(fano_same_pair_hist.items())}")
    print(f"all 128 Fano colourings blocked-line histogram: {sorted(fano_blocked_line_hist.items())}")
    print("proof role: arbitrary-label pigeonhole audit; independent of the Paley tournament quotient")
    print()
    print("Tournament/Fano audit")
    print("role: quotient perspective only; not the proof of the arbitrary-label Fano count")
    print("vertices: seven section residues F_7 (challenged alternatives: runners, teeth, seam events)")
    print("pair observable: chi7(y-x); switch: x -> y iff chi7(y-x)=+1")
    print(f"tie Hamiltonian path: {canonical_path} (no ties occur)")
    print(f"score histogram: {sorted(score_hist.items())}")
    print(f"directed 3-cycles: {len(cycles)}; Fano-line cycles: {len(fano_cycles)}/{len(fano_lines)}")
    print(f"SCCs: {sccs}")
    print(f"Hamiltonian-path count: {hamiltonian_paths}")
    print(f"edge flips under seam involution x -> -x: {negation_flips}/21")
    print(f"nonzero seam pairs mod 7: {sorted(seam_pairs_mod7)}")
    print("preserves: quadratic seam sign and Fano difference incidence")
    print("destroys: gcd quantum, tooth indices, cyclic adjacency, containment, and third support")
    print("scope: local arithmetic/owner law only; no seven-comb overlap constant and no LRC(14) closure")


if __name__ == "__main__":
    main()
