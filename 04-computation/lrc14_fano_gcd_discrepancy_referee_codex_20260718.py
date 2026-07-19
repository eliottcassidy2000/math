#!/usr/bin/env python3
"""Independent exact referee for THM-1166's seven-wall Fano/gcd debt.

The paper proof has four inputs, all replayed here with exact rationals:

1. the sharp one-period interval discrepancy ``rho(1-rho)/g``;
2. the actual-radius pair-density range ``1/91 <= rho(a,b) <= 1/14``;
3. the Fano decomposition of the optimal quadratic cover majorant; and
4. the resulting linewise gcd debt and common-dilate constant ``3523/36``.

Tournament-analysis declaration
--------------------------------
The faithful vertices are the seven Fano *lines* (three-speed proof
obligations), not runners.  Comparing the seven scalar line debts and
breaking ties by the displayed line order gives a transitive tournament,
whose tie Hamiltonian path is that order.  It preserves only priority.  It
destroys the point-line incidence, the three gcd labels on a line, and the
positioned interval discrepancy.  The weighted Fano hypergraph plus the
protected-interval sidecar is therefore retained as the proof carrier.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
from pathlib import Path


PAIR_SCAN = 500
DISCREPANCY_GRID = 84
FANO_LINES = (
    (0, 1, 2),
    (0, 3, 4),
    (0, 5, 6),
    (1, 3, 5),
    (1, 4, 6),
    (2, 3, 6),
    (2, 4, 5),
)


def fold(modulus: int, x: int) -> int:
    r = x % modulus
    return r * (modulus - r)


def pair_density(a: int, b: int) -> Fraction:
    """Haar measure of the two radius-1/14 danger combs (THM-965)."""
    assert 0 < a < b
    g = gcd(a, b)
    modulus = 14 * g
    numerator = (
        4 * a * b
        + fold(modulus, a + b)
        - fold(modulus, b - a)
    )
    return Fraction(numerator, 196 * a * b)


def q_majorant(c: int) -> Fraction:
    return Fraction(c, 1) - Fraction(2, 7) * Fraction(c * (c - 1), 2)


def line_term(k: int) -> Fraction:
    return Fraction(k, 3) - Fraction(2, 7) * Fraction(k * (k - 1), 2)


def main() -> None:
    # Sharp periodic-discrepancy envelope.  Normalize one period to length 1.
    # If a remainder interval has length r and the periodic set has mass rho,
    # its intersection x lies between max(0,r-(1-rho)) and min(r,rho).
    discrepancy_checks = 0
    for denominator in range(1, DISCREPANCY_GRID + 1):
        for occupied in range(denominator + 1):
            rho = Fraction(occupied, denominator)
            sharp = rho * (1 - rho)
            for remainder in range(denominator + 1):
                r = Fraction(remainder, denominator)
                upper_x = min(r, rho)
                lower_x = max(Fraction(0), r - (1 - rho))
                assert upper_x - rho * r <= sharp
                assert rho * r - lower_x <= sharp
                discrepancy_checks += 1

    single_density = Fraction(1, 7)
    single_error = single_density * (1 - single_density)
    pair_error = Fraction(1, 14) * (1 - Fraction(1, 14))
    assert single_error == Fraction(6, 49)
    assert pair_error == Fraction(13, 196)
    # Both constants are sharp for the density-only lemma: take a contiguous
    # occupied block and a remainder interval of the same length.

    # Exact pair range.  The fold formula is dilation-invariant, so primitive
    # A<B suffice.  The all-range proof reduces to tiny finite tails:
    #   lower: fold_plus-fold_minus >= -49, so AB>=27 is automatic;
    #   upper: fold_plus-fold_minus <= 49, so AB>=5 is automatic.
    lower = Fraction(1, 91)
    upper = Fraction(1, 14)
    independence = Fraction(1, 49)
    assert independence - Fraction(1, 4 * 27) >= lower
    assert Fraction(49) <= Fraction(10 * 5)

    lower_tail = []
    upper_tail = []
    for b in range(2, 27):
        for a in range(1, b):
            if gcd(a, b) == 1 and a * b <= 26:
                lower_tail.append((pair_density(a, b), a, b))
    for a, b in ((1, 2), (1, 3), (1, 4)):
        assert gcd(a, b) == 1 and a * b < 5
        upper_tail.append((pair_density(a, b), a, b))
    assert min(lower_tail) == (lower, 1, 13)
    assert max(upper_tail) == (upper, 1, 2)

    pair_checks = 0
    equality_low = []
    equality_high = []
    for b in range(2, PAIR_SCAN + 1):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            rho = pair_density(a, b)
            assert lower <= rho <= upper
            if rho == lower:
                equality_low.append((a, b))
            if rho == upper:
                equality_high.append((a, b))
            pair_checks += 1
    assert equality_low == [(1, 13)]
    assert equality_high == [(1, 2)]

    # Fano incidence and the exact quadratic decomposition, checked on all
    # 2^7 activation masks.  Each point occurs on three lines and every pair
    # occurs on exactly one line.
    point_degrees = Counter(v for line in FANO_LINES for v in line)
    pair_degrees = Counter(
        tuple(sorted(edge))
        for line in FANO_LINES
        for edge in combinations(line, 2)
    )
    assert point_degrees == Counter({v: 3 for v in range(7)})
    assert pair_degrees == Counter(combinations(range(7), 2))

    activation_checks = 0
    majorant_values = {}
    for mask in range(1 << 7):
        c = mask.bit_count()
        line_counts = [sum((mask >> v) & 1 for v in line) for line in FANO_LINES]
        assert sum(line_term(k) for k in line_counts) == q_majorant(c)
        if c:
            assert q_majorant(c) >= 1
        majorant_values[c] = q_majorant(c)
        activation_checks += 1
    assert majorant_values == {
        0: Fraction(0),
        1: Fraction(1),
        2: Fraction(12, 7),
        3: Fraction(15, 7),
        4: Fraction(16, 7),
        5: Fraction(15, 7),
        6: Fraction(12, 7),
        7: Fraction(1),
    }
    # Optimality within C-alpha*binom(C,2): C=7 forces alpha<=2/7.
    assert Fraction(7 - 1, 21) == Fraction(2, 7)

    # Constant ledger for the selected Fano line.  Coverage and the quadratic
    # majorant give one line functional >=L/7.  Its single terms are at most
    # L/7+(2/49)H_line; its three pair terms are at least
    # 3L/91-(13/196)G_line.  Rearrangement is exactly:
    # H_line+(13/28)G_line >= (3/13)L.
    h_coefficient = Fraction(2, 49)
    g_coefficient = Fraction(2, 7) * pair_error
    l_deficit = Fraction(2, 7) * 3 * lower
    assert g_coefficient == Fraction(13, 686)
    assert l_deficit == Fraction(6, 637)
    assert g_coefficient / h_coefficient == Fraction(13, 28)
    assert l_deficit / h_coefficient == Fraction(3, 13)

    # Lower-case LRC fattening supplies L>=1/(7m), so the selected line obeys
    # m[H_line+(13/28)G_line] >=3/91.
    protected_length_at_m_one = Fraction(1, 7)
    line_debt_at_m_one = Fraction(3, 13) * protected_length_at_m_one
    assert line_debt_at_m_one == Fraction(3, 91)

    # Common-dilate corollary.  If the seven deleted speeds are G*y_i with
    # distinct y_i, any line has harmonic mass <=H_3/G and three reciprocal
    # gcds <=3/G.  Combining with the necessary line debt gives G/m<=3523/36.
    harmonic_three = sum((Fraction(1, i) for i in range(1, 4)), Fraction(0))
    common_line_cap = harmonic_three + Fraction(13, 28) * 3
    common_dilate_bound = Fraction(91, 3) * common_line_cap
    assert harmonic_three == Fraction(11, 6)
    assert common_line_cap == Fraction(271, 84)
    assert common_dilate_bound == Fraction(3523, 36)

    # Tournament audit: scalar ordering of the seven line obligations is a
    # deliberately lossy transitive tournament.  The actual carrier is Fano.
    vertices = tuple(range(7))
    edges = {(i, j) for i in vertices for j in vertices if i < j}
    scores = Counter(sum((i, j) in edges for j in vertices if i != j) for i in vertices)
    directed_triangles = sum(
        (i, j) in edges and (j, k) in edges and (k, i) in edges
        for i, j, k in combinations(vertices, 3)
    )
    hamiltonian_paths = sum(
        all((path[i], path[i + 1]) in edges for i in range(6))
        for path in permutations(vertices)
    )
    assert scores == Counter({i: 1 for i in range(7)})
    assert directed_triangles == 0
    assert hamiltonian_paths == 1
    hamiltonian_path = vertices
    assert all((hamiltonian_path[i], hamiltonian_path[i + 1]) in edges for i in range(6))

    print("THM-1166 seven-wall Fano/gcd discrepancy referee")
    print(f"period-cell envelope checks: {discrepancy_checks}")
    print("sharp discrepancy: |mu(I cap pullback(A))-rho*|I|| <= rho(1-rho)/g")
    print(f"single-comb endpoint constant: {single_error}")
    print(f"distinct-pair endpoint constant: {pair_error}")
    print()
    print(f"primitive pair scan: {pair_checks} pairs with b<={PAIR_SCAN}")
    print(f"pair-density range: [{lower}, {upper}]")
    print(f"lower equality: {equality_low}; upper equality: {equality_high}")
    print("all-range tails: AB>=27 closes the lower bound; AB>=5 closes the upper bound")
    print()
    print(f"Fano activation masks checked: {activation_checks}")
    print(f"point degrees: {sorted(point_degrees.items())}")
    print(f"pair degrees: all {len(pair_degrees)} edges occur once")
    print(f"Q(C), C=0..7: {[majorant_values[c] for c in range(8)]}")
    print("optimal quadratic coefficient: 2/7 (forced by C=7)")
    print()
    print("necessary selected-line inequality:")
    print("  H_line + (13/28) sum_pairs(1/gcd) >= (3/13)|I|")
    print("with |I|>=1/(7m):")
    print("  m[H_line + (13/28) sum_pairs(1/gcd)] >= 3/91")
    print(f"common-dilate line cap: {common_line_cap}/G")
    print(f"common-dilate consequence: G/m <= {common_dilate_bound} = {float(common_dilate_bound):.12f}")
    print()
    print("Tournament/carrier audit")
    print("vertices: seven Fano-line proof obligations (not runners)")
    print("pair observable: scalar debt comparison; switch: larger debt first; ties: displayed line order")
    print(f"tie Hamiltonian path: {hamiltonian_path}")
    print(f"score histogram: {sorted(scores.items())}")
    print(f"directed cycles: 0; SCCs: 7 singletons; Hamiltonian paths: {hamiltonian_paths}")
    print("preserves: priority only")
    print("destroys: point-line incidence, three gcd labels, and interval position")
    print("faithful carrier: weighted Fano hypergraph + protected-interval endpoint sidecar")
    print()
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
