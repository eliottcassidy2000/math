#!/usr/bin/env python3
"""Scoped exact checks for the corrected boxeph-2026-07-21-S208.

This script verifies only the algebraic shadows that survived MISTAKE-223:

1. the braid-arrangement characteristic-polynomial region count;
2. a leading Vandermonde coalescence formula, distinguishing vanishing order
   from flat codimension;
3. det(I-x*M_g)=1-x-x^(g+1) for an explicit companion matrix; and
4. equality of its recurrence coefficients with the full-rank binomial
   gap-diagonal formula.

It does not verify a hyper-Bessel product, an Euler characteristic for bagel
cuts, or an identification of bagel and finite-shadow deficits.
"""

from fractions import Fraction as F
from itertools import permutations
from math import comb, factorial


def section(title):
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def vandermonde(values):
    out = F(1)
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            out *= values[j] - values[i]
    return out


section("CHECK 1  braid defining polynomial and REAL chamber count")
for n in range(2, 8):
    chi_minus_one = 1
    for j in range(1, n):
        chi_minus_one *= -1 - j
    regions = (-1) ** (n - 1) * chi_minus_one
    print(
        f"  n={n}: (-1)^(n-1)*chi(-1)={regions}; "
        f"n!={factorial(n)}; equal? {regions == factorial(n)}"
    )
print("  V(a)=prod_(i<j)(a_j-a_i) vanishes exactly on braid hyperplanes.")
print("  Scope: n! counts chambers of the REAL complement; the complex complement is connected.")


section("CHECK 2  leading coalescence; vanishing order is NOT codimension")
blocks = [[0, 1, 2], [3, 4]]
centres = [0, 10]
deltas = [[0, 1, 3], [0, 2]]
vanishing_order = sum(comb(len(block), 2) for block in blocks)
flat_codimension = sum(len(block) - 1 for block in blocks)
within = vandermonde([F(v) for v in deltas[0]]) * vandermonde(
    [F(v) for v in deltas[1]]
)
between = F(centres[1] - centres[0]) ** (len(blocks[0]) * len(blocks[1]))
for epsilon in (F(1, 10), F(1, 100), F(1, 1000)):
    values = []
    for centre, offsets in zip(centres, deltas):
        values.extend(F(centre) + epsilon * F(offset) for offset in offsets)
    leading = epsilon**vanishing_order * within * between
    ratio = vandermonde(values) / leading
    print(f"  epsilon={float(epsilon):.3g}: V/leading={float(ratio):.9f}")
print(f"  block sizes 3+2: vanishing order={vanishing_order}, codimension={flat_codimension}")
print("  The cross-block factor is a nonzero unit near this generic flat.")


def poly_add(left, right, scale=1):
    out = left[:] + [0] * max(0, len(right) - len(left))
    for i, value in enumerate(right):
        out[i] += scale * value
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def permutation_sign(perm):
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return -1 if inversions % 2 else 1


def companion_det_polynomial(g):
    """Return coefficients of det(I-x*M_g), low degree first."""
    size = g + 1
    matrix = [[0] * size for _ in range(size)]
    matrix[0][0] = 1
    matrix[0][g] = 1
    for row in range(1, size):
        matrix[row][row - 1] = 1

    polynomial_matrix = []
    for i in range(size):
        row = []
        for j in range(size):
            entry = [1 if i == j else 0]
            if matrix[i][j]:
                entry = poly_add(entry, [0, 1], scale=-matrix[i][j])
            row.append(entry)
        polynomial_matrix.append(row)

    determinant = [0]
    for perm in permutations(range(size)):
        term = [1]
        for row, column in enumerate(perm):
            term = poly_mul(term, polynomial_matrix[row][column])
        determinant = poly_add(determinant, term, scale=permutation_sign(perm))
    return determinant


def recurrence(g, length):
    values = [0] * length
    values[0] = 1
    for d in range(1, length):
        values[d] = values[d - 1] + (values[d - g - 1] if d >= g + 1 else 0)
    return values


def gap_binomial(g, d):
    return sum(comb(d - g * k, k) for k in range(d // (g + 1) + 1))


section("CHECK 3  explicit companion determinant and full-rank gap diagonal")
for g in range(1, 6):
    determinant = companion_det_polynomial(g)
    expected = [1, -1] + [0] * (g - 1) + [-1]
    rec = recurrence(g, 13)
    diagonal = [gap_binomial(g, d) for d in range(13)]
    print(
        f"  g={g}: det coeffs={determinant}; expected={expected}; "
        f"det PASS? {determinant == expected}; diagonal PASS? {rec == diagonal}"
    )
print("  The companion digraph has a loop; it is not a tournament adjacency matrix.")


section("NON-RESULTS MADE EXPLICIT")
print("  - No general NC2 iff follows from the special THM-2033 Vandermonde matrix.")
print("  - No hyper-Bessel block product follows from the coalescence formula.")
print("  - THM-2023 proves Phi_(p,q) Laguerre--Polya independently.")
print("  - No bagel-cut Euler complex or bagel/shadow valuation is computed here.")
print("  - bagel(n)=cake(n+1)-2 remains the strongest direct cutting identity.")
