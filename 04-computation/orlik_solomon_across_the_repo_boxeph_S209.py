#!/usr/bin/env python3
"""Exact controls for the corrected HYP-8830 arrangement wildcard.

This script deliberately separates four objects that S209 conflated:

* the ordinary braid and Shi arrangements;
* the Fourier annihilator L_v = {k in Z^n : k.v = 0};
* a finite toric arrangement and its layer poset; and
* the labelled phase-height walls v*t +/- delta in Z used by LRC.

It proves only the finite/combinatorial assertions printed below.  In
particular it computes no arithmetic Mobius formula for |G_delta|.
"""

from itertools import combinations, product
from math import factorial, gcd
from fractions import Fraction


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def stirling1_unsigned(n: int) -> list[int]:
    row = [1]
    for m in range(1, n + 1):
        nxt = [0] * (m + 1)
        for j in range(1, m + 1):
            nxt[j] = row[j - 1] + (m - 1) * (row[j] if j < len(row) else 0)
        row = nxt
    return row


def braid_poincare(n: int) -> list[int]:
    coefficients = [1]
    for k in range(1, n):
        coefficients = [
            (coefficients[i] if i < len(coefficients) else 0)
            + k * (coefficients[i - 1] if i else 0)
            for i in range(len(coefficients) + 1)
        ]
    return coefficients


section("P1  Classical braid-arrangement control")
for n in range(2, 8):
    poincare = braid_poincare(n)
    stirling = stirling1_unsigned(n)
    expected = [stirling[n - i] for i in range(n)]
    assert poincare == expected
    print(f"n={n}: Poincare={poincare}; top={(poincare[-1])}={(factorial(n - 1))}")
print("Scope: this is cohomology of one braid complement, not an invariant of each tournament.")


def count_braid_complement(n: int, q: int) -> int:
    return sum(len(set(point)) == n for point in product(range(q), repeat=n))


def count_shi_complement(n: int, q: int) -> int:
    total = 0
    for point in product(range(q), repeat=n):
        if all(
            (point[i] - point[j]) % q not in (0, 1)
            for i in range(n)
            for j in range(i + 1, n)
        ):
            total += 1
    return total


section("P2  Classical finite-field braid and Shi controls")
for n in (3, 4):
    for q in (11, 13):
        braid = count_braid_complement(n, q)
        braid_formula = 1
        for i in range(n):
            braid_formula *= q - i
        shi = count_shi_complement(n, q)
        shi_formula = q * (q - n) ** (n - 1)
        assert braid == braid_formula
        assert shi == shi_formula
        print(f"n={n}, q={q}: braid={braid}; Shi={shi}")
print("Scope: these textbook point counts do not identify Shi walls with LRC safety walls.")


section("P3  Minimal witness: a thickened safe set is not a standard toric complement")
delta = Fraction(1, 5)
safe_measure = 1 - 2 * delta
ordinary_complement_measure = Fraction(1)
assert safe_measure == Fraction(3, 5)
print(f"S={{1}}, delta={delta}: |G_delta|={safe_measure}")
print("ordinary complement of ker(z -> z) in the circle: Haar measure=1")
print("Therefore |G_delta| is not ordinary toric-complement volume.")


def det(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[1] - b[0] * a[1]


section("P4  Exact signed phase-height wall arithmetic")
speeds = (1, 2, 3)
for v in speeds:
    assert abs(det((v, 1), (v, -1))) == 2 * v
    print(f"self pair v={v}: |det((v,+1),(v,-1))|={2 * v}")
for v, w in combinations(speeds, 2):
    same = abs(det((v, 1), (w, 1)))
    opposite = abs(det((v, 1), (w, -1)))
    assert same == w - v
    assert opposite == v + w
    print(f"v={v}, w={w}: same-sign={same}=|v-w|; opposite-sign={opposite}=v+w")
print("Full-torus determinants expose pair-sum arithmetic; the top-slope proof selects the ruler.")


def relation_richness(v: tuple[int, int, int], bound: int = 2) -> int:
    count = sum(
        any(k) and sum(k[i] * v[i] for i in range(3)) == 0
        for k in product(range(-bound, bound + 1), repeat=3)
    )
    return count // 2


def roots_of_unity_union_size(speeds: tuple[int, ...]) -> int:
    roots = {
        Fraction(a, v) % 1
        for v in speeds
        for a in range(v)
    }
    return len(roots)


section("P5  Finite B=2 scout and a hostile standard-toric control")
triples = [
    (a, b, c)
    for a in range(1, 7)
    for b in range(a + 1, 8)
    for c in range(b + 1, 10)
    if gcd(gcd(a, b), c) == 1
]
assert len(triples) == 72
scored = sorted(((relation_richness(v), v) for v in triples), reverse=True)
leaders = [item for item in scored if item[0] == scored[0][0]]
assert leaders == [(4, (1, 2, 3))]
print(f"hard-coded universe: {len(triples)} primitive triples; unique B=2 leader={leaders[0]}")
ap_layers = roots_of_unity_union_size((1, 2, 3))
control_layers = roots_of_unity_union_size((2, 3, 4))
assert (ap_layers, control_layers) == (4, 6)
print(f"natural circle hypertori union: (1,2,3) has {ap_layers} points; (2,3,4) has {control_layers}")
print("Thus the B=2 relation count is not Betti/Mobius mass and is cutoff-dependent.")


section("SUMMARY")
print(
    "Classical braid/Shi formulas survive.  THM-1820's Fourier annihilator identity survives.\n"
    "The live LRC arrangement is the labelled, oriented phase-height cell complex for\n"
    "X_S={(v,+1),(v,-1)} with its selected inequality side and height functional.\n"
    "Ordinary complement cohomology is not known to recover the selected state;\n"
    "abstract relation lattices and bounded B=2 relation data are provably lossy."
)
