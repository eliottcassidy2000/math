#!/usr/bin/env python3
"""Exact referee for THM-2551.

The simultaneous C_13 action on ordered root pairs is free.  Its quotient is
the displacement d=b-h, so the semantic diagonal is one free orbit rather
than a rotation-fixed locus.  Endpoint swap descends to d -> -d.  On a
primitive integer displacement ledger, swap symmetry makes the parity of
the total quotient weight equal to the parity of the zero-displacement weight.
"""

from itertools import combinations

import sympy as sp


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def displacement(pair):
    h, b = pair
    return (b - h) % P


def translate(pair, u):
    h, b = pair
    return ((h + u) % P, (b + u) % P)


def orbit(pair):
    return {translate(pair, u) for u in range(P)}


print("== THM-2551: diagonal C_13 action and displacement quotient ==")
pairs = {(h, b) for h in range(P) for b in range(P)}
orbits = {frozenset(orbit(pair)) for pair in pairs}
require(len(orbits) == P, "wrong number of diagonal-translation orbits")
require(all(len(o) == P for o in orbits), "the action is not free")
require(
    {displacement(next(iter(o))) for o in orbits} == set(range(P)),
    "displacement does not classify the orbits",
)
for u in range(1, P):
    fixed = [pair for pair in pairs if translate(pair, u) == pair]
    require(not fixed, f"nontrivial translation {u} has a fixed pair")
diagonal = {(h, h) for h in range(P)}
require(diagonal in [set(o) for o in orbits], "diagonal is not one orbit")
print(f"  ordered pairs: {len(pairs)}")
print(f"  free orbits: {len(orbits)} of size {P}")
print("  orbit coordinate: d=b-h in F_13")
print("  semantic diagonal: d=0 is free, not rotation-fixed")


print("\n== invariant tables are exactly displacement ledgers ==")
sample_weights = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9]
table = {
    (h, b): sample_weights[(b - h) % P]
    for h in range(P)
    for b in range(P)
}
for u in range(P):
    require(
        all(table[translate(pair, u)] == table[pair] for pair in pairs),
        "displacement table is not invariant",
    )
row_sums = [sum(table[h, b] for b in range(P)) for h in range(P)]
column_sums = [sum(table[h, b] for h in range(P)) for b in range(P)]
require(len(set(row_sums)) == 1, "source margins are not uniform")
require(len(set(column_sums)) == 1, "later-root margins are not uniform")
require(row_sums[0] == sum(sample_weights), "wrong uniform margin")
require(sum(table[h, h] for h in range(P)) == P * sample_weights[0], "wrong hit")
print("  arbitrary 13-entry displacement ledger reconstructed the 169-entry table")
print(f"  common source/later-root margin: {row_sums[0]}")
print(f"  diagonal hit: 13*m_0 = {P * sample_weights[0]}")


print("\n== rotation-blind and Hall-perfect zero-arrival controls ==")
c = 1
single = {(h, (h + c) % P) for h in range(P)}
require(len(single) == P and not (single & diagonal), "single-offset hostile failed")
for h in range(P):
    require(sum((h, b) in single for b in range(P)) == 1, "bad source margin")
for b in range(P):
    require(sum((h, b) in single for h in range(P)) == 1, "bad target margin")

# The endpoint-swap-invariant sharp hostile uses the offsets +/-1.
double = single | {(h, (h - c) % P) for h in range(P)}
require(not (double & diagonal), "endpoint-swap hostile has a diagonal point")
require({(b, h) for h, b in double} == double, "hostile is not endpoint-swap invariant")

# Equal-mass sharp control: two copies of O_0 and O_1 union O_-1 have the
# same translation/swap symmetries and the same margins, but opposite arrival.
aligned = {(h, h): 2 for h in range(P)}
displaced = {pair: 1 for pair in double}
require(sum(aligned.values()) == sum(displaced.values()) == 2 * P, "mass mismatch")
for ledger in (aligned, displaced):
    require(
        [sum(ledger.get((h, b), 0) for b in range(P)) for h in range(P)]
        == [2] * P,
        "source-margin mismatch",
    )
    require(
        [sum(ledger.get((h, b), 0) for h in range(P)) for b in range(P)]
        == [2] * P,
        "target-margin mismatch",
    )
    require(
        all(
            ledger.get(translate(pair, u), 0) == ledger.get(pair, 0)
            for pair in pairs
            for u in range(P)
        ),
        "translation symmetry mismatch",
    )
    require(
        all(ledger.get((b, h), 0) == ledger.get((h, b), 0) for h, b in pairs),
        "endpoint-swap symmetry mismatch",
    )
require(sum(aligned[h, h] for h in range(P)) == 2 * P, "aligned hit mismatch")
require(sum(displaced.get((h, h), 0) for h in range(P)) == 0, "displaced hit mismatch")

hall_checks = 0
roots = list(range(P))
for size in range(P + 1):
    for subset in combinations(roots, size):
        S = set(subset)
        neighbours = {(h + eps * c) % P for h in S for eps in (-1, 1)}
        require(len(neighbours) >= len(S), "off-diagonal Hall inequality failed")
        hall_checks += 1
require(hall_checks == 2**P, "wrong Hall subset count")
print(f"  one-offset hostile: {len(single)} atoms, odd orbit weight 1, hit 0")
print(f"  endpoint-swap hostile: {len(double)} atoms, primitive orbit weight 2, hit 0")
print("  equal-mass symmetric controls: aligned hit 26, displaced hit 0")
print(f"  off-diagonal Hall checks: {hall_checks}/{2**P}")


print("\n== exact mixed-character control for the +/-1 hostile ==")
X = sp.symbols("X")
phi = sum(X**j for j in range(P))
nonzero_modes = 0
for a in range(1, P):
    mode = X**a + X**((-a) % P)
    remainder = sp.rem(mode, phi, domain=sp.QQ)
    require(remainder != 0, f"mixed character {a} vanished")
    nonzero_modes += 1
require(nonzero_modes == P - 1, "wrong nontrivial mode count")
print(f"  relative mixed modes nonzero: {nonzero_modes}/{P-1}")
print("  both one-root marginals are uniform; all their nontrivial modes vanish")


print("\n== primitive endpoint-swap parity ==")
parity_checks = 0
# A swap-invariant ledger is (n_0,n_1,...,n_6,n_6,...,n_1).
for n0 in range(4):
    # Exhaust every six-tuple with entries in {0,1,2}; 4*3^6=2916 controls.
    for code in range(3**6):
        q = code
        half = []
        for _ in range(6):
            half.append(q % 3)
            q //= 3
        total = n0 + 2 * sum(half)
        require(total % 2 == n0 % 2, "endpoint-swap parity identity failed")
        if total % 2 == 1:
            require(n0 > 0, "odd swap-invariant ledger has no zero displacement")
        parity_checks += 1
require(parity_checks == 4 * 3**6, "wrong parity control count")
print(f"  endpoint-symmetric primitive-ledger checks: {parity_checks}")
print("  total quotient weight mod 2 = zero-displacement weight mod 2")
print("  odd primitive weight => positive semantic diagonal")


print("\nall checks passed")
