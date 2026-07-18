#!/usr/bin/env python3
"""Exact referee for the corrected THM-1110 numerator-counting lemma.

The incoming count used only dangerous *unit* residues even when a speed was
not a unit modulo q.  This audit computes every unit-numerator kill set and
checks the gcd-stratified fibre formula exactly.
"""

from math import gcd


def phi(n: int) -> int:
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def units(q: int) -> tuple[int, ...]:
    return tuple(a for a in range(1, q) if gcd(a, q) == 1)


def forbidden(q: int) -> tuple[int, ...]:
    return tuple(r for r in range(q) if 14 * min(r, q - r) < q)


def killed(q: int, speed: int) -> tuple[int, ...]:
    bad = set(forbidden(q))
    return tuple(a for a in units(q) if (speed * a) % q in bad)


def gcd_stratum_count(q: int, g: int) -> int:
    """Exact |{a in U_q : v*a in W_q}| for gcd(v,q)=g<q."""
    quotient = q // g
    residue_count = sum(gcd(r, q) == g for r in forbidden(q))
    return phi(q) // phi(quotient) * residue_count


print("THM-1110 GCD-STRATIFIED NUMERATOR-COUNT REFEREE")

# Exhaustively verify that the kill count depends only on gcd(v,q), and equals
# the unit-reduction fibre formula, over a deterministic broad box.
rows = 0
for q in range(2, 241):
    for speed in range(1, q):
        g = gcd(speed, q)
        assert len(killed(q, speed)) == gcd_stratum_count(q, g)
        rows += 1
print(f"gcd-stratified fibre formula: PASS ({rows} exact (q,v) rows)")

q = 90
U = units(q)
W = forbidden(q)
print(f"q={q}: phi={len(U)}, forbidden={W}")

strata = tuple(
    (g, gcd_stratum_count(q, g))
    for g in range(1, q)
    if q % g == 0
)
print(f"proper gcd-stratum kill counts: {strata}")
assert dict(strata)[1] == 2
assert max(count for _, count in strata) == 8

# Three gcd-5 residues partition all unit numerators; adding the redundant unit
# speed 1 makes the actual speed set primitive without restoring a witness.
core = (5, 25, 35)
primitive = (1,) + core
kill_sets = {speed: killed(q, speed) for speed in primitive}
for speed in primitive:
    print(
        f"speed={speed:2d} gcd={gcd(speed,q):2d} "
        f"kills={kill_sets[speed]}"
    )

assert set().union(*(set(kill_sets[v]) for v in core)) == set(U)
assert gcd(*primitive) == 1
assert all(v % q for v in primitive)
assert set().union(*(set(kill_sets[v]) for v in primitive)) == set(U)
print("primitive counterexample: {1,5,25,35} blocks every unit numerator mod 90")

# The corrected union bound uses each speed's own gcd stratum.  Its unit-speed
# specialization recovers the incoming s*k_{q,1}<phi(q) statement.
unit_budget = 11 * gcd_stratum_count(q, 1)
assert unit_budget == 22 < phi(q)
arbitrary_three_budget = 3 * max(count for _, count in strata)
assert arbitrary_three_budget == phi(q)
print(f"all-unit specialization at q=90: 11*2={unit_budget}<24 (valid)")
print("arbitrary nonzero residues at q=90: 3*8=24 (sharp; no strict slack)")
print("ALL EXACT CHECKS PASS")
