#!/usr/bin/env python3
"""Exact audit for THM-1099 Route-B affine/lift guardrails.

All decisions use integer arithmetic or fractions.Fraction.  The exact M
evaluator uses the proved pair-sum-ruler fact: a global maximum occurs at a
reduced denominator dividing the sum of two active speeds.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt


def divisors(n: int) -> set[int]:
    out: set[int] = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.add(d)
            out.add(n // d)
    return out


def clearance_numerator(v: int, a: int, q: int) -> int:
    r = (v * a) % q
    return min(r, q - r)


def exact_M(values: tuple[int, ...]) -> tuple[F, tuple[int, int, int]]:
    rulers: set[int] = set()
    for i, x in enumerate(values):
        for y in values[i:]:
            rulers.update(divisors(x + y))

    best = F(0)
    owner = (0, 1, 0)
    for q in sorted(rulers):
        if q < 2:
            continue
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            val = min(clearance_numerator(v, a, q) for v in values)
            candidate = F(val, q)
            if candidate > best:
                best = candidate
                owner = (a, q, val)
    return best, owner


def covers(values: tuple[int, ...], top: int) -> bool:
    return all(any(v % d == 0 for v in values) for d in range(2, top + 1))


def sorted_residues(values: tuple[int, ...], a: int, q: int) -> tuple[int, ...]:
    return tuple(sorted((a * v) % q for v in values))


def internal_gaps(residues: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(y - x for x, y in zip(residues, residues[1:]))


def least_lift(base: int, q: int, d: int) -> int:
    """Least positive base+kq divisible by d for the positive bases used here."""
    period = d // gcd(q, d)
    candidates = [base + k * q for k in range(period) if (base + k * q) % d == 0]
    if not candidates:
        raise ValueError((base, q, d))
    return min(x for x in candidates if x > 0)


def least_owner_cost(bases: tuple[int, ...], q: int, d: int) -> int:
    possible = [
        least_lift(u, q, d)
        for u in bases
        if u % gcd(q, d) == 0
    ]
    return min(possible)


def tournament_telemetry(cost: dict[int, int]) -> tuple[tuple[int, ...], tuple[int, ...], int, int, int]:
    vertices = tuple(sorted(cost))

    def edge(x: int, y: int) -> bool:
        return (cost[x], x) < (cost[y], y)

    scores = tuple(sum(edge(x, y) for y in vertices if y != x) for x in vertices)
    path = tuple(sorted(vertices, key=lambda x: (cost[x], x)))
    cycles = sum(
        (edge(x, y) and edge(y, z) and edge(z, x))
        or (edge(x, z) and edge(z, y) and edge(y, x))
        for x, y, z in combinations(vertices, 3)
    )
    # A strict total-order tournament has singleton SCCs and one Hamiltonian path.
    scc_count = len(vertices) if cycles == 0 else -1
    hamiltonian_path_count = 1 if cycles == 0 else -1
    assert all(edge(path[i], path[i + 1]) for i in range(len(path) - 1))
    return scores, path, cycles, scc_count, hamiltonian_path_count


print("THM-1099 exact Route-B affine/lift guardrail audit")

# 1. Threshold trichotomy, exhaustively over a deterministic finite box.
threshold_rows = 0
for s in range(1, 81):
    for e in range(1, 2 * s + 1):
        q = 13 * s + e
        value = F(s, q)
        assert (F(1, 14) < value < F(1, 13)) == (0 < e < s)
        assert (value == F(1, 14)) == (e == s)
        assert (value < F(1, 14)) == (e > s)
        threshold_rows += 1
print(f"threshold trichotomy: PASS ({threshold_rows} exact (s,e) rows)")

# Signed observer+runner gap identity on representative legal packets.
packet_rows = 0
for s in range(2, 31):
    for e in range(1, s):
        q = 13 * s + e
        # A deterministic 13-point packet in [s,q-s].  It need not arise from
        # a speed family; the identity is representation-theoretic.
        residues = tuple(range(s, s + 13))
        assert residues[-1] <= q - s
        augmented = (0,) + residues + (q,)
        gaps = tuple(y - x for x, y in zip(augmented, augmented[1:]))
        assert sum(g - s for g in gaps) == q - 14 * s == e - s
        packet_rows += 1
print(f"fourteen-gap signed defect: PASS ({packet_rows} exact packets)")

# 2. The e=1 packing carry family.
packing_rows = 0
for s in range(2, 81):
    pts = tuple(i * s for i in range(1, 7)) + tuple(i * s + 1 for i in range(7, 13))
    gaps = internal_gaps(pts)
    assert pts[0] == s and pts[-1] == 12 * s + 1
    assert all(g >= s for g in gaps)
    assert gaps.count(s + 1) == 1
    assert pts != tuple(i * s for i in range(1, 13))
    packing_rows += 1
print(f"one-carry packing counterfamily: PASS ({packing_rows} exact scales)")

# 3. Coprime multiples and affine carry.
coset_rows = 0
carry_rows = 0
for s in range(2, 41):
    q = 13 * s + 1
    assert gcd(s, q) == 1
    assert len({(s * x) % q for x in range(q)}) == q
    coset_rows += 1

for s in range(2, 24):
    for e in range(1, s):
        q = 13 * s + e
        for x in range(q):
            for y in range(q):
                wrap = int(x + y >= q)
                z = (x + y) % q
                assert z % s == (x % s + y % s - wrap * e) % s
                carry_rows += 1
print(f"coprime multiple range = whole cyclic group: PASS ({coset_rows} scales)")
print(f"affine carry law: PASS ({carry_rows} exact additions)")

# 4-5. The exact death-star-S57 lift pair.
V0 = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12, 17, 19, 104)
V1 = (1, 2, 3, 5, 8, 9, 10, 11, 12, 17, 19, 104, 112)

M0, owner0 = exact_M(V0)
M1, owner1 = exact_M(V1)
assert M0 == F(8, 105)
assert M1 == F(3, 20)
assert covers(V0, 13) and not covers(V0, 14)
assert covers(V1, 14)

R0 = sorted_residues(V0, 8, 105)
R1 = sorted_residues(V1, 8, 105)
G0 = internal_gaps(R0)
assert R0 == R1
assert R0 == (8, 16, 24, 31, 40, 47, 56, 64, 72, 80, 88, 96, 97)
assert G0 == (8, 8, 7, 9, 7, 9, 8, 8, 8, 8, 8, 1)
assert (8 * 7) % 105 == (8 * 112) % 105 == 56

clearances = {v: clearance_numerator(v, 8, 105) for v in V1}
assert clearances[1] == clearances[104] == 8
assert min(c for v, c in clearances.items() if v not in (1, 104)) == 9

at_3_20 = tuple(clearance_numerator(v, 3, 20) for v in V1)
assert at_3_20 == (3, 6, 9, 5, 4, 7, 10, 7, 4, 9, 3, 8, 4)
assert min(at_3_20) == 3
assert F(3, 20) > F(1, 13)

assert gcd(105, 14) == 7 and 7 % 7 == 0
assert least_lift(7, 105, 14) == 112

print(f"V0 exact M: {M0}, maximizer a/q/val={owner0}")
print(f"V0 cover through 13 / 14: {covers(V0,13)} / {covers(V0,14)}")
print(f"V1 exact M: {M1}, maximizer a/q/val={owner1}")
print(f"V1 cover through 14: {covers(V1,14)}")
print(f"shared residues at 8/105: {R0}")
print(f"shared internal gaps: {G0}")
print("active-pair slack after 7->112: active=8, next=9 (numerators over 105)")
print(f"3/20 witness distance numerators: {at_3_20}")
print("CRT owner lift: gcd(105,14)=7 and least 14-carrier in 7+105Z is 112")

# A residue AP does not imply an integer AP.
assert (14 * 2) % 183 == (14 * 185) % 183 == 28
print("modular-lift guardrail: speeds 2 and 185 have the same residue 28 at 14/183")

# 7. Owner-cost tournament telemetry on alternate vertices d=2,...,14.
bases = V0  # all lie in [0,105), hence are canonical lift bases.
cost = {d: least_owner_cost(bases, 105, d) for d in range(2, 15)}
scores, path, cycles, scc_count, hp_count = tournament_telemetry(cost)
owners = tuple(range(2, 15))
score_hist = {x: scores.count(x) for x in sorted(set(scores))}
print(f"owner lift costs: {tuple((d,cost[d]) for d in owners)}")
print(f"owner-cost tournament path: {path}")
print(f"owner-cost score histogram: {score_hist}")
print(f"owner-cost directed 3-cycles / SCCs / Hamiltonian paths: {cycles} / {scc_count} / {hp_count}")

# The proof-bearing incidence remembers shared carriers, unlike the tournament.
selected_carrier = {
    d: min(v for v in V1 if v % d == 0)
    for d in range(2, 15)
}
shared_owner_pairs = sum(
    selected_carrier[d] == selected_carrier[e]
    for d, e in combinations(range(2, 15), 2)
)
print(f"V1 owner->least-carrier incidence: {tuple(sorted(selected_carrier.items()))}")
print(f"owner pairs sharing their least carrier: {shared_owner_pairs}")

print("ALL EXACT CHECKS PASS")
