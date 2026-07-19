#!/usr/bin/env python3
"""Exact referee for THM-1247 q=15 contracted-Fano guardrail.

Checks the sign-address mask classification, invariant Fano planes and common
contraction, CRT independence of speed chi_7, and the blocker-complete packet
``{1,2,7,8,18,19,20}`` using exact finite arithmetic.  Explicit failures are
used instead of ``assert``.
"""

from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: F) -> F:
    residue = x - (x.numerator // x.denominator)
    return min(residue, 1 - residue)


def mask(speed: int, q: int = 15) -> frozenset[int]:
    return frozenset(
        p for p in range(q)
        if 14 * min(speed * p % q, q - (speed * p % q)) < q
    )


def sign_address(residue: int) -> int:
    residue %= 15
    return min(residue, 15 - residue)


# Six masks and the seven-point sign-address refinement.
all_masks = {mask(s) for s in range(1, 15) if s % 15}
expected_masks = {
    frozenset((0, r, 15 - r)) for r in (1, 2, 4, 5, 7)
} | {frozenset((0, 3, 6, 9, 12))}
require(all_masks == expected_masks and len(all_masks) == 6, "six mask classes")

petal = {1: "A", 2: "B", 3: "C", 4: "D", 5: "E", 6: "C", 7: "F"}
cyclic_word = tuple(petal[sign_address(r)] for r in range(1, 15))
require(cyclic_word == tuple("ABCDE CFFC EDCBA".replace(" ", "")),
        "palindromic cyclic petal word")

# All labelled Fano planes and the two invariant under multiplication by two.
F0 = frozenset(
    frozenset(line) for line in ((1, 2, 3), (1, 4, 5), (1, 6, 7),
                                 (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 5, 6))
)
F1 = frozenset(
    frozenset(line) for line in ((1, 2, 6), (1, 3, 7), (1, 4, 5),
                                 (2, 3, 4), (2, 5, 7), (3, 5, 6), (4, 6, 7))
)
planes = set()
for perm in permutations(range(1, 8)):
    relabel = {i + 1: perm[i] for i in range(7)}
    planes.add(frozenset(frozenset(relabel[x] for x in line) for line in F0))
require(len(planes) == 30, "labelled Fano plane count")

T = {r: sign_address(2 * r) for r in range(1, 8)}


def transport(plane):
    return frozenset(frozenset(T[x] for x in line) for line in plane)


invariant = {plane for plane in planes if transport(plane) == plane}
require(invariant == {F0, F1}, "two multiplication-by-two invariant planes")
require(frozenset((3, 5, 6)) in F0 and frozenset((3, 5, 6)) in F1,
        "fixed negative line")

symbol = {1: "A", 2: "B", 3: "C", 4: "D", 5: "E", 6: "C", 7: "F"}


def contracted(plane):
    return sorted(tuple(sorted(symbol[x] for x in line)) for line in plane)


contracted_target = sorted((
    ("A", "B", "C"), ("A", "D", "E"), ("A", "C", "F"),
    ("B", "C", "D"), ("B", "E", "F"), ("C", "D", "F"),
    ("C", "C", "E"),
))
require(contracted(F0) == contracted(F1) == contracted_target,
        "common contracted Fano multihypergraph")


def chi7(speed: int) -> int:
    while speed % 7 == 0:
        speed //= 7
    residue = speed % 7
    return 1 if residue in (1, 2, 4) else -1


# Defining difference pair supports all four ordered chi_7 signs.
pair_base = {}
for signs in product((-1, 1), repeat=2):
    candidates = [a for a in range(1, 106)
                  if a % 15 == 1 and (chi7(a), chi7(a + 15)) == signs]
    require(candidates, f"CRT defining pair signs {signs}")
    pair_base[signs] = candidates[0]

other_residues = (5, 3, 2, 4, 7)
other_base = {}
for residue in other_residues:
    for sign in (-1, 1):
        candidates = [x for x in range(1, 106)
                      if x % 15 == residue and x % 7 != 0 and chi7(x) == sign]
        require(candidates, f"CRT mask residue {residue}, sign {sign}")
        other_base[(residue, sign)] = candidates[0]

colour_words = 0
for signs in product((-1, 1), repeat=7):
    a = pair_base[signs[:2]]
    speeds = [a, a + 15]
    for j, (residue, sign) in enumerate(zip(other_residues, signs[2:])):
        speeds.append(other_base[(residue, sign)] + 105 * (j + 2))
    require(len(set(speeds)) == 7, "CRT distinct speeds")
    require(tuple(chi7(v) for v in speeds) == signs, "CRT colour word")
    require(len({mask(v) for v in speeds}) == 6 and
            set().union(*(set(mask(v)) for v in speeds)) == set(range(15)),
            "CRT preserves sunflower")
    colour_words += 1

# Blocker-complete packet.
V = (1, 2, 7, 8, 18, 19, 20)
D = (2, 7, 8, 18, 19, 20)
q15_masks = tuple(mask(v) for v in V)
require(len(set(q15_masks)) == 6 and set().union(*(set(M) for M in q15_masks)) == set(range(15)),
        "witness full q15 clock")
require(mask(7) == mask(8), "defining pair mask equality")

beat_blocker = {2: 7, 13: 7, 3: 20, 6: 20, 9: 20, 12: 20,
                4: 19, 11: 19, 5: 18, 10: 18, 7: 2, 8: 2}
for p in range(2, 14):
    t = F(p, 15)
    blocker = beat_blocker[p]
    require(circle_norm(blocker * t) < F(1, 14), "q15 third blocker")
    if p not in (2, 13):
        require(circle_norm(7 * t) >= F(1, 14) and circle_norm(8 * t) >= F(1, 14),
                "defining pair safe beat")

invoice = tuple(sum(abs(di - dh) for di in D) for dh in D)
require(invoice == (62, 42, 40, 40, 42, 46), "all-pivot drift invoice")
require(sum(20 - di for di in D) == 46 and 14 * 46 > 20,
        "fastest pivot invoice")
require(all(210 * (D[i + 1] - D[i]) > 20 for i in range(5)),
        "macroscopic adjacent cuts")


def nearest_half_up(x: F) -> int:
    lo = x.numerator // x.denominator
    return lo if x - lo < F(1, 2) else lo + 1


spoke_times = {}
spoke_blockers = {}
for d in D:
    q = 1 + d
    p = nearest_half_up(F(q, 2))
    t = F(p, q)
    spoke_times[d] = t
    require(F(1, 14) < t < F(13, 14), "carrier spoke in gap")
    require(circle_norm(t) > F(1, 4) and circle_norm(d * t) > F(1, 4),
            "deep carrier spoke")
    blockers = frozenset(e for e in D if e != d and circle_norm(e * t) < F(1, 14))
    require(blockers, "carrier spoke blocked")
    spoke_blockers[d] = blockers

expected_spoke_times = {2: F(2, 3), 7: F(1, 2), 8: F(5, 9),
                        18: F(10, 19), 19: F(1, 2), 20: F(11, 21)}
expected_spoke_blockers = {
    2: frozenset((18,)), 7: frozenset((2, 8, 18, 20)),
    8: frozenset((18,)), 18: frozenset((2, 19)),
    19: frozenset((2, 8, 18, 20)), 20: frozenset((2, 19)),
}
require(spoke_times == expected_spoke_times, "carrier spoke phases")
require(spoke_blockers == expected_spoke_blockers, "carrier spoke blocker sets")

sum_beat_rows = (
    ((2, 7), 1, 18, 19), ((2, 8), 1, 20, 18), ((2, 18), 2, 20, 36),
    ((2, 19), 3, 7, 63), ((2, 20), 3, 7, 62), ((7, 8), 3, 20, 69),
    ((7, 18), 3, 8, 31), ((7, 19), 13, 2, 156), ((7, 20), 3, 18, 57),
    ((8, 18), 11, 7, 114), ((8, 19), 3, 18, 15), ((8, 20), 3, 19, 28),
    ((18, 19), 5, 7, 187), ((18, 20), 4, 19, 18), ((19, 20), 5, 8, 199),
)
for (x, y), p, blocker, expected_curvature in sum_beat_rows:
    q = x + y
    t = F(p, q)
    require(F(1, 14) < t < F(13, 14), "sum beat in gap")
    dx = min(x * p % q, q - (x * p % q))
    dy = min(y * p % q, q - (y * p % q))
    require(dx == dy and 14 * dx - q == expected_curvature > 0,
            "positive sum-beat curvature")
    require(circle_norm(blocker * t) < F(1, 14), "sum-beat third blocker")

lonely_distances = tuple(circle_norm(F(v, 12)) for v in V)
require(lonely_distances == (F(1, 12), F(1, 6), F(5, 12), F(1, 3),
                             F(1, 2), F(5, 12), F(1, 3)), "lonely ledger")
require(min(lonely_distances) == F(1, 12) > F(1, 14), "lonely witness")

require(all(((a + b) // gcd(a, b)) % 14 for a, b in combinations(V, 2)),
        "empty exact-seam graph")

# Cross-clearance tournament on carrier-spoke obligations.
edges = set()
for i, j in combinations(D, 2):
    cross_ij = circle_norm(j * spoke_times[i])
    cross_ji = circle_norm(i * spoke_times[j])
    winner = i if cross_ij < cross_ji or (cross_ij == cross_ji and i < j) else j
    edges.add((winner, j if winner == i else i))

scores = tuple(sorted(sum((i, j) in edges for j in D if j != i) for i in D))
triangles = sum(all((path[k], path[(k + 1) % 3]) in edges for k in range(3))
                for path in permutations(D, 3)) // 3
hamiltonian_paths = sum(all((path[k], path[k + 1]) in edges for k in range(5))
                        for path in permutations(D))
reach = {(i, j): i == j or (i, j) in edges for i in D for j in D}
for k in D:
    for i in D:
        for j in D:
            reach[i, j] = reach[i, j] or reach[i, k] and reach[k, j]
unseen = set(D)
scc = []
while unseen:
    root = min(unseen)
    component = {j for j in unseen if reach[root, j] and reach[j, root]}
    scc.append(len(component))
    unseen -= component
require(scores == (1, 1, 2, 3, 3, 5) and triangles == 3 and
        sorted(scc, reverse=True) == [5, 1] and hamiltonian_paths == 9,
        "cross-clearance tournament fingerprint")

print("THM-1247 q=15 CONTRACTED FANO / BLOCKER-COMPLETE GUARDRAIL EXACT AUDIT")
print("sign-address points = 7; mask classes = 6; contraction = P3~P6")
print(f"labelled Fano planes = {len(planes)}; multiplication-by-two invariant = {len(invariant)}")
print("common contracted negative line = C-E-C")
print(f"independent speed-chi7 colour words realized = {colour_words}")
print("guardrail packet = {1,2,7,8,18,19,20}")
print("q15 beat blockers = 12/12; carrier-spoke blockers = 6/6")
print(f"positive fast-fast sum-beat blockers = {len(sum_beat_rows)}/15")
print("explicit lonely witness = t=1/12, depth=1/12")
print(f"cross-clearance tournament = scores {scores}, triangles {triangles}, SCC {tuple(sorted(scc, reverse=True))}, Hamiltonian paths {hamiltonian_paths}")
print("RESULT: PASS")
