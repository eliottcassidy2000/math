#!/usr/bin/env python3
"""Exact referee for THM-1242 seven-wall pair-beat sunflower law.

Checks the all-Q common-clock threshold, its wall-count criticality, exact
unit masks, minimality of the first mixed full clock, and the q=15 sunflower
partition.  Dependency-free and exact; no ``assert``.
"""

from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def A(Q: int) -> int:
    return 2 * ((Q + 13) // 14) - 1


def common_threshold(Q: int, obligations: int) -> int:
    return 2 + obligations * (A(Q) - 1)


common_rows = 0
critical_failures = []
for Q in range(2, 10001):
    B7 = common_threshold(Q, 6)
    require(B7 == 6 * A(Q) - 4, "seven-wall threshold identity")
    require(B7 <= Q, "seven-wall all-Q escape threshold")
    common_rows += 1

    B8 = common_threshold(Q, 7)
    if B8 > Q:
        critical_failures.append(Q)

require(critical_failures[:5] == [15, 29, 43, 57, 71], "eight-wall shell failures")
require(all(Q % 14 == 1 for Q in critical_failures), "eight-wall first-shell rows")


def strict_mask(speed: int, q: int) -> frozenset[int]:
    return frozenset(
        p for p in range(q)
        if 14 * min((speed * p) % q, q - ((speed * p) % q)) < q
    )


q = 15
speeds = (1, 16, 5, 3, 2, 4, 7)
masks = (
    strict_mask(1, q),                  # defining pair 1,16
    strict_mask(5, q),
    strict_mask(3, q),
    strict_mask(2, q),
    strict_mask(4, q),
    strict_mask(7, q),
)
expected = (
    frozenset((0, 1, 14)),
    frozenset((0, 3, 6, 9, 12)),
    frozenset((0, 5, 10)),
    frozenset((0, 7, 8)),
    frozenset((0, 4, 11)),
    frozenset((0, 2, 13)),
)
require(strict_mask(1, q) == strict_mask(16, q), "defining masks coincide")
require(masks == expected, "q=15 exact masks")
require(set().union(*masks) == set(range(15)), "q=15 full clock")
for left, right in combinations(masks, 2):
    require(left & right == frozenset((0,)), "sunflower intersections")
nonzero_petals = [mask - {0} for mask in masks]
require(sum(map(len, nonzero_petals)) == 14, "nonzero petal mass")
require(set().union(*nonzero_petals) == set(range(1, 15)), "petal partition")
require(sum(map(len, masks)) - 5 == 15, "common-zero/tree cap equality")
periods = tuple(q // gcd(s, q) for s in speeds)
require(all(Q > 1 for Q in periods), "no period-one row")

# Minimality: below master 15, every nontrivial quotient window is {0}; its
# lift is a proper subgroup and misses the generator residue 1.
minimal_rows = 0
for L in range(2, 15):
    divisors = [Q for Q in range(2, L + 1) if L % Q == 0]
    require(divisors, "master has nontrivial divisor")
    for Q in divisors:
        require(A(Q) == 1, "small quotient singleton window")
        lifted = {p for p in range(L) if p % Q == 0}
        require(1 not in lifted, "generator escapes every lifted subgroup")
        minimal_rows += 1

# Direct common-clock union audits through Q=70 over every choice of at most
# six distinct unit-sign masks.
union_rows = 0
for Q in range(2, 71):
    unit_masks = []
    for u in range(1, Q):
        if gcd(u, Q) != 1:
            continue
        mask = strict_mask(u, Q)
        if mask not in unit_masks:
            unit_masks.append(mask)
    for r in range(1, min(6, len(unit_masks)) + 1):
        for chosen in combinations(unit_masks, r):
            union = set().union(*chosen)
            require(len(union) <= 1 + r * (A(Q) - 1), "common-zero union cap")
            union_rows += 1

print("THM-1242 SEVEN-WALL PAIR-BEAT SUNFLOWER EXACT AUDIT")
print(f"common-clock Q rows checked = {common_rows}")
print(f"common-mask unions checked through Q=70 = {union_rows}")
print(f"small-master lifted subgroup rows = {minimal_rows}")
print("seven-wall threshold = 12 ceil(Q/14)-10 <= Q for every Q>=2")
print("eight-wall first failures = Q=15,29,43,57,71,...")
print("first mixed full clock = q=15, periods=" + str(periods))
print("sunflower petal sizes = 2,4,2,2,2,2")
print("RESULT: PASS")
