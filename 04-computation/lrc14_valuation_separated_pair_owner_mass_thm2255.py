#!/usr/bin/env python3
"""Exact companion for THM-2255.

The load-bearing finite bank is deliberately tiny.  THM-1166 gives

  rho(a,b) = 1/49 + (F(a+b)-F(b-a))/(196ab)

for a coprime reduced pair a<b, where 0<=F<=49.  If exactly one of a,b is
divisible by 13, the analytic tail ab>=346 is below 25/1183.  This script
exhausts the complementary bank ab<=345, freezes the equality case, and
checks the two exclusive-owner deductions used in the theorem.
"""

from fractions import Fraction as Q
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fold14(value: int) -> int:
    residue = value % 14
    return residue * (14 - residue)


def overlap(a: int, b: int) -> Q:
    common = gcd(a, b)
    a //= common
    b //= common
    if a > b:
        a, b = b, a
    return Q(4 * a * b + fold14(a + b) - fold14(b - a), 196 * a * b)


CAP = Q(25, 1183)
TAIL_START = 346
TAIL_BOUND = Q(1, 49) + Q(1, 4 * TAIL_START)
DELTA = Q(961, 6930)
GENERIC_DISTINCT_CAP = Q(1, 14)

require(TAIL_BOUND < CAP, "analytic tail does not lie below the claimed cap")
require(overlap(1, 169) == CAP, "the ratio 1:169 is not an equality case")
require(overlap(1, 2) == GENERIC_DISTINCT_CAP, "generic pair cap control failed")

finite_bank: list[tuple[int, int, Q]] = []
for a in range(1, TAIL_START):
    for b in range(a + 1, TAIL_START):
        if a * b > TAIL_START - 1:
            break
        if gcd(a, b) != 1:
            continue
        if (a % 13 == 0) == (b % 13 == 0):
            continue
        finite_bank.append((a, b, overlap(a, b)))

require(len(finite_bank) == 70, "unexpected finite-bank size")
finite_max = max(value for _, _, value in finite_bank)
equalities = [(a, b) for a, b, value in finite_bank if value == finite_max]
require(finite_max == CAP, "finite bank exceeds or misses the claimed cap")
require(equalities == [(1, 169)], "unexpected equality class")
require(all(value <= CAP for _, _, value in finite_bank), "finite cap failure")

# Three distinct valuations: all three blocker pairs use the sharp ramified cap.
strict_unique_floor = DELTA - 3 * CAP
strict_labelled_floor = strict_unique_floor / 3
require(strict_unique_floor == Q(88159, 1171170), "strict owner floor mismatch")
require(strict_labelled_floor == Q(88159, 3513510), "strict labelled floor mismatch")
require(strict_unique_floor > 0, "strict unique-owner floor is not positive")

# Depth pattern (1,1,c): one same-depth pair uses 1/14 and two ramified pairs
# use 25/1183.
double_unique_floor = DELTA - GENERIC_DISTINCT_CAP - 2 * CAP
double_labelled_floor = double_unique_floor / 3
require(double_unique_floor == Q(14627, 585585), "double owner floor mismatch")
require(double_labelled_floor == Q(14627, 1756755), "double labelled floor mismatch")
require(double_unique_floor > 0, "double unique-owner floor is not positive")

# Hostile equality survives arbitrary common scaling.
for scale in (1, 2, 7, 13, 91, 1000):
    require(overlap(scale, 169 * scale) == CAP, "scaled equality control failed")

print("THM-2255 exact folded-formula companion")
print(f"analytic tail starts at ab={TAIL_START}")
print(f"tail bound at ab={TAIL_START}: {TAIL_BOUND} = {float(TAIL_BOUND):.15f}")
print(f"sharp ramified pair cap: {CAP} = {float(CAP):.15f}")
print(f"finite coprime bank size: {len(finite_bank)}")
print(f"finite maximum: {finite_max}")
print(f"reduced equality classes: {equalities}")
print(f"hostile equality rho(1,169): {overlap(1, 169)}")
print(f"strict-depth total exclusive-owner floor: {strict_unique_floor}")
print(f"strict-depth labelled pigeonhole floor: {strict_labelled_floor}")
print(f"double-depth total exclusive-owner floor: {double_unique_floor}")
print(f"double-depth labelled pigeonhole floor: {double_labelled_floor}")
print("all exact checks passed")
