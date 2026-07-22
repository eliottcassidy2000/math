#!/usr/bin/env python3
"""Exact referee for the triangle-packed pair bound in THM-2051.

The analytic Fejer/BV estimate is proved in THM-2051.  This script checks the
finite exceptional-ratio classification, all oriented two-exceptional-edge
triangles, the edge-to-triangle double count, and the final rational margin.

Tournament Analysis is not applicable.  The vertices here are speeds and the
observable is a symmetric covariance on unordered pairs.  Orienting it would
discard the signed rational value used by the triangle inequality; the useful
three-body object is a weighted triangle, not a tournament cycle.
"""

from fractions import Fraction as F
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fold14(n):
    residue = n % 14
    return residue * (14 - residue)


def pair_delta(a, b):
    """THM-965 centered covariance, reduced internally."""
    require(a > 0 and b > 0 and a != b, "speeds must be positive and distinct")
    divisor = gcd(a, b)
    a //= divisor
    b //= divisor
    if a > b:
        a, b = b, a
    return F(fold14(a + b) - fold14(b - a), 196 * a * b)


def ratio_delta(ratio):
    """Pair covariance for speeds in the positive rational ratio `ratio`."""
    require(ratio > 0 and ratio != 1, "ratio must be positive and nontrivial")
    if ratio < 1:
        ratio = 1 / ratio
    return pair_delta(ratio.denominator, ratio.numerator)


# Since fold14 lies in [0,49], delta(a,b)>=-1/(4ab) for a,b coprime.
# Thus delta<-1/220 is possible only when ab<=54.
threshold = F(-1, 220)
exceptional = [
    (a, b, pair_delta(a, b))
    for a in range(1, 55)
    for b in range(a + 1, 55)
    if a * b <= 54 and gcd(a, b) == 1 and pair_delta(a, b) < threshold
]
expected_exceptional = [
    (1, 10, F(-3, 490)),
    (1, 11, F(-4, 539)),
    (1, 12, F(-5, 588)),
    (1, 13, F(-6, 637)),
    (2, 11, F(-4, 539)),
    (3, 10, F(-3, 490)),
    (3, 11, F(-4, 539)),
]
require(exceptional == expected_exceptional, "exceptional ratio list failed")
require(F(-1, 4 * 55) == threshold, "large-product threshold failed")

# If a triangle has at most one exceptional edge, its covariance sum is at
# least -6/637-2/220.  If it has at least two, those edges share a vertex.
# Normalize that vertex to speed 1 and enumerate both orientations of the
# seven possible adjacent exceptional ratios.
triangle_cost = F(6, 637) + 2 * F(1, 220)
require(triangle_cost == F(1297, 70070), "coarse triangle cost failed")

exceptional_ratios = [F(13), F(12), F(11), F(11, 2), F(11, 3), F(10), F(10, 3)]
oriented_ratios = sorted(set(exceptional_ratios + [1 / r for r in exceptional_ratios]))
require(len(oriented_ratios) == 14, "unexpected oriented ratio count")

triangle_rows = []
for x in oriented_ratios:
    for y in oriented_ratios:
        if x == y:  # distinct speeds are required
            continue
        terms = (ratio_delta(x), ratio_delta(y), ratio_delta(x / y))
        triangle_rows.append((sum(terms, F(0)), x, y, *terms))

require(len(triangle_rows) == 182, "unexpected oriented triangle count")
exceptional_triangle_min = min(row[0] for row in triangle_rows)
require(exceptional_triangle_min == F(-150, 8281), "exceptional triangle minimum failed")
minimizers = [(row[1], row[2]) for row in triangle_rows if row[0] == exceptional_triangle_min]
require(minimizers == [(F(1, 13), F(13)), (F(13), F(1, 13))],
        "exceptional triangle equality classifier failed")
triangle_gap = triangle_cost + exceptional_triangle_min
require(triangle_gap == F(361, 910910) and triangle_gap > 0,
        "universal triangle comparison failed")

# In K_13 each edge occurs in eleven of the 286 triangles.
edge_sum_floor = -F(comb(13, 3), 11) * triangle_cost
require(edge_sum_floor == F(-1297, 2695), "triangle-to-edge double count failed")

base = F(2052, 16807)
pair_coefficient = F(24, 343)
pair_reserve = base + pair_coefficient * edge_sum_floor
require(pair_reserve == F(1668, 18865), "triangle-packed pair reserve failed")

K_high = F(10316592, 2401)
N = 2**20 + 1
fejer_tail = K_high * F(43, 2 * N)
require(fejer_tail == F(221806728, 2517633377), "Fejer tail failed")
final_margin = pair_reserve - fejer_tail
require(final_margin == F(43815012, 138469835735) and final_margin > 0,
        "final positive margin failed")

print("THM-2051 FEJER--BV TRIANGLE-PAIR REFEREE")
print(f"exceptional cutoff={threshold} list={exceptional}")
print(f"oriented exceptional triangles={len(triangle_rows)}")
print(f"exceptional triangle minimum={exceptional_triangle_min} minimizers={minimizers}")
print(f"universal triangle floor={-triangle_cost} comparison gap={triangle_gap}")
print(f"global pair-sum floor={edge_sum_floor} pair reserve={pair_reserve}")
print(f"H={2**20} higher-support tail<{fejer_tail} final margin>{final_margin}")
print("TOURNAMENT ANALYSIS=not applicable: symmetric weighted triangles are indispensable")
print("RESULT=PASS")
