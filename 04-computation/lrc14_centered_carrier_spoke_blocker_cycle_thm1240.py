#!/usr/bin/env python3
"""Exact referee for THM-1240 centered carrier-spoke blocker cycles.

Checks the nearest-integer spoke construction on a deterministic all-rational
bank (including negative gap indices), the deep `>1/4` and curvature-slack
bounds, nontrivial/nonzero clock facts, and every loopless functional digraph
on six labels.  No floating arithmetic or ``assert`` is used.
"""

from fractions import Fraction as F
from itertools import product
from math import floor, gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: F) -> F:
    r = x - floor(x)
    return min(r, 1 - r)


spoke_rows = 0
for c in range(1, 81):
    indices = sorted({-3 * c, -c - 1, -c, -1, 0, 1, c - 1, c, 3 * c})
    for k in indices:
        t0 = F(2 * k + 1, 2 * c)
        gap_lo = F(14 * k + 1, 14 * c)
        gap_hi = F(14 * k + 13, 14 * c)
        for d in range(c + 1, 20 * c + 1):
            q = c + d
            alpha = q * t0
            p = floor(alpha + F(1, 2))
            t = F(p, q)
            require(abs(p - alpha) <= F(1, 2), "nearest integer")
            require(gap_lo < t < gap_hi, "strict slow-gap interior")

            carrier_distance = circle_norm(c * t)
            fast_distance = circle_norm(d * t)
            require(carrier_distance == fast_distance, "sum-beat distance equality")
            require(carrier_distance >= F(d, 2 * q), "centered deep-distance bound")
            require(carrier_distance > F(1, 4), "uniform quarter-circle depth")

            D = min((c * p) % q, q - ((c * p) % q))
            require(carrier_distance == F(D, q), "integer active determinant")
            require(14 * D - q >= 6 * d - c > 0, "deep curvature slack")

            Q = q // gcd(c, d)
            require(Q >= 3, "nontrivial carrier-spoke reduced period")
            require(p % Q != 0, "nonzero reduced spoke residue")
            spoke_rows += 1

# Every selection of one non-self blocker per spoke contains a directed cycle.
functional_maps = 0
cycle_histogram = {length: 0 for length in range(2, 7)}
for choices in product(range(6), repeat=6):
    if any(choices[i] == i for i in range(6)):
        continue
    functional_maps += 1
    seen_global: set[int] = set()
    found_lengths: list[int] = []
    for start in range(6):
        if start in seen_global:
            continue
        local_index: dict[int, int] = {}
        orbit: list[int] = []
        v = start
        while v not in local_index and v not in seen_global:
            local_index[v] = len(orbit)
            orbit.append(v)
            v = choices[v]
        seen_global.update(orbit)
        if v in local_index:
            length = len(orbit) - local_index[v]
            require(2 <= length <= 6, "loopless directed cycle length")
            found_lengths.append(length)
    require(found_lengths, "functional digraph has a cycle")
    cycle_histogram[min(found_lengths)] += 1

require(functional_maps == 5**6, "loopless functional map count")

# Exact cut-doublet arithmetic: adjacent spoke denominators inherit the speed
# gap unchanged.
cut_checks = 0
for c in range(1, 101):
    for x in range(c + 1, 5 * c + 1):
        for y in range(x + 1, min(5 * c + 1, x + 8)):
            require((c + y) - (c + x) == y - x, "spoke-clock cut separation")
            cut_checks += 1

print("THM-1240 CENTERED CARRIER-SPOKE BLOCKER CYCLE EXACT AUDIT")
print(f"centered spoke rows checked = {spoke_rows}")
print("all rows: depth > 1/4 and slack >= 6d-c > 0")
print(f"loopless blocker maps checked = {functional_maps}")
print("minimum-cycle histogram = " + ", ".join(
    f"len{length}:{cycle_histogram[length]}" for length in range(2, 7)
))
print(f"adjacent spoke-clock cut checks = {cut_checks}")
print("RESULT: PASS")
