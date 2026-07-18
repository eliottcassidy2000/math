#!/usr/bin/env python3
"""Independent Fraction replay for THM-1094's frozen finite-bank facts.

Unlike the C++ certificate, core-safe sets are built by classifying the cells
of the complete rational breakpoint arrangement.  The script independently
checks the 66-core minimum, the finite-bank cardinality, and the exact global
extremal reported by the exhaustive C++ scan.
"""

from fractions import Fraction as F
from itertools import combinations


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def safe_set(core: tuple[int, ...]) -> list[tuple[F, F]]:
    breakpoints = {F(0), F(1)}
    for p in core:
        for j in range(p + 1):
            for sign in (-1, 1):
                point = F(j, p) + sign * F(1, 14 * p)
                if 0 <= point <= 1:
                    breakpoints.add(point)
    ordered = sorted(breakpoints)
    result: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(circle_distance(p * midpoint) >= F(1, 14) for p in core):
            if result and result[-1][1] == left:
                result[-1] = (result[-1][0], right)
            else:
                result.append((left, right))
    return result


def remove_bad(region: list[tuple[F, F]], speed: int) -> list[tuple[F, F]]:
    result: list[tuple[F, F]] = []
    radius = F(1, 14 * speed)
    for left, right in region:
        cursor = left
        j_lo = int(left * speed) - 1
        j_hi = int(right * speed) + 1
        for j in range(j_lo, j_hi + 1):
            bad_left = F(j, speed) - radius
            bad_right = F(j, speed) + radius
            if bad_right <= left or right <= bad_left:
                continue
            bad_left = max(left, bad_left)
            bad_right = min(right, bad_right)
            if cursor < bad_left:
                result.append((cursor, bad_left))
            cursor = max(cursor, bad_right)
        if cursor < right:
            result.append((cursor, right))
    return result


def longest(region: list[tuple[F, F]]) -> F:
    return max(right - left for left, right in region)


cores = list(combinations(range(1, 13), 10))
core_rows = [(longest(safe_set(core)), core) for core in cores]
min_ell = min(row[0] for row in core_rows)
min_cores = [core for ell, core in core_rows if ell == min_ell]
p0 = (1, 2, 3, 4, 7, 8, 9, 10, 11, 12)
scaled_non_p0 = [
    (ell * (13 * max(core) + 1), core)
    for ell, core in core_rows
    if core != p0
]
second_scaled_min = min(row[0] for row in scaled_non_p0)
second_scaled_cores = [
    core for value, core in scaled_non_p0 if value == second_scaled_min
]

# Independently count exactly the guarded finite domain used by C++.
finite_pairs = 0
for ell, core in core_rows:
    lo = 13 * max(core) + 1
    x_max = int(F(209, 7) / ell) + 1
    for k1 in range(lo, x_max + 1):
        threshold = F(k1) * (49 * ell * k1 + 185) / (56 * ell * k1 - 24)
        y_max = int(threshold) + 1
        finite_pairs += max(0, y_max - k1)

extremal_core = (1, 2, 3, 5, 6, 7, 8, 9, 10, 11)
k1, k2 = 153, 159
remainder = remove_bad(remove_bad(safe_set(extremal_core), k1), k2)
extremal_length = longest(remainder)
extremal_intervals = [
    interval for interval in remainder if interval[1] - interval[0] == extremal_length
]
metric = 3 * k2 * extremal_length
ratio = 1 / metric

assert len(cores) == 66
assert min_ell == F(1, 112)
assert min_cores == [p0]
assert second_scaled_min == F(1727, 1008)
assert second_scaled_cores == [
    (1, 2, 3, 5, 7, 8, 9, 10, 11, 12),
    (1, 2, 3, 6, 7, 8, 9, 10, 11, 12),
]
assert finite_pairs == 9_246_070
assert extremal_length == F(158, 56_763)
assert metric == F(158, 119)
assert ratio == F(119, 158)
assert extremal_intervals == [
    (F(505, 2142), F(177, 742)),
    (F(565, 742), F(1637, 2142)),
]

print("THM-1094 independent Fraction replay")
print(f"cores={len(cores)}")
print(f"finite_pairs_from_cutoff={finite_pairs}")
print(f"min_core_ell={min_ell}")
for core in min_cores:
    print(f"min_ell_core={list(core)}")
print(f"min_scaled_non_p0={second_scaled_min}")
for core in second_scaled_cores:
    print(f"min_scaled_non_p0_core={list(core)}")
print(f"finite_bank_extremal_core={list(extremal_core)}")
print(f"finite_bank_extremal_pair=({k1},{k2})")
print(f"finite_bank_extremal_longest_component={extremal_length}")
print(f"finite_bank_extremal_3k2L={metric}")
print(f"finite_bank_extremal_R={ratio}")
for left, right in extremal_intervals:
    print(f"finite_bank_extremal_interval=[{left},{right}]")
print("replay=PASS")
