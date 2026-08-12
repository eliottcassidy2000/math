#!/usr/bin/env python3
"""Exact arithmetic checks for the disconnected-low geometry note."""
from fractions import Fraction as F
from math import comb


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


shape_counts = {1: 1, 2: 8, 3: 94, 4: 1295, 5: 19389}

# Coefficient of x^6 in product_n (1-x^n)^(-c_n), by bounded multiset DP.
dp = [0] * 7
dp[0] = 1
for size, colours in shape_counts.items():
    nxt = [0] * 7
    for used in range(7):
        for multiplicity in range((6 - used) // size + 1):
            ways = comb(colours + multiplicity - 1, multiplicity)
            nxt[used + size * multiplicity] += dp[used] * ways
    dp = nxt

partition_rows = (
    ((5, 1), 19389),
    ((4, 2), 10360),
    ((4, 1, 1), 1295),
    ((3, 3), 4465),
    ((3, 2, 1), 752),
    ((3, 1, 1, 1), 94),
    ((2, 2, 2), 120),
    ((2, 2, 1, 1), 36),
    ((2, 1, 1, 1, 1), 8),
    ((1, 1, 1, 1, 1, 1), 1),
)
require(dp[6] == 36520, ("profile coefficient", dp))
require(sum(value for _, value in partition_rows) == 36520, partition_rows)


def multipartite_tree_count(parts):
    vertices = sum(parts)
    return vertices ** (len(parts) - 2) * __import__("math").prod(
        (vertices - part) ** (part - 1) for part in parts
    )


tree_counts = tuple((parts, multipartite_tree_count(parts)) for parts, _ in partition_rows)
expected_tree_counts = (1, 32, 48, 81, 216, 324, 384, 576, 864, 1296)
require(tuple(value for _, value in tree_counts) == expected_tree_counts, tree_counts)

debt_max = F(186636088362, 11773143757375)
homogeneous_gap = F(1, 21) - debt_max
require(homogeneous_gap == F(1121969414539, 35319431272125), homogeneous_gap)

far_ratio_floor = F(1, 49) - F(6 * 168, 49 * (8 * 168 - 14))
require(far_ratio_floor == F(23, 4655), far_ratio_floor)
require(far_ratio_floor - F(1, 294) == F(43, 27930), far_ratio_floor)
far_ratio_tree_gap = F(5, 294) - debt_max
require(far_ratio_tree_gap == F(570672686921, 494472037809750), far_ratio_tree_gap)

print("disconnected_profile_count", dp[6])
print("partition_rows", partition_rows)
print("multipartite_tree_counts", tree_counts)
print("homogeneous_gap", homogeneous_gap)
print("q_ge_8p_physical_floor", far_ratio_floor)
print("q_ge_8p_gap_over_1_294", far_ratio_floor - F(1, 294))
print("five_q_ge_8p_edges_minus_debt", far_ratio_tree_gap)
print("all_exact_checks=PASS")
