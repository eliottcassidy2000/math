#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2264.

The analytic input is the signed tensor minorant proved in THM-2085.
This script checks every proof-facing specialization and optimization used by
THM-2264.  It is deliberately dependency-free.
"""

from fractions import Fraction as F


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tensor_constant(alpha, heights):
    eps = [F(1, h + 1) for h in heights]
    uppers = [alpha + e for e in eps]
    budget = sum((2 * e / u for e, u in zip(eps, uppers)), F(0))
    product = F(1)
    for u in uppers:
        product *= u
    return product * (1 - budget), budget


def lrc14_budget(heights):
    return sum((F(14, 6 * h + 13) for h in heights), F(0))


alpha14 = F(6, 7)

# General equal-degree identity:
# k=n-1, H=2n+1, epsilon=1/(2n+2), and
# B=(alpha+epsilon)^(n-2)*(alpha-(2n-3)epsilon).
for n in range(5, 201):
    h = 2 * n + 1
    alpha = F(3, 5)  # an arbitrary exact control value
    direct, _ = tensor_constant(alpha, [h] * (n - 1))
    eps = F(1, 2 * n + 2)
    closed = (alpha + eps) ** (n - 2) * (
        alpha - (2 * n - 3) * eps
    )
    require(direct == closed, f"equal-degree identity failed at n={n}")

# LRC(14), H=29.
b29, budget29 = tensor_constant(alpha14, [29] * 13)
expected_b29 = F(187, 210) ** 12 * F(1, 42)
require(b29 == expected_b29, "H=29 constant mismatch")
require(budget29 == F(182, 187), "H=29 budget mismatch")
require(b29 > 0, "H=29 constant is not positive")

# H=28 is the first failed equal-degree certificate.
eps28 = F(1, 29)
last_factor28 = alpha14 - 25 * eps28
require(last_factor28 == -F(1, 203), "H=28 factor mismatch")
require(last_factor28 < 0, "H=28 should fail")

# The useful anisotropic profiles from the theorem.
profile_rows = [
    (1, 21, F(23, 25993)),
    (2, 25, F(13, 2771)),
    (3, 26, F(89, 31603)),
    (5, 27, F(1, 935)),
    (10, 28, F(65, 33847)),
]
for count, special_height, expected_margin in profile_rows:
    heights = [special_height] * count + [29] * (13 - count)
    budget = lrc14_budget(heights)
    require(1 - budget == expected_margin, "profile margin mismatch")
    require(budget < 1, "profile should pass")

# Convex balancing gives the optimal integer profile at fixed total height.
# These two adjacent balanced profiles certify the sharp total budget.
budget366 = lrc14_budget([28] * 11 + [29] * 2)
budget367 = lrc14_budget([28] * 10 + [29] * 3)
require(budget366 == F(33866, 33847), "budget 366 mismatch")
require(budget367 == F(33782, 33847), "budget 367 mismatch")
require(budget366 > 1 and budget367 < 1, "total-budget boundary mismatch")

# One-coordinate height 1, all other coordinates height 105.
rank_budget = lrc14_budget([1] + [105] * 12)
require(rank_budget == F(12194, 12217), "height-105 budget mismatch")
require(1 - rank_budget == F(23, 12217), "height-105 margin mismatch")

print("THM-2264 exact referee")
print("equal_degree_identity=PASS n=5..200")
print(f"lrc14_H29_budget={budget29}")
print(f"lrc14_H29_B={b29}")
print(f"lrc14_H28_last_factor={last_factor28}")
for count, special_height, expected_margin in profile_rows:
    print(
        "anisotropic_profile="
        f"{count}x{special_height}+(13-{count})x29 "
        f"margin={expected_margin}"
    )
print(f"balanced_total_366_budget={budget366}")
print(f"balanced_total_367_budget={budget367}")
print(f"height_1_105_budget={rank_budget}")
print(f"height_1_105_margin={1-rank_budget}")
print("ALL_CHECKS_PASS")
