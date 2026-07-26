#!/usr/bin/env python3
"""Exact companion for THM-2392.

The script uses only integer and rational arithmetic.  It checks the
thirteen-root toothpick combinatorics, DFT power sums through character
orthogonality, the blocker-cage ledgers, the 7-adic exact-overlap lemma,
all 150 strict profiles, and the final 124-ratio boundary bank.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fold14(value):
    residue = value % 14
    return residue * (14 - residue)


def defect(a, b):
    return fold14(a + b) - fold14(b - a)


def overlap(a, b):
    common = gcd(a, b)
    a //= common
    b //= common
    return F(1, 49) + F(defect(a, b), 196 * a * b)


def danger_intervals(speed):
    """Closed representatives of D_speed in [0,1]; endpoints have zero mass."""
    radius = F(1, 14 * speed)
    intervals = []
    for index in range(speed):
        center = F(index, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            intervals.append((F(0), right))
            intervals.append((1 + left, F(1)))
        elif right > 1:
            intervals.append((left, F(1)))
            intervals.append((F(0), right - 1))
        else:
            intervals.append((left, right))
    return sorted(intervals)


def interval_intersection_mass(first, second):
    total = F(0)
    for left_a, right_a in first:
        for left_b, right_b in second:
            left = max(left_a, left_b)
            right = min(right_a, right_b)
            if right > left:
                total += right - left
    return total


def gap_cap(gap):
    require(gap >= 1, "THM-2263 gap cap needs a positive gap")
    if gap % 2 == 0:
        return F(1, 49) + F(6, 49 * 13**gap)
    return F(1, 49) + F(5, 588 * 13**gap)


def character_sum_nonzero(exponent):
    """sum_{k=1}^{12} zeta^(k exponent) in Q."""
    return 12 if exponent % 13 == 0 else -1


print("THM-2392 CLEAN TOOTHPICK / CROSS-ANCESTOR CAGE -- exact companion")

# Same-line parent/child geometry.
same_line = interval_intersection_mass(
    danger_intervals(1), danger_intervals(13)
)
require(same_line == F(1, 91), "D_1 intersect D_13 mass is not 1/91")
require(overlap(1, 13) == same_line, "interval and folded overlap disagree")
print(f"same-line overlap rho(C,13C): {same_line}")

# A different 7-adic valuation means that exactly one reduced coefficient
# is divisible by seven.  Modulo fourteen, that coefficient equals its
# negative, so the folded defect vanishes.
seven_checks = 0
for a_residue in range(14):
    for b_residue in range(14):
        one_seven_multiple = (a_residue % 7 == 0) ^ (b_residue % 7 == 0)
        if one_seven_multiple:
            require(
                defect(a_residue, b_residue) == 0,
                f"7-adic folded cancellation failed at {a_residue},{b_residue}",
            )
            seven_checks += 1
require(seven_checks == 48, "wrong modulo-fourteen 7-adic test count")
print(f"different-nu7 exact law: {seven_checks} residue cells, rho=1/49")

# Clean-hole root profiles.  Labels are H,q1,...,q5.  The unique double
# pair is any two-subset.  An ordinary label outside it always exists.
labels = ("H", "q1", "q2", "q3", "q4", "q5")
ordinary = set(labels[1:])
pair_profiles = {}
selected_edges = {}
for pair in combinations(labels, 2):
    profile = {"H": 4, **{name: 2 for name in labels[1:]}}
    for name in pair:
        profile[name] -= 1
    outside = sorted(ordinary.difference(pair))
    require(len(outside) >= 3, f"too few ordinary labels outside {pair}")
    pair_profiles[pair] = tuple(profile[name] for name in labels)
    selected_edges[pair] = outside[0]

require(len(pair_profiles) == 15, "wrong double-pair type count")
guard_pairs = [pair for pair in pair_profiles if "H" in pair]
ordinary_pairs = [pair for pair in pair_profiles if "H" not in pair]
require(len(guard_pairs) == 5 and len(ordinary_pairs) == 10, "wrong pair split")
for pair in guard_pairs:
    require(
        sorted(pair_profiles[pair]) == [1, 2, 2, 2, 2, 3],
        f"wrong guard-pair singleton profile at {pair}",
    )
for pair in ordinary_pairs:
    require(
        sorted(pair_profiles[pair]) == [1, 1, 2, 2, 2, 4],
        f"wrong ordinary-pair singleton profile at {pair}",
    )

guard_word_count = 5 * factorial(12) // (factorial(3) * factorial(2) ** 4)
ordinary_word_count = (
    10 * factorial(12) // (factorial(4) * factorial(2) ** 3)
)
complete_words = 13 * (guard_word_count + ordinary_word_count)
require(complete_words == 648_648_000, "wrong complete root-word count")
labelled_cells = 2 * 15 * 13
exact_word_cells = 3 * 15 * 13
require(labelled_cells == 390 and exact_word_cells == 585, "wrong cell ledger")
print(
    "root words: 15 double pairs (5 guard, 10 ordinary), "
    f"{complete_words} complete words"
)
print(
    "deterministic charged cells: "
    f"owner/pair/edge={labelled_cells}, exact-word/pair/edge={exact_word_cells}"
)

# Two adjacent roots, after multiplication by a unit, are {0,1}.  Expand
# the second and fourth powers and sum each character over k != 0.
two_point = (0, 1)
square_sum_raw = 0
for left, right in product(two_point, repeat=2):
    square_sum_raw += character_sum_nonzero(right - left)

fourth_sum_raw = 0
for a, b, c, d in product(two_point, repeat=4):
    fourth_sum_raw += character_sum_nonzero(-a + b - c + d)

require(square_sum_raw == 22, "wrong unnormalized square sum")
require(fourth_sum_raw == 62, "wrong unnormalized fourth-power sum")
square_sum = F(square_sum_raw, 13**2)
fourth_sum = F(fourth_sum_raw, 13**4)
require(square_sum == F(22, 169), "wrong normalized square sum")
require(fourth_sum == F(62, 28561), "wrong normalized fourth-power sum")
print(
    "adjacent two-root DFT sums: "
    f"normalized ({square_sum}, {fourth_sum}); raw ({square_sum_raw}, {fourth_sum_raw})"
)

# General cage split after the two same-line and two nu_7-separated exact
# pieces.  One of the remaining two low cross pairs is large.
hole_floor = F(36, 343)
known_four = 2 * F(1, 91) + 2 * F(1, 49)
low_pair_baseline = (hole_floor - known_four) / 2
require(
    low_pair_baseline == F(1, 49) + F(3, 4459),
    "wrong low-cross baseline",
)
delta_threshold = 2 * (low_pair_baseline - F(1, 49))
require(delta_threshold == F(6, 4459), "wrong delta threshold")

theta = F(3, 4459)
cell_theta = theta / 390
product_theta = F(1, 2 * (delta_threshold - theta))
product_zero = F(1, 2 * delta_threshold)
require(cell_theta == F(1, 579670), "wrong theta cell floor")
require(743 < product_theta < 744, "wrong theta product threshold")
require(371 < product_zero < 372, "wrong zero-delta product threshold")
require(low_pair_baseline == F(94, 4459), "wrong zero-delta overlap")
print(
    "general dichotomy: delta<6/4459 -> "
    "ab<=1/[2(6/4459-delta)]"
)
print(
    "theta=3/4459: cell>=1/579670 or ab<=743; "
    "delta=0: rho>=94/4459 and ab<=371"
)

# Repeated-first profiles.  The high pairs are nu_7-separated, hence
# exactly 1/49 rather than merely receiving a thirteen-gap cap.
repeated_cage = 2 * F(1, 91) + 2 * F(23, 1092) + 2 * F(1, 49)
repeated_delta = hole_floor - repeated_cage
repeated_cell = repeated_delta / 390
require(repeated_cage == F(401, 3822), "wrong repeated-first cage cap")
require(repeated_delta == F(1, 26754), "wrong repeated-first clean floor")
require(repeated_cell == F(1, 10434060), "wrong repeated-first cell floor")
require(
    2 * repeated_cell / 13 == F(1, 67821390),
    "wrong repeated-first coefficient floor",
)
require(
    4 * repeated_cell / 169 == F(1, 440839035),
    "wrong repeated-first energy floor",
)
require(
    repeated_cell * square_sum == F(11, 881678070),
    "wrong repeated-first summed square floor",
)
require(
    repeated_cell * fourth_sum == F(31, 149003593830),
    "wrong repeated-first summed fourth floor",
)
print(
    "repeated-first: cage<=401/3822, delta>=1/26754, "
    "charged cell>=1/10434060"
)

# All 150 strict profiles.  The high C3 pairs are again exact 1/49.
strict_rows = []
for c_depth in range(5, 20):
    for middle_depth in range(2, c_depth):
        if middle_depth == 2:
            cage_cap = None
            clean_floor = None
        else:
            cage_cap = (
                2 * F(1, 91)
                + gap_cap(middle_depth)
                + gap_cap(middle_depth - 2)
                + 2 * F(1, 49)
            )
            clean_floor = hole_floor - cage_cap
        strict_rows.append((middle_depth, c_depth, cage_cap, clean_floor))

require(len(strict_rows) == 150, "wrong strict-profile count")
passing = [row for row in strict_rows if row[3] is not None and row[3] > 0]
require(len(passing) == 135, "wrong uniformly positive strict count")
minimum_clean = min(row[3] for row in passing)
minimizers = [(row[0], row[1]) for row in passing if row[3] == minimum_clean]
require(minimum_clean == F(6042, 9796423), "wrong strict clean floor")
require(
    minimizers == [(4, c_depth) for c_depth in range(5, 20)],
    "wrong strict worst-profile tie family",
)
strict_cell = minimum_clean / 390
require(strict_cell == F(1007, 636767495), "wrong strict cell floor")

even_correction = F(6, 49) * (F(1, 13**4) + F(1, 13**2))
odd_correction = F(5, 588) * (F(1, 13**3) + F(1, 13))
require(
    even_correction - odd_correction == F(85, 1199562),
    "wrong parity-boundary comparison",
)
require(
    max(row[2] for row in passing) == F(146022, 1399489),
    "wrong maximum strict cage cap",
)
print(
    "strict scan: 135/150 uniformly clean; worst middle depth b=4 "
    "(15 tied c-values)"
)
print(
    "strict uniform floors: delta>=6042/9796423, "
    "charged cell>=1007/636767495"
)

# At b=2, different low nu_7 values close the remaining zero-gap pair.
b2_cage_separated = (
    2 * F(1, 91) + F(1, 49) + gap_cap(2) + 2 * F(1, 49)
)
b2_delta_separated = hole_floor - b2_cage_separated
b2_cell_separated = b2_delta_separated / 390
require(b2_cage_separated == F(864, 8281), "wrong separated b=2 cage")
require(b2_delta_separated == F(36, 57967), "wrong separated b=2 clean floor")
require(b2_cell_separated == F(6, 3767855), "wrong separated b=2 cell")
print(
    "b=2 with different low nu7: delta>=36/57967, "
    "charged cell>=6/3767855"
)

# Same-(nu_13,nu_7) b=2 boundary.  If delta=0, the unresolved pair must
# exceed a fixed folded-defect threshold.  Enumerate exactly after the
# coarse product cutoff.
other_five = 2 * F(1, 91) + gap_cap(2) + 2 * F(1, 49)
required_overlap = hole_floor - other_five
required_excess = required_overlap - F(1, 49)
required_defect_ratio = 196 * required_excess
require(other_five == F(695, 8281), "wrong five-piece cage cap")
require(required_overlap == F(1219, 57967), "wrong residual overlap threshold")
require(required_excess == F(36, 57967), "wrong residual excess")
require(required_defect_ratio == F(144, 1183), "wrong defect-ratio threshold")
coarse_product = F(1, 4 * required_excess)
require(402 < coarse_product < 403, "wrong residual coarse product cutoff")

ratio_bank = []
for a in range(1, 403):
    for b in range(a, 403):
        if a * b > 402:
            break
        if gcd(a, b) != 1 or gcd(a * b, 91) != 1:
            continue
        ratio = F(defect(a, b), a * b)
        if ratio >= required_defect_ratio:
            ratio_bank.append((a, b, a * b, ratio))

require(len(ratio_bank) == 124, "wrong residual ratio-bank size")
max_product = max(row[2] for row in ratio_bank)
max_small_coordinate = max(row[0] for row in ratio_bank)
max_large_coordinate = max(row[1] for row in ratio_bank)
require(max_product == 345, "wrong exact residual product cap")
require(max_small_coordinate == 11, "wrong small-coordinate cap")
require(max_large_coordinate == 197, "wrong large-coordinate cap")
require(
    [(a, b) for a, b, prod, _ in ratio_bank if prod == max_product]
    == [(3, 115)],
    "wrong product-maximizing residual ratio",
)
require(
    [(a, b) for a, b, _, _ in ratio_bank if b == max_large_coordinate]
    == [(1, 197)],
    "wrong coordinate-maximizing residual ratio",
)
print(
    "same-(nu13,nu7) b=2 zero-clean bank: "
    "124 unordered ratios; ab<=345, min<=11, max<=197"
)
print("bank witnesses: max product (3,115); max coordinate (1,197)")
print(
    "VERDICT: 15 repeated + 135 strict profiles force a positive toothpick; "
    "the last 15 reduce to the explicit b=2 ratio bank"
)
