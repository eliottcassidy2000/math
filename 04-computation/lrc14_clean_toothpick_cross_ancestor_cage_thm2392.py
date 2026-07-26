#!/usr/bin/env python3
"""Exact companion for THM-2392.

The script uses only integer and rational arithmetic.  It checks the
thirteen-root toothpick combinatorics, DFT power sums through character
orthogonality, the blocker-cage ledgers, the 7-adic exact-overlap lemma,
all 150 strict profiles, and the compatible ten-orientation boundary bank.
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


def interval_union(*families):
    merged = []
    for left, right in sorted(interval for family in families for interval in family):
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return [(left, right) for left, right in merged]


def danger_intervals_scaled(speed, denominator):
    require(denominator % (14 * speed) == 0, "incompatible exact grid")
    radius = denominator // (14 * speed)
    step = denominator // speed
    intervals = []
    for index in range(speed):
        center = index * step
        left = center - radius
        right = center + radius
        if left < 0:
            intervals.append((0, right))
            intervals.append((denominator + left, denominator))
        elif right > denominator:
            intervals.append((left, denominator))
            intervals.append((0, right - denominator))
        else:
            intervals.append((left, right))
    return sorted(intervals)


def integer_interval_union(*families):
    merged = []
    for left, right in sorted(interval for family in families for interval in family):
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return merged


def integer_intersection_length(first, second):
    i = 0
    j = 0
    total = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if right > left:
            total += right - left
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def compatible_low_union_mass(a, b):
    denominator = 14 * 169 * a * b
    quotient_union = integer_interval_union(
        danger_intervals_scaled(b, denominator),
        danger_intervals_scaled(13 * a, denominator),
    )
    original_union = integer_interval_union(
        danger_intervals_scaled(13 * b, denominator),
        danger_intervals_scaled(169 * a, denominator),
    )
    return F(integer_intersection_length(quotient_union, original_union), denominator)


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
for pair in combinations(labels, 2):
    profile = {"H": 4, **{name: 2 for name in labels[1:]}}
    for name in pair:
        profile[name] -= 1
    outside = sorted(ordinary.difference(pair))
    require(len(outside) >= 3, f"too few ordinary labels outside {pair}")
    pair_profiles[pair] = tuple(profile[name] for name in labels)

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
q_star = "q1"
q_star_inside = [pair for pair in pair_profiles if q_star in pair]
q_star_outside = [pair for pair in pair_profiles if q_star not in pair]
require(
    len(q_star_inside) == 5 and len(q_star_outside) == 10,
    "wrong q-star status split",
)
labelled_cells = 2 * 2 * 13
exact_word_cells = 3 * 2 * 13
require(labelled_cells == 52 and exact_word_cells == 78, "wrong cell ledger")
print(
    "root words: 15 double pairs (5 guard, 10 ordinary), "
    f"{complete_words} complete words"
)
print(
    "top-labelled charged cells: "
    f"owner/status/support={labelled_cells}, exact-word/status/support={exact_word_cells}"
)

# A singleton is {0}; two adjacent roots, after multiplication by a unit,
# are {0,1}.  Expand the second and fourth powers and sum each character.
one_point = (0,)
two_point = (0, 1)
singleton_square_raw = 0
for left, right in product(one_point, repeat=2):
    singleton_square_raw += character_sum_nonzero(right - left)
singleton_fourth_raw = 0
for a, b, c, d in product(one_point, repeat=4):
    singleton_fourth_raw += character_sum_nonzero(-a + b - c + d)
square_sum_raw = 0
for left, right in product(two_point, repeat=2):
    square_sum_raw += character_sum_nonzero(right - left)

fourth_sum_raw = 0
for a, b, c, d in product(two_point, repeat=4):
    fourth_sum_raw += character_sum_nonzero(-a + b - c + d)

require(square_sum_raw == 22, "wrong unnormalized square sum")
require(fourth_sum_raw == 62, "wrong unnormalized fourth-power sum")
require(singleton_square_raw == 12, "wrong singleton square sum")
require(singleton_fourth_raw == 12, "wrong singleton fourth-power sum")
singleton_square_sum = F(singleton_square_raw, 13**2)
singleton_fourth_sum = F(singleton_fourth_raw, 13**4)
square_sum = F(square_sum_raw, 13**2)
fourth_sum = F(fourth_sum_raw, 13**4)
require(singleton_square_sum == F(12, 169), "wrong normalized singleton square")
require(
    singleton_fourth_sum == F(12, 28561),
    "wrong normalized singleton fourth power",
)
require(square_sum == F(22, 169), "wrong normalized square sum")
require(fourth_sum == F(62, 28561), "wrong normalized fourth-power sum")
print(
    "top-labelled DFT sums singleton/adjacent: "
    f"({singleton_square_sum},{singleton_fourth_sum}) / ({square_sum},{fourth_sum})"
)

# THM-2391 supplies a unique excess among seven siblings.  Its location
# resolves simultaneous ownership at d=0 and two exclusive owners at
# each nonzero d.  Cross with the 26 singleton/edge supports.
septimal_owner_categories = 1 + 6 * 2
target_support_categories = 2 * 13
joint_tensor_cells = septimal_owner_categories * target_support_categories
require(septimal_owner_categories == 13, "wrong septimal owner-category count")
require(target_support_categories == 26, "wrong target support count")
require(joint_tensor_cells == 338, "wrong C7 x C13 cell ledger")
for d in range(7):
    # The normalized DFT of 1_{d} is one seventh times a root of unity,
    # hence every one of its seven colours has exact squared magnitude 1/49.
    septimal_magnitude_squared = F(1, 49)
    require(
        septimal_magnitude_squared > 0,
        f"septimal singleton lost a Fourier colour at d={d}",
    )
print(
    "same-parent tensor cells: 13 owner-resolved septimal categories "
    "* 26 target supports = 338"
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
cell_theta = theta / 52
product_theta = F(1, 2 * (delta_threshold - theta))
product_zero = F(1, 2 * delta_threshold)
require(cell_theta == F(3, 231868), "wrong theta cell floor")
require(743 < product_theta < 744, "wrong theta product threshold")
require(371 < product_zero < 372, "wrong zero-delta product threshold")
require(low_pair_baseline == F(94, 4459), "wrong zero-delta overlap")
print(
    "general dichotomy: delta<6/4459 -> "
    "ab<=1/[2(6/4459-delta)]"
)
print(
    "theta=3/4459: cell>=3/231868 or ab<=743; "
    "delta=0: rho>=94/4459 and ab<=371"
)

# Repeated-first profiles.  The high pairs are nu_7-separated, hence
# exactly 1/49 rather than merely receiving a thirteen-gap cap.
repeated_cage = 2 * F(1, 91) + 2 * F(23, 1092) + 2 * F(1, 49)
repeated_delta = hole_floor - repeated_cage
repeated_cell = repeated_delta / 52
repeated_tensor_cell = repeated_delta / joint_tensor_cells
require(repeated_cage == F(401, 3822), "wrong repeated-first cage cap")
require(repeated_delta == F(1, 26754), "wrong repeated-first clean floor")
require(repeated_cell == F(1, 1391208), "wrong repeated-first cell floor")
require(
    repeated_tensor_cell == F(1, 9042852),
    "wrong repeated-first tensor-cell floor",
)
require(
    2 * repeated_cell / 13 == F(1, 9042852),
    "wrong repeated-first coefficient floor",
)
require(
    4 * repeated_cell / 169 == F(1, 58778538),
    "wrong repeated-first energy floor",
)
require(
    repeated_cell * singleton_square_sum == F(1, 19592846),
    "wrong repeated-first uniform summed square floor",
)
require(
    repeated_cell * singleton_fourth_sum == F(1, 3311190974),
    "wrong repeated-first uniform summed fourth floor",
)
print(
    "repeated-first: cage<=401/3822, delta>=1/26754, "
    "top cell>=1/1391208, tensor cell>=1/9042852"
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
strict_cell = minimum_clean / 52
strict_tensor_cell = minimum_clean / joint_tensor_cells
require(strict_cell == F(3021, 254706998), "wrong strict cell floor")
require(
    strict_tensor_cell == F(3021, 1655595487),
    "wrong strict tensor-cell floor",
)

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
    "top cell>=3021/254706998, tensor cell>=3021/1655595487"
)

# At b=2, different low nu_7 values close the remaining zero-gap pair.
b2_cage_separated = (
    2 * F(1, 91) + F(1, 49) + gap_cap(2) + 2 * F(1, 49)
)
b2_delta_separated = hole_floor - b2_cage_separated
b2_cell_separated = b2_delta_separated / 52
b2_tensor_separated = b2_delta_separated / joint_tensor_cells
require(b2_cage_separated == F(864, 8281), "wrong separated b=2 cage")
require(b2_delta_separated == F(36, 57967), "wrong separated b=2 clean floor")
require(b2_cell_separated == F(9, 753571), "wrong separated b=2 cell")
require(
    b2_tensor_separated == F(18, 9796423),
    "wrong separated b=2 tensor cell",
)
print(
    "b=2 with different low nu7: delta>=36/57967, "
    "top cell>=9/753571, tensor cell>=18/9796423"
)

# Same-(nu_13,nu_7) b=2 boundary.  Orient
# (C2,c1)=13h(a,b), so (C1,C2,c1,c2)=h(b,13a,13b,169a).
# The two cross-pair defects agree because 169 == 1 modulo 14.  Their
# denominators differ by 169.
paired_base = 2 * F(1, 91) + 4 * F(1, 49)
paired_excess = hole_floor - paired_base
paired_defect_threshold = paired_excess * 196 * F(169, 170)
require(paired_base == F(66, 637), "wrong paired blocker base invoice")
require(paired_excess == F(6, 4459), "wrong paired excess invoice")
require(
    paired_defect_threshold == F(156, 595),
    "wrong paired defect-ratio threshold",
)
coarse_product = F(49, 1) / paired_defect_threshold
require(186 < coarse_product < 187, "wrong paired coarse product cutoff")

preliminary_bank = []
for a in range(1, 187):
    for b in range(a, 187):
        if a * b > 186:
            break
        if gcd(a, b) != 1 or gcd(a * b, 91) != 1:
            continue
        ratio = F(defect(a, b), a * b)
        if ratio >= paired_defect_threshold:
            preliminary_bank.append((a, b, a * b, ratio))

require(len(preliminary_bank) == 52, "wrong preliminary paired-bank size")
require(max(row[2] for row in preliminary_bank) == 177, "wrong paired product cap")
require(max(row[0] for row in preliminary_bank) == 10, "wrong paired small cap")
require(max(row[1] for row in preliminary_bank) == 85, "wrong paired large cap")
small_census = {}
for a, _, _, _ in preliminary_bank:
    small_census[a] = small_census.get(a, 0) + 1
require(
    small_census == {1: 24, 2: 9, 3: 8, 4: 4, 5: 4, 8: 1, 9: 1, 10: 1},
    "wrong preliminary-bank small-coordinate census",
)

# Use the complete compatible low union, not a sum of separately sharp
# pair caps.  No-clean absorption needs mu(U_ab) >= 22/343.
compatible_threshold = hole_floor - 2 * F(1, 49)
require(compatible_threshold == F(22, 343), "wrong compatible-union threshold")
compatible_bank = []
for left, right, _, _ in preliminary_bank:
    orientations = [(left, right)] if left == right else [(left, right), (right, left)]
    for a, b in orientations:
        mass = compatible_low_union_mass(a, b)
        if mass >= compatible_threshold:
            compatible_bank.append((a, b, mass))

expected_compatible = [
    (1, 1, F(193, 1183)),
    (1, 2, F(114, 1183)),
    (2, 1, F(239, 2366)),
    (1, 3, F(263, 3549)),
    (3, 1, F(95, 1183)),
    (4, 1, F(331, 4732)),
    (2, 3, F(43, 546)),
    (3, 2, F(95, 1183)),
    (3, 4, F(491, 7098)),
    (4, 3, F(331, 4732)),
]
require(
    compatible_bank == expected_compatible,
    f"wrong oriented compatible bank: {compatible_bank}",
)
unordered_compatible = sorted(
    {tuple(sorted((a, b))) for a, b, _ in compatible_bank}
)
require(
    unordered_compatible == [(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)],
    "wrong six-ratio compatible bank",
)
print(
    "same-(nu13,nu7) b=2 preliminary bank: "
    "52 unordered; ab<=177, min<=10, max<=85"
)
print(
    "compatible low-union bank: 10 oriented / 6 unordered ratios "
    "{1:1,1:2,1:3,1:4,2:3,3:4}"
)

# THM-2391: at M>=2 both blocker steps have the same absolute value,
# hence C2/C1=13a/b is +/-1 modulo 7^M.
def address_survivors(modulus):
    return [
        (a, b)
        for a, b, _ in compatible_bank
        if (13 * a - b) % modulus == 0 or (13 * a + b) % modulus == 0
    ]


mod49 = address_survivors(49)
mod343 = address_survivors(343)
require(mod49 == [(4, 3)], f"wrong M>=2 address survivor: {mod49}")
require(mod343 == [], f"wrong M>=3 address survivor: {mod343}")
print("THM-2391 address filter: mod49 -> {(4,3)}; mod343 -> empty")

expected_digits = [
    (1, 1, -1, 2),
    (1, 2, 3, 1),
    (2, 1, -2, 4),
    (1, 3, 2, 1),
    (3, 1, -3, 6),
    (4, 1, 3, 7),
    (2, 3, -3, 5),
    (3, 2, 2, 5),
    (3, 4, 1, 5),
    (4, 3, 1, 7),
]
digit_rows = []
for a, b, _ in compatible_bank:
    balanced = [sigma for sigma in (-3, -2, -1, 1, 2, 3) if (13 * a - sigma * b) % 7 == 0]
    require(len(balanced) == 1, f"nonunique balanced digit for {(a, b)}")
    sigma = balanced[0]
    d = (13 * a - sigma * b) // 7
    digit_rows.append((a, b, sigma, d))
require(digit_rows == expected_digits, f"wrong first-digit address table: {digit_rows}")
require({row[3] for row in digit_rows} == {1, 2, 4, 5, 6, 7}, "wrong descendant set")
require(
    [(a, b) for a, b, _, d in digit_rows if d == a] == [(1, 2), (1, 3)],
    "wrong exact-C2 descendant rows",
)
require(
    [(a, b) for a, b, _, d in digit_rows if d == 7] == [(4, 1), (4, 3)],
    "wrong septimal-climb rows",
)
print("balanced first-digit descendants d: {1,2,4,5,6,7}; exact C2 rows 2; climbs 2")
print(
    "VERDICT: 15 repeated + 135 strict profiles force a positive toothpick; "
    "b=2 no-clean has 10 orientations at M=1, one at M=2, none at M>=3"
)
