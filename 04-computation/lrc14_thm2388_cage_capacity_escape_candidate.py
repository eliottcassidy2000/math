#!/usr/bin/env python3
"""Exact referee for the THM-2388 blocker-cage escape candidate.

THM-2388 and THM-2263 are proved and independently audited.  This script
verifies the arithmetic implication from THM-2388's hole-mass/cage
inclusions and THM-2263's pair caps.
"""

from fractions import Fraction
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HOLE_MASS = Fraction(36, 343)
SAME_OWNER = Fraction(1, 91)
SHALLOW_GAP_ONE = Fraction(23, 1092)


def gap_cap(gap):
    """THM-2263 sharp upper cap at a positive 13-adic valuation gap."""
    require(gap >= 1, "gap cap called at a nonpositive gap")
    power = 13 ** gap
    if gap % 2 == 0:
        return Fraction(1, 49) + Fraction(6, 49 * power)
    return Fraction(1, 49) + Fraction(5, 588 * power)


def fold14(value):
    residue = value % 14
    return residue * (14 - residue)


def overlap(speed_a, speed_b):
    """THM-1166/2263 exact danger-pair overlap."""
    common = gcd(speed_a, speed_b)
    reduced_a = speed_a // common
    reduced_b = speed_b // common
    return Fraction(
        4 * reduced_a * reduced_b
        + fold14(reduced_a + reduced_b)
        - fold14(reduced_b - reduced_a),
        196 * reduced_a * reduced_b,
    )


require(gap_cap(1) == SHALLOW_GAP_ONE, "gap-one cap changed")
require(gap_cap(2) == Fraction(25, 1183), "gap-two cap changed")
require(gap_cap(3) == Fraction(3767, 184548), "gap-three cap changed")
require(gap_cap(4) == Fraction(583, 28561), "gap-four cap changed")

# Among gaps 3,...,18, the odd and even defects decrease separately.
# The even d=4 defect is larger than the odd d=3 defect by the exact
# ratio 72/65, so d=4 is the global maximum.
deep_gap_caps = tuple((gap, gap_cap(gap)) for gap in range(3, 19))
deep_worst_gap, deep_worst_cap = max(
    deep_gap_caps,
    key=lambda item: item[1],
)
require(
    (deep_worst_gap, deep_worst_cap) == (4, Fraction(583, 28561)),
    "deep-gap maximum changed",
)


def cage_bound(profile):
    """Six-pair union cap, or None if a cross pair has valuation gap zero."""
    blocker_valuations = profile
    quotient_valuations = tuple(value - 1 for value in profile)
    total = Fraction(0)
    terms = []
    for quotient_index in range(3):
        for original_index in range(2):
            if quotient_index == original_index:
                # D_(C_j) intersect D_(13 C_j) has exact mass 1/91.
                value = SAME_OWNER
                kind = "same-owner"
            else:
                gap = abs(
                    quotient_valuations[quotient_index]
                    - blocker_valuations[original_index]
                )
                if gap == 0:
                    return None, ()
                value = gap_cap(gap)
                kind = f"gap-{gap}"
            total += value
            terms.append((quotient_index + 1, original_index + 1, kind, value))
    return total, tuple(terms)


# All fifteen repeated-first profiles escape.  Their uniform worst case
# has deep gap four, i.e. c=6.
repeated_rows = []
for deepest in range(5, 20):
    profile = (1, 1, deepest)
    bound, terms = cage_bound(profile)
    require(bound is not None, "repeated-first profile gained a zero gap")
    require(bound < HOLE_MASS, "repeated-first cage reached the hole mass")
    repeated_rows.append((profile, bound, HOLE_MASS - bound, terms))

repeated_worst = min(repeated_rows, key=lambda row: row[2])
require(repeated_worst[0] == (1, 1, 6), "repeated worst profile changed")
require(
    repeated_worst[1] == Fraction(17981, 171366),
    "repeated cage cap changed",
)
REPEATED_MARGIN = repeated_worst[2]
require(
    REPEATED_MARGIN == Fraction(1693, 58778538),
    "repeated off-cage margin changed",
)

# Strict profiles: positive-gap pair caps already settle 117/150.
# Thirty profiles have a cross quotient/original pair at the same
# valuation (b=2 or c=b+1); three further rows miss this coarse union
# threshold despite all cross gaps being positive.
strict_pass = []
strict_zero_gap = []
strict_coarse_fail = []
for deepest in range(5, 20):
    for middle in range(2, deepest):
        profile = (1, middle, deepest)
        bound, terms = cage_bound(profile)
        if bound is None:
            strict_zero_gap.append(profile)
        elif bound < HOLE_MASS:
            strict_pass.append((profile, bound, HOLE_MASS - bound, terms))
        else:
            strict_coarse_fail.append((profile, bound, bound - HOLE_MASS, terms))

require(len(strict_pass) == 117, "strict passing profile count changed")
require(len(strict_zero_gap) == 30, "strict zero-gap profile count changed")
require(
    tuple(profile for profile, *_ in strict_coarse_fail)
    == ((1, 3, 6), (1, 4, 6), (1, 4, 7)),
    "strict positive-gap coarse failures changed",
)
require(
    set(strict_zero_gap)
    == (
        {(1, 2, deepest) for deepest in range(5, 20)}
        | {(1, middle, middle + 1) for middle in range(4, 19)}
    ),
    "strict zero-gap family changed",
)

strict_worst = min(strict_pass, key=lambda row: row[2])
require(strict_worst[0] == (1, 3, 5), "strict worst passing profile changed")
require(
    strict_worst[2] == Fraction(67, 2260713),
    "strict passing margin changed",
)

# Compatibility of the four cross pairs does not by itself settle the
# three positive-gap coarse failures.  These common unit-part ladders
# realize actual pair sums whose six-pair union invoices still exceed
# the hole floor.  They are hostiles to a proposed "simultaneous maxima
# are incompatible, therefore the threshold clears" shortcut; they say
# nothing about the actual cage union, where inclusion-exclusion or X
# can still save the argument.
compatibility_hostiles = {
    (1, 3, 6): (
        (1, 12, 12),
        Fraction(251891, 2399124),
        Fraction(4307, 117557076),
    ),
    (1, 4, 6): (
        (1, 1, 12),
        Fraction(20991, 199927),
        Fraction(363, 9796423),
    ),
    (1, 4, 7): (
        (1, 1, 1),
        Fraction(273066, 2599051),
        Fraction(13686, 127353499),
    ),
}
for profile, (unit_parts, expected_invoice, expected_excess) in (
    compatibility_hostiles.items()
):
    _, middle, deepest = profile
    unit_1, unit_2, unit_3 = unit_parts
    cross_terms = (
        overlap(unit_1, 13**middle * unit_2),
        overlap(13 ** (middle - 2) * unit_2, unit_1),
        overlap(13 ** (deepest - 2) * unit_3, unit_1),
        overlap(13 ** (deepest - middle - 1) * unit_3, unit_2),
    )
    invoice = 2 * SAME_OWNER + sum(cross_terms, Fraction())
    require(invoice == expected_invoice, f"compatible hostile changed at {profile}")
    require(
        invoice - HOLE_MASS == expected_excess,
        f"compatible hostile excess changed at {profile}",
    )

pair_only_escaped_profile_count = len(repeated_rows) + len(strict_pass)
require(pair_only_escaped_profile_count == 132, "pair-only profile count changed")
PAIR_ONLY_UNIFORM_MARGIN = min(
    REPEATED_MARGIN,
    strict_worst[2],
)
require(
    PAIR_ONLY_UNIFORM_MARGIN == REPEATED_MARGIN,
    "pair-only uniform margin changed",
)

# The live THM-2388 lane has additional septimal typing. If two speeds
# have distinct 7-adic valuations, divide by their gcd. One reduced
# coefficient is 0 or 7 modulo 14, and the other is a 7-unit. The folded
# numerator in THM-1166 is then identically zero, so the overlap is
# exactly 1/49.
SEPTIMALLY_SEPARATED = Fraction(1, 49)
for high_residue in (0, 7):
    for low_residue in range(14):
        if low_residue % 7 == 0:
            continue
        require(
            fold14(high_residue + low_residue)
            == fold14(low_residue - high_residue),
            "septimal folded cancellation changed",
        )

# In THM-2388, nu_7(C_3)=nu_7(c_3)>M while c_1,c_2 are below M.
# Hence both C_3/c_j cross pairs contribute exactly 1/49, not merely
# their 13-adic caps.
live_repeated_bound = (
    2 * SAME_OWNER
    + 2 * SHALLOW_GAP_ONE
    + 2 * SEPTIMALLY_SEPARATED
)
LIVE_REPEATED_MARGIN = HOLE_MASS - live_repeated_bound
require(live_repeated_bound == Fraction(401, 3822), "live repeated cap changed")
require(
    LIVE_REPEATED_MARGIN == Fraction(1, 26754),
    "live repeated margin changed",
)

live_strict_rows = []
same_septimal_residual = []
for deepest in range(5, 20):
    for middle in range(2, deepest):
        profile = (1, middle, deepest)
        if middle == 2:
            same_septimal_residual.append(profile)
            continue
        bound = (
            2 * SAME_OWNER
            + gap_cap(middle)
            + gap_cap(middle - 2)
            + 2 * SEPTIMALLY_SEPARATED
        )
        require(bound < HOLE_MASS, "live strict cage reached the hole mass")
        live_strict_rows.append((profile, bound, HOLE_MASS - bound))

require(len(live_strict_rows) == 135, "live strict profile count changed")
require(len(same_septimal_residual) == 15, "same-septimal residual count changed")
live_strict_margin = min(row[2] for row in live_strict_rows)
live_strict_worst = {
    profile for profile, _, margin in live_strict_rows
    if margin == live_strict_margin
}
require(
    live_strict_margin == Fraction(6042, 9796423),
    "live strict margin changed",
)
require(
    live_strict_worst
    == {(1, 4, deepest) for deepest in range(5, 20)},
    "live strict worst family changed",
)

# In the remaining b=2 family, C_2 and c_1 have the same 13-adic
# valuation. If their 7-adic valuations differ, the same exact septimal
# cancellation closes the pair invoice. Thus the genuine residue is
# b=2 together with nu_7(C_2)=nu_7(c_1).
b2_septimally_distinct_bound = (
    2 * SAME_OWNER
    + gap_cap(2)
    + 3 * SEPTIMALLY_SEPARATED
)
B2_SEPTIMALLY_DISTINCT_MARGIN = HOLE_MASS - b2_septimally_distinct_bound
require(
    b2_septimally_distinct_bound == Fraction(864, 8281),
    "b=2 septimally distinct cap changed",
)
require(
    B2_SEPTIMALLY_DISTINCT_MARGIN == Fraction(36, 57967),
    "b=2 septimally distinct margin changed",
)

# If the b=2 pair also has equal septimal valuation, a large enough
# actual overlap is necessary for the cage invoice to absorb the hole
# floor. Reduce (C_2,c_1)=g(a,b), a<=b. THM-1166 turns the threshold
# into a finite folded-defect bank.
b2_base_without_cross_ancestor = (
    2 * SAME_OWNER
    + gap_cap(2)
    + 2 * SEPTIMALLY_SEPARATED
)
CROSS_ANCESTOR_THRESHOLD = HOLE_MASS - b2_base_without_cross_ancestor
require(
    b2_base_without_cross_ancestor == Fraction(695, 8281),
    "b=2 base invoice changed",
)
require(
    CROSS_ANCESTOR_THRESHOLD == Fraction(1219, 57967),
    "cross-ancestor threshold changed",
)
require(
    196 * (CROSS_ANCESTOR_THRESHOLD - Fraction(1, 49))
    == Fraction(144, 1183),
    "folded cross-ancestor threshold changed",
)

# The folded numerator is at most 49. Hence every survivor has
# ab<=floor(49*1183/144)=402. Exact enumeration sharpens the bank.
CROSS_ANCESTOR_PRODUCT_TAIL = (49 * 1183) // 144
require(CROSS_ANCESTOR_PRODUCT_TAIL == 402, "cross-ancestor tail changed")
cross_ancestor_bank = []
for reduced_a in range(1, CROSS_ANCESTOR_PRODUCT_TAIL + 1):
    for reduced_b in range(
        reduced_a,
        CROSS_ANCESTOR_PRODUCT_TAIL // reduced_a + 1,
    ):
        if gcd(reduced_a, reduced_b) != 1:
            continue
        if reduced_a % 7 == 0 or reduced_b % 7 == 0:
            continue
        if reduced_a % 13 == 0 or reduced_b % 13 == 0:
            continue
        folded_defect = (
            fold14(reduced_a + reduced_b)
            - fold14(reduced_b - reduced_a)
        )
        folded_survives = (
            1183 * folded_defect
            >= 144 * reduced_a * reduced_b
        )
        overlap_survives = (
            overlap(reduced_a, reduced_b)
            >= CROSS_ANCESTOR_THRESHOLD
        )
        require(
            folded_survives == overlap_survives,
            "folded/overlap bank equivalence changed",
        )
        if folded_survives:
            cross_ancestor_bank.append(
                (reduced_a, reduced_b, folded_defect)
            )

require(len(cross_ancestor_bank) == 124, "cross-ancestor bank size changed")
require(
    max(a * b for a, b, _ in cross_ancestor_bank) == 345,
    "cross-ancestor maximum product changed",
)
require(
    max(a for a, _, _ in cross_ancestor_bank) == 11,
    "cross-ancestor minimum-coordinate ceiling changed",
)
require(
    max(b for _, b, _ in cross_ancestor_bank) == 197,
    "cross-ancestor maximum-coordinate ceiling changed",
)
cross_ancestor_counts = {
    reduced_a: sum(1 for a, _, _ in cross_ancestor_bank if a == reduced_a)
    for reduced_a in sorted({a for a, _, _ in cross_ancestor_bank})
}
require(
    cross_ancestor_counts
    == {1: 49, 2: 20, 3: 22, 4: 11, 5: 11,
        6: 1, 8: 2, 9: 3, 10: 2, 11: 3},
    "cross-ancestor first-coordinate census changed",
)

# The gap-two cross pair is algebraically coupled to the same ancestor
# ratio. With (C_2,c_1)=13h(a,b), the other pair is (C_1,c_2)=(b,169a)
# after common scaling. Because 169=1 mod 14, its folded defect is the
# same Delta(a,b), divided by 169. The two cross corrections therefore
# occur in the exact ratio 1:1/169.
COUPLED_FOLDED_THRESHOLD = Fraction(156, 595)
COUPLED_PRODUCT_TAIL = (49 * 595) // 156
require(COUPLED_PRODUCT_TAIL == 186, "coupled product tail changed")
coupled_cross_ancestor_bank = []
for reduced_a, reduced_b, folded_defect in cross_ancestor_bank:
    coupled_invoice = (
        2 * SAME_OWNER
        + 2 * SEPTIMALLY_SEPARATED
        + overlap(reduced_a, reduced_b)
        + overlap(reduced_b, 169 * reduced_a)
    )
    folded_survives = (
        Fraction(folded_defect, reduced_a * reduced_b)
        >= COUPLED_FOLDED_THRESHOLD
    )
    invoice_survives = coupled_invoice >= HOLE_MASS
    require(
        folded_survives == invoice_survives,
        "coupled folded/invoice equivalence changed",
    )
    if folded_survives:
        coupled_cross_ancestor_bank.append(
            (reduced_a, reduced_b, folded_defect)
        )

require(
    len(coupled_cross_ancestor_bank) == 52,
    "coupled cross-ancestor bank size changed",
)
require(
    max(a * b for a, b, _ in coupled_cross_ancestor_bank) == 177,
    "coupled cross-ancestor maximum product changed",
)
require(
    max(a for a, _, _ in coupled_cross_ancestor_bank) == 10,
    "coupled cross-ancestor minimum-coordinate ceiling changed",
)
require(
    max(b for _, b, _ in coupled_cross_ancestor_bank) == 85,
    "coupled cross-ancestor maximum-coordinate ceiling changed",
)


def merged_intervals(intervals):
    """Merge half-open integer intervals on one common endpoint grid."""
    merged = []
    for left, right in sorted(intervals):
        require(left < right, "empty interval entered merger")
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return tuple((left, right) for left, right in merged)


def danger_intervals(speed, modulus):
    """D_speed as half-open cells on Z/modulus; endpoints have zero mass."""
    require(modulus % (14 * speed) == 0, "incompatible endpoint grid")
    step = modulus // speed
    radius = modulus // (14 * speed)
    intervals = []
    for tooth in range(speed):
        center = tooth * step
        left = center - radius
        right = center + radius
        if left < 0:
            intervals.append((0, right))
            intervals.append((left + modulus, modulus))
        elif right > modulus:
            intervals.append((left, modulus))
            intervals.append((0, right - modulus))
        else:
            intervals.append((left, right))
    return merged_intervals(intervals)


def union_intersection_length(first, second):
    """Length of the intersection of two merged interval unions."""
    first_index = 0
    second_index = 0
    total = 0
    while first_index < len(first) and second_index < len(second):
        left = max(first[first_index][0], second[second_index][0])
        right = min(first[first_index][1], second[second_index][1])
        if left < right:
            total += right - left
        if first[first_index][1] <= second[second_index][1]:
            first_index += 1
        else:
            second_index += 1
    return total


def intersection_intervals(first, second):
    """Return the merged half-open intervals in the intersection."""
    first_index = 0
    second_index = 0
    result = []
    while first_index < len(first) and second_index < len(second):
        left = max(first[first_index][0], second[second_index][0])
        right = min(first[first_index][1], second[second_index][1])
        if left < right:
            result.append((left, right))
        if first[first_index][1] <= second[second_index][1]:
            first_index += 1
        else:
            second_index += 1
    return merged_intervals(result)


def scaled_low_cage_intervals(scale, modulus):
    """Intervals of U_(4,3)(scale*x) on one compatible grid."""
    quotient_union = merged_intervals(
        danger_intervals(3 * scale, modulus)
        + danger_intervals(52 * scale, modulus)
    )
    original_union = merged_intervals(
        danger_intervals(39 * scale, modulus)
        + danger_intervals(676 * scale, modulus)
    )
    return intersection_intervals(quotient_union, original_union)


def low_cage_union(reduced_c2, reduced_c1):
    """Measure of (D_C1 union D_C2) intersect (D_c1 union D_c2)."""
    modulus = 14 * 169 * reduced_c2 * reduced_c1
    quotient_union = merged_intervals(
        danger_intervals(reduced_c1, modulus)
        + danger_intervals(13 * reduced_c2, modulus)
    )
    original_union = merged_intervals(
        danger_intervals(13 * reduced_c1, modulus)
        + danger_intervals(169 * reduced_c2, modulus)
    )
    return Fraction(
        union_intersection_length(quotient_union, original_union),
        modulus,
    )


# The C_3 part of the cage costs at most 2/49, so absorption requires
# the exact low-low union to have mass at least 22/343. Both orientations
# of an unordered ancestor ratio must be retained.
LOW_CAGE_THRESHOLD = HOLE_MASS - 2 * SEPTIMALLY_SEPARATED
require(LOW_CAGE_THRESHOLD == Fraction(22, 343), "low-cage threshold changed")
oriented_low_cage_bank = {}
for reduced_a, reduced_b, _ in coupled_cross_ancestor_bank:
    orientations = {(reduced_a, reduced_b), (reduced_b, reduced_a)}
    for reduced_c2, reduced_c1 in sorted(orientations):
        union_mass = low_cage_union(reduced_c2, reduced_c1)
        if union_mass >= LOW_CAGE_THRESHOLD:
            oriented_low_cage_bank[(reduced_c2, reduced_c1)] = union_mass

expected_oriented_low_cage_bank = {
    (1, 1): Fraction(193, 1183),
    (1, 2): Fraction(114, 1183),
    (2, 1): Fraction(239, 2366),
    (1, 3): Fraction(263, 3549),
    (3, 1): Fraction(95, 1183),
    (4, 1): Fraction(331, 4732),
    (2, 3): Fraction(43, 546),
    (3, 2): Fraction(95, 1183),
    (3, 4): Fraction(491, 7098),
    (4, 3): Fraction(331, 4732),
}
require(
    oriented_low_cage_bank == expected_oriented_low_cage_bank,
    "oriented low-cage bank changed",
)
unordered_low_cage_bank = {
    tuple(sorted(pair)) for pair in oriented_low_cage_bank
}
require(
    unordered_low_cage_bank
    == {(1, 1), (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)},
    "unordered low-cage bank changed",
)

# Complete the quantitative dichotomy. Directly scan every admissible
# orientation with ab<=186, not only the 52 pair-sum survivors. Among
# the orientations below the low-cage threshold, the smallest positive
# margin is 9/22295. For ab>186, the coupled folded bound applies.
small_oriented_low_cage_bank = {}
small_oriented_margins = []
for reduced_a in range(1, COUPLED_PRODUCT_TAIL + 1):
    for reduced_b in range(
        reduced_a,
        COUPLED_PRODUCT_TAIL // reduced_a + 1,
    ):
        if gcd(reduced_a, reduced_b) != 1:
            continue
        if reduced_a % 7 == 0 or reduced_b % 7 == 0:
            continue
        if reduced_a % 13 == 0 or reduced_b % 13 == 0:
            continue
        for reduced_c2, reduced_c1 in {
            (reduced_a, reduced_b),
            (reduced_b, reduced_a),
        }:
            union_mass = low_cage_union(reduced_c2, reduced_c1)
            if union_mass >= LOW_CAGE_THRESHOLD:
                small_oriented_low_cage_bank[
                    (reduced_c2, reduced_c1)
                ] = union_mass
            else:
                small_oriented_margins.append(
                    LOW_CAGE_THRESHOLD - union_mass
                )

require(
    small_oriented_low_cage_bank == expected_oriented_low_cage_bank,
    "complete small low-cage bank changed",
)
SMALL_NONBANK_MARGIN = min(small_oriented_margins)
require(
    SMALL_NONBANK_MARGIN == Fraction(9, 22295),
    "small nonbank margin changed",
)

# The best folded quotient after the product-186 cutoff occurs at
# (a,b)=(3,73): Delta/(ab)=48/219=16/73. It is enough to inspect
# products through 223, because 49/224<16/73 controls the tail.
tail_folded_rows = []
for reduced_a in range(1, 224):
    for reduced_b in range(reduced_a, 224):
        product = reduced_a * reduced_b
        if not 187 <= product <= 223:
            continue
        if gcd(reduced_a, reduced_b) != 1:
            continue
        if reduced_a % 7 == 0 or reduced_b % 7 == 0:
            continue
        if reduced_a % 13 == 0 or reduced_b % 13 == 0:
            continue
        tail_folded_rows.append(
            (
                Fraction(
                    fold14(reduced_a + reduced_b)
                    - fold14(reduced_b - reduced_a),
                    product,
                ),
                reduced_a,
                reduced_b,
            )
        )
tail_folded_max = max(row[0] for row in tail_folded_rows)
require(tail_folded_max == Fraction(16, 73), "folded tail maximum changed")
require(
    {(a, b) for value, a, b in tail_folded_rows if value == tail_folded_max}
    == {(3, 73)},
    "folded tail equality locus changed",
)
require(Fraction(49, 224) < tail_folded_max, "folded analytic tail failed")
TAIL_NONBANK_MARGIN = (
    Fraction(170, 169 * 196)
    * (COUPLED_FOLDED_THRESHOLD - tail_folded_max)
)
require(
    TAIL_NONBANK_MARGIN == Fraction(934, 4231591),
    "tail nonbank margin changed",
)
B2_NONBANK_MARGIN = min(SMALL_NONBANK_MARGIN, TAIL_NONBANK_MARGIN)
require(
    B2_NONBANK_MARGIN == TAIL_NONBANK_MARGIN,
    "b=2 nonbank margin changed",
)
require(
    B2_NONBANK_MARGIN > LIVE_REPEATED_MARGIN,
    "b=2 nonbank fell below the 150-profile uniform floor",
)

# THM-2377's balanced septimal descendant becomes a ten-row exact bank.
# With c1=13hb and c2=169ha, choose balanced rho with 13a=rho*b mod 7;
# then (c2-rho*c1)/7=13h*d.
def balanced_seven(residue):
    residue %= 7
    require(residue != 0, "balanced seven called at zero")
    return residue if residue <= 3 else residue - 7


bockstein_descendants = {}
for reduced_c2, reduced_c1 in oriented_low_cage_bank:
    inverse_c1 = pow(reduced_c1, -1, 7)
    rho = balanced_seven(13 * reduced_c2 * inverse_c1)
    numerator = 13 * reduced_c2 - rho * reduced_c1
    require(numerator % 7 == 0, "Bockstein descendant lost integrality")
    bockstein_descendants[(reduced_c2, reduced_c1)] = (
        rho,
        numerator // 7,
    )
require(
    bockstein_descendants
    == {
        (1, 1): (-1, 2),
        (1, 2): (3, 1),
        (2, 1): (-2, 4),
        (1, 3): (2, 1),
        (3, 1): (-3, 6),
        (4, 1): (3, 7),
        (2, 3): (-3, 5),
        (3, 2): (2, 5),
        (3, 4): (1, 5),
        (4, 3): (1, 7),
    },
    "Bockstein descendant bank changed",
)

# At septimal depth M>=2, THM-2391 makes the two low-blocker steps
# have one common absolute value. Since C_2/C_1=13a/b, the ten-bank
# address must satisfy 13a=+/-b modulo 7^M. Only (4,3) survives
# modulo 49, and no address survives modulo 343.
def thm2391_address_survivors(modulus):
    return tuple(
        (reduced_c2, reduced_c1)
        for reduced_c2, reduced_c1 in oriented_low_cage_bank
        if (
            (13 * reduced_c2 - reduced_c1) % modulus == 0
            or (13 * reduced_c2 + reduced_c1) % modulus == 0
        )
    )


M2_ADDRESS_BANK = thm2391_address_survivors(49)
M3_ADDRESS_BANK = thm2391_address_survivors(343)
require(M2_ADDRESS_BANK == ((4, 3),), "M=2 address bank changed")
require(M3_ADDRESS_BANK == (), "M>=3 address bank changed")

# THM-2391 reduces the M=2 no-clean bank to the sole orientation
# (C_2:c_1)=(4:3). Write its common low scale and top speed, after
# dividing their gcd, as coprime p,n. Then p is a 7*13-unit, while n
# is divisible by 49 and is a 13-unit. The compatible low union is
#
#   U_p=(D_(3p) union D_(52p))
#       intersection (D_(39p) union D_(676p)).
#
# Excluding the top comb removes U_p intersection D_n. The unexcluded
# low union plus the two high cross pieces exceeds the hole floor by
# only 1347/231868, so any larger top intersection forces a clean hole.
M2_LOW_UNION = expected_oriented_low_cage_bank[(4, 3)]
M2_CAGE_EXCESS = (
    M2_LOW_UNION + 2 * SEPTIMALLY_SEPARATED - HOLE_MASS
)
require(M2_LOW_UNION == Fraction(331, 4732), "M=2 low union changed")
require(
    M2_CAGE_EXCESS == Fraction(1347, 231868),
    "M=2 cage excess changed",
)

# On the base endpoint grid, the linear interval list has 139 pieces;
# its first and last pieces join through the circle endpoint. Thus the
# whole union has 138 circular components and total variation 276.
M2_BASE_MODULUS = 14 * 2028
m2_base_intervals = scaled_low_cage_intervals(1, M2_BASE_MODULUS)
require(len(m2_base_intervals) == 139, "M=2 linear interval count changed")
require(
    m2_base_intervals[0][0] == 0
    and m2_base_intervals[-1][1] == M2_BASE_MODULUS,
    "M=2 circular boundary component changed",
)
M2_CIRCULAR_COMPONENTS = len(m2_base_intervals) - 1
M2_VARIATION = 2 * M2_CIRCULAR_COMPONENTS
require(M2_CIRCULAR_COMPONENTS == 138, "M=2 component count changed")
require(M2_VARIATION == 276, "M=2 variation changed")
require(
    Fraction(
        sum(right - left for left, right in m2_base_intervals),
        M2_BASE_MODULUS,
    )
    == M2_LOW_UNION,
    "M=2 base interval mass changed",
)

# Fourier gives
#
# |mu(U_p intersect D_n)-mu(U)/7|
#  <= sum_(t!=0) [138/(pi*n*|t|)] [1/(pi*p*|t|)]
#  =46/(np).
#
# The equilibrium share beats the cage excess by 485/115934.
M2_EQUILIBRIUM_RESERVE = M2_LOW_UNION / 7 - M2_CAGE_EXCESS
require(
    M2_EQUILIBRIUM_RESERVE == Fraction(485, 115934),
    "M=2 equilibrium reserve changed",
)
M2_TAIL_START = 10996
require(
    M2_EQUILIBRIUM_RESERVE
    - Fraction(46, M2_TAIL_START)
    == Fraction(12, 159351283),
    "M=2 BV tail reserve changed",
)
require(
    M2_EQUILIBRIUM_RESERVE
    - Fraction(46, M2_TAIL_START - 1)
    < 0,
    "M=2 BV cutoff is no longer first",
)


def scaled_m2_base_intervals(reduced_low_scale, modulus):
    """Lift the already merged base union through x -> p*x."""
    denominator = M2_BASE_MODULUS * reduced_low_scale
    require(modulus % denominator == 0, "M=2 lifted grid changed")
    factor = modulus // denominator
    return tuple(
        (
            (sheet * M2_BASE_MODULUS + left) * factor,
            (sheet * M2_BASE_MODULUS + right) * factor,
        )
        for sheet in range(reduced_low_scale)
        for left, right in m2_base_intervals
    )


def m2_top_intersection(reduced_low_scale, reduced_top_speed):
    """Exact mu(U_p intersection D_n) for one coprime primitive pair."""
    modulus = 14 * lcm(
        2028 * reduced_low_scale,
        reduced_top_speed,
    )
    low_union = scaled_m2_base_intervals(reduced_low_scale, modulus)
    top_comb = danger_intervals(reduced_top_speed, modulus)
    return Fraction(
        union_intersection_length(low_union, top_comb),
        modulus,
    )


# Since n=49m, the BV-open core np<10996 is pm<=224. Exhaust the exact
# primitive universe, including the stronger higher-seven-depth controls.
m2_finite_rows = []
for reduced_low_scale in range(1, 225):
    if reduced_low_scale % 7 == 0 or reduced_low_scale % 13 == 0:
        continue
    for top_multiplier in range(
        1,
        224 // reduced_low_scale + 1,
    ):
        reduced_top_speed = 49 * top_multiplier
        if reduced_top_speed % 13 == 0:
            continue
        if gcd(reduced_low_scale, reduced_top_speed) != 1:
            continue
        m2_finite_rows.append(
            (
                m2_top_intersection(
                    reduced_low_scale,
                    reduced_top_speed,
                ),
                reduced_low_scale,
                reduced_top_speed,
            )
        )

require(len(m2_finite_rows) == 753, "M=2 finite pair count changed")
M2_TOP_INTERSECTION, m2_min_p, m2_min_n = min(m2_finite_rows)
require(
    (M2_TOP_INTERSECTION, m2_min_p, m2_min_n)
    == (Fraction(1849, 231868), 1, 49),
    "M=2 exact top-intersection minimum changed",
)
require(
    {
        (reduced_low_scale, reduced_top_speed)
        for value, reduced_low_scale, reduced_top_speed in m2_finite_rows
        if value == M2_TOP_INTERSECTION
    }
    == {(1, 49)},
    "M=2 minimum is no longer unique",
)
M2_CLEAN_MARGIN = M2_TOP_INTERSECTION - M2_CAGE_EXCESS
require(
    M2_CLEAN_MARGIN == Fraction(251, 115934),
    "M=2 clean-hole margin changed",
)

LIVE_ESCAPED_PROFILE_COUNT = len(repeated_rows) + len(live_strict_rows)
require(LIVE_ESCAPED_PROFILE_COUNT == 150, "live escaped profile count changed")
LIVE_UNIFORM_MARGIN = min(LIVE_REPEATED_MARGIN, live_strict_margin)
require(
    LIVE_UNIFORM_MARGIN == LIVE_REPEATED_MARGIN,
    "live uniform margin changed",
)

# A general deterministic ordinary-label deletion partitions by a
# lower owner, the double pair, and the adjacent word translate.
FALLBACK_PATTERN_COUNT = 2 * 15 * 13
FALLBACK_CELL_MASS = LIVE_UNIFORM_MARGIN / FALLBACK_PATTERN_COUNT
require(FALLBACK_PATTERN_COUNT == 390, "fallback pattern count changed")
require(
    FALLBACK_CELL_MASS == Fraction(1, 10434060),
    "fallback cell mass floor changed",
)

# The canonical sharpening deletes the distinguished top label q_*.
# If q_* is outside the unique double pair, the complement of the other
# five unit masks is its adjacent two-root word. If q_* belongs to the
# double pair, that complement is the one exclusive q_* root. Thus the
# only cell data are a deterministic low owner (two choices), the
# singleton/adjacent status (two choices), and its support (thirteen
# choices in either status).
TOP_PATTERN_COUNT = 2 * 2 * 13
TOP_CELL_MASS = LIVE_UNIFORM_MARGIN / TOP_PATTERN_COUNT
require(TOP_PATTERN_COUNT == 52, "top-labelled pattern count changed")
require(
    TOP_CELL_MASS == Fraction(1, 1391208),
    "top-labelled cell mass floor changed",
)

# THM-2391 supplies on the same parent a seven-root lower-load word
# L_7-1=1_{s=d}.  The exact owner categories are: d=0 with both low
# blockers active, or d in F_7^* together with one of two exclusive
# owners. Combining these thirteen categories with the 26 possible
# singleton/adjacent q_* supports gives an owner-resolved C_7 x C_13
# tensor cell.
OWNER_ADDRESS_COUNT = 1 + 6 * 2
TOP_SUPPORT_STATUS_COUNT = 2 * 13
TENSOR_PATTERN_COUNT = OWNER_ADDRESS_COUNT * TOP_SUPPORT_STATUS_COUNT
TENSOR_CELL_MASS = LIVE_UNIFORM_MARGIN / TENSOR_PATTERN_COUNT
require(OWNER_ADDRESS_COUNT == 13, "owner/address count changed")
require(TOP_SUPPORT_STATUS_COUNT == 26, "top support/status count changed")
require(TENSOR_PATTERN_COUNT == 338, "tensor pattern count changed")
require(
    TENSOR_CELL_MASS == Fraction(1, 9042852),
    "owner-resolved tensor cell mass floor changed",
)

# For normalized C_13 Fourier transform, a singleton word has twelve
# equal nonzero modes. A two-adjacent-root word has the displayed
# Parseval and fourth-power ledgers. Both are nonzero in every nonzero
# target colour.
SINGLETON_ENERGY = Fraction(12, 169)
SINGLETON_FOURTH = Fraction(12, 28561)
TWO_ROOT_ENERGY = Fraction(22, 169)
TWO_ROOT_FOURTH = Fraction(62, 28561)
TOP_CELL_ENERGY_FLOOR = TOP_CELL_MASS ** 2 * SINGLETON_ENERGY
TOP_CELL_FOURTH_FLOOR = TOP_CELL_MASS ** 4 * SINGLETON_FOURTH

# With normalized transforms in both cyclic variables, the fixed
# tensor is rho*(delta_d tensor A). Every one of its 7*12 coefficients
# with arbitrary septimal colour and nonzero target colour survives.
# Summing their squared moduli contributes rho^2 E(A)/7.
TENSOR_ENERGY_FLOOR = (
    TENSOR_CELL_MASS ** 2 * SINGLETON_ENERGY / 7
)
TENSOR_FOURTH_FLOOR = (
    TENSOR_CELL_MASS ** 4 * SINGLETON_FOURTH / (7 ** 3)
)

# After promoted THM-2393, the only no-clean boundary is the M=1
# common-core chain. Its high-safe base has mass 396/637, while the
# one-address low-union base has mass 145/169. On their intersection,
# total lower-unit incidence six forces the unique U-address to be the
# sole K-hole. Root disintegration divides its base mass by seven.
COMMON_CORE_HIGH_SAFE = Fraction(396, 637)
COMMON_CORE_ONE_ADDRESS = Fraction(145, 169)
COMMON_CORE_FORCED_BASE = (
    COMMON_CORE_HIGH_SAFE + COMMON_CORE_ONE_ADDRESS - 1
)
COMMON_CORE_FORCED_ROOT = COMMON_CORE_FORCED_BASE / 7
require(
    COMMON_CORE_FORCED_BASE == Fraction(3972, 8281),
    "common-core forced base mass changed",
)
require(
    COMMON_CORE_FORCED_ROOT == Fraction(3972, 57967),
    "common-core forced root mass changed",
)

print("LRC14 THM-2388 blocker-cage capacity escape -- exact candidate")
print("dependency_status: THM-2388_PROVED_AUDITED")
print(f"hole_mass_floor: {HOLE_MASS}")
print(f"deep_gap_max: d={deep_worst_gap},cap={deep_worst_cap}")
print(f"repeated_profiles_escaped: {len(repeated_rows)}/15")
print(f"repeated_worst_profile: {repeated_worst[0]}")
print(f"repeated_cage_cap: {repeated_worst[1]}")
print(f"repeated_off_cage_margin: {REPEATED_MARGIN}")
print(f"strict_profiles_escaped_by_pair_caps: {len(strict_pass)}/150")
print(f"strict_zero_gap_profiles_open: {len(strict_zero_gap)}")
print(
    "strict_positive_gap_coarse_failures: "
    + ",".join(str(profile) for profile, *_ in strict_coarse_fail)
)
print(
    "compatible_pair_sum_hostiles: "
    + ",".join(
        f"{profile}:{invoice}"
        for profile, (_, invoice, _) in compatibility_hostiles.items()
    )
)
print(
    "thirteen_adic_pair_only_escaped_profiles: "
    f"{pair_only_escaped_profile_count}/165"
)
print(f"septimally_separated_pair_overlap: {SEPTIMALLY_SEPARATED}")
print(f"live_repeated_off_cage_margin: {LIVE_REPEATED_MARGIN}")
print(f"live_strict_profiles_escaped: {len(live_strict_rows)}/150")
print(f"live_strict_worst_b_family: b=4,count={len(live_strict_worst)}")
print(f"live_strict_off_cage_margin: {live_strict_margin}")
print(f"live_escaped_profiles_total: {LIVE_ESCAPED_PROFILE_COUNT}/165")
print(f"live_uniform_off_cage_margin: {LIVE_UNIFORM_MARGIN}")
print(f"b2_septimally_distinct_margin: {B2_SEPTIMALLY_DISTINCT_MARGIN}")
print("residual_lane: b=2_and_nu7(C2)=nu7(c1)")
print(f"cross_ancestor_overlap_threshold: {CROSS_ANCESTOR_THRESHOLD}")
print(f"cross_ancestor_product_tail: {CROSS_ANCESTOR_PRODUCT_TAIL}")
print(f"cross_ancestor_unordered_ratio_bank: {len(cross_ancestor_bank)}")
print("cross_ancestor_ratio_bounds: min_coordinate<=11,max_coordinate<=197")
print(
    "cross_ancestor_first_coordinate_census: "
    + ",".join(f"{key}:{value}" for key, value in cross_ancestor_counts.items())
)
print(f"coupled_cross_ancestor_folded_threshold: {COUPLED_FOLDED_THRESHOLD}")
print(f"coupled_cross_ancestor_product_tail: {COUPLED_PRODUCT_TAIL}")
print(
    "coupled_cross_ancestor_unordered_ratio_bank: "
    f"{len(coupled_cross_ancestor_bank)}"
)
print("coupled_cross_ancestor_ratio_bounds: product<=177,a<=10,b<=85")
print(f"low_cage_union_threshold: {LOW_CAGE_THRESHOLD}")
print(f"low_cage_oriented_ratio_bank: {len(oriented_low_cage_bank)}")
print(f"low_cage_unordered_ratio_bank: {len(unordered_low_cage_bank)}")
print(
    "low_cage_oriented_survivors: "
    + ",".join(
        f"{pair}:{mass}"
        for pair, mass in oriented_low_cage_bank.items()
    )
)
print(f"small_nonbank_off_cage_margin: {SMALL_NONBANK_MARGIN}")
print(f"folded_tail_max: {tail_folded_max}@3:73")
print(f"tail_nonbank_off_cage_margin: {TAIL_NONBANK_MARGIN}")
print(f"b2_nonbank_uniform_off_cage_margin: {B2_NONBANK_MARGIN}")
print(f"all_nonexceptional_uniform_off_cage_margin: {LIVE_UNIFORM_MARGIN}")
print(
    "oriented_bockstein_descendants: "
    + ",".join(
        f"{pair}:rho={rho},d={descendant}"
        for pair, (rho, descendant) in bockstein_descendants.items()
    )
)
print(f"m2_address_bank: {M2_ADDRESS_BANK}")
print(f"m3_address_bank: {M3_ADDRESS_BANK}")
print(f"m2_low_union: {M2_LOW_UNION}")
print(f"m2_cage_excess: {M2_CAGE_EXCESS}")
print(f"m2_circular_components: {M2_CIRCULAR_COMPONENTS}")
print(f"m2_whole_union_variation: {M2_VARIATION}")
print(f"m2_bv_tail_start: {M2_TAIL_START}")
print(f"m2_bv_tail_reserve: {Fraction(12, 159351283)}")
print(f"m2_finite_primitive_pairs: {len(m2_finite_rows)}")
print(f"m2_min_top_intersection: {M2_TOP_INTERSECTION}@1:49")
print(f"m2_uniform_clean_margin: {M2_CLEAN_MARGIN}")
print(f"fallback_pattern_count: {FALLBACK_PATTERN_COUNT}")
print(f"fallback_two_root_cell_mass: {FALLBACK_CELL_MASS}")
print(f"top_labelled_pattern_count: {TOP_PATTERN_COUNT}")
print(f"top_labelled_cell_mass: {TOP_CELL_MASS}")
print(f"owner_address_count: {OWNER_ADDRESS_COUNT}")
print(f"owner_resolved_tensor_pattern_count: {TENSOR_PATTERN_COUNT}")
print(f"owner_resolved_tensor_cell_mass: {TENSOR_CELL_MASS}")
print(f"singleton_nonzero_target_energy: {SINGLETON_ENERGY}")
print(f"singleton_nonzero_target_fourth_energy: {SINGLETON_FOURTH}")
print(f"two_root_nonzero_target_energy: {TWO_ROOT_ENERGY}")
print(f"two_root_nonzero_target_fourth_energy: {TWO_ROOT_FOURTH}")
print(f"top_cell_nonzero_target_energy_floor: {TOP_CELL_ENERGY_FLOOR}")
print(f"top_cell_nonzero_target_fourth_floor: {TOP_CELL_FOURTH_FLOOR}")
print(f"tensor_nonzero_target_energy_floor: {TENSOR_ENERGY_FLOOR}")
print(f"tensor_nonzero_target_fourth_floor: {TENSOR_FOURTH_FLOOR}")
print(f"common_core_high_safe_base: {COMMON_CORE_HIGH_SAFE}")
print(f"common_core_one_address_base: {COMMON_CORE_ONE_ADDRESS}")
print(f"common_core_forced_middle_base: {COMMON_CORE_FORCED_BASE}")
print(f"common_core_forced_middle_root: {COMMON_CORE_FORCED_ROOT}")
print("common_core_address_transport: D_13h_to_D_h")
print("terminal_owner_current_transport: OPEN")
print("all_checks: PASS")
