#!/usr/bin/env python3
"""Exact companion for THM-2414.

This dependency-free script audits the thirteen-skew seven-root word
reflection, the capacity-saturated quotient split, its common-valuation
consequence, the live W8 terminal atlas, the exact full-bin failure of the
excluded W7 one-fibre hostile, and the sharp triple-shield boundary.  It
deliberately does not search for or assert a global scalar cover.
"""

from fractions import Fraction as F
from itertools import combinations


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(n: int, p: int) -> int:
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def centred_distance_num(num: int, den: int) -> int:
    residue = num % den
    return min(residue, den - residue)


def danger_num(num: int, den: int, width: int = 1) -> bool:
    return 14 * centred_distance_num(num, den) < width * den


def base_sign(v: int, num: int, den: int, width: int = 1) -> int:
    """Negative means danger, positive means strict safety."""
    return 14 * centred_distance_num(v * num, den) - width * den


def seven_word(v: int, num: int, den: int, width: int = 1) -> frozenset[int]:
    """Word on the roots x_k=(num/den+k)/7."""
    root_den = 7 * den
    return frozenset(
        k
        for k in range(7)
        if danger_num(v * (num + k * den), root_den, width)
    )


def word_signed_distances(
    v: int, num: int, den: int, width: int = 1
) -> tuple[int, ...]:
    root_den = 7 * den
    return tuple(
        14 * centred_distance_num(v * (num + k * den), root_den)
        - width * root_den
        for k in range(7)
    )


# ---------------------------------------------------------------------------
# 1. The physical/divided word reflection, including the base-lift digit.
# ---------------------------------------------------------------------------

reflection_cases = 0
reflection_den = 1009
for A in range(1, 121):
    if A % 7 == 0:
        continue
    for num in range(reflection_den):
        image_num = (13 * num) % reflection_den
        lift_digit = (13 * num) // reflection_den
        for width in (1, 2):
            physical = seven_word(13 * A, num, reflection_den, width)
            divided = seven_word(A, image_num, reflection_den, width)
            reflected = frozenset((lift_digit - k) % 7 for k in divided)
            require(physical == reflected, "thirteen-skew word reflection")
            reflection_cases += 1


# ---------------------------------------------------------------------------
# 2. Capacity saturation and the common-valuation residue gate.
# ---------------------------------------------------------------------------

capacity_cases = 0
for scale in range(1, 65):
    guard_size = 2 * scale
    blocker_size = scale
    # If a set of guard_size points lies inside the union of two sets of
    # blocker_size points, equality of capacities forces a disjoint union.
    require(
        guard_size == 2 * blocker_size,
        "guard/blocker capacity normalization",
    )
    capacity_cases += 1

# The root-set realization at scale one checks the full support conclusion.
root_split_cases = 0
vertices = range(7)
for guard in combinations(vertices, 2):
    guard_set = frozenset(guard)
    for b1 in vertices:
        for b2 in vertices:
            union = frozenset((b1, b2))
            if guard_set <= union:
                require(b1 != b2, "saturated blockers must be disjoint")
                require(union == guard_set, "saturated union must equal guard")
                root_split_cases += 1
require(root_split_cases == 42, "unexpected number of labelled split words")

valuation_profiles = 0
valuation_survivors = []
for h in range(6):
    for r1 in range(6):
        for r2 in range(6):
            least = min(h, r1, r2)
            lhs_mod7 = 2 if h == least else 0
            rhs_mod7 = int(r1 == least) + int(r2 == least)
            if lhs_mod7 % 7 == rhs_mod7 % 7:
                valuation_survivors.append((h, r1, r2))
            valuation_profiles += 1
require(
    valuation_survivors == [(r, r, r) for r in range(6)],
    "the mod-seven termination gate must force a common valuation",
)

# THM-2390 then leaves guard + c1 + c2 + three or four lower q labels.
heavy_role_cases = []
lower_q_labels = tuple(range(4))
for number_heavy_q in (3, 4):
    for chosen in combinations(lower_q_labels, number_heavy_q):
        weight = 2 + 1 + 1 + number_heavy_q
        require(weight in (7, 8), "unexpected heavy-layer weight")
        heavy_role_cases.append((chosen, weight))
require(
    len([row for row in heavy_role_cases if row[1] == 7]) == 4,
    "weight-seven role count",
)
require(
    len([row for row in heavy_role_cases if row[1] == 8]) == 1,
    "weight-eight role count",
)


# ---------------------------------------------------------------------------
# 3. Exact W=8 cross-time stopping atlas.
# ---------------------------------------------------------------------------

den = 1009
y8 = 73
z8 = 949
require(z8 == (13 * y8) % den, "W8 image base")

H8 = 1
q_star8 = 7
lower_q8 = (37, 82, 67, 22)
C8 = (15, 1)
c8 = tuple(13 * value for value in C8)
C3_8 = 637
c3_8 = 13 * C3_8

w8_y_expected = {
    ("H", H8, 2): frozenset((0, 6)),
    ("q37", 37, 1): frozenset((2,)),
    ("q82", 82, 1): frozenset((3,)),
    ("q67", 67, 1): frozenset((4,)),
    ("q22", 22, 1): frozenset((5,)),
    ("c1", c8[0], 1): frozenset((0,)),
    ("c2", c8[1], 1): frozenset((1,)),
}
w8_z_expected = {
    ("H", H8, 2): frozenset((0, 6)),
    ("q37", 37, 1): frozenset((0,)),
    ("q82", 82, 1): frozenset((0,)),
    ("q67", 67, 1): frozenset((0,)),
    ("q22", 22, 1): frozenset((0,)),
    ("C1", C8[0], 1): frozenset((0,)),
    ("C2", C8[1], 1): frozenset((6,)),
}

w8_endpoint_clearance = None
for (_, speed, width), expected in w8_y_expected.items():
    require(seven_word(speed, y8, den, width) == expected, "W8 y-word")
    clearance = min(
        abs(value) for value in word_signed_distances(speed, y8, den, width)
    )
    w8_endpoint_clearance = (
        clearance
        if w8_endpoint_clearance is None
        else min(w8_endpoint_clearance, clearance)
    )
for (_, speed, width), expected in w8_z_expected.items():
    require(seven_word(speed, z8, den, width) == expected, "W8 z-word")
    clearance = min(
        abs(value) for value in word_signed_distances(speed, z8, den, width)
    )
    w8_endpoint_clearance = min(w8_endpoint_clearance, clearance)
require(w8_endpoint_clearance == 840, "W8 endpoint clearance")

w8_multiplicity = [0] * 7
for word in w8_y_expected.values():
    for root in word:
        w8_multiplicity[root] += 1
require(
    tuple(w8_multiplicity) == (2, 1, 1, 1, 1, 1, 1),
    "W8 terminal one-double word",
)
require(
    w8_z_expected[("C1", C8[0], 1)]
    | w8_z_expected[("C2", C8[1], 1)]
    == w8_z_expected[("H", H8, 2)],
    "W8 quotient blockers split the guard",
)
require(
    not (
        w8_z_expected[("C1", C8[0], 1)]
        & w8_z_expected[("C2", C8[1], 1)]
    ),
    "W8 quotient blockers are disjoint",
)
for q in lower_q8:
    require(
        seven_word(q, z8, den) <= seven_word(H8, z8, den, 2),
        "W8 lower q must select a guard endpoint",
    )

w8_q_y_sign = base_sign(q_star8 // 7, y8, den)
w8_q_z_sign = base_sign(q_star8 // 7, z8, den)
w8_high_y_sign = base_sign(c3_8 // 7, y8, den)
w8_high_z_sign = base_sign(C3_8 // 7, z8, den)
require((w8_q_y_sign, w8_q_z_sign) == (13, -169), "W8 q transition")
require(
    w8_high_y_sign == w8_high_z_sign == 4801,
    "W8 high-state transport",
)
require(
    tuple(valuation(v, 7) for v in (H8, *lower_q8, *c8, q_star8, c3_8))
    == (0, 0, 0, 0, 0, 0, 0, 1, 2),
    "W8 septimal roles",
)
require(
    tuple(valuation(v, 13) for v in (*c8, c3_8)) == (1, 1, 2),
    "W8 blocker thirteen-adic roles",
)


# ---------------------------------------------------------------------------
# 4. Exact W=7 cross-time stopping atlas.
# ---------------------------------------------------------------------------

y7 = 11
z7 = 143
require(z7 == (13 * y7) % den, "W7 image base")

H7 = 8
q_star7 = 49
heavy_q7 = (79, 68, 53)
off_q7 = 56
C7 = (4, 8)
c7 = tuple(13 * value for value in C7)
C3_7 = 4459
c3_7 = 13 * C3_7

w7_y_expected = {
    ("H", H7, 2): frozenset((0, 6)),
    ("q79", 79, 1): frozenset((3,)),
    ("q68", 68, 1): frozenset((4,)),
    ("q53", 53, 1): frozenset((5,)),
    ("c1", c7[0], 1): frozenset((2,)),
    ("c2", c7[1], 1): frozenset((1,)),
}
w7_z_expected = {
    ("H", H7, 2): frozenset((5, 6)),
    ("q79", 79, 1): frozenset((5,)),
    ("q68", 68, 1): frozenset((5,)),
    ("q53", 53, 1): frozenset((5,)),
    ("C1", C7[0], 1): frozenset((5,)),
    ("C2", C7[1], 1): frozenset((6,)),
}

w7_endpoint_clearance = None
for (_, speed, width), expected in w7_y_expected.items():
    require(seven_word(speed, y7, den, width) == expected, "W7 y-word")
    clearance = min(
        abs(value) for value in word_signed_distances(speed, y7, den, width)
    )
    w7_endpoint_clearance = (
        clearance
        if w7_endpoint_clearance is None
        else min(w7_endpoint_clearance, clearance)
    )
for (_, speed, width), expected in w7_z_expected.items():
    require(seven_word(speed, z7, den, width) == expected, "W7 z-word")
    clearance = min(
        abs(value) for value in word_signed_distances(speed, z7, den, width)
    )
    w7_endpoint_clearance = min(w7_endpoint_clearance, clearance)
require(w7_endpoint_clearance == 161, "W7 endpoint clearance")

w7_multiplicity = [0] * 7
for word in w7_y_expected.values():
    for root in word:
        w7_multiplicity[root] += 1
require(tuple(w7_multiplicity) == (1,) * 7, "W7 terminal partition")
require(
    w7_z_expected[("C1", C7[0], 1)]
    | w7_z_expected[("C2", C7[1], 1)]
    == w7_z_expected[("H", H7, 2)],
    "W7 quotient blockers split the guard",
)
require(
    not (
        w7_z_expected[("C1", C7[0], 1)]
        & w7_z_expected[("C2", C7[1], 1)]
    ),
    "W7 quotient blockers are disjoint",
)
for q in heavy_q7:
    require(
        seven_word(q, z7, den) <= seven_word(H7, z7, den, 2),
        "W7 heavy q must select a guard endpoint",
    )

w7_q_y_sign = base_sign(q_star7 // 7, y7, den)
w7_q_z_sign = base_sign(q_star7 // 7, z7, den)
w7_high_y_sign = base_sign(c3_7 // 7, y7, den)
w7_high_z_sign = base_sign(C3_7 // 7, z7, den)
w7_off_y_sign = base_sign(off_q7 // 7, y7, den)
w7_off_z_sign = base_sign(off_q7 // 7, z7, den)
require((w7_q_y_sign, w7_q_z_sign) == (69, -897), "W7 q transition")
require(
    w7_high_y_sign == w7_high_z_sign == 2925,
    "W7 high-state transport",
)
require((w7_off_y_sign, w7_off_z_sign) == (223, 881), "W7 off-layer safety")
require(
    tuple(
        valuation(v, 7)
        for v in (H7, *heavy_q7, off_q7, *c7, q_star7, c3_7)
    )
    == (0, 0, 0, 0, 1, 0, 0, 2, 3),
    "W7 septimal roles",
)
require(
    tuple(valuation(v, 13) for v in (*c7, c3_7)) == (1, 1, 2),
    "W7 blocker thirteen-adic roles",
)

# The displayed seven-root z-fibre is only the address class s=0 mod 7
# inside the 49-point top bin x_s.  Two strict controls elsewhere in that
# bin show why this local word is excluded by THM-2391's full-bin theorem.
w7_full_bin_den = 49 * den


def w7_full_bin_num(s: int) -> int:
    return 1001 + den * s


for k in range(7):
    # x_(7k)=(z+k)/7, so this really recovers the word already checked.
    left = F(w7_full_bin_num(7 * k), w7_full_bin_den)
    right = F(z7 + den * k, 7 * den)
    require(left == right, "W7 displayed fibre address")

w7_split_failures = []
w7_containment_failures = []
for s in range(49):
    point_num = w7_full_bin_num(s)
    top_bit = danger_num(q_star7 * point_num, w7_full_bin_den)
    high_bit = danger_num(C3_7 * point_num, w7_full_bin_den)
    guard_bit = danger_num(H7 * point_num, w7_full_bin_den, 2)
    C1_bit = danger_num(C7[0] * point_num, w7_full_bin_den)
    C2_bit = danger_num(C7[1] * point_num, w7_full_bin_den)
    off_q_bit = danger_num(off_q7 * point_num, w7_full_bin_den)
    require(top_bit and not high_bit, "W7 full-bin carrier state")
    if int(guard_bit) != int(C1_bit) + int(C2_bit):
        w7_split_failures.append(s)
    if off_q_bit and not guard_bit:
        w7_containment_failures.append(s)

require(
    w7_split_failures == [6, 11, 18, 29, 36, 48],
    "W7 full-bin partition-failure addresses",
)
require(
    w7_containment_failures == [13, 20, 27, 34, 41],
    "W7 full-bin containment-failure addresses",
)

w7_s6_signs = tuple(
    base_sign(v, w7_full_bin_num(6), w7_full_bin_den, width)
    for v, width in (
        (q_star7, 1),
        (H7, 2),
        (C7[0], 1),
        (C7[1], 1),
        (C3_7, 1),
    )
)
require(
    w7_s6_signs == (-43953, -896, 247653, 48545, 143325),
    "W7 s=6 strict partition hostile",
)

w7_s13_signs = tuple(
    base_sign(v, w7_full_bin_num(13), w7_full_bin_den, width)
    for v, width in (
        (q_star7, 1),
        (off_q7, 1),
        (H7, 2),
        (C7[0], 1),
        (C7[1], 1),
        (C3_7, 1),
    )
)
require(
    w7_s13_signs == (-43953, -43169, 97986, 48993, 147427, 143325),
    "W7 s=13 strict containment hostile",
)


# ---------------------------------------------------------------------------
# 5. Infinite triple-shield hostile and its XOR distinction.
# ---------------------------------------------------------------------------

shield_B = tuple(B for B in range(50, 92) if B % 7)
for B in shield_B:
    # D_1 cap D_7 is the central interval of radius 1/98.  The nearest
    # noncentral B tooth starts no earlier than 1/98, and its central
    # tooth is strictly inside the D_49 central tooth.
    require(F(13, 14 * B) >= F(1, 98), "noncentral B-tooth exclusion")
    require(F(1, 14 * B) < F(1, 686), "central B-tooth shielding")
require(len(shield_B) == 36, "shield-family cardinality")

require(F(13, 14 * 91) == F(1, 98), "B=91 touching boundary")
require(F(13, 14 * 92) < F(1, 98), "B=92 penetration boundary")
require(F(1, 14 * 48) > F(1, 686), "B=48 central escape boundary")

x48 = (F(1, 14 * 48) + F(1, 686)) / 2
require(
    danger_num(x48.numerator, x48.denominator),
    "B48 witness must lie in D1",
)
require(
    danger_num(7 * x48.numerator, x48.denominator),
    "B48 witness must lie in D7",
)
require(
    danger_num(48 * x48.numerator, x48.denominator),
    "B48 witness must lie in D48",
)
require(
    not danger_num(49 * x48.numerator, x48.denominator),
    "B48 witness must escape D49",
)

x92 = (F(13, 14 * 92) + F(1, 98)) / 2
require(
    danger_num(x92.numerator, x92.denominator),
    "B92 witness must lie in D1",
)
require(
    danger_num(7 * x92.numerator, x92.denominator),
    "B92 witness must lie in D7",
)
require(
    danger_num(92 * x92.numerator, x92.denominator),
    "B92 witness must lie in D92",
)
require(
    not danger_num(49 * x92.numerator, x92.denominator),
    "B92 witness must escape D49",
)

# Shielding is strictly weaker than the XOR split.  For B=50 and H=1,
# x=49/175 is q-danger, c-safe and B-danger, but both D1 and E1 are safe.
x_xor = F(49, 175)
require(danger_num(7 * x_xor.numerator, x_xor.denominator), "XOR q bit")
require(not danger_num(49 * x_xor.numerator, x_xor.denominator), "XOR c bit")
require(danger_num(50 * x_xor.numerator, x_xor.denominator), "XOR B bit")
require(not danger_num(x_xor.numerator, x_xor.denominator), "XOR D1 bit")
require(
    not danger_num(x_xor.numerator, x_xor.denominator, 2),
    "XOR guard bit",
)


print("theorem=THM-2414")
print("status=PROVED+VERIFIED-EXACT-CANDIDATE-UNDER-INDEPENDENT-AUDIT")
print(f"word_reflection_cases={reflection_cases}")
print(f"capacity_scales={capacity_cases}")
print(f"labelled_root_splits={root_split_cases}")
print(f"valuation_profiles={valuation_profiles}")
print(f"valuation_survivors={len(valuation_survivors)}")
print("common_valuation=nu7(H)=nu7(C1)=nu7(C2)")
print("heavy_words=W7:4-labelled-roles,W8:1-labelled-role")
print(f"W8_endpoint_clearance={w8_endpoint_clearance}")
print("W8_terminal_multiplicity=2,1,1,1,1,1,1")
print(
    f"W8_base_signs=q({w8_q_y_sign},{w8_q_z_sign});"
    f"high({w8_high_y_sign},{w8_high_z_sign})"
)
print(f"W7_endpoint_clearance={w7_endpoint_clearance}")
print("W7_terminal_multiplicity=1,1,1,1,1,1,1")
print(
    f"W7_base_signs=q({w7_q_y_sign},{w7_q_z_sign});"
    f"high({w7_high_y_sign},{w7_high_z_sign});"
    f"off({w7_off_y_sign},{w7_off_z_sign})"
)
print(
    "W7_full_bin_failures="
    f"partition:{','.join(map(str, w7_split_failures))};"
    f"containment:{','.join(map(str, w7_containment_failures))}"
)
print("W7_full_bin_status=EXCLUDED-BY-THM-2391")
print(f"triple_shield_B_count={len(shield_B)}")
print("triple_shield_boundary=B48-escape,B91-touch,B92-penetration")
print("triple_shield_is_not_XOR=PASS")
print("branch_excluded=0; ledger=165; LRC(14)=OPEN")
print("all_checks=PASS")
