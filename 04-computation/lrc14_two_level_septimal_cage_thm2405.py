#!/usr/bin/env python3
"""Exact companion for THM-2405.

The script independently reconstructs the two septimal one-sheet
identities, the high-band cap, the complete THM-2392 b=2 compatible
union table, and the resulting residual ratio/floor table.  All
measure calculations use exact rational interval arrangements.
"""

from fractions import Fraction as F
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac_part(x: F) -> F:
    return x - x.numerator // x.denominator


def circle_norm(x: F) -> F:
    residue = frac_part(x)
    return min(residue, 1 - residue)


def danger(speed: int, x: F, width: int = 1) -> bool:
    return circle_norm(speed * x) < F(width, 14)


def danger_boundaries(speed: int, width: int = 1) -> set[F]:
    denominator = 14 * speed
    return {
        F((14 * k + sign * width) % denominator, denominator)
        for k in range(speed)
        for sign in (-1, 1)
    }


def interval_measure(
    predicate,
    speed_widths: tuple[tuple[int, int], ...],
) -> F:
    points = {F(0), F(1)}
    for speed, width in speed_widths:
        points.update(danger_boundaries(speed, width))
    ordered = sorted(points)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        if left < right and predicate((left + right) / 2):
            total += right - left
    return total


def set_measure(speed: int, width: int = 1) -> F:
    return interval_measure(
        lambda x: danger(speed, x, width),
        ((speed, width),),
    )


# ---------------------------------------------------------------------------
# 1. The strict-open seven-sheet count.
# ---------------------------------------------------------------------------

sheet_checks = 0
endpoint_skips = 0
for speed in range(1, 43):
    if speed % 7 == 0:
        continue
    for denominator in (101, 103):
        for numerator in range(denominator):
            x = F(numerator, denominator)
            norms = tuple(
                circle_norm(speed * (x + F(sheet, 7)))
                for sheet in range(7)
            )
            if F(1, 14) in norms:
                endpoint_skips += 1
                continue
            require(
                sum(value < F(1, 14) for value in norms) == 1,
                "a seven-unit did not occupy exactly one strict-open sheet",
            )
            sheet_checks += 1


# ---------------------------------------------------------------------------
# 2. Exact high-band and two-level independence controls.
# ---------------------------------------------------------------------------

require(set_measure(1) == F(1, 7), "ordinary danger mass")
same_line = interval_measure(
    lambda x: danger(1, x) and danger(13, x),
    ((1, 1), (13, 1)),
)
require(same_line == F(1, 91), "same-line 1:13 overlap")
high_band_mass = F(1, 7) - same_line
require(high_band_mass == F(12, 91), "wrong high-band mass")

high_band_cases = (
    # (septimal depth M, q_star, C_3)
    (1, 7, 49),
    (1, 21, 98),
    (2, 49, 343),
    (2, 98, 686),
)

low_units = (1, 2, 3, 5, 6, 13, 27, 29, 39, 52, 156, 676)
high_band_intersection_checks = 0

for depth, q_star, C3 in high_band_cases:
    c3 = 13 * C3
    require(q_star % (7**depth) == 0, "q-star depth lower bound")
    require(q_star % (7 ** (depth + 1)) != 0, "q-star depth exactness")
    require(C3 % (7 ** (depth + 1)) == 0, "C3 is not above q-star")

    high_safe = lambda x, C3=C3, c3=c3, q_star=q_star: (
        danger(C3, x)
        and not danger(c3, x)
        and not danger(q_star, x)
    )
    high_safe_specs = ((C3, 1), (c3, 1), (q_star, 1))
    high_safe_mass = interval_measure(high_safe, high_safe_specs)
    require(high_safe_mass == F(72, 637), "wrong 72/637 high-safe band")

    for low in low_units:
        if low % 7 == 0:
            continue
        intersection = interval_measure(
            lambda x, low=low: high_safe(x) and danger(low, x),
            (*high_safe_specs, (low, 1)),
        )
        require(
            intersection == F(72, 4459),
            "one primitive low sheet did not take 1/7 of the band",
        )
        high_band_intersection_checks += 1


# The two-sheet cap is attained at the level of the abstract hypotheses.
C3_sharp = 49
c3_sharp = 13 * C3_sharp
q_sharp = 7
v_sharp, w_sharp = 27, 29
sharp_union = interval_measure(
    lambda x: (
        danger(C3_sharp, x)
        and not danger(c3_sharp, x)
        and not danger(q_sharp, x)
        and (danger(v_sharp, x) or danger(w_sharp, x))
    ),
    (
        (C3_sharp, 1),
        (c3_sharp, 1),
        (q_sharp, 1),
        (v_sharp, 1),
        (w_sharp, 1),
    ),
)
require(sharp_union == F(144, 4459), "abstract two-sheet cap not attained")


# Hypothesis hostiles.
nonunit_invariant = interval_measure(
    lambda x: danger(7, x) and danger(7, x),
    ((7, 1),),
)
require(nonunit_invariant == F(1, 7), "nonunit hostile mass")
require(
    nonunit_invariant != F(1, 49),
    "the 7-unit hypothesis was accidentally unnecessary",
)

same_depth_band = interval_measure(
    lambda x: danger(49, x) and not danger(13 * 49, x) and danger(49, x),
    ((49, 1), (13 * 49, 1)),
)
require(same_depth_band == F(12, 91), "same-depth hostile mass")
require(
    same_depth_band != F(12, 637),
    "the strict septimal depth gap was accidentally unnecessary",
)


# ---------------------------------------------------------------------------
# 3. Reconstruct the complete compatible low-cage table.
# ---------------------------------------------------------------------------

expected_unions = {
    (1, 1): F(193, 1183),
    (1, 2): F(114, 1183),
    (2, 1): F(239, 2366),
    (1, 3): F(263, 3549),
    (3, 1): F(95, 1183),
    (4, 1): F(331, 4732),
    (2, 3): F(43, 546),
    (3, 2): F(95, 1183),
    (3, 4): F(491, 7098),
    (4, 3): F(331, 4732),
}

computed_unions: dict[tuple[int, int], F] = {}
for (a, b), expected in expected_unions.items():
    require(gcd(a, b) == 1 and gcd(a * b, 91) == 1, "invalid ratio row")
    C1, C2, c1, c2 = b, 13 * a, 13 * b, 169 * a
    union = interval_measure(
        lambda x, C1=C1, C2=C2, c1=c1, c2=c2: (
            (danger(C1, x) or danger(C2, x))
            and (danger(c1, x) or danger(c2, x))
        ),
        ((C1, 1), (C2, 1), (c1, 1), (c2, 1)),
    )
    require(union == expected, "compatible low-cage table mismatch")
    computed_unions[(a, b)] = union


# ---------------------------------------------------------------------------
# 4. Exact floor and residual tables.
# ---------------------------------------------------------------------------

hole_mass = F(36, 343)
high_cap = F(144, 4459)
required_low_mass = hole_mass - high_cap
require(required_low_mass == F(324, 4459), "wrong residual low threshold")

floors = {
    ratio: required_low_mass - union
    for ratio, union in computed_unions.items()
    if required_low_mass > union
}
expected_floors = {
    (4, 1): F(629, 231868),
    (3, 4): F(1213, 347802),
    (4, 3): F(629, 231868),
}
require(floors == expected_floors, "wrong positive clean-mass floors")

surviving = tuple(
    ratio for ratio in expected_unions if ratio not in expected_floors
)
require(
    surviving
    == ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)),
    "wrong seven-ratio residual",
)

charged_floors = {ratio: floor / 52 for ratio, floor in floors.items()}
tensor_floors = {ratio: floor / 338 for ratio, floor in floors.items()}
require(
    charged_floors
    == {
        (4, 1): F(629, 12057136),
        (3, 4): F(1213, 18085704),
        (4, 3): F(629, 12057136),
    },
    "wrong charged-cell floors",
)
require(
    tensor_floors
    == {
        (4, 1): F(629, 78371384),
        (3, 4): F(1213, 117557076),
        (4, 3): F(629, 78371384),
    },
    "wrong C7xC13 tensor floors",
)

# THM-2392 leaves only (4,3) at M=2, and it is now in the positive table.
m2_bank = ((4, 3),)
require(all(ratio in floors for ratio in m2_bank), "M=2 was not eliminated")


print("theorem=THM-2405")
print("status=PROVED-CANDIDATE+VERIFIED-EXACT;INDEPENDENT-AUDIT-PENDING")
print(f"seven_sheet_checks={sheet_checks};endpoint_skips={endpoint_skips}")
print(f"same_line={same_line};high_band={high_band_mass}")
print("high_safe_band=72/637")
print(f"high_band_intersection_checks={high_band_intersection_checks}")
print(f"one_low_sheet=72/4459;two_low_cap={high_cap}")
print(f"abstract_cap_sharp={sharp_union};witness=27,29|q=7,C3=49")
print("hostiles=nonunit-low:1/7!=1/49;same-depth-q:12/91!=12/637")
print(f"compatible_union_rows={len(computed_unions)}")
print(f"required_low_mass={required_low_mass}")
print("eliminated_oriented_ratios=4:1,3:4,4:3")
print("surviving_oriented_ratios=1:1,1:2,2:1,1:3,3:1,2:3,3:2")
print("clean_floors=4:1=629/231868,3:4=1213/347802,4:3=629/231868")
print(
    "charged_floors="
    "4:1=629/12057136,3:4=1213/18085704,4:3=629/12057136"
)
print(
    "tensor_floors="
    "4:1=629/78371384,3:4=1213/117557076,4:3=629/78371384"
)
print("M2_no_clean=EMPTY;remaining_no_clean=M1,seven-oriented-ratios")
print("endpoint_scope=a.e.-strict-open")
print("row_decrement=0;canonical_target_landing=0;LRC(14)=OPEN")
print("all_checks=PASS")
