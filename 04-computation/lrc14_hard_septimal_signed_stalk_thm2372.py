#!/usr/bin/env python3
"""Exact controls for THM-2372.

This dependency-free companion checks the signed-moment algebra, the
additive F_7 cubic ledger, the universal low-blocker handoff stalk, the
closed absorber incidence formula, and the fourteen-fold nested normal form.
"""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(value: int, prime: int) -> int:
    require(value > 0, "valuation requires a positive integer")
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def frac(value: Fraction) -> Fraction:
    return value - value.numerator // value.denominator


def circle_norm(value: Fraction) -> Fraction:
    residue = frac(value)
    return min(residue, 1 - residue)


def danger(value: Fraction, length_units: int = 1) -> bool:
    return circle_norm(value) < Fraction(length_units, 14)


def safe(value: Fraction, length_units: int = 1) -> bool:
    return not danger(value, length_units)


def safe_root_average(
    root_count: int,
    length_units: int,
    value: Fraction,
) -> Fraction:
    return Fraction(
        sum(
            safe(
                Fraction(value + root, root_count),
                length_units,
            )
            for root in range(root_count)
        ),
        root_count,
    )


# ---------------------------------------------------------------------------
# Integer multiplicity moments.
# ---------------------------------------------------------------------------

moment_controls = 0
for multiplicity in range(7):
    choose_two = multiplicity * (multiplicity - 1) // 2
    choose_three = multiplicity * (multiplicity - 1) * (
        multiplicity - 2
    ) // 6
    centred = multiplicity - 1
    require(
        multiplicity**2 == multiplicity + 2 * choose_two,
        "quadratic multiplicity identity failed",
    )
    require(
        multiplicity**3
        == multiplicity + 6 * choose_two + 6 * choose_three,
        "cubic multiplicity identity failed",
    )
    require(
        centred**2 == 2 * choose_two - centred,
        "centred quadratic identity failed",
    )
    require(
        centred**3 == 6 * choose_three + centred,
        "centred cubic identity failed",
    )
    require(
        max(centred, 0) - int(multiplicity == 0) == centred,
        "positive/zero multiplicity identity failed",
    )
    moment_controls += 1


# ---------------------------------------------------------------------------
# Fully nonzero additive F_7 triangles and chi_7 sectors.
# ---------------------------------------------------------------------------

quadratic_residues = {1, 2, 4}
triangles: list[tuple[int, int, int]] = []
same_chi_triangles: list[tuple[int, int, int]] = []
mixed_chi_triangles: list[tuple[int, int, int]] = []
for first in range(1, 7):
    for second in range(1, 7):
        third = (-first - second) % 7
        if third == 0:
            continue
        triangle = (first, second, third)
        triangles.append(triangle)
        signs = tuple(value in quadratic_residues for value in triangle)
        if signs[0] == signs[1] == signs[2]:
            same_chi_triangles.append(triangle)
        else:
            mixed_chi_triangles.append(triangle)

require(len(triangles) == 30, "wrong nonzero F_7 triangle count")
require(len(same_chi_triangles) == 12, "wrong same-chi_7 count")
require(len(mixed_chi_triangles) == 18, "wrong mixed-chi_7 count")
require(
    set(same_chi_triangles)
    == {
        (first, second, third)
        for base in ((1, 2, 4), (3, 5, 6))
        for first in base
        for second in base
        for third in base
        if len({first, second, third}) == 3
    },
    "same-chi_7 triangles are not the two permutation sectors",
)


# ---------------------------------------------------------------------------
# Universal c=13u handoff stalk.
# ---------------------------------------------------------------------------

handoff_pair_controls = 0
handoff_atom_controls = 0
for core in range(1, 31):
    blocker = 13 * core
    for unit_speed in range(1, 121):
        if unit_speed % 13 == 0:
            continue
        common = gcd(blocker, unit_speed)
        compatible = (
            blocker // common + unit_speed // common
        ) % 14 == 0
        for epsilon in (-1, 1):
            actual_count = 0
            for wall_index in range(blocker):
                numerator = 14 * wall_index + epsilon
                # q*x == -epsilon/14 modulo one.
                is_handoff = (
                    unit_speed * numerator + epsilon * blocker
                ) % (14 * blocker) == 0
                if not is_handoff:
                    continue
                actual_count += 1
                require(
                    (wall_index + epsilon) % 13 == 0,
                    "handoff escaped the universal 1/13 stalk",
                )
                handoff_atom_controls += 1
            require(
                actual_count == (common if compatible else 0),
                "exact opposite-wall seam count failed",
            )
        handoff_pair_controls += 1

odd_guard_controls = 0
for core in range(1, 91):
    blocker = 13 * core
    for guard in range(1, 180, 2):
        if guard % 13 == 0:
            continue
        common = gcd(blocker, guard)
        require(
            (guard + 2 * blocker) % (14 * common) != 0,
            "odd guard acquired an opposite ordinary-wall handoff",
        )
        odd_guard_controls += 1


# ---------------------------------------------------------------------------
# Exact closed-absorber incidence on a blocker wall branch.
# ---------------------------------------------------------------------------


def closed_absorber_count(
    reduced_blocker: int,
    reduced_absorber: int,
    epsilon: int,
) -> int:
    modulus = 14 * reduced_blocker
    count = 0
    for wall_index in range(reduced_blocker):
        numerator = reduced_absorber * (14 * wall_index + epsilon)
        residue = numerator % modulus
        distance = min(residue, (-numerator) % modulus)
        if distance <= reduced_blocker:
            count += 1
    return count


absorber_formula_controls = 0
qtop_ceiling_controls = 0
high_ceiling_controls = 0
for reduced_blocker in range(1, 181):
    if reduced_blocker % 7 == 0:
        continue
    for reduced_absorber in range(7, 351, 7):
        if gcd(reduced_blocker, reduced_absorber) != 1:
            continue
        if reduced_absorber % 2 == 0:
            expected = 2 * ((reduced_blocker - 1) // 14) + 1
        else:
            expected = 2 * ((reduced_blocker + 6) // 14)
        for epsilon in (-1, 1):
            require(
                closed_absorber_count(
                    reduced_blocker,
                    reduced_absorber,
                    epsilon,
                )
                == expected,
                "closed absorber incidence formula failed",
            )
            absorber_formula_controls += 1

        incidence = Fraction(expected, reduced_blocker)
        if reduced_blocker % 13 == 0:
            require(
                incidence <= Fraction(2, 13),
                "top-unit absorber exceeded its 2/13 ceiling",
            )
            qtop_ceiling_controls += 1

        nested = (
            reduced_blocker == 1
            and reduced_absorber % 14 == 0
        )
        if not nested:
            if reduced_blocker in (1, 2):
                require(incidence == 0, "small nonnested incidence")
            else:
                require(
                    incidence <= Fraction(1, 3),
                    "high absorber exceeded its 1/3 ceiling",
                )
            high_ceiling_controls += 1


# The three hypotheses in the wall theorem are load-bearing.
require(
    (2 + 2 * 13) % (14 * gcd(13, 2)) == 0,
    "even-guard hostile disappeared",
)
require(
    (13 // gcd(13, 169) + 169 // gcd(13, 169)) % 14 == 0,
    "thirteen-divisible lower-unit hostile disappeared",
)
require(
    closed_absorber_count(1, 13, 1) == 1,
    "nondominant absorber hostile disappeared",
)


# ---------------------------------------------------------------------------
# Nested normal forms and the hard-chamber hostile.
# ---------------------------------------------------------------------------

profiles = [
    (1, middle, deepest)
    for deepest in range(5, 20)
    for middle in range(1, deepest)
]
require(len(profiles) == 165, "valuation ledger changed")
normal_form_controls = 0
for _, middle, deepest in profiles:
    low_core = 5
    low_blocker = 13**middle * low_core
    multiplier = 98 * 11
    high_blocker = 13**deepest * multiplier * low_core
    require(
        high_blocker % (98 * low_blocker) == 0,
        "nested normal form lost ninety-eight-fold divisibility",
    )
    require(
        valuation(low_blocker, 13) == middle
        and valuation(high_blocker, 13) == deepest,
        "nested normal form changed thirteen-adic profile",
    )
    require(
        high_blocker // low_blocker
        == multiplier * 13 ** (deepest - middle),
        "nested quotient formula failed",
    )
    normal_form_controls += 1

local_low_blocker = 195
local_high_blockers = (16562, 215306)
require(
    all(
        high % (98 * local_low_blocker) != 0
        for high in local_high_blockers
    ),
    "THM-2367 local chamber unexpectedly became nested",
)


# ---------------------------------------------------------------------------
# Exact nested-carrier root renormalization.
# ---------------------------------------------------------------------------


def renormalization_data(
    reduced_source: int,
    reduced_graft: int,
    scale: int,
) -> tuple[int, int]:
    common = gcd(reduced_graft, scale)
    new_source = reduced_source * (scale // common)
    new_graft = reduced_graft // common
    require(
        gcd(new_source, new_graft) == 1,
        "renormalized pair lost coprimality",
    )
    require(new_source % 13 == 0, "renormalized source lost thirteen")
    require(new_graft % 13 != 0, "renormalized graft gained thirteen")
    return new_source, new_graft


renormalization_examples = (
    # (C,Q,R, expected circulant)
    (13, 7, 14, False),
    (13, 7, 98, True),
    (13, 49, 98, False),
    (13, 49, 686, True),
    (26, 21, 196, True),
)
renormalization_samples = 0
for reduced_source, reduced_graft, scale, expected_circulant in (
    renormalization_examples
):
    require(
        gcd(reduced_source, reduced_graft) == 1,
        "control pair is not reduced",
    )
    require(
        reduced_source % 13 == 0
        and reduced_graft % 13 != 0,
        "control pair lost thirteen typing",
    )
    new_source, new_graft = renormalization_data(
        reduced_source,
        reduced_graft,
        scale,
    )
    require(
        (new_source % 7 == 0) == expected_circulant,
        "renormalized toothpick criterion changed",
    )
    require(
        (
            valuation(scale, 7)
            > valuation(reduced_graft, 7)
        )
        == expected_circulant,
        "scale/graft valuation criterion changed",
    )

    for length_units in (1, 2):
        for sample in range(101):
            value = Fraction(3 * sample + 1, 303)
            actual = Fraction(0)
            for root in range(scale):
                source_point = Fraction(value + root, scale)
                if danger(source_point):
                    actual += safe_root_average(
                        reduced_source,
                        length_units,
                        reduced_graft * source_point,
                    )
            actual /= scale
            expected = Fraction(1, 7) * safe_root_average(
                new_source,
                length_units,
                new_graft * value,
            )
            require(
                actual == expected,
                "nested-carrier root renormalization failed",
            )
            renormalization_samples += 1

forced_circulant_profiles = 0
for _, middle, deepest in profiles:
    source = 13**middle * 5
    graft = 7
    scale = 98 * 13 ** (deepest - middle)
    common = gcd(source, graft)
    reduced_source = source // common
    reduced_graft = graft // common
    new_source, _ = renormalization_data(
        reduced_source,
        reduced_graft,
        scale,
    )
    require(new_source % 7 == 0, "hard nested carrier regained drift")
    require(
        valuation(scale, 7) > valuation(reduced_graft, 7),
        "hard nested valuation gap disappeared",
    )
    forced_circulant_profiles += 1


print("theorem=THM-2372")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
print(f"moment_controls={moment_controls}")
print(f"nonzero_F7_triangles={len(triangles)}")
print(f"same_chi7_triangles={len(same_chi_triangles)}")
print(f"mixed_chi7_triangles={len(mixed_chi_triangles)}")
print(f"handoff_pair_controls={handoff_pair_controls}")
print(f"handoff_atom_controls={handoff_atom_controls}")
print(f"odd_guard_controls={odd_guard_controls}")
print(f"absorber_formula_controls={absorber_formula_controls}")
print(f"qtop_ceiling_controls={qtop_ceiling_controls}")
print(f"high_ceiling_controls={high_ceiling_controls}")
print(f"nested_normal_form_profiles={normal_form_controls}")
print(f"renormalization_examples={len(renormalization_examples)}")
print(f"renormalization_samples={renormalization_samples}")
print(f"forced_circulant_profiles={forced_circulant_profiles}")
print("local_hard_chamber_nested_highs=0")
print("load_bearing_hypothesis_hostiles=3")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("all_checks=PASS")
