#!/usr/bin/env python3
"""Exact referee for THM-2275.

The script uses only integer arithmetic and fractions.  It checks:

* the degree-19/20 mixed Selberg tensor boundary;
* the N=231/232 symmetric mixed-crossing boundary;
* the exact 7-unit Fourier-support condition through height 462; and
* the primitive THM-2266 dependency-ray classification.
"""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ALPHA_GUARD = Fraction(5, 7)
ALPHA_SAFE = Fraction(6, 7)
DELTA_FIVE = Fraction(961, 6930)
DEEP_SAFE_FLOOR = Fraction(55, 91)
PAIR_PRODUCT_MAX = 757


def tensor_budget(degree: int) -> Fraction:
    epsilon = Fraction(1, degree + 1)
    return (
        2 * epsilon / (ALPHA_GUARD + epsilon)
        + 8 * 2 * epsilon / (ALPHA_SAFE + epsilon)
    )


def tensor_constant(degree: int) -> Fraction:
    epsilon = Fraction(1, degree + 1)
    return (
        (ALPHA_GUARD + epsilon)
        * (ALPHA_SAFE + epsilon) ** 8
        * (1 - tensor_budget(degree))
    )


def crossing_data(bandwidth: int) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    unit_lower = DELTA_FIVE - Fraction(9, bandwidth)
    deep_lower = DEEP_SAFE_FLOOR - Fraction(9, 2 * bandwidth)
    joint_error = Fraction(27, 2 * bandwidth)
    margin = unit_lower * deep_lower - joint_error
    return unit_lower, deep_lower, joint_error, margin


def coprime_atlas(require_thirteen_units: bool) -> tuple[tuple[int, int], ...]:
    return tuple(
        (left, right)
        for left in range(1, PAIR_PRODUCT_MAX + 1)
        for right in range(1, PAIR_PRODUCT_MAX // left + 1)
        if gcd(left, right) == 1
        and (
            not require_thirteen_units
            or (left % 13 != 0 and right % 13 != 0)
        )
    )


def main() -> None:
    budget_19 = tensor_budget(19)
    budget_20 = tensor_budget(20)
    constant_20 = tensor_constant(20)

    require(budget_19 == Fraction(13762, 13589), "degree-19 budget drift")
    require(1 - budget_19 == Fraction(-173, 13589), "degree-19 bracket drift")
    require(budget_20 == Fraction(147, 152), "degree-20 budget drift")
    require(1 - budget_20 == Fraction(5, 152), "degree-20 bracket drift")
    require(
        constant_20 == Fraction(8938717390, 794280046581),
        "degree-20 tensor constant drift",
    )
    require(constant_20 > 0, "degree-20 tensor must be positive")

    # Every coordinate budget is decreasing in its degree.  Hence the
    # all-degree-19 row is the best row with maximum degree at most 19.
    require(budget_19 > 1, "degree-19 certificate should fail")
    require(budget_20 < 1, "degree-20 certificate should pass")

    unit_231, deep_231, error_231, margin_231 = crossing_data(231)
    unit_232, deep_232, error_232, margin_232 = crossing_data(232)
    require(unit_231 == Fraction(691, 6930), "N=231 unit lower drift")
    require(deep_231 == Fraction(1171, 2002), "N=231 deep lower drift")
    require(error_231 == Fraction(9, 154), "N=231 error drift")
    require(margin_231 == Fraction(-1649, 13873860), "N=231 margin drift")
    require(unit_232 == Fraction(80291, 803880), "N=232 unit lower drift")
    require(deep_232 == Fraction(24701, 42224), "N=232 deep lower drift")
    require(error_232 == Fraction(27, 464), "N=232 error drift")
    require(
        margin_232 == Fraction(8134831, 33943029120),
        "N=232 margin drift",
    )
    require(margin_231 < 0 < margin_232, "crossing boundary sign drift")

    first_valid = None
    for bandwidth in range(2, 233):
        unit_lower, deep_lower, _, margin = crossing_data(bandwidth)
        if unit_lower > 0 and deep_lower > 0 and margin > 0:
            first_valid = bandwidth
            break
    require(first_valid == 232, "first symmetric bandwidth drift")

    height = 2 * first_valid - 2
    require(height == 462, "crossing height drift")
    for interval_numerator in (5, 6):
        for mode in range(-height, height + 1):
            if mode == 0:
                continue
            interval_coefficient_nonzero = (interval_numerator * mode) % 7 != 0
            require(
                interval_coefficient_nonzero == (mode % 7 != 0),
                "7-unit Fourier support drift",
            )
    positive_seven_unit_modes = sum(
        1 for mode in range(1, height + 1) if mode % 7 != 0
    )
    require(positive_seven_unit_modes == 396, "7-unit mode census drift")

    generic_atlas = coprime_atlas(False)
    unit_atlas = coprime_atlas(True)
    require(len(generic_atlas) == 3643, "generic THM-2266 atlas drift")
    require(
        len({tuple(sorted(pair)) for pair in generic_atlas}) == 1822,
        "generic unordered THM-2266 atlas drift",
    )
    require(len(unit_atlas) == 3279, "13-unit oriented atlas drift")
    require(
        len({tuple(sorted(pair)) for pair in unit_atlas}) == 1640,
        "13-unit unordered atlas drift",
    )

    for left, right in unit_atlas:
        require(gcd(13 * right, left) == 1, "H-owner lift not primitive")
        require(gcd(right, 13 * left) == 1, "owner-q lift not primitive")

    h_dependency_before_parity = tuple(
        (left, right)
        for left, right in unit_atlas
        if max(13 * right, left) <= 20
    )
    h_dependency = tuple(
        (left, right)
        for left, right in h_dependency_before_parity
        if left % 2 == 1
    )
    q_dependency = tuple(
        (left, right)
        for left, right in unit_atlas
        if max(right, 13 * left) <= 20
    )

    expected_h_dependency = tuple(
        (value, 1) for value in (1, 3, 5, 7, 9, 11, 15, 17, 19)
    )
    expected_q_dependency = tuple(
        (1, value)
        for value in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20)
    )
    require(
        h_dependency == expected_h_dependency,
        "odd guard-owner dependency rays drift",
    )
    require(
        q_dependency == expected_q_dependency,
        "owner-q dependency rays drift",
    )
    require(len(h_dependency) == 9, "odd H-ray count drift")
    require(len(q_dependency) == 19, "q-ray count drift")

    max_h_owner_lift = max(max(13 * right, left) for left, right in unit_atlas)
    max_owner_q_lift = max(max(right, 13 * left) for left, right in unit_atlas)
    require(max_h_owner_lift == 9841, "H-owner lifted height drift")
    require(max_owner_q_lift == 9841, "owner-q lifted height drift")

    print("THM-2275 exact referee")
    print(
        "mixed_tensor_degree_19 "
        f"budget={budget_19} bracket={1 - budget_19}"
    )
    print(
        "mixed_tensor_degree_20 "
        f"budget={budget_20} bracket={1 - budget_20} "
        f"constant={constant_20}"
    )
    print(
        "crossing_N_231 "
        f"unit_lower={unit_231} deep_lower={deep_231} "
        f"joint_error={error_231} margin={margin_231}"
    )
    print(
        "crossing_N_232 "
        f"unit_lower={unit_232} deep_lower={deep_232} "
        f"joint_error={error_232} margin={margin_232}"
    )
    print(
        f"first_symmetric_bandwidth={first_valid} "
        f"scalar_crossing_height={height} "
        f"positive_7_unit_modes={positive_seven_unit_modes}"
    )
    print(
        f"generic_pair_atlas={len(generic_atlas)} "
        f"generic_unordered={len({tuple(sorted(pair)) for pair in generic_atlas})}"
    )
    print(
        f"thirteen_unit_pair_atlas={len(unit_atlas)} "
        f"thirteen_unit_unordered={len({tuple(sorted(pair)) for pair in unit_atlas})}"
    )
    print(f"odd_H_dependency_rays={tuple(left for left, _ in h_dependency)}")
    print(f"q_dependency_rays={tuple(right for _, right in q_dependency)}")
    print(
        f"max_scalar_pair_lift_height="
        f"{max(max_h_owner_lift, max_owner_q_lift)}"
    )
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
