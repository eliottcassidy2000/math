#!/usr/bin/env python3
"""Exact referee for THM-2275.

The script uses only integer arithmetic and fractions.  It checks:

* the degree-19/20 mixed Selberg tensor boundary;
* the N=231/232 symmetric mixed-crossing boundary;
* the N=354/355 and N=1058/1059 adaptive-cut boundaries;
* the exact 7-unit Fourier-support condition through height 462; and
* the primitive THM-2266 dependency rays and all fixed-section lift heights.
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
GUARD_OWNER_PAIR_FLOOR = Fraction(4, 7)
SEVEN_TERMINAL_SAFE_FLOOR = Fraction(15, 154)
OWNER_UNIT_PAIR_FLOOR = Fraction(66, 91)
DELTA_SIX = Fraction(191, 6930)
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


def adaptive_crossing_data(
    bandwidth: int, left_floor: Fraction, right_floor: Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    left_lower = left_floor - Fraction(3, bandwidth)
    right_lower = right_floor - Fraction(21, 2 * bandwidth)
    joint_error = Fraction(27, 2 * bandwidth)
    margin = left_lower * right_lower - joint_error
    return left_lower, right_lower, joint_error, margin


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

    guard_owner_354 = adaptive_crossing_data(
        354, GUARD_OWNER_PAIR_FLOOR, SEVEN_TERMINAL_SAFE_FLOOR
    )
    guard_owner_355 = adaptive_crossing_data(
        355, GUARD_OWNER_PAIR_FLOOR, SEVEN_TERMINAL_SAFE_FLOOR
    )
    require(
        guard_owner_354[3] == Fraction(-3, 15010072),
        "N=354 guard-owner margin drift",
    )
    require(
        guard_owner_355[3] == Fraction(21177, 135854950),
        "N=355 guard-owner margin drift",
    )
    require(
        guard_owner_354[3] < 0 < guard_owner_355[3],
        "guard-owner adaptive boundary sign drift",
    )
    first_guard_owner = None
    for bandwidth in range(2, 356):
        left_lower, right_lower, _, margin = adaptive_crossing_data(
            bandwidth, GUARD_OWNER_PAIR_FLOOR, SEVEN_TERMINAL_SAFE_FLOOR
        )
        if left_lower > 0 and right_lower > 0 and margin > 0:
            first_guard_owner = bandwidth
            break
    require(first_guard_owner == 355, "first guard-owner bandwidth drift")
    guard_owner_height = 2 * first_guard_owner - 2
    require(guard_owner_height == 708, "guard-owner crossing height drift")

    owner_unit_1058 = adaptive_crossing_data(
        1058, OWNER_UNIT_PAIR_FLOOR, DELTA_SIX
    )
    owner_unit_1059 = adaptive_crossing_data(
        1059, OWNER_UNIT_PAIR_FLOOR, DELTA_SIX
    )
    require(
        owner_unit_1058[3] == Fraction(-861505, 47060301288),
        "N=1058 owner-unit margin drift",
    )
    require(
        owner_unit_1059[3] == Fraction(44021, 78582173670),
        "N=1059 owner-unit margin drift",
    )
    require(
        owner_unit_1058[3] < 0 < owner_unit_1059[3],
        "owner-unit adaptive boundary sign drift",
    )
    first_owner_unit = None
    for bandwidth in range(2, 1060):
        left_lower, right_lower, _, margin = adaptive_crossing_data(
            bandwidth, OWNER_UNIT_PAIR_FLOOR, DELTA_SIX
        )
        if left_lower > 0 and right_lower > 0 and margin > 0:
            first_owner_unit = bandwidth
            break
    require(first_owner_unit == 1059, "first owner-unit bandwidth drift")
    owner_unit_height = 2 * first_owner_unit - 2
    require(owner_unit_height == 2116, "owner-unit crossing height drift")

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

    small_original_height = 2 * 20
    natural_crossing_original_height = 2 * height
    guard_owner_crossing_original_height = 2 * guard_owner_height
    owner_unit_crossing_original_height = 2 * owner_unit_height
    guard_owner_pair_original_height = 2 * max_h_owner_lift
    owner_unit_pair_original_height = max_owner_q_lift
    uniform_scalar_rank_height = max(
        max_h_owner_lift,
        max_owner_q_lift,
        guard_owner_height,
        owner_unit_height,
    )
    uniform_original_rank_height = max(
        guard_owner_pair_original_height,
        owner_unit_pair_original_height,
        guard_owner_crossing_original_height,
        owner_unit_crossing_original_height,
    )
    require(small_original_height == 40, "small fixed-section lift drift")
    require(
        natural_crossing_original_height == 924,
        "natural crossing fixed-section lift drift",
    )
    require(
        guard_owner_crossing_original_height == 1416,
        "guard-owner crossing fixed-section lift drift",
    )
    require(
        owner_unit_crossing_original_height == 4232,
        "owner-unit crossing fixed-section lift drift",
    )
    require(
        guard_owner_pair_original_height == 19682,
        "guard-owner pair fixed-section lift drift",
    )
    require(
        owner_unit_pair_original_height == 9841,
        "owner-unit pair fixed-section lift drift",
    )
    require(uniform_scalar_rank_height == 9841, "scalar rank height drift")
    require(
        uniform_original_rank_height == 19682,
        "original rank height drift",
    )

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
        "guard_owner_N_354 "
        f"left_lower={guard_owner_354[0]} "
        f"right_lower={guard_owner_354[1]} "
        f"joint_error={guard_owner_354[2]} "
        f"margin={guard_owner_354[3]}"
    )
    print(
        "guard_owner_N_355 "
        f"left_lower={guard_owner_355[0]} "
        f"right_lower={guard_owner_355[1]} "
        f"joint_error={guard_owner_355[2]} "
        f"margin={guard_owner_355[3]}"
    )
    print(
        f"first_guard_owner_bandwidth={first_guard_owner} "
        f"scalar_crossing_height={guard_owner_height} "
        f"original_crossing_height={guard_owner_crossing_original_height}"
    )
    print(
        "owner_unit_N_1058 "
        f"left_lower={owner_unit_1058[0]} "
        f"right_lower={owner_unit_1058[1]} "
        f"joint_error={owner_unit_1058[2]} "
        f"margin={owner_unit_1058[3]}"
    )
    print(
        "owner_unit_N_1059 "
        f"left_lower={owner_unit_1059[0]} "
        f"right_lower={owner_unit_1059[1]} "
        f"joint_error={owner_unit_1059[2]} "
        f"margin={owner_unit_1059[3]}"
    )
    print(
        f"first_owner_unit_bandwidth={first_owner_unit} "
        f"scalar_crossing_height={owner_unit_height} "
        f"original_crossing_height={owner_unit_crossing_original_height}"
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
    print(
        "fixed_section_heights "
        f"small={small_original_height} "
        f"natural_crossing={natural_crossing_original_height} "
        f"guard_owner_pair={guard_owner_pair_original_height} "
        f"owner_unit_pair={owner_unit_pair_original_height} "
        f"uniform_scalar_rank={uniform_scalar_rank_height} "
        f"uniform_original_rank={uniform_original_rank_height}"
    )
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
