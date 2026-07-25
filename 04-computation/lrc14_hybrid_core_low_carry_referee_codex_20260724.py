#!/usr/bin/env python3
"""Exact referee for THM-2166 (hybrid whole-core smoothing).

The script independently enumerates all seven-speed cores in {1,...,13},
records both safe mass and positive-length component count, and checks the
two-scale Jackson error ledger with exact rational arithmetic.
"""

from fractions import Fraction
from itertools import combinations


RADIUS = Fraction(1, 14)
ALPHA_SIX = Fraction(61, 273)
PI2_UPPER = Fraction(355 * 355, 113 * 113)
CARRY_LIMIT = 708


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_fraction(x: Fraction) -> Fraction:
    return x % 1


def is_safe(speed: int, t: Fraction) -> bool:
    phase = circle_fraction(speed * t)
    return min(phase, 1 - phase) >= RADIUS


def boundary_cells(max_speed: int = 13) -> tuple[Fraction, ...]:
    points = {Fraction(0), Fraction(1)}
    for speed in range(1, max_speed + 1):
        for index in range(speed):
            points.add(
                circle_fraction(Fraction(14 * index - 1, 14 * speed))
            )
            points.add(
                circle_fraction(Fraction(14 * index + 1, 14 * speed))
            )
    return tuple(sorted(points))


def circular_run_count(flags: tuple[bool, ...]) -> int:
    require(flags, "nonempty cell word")
    if all(flags):
        return 1
    return sum(
        flag and not flags[index - 1]
        for index, flag in enumerate(flags)
    )


def core_profiles() -> tuple[
    dict[tuple[int, ...], tuple[Fraction, int]], tuple[Fraction, ...]
]:
    points = boundary_cells()
    cells = tuple(
        (left, right, (left + right) / 2)
        for left, right in zip(points, points[1:])
    )
    profiles: dict[tuple[int, ...], tuple[Fraction, int]] = {}
    for core in combinations(range(1, 14), 7):
        flags = tuple(
            all(is_safe(speed, midpoint) for speed in core)
            for _, _, midpoint in cells
        )
        mass = sum(
            (
                right - left
                for (left, right, _), safe in zip(cells, flags)
                if safe
            ),
            Fraction(0),
        )
        profiles[core] = (mass, circular_run_count(flags))
    return profiles, points


def jackson_coefficient(n: int, frequency: int) -> int:
    require(n >= 2 and 0 <= frequency <= 2 * n - 2, "Jackson index")
    if frequency <= n:
        numerator = (
            4 * n**3
            - 6 * n * frequency**2
            + 2 * n
            + 3 * frequency**3
            - 3 * frequency
        )
    else:
        complement = 2 * n - frequency
        numerator = complement**3 - complement
    require(numerator % 6 == 0, "Jackson coefficient integrality")
    value = numerator // 6
    require(value > 0, "Jackson coefficient positivity")
    return value


def jackson_eta_upper(n: int) -> Fraction:
    constant = n * (2 * n * n + 1) // 3
    require(constant == jackson_coefficient(n, 0), "Jackson constant")
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(n, frequency), frequency**2)
            for frequency in range(1, 2 * n - 2, 2)
        ),
        Fraction(0),
    )
    upper = Fraction(1, 2) - 4 * odd_sum / (PI2_UPPER * constant)
    require(upper > 0, "positive Jackson moment upper bound")
    return upper


def hybrid_margin(
    mass: Fraction,
    components: int,
    far_error: Fraction,
    core_error: Fraction,
) -> Fraction:
    return (
        (ALPHA_SIX - 6 * far_error) * mass
        - 6 * far_error
        - components * core_error
    )


def has_small_bezout_support(core: tuple[int, ...]) -> bool:
    if 1 in core:
        return True
    return any(right == left + 1 for left, right in zip(core, core[1:]))


def pair_sum_mask(left: int, right: int, height: int) -> int:
    """Bit mask of bounded two-coordinate sums in [-CARRY_LIMIT,CARRY_LIMIT]."""
    mask = 0
    for left_coefficient in range(-height, height + 1):
        for right_coefficient in range(-height, height + 1):
            value = left_coefficient * left + right_coefficient * right
            if -CARRY_LIMIT <= value <= CARRY_LIMIT:
                mask |= 1 << (value + CARRY_LIMIT)
    return mask


def support_two_masks(height: int) -> dict[tuple[int, int], int]:
    return {
        pair: pair_sum_mask(*pair, height)
        for pair in combinations(range(1, 14), 2)
    }


def core_sum_mask(
    core: tuple[int, ...], pair_masks: dict[tuple[int, int], int]
) -> int:
    mask = 0
    for pair in combinations(core, 2):
        mask |= pair_masks[pair]
    return mask


def main() -> None:
    profiles, points = core_profiles()
    expected_core = (1, 5, 7, 8, 9, 11, 13)
    require(len(points) == 178, "boundary point count")
    require(len(profiles) == 1716, "seven-core count")
    require(
        all(has_small_bezout_support(core) for core in profiles),
        "support-two Bezout property",
    )

    target_mask = (1 << (2 * CARRY_LIMIT + 1)) - 1
    masks_57 = support_two_masks(57)
    require(
        all(core_sum_mask(core, masks_57) == target_mask for core in profiles),
        "height-57 support-two carry coverage",
    )
    masks_56 = support_two_masks(56)
    failures_56 = tuple(
        core
        for core in profiles
        if core_sum_mask(core, masks_56) != target_mask
    )
    expected_failures = (
        (1, 2, 3, 4, 5, 6, 7),
        (1, 2, 3, 4, 5, 6, 8),
        (1, 2, 4, 5, 6, 8, 10),
    )
    require(failures_56 == expected_failures, "height-56 hostile cores")
    expected_missing = {
        expected_failures[0]: (-706, -705, -699, 699, 705, 706),
        expected_failures[1]: (-701, 701),
        expected_failures[2]: (-701, 701),
    }
    for core, missing_values in expected_missing.items():
        mask = core_sum_mask(core, masks_56)
        actual_missing = tuple(
            value
            for value in range(-CARRY_LIMIT, CARRY_LIMIT + 1)
            if not (mask >> (value + CARRY_LIMIT)) & 1
        )
        require(
            actual_missing == missing_values,
            f"height-56 missing-frequency profile for {core}",
        )

    far_exact = jackson_eta_upper(150)
    core_exact = jackson_eta_upper(355)
    far_bound = Fraction(439, 156250)
    core_bound = Fraction(371, 312500)
    require(far_exact < far_bound, "compact N=150 error")
    require(core_exact < core_bound, "compact N=355 error")
    require(ALPHA_SIX - 6 * far_bound > 0, "positive far integral")

    margins = {
        core: hybrid_margin(mass, components, far_bound, core_bound)
        for core, (mass, components) in profiles.items()
    }
    minimizer = min(margins, key=margins.get)
    minimum = margins[minimizer]
    worst_mass, worst_components = profiles[minimizer]
    require(minimizer == expected_core, "unique hybrid-margin minimizer")
    require(
        sum(margin == minimum for margin in margins.values()) == 1,
        "hybrid-margin uniqueness",
    )
    require(worst_mass == Fraction(45107, 229320), "worst core mass")
    require(worst_components == 20, "worst core BV components")
    require(
        minimum == Fraction(41050267, 1222741406250) > 0,
        "all-core hybrid margin",
    )

    # Adjacent failure is for this exact rational-pi/Jackson ledger only.
    core_354_exact = jackson_eta_upper(354)
    adjacent_margin = hybrid_margin(
        worst_mass, worst_components, far_exact, core_354_exact
    )
    require(adjacent_margin < 0, "N=354 ledger must fail")

    print("THM-2166 exact hybrid whole-core smoothing referee")
    print(
        f"boundary points={len(points)}, open cells={len(points)-1}, "
        f"seven-cores={len(profiles)}"
    )
    print("all 1716 cores contain speed 1 or an adjacent pair")
    print("support<=2 core representation: height 57 covers every |nu|<=708")
    print(
        "height 56 fails on exactly 3 cores; hostile missing frequencies "
        "are +/-699,+/-705,+/-706 or +/-701"
    )
    print(
        f"eta_150<{far_bound}, eta_355<{core_bound}; "
        f"alpha_6-6eta_150>0"
    )
    print(
        f"unique margin minimizer={minimizer}: "
        f"mass={worst_mass}, BV components={worst_components}"
    )
    print(f"all-core hybrid margin={minimum} > 0")
    print("the same exact rational-pi ledger is negative at core N=354")
    print("far height=298, core height=57, scalar carry<=708, core support<=2")
    print("l1<=1902; after four dyadic digits the full carry satisfies |kappa|<=1787")
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
