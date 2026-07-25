#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2298.

The script checks three independent parts of the theorem:

* the complete finite local-core enumeration for the three weighted facet
  geometries, using the generalized Irwin--Hall truncated-power formula;
* the two weighted critical-cancellation ledgers; and
* the inherited Machin/Jackson comparison, including the first passing
  bandwidth N=264.

All validity gates use ``require`` so ordinary and optimized Python execute
the same checks.
"""

from fractions import Fraction
from itertools import product
from math import factorial, gcd
import sys


ORDINARY_DANGER = Fraction(1, 7)
GUARD_DANGER = Fraction(2, 7)
ORDINARY_SAFE = Fraction(6, 7)
GUARD_SAFE = Fraction(5, 7)

ORDINARY_DANGER_RADIUS = Fraction(1, 14)
GUARD_DANGER_RADIUS = Fraction(1, 7)
ORDINARY_SAFE_RADIUS = Fraction(3, 7)

SHARP_PI_UPPER_NUMERATOR = 104_348
SHARP_PI_UPPER_DENOMINATOR = 33_215
TOTAL_COORDINATES = 9
TELESCOPE_FACTOR = 2 * TOTAL_COORDINATES
JACKSON_N = 264
SCALAR_HEIGHT = 2 * JACKSON_N - 2
FIXED_SECTION_HEIGHT = 2 * SCALAR_HEIGHT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -floor_fraction(-value)


def truncated_box_cdf(
    endpoint: Fraction,
    widths: tuple[Fraction, ...],
) -> Fraction:
    """Volume of {0<=s_i<=width_i, sum s_i<=endpoint}."""

    dimension = len(widths)
    require(dimension >= 1, "truncated box must have positive dimension")
    total = Fraction(0)
    for mask in range(1 << dimension):
        shift = sum(
            (
                widths[index]
                for index in range(dimension)
                if (mask >> index) & 1
            ),
            Fraction(0),
        )
        remainder = endpoint - shift
        if remainder > 0:
            sign = -1 if mask.bit_count() % 2 else 1
            total += sign * remainder**dimension
    return total / factorial(dimension)


def periodic_strip_volume(
    radii: tuple[Fraction, ...],
    coefficients: tuple[int, ...],
    center: Fraction,
    target_radius: Fraction,
) -> Fraction:
    """Exact box volume where center+sum b_i*y_i is target-close to Z."""

    require(
        len(radii) == len(coefficients) == 4,
        "facet geometry must have four coordinates",
    )
    require(
        Fraction(0) < target_radius < Fraction(1, 2),
        "target strips must be disjoint",
    )
    require(coefficients[0] > 0, "facet coefficient must be nonzero")
    require(
        all(0 <= coefficient <= 6 for coefficient in coefficients),
        "local coefficient outside 0..6",
    )

    active = tuple(
        index
        for index, coefficient in enumerate(coefficients)
        if coefficient
    )
    inactive_volume = Fraction(1)
    for index, coefficient in enumerate(coefficients):
        if not coefficient:
            inactive_volume *= 2 * radii[index]

    widths = tuple(
        2 * coefficients[index] * radii[index]
        for index in active
    )
    jacobian = Fraction(1)
    for index in active:
        jacobian /= coefficients[index]

    offset = center - sum(
        (
            coefficients[index] * radii[index]
            for index in active
        ),
        Fraction(0),
    )
    lower_range = offset
    upper_range = offset + sum(widths, Fraction(0))

    first_integer = floor_fraction(lower_range - target_radius) - 1
    last_integer = ceil_fraction(upper_range + target_radius) + 1
    strip_volume = Fraction(0)
    for integer in range(first_integer, last_integer + 1):
        upper = Fraction(integer) + target_radius - offset
        lower = Fraction(integer) - target_radius - offset
        strip_volume += (
            truncated_box_cdf(upper, widths)
            - truncated_box_cdf(lower, widths)
        )

    volume = inactive_volume * jacobian * strip_volume
    box_volume = product_fraction(2 * radius for radius in radii)
    require(
        Fraction(0) <= volume <= box_volume,
        "periodic strip volume escaped its box",
    )
    return volume


def product_fraction(values) -> Fraction:
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def local_core_minimum(
    facet_radius: Fraction,
    target_multiplier: int,
    maximum_h: int,
) -> tuple[Fraction, tuple[tuple[int, tuple[int, ...]], ...], int]:
    """Enumerate every raw absolute local tuple, caching safe permutations."""

    require(target_multiplier in (1, 2), "unknown target multiplier")
    require(1 <= maximum_h <= 6, "local h bound outside 1..6")
    radii = (
        facet_radius,
        ORDINARY_SAFE_RADIUS,
        ORDINARY_SAFE_RADIUS,
        ORDINARY_SAFE_RADIUS,
    )
    cache: dict[tuple[int, int, tuple[int, int, int]], Fraction] = {}
    minimum: Fraction | None = None
    witnesses: list[tuple[int, tuple[int, ...]]] = []
    raw_cases = 0

    for h in range(1, maximum_h + 1):
        target_radius = Fraction(target_multiplier * h, 14)
        for facet_coefficient in range(1, 7):
            for safe_coefficients in product(range(7), repeat=3):
                common_divisor = h
                for coefficient in (
                    facet_coefficient,
                    *safe_coefficients,
                ):
                    common_divisor = gcd(common_divisor, coefficient)
                if common_divisor != 1:
                    continue

                raw_cases += 1
                sorted_safe = tuple(sorted(safe_coefficients))
                cache_key = (h, facet_coefficient, sorted_safe)
                if cache_key not in cache:
                    coefficients = (facet_coefficient, *sorted_safe)
                    center = Fraction(sum(sorted_safe) % 2, 2)
                    cache[cache_key] = periodic_strip_volume(
                        radii,
                        coefficients,
                        center,
                        target_radius,
                    ) / h
                deficit = cache[cache_key]
                raw_witness = (
                    h,
                    (facet_coefficient, *safe_coefficients),
                )
                if minimum is None or deficit < minimum:
                    minimum = deficit
                    witnesses = [raw_witness]
                elif deficit == minimum:
                    witnesses.append(raw_witness)

    require(minimum is not None, "local enumeration was empty")
    return minimum, tuple(witnesses), raw_cases


def alternating_arctangent_partial(
    denominator: int,
    last_index: int,
) -> Fraction:
    require(denominator > 1, "arctangent denominator must exceed one")
    require(last_index >= 0, "arctangent index must be nonnegative")
    terms = tuple(
        Fraction(
            1,
            (2 * index + 1) * denominator ** (2 * index + 1),
        )
        for index in range(last_index + 1)
    )
    require(
        all(left > right for left, right in zip(terms, terms[1:])),
        "arctangent terms must strictly decrease",
    )
    return sum(
        (
            term if index % 2 == 0 else -term
            for index, term in enumerate(terms)
        ),
        Fraction(0),
    )


def jackson_coefficient(bandwidth: int, frequency: int) -> int:
    require(bandwidth >= 2, "Jackson bandwidth must be at least two")
    require(
        0 <= frequency <= 2 * bandwidth - 2,
        "Jackson frequency outside support",
    )
    if frequency <= bandwidth:
        return (
            4 * bandwidth**3
            - 6 * bandwidth * frequency**2
            + 2 * bandwidth
            + 3 * frequency**3
            - 3 * frequency
        ) // 6
    reflected = 2 * bandwidth - frequency
    return (reflected**3 - reflected) // 6


def jackson_coefficient_by_convolution(
    bandwidth: int,
    frequency: int,
) -> int:
    return sum(
        (bandwidth - abs(index))
        * (bandwidth - abs(frequency - index))
        for index in range(-(bandwidth - 1), bandwidth)
        if abs(frequency - index) < bandwidth
    )


def jackson_eta_cap(bandwidth: int) -> Fraction:
    c_zero = bandwidth * (2 * bandwidth * bandwidth + 1) // 3
    odd_sum = sum(
        (
            Fraction(
                jackson_coefficient(bandwidth, frequency),
                frequency**2,
            )
            for frequency in range(1, 2 * bandwidth - 2, 2)
        ),
        Fraction(0),
    )
    return Fraction(1, 2) - Fraction(
        4 * SHARP_PI_UPPER_DENOMINATOR**2,
        SHARP_PI_UPPER_NUMERATOR**2 * c_zero,
    ) * odd_sum


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    ordinary_minimum, ordinary_witnesses, ordinary_cases = (
        local_core_minimum(
            ORDINARY_DANGER_RADIUS,
            target_multiplier=1,
            maximum_h=6,
        )
    )
    guard_facet_minimum, guard_facet_witnesses, guard_facet_cases = (
        local_core_minimum(
            GUARD_DANGER_RADIUS,
            target_multiplier=1,
            maximum_h=6,
        )
    )
    guard_target_minimum, guard_target_witnesses, guard_target_cases = (
        local_core_minimum(
            ORDINARY_DANGER_RADIUS,
            target_multiplier=2,
            maximum_h=3,
        )
    )

    require(ordinary_cases == 11_664, "ordinary raw case count changed")
    require(
        guard_facet_cases == 11_664,
        "guard-facet raw case count changed",
    )
    require(
        guard_target_cases == 5_928,
        "guard-target raw case count changed",
    )
    require(
        ordinary_minimum == Fraction(9, 2_401),
        "ordinary local minimum changed",
    )
    require(
        guard_facet_minimum == Fraction(36, 2_401),
        "guard-facet local minimum changed",
    )
    require(
        guard_target_minimum == Fraction(36, 2_401),
        "guard-target local minimum changed",
    )

    ordinary_canonical = tuple(
        witness
        for witness in ordinary_witnesses
        if witness == (1, (1, 0, 0, 1))
    )
    require(
        ordinary_canonical == ((1, (1, 0, 0, 1)),),
        "canonical ordinary minimizer disappeared",
    )
    require(
        all(
            h == 1
            and coefficient[0] == 1
            and sorted(coefficient[1:]) == [0, 0, 1]
            for h, coefficient in ordinary_witnesses
        ),
        "ordinary minimizer orbit changed",
    )
    require(
        all(
            h == 1
            and (
                (
                    coefficient[0] == 1
                    and sorted(coefficient[1:]) == [0, 0, 1]
                )
                or coefficient == (6, 0, 0, 0)
            )
            for h, coefficient in guard_facet_witnesses
        ),
        "guard-facet minimizer orbits changed",
    )
    require(
        all(
            h == 1
            and coefficient[0] == 1
            and sorted(coefficient[1:]) == [0, 0, 1]
            for h, coefficient in guard_target_witnesses
        ),
        "guard-target minimizer orbit changed",
    )

    # Strong residue/slice regimes. These finite sweeps certify the elementary
    # floor inequalities at the only adjacent boundaries; the theorem proves
    # their all-integer continuations.
    for h in range(7, 100_001):
        require(
            Fraction(h // 7, h) >= Fraction(1, 13),
            f"ordinary residue floor failed at h={h}",
        )
    for h in range(4, 100_001):
        require(
            Fraction((2 * h) // 7, h) >= Fraction(1, 6),
            f"guard residue floor failed at h={h}",
        )
    for coefficient in range(7, 100_001):
        require(
            Fraction(coefficient // 7, coefficient)
            >= Fraction(1, 13),
            f"ordinary danger slice failed at B={coefficient}",
        )
        require(
            Fraction((2 * coefficient) // 7, coefficient)
            >= Fraction(2, 13),
            f"guard danger slice failed at B={coefficient}",
        )
        require(
            Fraction((6 * coefficient) // 7, coefficient)
            >= Fraction(6, 13),
            f"ordinary safe slice failed at B={coefficient}",
        )

    ordinary_uniform_deficit = min(
        ordinary_minimum,
        ORDINARY_DANGER * ORDINARY_SAFE**3 / 13,
    )
    guard_facet_uniform_deficit = min(
        guard_facet_minimum,
        GUARD_DANGER * ORDINARY_SAFE**3 / 13,
    )
    guard_target_uniform_deficit = min(
        guard_target_minimum,
        GUARD_DANGER * ORDINARY_SAFE**3 / 13,
    )
    require(
        ordinary_uniform_deficit == Fraction(117, 31_213),
        "ordinary uniform deficit changed",
    )
    require(
        guard_facet_uniform_deficit == Fraction(432, 31_213),
        "guard-facet uniform deficit changed",
    )
    require(
        guard_target_uniform_deficit == Fraction(432, 31_213),
        "guard-target uniform deficit changed",
    )

    guard_basis_floor = 5 * guard_facet_uniform_deficit
    guard_extra_floor = (
        4 * ordinary_uniform_deficit
        + guard_target_uniform_deficit
    )
    uniform_rank_five_floor = min(
        guard_basis_floor,
        guard_extra_floor,
    )
    require(
        guard_basis_floor == Fraction(2_160, 31_213),
        "guard-basis floor changed",
    )
    require(
        guard_extra_floor == Fraction(900, 31_213),
        "guard-extra floor changed",
    )
    require(
        uniform_rank_five_floor == Fraction(900, 31_213),
        "uniform rank-five floor changed",
    )
    require(
        Fraction(1_296, 16_807) > uniform_rank_five_floor,
        "THM-2295 lower-rank floor no longer dominates",
    )

    tangent_one_fifth = Fraction(1, 5)
    tangent_double = (
        2 * tangent_one_fifth / (1 - tangent_one_fifth**2)
    )
    tangent_quadruple = (
        2 * tangent_double / (1 - tangent_double**2)
    )
    tangent_one_239 = Fraction(1, 239)
    tangent_difference = (
        (tangent_quadruple - tangent_one_239)
        / (1 + tangent_quadruple * tangent_one_239)
    )
    require(tangent_double == Fraction(5, 12), "Machin double changed")
    require(
        tangent_quadruple == Fraction(120, 119),
        "Machin quadruple changed",
    )
    require(tangent_difference == 1, "Machin tangent identity failed")
    machin_pi_upper = 4 * (
        4 * alternating_arctangent_partial(5, 6)
        - alternating_arctangent_partial(239, 1)
    )
    require(
        Fraction(0)
        < machin_pi_upper
        < Fraction(
            SHARP_PI_UPPER_NUMERATOR,
            SHARP_PI_UPPER_DENOMINATOR,
        ),
        "Machin rational upper bound for pi failed",
    )

    for bandwidth in (263, 264):
        for frequency in range(2 * bandwidth - 1):
            require(
                jackson_coefficient(bandwidth, frequency)
                == jackson_coefficient_by_convolution(
                    bandwidth,
                    frequency,
                ),
                "Jackson closed/convolution coefficient mismatch "
                f"at N={bandwidth}, k={frequency}",
            )

    accepted: list[int] = []
    selected: dict[int, Fraction] = {}
    for bandwidth in range(2, JACKSON_N + 1):
        eta_cap = jackson_eta_cap(bandwidth)
        require(eta_cap > 0, f"eta cap nonpositive at N={bandwidth}")
        gap = uniform_rank_five_floor - TELESCOPE_FACTOR * eta_cap
        if gap > 0:
            accepted.append(bandwidth)
        if bandwidth in (263, 264):
            selected[bandwidth] = gap
    require(
        accepted == [JACKSON_N],
        "N=264 was not the first passing exact certificate",
    )
    require(
        Fraction(-7, 1_000_000)
        < selected[263]
        < Fraction(-6, 1_000_000)
        < 0,
        "N=263 failure wrapper changed",
    )
    require(
        0
        < Fraction(102, 1_000_000)
        < selected[264]
        < Fraction(103, 1_000_000),
        "N=264 success wrapper changed",
    )
    require(
        selected[264] > Fraction(1, 10_000),
        "N=264 exact margin no longer exceeds 1/10000",
    )
    require(SCALAR_HEIGHT == 526, "scalar height changed")
    require(FIXED_SECTION_HEIGHT == 1_052, "fixed height changed")

    print("THM-2298 weighted rank-five facet and rank-six referee")
    print(
        "ordinary_local_raw_cases=11664,min=9/2401,"
        "orbit=h1:(1;0,0,1)"
    )
    print(
        "guardfacet_ordinary_local_raw_cases=11664,min=36/2401,"
        "orbits=h1:(1;0,0,1)|(6;0,0,0)"
    )
    print(
        "ordinaryfacet_guard_local_raw_cases=5928,min=36/2401,"
        "orbit=h1:(1;0,0,1)"
    )
    print(
        "uniform_deficits=ordinary:117/31213,"
        "guardfacet:432/31213,guardtarget:432/31213"
    )
    print(
        f"guard_basis_floor={guard_basis_floor},"
        f"guard_extra_floor={guard_extra_floor}"
    )
    print(f"uniform_rank5_safe_floor={uniform_rank_five_floor}")
    print("Machin_pi_upper=104348/33215")
    print("exact_Jackson_scan_N=2..264:first_accepted=264")
    print("N263_FAIL:-7/1000000<gap<-6/1000000<0")
    print("N264_PASS:1/10000<102/1000000<gap<103/1000000")
    print(
        f"scalar_rank>=6_height={SCALAR_HEIGHT},"
        f"fixed_section_height={FIXED_SECTION_HEIGHT}"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
