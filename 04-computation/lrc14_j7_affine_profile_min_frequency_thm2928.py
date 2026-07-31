#!/usr/bin/env python3
"""Dependency-free exact referee for critical seven-tail affine-profile rays."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import lcm


RADIUS = F(1, 14)
SEVENTH = F(1, 7)
RULER = 5_045_040


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(value: F) -> F:
    return value - value.numerator // value.denominator


def circle_norm(value: F) -> F:
    remainder = fractional_part(value)
    return min(remainder, 1 - remainder)


def merge_intervals(
    intervals: list[tuple[F, F]],
) -> tuple[tuple[F, F], ...]:
    rows = sorted((left, right) for left, right in intervals if left < right)
    if not rows:
        return ()
    merged: list[list[F]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def profile_by_union(residues: tuple[int, ...], time: F) -> F:
    intervals: list[tuple[F, F]] = []
    for residue in residues:
        center = fractional_part(-residue * time)
        left = center - RADIUS
        right = center + RADIUS
        if left < 0:
            intervals.append((left + 1, F(1)))
            intervals.append((F(0), right))
        elif right > 1:
            intervals.append((left, F(1)))
            intervals.append((F(0), right - 1))
        else:
            intervals.append((left, right))
    danger = merge_intervals(intervals)
    return 1 - sum((right - left for left, right in danger), F(0))


def shifted_comb_intervals(
    speed: int,
    phase: F,
) -> tuple[tuple[F, F], ...]:
    phase = fractional_part(phase)
    radius = RADIUS / speed
    rows: list[tuple[F, F]] = []
    for index in range(speed):
        center = fractional_part(F(index, speed) - phase / speed)
        left = center - radius
        right = center + radius
        if left < 0:
            rows.append((left + 1, F(1)))
            rows.append((F(0), right))
        elif right > 1:
            rows.append((left, F(1)))
            rows.append((F(0), right - 1))
        else:
            rows.append((left, right))
    return tuple(sorted(rows))


def general_profile_and_l1(
    slopes: tuple[int, ...],
    residues: tuple[int, ...],
    time: F,
) -> tuple[F, F]:
    event_intervals = tuple(
        shifted_comb_intervals(slope, residue * time)
        for slope, residue in zip(slopes, residues)
    )
    all_intervals = [
        interval for event in event_intervals for interval in event
    ]
    union = merge_intervals(all_intervals)
    profile = 1 - sum((right - left for left, right in union), F(0))
    points = sorted(
        {
            F(0),
            F(1),
            *(
                endpoint
                for left, right in all_intervals
                for endpoint in (left, right)
            ),
        }
    )
    l1 = F(0)
    total_multiplicity = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        multiplicity = sum(
            any(a < midpoint < b for a, b in event)
            for event in event_intervals
        )
        width = right - left
        l1 += width * abs(multiplicity - 1)
        total_multiplicity += width * multiplicity
    require(total_multiplicity == 1, "seven danger masses did not sum to one")
    require(l1 == 2 * profile, "multiplicity L1/profile identity failed")
    return profile, l1


def polynomial_div_exact(
    numerator: tuple[int, ...],
    denominator: tuple[int, ...],
) -> tuple[int, ...]:
    rows = list(numerator)
    while len(rows) >= len(denominator):
        lead = rows[-1]
        require(
            lead % denominator[-1] == 0,
            "cyclotomic polynomial division was not exact",
        )
        quotient = lead // denominator[-1]
        offset = len(rows) - len(denominator)
        for index, value in enumerate(denominator):
            rows[offset + index] -= quotient * value
        while rows and rows[-1] == 0:
            rows.pop()
    require(not rows, "cyclotomic polynomial left a remainder")
    quotient_degree = len(numerator) - len(denominator)
    # Re-run long division while retaining the quotient coefficients.
    rows = list(numerator)
    result = [0] * (quotient_degree + 1)
    while len(rows) >= len(denominator):
        quotient = rows[-1] // denominator[-1]
        offset = len(rows) - len(denominator)
        result[offset] = quotient
        for index, value in enumerate(denominator):
            rows[offset + index] -= quotient * value
        while rows and rows[-1] == 0:
            rows.pop()
    return tuple(result)


@lru_cache(maxsize=None)
def cyclotomic(order: int) -> tuple[int, ...]:
    polynomial = [-1] + [0] * (order - 1) + [1]
    for divisor in range(1, order):
        if order % divisor == 0:
            polynomial = list(
                polynomial_div_exact(tuple(polynomial), cyclotomic(divisor))
            )
    return tuple(polynomial)


def polynomial_remainder(
    numerator: tuple[int, ...],
    denominator: tuple[int, ...],
) -> tuple[int, ...]:
    rows = list(numerator)
    while len(rows) >= len(denominator):
        quotient = rows[-1] // denominator[-1]
        offset = len(rows) - len(denominator)
        for index, value in enumerate(denominator):
            rows[offset + index] -= quotient * value
        while rows and rows[-1] == 0:
            rows.pop()
    return tuple(rows)


def minimum_frequency_sum_is_zero(
    minimum_residues: tuple[int, ...],
    time: F,
) -> bool:
    order = time.denominator
    if order == 1:
        return False
    coefficients = [0] * order
    for residue in minimum_residues:
        coefficients[(residue * time.numerator) % order] += 1
    return not polynomial_remainder(tuple(coefficients), cyclotomic(order))


def center_gaps(
    residues: tuple[int, ...],
    time: F,
) -> tuple[F, ...]:
    centers = tuple(sorted(fractional_part(residue * time) for residue in residues))
    return tuple(
        (
            centers[(index + 1) % len(centers)] - centers[index]
            if index + 1 < len(centers)
            else centers[0] + 1 - centers[index]
        )
        for index in range(len(centers))
    )


def profile_by_gaps(residues: tuple[int, ...], time: F) -> F:
    return sum(
        (max(F(0), gap - SEVENTH) for gap in center_gaps(residues, time)),
        F(0),
    )


def is_seventh_coset(residues: tuple[int, ...], time: F) -> bool:
    return all(gap == SEVENTH for gap in center_gaps(residues, time))


def seventh_grid_distance(value: F) -> F:
    return circle_norm(7 * value) / 7


def profile_breakpoints(residues: tuple[int, ...]) -> tuple[F, ...]:
    points = {F(0), F(1)}
    for first, second in combinations(residues, 2):
        difference = abs(first - second)
        for numerator in range(7 * difference + 1):
            if numerator % 7 in (0, 1, 6):
                points.add(F(numerator, 7 * difference))
    return tuple(sorted(points))


def integrate_profile(
    residues: tuple[int, ...],
    left: F,
    right: F,
) -> F:
    points = tuple(
        sorted(
            {
                left,
                right,
                *(
                    point
                    for point in profile_breakpoints(residues)
                    if left < point < right
                ),
            }
        )
    )
    total = F(0)
    for first, second in zip(points, points[1:]):
        first_value = profile_by_gaps(residues, first)
        second_value = profile_by_gaps(residues, second)
        midpoint = (first + second) / 2
        require(
            2 * profile_by_gaps(residues, midpoint)
            == first_value + second_value,
            "profile breakpoint list missed an affine change",
        )
        total += (second - first) * (first_value + second_value) / 2
    return total


def norm_primitive(value: F) -> F:
    cycles = value.numerator // value.denominator
    remainder = value - cycles
    if remainder <= F(1, 2):
        partial = remainder * remainder / 2
    else:
        partial = F(1, 8) + (remainder - F(1, 2)) / 2 - (
            remainder - F(1, 2)
        ) ** 2 / 2
    return F(cycles, 4) + partial


def integrate_norm(
    frequency: int,
    left: F,
    length: F,
) -> F:
    return (
        norm_primitive(frequency * (left + length))
        - norm_primitive(frequency * left)
    ) / frequency


def minimum_norm_integral(frequency: int, length: F) -> F:
    scaled = frequency * length
    whole = scaled.numerator // scaled.denominator
    remainder = scaled - whole
    return (whole + remainder * remainder) / (4 * frequency)


def body_safe_ranges(body: tuple[int, ...]) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Exact open-cell ranges of the literal six-body carrier."""
    ruler = lcm(*(14 * speed for speed in body))
    danger: list[tuple[int, int]] = []
    for speed in body:
        half_tooth = ruler // (14 * speed)
        period = ruler // speed
        for tooth in range(speed + 1):
            center = tooth * period
            danger.append(
                (
                    max(0, center - half_tooth),
                    min(ruler, center + half_tooth),
                )
            )
    merged: list[list[int]] = []
    for left, right in sorted(danger):
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    safe: list[tuple[int, int]] = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < ruler:
        safe.append((cursor, ruler))
    require(safe, f"empty body carrier: {body}")
    return ruler, tuple(safe)


def minimum_longest_body_component() -> tuple[F, tuple[int, ...]]:
    """Independent all-3,003-body longest-component census."""
    minimum: tuple[F, tuple[int, ...]] | None = None
    for body in combinations(range(1, 15), 6):
        ruler, safe = body_safe_ranges(body)
        longest = F(max(right - left for left, right in safe), ruler)
        candidate = (longest, body)
        if minimum is None or candidate < minimum:
            minimum = candidate
    require(minimum is not None, "body census was empty")
    return minimum


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    packets = (
        (0, 1, 2, 3, 4, 5, 6),
        (0, 2, 4, 6, 8, 10, 12),
        (1, 4, 7, 11, 16, 22, 29),
        (17, 22, 23, 26, 50, 77, 108),
        (16, 18, 22, 26, 44, 52, 182),
    )
    times = {
        F(0),
        F(1, 7),
        F(1, 14),
        F(3, 28),
        F(13, 91),
        F(22, 147),
        F(1, 2),
    }
    for denominator in range(2, 43):
        for numerator in range(denominator):
            times.add(F(numerator, denominator))

    point_checks = 0
    zero_checks = 0
    pair_checks = 0
    for residues in packets:
        require(len(residues) == len(set(residues)) == 7, "bad residue packet")
        for time in sorted(times):
            profile = profile_by_union(residues, time)
            gaps = center_gaps(residues, time)
            require(profile == profile_by_gaps(residues, time), "gap law failed")
            require(
                sum(abs(gap - SEVENTH) for gap in gaps) == 2 * profile,
                "L1 gap identity failed",
            )
            coset = is_seventh_coset(residues, time)
            require((profile == 0) == coset, "zero/coset equivalence failed")
            point_checks += 1
            zero_checks += profile == 0
            for first, second in combinations(residues, 2):
                require(
                    profile
                    >= seventh_grid_distance((first - second) * time) / 2,
                    "pair-distance profile floor failed",
                )
                pair_checks += 1

    affine_packets = (
        (
            (1, 1, 1, 1, 1, 1, 1),
            (0, 1, 2, 3, 4, 5, 6),
        ),
        (
            (1, 2, 3, 4, 5, 6, 7),
            (0, 1, 3, 6, 10, 15, 21),
        ),
        (
            (1, 2, 2, 3, 3, 4, 5),
            (4, 0, 7, 2, 11, 5, 13),
        ),
        (
            (1, 1, 2, 2, 3, 3, 4),
            (0, 7, 1, 8, 2, 9, 3),
        ),
        (
            (2, 2, 3, 3, 5, 5, 8),
            (1, 8, 2, 9, 3, 10, 4),
        ),
    )
    affine_checks = 0
    affine_zero_checks = 0
    unique_minimum_positive_checks = 0
    for slopes, residues in affine_packets:
        minimum_slope = min(slopes)
        minimum_residues = tuple(
            residue
            for slope, residue in zip(slopes, residues)
            if slope == minimum_slope
        )
        formal_coefficients = Counter(minimum_residues)
        require(
            formal_coefficients
            and all(value > 0 for value in formal_coefficients.values()),
            "minimum-frequency exponential polynomial vanished formally",
        )
        for time in sorted(times):
            profile, _ = general_profile_and_l1(slopes, residues, time)
            if profile == 0:
                require(
                    minimum_frequency_sum_is_zero(minimum_residues, time),
                    "zero profile violated the minimum-frequency obstruction",
                )
                affine_zero_checks += 1
            if len(minimum_residues) == 1:
                require(
                    profile > 0,
                    "unique minimum slope lost its pointwise profile floor",
                )
                unique_minimum_positive_checks += 1
            affine_checks += 1

    interval_checks = 0
    equality_checks = 0
    for multiplier in range(1, 41):
        frequency = 7 * multiplier
        for denominator in range(7, 85, 7):
            for numerator in range(1, denominator // 7 + 1):
                length = F(numerator, denominator)
                if length > SEVENTH:
                    continue
                exact_minimum = minimum_norm_integral(frequency, length)
                require(
                    exact_minimum >= 7 * length * length / 4,
                    "multiple-of-seven interval floor failed",
                )
                scaled = frequency * length
                remainder = fractional_part(scaled)
                minimizing_left = fractional_part(-remainder / (2 * frequency))
                require(
                    integrate_norm(frequency, minimizing_left, length)
                    == exact_minimum,
                    "claimed norm-integral minimum was not attained",
                )
                for left in (
                    F(0),
                    F(1, 17),
                    F(5, 43),
                    minimizing_left,
                ):
                    require(
                        integrate_norm(frequency, left, length) >= exact_minimum,
                        "hostile interval beat the exact norm minimum",
                    )
                    interval_checks += 1
                equality_checks += exact_minimum == 7 * length * length / 4

    profile_integral_checks = 0
    interval_bank = (
        (F(1, 100), F(1, 14)),
        (F(23, 1092), F(23, 1092)),
        (F(71, 980), F(3, 98)),
        (F(117, 1001), F(1, 13)),
    )
    for residues in packets:
        difference = abs(residues[1] - residues[0])
        for left, length in interval_bank:
            right = left + length
            exact_profile = integrate_profile(residues, left, right)
            exact_pair_floor = (
                integrate_norm(7 * difference, left, length) / 14
            )
            require(
                exact_profile >= exact_pair_floor >= length * length / 8,
                "integrated profile floor failed",
            )
            profile_integral_checks += 1

    longest, longest_body = minimum_longest_body_component()
    require(
        (longest, longest_body) == (F(23, 1092), (1, 6, 7, 8, 10, 13)),
        "all-body longest-component minimum changed",
    )
    profile_floor = longest * longest / 8
    require(profile_floor == F(529, 9_539_712), "uniform profile floor changed")
    threshold_ratio = F(35, 2) / profile_floor
    threshold = threshold_ratio.numerator // threshold_ratio.denominator + 1
    require(
        threshold == 315_586
        and F(35, 2 * threshold) < profile_floor
        and F(35, 2 * (threshold - 1)) >= profile_floor,
        "uniform common-block threshold changed",
    )

    print("LRC14 j7 affine-profile minimum-frequency exact referee")
    print(f"packets={len(packets)}")
    print(f"point_checks={point_checks}")
    print(f"zero_profile_checks={zero_checks}")
    print(f"pair_distance_checks={pair_checks}")
    print(f"general_affine_checks={affine_checks}")
    print(f"general_affine_zero_checks={affine_zero_checks}")
    print(
        "unique_minimum_positive_checks="
        f"{unique_minimum_positive_checks}"
    )
    print(f"interval_formula_checks={interval_checks}")
    print(f"interval_equality_checks={equality_checks}")
    print(f"profile_integral_checks={profile_integral_checks}")
    print(
        "zero_law=Phi_b(t)=0 iff centers form a coset of "
        "(1/7)Z/Z"
    )
    print(
        "pointwise_floor=Phi_b(t)>=dist((b_i-b_j)t,"
        "(1/7)Z/Z)/2"
    )
    print(
        "general_affine_zero_obstruction=at x-frequency min(c), "
        "hat_d(1)*sum_(c_i=min(c)) exp(2*pi*i*r_i*t)=0"
    )
    print(
        "general_affine_consequence=no fixed seven-tail affine profile "
        "vanishes on a positive interval"
    )
    print(
        "interval_formula=inf_I integral_I ||n*t|| dt="
        "(floor(n*ell)+frac(n*ell)^2)/(4*n)"
    )
    print(f"longest_component={ftext(longest)}")
    print(f"longest_component_body={longest_body}")
    print(f"uniform_profile_floor={ftext(profile_floor)}")
    print(f"uniform_common_block_k_min={threshold}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
