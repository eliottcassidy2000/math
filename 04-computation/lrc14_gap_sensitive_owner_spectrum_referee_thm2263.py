#!/usr/bin/env python3
"""Independent rational interval-intersection referee for THM-2263."""

from fractions import Fraction
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def danger_intervals(speed):
    """Disjoint half-open representatives; endpoint choices do not affect mass."""
    pieces = []
    for index in range(speed):
        left = Fraction(14 * index - 1, 14 * speed)
        right = Fraction(14 * index + 1, 14 * speed)
        if left < 0:
            pieces.append((Fraction(), right))
            pieces.append((left + 1, Fraction(1)))
        elif right > 1:
            pieces.append((left, Fraction(1)))
            pieces.append((Fraction(), right - 1))
        else:
            pieces.append((left, right))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return tuple(merged)


def direct_overlap(left_speed, right_speed):
    left = danger_intervals(left_speed)
    right = danger_intervals(right_speed)
    i = 0
    j = 0
    total = Fraction()
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        stop = min(left[i][1], right[j][1])
        if start < stop:
            total += stop - start
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def gap_upper(depth):
    power = P**depth
    if depth % 2 == 0:
        return Fraction(1, 49) + Fraction(6, 49 * power)
    return Fraction(1, 49) + Fraction(5, 588 * power)


def gap_lower(depth):
    power = P**depth
    if depth % 2 == 0:
        return Fraction(1, 49) - Fraction(5, 588 * power)
    return Fraction(1, 49) - Fraction(6, 49 * power)


def main():
    endpoint_checks = 0
    for depth in range(1, 5):
        power = P**depth
        upper_pair = (1, power) if depth % 2 == 0 else (12, power)
        lower_pair = (12, power) if depth % 2 == 0 else (1, power)
        require(
            direct_overlap(*upper_pair) == gap_upper(depth),
            f"direct upper endpoint failed at gap {depth}",
        )
        require(
            direct_overlap(*lower_pair) == gap_lower(depth),
            f"direct lower endpoint failed at gap {depth}",
        )
        endpoint_checks += 2

    # Directly recover the universal maximum and second maximum among the
    # only reduced pairs that can beat the analytic ab>=10 tail.
    small_pairs = [
        (direct_overlap(a, b), a, b)
        for a in range(1, 10)
        for b in range(a + 1, 10)
        if a * b <= 9 and gcd(a, b) == 1
    ]
    require(max(small_pairs) == (Fraction(1, 14), 1, 2), "pair maximum drift")
    nonmax = [row for row in small_pairs if row[1:] != (1, 2)]
    require(
        max(row[0] for row in nonmax) == Fraction(1, 21),
        "pair second maximum drift",
    )

    # Independent direct measures for both simultaneous equality carriers.
    strict = (1, P**2, P**4)
    strict_sum = sum(
        (
            direct_overlap(strict[i], strict[j])
            for i in range(3)
            for j in range(i + 1, 3)
        ),
        Fraction(),
    )
    require(strict_sum == Fraction(12531, 199927), "strict carrier drift")

    repeated = (1, 2, P**4)
    repeated_sum = sum(
        (
            direct_overlap(repeated[i], repeated[j])
            for i in range(3)
            for j in range(i + 1, 3)
        ),
        Fraction(),
    )
    require(repeated_sum == Fraction(3206, 28561), "repeated carrier drift")

    print("THM-2263 INDEPENDENT INTERVAL REFEREE")
    print(f"direct endpoint checks: {endpoint_checks}")
    print(f"small same-valuation pairs: {len(small_pairs)}")
    print(f"strict equality carrier sum: {strict_sum}")
    print(f"repeated equality carrier sum: {repeated_sum}")
    print("all independent checks passed")


if __name__ == "__main__":
    main()
