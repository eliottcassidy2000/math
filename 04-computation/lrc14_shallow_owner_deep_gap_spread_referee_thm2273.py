#!/usr/bin/env python3
"""Independent interval/root referee for THM-2273.

This path never calls the folded formulas in the primary companion.  It
constructs centered combs as exact rational interval unions, checks the two
finite guard/deep endpoint banks directly, verifies the relevant danger-pair
equality, audits the three root-occupancy caps on a hostile torsion grid, and
checks the geometry of deepest-successor safe gaps.
"""

from fractions import Fraction as Q
from math import gcd


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def merge(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    result: list[tuple[Q, Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not result or left > result[-1][1]:
            result.append((left, right))
        elif right > result[-1][1]:
            result[-1] = (result[-1][0], right)
    return result


def centered_comb(speed: int, denominator: int) -> list[tuple[Q, Q]]:
    """Intervals of ||speed*x||<1/denominator in [0,1], endpoints ignored."""
    radius = Q(1, denominator * speed)
    pieces: list[tuple[Q, Q]] = []
    for index in range(speed):
        center = Q(index, speed)
        left = center - radius
        right = center + radius
        if left < 0:
            pieces.append((Q(0), right))
            pieces.append((Q(1) + left, Q(1)))
        elif right > 1:
            pieces.append((left, Q(1)))
            pieces.append((Q(0), right - 1))
        else:
            pieces.append((left, right))
    return merge(pieces)


def complement(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    output: list[tuple[Q, Q]] = []
    cursor = Q(0)
    for left, right in intervals:
        if cursor < left:
            output.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        output.append((cursor, Q(1)))
    return output


def intersection_mass(
    left_intervals: list[tuple[Q, Q]],
    right_intervals: list[tuple[Q, Q]],
) -> Q:
    i = 0
    j = 0
    total = Q(0)
    while i < len(left_intervals) and j < len(right_intervals):
        left = max(left_intervals[i][0], right_intervals[j][0])
        right = min(left_intervals[i][1], right_intervals[j][1])
        if left < right:
            total += right - left
        if left_intervals[i][1] < right_intervals[j][1]:
            i += 1
        elif right_intervals[j][1] < left_intervals[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


def direct_guard_deep(guard: int, deep: int) -> Q:
    guard_safe = complement(centered_comb(guard, 7))
    danger = centered_comb(deep, 14)
    return intersection_mass(guard_safe, danger)


def direct_danger_pair(left: int, right: int) -> Q:
    return intersection_mass(centered_comb(left, 14), centered_comb(right, 14))


def norm_numerator(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def in_danger(speed: int, numerator: int, modulus: int) -> bool:
    return 14 * norm_numerator(speed * numerator, modulus) < modulus


def in_guard(speed: int, numerator: int, modulus: int) -> bool:
    return 7 * norm_numerator(speed * numerator, modulus) > modulus


def root_count(
    predicate,
    speed: int,
    parent_numerator: int,
    parent_modulus: int,
) -> int:
    child_modulus = P * parent_modulus
    return sum(
        predicate(
            speed,
            parent_numerator + sheet * parent_modulus,
            child_modulus,
        )
        for sheet in range(P)
    )


def main() -> None:
    # Direct finite endpoint banks at the first odd and even depths.
    odd_rows = []
    for a in range(1, 3, 2):
        for k in range(1, 3):
            if a * k > 2 or gcd(a, k) != 1 or a % P == 0 or k % P == 0:
                continue
            odd_rows.append((direct_guard_deep(a, P * k), a, k))
    require(len(odd_rows) == 2, "direct odd bank size drift")
    odd_max = max(value for value, _, _ in odd_rows)
    odd_equal = {(a, k) for value, a, k in odd_rows if value == odd_max}
    require(
        odd_max == Q(5, 49) + Q(5, 49 * P),
        "direct odd guard/deep cap drift",
    )
    require(odd_equal == {(1, 1)}, "direct odd equality drift")

    even_rows = []
    for a in range(1, 15, 2):
        for k in range(1, 15 // a + 1):
            if a * k > 14 or gcd(a, k) != 1 or a % P == 0 or k % P == 0:
                continue
            even_rows.append((direct_guard_deep(a, P**2 * k), a, k))
    require(len(even_rows) == 22, "direct even bank size drift")
    even_max = max(value for value, _, _ in even_rows)
    even_equal = {(a, k) for value, a, k in even_rows if value == even_max}
    require(
        even_max == Q(5, 49) + Q(5, 294 * P**2),
        "direct even guard/deep cap drift",
    )
    require(even_equal == {(1, 6)}, "direct even equality drift")

    # An additional odd-depth equality and the shallow gap-two pair equality.
    require(
        direct_guard_deep(1, P**3)
        == Q(5, 49) + Q(5, 49 * P**3),
        "direct depth-three guard/deep equality drift",
    )
    require(
        direct_danger_pair(1, P**2) == Q(25, 1183),
        "direct gap-two danger equality drift",
    )

    # Independent exact root counts. Powers of thirteen avoid every 7/14
    # endpoint, so all strict predicates are unambiguous.
    modulus = P**2
    for parent in range(modulus):
        guard_parent = int(in_guard(1, parent, modulus))
        danger_parent = int(in_danger(1, parent, modulus))
        guard_roots = root_count(in_guard, 1, parent, modulus)
        danger_roots = root_count(in_danger, 1, parent, modulus)
        require(
            guard_roots == 10 - guard_parent,
            f"guard root law drift at parent {parent}",
        )
        require(
            danger_roots == 2 - danger_parent,
            f"danger root law drift at parent {parent}",
        )
        require(
            P - danger_roots == 11 + danger_parent,
            f"danger-complement root law drift at parent {parent}",
        )

    # Direct safe-gap geometry of D_w^c.
    for speed in (1, P, P**2):
        safe = complement(centered_comb(speed, 14))
        require(len(safe) == speed, f"safe-gap component count drift at {speed}")
        require(
            all(right - left == Q(6, 7 * speed) for left, right in safe),
            f"safe-gap length drift at {speed}",
        )
        require(
            sum((right - left for left, right in safe), Q(0)) == Q(6, 7),
            f"safe-gap total mass drift at {speed}",
        )

    # The common-time optimization is checked independently as a minimax
    # identity on the exact rational transport factors.
    a_factor = Q(169, 20)
    c_factor = Q(169, 120)
    b_factor = Q(2197, 240)
    optimum = Q(2197, 460)
    require(
        a_factor * Q(11, 23) + c_factor * Q(12, 23) == optimum,
        "direct minimax first branch drift",
    )
    require(
        b_factor * Q(12, 23) == optimum,
        "direct minimax second branch drift",
    )

    print("THM-2273 independent interval/root referee")
    print(f"direct_odd_bank={len(odd_rows)} max={odd_max} equality={sorted(odd_equal)}")
    print(
        f"direct_even_bank={len(even_rows)} "
        f"max={even_max} equality={sorted(even_equal)}"
    )
    print(f"direct_gap_two_danger={direct_danger_pair(1, P**2)}")
    print(f"root_laws_checked_parents={modulus}")
    print("safe_gap_speeds_checked=(1,13,169)")
    print(f"common_time_minimax_factor={optimum}")
    print("status=INDEPENDENT_DIRECT_INTERVAL_ROOT_PASS")


if __name__ == "__main__":
    main()
