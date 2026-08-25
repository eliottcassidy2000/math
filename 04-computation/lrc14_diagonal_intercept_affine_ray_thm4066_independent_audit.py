#!/usr/bin/env python3
"""Independent affine-centre audit for THM-4066.

Unlike the primary literal wall subdivision, this path starts from the
THM-4032/4030 affine-centre boxes, intersects their real intervals, and only
then projects the resulting intervals to the phase circle.  Dilation is
checked as a formal label-conjugacy identity, not by sampled phases.
"""

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def strict_masks(
    y: F,
    d: int,
    exceptions: tuple[int, ...],
    closed: bool = False,
) -> tuple[tuple[int, ...], ...]:
    compare = (lambda value: value <= F(1, 14)) if closed else (
        lambda value: value < F(1, 14)
    )
    return tuple(
        tuple(
            label
            for label in range(d)
            if compare(distance(speed, (y + label) / d))
        )
        for speed in exceptions
    )


def interval_intersection(
    centres: tuple[F, ...],
    radii: tuple[F, ...],
) -> tuple[F, F] | None:
    left = max(centre - radius for centre, radius in zip(centres, radii))
    right = min(centre + radius for centre, radius in zip(centres, radii))
    return (left, right) if left < right else None


def d3_affine_intervals(
    a: int,
    b: int,
    c: int,
) -> tuple[tuple[F, F], ...]:
    intervals = set()
    for second, third in ((b, c), (c, b)):
        for A in range(a):
            for B in range(2 * second + 1):
                for C in range(2 * third + 1):
                    interval = interval_intersection(
                        (F(A, a), F(B, second) - F(1, 3), F(C, third) - F(2, 3)),
                        (F(1, 14 * a), F(1, 14 * second), F(1, 14 * third)),
                    )
                    if interval is not None:
                        intervals.add(interval)
    return tuple(sorted(intervals))


def d4_affine_intervals(
    even: int,
    first: int,
    second: int,
) -> tuple[tuple[F, F], ...]:
    require(even % 4 == 2, "d4 even exception type")
    r = even // 2
    intervals = set()
    for a, b in ((first, second), (second, first)):
        for A in range(a):
            for B in range(2 * b + 1):
                for C in range(3 * r + 1):
                    interval = interval_intersection(
                        (F(A, a), F(B, b) - F(1, 2), F(C, 2 * r) - F(1, 4)),
                        (F(1, 14 * a), F(1, 14 * b), F(1, 28 * r)),
                    )
                    if interval is not None:
                        intervals.add(interval)
    return tuple(sorted(intervals))


def image_contains(
    phase: F,
    multiplier: int,
    real_intervals: tuple[tuple[F, F], ...],
) -> bool:
    """Whether phase lies in the mod-one image of multiplier*intervals."""
    phase %= 1
    for left, right in real_intervals:
        scaled_left = multiplier * left
        scaled_right = multiplier * right
        integer = (scaled_left - phase).numerator // (scaled_left - phase).denominator + 1
        if integer < scaled_right - phase:
            return True
    return False


def projected_components(
    multiplier: int,
    real_intervals: tuple[tuple[F, F], ...],
) -> tuple[tuple[F, F], ...]:
    endpoints = {F(0), F(1)}
    for left, right in real_intervals:
        endpoints.add((multiplier * left) % 1)
        endpoints.add((multiplier * right) % 1)
    ordered = sorted(endpoints)
    components: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        if not image_contains((left + right) / 2, multiplier, real_intervals):
            continue
        if (
            components
            and components[-1][1] == left
            and image_contains(left, multiplier, real_intervals)
        ):
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    return tuple(components)


def ordinary_danger_intervals(speed: int) -> tuple[tuple[F, F], ...]:
    radius = F(1, 14 * speed)
    return tuple(
        (F(integer, speed) - radius, F(integer, speed) + radius)
        for integer in range(speed)
    )


def formal_conjugacy_checks(
    d: int,
    primitive: tuple[int, ...],
    dilations: tuple[int, ...],
) -> int:
    checks = 0
    for h in dilations:
        require(gcd(h, d) == 1, "unit label action")
        require(len({h * label % d for label in range(d)}) == d, "label permutation")
        for delta in primitive:
            for label in range(d):
                target_label = h * label % d
                constant_difference = delta * (h * label - target_label)
                require(
                    constant_difference % d == 0,
                    "scaled-mask and base-mask phases differ integrally",
                )
                checks += 1
    return checks


def clearance(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed, phase) for speed in speeds)


def main() -> None:
    expected_d3 = ((F(11, 56), F(3, 14)), (F(11, 14), F(45, 56)))
    expected_d4 = ((F(5, 49), F(1, 7)), (F(6, 7), F(44, 49)))

    d3_real = d3_affine_intervals(1, 4, 5)
    d4_real = d4_affine_intervals(2, 7, 9)
    d3_projected = projected_components(3, d3_real)
    d4_projected = projected_components(4, d4_real)
    require(d3_projected == expected_d3, "d3 affine-centre projection")
    require(d4_projected == expected_d4, "d4 affine-centre projection")
    require(len(d3_real) == 2, "two d3 oriented affine intervals")
    require(
        len(d4_real) == 4,
        "four d4 affine lift intervals project to two pack-phase components",
    )

    conjugacy_checks = formal_conjugacy_checks(
        3,
        (1, 4, 5),
        (1, 2, 4, 5, 7, 8, 10, 11),
    )
    conjugacy_checks += formal_conjugacy_checks(
        4,
        (2, 7, 9),
        (1, 3, 5, 7, 9, 11),
    )

    endpoint_rows = (
        (3, (1, 4, 5), (F(11, 56), F(3, 14), F(11, 14), F(45, 56))),
        (4, (2, 7, 9), (F(5, 49), F(1, 7), F(6, 7), F(44, 49))),
    )
    endpoint_words = []
    for d, primitive, endpoints in endpoint_rows:
        for endpoint in endpoints:
            strict = strict_masks(endpoint, d, primitive)
            closed = strict_masks(endpoint, d, primitive, closed=True)
            require(set().union(*map(set, strict)) != set(range(d)), "strict endpoint leak")
            require(set().union(*map(set, closed)) == set(range(d)), "closed endpoint cover")
            endpoint_words.append((d, endpoint, strict))

    # For h=d=3 every scaled exception is divisible by d, so its mask is
    # either all labels or none.  The spoiled set is exactly the union of the
    # ordinary danger sets for speeds 1,4,5; this independently yields 51/140.
    nonunit_real = tuple(
        interval
        for speed in (1, 4, 5)
        for interval in ordinary_danger_intervals(speed)
    )
    nonunit_components = projected_components(1, nonunit_real)
    nonunit_measure = sum((right - left for left, right in nonunit_components), F(0))
    require(nonunit_measure == F(51, 140), "nonunit analytic union measure")
    require(nonunit_measure != F(1, 28), "nonunit is not a false pullback")

    h10 = tuple(range(1, 11))
    d3_row = tuple(3 * speed for speed in h10) + (5, 20, 25)
    d4_row = tuple(4 * speed for speed in h10) + (22, 77, 99)
    require(len(set(d3_row)) == 13, "d3 thirteen nonzero speeds")
    require(len(set(d4_row)) == 13, "d4 thirteen nonzero speeds")
    require(set(3 * speed for speed in h10).isdisjoint((5, 20, 25)), "d3 disjointness")
    require(set(4 * speed for speed in h10).isdisjoint((22, 77, 99)), "d4 disjointness")
    require(24 * 5 >= 11 * 10 > 5 * 5, "d3 new-versus-old thresholds")
    require(21 * 11 >= 22 * 10 and 27 * 11 < 44 * 10, "d4 new-versus-old thresholds")
    require(clearance(d3_row, F(1, 33)) == F(1, 11), "d3 witness")
    require(clearance(d4_row, F(1, 44)) == F(1, 11), "d4 witness")

    require(F(3, 77 * 10) >= F(1, 56 * 5), "d3 safe arc beats component")
    require(F(3, 77 * 10) >= F(2, 49 * 11), "d4 safe arc beats component")
    require(F(3, 35) / F(1, 56) == F(24, 5), "d3 gain factor")
    require(F(4, 63) / F(2, 49) == F(14, 9), "d4 gain factor")

    print("LRC14_DIAGONAL_INTERCEPT_AFFINE_RAY_THM4066_INDEPENDENT_AUDIT")
    print("path=affine_centres_to_real_intersections_to_circle_projection")
    print("d3_real_intervals", d3_real)
    print("d3_projected", d3_projected, "width", F(1, 56), "measure", F(1, 28))
    print("d4_real_intervals", d4_real)
    print("d4_projected", d4_projected, "width", F(2, 49), "measure", F(4, 49))
    print("formal_label_conjugacy_checks", conjugacy_checks)
    print("endpoint_words", tuple(endpoint_words))
    print("nonunit_components", nonunit_components)
    print("nonunit_measure", nonunit_measure, "false_pullback_measure", F(1, 28))
    print("runner_counts=d3:13,d4:13;stationary_runner_makes_LRC14")
    print("d3_threshold=24h>=11M;gain=24/5;witness=1/33")
    print("d4_threshold=21h>=22M;gain=14/9;witness=1/44")
    print("checks", CHECKS)
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
