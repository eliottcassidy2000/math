#!/usr/bin/env python3
"""Exact hostile checks for the proposed THM-4088.

This is not a substitute for the elementary proof.  It checks the two sharp
boundaries (strictness and the tight initial segment) and exhausts the
Lipschitz implication on a small rational universe.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def dist_z(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def circle_dist(x: Fraction, y: Fraction) -> Fraction:
    return dist_z(x - y)


def lonely_profile(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(dist_z(v * t) for v in speeds)


def require(condition: bool, payload) -> None:
    if not condition:
        raise RuntimeError(payload)


def all_nonempty_subsets(n: int):
    base = range(1, n + 1)
    for size in range(1, n + 1):
        yield from combinations(base, size)


def tight_initial_segment_audit() -> int:
    checked = 0
    for n in range(2, 15):
        denominator = 60 * n
        equality = set()
        strict = set()
        speeds = tuple(range(1, n))
        threshold = Fraction(1, n)
        for k in range(denominator):
            value = lonely_profile(speeds, Fraction(k, denominator))
            checked += 1
            if value == threshold:
                equality.add(k)
            elif value > threshold:
                strict.add(k)
        expected = {60 * a for a in range(1, n) if gcd(a, n) == 1}
        require(equality == expected, (n, equality, expected))
        require(not strict, (n, strict))
    return checked


def lipschitz_and_margin_audit() -> tuple[int, int]:
    grid_denominator = 84
    grid = [Fraction(k, grid_denominator) for k in range(grid_denominator)]
    threshold = Fraction(1, 14)
    lipschitz_checks = 0
    margin_checks = 0
    for speeds in all_nonempty_subsets(7):
        bound = max(speeds)
        profile = [lonely_profile(speeds, t) for t in grid]
        for i, t in enumerate(grid):
            for j, u in enumerate(grid):
                require(
                    abs(profile[i] - profile[j]) <= bound * circle_dist(t, u),
                    (speeds, t, u, profile[i], profile[j]),
                )
                lipschitz_checks += 1
                margin = profile[i] - threshold
                if margin > 0 and circle_dist(t, u) < margin / bound:
                    require(profile[j] > threshold, (speeds, t, u, margin))
                    margin_checks += 1
    return lipschitz_checks, margin_checks


def open_radius_boundary_audit() -> None:
    speeds = (1,)
    threshold = Fraction(1, 4)
    center = Fraction(1, 2)
    margin = lonely_profile(speeds, center) - threshold
    radius = margin / max(speeds)
    require(radius == Fraction(1, 4), radius)
    require(lonely_profile(speeds, center - radius) == threshold, "left endpoint")
    require(lonely_profile(speeds, center + radius) == threshold, "right endpoint")


def main() -> None:
    tight_checks = tight_initial_segment_audit()
    lipschitz_checks, margin_checks = lipschitz_and_margin_audit()
    open_radius_boundary_audit()
    print("THM-4088 independent exact probe")
    print(f"tight-initial-segment grid values: {tight_checks}")
    print(f"circle-Lipschitz pairs: {lipschitz_checks}")
    print(f"strict-margin implications: {margin_checks}")
    print("open-radius endpoint hostile: PASS")
    print("PASS")


if __name__ == "__main__":
    main()
