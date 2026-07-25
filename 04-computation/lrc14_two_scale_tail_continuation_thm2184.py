#!/usr/bin/env python3
"""Exact one-dimensional controls for THM-2184."""

from fractions import Fraction
from math import lcm


RADIUS = Fraction(1, 14)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def danger_arcs(speeds: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    arcs: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        half_width = Fraction(1, 14 * speed)
        for index in range(speed):
            left = (Fraction(index, speed) - half_width) % 1
            right = left + 2 * half_width
            if right <= 1:
                arcs.append((left, right))
            else:
                arcs.append((left, Fraction(1)))
                arcs.append((Fraction(0), right - 1))
    return arcs


def safe_measure(speeds: tuple[int, ...]) -> Fraction:
    arcs = sorted(danger_arcs(speeds))
    if not arcs:
        return Fraction(1)
    danger = Fraction(0)
    current_left, current_right = arcs[0]
    for left, right in arcs[1:]:
        if left <= current_right:
            current_right = max(current_right, right)
        else:
            danger += current_right - current_left
            current_left, current_right = left, right
    danger += current_right - current_left
    return 1 - danger


def pair_by_boundary_cells(first: int, second: int) -> Fraction:
    points = {Fraction(0), Fraction(1)}
    for speed in (first, second):
        for index in range(speed):
            points.add(Fraction(14 * index - 1, 14 * speed) % 1)
            points.add(Fraction(14 * index + 1, 14 * speed) % 1)
    ordered = sorted(points)
    total = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(
            min((speed * midpoint) % 1, 1 - (speed * midpoint) % 1)
            >= RADIUS
            for speed in (first, second)
        ):
            total += right - left
    return total


def main() -> None:
    profile_pair = Fraction(36, 49)
    profile_at_zero = safe_measure((1, 1))
    require(profile_at_zero == Fraction(6, 7), "profile value at zero")
    require(profile_at_zero != profile_pair, "control profile became constant")
    expected_pair_deviations = (
        -Fraction(2, 3185),
        -Fraction(2, 12789),
        -Fraction(2, 80017),
        -Fraction(2, 320117),
    )
    for index, level in enumerate((1, 2, 5, 10)):
        row = (14 * level - 1, 14 * level + 1)
        union_value = safe_measure(row)
        cell_value = pair_by_boundary_cells(*row)
        require(union_value == cell_value, "pair evaluators disagree")
        deviation = union_value - profile_pair
        require(deviation == expected_pair_deviations[index], "pair deviation")
        bound = Fraction(5 * 2, 2 * 14 * level)
        require(abs(deviation) <= bound, "pair escaped theorem bound")

    core = (1, 2, 3, 4, 5, 6, 8)
    normalized_tail = (3, 12)
    endpoint_modulus = 1
    for speed in core:
        endpoint_modulus = lcm(endpoint_modulus, 14 * speed)
    require(endpoint_modulus == 1680, "endpoint modulus")

    core_mass = safe_measure(core)
    tail_mass = safe_measure(normalized_tail)
    product_profile = core_mass * tail_mass
    require(core_mass == Fraction(27, 70), "core mass")
    require(tail_mass == Fraction(3, 4), "tail mass")
    require(product_profile == Fraction(81, 280), "product profile")

    expected_proportional_deviations = (
        Fraction(47, 1412040),
        Fraction(47, 2823240),
        Fraction(47, 7056840),
    )
    for index, level in enumerate((1, 2, 5)):
        scale = level * endpoint_modulus + 1
        row = core + tuple(scale * speed for speed in normalized_tail)
        deviation = safe_measure(row) - product_profile
        require(
            deviation == expected_proportional_deviations[index],
            "proportional deviation",
        )
        bound = Fraction(5 * sum(normalized_tail), 2 * level * endpoint_modulus)
        require(abs(deviation) <= bound, "proportional bound")

    aligned_row = core + tuple(
        endpoint_modulus * speed for speed in normalized_tail
    )
    require(safe_measure(aligned_row) == product_profile, "exact aligned product")

    print("THM-2184 exact two-scale continuation controls")
    print(f"nonconstant_profile_at_zero={profile_at_zero}")
    print("nonconstant_profile=36/49")
    print(f"nonconstant_levels={(1, 2, 5, 10)}")
    print(f"nonconstant_deviations={expected_pair_deviations}")
    print(f"core={core}")
    print(f"endpoint_modulus={endpoint_modulus}")
    print(f"core_mass={core_mass}")
    print(f"normalized_tail={normalized_tail}")
    print(f"normalized_tail_mass={tail_mass}")
    print(f"proportional_profile={product_profile}")
    print(f"proportional_levels={(1, 2, 5)}")
    print(f"proportional_deviations={expected_proportional_deviations}")
    print("aligned_deviation=0")
    print("independent_pair_evaluators=interval_union,boundary_cells")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
