#!/usr/bin/env python3
"""Exact companion for THM-2440.

All coverage radii and endpoint decisions use fractions.  Touching open
teeth are merged because the theorem concerns almost-everywhere coverage.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


SHARP = Fraction(15, 182)
TARGET = Fraction(8, 91)
HALF = Fraction(1, 2)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def distance_to_integer(x: Fraction) -> Fraction:
    residue = x % 1
    return min(residue, 1 - residue)


def in_danger(speed: int, y: Fraction) -> bool:
    return distance_to_integer(speed * y) < Fraction(1, 14)


def closed_center_count(u: int, n: int) -> int:
    d = gcd(u, n)
    order = n // d
    return d * (2 * (order // 14) + 1)


def half_circle_teeth(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    """Closed supports on [0,1/2], used only for a.e. seam merging."""

    intervals: list[tuple[Fraction, Fraction]] = []
    for tooth in range(speed // 2 + 2):
        left = Fraction(14 * tooth - 1, 14 * speed)
        right = Fraction(14 * tooth + 1, 14 * speed)
        left = max(Fraction(0), left)
        right = min(HALF, right)
        if left <= right and right >= 0 and left <= HALF:
            intervals.append((left, right))
    return tuple(intervals)


def centred_ae_radius(a: int, b: int) -> Fraction:
    """Right radius of the a.e.-connected union component at zero."""

    intervals = sorted(half_circle_teeth(a) + half_circle_teeth(b))
    require(intervals and intervals[0][0] == 0, "zero tooth disappeared")
    right = Fraction(0)
    for left, endpoint in intervals:
        if left > right:
            break
        right = max(right, endpoint)
    return right


def analytic_bound(a: int, b: int) -> Fraction:
    a, b = sorted((a, b))
    if b < 13 * a:
        return Fraction(1, 14 * a)
    return Fraction(1, 14 * a) + Fraction(1, 7 * b)


# 1. Centre-mask formula and the proper-mask obstruction.
center_orders_checked = 0
for order in range(1, 1001):
    direct = sum(
        min(root, order - root) * 14 <= order
        for root in range(order)
    )
    formula = 2 * (order // 14) + 1
    require(direct == formula, f"centre formula failed at order {order}")
    if order >= 3:
        require(Fraction(formula, order) <= Fraction(1, 3), "one-third bound failed")
    center_orders_checked += 1

proper_order_pairs_checked = 0
for first in range(2, 1001):
    first_count = 2 * (first // 14) + 1
    for second in range(2, 1001):
        second_count = 2 * (second // 14) + 1
        if first == second == 2:
            union_bound = Fraction(1, 2)
        else:
            union_bound = Fraction(first_count, first) + Fraction(second_count, second)
        require(union_bound < 1, "two proper centre masks reached full density")
        proper_order_pairs_checked += 1

for n in range(1, 121):
    for u in range(1, 121):
        order = n // gcd(u, n)
        require(
            closed_center_count(u, n)
            == gcd(u, n) * (2 * (order // 14) + 1),
            "gcd centre count drifted",
        )


# 2. Exact reduced radii and equality.
reduced_pairs_checked = 0
sharp_pairs: list[tuple[int, int]] = []
second_radius = Fraction(0)
second_pairs: list[tuple[int, int]] = []
for a in range(1, 121):
    for b in range(a, 121):
        radius = centred_ae_radius(a, b)
        bound = analytic_bound(a, b)
        require(radius <= bound, f"analytic radius bound failed at {(a, b)}")
        require(bound <= SHARP, f"universal sharp bound failed at {(a, b)}")
        if radius == SHARP:
            sharp_pairs.append((a, b))
        elif radius > second_radius:
            second_radius = radius
            second_pairs = [(a, b)]
        elif radius == second_radius:
            second_pairs.append((a, b))
        reduced_pairs_checked += 1

require(sharp_pairs == [(1, 13)], "sharp reduced pair changed")
require(centred_ae_radius(1, 13) == SHARP, "sharp seam radius changed")
require(centred_ae_radius(1, 13) < TARGET, "target ceased to be strict")


# 3. Bounded scaled census, using the proved divisibility gate followed by
# an independent exact interval merge on the reduced pair.
triple_checks = 0
target_covers: list[tuple[int, int, int]] = []
sharp_scaled_triples: list[tuple[int, int, int]] = []
for n in range(1, 121):
    for a_speed in range(1, 121):
        for b_speed in range(a_speed + 1, 121):
            triple_checks += 1
            if a_speed % n or b_speed % n:
                continue
            radius = centred_ae_radius(a_speed // n, b_speed // n)
            if radius >= TARGET:
                target_covers.append((a_speed, b_speed, n))
            if radius >= SHARP:
                sharp_scaled_triples.append((a_speed, b_speed, n))

expected_scaled = [(n, 13 * n, n) for n in range(1, 10)]
require(triple_checks == 856800, "bounded triple universe changed")
require(target_covers == [], "a bounded target cover appeared")
require(sharp_scaled_triples == expected_scaled, "scaled equality census changed")


# 4. The load-bearing open-seam hostile from THM-2434.
seam = Fraction(5, 14)
epsilon = Fraction(1, 46200)
require(not in_danger(3, seam), "speed 3 caught the hostile seam")
require(not in_danger(11, seam), "speed 11 caught the hostile seam")
left_owners = tuple(u for u in (3, 11) if in_danger(u, seam - epsilon))
right_owners = tuple(u for u in (3, 11) if in_danger(u, seam + epsilon))
require(left_owners == (3,), "left seam owner changed")
require(right_owners == (11,), "right seam owner changed")


# 5. Explicit sharp seam and a strict point beyond it.
sharp_seam = Fraction(13, 182)
beyond = Fraction(31, 364)  # midpoint of 15/182 and 16/182
require(not in_danger(1, sharp_seam), "D1 caught the sharp seam")
require(not in_danger(13, sharp_seam), "D13 caught the sharp seam")
require(in_danger(1, sharp_seam - epsilon), "D1 lost its sharp-side tooth")
require(in_danger(13, sharp_seam + epsilon), "D13 lost its sharp-side tooth")
require(not in_danger(1, beyond), "D1 reached beyond the sharp radius")
require(not in_danger(13, beyond), "D13 reached beyond the sharp radius")


print("THM-2440 exact companion")
print(f"center_orders_checked={center_orders_checked}")
print(f"proper_order_pairs_checked={proper_order_pairs_checked}")
print(f"reduced_pairs_checked={reduced_pairs_checked}")
print(f"sharp_radius={SHARP}")
print("sharp_reduced_pairs=" + ";".join(f"{a},{b}" for a, b in sharp_pairs))
print(f"second_radius_through_120={second_radius}")
print(
    "second_pairs_through_120="
    + ";".join(f"{a},{b}" for a, b in second_pairs)
)
print(f"bounded_triple_checks={triple_checks}")
print(f"target_radius={TARGET}")
print(f"target_covers={len(target_covers)}")
print(f"sharp_scaled_triples={len(sharp_scaled_triples)}")
print("sharp_scaled_n=" + ",".join(str(n) for _, _, n in sharp_scaled_triples))
print(f"hostile_seam={seam}")
print("hostile_seam_owners=left:3,right:11,point:none")
print(f"sharp_internal_seam={sharp_seam}")
print(f"strict_beyond_point={beyond}")
print("ALL CHECKS PASSED")
