#!/usr/bin/env python3
"""Exact referee for THM-2080's unequal comb-overlap law.

The theorem is analytic; this finite program is an independent arithmetic
audit of its exact formula, small-product equality ledger, and Hunter-star
consequence.  It deliberately uses explicit checks rather than assertions so
that ``python -O`` runs the same referee.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def circle_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def centered_intersection_exact(h: int, q: int) -> F:
    """Exact measure of {||ht||<1/7} intersect {||qt||<1/14}."""
    boundaries = {F(0), F(1)}
    for speed, radius in ((h, F(1, 7)), (q, F(1, 14))):
        for k in range(speed + 1):
            for sign in (-1, 1):
                point = F(k, speed) + sign * radius / speed
                if 0 <= point <= 1:
                    boundaries.add(point)
    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if (
            circle_distance(h * midpoint) < F(1, 7)
            and circle_distance(q * midpoint) < F(1, 14)
        ):
            measure += right - left
    return measure


def overlap_formula(h: int, q: int) -> F:
    g = gcd(h, q)
    a, b = h // g, q // g
    x, y = F(a % 14, 14), F(b % 7, 7)
    correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
    return F(2, 49) + F(2, a * b) * correction


def comb_union_measure(speeds_and_radii: tuple[tuple[int, F], ...]) -> F:
    boundaries = {F(0), F(1)}
    for speed, radius in speeds_and_radii:
        for k in range(speed + 1):
            for sign in (-1, 1):
                point = F(k, speed) + sign * radius / speed
                if 0 <= point <= 1:
                    boundaries.add(point)
    points = sorted(boundaries)
    measure = F(0)
    for left, right in zip(points, points[1:]):
        midpoint = (left + right) / 2
        if any(circle_distance(v * midpoint) < radius for v, radius in speeds_and_radii):
            measure += right - left
    return measure


def main() -> None:
    print("THM-2080 UNEQUAL COMB OVERLAP REFEREE")

    pair_checks = 0
    equality_pairs: list[tuple[int, int]] = []
    minimum = (F(1), 0, 0)
    for h in range(1, 64, 2):
        for q in range(1, 127):
            direct = centered_intersection_exact(h, q)
            formula = overlap_formula(h, q)
            require(direct == formula, f"pair formula mismatch h={h}, q={q}")
            require(formula >= F(1, 42), f"floor failure h={h}, q={q}: {formula}")
            require(
                (formula == F(1, 42)) == (q == 6 * h),
                f"equality classification failure h={h}, q={q}: {formula}",
            )
            if formula == F(1, 42):
                equality_pairs.append((h, q))
            if (formula, h, q) < minimum:
                minimum = (formula, h, q)
            pair_checks += 1

    negative_small: list[tuple[int, int, F]] = []
    small_checks = 0
    for a in range(1, 15):
        for b in range(1, 15):
            if gcd(a, b) != 1 or a * b > 14:
                continue
            x, y = F(a % 14, 14), F(b % 7, 7)
            correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
            value = F(2, 49) + F(2, a * b) * correction
            require(value >= F(1, 42), f"small table floor failure {(a, b)}")
            if correction < 0:
                negative_small.append((a, b, value))
            small_checks += 1
    require(
        [(a, b) for a, b, value in negative_small if value == F(1, 42)]
        == [(1, 6), (12, 1)],
        "small-product equality ledger changed",
    )

    hunter_bounds = {s: F(12 + 5 * s, 42) for s in range(1, 7)}
    require(all(hunter_bounds[s] < 1 for s in range(1, 6)), "Hunter s<=5 arithmetic")
    require(hunter_bounds[6] == 1, "Hunter six-event boundary")

    negative_residue_pairs: list[tuple[int, int]] = []
    neutral_residue_pairs: list[tuple[int, int]] = []
    for a_mod in range(1, 14, 2):
        for b_mod in range(7):
            x, y = F(a_mod, 14), F(b_mod, 7)
            correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
            predicted_negative = (
                a_mod in {1, 3, 5} and b_mod in {4, 5, 6}
            ) or (
                a_mod in {9, 11, 13} and b_mod in {1, 2, 3}
            )
            require((correction < 0) == predicted_negative, "depth-four sign router")
            require(
                (correction == 0) == (a_mod == 7 or b_mod == 0),
                "depth-four neutral router",
            )
            if correction < 0:
                negative_residue_pairs.append((a_mod, b_mod))
            if correction == 0:
                neutral_residue_pairs.append((a_mod, b_mod))

    depth_four_control = (4, 5, 6, 11, 12, 13, 360360)
    depth_four_invoice = F(0)
    for q in depth_four_control:
        a, b = 1, q
        x, y = F(a % 14, 14), F(b % 7, 7)
        correction = min(x, y) + max(F(0), x + y - 1) - 2 * x * y
        depth_four_invoice += correction / (a * b)
    require(depth_four_invoice == F(-5167, 210210), "depth-four hostile invoice")

    # Hostile direct-union audit: every six-subset through 12 and every odd
    # guard through 23 leaves positive complement.  This is not used as proof.
    containment_checks = 0
    smallest_complement = (F(1), (), 0)
    for Q in combinations(range(1, 13), 6):
        for h in range(1, 24, 2):
            union = comb_union_measure(
                ((h, F(1, 7)),) + tuple((q, F(1, 14)) for q in Q)
            )
            complement = 1 - union
            require(complement > 0, f"hostile containment witness Q={Q}, h={h}")
            if (complement, Q, h) < smallest_complement:
                smallest_complement = (complement, Q, h)
            containment_checks += 1

    print(f"exact pair/formula checks: {pair_checks}")
    print(f"minimum pair overlap: {minimum[0]} at h={minimum[1]}, q={minimum[2]}")
    print(f"equality pairs in audit range: {equality_pairs}")
    print(f"coprime small-product checks: {small_checks}")
    print("negative small-product table:")
    for a, b, value in negative_small:
        print(f"  ({a:2d},{b:2d}) -> {value}")
    print(f"Hunter bounds s=1..6: {hunter_bounds}")
    print(f"depth-four negative residue pairs: {negative_residue_pairs}")
    print(f"depth-four neutral residue pairs: {neutral_residue_pairs}")
    print(
        f"depth-four hostile invoice: {depth_four_invoice} "
        f"at Q={depth_four_control}, h=1"
    )
    print(f"hostile direct-containment checks: {containment_checks}")
    print(
        "smallest hostile complement: "
        f"{smallest_complement[0]} at Q={smallest_complement[1]}, h={smallest_complement[2]}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
