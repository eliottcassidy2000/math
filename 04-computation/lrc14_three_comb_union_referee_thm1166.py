#!/usr/bin/env python3
"""Exact referee for the sharp three-comb union addendum to THM-1166.

For D_s={t in R/Z: ||st||<1/14}, the addendum proves

    measure(D_a union D_b union D_c) <= 36/91

for all distinct positive integer speeds, with equality only for a scaled
and permuted copy of (1,12,13).

The infinite tail is discharged by THM-1166's folded pair-overlap bound.
Only configurations whose two smallest reduced ratio-products are at most
63 remain.  This script evaluates that complete rational bank and checks
triple intersections independently by interval intersection and by a full
breakpoint arrangement.  It also checks the scalar 5+3 invoice arithmetic
used in the THM-2168 addendum.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


PAIR_BASE = F(1, 49)
TREE_TARGET = F(3, 91)
UNION_TARGET = F(36, 91)
DELTA_5 = F(961, 6930)
INVOICE_CONSTANT = F(12493, 35640)
PRODUCT_CUTOFF = 63


def require(condition: bool, message: str) -> None:
    """Optimization-stable certificate check."""
    if not condition:
        raise RuntimeError(message)


def fold(residue: int, modulus: int) -> int:
    residue %= modulus
    return residue * (modulus - residue)


def rho_fold(a: int, b: int) -> F:
    """THM-1166's exact measure of D_a intersect D_b."""
    require(a > 0 and b > 0 and a != b, "pair speeds must be distinct positive")
    if a > b:
        a, b = b, a
    common = gcd(a, b)
    modulus = 14 * common
    numerator = (
        4 * a * b
        + fold(a + b, modulus)
        - fold(b - a, modulus)
    )
    return F(numerator, 196 * a * b)


def oriented_reduced_ratios(product_cap: int) -> list[F]:
    """All oriented nonunit reduced ratios with numerator*denominator <= cap."""
    return sorted(
        {
            F(numerator, denominator)
            for numerator in range(1, product_cap + 1)
            for denominator in range(1, product_cap + 1)
            if numerator != denominator
            and numerator * denominator <= product_cap
            and gcd(numerator, denominator) == 1
        }
    )


def danger_intervals(speed: int) -> list[tuple[F, F]]:
    """The components of D_speed in [0,1], with endpoint choices immaterial."""
    radius = F(1, 14 * speed)
    pieces: list[tuple[F, F]] = []
    for index in range(speed):
        center = F(index, speed)
        left, right = center - radius, center + radius
        if left < 0:
            pieces.extend(((F(0), right), (F(1) + left, F(1))))
        elif right > 1:
            pieces.extend(((left, F(1)), (F(0), right - F(1))))
        else:
            pieces.append((left, right))
    pieces.sort()
    merged: list[tuple[F, F]] = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return merged


def intersect_intervals(
    lefts: list[tuple[F, F]], rights: list[tuple[F, F]]
) -> list[tuple[F, F]]:
    """Exact intersection of two sorted interval unions."""
    out: list[tuple[F, F]] = []
    i = j = 0
    while i < len(lefts) and j < len(rights):
        left = max(lefts[i][0], rights[j][0])
        right = min(lefts[i][1], rights[j][1])
        if left < right:
            out.append((left, right))
        if lefts[i][1] < rights[j][1]:
            i += 1
        else:
            j += 1
    return out


def interval_measure(parts: list[tuple[F, F]]) -> F:
    return sum((right - left for left, right in parts), F(0))


def triple_mass_intervals(speeds: tuple[int, int, int]) -> F:
    """First exact route: intersect the three individual comb decompositions."""
    common = danger_intervals(speeds[0])
    for speed in speeds[1:]:
        common = intersect_intervals(common, danger_intervals(speed))
    return interval_measure(common)


def circular_norm(value: F) -> F:
    fractional = value - (value.numerator // value.denominator)
    return min(fractional, 1 - fractional)


def triple_mass_arrangement(speeds: tuple[int, int, int]) -> F:
    """Independent exact route: classify the complete endpoint arrangement."""
    breakpoints = {F(0), F(1)}
    for speed in speeds:
        for index in range(speed + 1):
            for sign in (-1, 1):
                point = F(index, speed) + sign * F(1, 14 * speed)
                if 0 <= point <= 1:
                    breakpoints.add(point)
    ordered = sorted(breakpoints)
    total = F(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(circular_norm(speed * midpoint) < F(1, 14) for speed in speeds):
            total += right - left
    return total


def primitive_triple(first: F, second: F) -> tuple[int, int, int]:
    """Normalize the rational speed triple (1,first,second)."""
    scale = lcm(first.denominator, second.denominator)
    speeds = [
        scale,
        first.numerator * (scale // first.denominator),
        second.numerator * (scale // second.denominator),
    ]
    common = gcd(gcd(speeds[0], speeds[1]), speeds[2])
    return tuple(sorted(speed // common for speed in speeds))


EXPECTED_NET = {
    (1, 12, 13): F(3, 91),
    (1, 12, 27): F(29, 756),
    (1, 12, 40): F(17, 420),
    (1, 12, 66): F(13, 308),
    (1, 13, 169): F(50, 1183),
    (4, 9, 108): F(8, 189),
    (1, 10, 55): F(3, 70),
    (3, 10, 120): F(3, 70),
    (1, 12, 156): F(47, 1092),
    (1, 13, 156): F(47, 1092),
    (12, 13, 156): F(47, 1092),
    (2, 11, 132): F(10, 231),
    (1, 12, 144): F(11, 252),
    (2, 11, 110): F(17, 385),
}


def main() -> None:
    print("THM-1166 sharp three-comb union referee")

    tail_margin = (
        2 * (PAIR_BASE - F(1, 4 * (PRODUCT_CUTOFF + 1))) - TREE_TARGET
    )
    require(tail_margin == F(3, 81536), "wrong product-tail margin")
    print("analytic reduction: p_2>=64 gives tree excess 3/81536")

    ratios = oriented_reduced_ratios(PRODUCT_CUTOFF)
    require(len(ratios) == 208, "wrong oriented ratio-bank size")

    candidate_rows: dict[
        tuple[int, int, int], tuple[F, F, tuple[F, F, F]]
    ] = {}
    checked = 0
    for first in ratios:
        for second in ratios:
            if first == second:
                continue
            checked += 1
            quotient = first / second
            weights = (
                rho_fold(first.numerator, first.denominator),
                rho_fold(second.numerator, second.denominator),
                rho_fold(quotient.numerator, quotient.denominator),
            )
            tree_weight = sum(sorted(weights)[-2:], F(0))
            if tree_weight > TREE_TARGET:
                continue
            speeds = primitive_triple(first, second)
            if speeds in candidate_rows:
                continue

            triple_a = triple_mass_intervals(speeds)
            triple_b = triple_mass_arrangement(speeds)
            require(triple_a == triple_b, f"triple-route mismatch at {speeds}")

            direct_pair_weights = tuple(
                rho_fold(a, b) for a, b in combinations(speeds, 2)
            )
            for a, b in combinations(speeds, 2):
                interval_pair = interval_measure(
                    intersect_intervals(danger_intervals(a), danger_intervals(b))
                )
                require(
                    rho_fold(a, b) == interval_pair,
                    f"folded/interval pair mismatch at {(a, b)}",
                )
            net_redundancy = sum(direct_pair_weights, F(0)) - triple_a
            candidate_rows[speeds] = (
                net_redundancy,
                triple_a,
                direct_pair_weights,
            )

    require(checked == 43056, "wrong ordered configuration count")
    require(len(candidate_rows) == 14, "wrong low-tree primitive census")
    require(
        {speeds: row[0] for speeds, row in candidate_rows.items()} == EXPECTED_NET,
        "low-tree net-redundancy table mismatch",
    )

    minimum = min(row[0] for row in candidate_rows.values())
    minimizers = {
        speeds for speeds, row in candidate_rows.items() if row[0] == minimum
    }
    require(minimum == TREE_TARGET, "wrong global net-redundancy floor")
    require(minimizers == {(1, 12, 13)}, "wrong equality triple")
    require(F(3, 7) - minimum == UNION_TARGET, "wrong union constant")

    print(f"oriented reduced ratios with product<=63: {len(ratios)}")
    print(f"ordered incident-ratio configurations: {checked}")
    print(f"primitive triples with tree weight<=3/91: {len(candidate_rows)}")
    print("exact interval/breakpoint and folded/interval cross-checks: PASS")
    print("minimum net pair-minus-triple redundancy: 3/91 at (1,12,13)")
    print("next net redundancy: 29/756 at (1,12,27)")
    print("global three-comb union maximum: 36/91")
    print("equality: scaled and permuted copies of (1,12,13)")

    require(
        DELTA_5 / UNION_TARGET == INVOICE_CONSTANT,
        "wrong scalar 5+3 invoice constant",
    )
    scaled_invoices = [
        INVOICE_CONSTANT * 13**depth for depth in range(1, 7)
    ]
    floors = [
        (value.numerator + value.denominator - 1) // value.denominator
        for value in scaled_invoices
    ]
    require(
        floors == [5, 60, 771, 10012, 130151, 1691957],
        "wrong scalar 5+3 integer floors",
    )
    print("scalar 5+3 invoice constant: 12493/35640")
    print("terminal depths 1..6 force B floors: 5,60,771,10012,130151,1691957")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"source_sha256={source_hash}")


if __name__ == "__main__":
    main()
