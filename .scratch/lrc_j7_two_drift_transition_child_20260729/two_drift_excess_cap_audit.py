#!/usr/bin/env python3
"""Independent exact audit of the k-aligned absolute first-drift caps.

This is a scratch verifier, not a canonical theorem artifact.

For every six-body root E subset {1,...,14}, it reconstructs the literal
carrier twice: once on a common integer ruler and once with a separate
Fraction endpoint sweep.  It then evaluates, for ``2 <= k <= 6`` and
``d=7-k``,

    Delta_E >= eta_k h_E,
    delta_E(z) <= 6 r_E/(49z),
    z_small <= B_E(k) := 6 d r_E/(49 eta_k h_E).

For ``d>=2`` the last comparison is strict because the distinct drift labels
have reciprocal sum strictly below ``d/z_small``.  For ``d=1`` it is
non-strict.  Every global maximum below is nonintegral, so both cases give
the displayed floor as the integral cap.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations

import exact_carrier as X


RULER = 14 * 360_360
ETA = {
    2: F(1, 91),
    3: F(3, 91),
    4: F(51, 1183),
    5: F(88, 1365),
    6: F(22, 273),
}
EXPECTED = {
    2: (F(222_522_300, 32_029), 6_947),
    3: (F(59_339_280, 32_029), 1_852),
    4: (F(578_557_980, 544_493), 1_062),
    5: (F(15_171_975, 32_029), 473),
    6: (F(6_068_790, 32_029), 189),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def merge_integer_intervals(
    rows: list[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    ordered = sorted((left, right) for left, right in rows if left < right)
    if not ordered:
        return ()
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for left, right in ordered[1:]:
        if left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1][1] = right
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def direct_integer_carrier(body: tuple[int, ...]) -> tuple[F, int]:
    """Reconstruct C_E on the common integer endpoint ruler."""
    danger: list[tuple[int, int]] = []
    for speed in body:
        require(RULER % (14 * speed) == 0, "inexact integer ruler")
        radius = RULER // (14 * speed)
        step = RULER // speed
        danger.append((0, radius))
        danger.extend(
            (index * step - radius, index * step + radius)
            for index in range(1, speed)
        )
        danger.append((RULER - radius, RULER))
    union = merge_integer_intervals(danger)
    carrier: list[tuple[int, int]] = []
    cursor = 0
    for left, right in union:
        if cursor < left:
            carrier.append((cursor, left))
        if right > cursor:
            cursor = right
    if cursor < RULER:
        carrier.append((cursor, RULER))
    return (
        F(sum(right - left for left, right in carrier), RULER),
        len(carrier),
    )


def main() -> None:
    rows: list[tuple[F, tuple[int, ...], F, int]] = []
    for body in combinations(range(1, 15), 6):
        mass_integer, components_integer = direct_integer_carrier(body)
        carrier_fraction = X.carrier(body)
        mass_fraction = X.mass(carrier_fraction)
        components_fraction = len(carrier_fraction)
        require(
            (mass_integer, components_integer)
            == (mass_fraction, components_fraction),
            f"{body}: independent carrier reconstructions disagree",
        )
        rows.append(
            (
                F(components_integer, 1) / mass_integer,
                body,
                mass_integer,
                components_integer,
            )
        )

    require(len(rows) == 3_003, "six-body universe changed")
    rows.sort(key=lambda row: (-row[0], row[1]))
    maximum = rows[0]
    require(
        maximum
        == (
            F(3_993_990, 32_029),
            (1, 10, 11, 12, 13, 14),
            F(32_029, 105_105),
            38,
        ),
        "maximum component/mass row changed",
    )
    require(
        sum(row[0] == maximum[0] for row in rows) == 1,
        "maximum component/mass ratio is no longer unique",
    )

    print("LRC14 k-aligned absolute first-drift cap scratch audit")
    print("roots=3003;carrier_reconstructions=integer-ruler+Fraction-endpoint")
    print(
        "maximum_component_mass_ratio="
        f"{maximum[0]};body={maximum[1]};h={maximum[2]};r={maximum[3]}"
    )
    print("maximum_component_mass_ratio_unique=TRUE")
    for k in range(2, 7):
        d = 7 - k
        bounds = tuple(
            (
                F(6 * d, 49) * ratio / ETA[k],
                body,
                mass,
                components,
            )
            for ratio, body, mass, components in rows
        )
        bound = max(row[0] for row in bounds)
        cap = bound.numerator // bound.denominator
        require((bound, cap) == EXPECTED[k], f"k={k}: cap changed")
        require(bound.denominator != 1, f"k={k}: equality convention became live")
        require(
            sum(row[0] == bound for row in bounds) == 1,
            f"k={k}: maximum rational bound is not unique",
        )
        require(
            sum(row[0].numerator // row[0].denominator == cap for row in bounds)
            == 1,
            f"k={k}: maximum integer cap is not unique",
        )
        print(
            f"k={k};d={d};eta={ETA[k]};bound={bound};"
            f"decimal={float(bound):.15f};integral_cap={cap};"
            f"strict={'TRUE' if d >= 2 else 'FALSE'};"
            "maximum_bound_unique=TRUE;maximum_cap_unique=TRUE"
        )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
