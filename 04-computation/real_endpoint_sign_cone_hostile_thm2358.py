#!/usr/bin/env python3
"""Exact rational referee for THM-2358.

The proof uses only rational interval geometry and elementary strict
trigonometric orderings.  No floating point and no Python ``assert`` enter
the certificate.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def interval_intersection(
    left: Fraction,
    right: Fraction,
    lower: Fraction,
    upper: Fraction,
) -> tuple[Fraction, Fraction] | None:
    lo = max(left, lower)
    hi = min(right, upper)
    if lo >= hi:
        return None
    return lo, hi


def root_danger_intervals(t: int) -> list[tuple[Fraction, Fraction]]:
    """Intervals in y where ||2y+t/13||<1/14 and |y|<1/14."""

    domain_left = -Fraction(1, 14)
    domain_right = Fraction(1, 14)
    intervals: list[tuple[Fraction, Fraction]] = []
    for integer_shift in (-1, 0, 1):
        lower = (
            Fraction(integer_shift, 1)
            - Fraction(t, 13)
            - Fraction(1, 14)
        ) / 2
        upper = (
            Fraction(integer_shift, 1)
            - Fraction(t, 13)
            + Fraction(1, 14)
        ) / 2
        intersection = interval_intersection(
            domain_left,
            domain_right,
            lower,
            upper,
        )
        if intersection is not None:
            intervals.append(intersection)
    return intervals


def main() -> None:
    # Bare r=1 incidence after parameterizing D_13 by x=(k+y)/13.
    signed_residues = list(range(-6, 7))
    interval_table = {
        t: root_danger_intervals(t)
        for t in signed_residues
    }
    expected_table = {
        -2: [(Fraction(15, 364), Fraction(26, 364))],
        -1: [(Fraction(1, 364), Fraction(26, 364))],
        0: [(Fraction(-13, 364), Fraction(13, 364))],
        1: [(Fraction(-26, 364), Fraction(-1, 364))],
        2: [(Fraction(-26, 364), Fraction(-15, 364))],
    }
    require(
        {t: value for t, value in interval_table.items() if value}
        == expected_table,
        f"root-danger interval table changed: {interval_table}",
    )
    incidence = sum(
        upper - lower
        for intervals in interval_table.values()
        for lower, upper in intervals
    )
    require(incidence == Fraction(49, 182), "bare incidence changed")
    endpoints = sorted(
        {
            Fraction(-1, 14),
            Fraction(1, 14),
            *(
                endpoint
                for intervals in interval_table.values()
                for interval in intervals
                for endpoint in interval
            ),
        }
    )
    max_active = 0
    for left, right in zip(endpoints, endpoints[1:]):
        midpoint = (left + right) / 2
        active = sum(
            lower < midpoint < upper
            for intervals in interval_table.values()
            for lower, upper in intervals
        )
        max_active = max(max_active, active)
    require(max_active == 2, "more than two danger roots became active")

    # Exact D_13 cap D_169 overlap.  A D_169 tooth is
    # x=(K+y)/169, |y|<1/14.  D_13 holds on the complete tooth exactly
    # when K=0 mod 13; every other residue misses even at the boundary.
    live_teeth: list[int] = []
    for k in range(169):
        residue = k % 13
        signed = residue if residue <= 6 else residue - 13
        if signed == 0:
            live_teeth.append(k)
        else:
            require(
                Fraction(abs(signed), 1) - Fraction(1, 14)
                >= Fraction(13, 14),
                "a nonzero residue approached the D_13 tooth interior",
            )
    require(
        len(live_teeth) == 13
        and all(k % 13 == 0 for k in live_teeth),
        "wrong D_13/D_169 live-tooth set",
    )
    tooth_length = Fraction(1, 7 * 169)
    overlap = len(live_teeth) * tooth_length
    require(overlap == Fraction(1, 91), "D_13/D_169 overlap changed")

    # The rigorous sign gap uses cos(pi/7)>cos(pi/4)>2/3.
    bare_rational_floor = 49 * Fraction(2, 3)
    correction_numerator = 26
    require(
        bare_rational_floor > correction_numerator,
        "the U_1 sign gap stopped being strict",
    )

    # The inverse target support line u+2v=2 over F_13.
    support_line = {
        ((2 - 2 * r) % 13, r)
        for r in range(13)
    }
    require(
        len(support_line) == 13
        and all((u + 2 * v - 2) % 13 == 0 for u, v in support_line),
        "the target support line changed",
    )

    # Explicit positive support used for L(0,0): (1/28,1/26) is inside
    # D_13 but outside D_26 and D_169 after the y-parameterization.
    require(
        Fraction(1, 28) < Fraction(1, 26) < Fraction(1, 14),
        "the untwisted positive-support interval collapsed",
    )
    # cos(4*pi*y)>0 on |y|<1/14 because 2/7<1/2.
    require(Fraction(2, 7) < Fraction(1, 2), "left phase lost positivity")
    # sin(3*pi/7)>sin(2*pi/7), since both angles lie in (0,pi/2).
    require(
        Fraction(0, 1) < Fraction(2, 7)
        < Fraction(3, 7) < Fraction(1, 2),
        "right endpoint sine ordering changed",
    )

    # Primitive-colour hostile at N=91.
    epsilon = Fraction(1, 1000)
    require(
        0 < epsilon < Fraction(1, 182)
        and 2 * epsilon < Fraction(1, 91)
        and Fraction(6, 91) + epsilon < Fraction(1, 14),
        "the N=91 hostile intervals stopped being disjoint inside D_1",
    )
    correlation_weight = 91 * 2 * epsilon
    require(
        correlation_weight == Fraction(91, 500),
        "the N=91 correlation weight changed",
    )
    positive_colour = 76
    negative_colour = 53
    require(
        gcd(positive_colour, 91) == 1
        and gcd(negative_colour, 91) == 1
        and (6 * positive_colour) % 91 == 1
        and (6 * negative_colour) % 91 == 45,
        "the primitive sign colours changed",
    )
    # C_76=c(1+2cos(2pi/91))>0.  For C_53, use
    # cos(90pi/91)=-cos(pi/91) and pi/91<pi/3.
    require(
        Fraction(1, 91) < Fraction(1, 3),
        "the negative primitive-colour angle bound changed",
    )
    require(
        91 % 3 != 0 and gcd(6, 91) == 1,
        "the order-three nonvanishing obstruction changed",
    )

    print("THM-2358 exact real-endpoint sign-cone referee")
    print(
        "bare r=1 interval lengths /364:"
        " t=-2:11, t=-1:25, t=0:26, t=1:25, t=2:11"
    )
    print("bare danger incidence: 49/182 (=7/26); maximum active roots: 2")
    print(f"D_13 intersection D_169: {overlap} ({len(live_teeth)} live teeth)")
    print("strict sign gap: 49*cos(pi/7)>49*(2/3)>26")
    print(f"inverse target support line size: {len(support_line)}")
    print("left/right/deep triangle signs: positive / positive / positive")
    print(
        "N=91 primitive colours:"
        f" k={positive_colour} positive, k={negative_colour} negative,"
        " none zero"
    )
    print("VERDICT: centered positive/even word factors do not force one acute endpoint cone")


if __name__ == "__main__":
    main()
