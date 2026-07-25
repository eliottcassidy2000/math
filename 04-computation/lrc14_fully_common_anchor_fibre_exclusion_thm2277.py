#!/usr/bin/env python3
"""Exact referee for THM-2277.

The proof is structural.  This companion checks the finite residue-class
certificate behind the all-n fibre cap, the complete height-20 ray packet,
and the exact THM-2137 boundary-complexity margins.
"""

from fractions import Fraction
from itertools import combinations
from math import factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_seventh(n: int) -> int:
    if n < 1:
        raise ValueError("n must be positive")
    return (n + 6) // 7


def verify_all_n_fibre_caps() -> None:
    """Verify the cap by the seven residue classes of n, not by a cutoff."""

    # The q=0 representatives n=1,...,6 give every exceptional equality.
    small = tuple(Fraction(ceil_seventh(n), n) for n in range(1, 7))
    require(
        small
        == (
            Fraction(1, 1),
            Fraction(1, 2),
            Fraction(1, 3),
            Fraction(1, 4),
            Fraction(1, 5),
            Fraction(1, 6),
        ),
        "small fibre-ratio table changed",
    )

    # Write n=7q+r with q>=1 and 0<=r<=6.  If r=0, ceil(n/7)=q;
    # otherwise it is q+1.  The two slacks below are affine increasing
    # functions of q, so checking q=1 proves the claims for every q>=1.
    for residue in range(7):
        extra = int(residue != 0)
        half_slope = 7 - 2
        half_intercept = residue - 2 * extra
        third_slope = 7 - 3
        third_intercept = residue - 3 * extra
        require(
            half_slope > 0 and half_slope + half_intercept > 0,
            f"half-cap residue {residue} failed",
        )
        require(
            third_slope > 0 and third_slope + third_intercept > 0,
            f"third-cap residue {residue} failed",
        )

    # Therefore ceil(n/7)/n <=1/2 for n>=2, with equality only n=2,
    # and <=1/3 for n>=3, with equality only n=3.
    require(2 * ceil_seventh(2) == 2, "n=2 equality missing")
    require(3 * ceil_seventh(3) == 3, "n=3 equality missing")


def fraction_ceiling(value: Fraction) -> int:
    return (value.numerator + value.denominator - 1) // value.denominator


verify_all_n_fibre_caps()

guard_rays = tuple(a for a in range(1, 21) if a % 2 == 1 and a % 13 != 0)
unit_rays = tuple(b for b in range(1, 21) if b % 13 != 0)

require(guard_rays == (1, 3, 5, 7, 9, 11, 15, 17, 19), "guard rays changed")
require(
    unit_rays
    == (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20),
    "unit rays changed",
)

unit_packets = tuple(combinations(unit_rays, 5))
unlabelled_packets = len(guard_rays) * len(unit_packets)
labelled_packets = unlabelled_packets * factorial(5)
boundary_budgets = tuple(
    a + 13 + sum(packet) for a in guard_rays for packet in unit_packets
)

require(unlabelled_packets == 104_652, "unlabelled packet count changed")
require(labelled_packets == 12_558_240, "labelled packet count changed")
require(max(boundary_budgets) == 122, "boundary-budget maximum changed")

private_mass = Fraction(2593, 90090)
deep_factor = 13**3
one_comb_sigma = Fraction(1, 7)
two_comb_sigma_cap = Fraction(25, 91)

one_comb_floor = deep_factor * private_mass / one_comb_sigma
two_comb_floor = deep_factor * private_mass / two_comb_sigma_cap
one_comb_margin = one_comb_floor - 122
two_comb_margin = two_comb_floor - 122

require(one_comb_floor == Fraction(438217, 990), "one-comb floor changed")
require(one_comb_margin == Fraction(317437, 990), "one-comb margin changed")
require(two_comb_floor == Fraction(5696821, 24750), "two-comb floor changed")
require(two_comb_margin == Fraction(2677321, 24750), "two-comb margin changed")
require(one_comb_margin > 0 and two_comb_margin > 0, "boundary contradiction failed")
require(fraction_ceiling(one_comb_floor) == 443, "one-comb integer floor changed")
require(fraction_ceiling(two_comb_floor) == 231, "two-comb integer floor changed")

print("THM-2277 exact referee")
print("fibre_ratio: n=1 -> 1")
print("fibre_ratio: n>=2 -> <=1/2, equality iff n=2")
print("fibre_ratio: n>=3 -> <=1/3, equality iff n=3")
print(f"guard_rays={len(guard_rays)} unit_rays={len(unit_rays)}")
print(f"unlabelled_multiplier_packets={unlabelled_packets}")
print(f"labelled_multiplier_packets={labelled_packets}")
print(f"max_boundary_budget={max(boundary_budgets)}")
print(f"private_mass={private_mass} deep_factor={deep_factor}")
print(
    f"one_comb_floor={one_comb_floor} "
    f"integer_requirement={fraction_ceiling(one_comb_floor)} "
    f"margin_over_122={one_comb_margin}"
)
print(
    f"two_comb_sigma_cap={two_comb_sigma_cap} "
    f"two_comb_floor={two_comb_floor} "
    f"integer_requirement={fraction_ceiling(two_comb_floor)} "
    f"margin_over_122={two_comb_margin}"
)
print("verdict=PASS")
