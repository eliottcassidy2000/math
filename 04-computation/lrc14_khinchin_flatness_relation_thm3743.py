#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3743.

The cited input is Khinchin's flatness theorem.  This file verifies the full
quotient-width algebra, the explicit d^(5/2) numerical specialization, and
hostile comparisons with existing LRC(14) relation mechanisms.  It does not
attempt to computationally prove the flatness theorem itself.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import gcd, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right, strict=True))


def l1(vector: tuple[int, ...]) -> int:
    return sum(abs(value) for value in vector)


def exact_corner_width(vector: tuple[int, ...]) -> Fraction:
    low = Fraction(1, 14)
    high = Fraction(13, 14)
    values = []
    for corner in product((low, high), repeat=len(vector)):
        values.append(sum(Fraction(coefficient) * coordinate for coefficient, coordinate in zip(vector, corner, strict=True)))
    return max(values) - min(values)


checks = 0


print("=== Quotient-lattice normalization ===")
dimension = 12
box_low = Fraction(1, 14)
box_high = Fraction(13, 14)
box_length = box_high - box_low
require(box_length == Fraction(6, 7), "LRC box side length")
checks += 1
print("B=[1/14,13/14]^13 has side length 6/7; Q=R^13/Rn has dimension 12.")
print("For primitive n, Lambda=pi(Z^13) and Lambda*={a in Z^13:a.n=0}.")


print("\n=== Exhaustive corner checks for the width formula ===")
controls = (
    (1, -1) + (0,) * 11,
    (2, -1) + (0,) * 11,
    (3, -2, 1) + (0,) * 10,
    tuple((-1) ** index * (index + 1) for index in range(13)),
    (29,) * 3 + (-28,) * 10,
)
for vector in controls:
    width = exact_corner_width(vector)
    require(width == Fraction(6, 7) * l1(vector), "box width equals (6/7)||a||_1")
    checks += 1
    print(f"a_l1={l1(vector):3d}  width={width}")
print("The identity follows coordinatewise for every real box point; corners are hostile replays.")


print("\n=== Flatness consequence and explicit classical bound ===")
print("A hypothetical counterexample makes pi(B) a full-dimensional Lambda-free zonotope.")
print("Khinchin flatness therefore selects primitive 0!=a in Lambda* with")
print("  (6/7)||a||_1 <= Flt(12), hence ||a||_1 <= (7/6)Flt(12).")

# The cited explicit inequality Flt(d)<=d^(5/2) gives
# (7/6)*12^(5/2)=336*sqrt(3).  Its floor is checked without floats.
require(581 * 581 < 336 * 336 * 3 < 582 * 582, "floor(336 sqrt(3))=581")
explicit_l1_cap = 581
checks += 1
print("Using the explicit Flt(d)<=d^(5/2) input: floor(336*sqrt(3))=581.")


print("\n=== Existing relation-canon comparison ===")
ap = tuple(range(1, 14))
lifted_ap = tuple(range(1, 13)) + (5460,)
short_relation = (2, -1) + (0,) * 11
require(gcd(*ap) == 1 and gcd(*lifted_ap) == 1, "primitive controls")
require(dot(short_relation, ap) == 0, "AP short relation")
require(dot(short_relation, lifted_ap) == 0, "lifted AP short relation")
require(l1(short_relation) == 3, "short relation l1")
checks += 4

selberg_caps = (29,) * 3 + (28,) * 10
require(sum(selberg_caps) == 367, "THM-2144 total cap")
require(sum(selberg_caps) < explicit_l1_cap, "existing total cap beats simple flatness numerical cap")
checks += 2
print("AP {1,...,13} and {1,...,12,5460} already have (2,-1,0,...), l1=3.")
print("THM-2144 forces a relation with l1<=3*29+10*28=367 on every zero-safe row.")
print("Thus the simple explicit flatness cap 581 is numerically weaker; its geometric direction is different.")
print("THM-2051 instead forces support 3..5 with coefficient height <=2^20: sparse but much taller.")


print("\n=== Restricted-zonotope sidecar ===")
center = (Fraction(1, 2),) * 13
twice_center = tuple(2 * coordinate for coordinate in center)
require(all(value.denominator == 1 for value in twice_center), "center is a two-torsion lattice class")
require(Fraction(3, 7) * 2 == box_length, "centered generator half-length")
checks += 2
print("pi(B) has 13 segment generators in dimension 12, is centrally symmetric,")
print("and its center class has order dividing two because 2*(1/2,...,1/2) is integral.")
print("The useful open constant is flatness for this restricted 13-generator, two-torsion-centered class.")


print("\n=== Semantic boundary ===")
print("The relation direction alone carries no active facets, slice intercept, root, owner, or temporal word.")
print("A recursive LRC gain needs those sidecars; neither relation existence nor small l1 proves loneliness.")
print(f"CHECKS={checks}")
