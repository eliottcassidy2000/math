#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3743.

The cited input is Khinchin's flatness theorem.  This file verifies the full
quotient-width algebra, the explicit general numerical specialization, the
Graver split, and hostile comparisons with existing LRC(14) relation
mechanisms.  It does not attempt to computationally prove the flatness theorem
itself.
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

# Averkov--Hofscheier--Nill record the explicit classical inequality
# Flt(d)<=sqrt((d+1)(2d+1)/6)*d^(3/2).  At d=12 this squares to 93600;
# after multiplying by 7/6 the relation bound is 70*sqrt(26).
flatness_bound_squared = Fraction((dimension + 1) * (2 * dimension + 1), 6) * dimension**3
relation_bound_squared = Fraction(7, 6) ** 2 * flatness_bound_squared
require(flatness_bound_squared == 93600, "dimension-twelve flatness square")
require(relation_bound_squared == 127400, "dimension-twelve relation square")
explicit_l1_cap = isqrt(relation_bound_squared.numerator // relation_bound_squared.denominator)
require(explicit_l1_cap == 356, "floor(70 sqrt(26))=356")
require(explicit_l1_cap**2 < relation_bound_squared < (explicit_l1_cap + 1) ** 2, "sharp integer floor")
checks += 4
print("Using Flt(d)<=sqrt((d+1)(2d+1)/6)*d^(3/2):")
print("  Flt(12)<=60*sqrt(26), so ||a||_1<=70*sqrt(26)<357.")
print("Exact integer consequence: ||a||_1<=356.")


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
require(explicit_l1_cap < sum(selberg_caps), "flatness improves the existing total cap")
checks += 2
print("AP {1,...,13} and {1,...,12,5460} already have (2,-1,0,...), l1=3.")
print("THM-2144 forces a relation with l1<=3*29+10*28=367 on every zero-safe row.")
print("The flatness cap 356 improves that total-coefficient number by 11; its geometric direction is different.")
print("THM-2051 instead forces support 3..5 with coefficient height <=2^20: sparse but much taller.")


print("\n=== Minimal-width Graver split ===")
coprime_pair_ratios = sum(
    1
    for first in range(1, explicit_l1_cap)
    for second in range(first + 1, explicit_l1_cap + 1 - first)
    if gcd(first, second) == 1
)
require(coprime_pair_ratios == 19314, "bounded coprime pair-ratio atlas")
triple_speeds = (3, 4, 5)
triple_relation = (1, -2, 1)
pair_norms = []
for first in range(3):
    for second in range(first + 1, 3):
        common = gcd(triple_speeds[first], triple_speeds[second])
        pair_norms.append((triple_speeds[first] + triple_speeds[second]) // common)
require(dot(triple_relation, triple_speeds) == 0, "genuine triple relation")
require(l1(triple_relation) == 4 < min(pair_norms), "triple can beat every pair relation")
checks += 3
print("An l1-minimal relation is conformally indecomposable, hence a Graver element.")
print("Support two gives a reduced speed ratio with numerator+denominator<=356;")
print(f"there are {coprime_pair_ratios} unordered distinct coprime ratios in that atlas.")
print("The higher branch is genuine: (1,-2,1) on speeds (3,4,5) has l1=4 and beats every pair relation.")


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
