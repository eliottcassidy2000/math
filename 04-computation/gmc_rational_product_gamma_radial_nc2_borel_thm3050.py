#!/usr/bin/env python3
"""Exact rational referee for THM-3050."""

from fractions import Fraction
from math import factorial
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jdump(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def rising(alpha, n):
    out = Fraction(1)
    for j in range(n):
        out *= alpha + j
    return out


def weight(scale, shapes, n):
    out = scale ** n
    for alpha in shapes:
        out *= rising(alpha, n)
    return out


def vp_integer(n, p):
    n = abs(n)
    if n == 0:
        return 10**9
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def vp_fraction(x, p):
    return vp_integer(x.numerator, p) - vp_integer(x.denominator, p)


def harmonic(k):
    return sum(Fraction(1, j) for j in range(1, k + 1))


def character(k):
    kfac = factorial(k)
    a = kfac * (harmonic(k) - 1)
    b = kfac * (k + 1 - 2 * harmonic(k))
    require(a.denominator == b.denominator == 1, "character not integral")
    return int(a), int(b), int(a + b)


def width_flag(k, t, n):
    _, b, interior = character(k)
    out = Fraction(1)
    for s in range(1, n):
        out *= (1 + s * t) ** interior
    if n:
        out *= (1 + n * t) ** b
    return out


print("THM-3050 RATIONAL PRODUCT-GAMMA RADIAL NC2 AND BOREL ORDER")

families = (
    (Fraction(2, 3), (Fraction(1, 2),)),
    (Fraction(3, 5), (Fraction(1, 2), Fraction(2, 3))),
    (Fraction(2), (Fraction(1, 3), Fraction(5, 4), Fraction(7, 6))),
    (Fraction(7, 11), (Fraction(2, 5), Fraction(4, 7), Fraction(9, 8), Fraction(11, 6))),
)

prime_block_cells = 0
strict_block_cells = 0
prime_block_digest = []
for scale, shapes in families:
    bad_denominators = [scale.denominator, scale.numerator] + [a.denominator for a in shapes]
    for p in (5, 7, 11, 13, 17):
        if any(value % p == 0 for value in bad_denominators):
            continue
        for a0 in range(4):
            base = weight(scale, shapes, p * a0)
            for n in range(p * a0, p * (a0 + 3) + 1):
                require(vp_fraction(weight(scale, shapes, n) / base, p) >= 0, "side ratio not p-integral")
                prime_block_cells += 1
            for aprime in range(a0 + 1, a0 + 4):
                valuation = vp_fraction(weight(scale, shapes, p * aprime) / base, p)
                lower = len(shapes) * (aprime - a0)
                require(valuation >= lower, "product-Gamma strict block lost divisibility")
                strict_block_cells += 1
                if len(prime_block_digest) < 12:
                    prime_block_digest.append({"A0": a0, "Aprime": aprime, "J": len(shapes), "p": p, "v": valuation})
print(jdump({"prime_block_cells": prime_block_cells, "strict_block_cells": strict_block_cells, "strict_digest": prime_block_digest}))

borel_cells = 0
borel_digest = []
for scale, shapes in families:
    order = len(shapes)
    for n in range(13):
        coeff = weight(scale, shapes, n) / (factorial(n) ** order)
        next_coeff = weight(scale, shapes, n + 1) / (factorial(n + 1) ** order)
        ratio = scale
        for alpha in shapes:
            ratio *= alpha + n
        ratio /= (n + 1) ** order
        require(next_coeff == coeff * ratio, "critical hypergeometric recurrence failed")
        borel_cells += 1
        if len(borel_digest) < 8:
            borel_digest.append({"J": order, "c": str(scale), "coefficient": str(coeff), "n": n})
print(jdump({"borel_critical_cells": borel_cells, "digest": borel_digest, "radius_law":"0 if sigma<J; 1/c if sigma=J; infinity if sigma>J"}))

corner_cells = 0
for k in range(2, 6):
    a_count, b_count, interior = character(k)
    for t in (Fraction(1, 2), Fraction(2, 3), Fraction(3, 2)):
        shapes = (Fraction(1, 1) / t,) * a_count + (Fraction(1, 1) / t + 1,) * b_count
        scale = t ** interior
        for n in range(7):
            require(weight(scale, shapes, n) == width_flag(k, t, n), "THM-3047 specialization failed")
            corner_cells += 1
print(jdump({"thm3047_specialization_cells": corner_cells, "k":"2..5", "rational_t":["1/2","2/3","3/2"]}))

controls = []
for scale, shapes in families:
    w1 = weight(scale, shapes, 1)
    require(2 * w1 > 0, "balanced positive control failed")
    controls.append({"J": len(shapes), "E[(Z+Zbar)^2]": str(2 * w1), "one_sided_Z3_all_moments": 0})
print(jdump({"nullcone_controls": controls, "ordinary_exp_hostile":"P=Z^d,d>=3: polynomial moments vanish but absolute exponential moment diverges"}))

print("all_exact_checks=PASS")
