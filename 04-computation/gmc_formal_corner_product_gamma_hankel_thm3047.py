#!/usr/bin/env python3
"""Exact referee for THM-3047 using only rational arithmetic."""

from fractions import Fraction
from itertools import combinations
from math import comb, factorial
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def jdump(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def character(k):
    kfac = factorial(k)
    harmonic = sum(Fraction(1, r) for r in range(1, k + 1))
    a = kfac * (harmonic - 1)
    b = kfac * (k + 1 - 2 * harmonic)
    interior = a + b
    require(a.denominator == b.denominator == interior.denominator == 1, "character lost integrality")
    return int(a), int(b), int(interior)


def rising(x, m):
    out = Fraction(1)
    for j in range(m):
        out *= x + j
    return out


def width_flag(k, t, m):
    a, b, interior = character(k)
    out = Fraction(1)
    for s in range(1, m):
        out *= (1 + s * t) ** interior
    if m:
        out *= (1 + m * t) ** b
    return out


def gamma_product_moment(k, t, m):
    a_count, b_count, interior = character(k)
    shape = 1 / t
    return (t ** (interior * m)) * (rising(shape, m) ** a_count) * (rising(shape + 1, m) ** b_count)


def ratio(k, t, m):
    a, b, _ = character(k)
    return (1 + m * t) ** a * (1 + (m + 1) * t) ** b


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    det = Fraction(1)
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        value = a[col][col]
        det *= value
        for j in range(col, n):
            a[col][j] /= value
        for row in range(col + 1, n):
            scale = a[row][col]
            if scale:
                for j in range(col, n):
                    a[row][j] -= scale * a[col][j]
    return det


print("THM-3047 FORMAL CORNER PRODUCT-GAMMA WIDTH MOMENTS")

characters = []
for k in range(2, 13):
    a, b, interior = character(k)
    require(a >= 1 and b >= 0 and interior == a + b, "invalid character signs")
    characters.append({"A": a, "B": b, "I": interior, "k": k})
print(jdump({"characters": characters}))

def width_factor_profile(k, m):
    """Exponent ledger for the formal factors (1+s*t), without giant integers."""
    _, b, interior = character(k)
    return {s: (interior if s < m else b) for s in range(1, m + 1)}


def gamma_factor_profile(k, m):
    """The same ledger after cancelling t^(I*m) in the Gamma moments."""
    a, b, _ = character(k)
    profile = {}
    for s in range(1, m):
        profile[s] = profile.get(s, 0) + a
    for s in range(1, m + 1):
        profile[s] = profile.get(s, 0) + b
    return profile


symbolic_profile_cells = 0
for k in range(2, 13):
    for m in range(17):
        require(width_factor_profile(k, m) == gamma_factor_profile(k, m), "formal Gamma factor profile failed")
        symbolic_profile_cells += 1

positive_t = (Fraction(1, 7), Fraction(2, 5), Fraction(1), Fraction(3, 2))
moment_cells = 0
ratio_cells = 0
# Materialized rational checks stop at k=5; the all-k test above works with
# exponent ledgers so k! exponents never create hundred-million-digit integers.
for k in range(2, 6):
    for t in positive_t:
        for m in range(17):
            require(width_flag(k, t, m) == gamma_product_moment(k, t, m), "Gamma moment identity failed")
            moment_cells += 1
        for m in range(1, 16):
            lhs = width_flag(k, t, m - 1) * width_flag(k, t, m + 1)
            rhs = width_flag(k, t, m) ** 2
            require(lhs > rhs, "strict adjacent log-convexity failed")
            require(lhs / rhs == ratio(k, t, m) / ratio(k, t, m - 1), "exact curvature ratio failed")
            ratio_cells += 1
print(jdump({"moment_identity":{"formal_profile_cells":symbolic_profile_cells,"rational_cells":moment_cells,"rational_k":"2..5","m":"0..16","positive_t":["1/7","2/5","1","3/2"]},"adjacent_curvature_cells":ratio_cells}))

# Exact multiplicative form of the alternating finite-difference hierarchy.
finite_difference_cells = 0
for k in range(2, 6):
    for t in positive_t:
        for m in range(11):
            for order in range(1, 7):
                numerator = Fraction(1)
                denominator = Fraction(1)
                for j in range(order + 1):
                    coefficient = (-1) ** (order - j) * comb(order, j)
                    value = ratio(k, t, m + j)
                    if coefficient > 0:
                        numerator *= value ** coefficient
                    elif coefficient < 0:
                        denominator *= value ** (-coefficient)
                quotient = numerator / denominator
                if order % 2:
                    require(quotient > 1, "odd finite log difference lost positivity")
                else:
                    require(quotient < 1, "even finite log difference lost negativity")
                finite_difference_cells += 1
print(jdump({"alternating_log_difference_cells":finite_difference_cells,"orders":"1..6"}))

# Strict total positivity: every tested generalized Hankel minor is positive.
hankel_cells = 0
hankel_digest_data = []
for k in range(2, 6):
    for t in (Fraction(1, 2), Fraction(1)):
        moments = [width_flag(k, t, m) for m in range(13)]
        for size in range(1, 4):
            index_sets = list(combinations(range(7), size))
            for rows in index_sets:
                for cols in index_sets:
                    det = determinant([[moments[i + j] for j in cols] for i in rows])
                    require(det > 0, "generalized Hankel minor is not strictly positive")
                    hankel_cells += 1
                    if len(hankel_digest_data) < 12:
                        hankel_digest_data.append([k, str(t), list(rows), list(cols), str(det)])
print(jdump({"strict_hankel_minors":{"cells":hankel_cells,"index_window":"0..6","k":"2..5","sizes":"1..3","t":["1/2","1"],"positive_controls":hankel_digest_data}}))

# Sharp sign/scope boundaries.
t_zero = [Fraction(1) for _ in range(5)]
zero_hankel = determinant([[t_zero[i + j] for j in range(2)] for i in range(2)])
negative_control = width_flag(2, Fraction(-2), 2)
require(zero_hankel == 0, "t=0 should be rank one")
require(negative_control == -1, "negative-t hostile changed")
print(jdump({"boundaries":{"moving_lower_full_corner":"uncontrolled extra low-resultant transport","negative_t_k2_F2":str(negative_control),"physical_width":False,"raw_chart":False,"t0_hankel2":str(zero_hankel)}}))

print("all_exact_checks=PASS")
