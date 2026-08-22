#!/usr/bin/env python3
"""Exact controls for THM-3183's Hecke square and wedge continuant."""

from itertools import combinations
from math import comb, factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compound(matrix):
    """Second compound in the ordered pairs (01,02,12)."""

    pairs = ((0, 1), (0, 2), (1, 2))
    return sp.Matrix([
        [sp.det(matrix.extract(rows, columns)) for columns in pairs]
        for rows in pairs
    ])


def vp_rational(value, prime):
    value = sp.Rational(value)
    require(value != 0, "valuation of zero")
    numerator = abs(int(value.p))
    denominator = abs(int(value.q))
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def smith_valuations(matrix, prime):
    """Invariant-factor valuations from determinantal divisors."""

    one = min(vp_rational(entry, prime)
              for entry in matrix if entry != 0)
    two = min(
        vp_rational(sp.det(matrix.extract(rows, columns)), prime)
        for rows in combinations(range(3), 2)
        for columns in combinations(range(3), 2)
        if sp.det(matrix.extract(rows, columns)) != 0
    )
    three = vp_rational(sp.det(matrix), prime)
    return one, two - one, three - two


n, d, v = sp.symbols("n d v")
Delta = 1 - 4 * d * v


def weighted_transfer(index=n, offset=d, parameter=v):
    return sp.Matrix([
        [-(index + 1), 2 * parameter * (index + 1), offset],
        [-(index + 1) * (2 * index + 3 + 2 * offset),
         (index + 1) * (1 + 2 * parameter * (2 * index + 3)),
         offset * (2 * index + 3)],
        [0, 0, offset],
    ])


def scalar_transfer(index=n, offset=d, parameter=v):
    discriminant = 1 - 4 * offset * parameter
    return sp.Matrix([
        [2 * (index + 1) * (2 * index + 1) * parameter,
         index * (index + 1) * discriminant,
         offset - index - 1],
        [1, 0, 0],
        [0, 0, offset],
    ])


def gauge(index=n, offset=d, parameter=v):
    discriminant = 1 - 4 * offset * parameter
    return sp.Matrix([
        [1, 0, 0],
        [(1 + 2 * (2 * index + 1) * parameter) / (2 * parameter),
         index * discriminant / (2 * parameter),
         -1 / (2 * parameter)],
        [0, 0, 1],
    ])


G = weighted_transfer()
S = scalar_transfer()
P_in = gauge()
P_out = gauge(n + 1, d, v)

require((P_out * S - G * P_in).applyfunc(sp.factor) == sp.zeros(3),
        "Hecke square does not commute")
require(sp.factor(sp.det(P_in) - n * Delta / (2 * v)) == 0,
        "input gauge determinant")
require(sp.factor(sp.det(P_out) - (n + 1) * Delta / (2 * v)) == 0,
        "output gauge determinant")
require(sp.factor(sp.det(S) + d * n * (n + 1) * Delta) == 0,
        "scalar determinant")
require(sp.factor(sp.det(G) + d * (n + 1) ** 2 * Delta) == 0,
        "weighted determinant")

# Exact scalar exterior transfer and its two-step hidden-to-visible pivot.
a = 2 * (n + 1) * (2 * n + 1) * v
b = n * (n + 1) * Delta
c = d - n - 1
expected_wedge = sp.Matrix([
    [-b, -c, 0],
    [0, a * d, b * d],
    [0, d, 0],
])
W = compound(S)
require((W - expected_wedge).applyfunc(sp.factor) == sp.zeros(3),
        "scalar exterior matrix")
W_next = compound(scalar_transfer(n + 1, d, v))
source = sp.Matrix([0, 0, 1])
two_step = (W_next * W * source).applyfunc(sp.factor)
expected_pivot = sp.factor((n + 2 - d) * n * (n + 1) * Delta * d)
require(sp.factor(two_step[0] - expected_pivot) == 0,
        "two-step visible pivot")
require(sp.factor(two_step[1]
                  - b * d * (2 * (n + 2) * (2 * n + 3) * v) * d) == 0,
        "two-step transverse coordinate")
require(sp.factor(two_step[2] - b * d**2) == 0,
        "two-step returning coordinate")

# Finite exact determinantal-divisor controls of all four reset maps and
# their exterior squares.  The fixtures cover p=5,7,11,13 and distinct s,v.
fixtures = ((5, 1, 1), (7, 2, 2), (11, 3, 2), (13, 4, 3))
for prime, offset_small, parameter in fixtures:
    offset = prime + offset_small
    discriminant = 1 - 4 * offset * parameter
    require((2 * offset_small * parameter * discriminant) % prime != 0,
            "nonunit Smith fixture")
    prime_exact = sp.Integer(prime)
    offset_exact = sp.Integer(offset)
    parameter_exact = sp.Integer(parameter)
    maps = (
        gauge(prime_exact - 1, offset_exact, parameter_exact),
        scalar_transfer(prime_exact - 1, offset_exact, parameter_exact),
        gauge(prime_exact, offset_exact, parameter_exact),
        weighted_transfer(prime_exact - 1, offset_exact, parameter_exact),
    )
    expected = ((0, 0, 0), (0, 0, 1), (0, 0, 1), (0, 1, 1))
    expected_exterior = ((0, 0, 0), (0, 1, 1),
                         (0, 1, 1), (1, 1, 2))
    actual = tuple(smith_valuations(item, prime) for item in maps)
    actual_exterior = tuple(
        smith_valuations(compound(item), prime) for item in maps)
    require(actual == expected, ("reset Smith profile", prime, actual))
    require(actual_exterior == expected_exterior,
            ("exterior reset Smith profile", prime, actual_exterior))

# Sharp same-Smith projection hostile.
pi = sp.symbols("pi")
reset = sp.diag(1, pi, pi)
first_layer = compound(reset) / pi
identity = sp.eye(3)
unit_tail = sp.eye(3)
unit_tail[0, 2] = 1
deep_tail = sp.eye(3)
deep_tail[0, 2] = pi
hostile = tuple(sp.factor((compound(tail) * first_layer * source)[0])
                for tail in (identity, unit_tail, deep_tail))
require(hostile == (0, -pi, -pi**2), "same-Smith hostile")
require(all(sp.factor(sp.det(tail) - 1) == 0
            for tail in (identity, unit_tail, deep_tail)),
        "hostile tail is not unimodular")

# Offset-six coefficient-degree PRS pivots, normalized by (2p)!.
p = sp.symbols("p", integer=True, positive=True)
offset6 = p + 6


def factorial_ratio(argument):
    shift = sp.expand(argument - 2 * p)
    require(shift.is_Integer, "invalid top-factorial shift")
    shift = int(shift)
    if shift >= 0:
        return sp.prod(2 * p + index for index in range(1, shift + 1))
    return sp.cancel(1 / sp.prod(2 * p - index
                                for index in range(0, -shift)))


def top_coefficient(degree_shift, coefficient_shift):
    """[v^(p+coefficient_shift)] M_(p+degree_shift)/(2p)! ."""

    codimension = degree_shift - coefficient_shift
    require(codimension >= 0, "negative top codimension")
    choose = (sp.prod(p + degree_shift - index
                      for index in range(codimension))
              * sp.Rational(1, factorial(codimension)))
    total = sum(
        sp.binomial(codimension, ell)
        * offset6 ** (codimension - ell)
        * (-1) ** ell
        * factorial_ratio(2 * (p + coefficient_shift) + ell)
        for ell in range(codimension + 1)
    )
    return sp.factor(sp.cancel(choose * total))


A = {shift: top_coefficient(4, shift) for shift in range(-2, 5)}
B = {shift: top_coefficient(5, shift) for shift in range(-1, 6)}
C = (2 * p + 10) * (2 * p + 9)
R = {
    shift: sp.factor(
        (2 * p + 7) * B[shift]
        - (2 * p + 7) * C * A[shift - 1]
        + 2 * (p + 5) * (p + 6) * A[shift]
    )
    for shift in range(-1, 4)
}
H = 24 * p + 109
D2 = (p + 5) * (2 * p + 3) * H**2
N1 = -2 * (2 * p + 5) * (2 * p + 7) * (2 * p + 3) * H
N0 = 4 * (p + 6) * (2 * p + 5) * (28 * p + 123)
S_prs = {
    shift: sp.factor(D2 * A[shift]
                     - N1 * R[shift - 1]
                     - N0 * R[shift])
    for shift in range(3)
}
J = (256 * p**4 - 27648 * p**3 - 365600 * p**2
     - 1528800 * p - 2096649)
expected_R = (-8 * sp.prod(p + index for index in range(1, 6))
              * (2 * p + 1) * (2 * p + 3) * H)
expected_S = (4 * sp.prod(p + index for index in range(1, 6))
              * (2 * p + 7) * J)
require(sp.factor(R[3] - expected_R) == 0, "offset-six H pivot")
require(sp.factor(S_prs[2] - expected_S) == 0, "offset-six J pivot")


def direct_moment_coefficient(order, scalar_d, coefficient_index):
    return comb(order, coefficient_index) * sum(
        comb(order - coefficient_index, ell)
        * scalar_d ** (order - coefficient_index - ell)
        * (-1) ** ell
        * factorial(2 * coefficient_index + ell)
        for ell in range(order - coefficient_index + 1)
    )


DIRECT_TOP_JET_CHECKS = 0
for prime in (5, 7, 11):
    denominator = factorial(2 * prime)
    scalar_d = prime + 6
    direct_a = {
        shift: direct_moment_coefficient(prime + 4, scalar_d,
                                         prime + shift)
        for shift in A
    }
    direct_b = {
        shift: direct_moment_coefficient(prime + 5, scalar_d,
                                         prime + shift)
        for shift in B
    }
    direct_r = {
        shift: ((2 * prime + 7) * direct_b[shift]
                - (2 * prime + 7) * ((2 * prime + 10) * (2 * prime + 9))
                  * direct_a[shift - 1]
                + 2 * (prime + 5) * (prime + 6) * direct_a[shift])
        for shift in R
    }
    h_value = 24 * prime + 109
    d2_value = (prime + 5) * (2 * prime + 3) * h_value**2
    n1_value = (-2 * (2 * prime + 5) * (2 * prime + 7)
                * (2 * prime + 3) * h_value)
    n0_value = (4 * (prime + 6) * (2 * prime + 5)
                * (28 * prime + 123))
    direct_s = {
        shift: (d2_value * direct_a[shift]
                - n1_value * direct_r[shift - 1]
                - n0_value * direct_r[shift])
        for shift in S_prs
    }
    for symbolic, direct, label in ((A, direct_a, "A"),
                                    (B, direct_b, "B"),
                                    (R, direct_r, "R"),
                                    (S_prs, direct_s, "S")):
        for shift, expression in symbolic.items():
            require(expression.subs(p, prime)
                    == sp.Rational(direct[shift], denominator),
                    ("direct top jet", label, prime, shift))
            DIRECT_TOP_JET_CHECKS += 1

j = sp.symbols("j", integer=True)
bare_pivot = sp.factor(expected_pivot.subs({n: p + j, d: p + 6}))
expected_bare = sp.factor(
    (j - 4) * (p + j) * (p + j + 1) * (p + 6)
    * (1 - 4 * (p + 6) * v)
)
require(sp.factor(bare_pivot - expected_bare) == 0,
        "offset-six bare continuant pivot")

print("THM-3183 FACTORIAL HECKE LATTICE / WEDGE CONTINUANT EXACT CONTROL")
print("hecke_square=P_(n+1)S_n=G_nP_n")
print("reset_smith_profiles=P_in:111,S:11p,P_out:11p,G:1pp")
print("reset_length_balance=2=1+1-0")
print("exterior_profiles=P_in:111,S:1pp,P_out:1pp,G:ppp2")
print("exterior_length_balance=4=2+2-0")
print("smith_fixture_count=" + str(len(fixtures)))
print("two_step_pivot=(n+2-d)n(n+1)Delta*d")
print("same_smith_hostile_observed_valuations=(infinity,1,2)")
print("offset6_prs_leading_pivots=(H,J)")
print("offset6_direct_top_jet_checks=" + str(DIRECT_TOP_JET_CHECKS))
print("offset6_bare_wall=(j-4)(p+j)(p+j+1)(p+6)Delta")
print("scope=no_floor_s_over_2_or_arbitrary_offset_claim")
print("ALL EXACT CHECKS PASSED")
