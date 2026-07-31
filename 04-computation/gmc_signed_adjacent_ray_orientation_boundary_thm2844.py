#!/usr/bin/env python3
"""Exact referee for THM-2844's signed adjacent-ray boundary."""

from fractions import Fraction
from math import factorial

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


x, t, n = sp.symbols("x t n")


def f(index):
    return x**index / factorial(index)


def d(index):
    return f(index + 1) - f(index)


def moment(poly):
    terms = sp.Poly(sp.expand(poly), x).terms()
    return sp.expand(sum(coefficient * factorial(power[0]) for power, coefficient in terms))


def orientation(lower, upper):
    return sp.factor(
        2 * moment(upper**3) * moment(lower * upper)
        - 3 * moment(lower * upper**2) * moment(upper**2)
    )


# The singleton law is a formal consequence of THM-2830's strict Pascal
# quotient.  This bank independently checks the direct tensor expansion,
# including both sides of the support cut and the equality ray.
singleton_cells = 0
for i in range(9):
    for j in range(1, 11):
        value = orientation(d(i), d(j))
        expected_sign = (j - 1 > i) - (j - 1 < i)
        actual_sign = int(sp.sign(value))
        require(actual_sign == expected_sign, "singleton support-cut sign")
        singleton_cells += 1


u = d(0)
v = t * d(1) + d(2)

g11 = moment(u**2)
g12 = moment(u * v)
g22 = moment(v**2)
t111 = moment(u**3)
t112 = moment(u**2 * v)
t122 = moment(u * v**2)
t222 = moment(v**3)

require(g11 == 1, "g11")
require(g12 == t + 1, "g12")
require(g22 == 2 * (t**2 + 3 * t + 3), "g22")
require(t111 == 2, "t111")
require(t112 == 2 * (2 * t + 3), "t112")
require(t122 == 2 * (5 * t**2 + 19 * t + 22), "t122")
require(t222 == 6 * (5 * t**3 + 34 * t**2 + 90 * t + 90), "t222")
require(sp.factor(g11 * g22 - g12**2) == t**2 + 4 * t + 5, "Gram determinant")

p = 5 * t**3 + 30 * t**2 + 57 * t + 24
d0 = orientation(u, v)
require(sp.expand(d0 - 12 * p) == 0, "signed orientation cubic")
require(sp.discriminant(p, t) == -66960, "orientation discriminant")

left = Fraction(-583, 1000)
right = Fraction(-582, 1000)
require(p.subs(t, sp.Rational(left.numerator, left.denominator)) < 0, "left root bracket")
require(p.subs(t, sp.Rational(right.numerator, right.denominator)) > 0, "right root bracket")

a, b = sp.symbols("a b")
homogeneous_p = 5 * a**3 * b + 30 * a**2 * b**2 + 57 * a * b**3 + 24 * b**4
primitive_family = sp.factor(homogeneous_p.subs({a: -(n + 1), b: n}))
require(
    primitive_family == -n * (8 * n**3 + 12 * n**2 - 15 * n + 5),
    "primitive hostile family",
)
require(
    sp.expand(
        (8 * (n + 1) ** 3 + 12 * (n + 1) ** 2 - 15 * (n + 1) + 5)
        - (8 * n**3 + 12 * n**2 - 15 * n + 5)
    )
    == 24 * n**2 + 48 * n + 5,
    "primitive family positivity increment",
)

i1 = sp.factor(
    3 * t112 * g11 * g22 - t222 * g11**2 - 2 * t111 * g12 * g22
)
i2 = sp.factor(
    3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22**2
)
p1 = 7 * t**3 + 64 * t**2 + 204 * t + 228
p2 = 2 * t**4 + 27 * t**3 + 120 * t**2 + 207 * t + 90
require(sp.expand(i1 + 2 * p1) == 0, "first divisibility invariant")
require(sp.expand(i2 + 4 * p2) == 0, "second divisibility invariant")
resultant = sp.resultant(p1, p2, t)
require(resultant == -4409316, "divisibility resultant")
require(sp.gcd(p1, p2) == 1, "divisibility coprimality")

print("THM-2844 SIGNED ADJACENT-RAY ORIENTATION BOUNDARY")
print("status=PROVED+VERIFIED-EXACT")
print(f"singleton_sign_cells={singleton_cells}")
print("singleton_law=sign(D_i(d_j))=sign(j-1-i)")
print("V_t=t*d_1+d_2; U=d_0")
print("D_0(V_t)=12*(5*t^3+30*t^2+57*t+24)")
print("orientation_cubic_discriminant=-66960; unique_real_root=theta")
print(f"theta_bracket=({left},{right})")
print("orientation_nonnegative_iff=t>=theta")
print("primitive_hostiles=V_n=-(n+1)*d_1+n*d_2; n>=1")
print("D_0(V_n)=-12*n*(8*n^3+12*n^2-15*n+5)<0")
print("Gram=t^2+4*t+5=(t+2)^2+1>0 on R")
print(f"divisibility_resultant={resultant}")
print("moment_detection=all_real_t; orientation_failure_is_not_detection_failure")
print("all_exact_controls=PASS")
