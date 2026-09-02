#!/usr/bin/env python3
"""Hostile-referee exact probe for THM-4339's cubic extinction chart.

This does not duplicate the inherited global hull enumeration; it targets the
local strata easiest to misclassify: the persistent double node, the general
(not merely linear) simple-zero triple face, and the differential orders on
all balanced carriers.
"""

from sympy import Poly, Rational, diff, expand, factor, limit, simplify, symbols
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def require(value, label):
    if value is not True:
        raise AssertionError(label)


s, d, b, x, P, X, Y = symbols("s d b x P X Y")
a, c, W = symbols("a c W", nonzero=True)
B0, B1, B2 = symbols("B0 B1 B2", nonzero=True)
C0, C1 = symbols("C0 C1")


def local_equation(A, B, C):
    """The exact reciprocal equation with P=a+x and delta=d."""
    p = a + x
    return expand(
        (1 - d**2 * p * b**2)
        * (b**2 - p**2 * A - d * b * B - d**2 * b**2 * C)
        - d**2 * b**2 / 2
    )


# Double root, B(a)=0.  The quadratic tangent cone is a node uniformly over
# the base because its discriminant has a nonzero constant term.
A_double = W * x**2 * (a + x - c)
B_double_zero = B1 * x + B2 * x**2
F_double_zero = local_equation(A_double, B_double_zero, C0 + C1 * x)

def homogeneous_piece(poly, degree):
    out = 0
    for term in Poly(expand(poly), x, b).terms():
        (ix, ib), coefficient = term
        if ix + ib == degree:
            out += coefficient * x**ix * b**ib
    return expand(out)


double_tangent = homogeneous_piece(F_double_zero, 2)
q0 = 1 - d**2 * (C0 + Rational(1, 2))
expected_double_tangent = q0 * b**2 - d * B1 * x * b - a**2 * W * (a - c) * x**2
require(simplify(double_tangent - expected_double_tangent) == 0,
        "persistent double tangent cone")
double_tangent_disc = factor(Poly(double_tangent, b).discriminant())
require(simplify(double_tangent_disc / x**2
                 - (d**2 * B1**2 + 4 * q0 * a**2 * W * (a - c))) == 0,
        "persistent double node unit")


# Triple root and a GENERAL simple zero B=x(B1+B2*x), not the restricted
# linear-J specialization.  B2 is later and the weighted face is unchanged.
A_triple = W * x**3
B_triple_simple = B1 * x + B2 * x**2
F_triple_simple = local_equation(A_triple, B_triple_simple, C0 + C1 * x)
weighted_triple_simple = expand(F_triple_simple.subs({x: s**12 * X,
                                                       b: s**18 * Y,
                                                       d: s**6}))
triple_simple_face = factor(limit(weighted_triple_simple / s**36, s, 0))
expected_triple_simple_face = Y**2 - B1 * X * Y - a**2 * W * X**3
require(simplify(triple_simple_face - expected_triple_simple_face) == 0,
        "general simple-zero triple face")
require(all(simplify(diff(triple_simple_face, q).subs({X: 0, Y: 0})) == 0
            for q in (X, Y)), "simple-zero cubic singular")
require(simplify(
    (Y - B1 * X / 2)**2
    - X**2 * (B1**2 / 4 + a**2 * W * X)
    - triple_simple_face
) == 0, "simple-zero cubic rational normalization")


# Weighted denominator orders are D-r provided d(face)/dX is nonzero in the
# function field.  Check it is not the zero polynomial in every carrier.
double_face = Y**2 - B0 * Y - a**2 * W * (a - c) * X**2
triple_face = Y**2 - B0 * Y - a**2 * W * X**3
for face, label in ((double_face, "double"),
                    (triple_face, "triple"),
                    (triple_simple_face, "triple-simple-zero")):
    require(diff(face, X) != 0, label + " generic denominator")

orders = {
    "double_conic": 16 + 3 * 6 - (12 - 6),
    "triple_elliptic": 16 + 3 * 6 - (12 - 4),
    "triple_simple_zero_cubic": 16 + 3 * 18 - (36 - 12),
}
require(orders == {"double_conic": 28,
                   "triple_elliptic": 26,
                   "triple_simple_zero_cubic": 46},
        "differential exponent ledger")

print("THM4339_HOSTILE_REFEREE_LOCAL_PROBE=PASS")
print("double_tangent=" + str(double_tangent))
print("double_tangent_discriminant=" + str(double_tangent_disc))
print("general_triple_simple_zero_face=" + str(triple_simple_face))
print("orders=" + str(orders))
