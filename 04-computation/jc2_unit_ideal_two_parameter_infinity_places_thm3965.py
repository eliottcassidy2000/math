#!/usr/bin/env python3
"""Exact companion for THM-3965's two-parameter infinity-place gate."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.expand(expression) == 0, message)


A, C, U, V, t, g, h = sp.symbols("A C U V t g h")


# ---------------------------------------------------------------------------
# The unit-ideal binary cubic and its finite-root incidence.
# ---------------------------------------------------------------------------

Phi = A * U**3 + (C + g * A) * U**2 * V + (A * C - 1) * U * V**2 + h * A * V**3
psi = sp.expand(Phi.subs({U: t, V: 1}))
psi_t = sp.diff(psi, t)

zero(A * C - (A * C - 1) - 1, "coefficient ideal is one")
zero(Phi.subs(A, 0) - U * V * (C * U - V),
     "A=0 specialization is the three-factor scalar-unit hostile")
gate(sp.factor(Phi.subs(h, 0)).has(U),
     "h=0 binary cubic has the root U=0")

# Over k(A), linearity in C reduces irreducibility to a gcd.  The two roots
# of the C coefficient t(t+A) give the exact nonzero constant rows.
C_coefficient = sp.diff(psi, C)
C_free = sp.expand(psi.subs(C, 0))
zero(C_coefficient - t * (t + A), "linear-C coefficient")
zero(C_free.subs(t, 0) - h * A, "constant row at t=0")
zero(C_free.subs(t, -A) - A * (-A**3 + g * A**2 + 1 + h),
     "constant row at t=-A")

alpha = 2 * t**3 + g * t**2 - h
beta = t**4 - 2 * h * t
G = sp.expand(alpha * A**2 + beta * A + t**2)
zero(sp.resultant(psi, psi_t, C) - G,
     "linear-color resultant is the quadratic incidence")

color_equation = sp.expand((A + 2 * t) * C + 2 * A * g * t + 3 * A * t**2 - 1)
zero(color_equation - psi_t, "generic color recovery is the derivative equation")

# The only simultaneous zero of the two C coefficients is the raw origin.
zero(C_coefficient.subs(A, -2 * t) + t**2,
     "simultaneous color-coefficient row")
gate(psi_t.subs({A: 0, t: 0}) == -1,
     "raw incidence origin is extraneous, not a repeated root")


# ---------------------------------------------------------------------------
# Irreducibility and the unique reducible seam.
# ---------------------------------------------------------------------------

D = sp.expand(t**6 - 4 * (h + 2) * t**3 - 4 * g * t**2 + 4 * h * (h + 1))
zero(sp.discriminant(G, A) - t**2 * D,
     "quadratic-incidence discriminant")
gate(alpha.subs(t, 0) == -h,
     "h nonzero makes the quadratic incidence primitive")

# A rational square of the monic polynomial D is a polynomial square.  The
# monic cubic square comparison has a unique solution.
a2, b1, c0 = sp.symbols("a2 b1 c0")
Q = t**3 + a2 * t**2 + b1 * t + c0
square_rows = sp.Poly(sp.expand(D - Q**2), t).all_coeffs()
square_gb = sp.groebner(square_rows, a2, b1, c0, g, h, order="lex")
gate([sp.factor(p.as_expr()) for p in square_gb.polys] ==
     [a2, b1, 3 * c0 + 4, g, 3 * h + 4],
     "unique square-discriminant seam")
zero(D.subs({g: 0, h: -sp.Rational(4, 3)}) -
     (t**3 - sp.Rational(4, 3))**2,
     "exceptional D is a square")

G_exceptional = sp.factor(G.subs({g: 0, h: -sp.Rational(4, 3)}))
gate(G_exceptional == (2 * A + t) * (3 * A * t**3 + 2 * A + 3 * t) / 3,
     "exceptional incidence has exactly two factors")


# ---------------------------------------------------------------------------
# The projective A-infinity packet.
# ---------------------------------------------------------------------------

alpha_disc = sp.factor(sp.discriminant(alpha, t))
gate(alpha_disc == 4 * h * (g**3 - 27 * h),
     "alpha repeated-root discriminant")
zero(alpha.subs(h, g**3 / 27) -
     2 * (t + g / 3)**2 * (t - g / 6),
     "repeated-alpha boundary still has two distinct roots")

# With z=1/A, the recovered color is regular at every A-infinity point.
z = sp.symbols("z")
C_at_infinity_chart = -(2 * g * t + 3 * t**2 - z) / (1 + 2 * t * z)
zero(sp.limit(C_at_infinity_chart, z, 0) + 2 * g * t + 3 * t**2,
     "color remains finite at A infinity")

Delta = sp.expand(sp.discriminant(psi, t))
Delta_poly = sp.Poly(Delta, A, C)
gate(Delta_poly.total_degree() == 7,
     "target discriminant has fixed degree seven")
Delta_top = sum(coeff * A**monom[0] * C**monom[1]
                for monom, coeff in Delta_poly.terms()
                if sum(monom) == 7)
zero(Delta_top + 4 * A**4 * C**3,
     "fixed projective infinity support")
gate(Delta_poly.coeff_monomial(A) == 4,
     "linear A term makes the irreducible discriminant reduced")


# ---------------------------------------------------------------------------
# Unique reducible seam: both components still have multiple places.
# ---------------------------------------------------------------------------

Delta_exceptional = sp.factor(Delta.subs({g: 0, h: -sp.Rational(4, 3)}))
E1 = 12 * A**3 - 3 * A * C - 1
E2 = A * C**3 + 12 * A + 3 * C**2
zero(Delta_exceptional + E1 * E2 / 3,
     "exceptional target discriminant splits into two reduced factors")

# First factor: parameter A on G_m, with poles at A=0 and A=infinity.
C1 = 4 * A**2 - 1 / (3 * A)
zero(E1.subs(C, C1), "first exceptional component parametrization")

# Second factor: parameter t; t=0 and the three roots of 3t^3+2 are poles.
A2 = -3 * t / (3 * t**3 + 2)
C2 = 2 / t
zero(sp.together(E2.subs({A: A2, C: C2})).as_numer_denom()[0],
     "second exceptional component parametrization")
gate(sp.degree(3 * t**3 + 2, t) == 3 and
     sp.discriminant(3 * t**3 + 2, t) != 0,
     "second component has three distinct A-pole points besides t=0")


# ---------------------------------------------------------------------------
# Hostile genus-zero and repeated-alpha degenerations.
# ---------------------------------------------------------------------------

D_rational = sp.factor(D.subs({g: -sp.Rational(3, 4), h: -1}))
gate(D_rational == t**2 * (t - 1)**2 * (t**2 + 2 * t + 3),
     "rational normalization hostile has a residual conic")
alpha_rational = sp.factor(alpha.subs({g: -sp.Rational(3, 4), h: -1}))
gate(sp.discriminant(alpha_rational, t) == -sp.Rational(1701, 16),
     "rational hostile still has three distinct A-infinity points")

alpha_double = sp.factor(alpha.subs({g: 3, h: 1}))
gate(alpha_double == (t + 1)**2 * (2 * t - 1),
     "double-alpha hostile realizes the sharp two-support boundary")
gate(sp.factor(D.subs({g: 3, h: 1})) ==
     t**6 - 12 * t**3 - 12 * t**2 + 8,
     "double-alpha hostile remains outside the reducible seam")


summary = {
    "checks": CHECKS,
    "family": "A U3+(C+gA)U2V+(AC-1)UV2+hAV3",
    "domain": "h=0 reducible; h nonzero irreducible unit-ideal nonmonogenic",
    "incidence": "G=alpha A2+beta A+t2; alpha=2t3+gt2-h",
    "seam": "G reducible iff (g,h)=(0,-4/3)",
    "places": "alpha has at least two distinct roots; each gives A-infinity normalization place",
    "exception": "two discriminant components have respectively two and four infinity places",
    "hostile": "rational residual conic and repeated-alpha boundary still fail one-place",
    "scope": "constant two-parameter THM3907 deformation; bivariate coefficients and JC2 open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3965 unit-ideal two-parameter infinity-place companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=A_U3_PLUS_C_PLUS_GA_U2V_PLUS_AC_MINUS1_UV2_PLUS_HA_V3")
print("DOMAIN=H_ZERO_REDUCIBLE;H_NONZERO_IRREDUCIBLE_UNIT_IDEAL_NONMONOGENIC")
print("INCIDENCE=G_ALPHA_A2_PLUS_BETA_A_PLUS_T2;ALPHA_2T3_PLUS_GT2_MINUS_H")
print("REDUCIBLE_SEAM=G_ZERO_H_MINUS4_OVER3;TWO_TARGET_COMPONENTS")
print("INFINITY=ALPHA_HAS_AT_LEAST_TWO_DISTINCT_ROOTS;EACH_GIVES_A_PLACE")
print("EXCEPTIONAL_PLACES=FIRST_COMPONENT_2;SECOND_COMPONENT_4")
print("HOSTILES=RATIONAL_CONIC_STILL_3_ALPHA_PLACES;DOUBLE_ALPHA_HAS_2_DISTINCT_ALPHA_ROOTS")
print("CONCLUSION=NO_ONE_PLACE_DISCRIMINANT_MEMBER_IN_THIS_TWO_PARAMETER_FAMILY")
print("SCOPE=CONSTANT_G_H_DEFORMATION_ONLY;BIVARIATE_COEFFICIENTS_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
