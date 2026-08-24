#!/usr/bin/env python3
"""Exact companion for THM-3945's non-simple J-sextic line gate.

The cited input is Degtyarev's four-cuspidal trigonal model and its
stereographic construction of the J_{2,0}+4A2 and J_{2,3}+3A2 sextics.
This script freezes the repo-derived part: the normalization, cone embedding,
both center charts, and the coefficient contradictions excluding a sixth
power among all projected line pullbacks.
"""

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
    gate(sp.cancel(expression) == 0, message)


t, x, y = sp.symbols("t x y")
T, S = sp.symbols("T S")
a, b = sp.symbols("a b")
A, B, C = sp.symbols("A B C")


# ---------------------------------------------------------------------------
# The four-cuspidal trigonal model and its normalization.
# ---------------------------------------------------------------------------

f = 4 * y**3 - (24 * x**3 + 3) * y + (8 * x**6 + 20 * x**3 - 1)
d = t**3 + 2
xt = 3 * t / d
yt = -(t**6 - 20 * t**3 - 8) / (2 * d**2)
zero(f.subs({x: xt, y: yt}), "trigonal parametrization")

# A rational inverse on a dense open set proves that the parametrization is
# birational, rather than merely a multiple cover.
inverse_numerator = 4 * x * y - 6 * x
inverse_denominator = 4 * x**3 - 2 * y - 1
zero((inverse_denominator * t - inverse_numerator).subs({x: xt, y: yt}),
     "trigonal rational inverse")

# The contraction Sigma_2 -> Q embeds the affine O(2)-chart as
# [1:x:x^2:y] in the quadric cone Z1^2=Z0*Z2.  On the normalization the
# four coordinates become binary sextics.
D = T**3 + 2 * S**3
Z0 = D**2
Z1 = 3 * T * S**2 * D
Z2 = 9 * T**2 * S**4
Z3 = -(T**6 - 20 * T**3 * S**3 - 8 * S**6) / 2
zero(Z1**2 - Z0 * Z2, "quadric-cone equation")
gate(all(sp.Poly(z, T, S).total_degree() == 6 for z in (Z0, Z1, Z2, Z3)),
     "cone coordinates are binary sextics")
gate(sp.gcd(sp.gcd(sp.gcd(Z0, Z1), Z2), Z3) == 1,
     "cone normalization map is basepoint free")
zero(Z1.subs({T: t, S: 1}) / Z0.subs({T: t, S: 1}) - xt,
     "cone x-coordinate")
zero(Z3.subs({T: t, S: 1}) / Z0.subs({T: t, S: 1}) - yt,
     "cone y-coordinate")


def coefficient(poly: sp.Expr, degree: int) -> sp.Expr:
    """Coefficient of t^degree after dehomogenizing S=1,T=t."""
    return sp.Poly(sp.expand(poly.subs({T: t, S: 1})), t).coeff_monomial(t**degree)


# ---------------------------------------------------------------------------
# Finite-base centers o=[1:a:a^2:b].
# ---------------------------------------------------------------------------

Uf = Z1 - a * Z0
Vf = Z2 - a**2 * Z0
Wf = Z3 - b * Z0
Ff = sp.expand(A * Uf + B * Vf + C * Wf)

finite_expected = {
    5: 0,
    4: 3 * A,
    3: -4 * a * A - 4 * a**2 * B + (10 - 4 * b) * C,
    2: 9 * B,
    1: 6 * A,
    0: -4 * a * A - 4 * a**2 * B + (4 - 4 * b) * C,
}
for degree, expected in finite_expected.items():
    zero(coefficient(Ff, degree) - expected,
         f"finite-center t^{degree} coefficient")

# A sixth power (T-rS)^6 with finite root has t^5 coefficient -6r.
# Since every Ff has t^5 coefficient zero, r=0.  Then t^4 and t^2 force
# A=B=0, while t^3 and t^0 would force b=5/2 and b=1 respectively.
gate(sp.solve([10 - 4 * b, 4 - 4 * b], [b], dict=True) == [],
     "finite center cannot support a sixth power at a finite point")

# A sixth power at t=infinity is constant after dehomogenization.  Again
# t^4,t^2 force A=B=0, but Wf cannot have both t^6 and t^3 zero.
gate(sp.solve([-sp.Rational(1, 2) - b, 10 - 4 * b], [b], dict=True) == [],
     "finite center cannot support a sixth power at infinity")


# ---------------------------------------------------------------------------
# Centers over the base point infinity: o=[0:0:1:b].
# ---------------------------------------------------------------------------

Ui = Z0
Vi = Z1
Wi = Z3 - b * Z2
Fi = sp.expand(A * Ui + B * Vi + C * Wi)

infinity_expected = {
    6: A - sp.Rational(1, 2) * C,
    5: 0,
    4: 3 * B,
    3: 4 * A + 10 * C,
    2: -9 * b * C,
    1: 6 * B,
    0: 4 * A + 4 * C,
}
for degree, expected in infinity_expected.items():
    zero(coefficient(Fi, degree) - expected,
         f"infinite-base-center t^{degree} coefficient")

# For a finite sixth-power root, r=0 as above.  The t^4,t^3,t^0 rows
# already have full rank in A,B,C, independently of b.
finite_root_minor = sp.Matrix([
    [0, 3, 0],
    [4, 0, 10],
    [4, 0, 4],
]).det()
gate(finite_root_minor == 72,
     "infinite-base center finite-root coefficient rank")

# For a sixth-power root at infinity, the t^6,t^4,t^3 rows have full rank.
infinite_root_minor = sp.Matrix([
    [1, 0, -sp.Rational(1, 2)],
    [0, 3, 0],
    [4, 0, 10],
]).det()
gate(infinite_root_minor == 36,
     "infinite-base center infinity-root coefficient rank")


# Every nonvertex point of Q is in exactly one of the two center charts:
# Z0!=0 scales to [1:a:a^2:b]; if Z0=0, the cone equation gives Z1=0 and
# nonvertex means Z2!=0, scaling to [0:0:1:b].  The vertex is the contracted
# exceptional section and is excluded in the stereographic construction.
center_charts = ["[1:a:a^2:b]", "[0:0:1:b]"]

summary = {
    "checks": CHECKS,
    "normalization": "P1",
    "cone": "Z1^2=Z0*Z2",
    "center_charts": center_charts,
    "finite_center_sixth_power": False,
    "infinite_base_center_sixth_power": False,
    "one_place_line": False,
    "scope": ["J_{2,0}+4A2", "J_{2,3}+3A2"],
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3945 non-simple weight-eight J-sextic one-place companion")
print(f"CHECKS={CHECKS}")
print("TRIGONAL_NORMALIZATION=P1;CONE=Z1^2-Z0*Z2")
print("CENTER_CHARTS=FINITE_BASE,INFINITE_BASE")
print("FINITE_CENTER_SIXTH_POWER=NONE")
print("INFINITE_BASE_CENTER_SIXTH_POWER=NONE")
print("J20_PLUS_4A2_ONE_PLACE_LINE=NO")
print("J23_PLUS_3A2_ONE_PLACE_LINE=NO")
print(f"SEMANTIC_SHA256={semantic}")
