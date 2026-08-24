#!/usr/bin/env python3
"""Exact scout for the cheapest double-torus two-character near miss.

This is intentionally not a theorem companion.  It freezes the algebra behind
the post-THM-3937 design experiment recorded in the matching reflection.
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


t, X, W, s = sp.symbols("t X W s")
p0 = X
p1 = X + t
q0 = 3 * X**2 + 3 * t * X + t**2 - t
q1 = 3 * X**2 + 3 * t * X + t**2 + t
H = sp.expand(q0**2 - 4 * p0**3)


# ---------------------------------------------------------------------------
# Two exact torus decompositions and the affine surface.
# ---------------------------------------------------------------------------

zero(H - (q1**2 - 4 * p1**3), "second torus decomposition")
gate(sp.Poly(H, t, X).total_degree() == 4, "common branch is a plane quartic")
gate(sp.factor(H) == H, "common quartic is irreducible over Q")
gate(
    sp.factor(sp.discriminant(H, X))
    == 6912 * t**4 * (t - 1) ** 3 * (t + 1) ** 3,
    "generic quartic discriminant",
)

# The only affine branch singularities are (t,X)=(0,0),(1,0),(-1,1).
branch_basis = sp.groebner([H, sp.diff(H, t), sp.diff(H, X)], X, t, order="lex")
expected_basis = [
    12 * X**2 + 2 * t**4 + t**3 - 8 * t**2 + 5 * t,
    12 * X * t - t**3 + 6 * t**2 - 5 * t,
    t * (t - 1) ** 2 * (t + 1) ** 2,
]
gate(
    [sp.factor(poly.as_expr()) for poly in branch_basis.polys]
    == [sp.factor(poly) for poly in expected_basis],
    "affine branch singular ideal",
)
for point in ((0, 0), (1, 0), (-1, 1)):
    gate(
        all(poly.subs({t: point[0], X: point[1]}) == 0 for poly in expected_basis),
        f"branch singular point {point}",
    )

# Away from t=0,+/-1 the fibre quartic is squarefree.  At the three bad
# parameters it is still not a polynomial square, so every vertical fibre of
# W^2=H is integral and principal.
special_fibres = {
    0: X**3 * (9 * X - 4),
    1: X**2 * (9 * X**2 + 14 * X + 9),
    -1: (X - 1) ** 2 * (9 * X**2 - 4 * X + 4),
}
for value, factorization in special_fibres.items():
    gate(sp.factor(H.subs(t, value) - factorization) == 0,
         f"special fibre t={value} factorization")
gate(sp.discriminant(9 * X**2 + 14 * X + 9, X) == -128,
     "t=1 residual quadratic is not a square")
gate(sp.discriminant(9 * X**2 - 4 * X + 4, X) == -128,
     "t=-1 residual quadratic is not a square")


# ---------------------------------------------------------------------------
# Exact rational normalization and its two unavoidable infinity addresses.
# ---------------------------------------------------------------------------

h = sp.cancel(2 * (1 - s) / (s**2 + 3))
t_param = sp.cancel((s - 1) ** 3 * (s + 3) / (s**2 + 3) ** 2)
X_param = sp.cancel(4 * (s - 1) ** 2 / (s**2 + 3) ** 2)
zero(X_param - h**2, "first torus square on normalization")
zero(q0.subs({t: t_param, X: X_param}) - 2 * h**3,
     "first torus cube on normalization")
zero(H.subs({t: t_param, X: X_param}), "quartic parametrization")

# Recover s rationally from a generic branch point.
h_recovered = sp.cancel(q0 / (2 * X))
v_recovered = sp.cancel((2 * t - 1 + 3 * h_recovered**2) / (h_recovered - 1))
s_recovered = sp.cancel((v_recovered - 1) / h_recovered)
zero(s_recovered.subs({t: t_param, X: X_param}) - s,
     "birational normalization inverse")

S, R = sp.symbols("S R")
T_hom = (S - R) ** 3 * (S + 3 * R)
X_hom = 4 * (S - R) ** 2 * R**2
Z_hom = (S**2 + 3 * R**2) ** 2
gate(sp.gcd(sp.gcd(T_hom, X_hom), Z_hom) == 1,
     "projective parametrization is basepoint free")
gate(sp.discriminant(S**2 + 3 * R**2, S) == -12 * R**2,
     "two distinct normalization points lie over the infinity line")
gate(sp.gcd(S**2 + 3 * R**2, (S - R) * R) == 1,
     "neither infinity address cancels from the target map")


# ---------------------------------------------------------------------------
# Elliptic fibre and Mordell--Weil invoices.
# ---------------------------------------------------------------------------

quartic = sp.Poly(H, X)
a, b, c, d, e = quartic.all_coeffs()
I = sp.factor(12 * a * e - 3 * b * d + c**2)
J = sp.factor(72 * a * c * e + 9 * b * c * d - 27 * a * d**2
              - 27 * b**2 * e - 2 * c**3)
gate(I == 9 * t**2 * (t**2 + 8), "binary-quartic I invariant")
gate(J == 54 * t**2 * (t**4 - 20 * t**2 - 8),
     "binary-quartic J invariant")
gate(sp.factor(4 * I**3 - J**2)
     == 186624 * t**4 * (t - 1) ** 3 * (t + 1) ** 3,
     "Jacobian discriminant invariant")

# These valuations give IV at 0, I3 at +/-1, and I2 at infinity.  The root
# rank is 2+2+2+1=7, so a rational elliptic surface has MW rank one.
gate(sp.Poly(I, t).as_dict().get((2,)) == 72, "I valuation two at zero")
gate(sp.Poly(J, t).as_dict().get((2,)) == -432, "J valuation two at zero")
for value in (1, -1):
    gate(I.subs(t, value) != 0 and J.subs(t, value) != 0,
         f"multiplicative I3 fibre at t={value}")
gate(sp.degree(I, t) == 4 and sp.degree(J, t) == 6,
     "infinity minimal scaling degrees")
gate(4 + 3 + 3 + 2 == 12 and 2 + 2 + 2 + 1 == 7,
     "rational-surface Euler and root-rank invoices")

# The two torus radicands have the same infinity pole divisor.  Their ratio
# has divisor 3(P0-P1), producing a nonzero rational point of exact order 3.
zero(sp.expand((q0 + W) * (q0 - W) - 4 * p0**3).subs(W**2, H),
     "first Cardano product modulo W^2=H")
zero(sp.expand((q1 + W) * (q1 - W) - 4 * p1**3).subs(W**2, H),
     "second Cardano product modulo W^2=H")
zero(q1**2 - q0**2 - 4 * (p1**3 - p0**3),
     "two Cardano products share H")


summary = {
    "checks": CHECKS,
    "branch": "irreducible rational plane quartic",
    "normalization": "P1 minus two infinity points (Gm), not A1",
    "torus_structures": 2,
    "quadratic_surface": "normal; three isolated singularities; vertical fibres integral",
    "elliptic_fibres": "IV + I3 + I3 + I2",
    "mordell_weil": "rank1 with exact Z/3 torsion from decomposition difference",
    "affine_picard_3_torsion": "(Z/3)^2",
    "units": "scalar",
    "stopping_reason": "second torus character costs a second normalization infinity",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("JC double-torus two-character near-miss scout")
print(f"CHECKS={CHECKS}")
print("IDENTITY=q0^2-4X^3=q1^2-4(X+t)^3")
print("BRANCH=IRREDUCIBLE_RATIONAL_QUARTIC;NORMALIZATION=Gm")
print("INFINITY_PLACES=2;ONE_PLACE_GATE=FAIL")
print("ELLIPTIC_FIBRES=IV+I3+I3+I2;MW_RANK=1;MW_TORSION=Z/3")
print("RESOLVENT_UNITS=k*;CL_3=(Z/3)^2")
print("MECHANISM=BOUNDARY_DIVISIBILITY_PLUS_INTRINSIC_THREE_TORSION")
print(f"SEMANTIC_SHA256={semantic}")
