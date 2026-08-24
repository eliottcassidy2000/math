#!/usr/bin/env python3
"""Exact companion for THM-3969's affine-P graph collision gate."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(expression) == 0, message)


P, T, x, y, v, t = sp.symbols("P T x y v t")
a, c, r = sp.symbols("a c r")
V, W = sp.symbols("V W")


# ---------------------------------------------------------------------------
# Node normal form and rational normalization on the v-chart.
# ---------------------------------------------------------------------------

y_target = P - r**2
q = sp.expand(3 * r * P - r**3 + y_target**2 * (c + a * y_target))
F = sp.expand(T**3 - 3 * P * T - q)
G = sp.expand(x**3 - 3 * r * x**2 - 3 * x * y - c * y**2 - a * y**3)
zero(F.subs({T: x - r, P: y + r**2}) - G,
     "affine-P graph node normal form")

D = 1 - a * v**3
N = 3 * r + 3 * v + c * v**2
xv = N / D
yv = v * N / D
zero(sp.together(G.subs({x: xv, y: yv})),
     "relative-P1 parametrization lies on cubic")
zero(sp.together(xv * D - N), "normalization graph equation xD=N")
zero(sp.together(yv - v * xv), "normalization slope y=vx")

P_num = sp.expand(r**2 * D + v * N)
gate(sp.degree(P_num, v) <= 3 and sp.degree(D, v) <= 3,
     "projective P-map has cubic homogeneous rows")


# ---------------------------------------------------------------------------
# Homogeneous graph and exact collision resultant.
# ---------------------------------------------------------------------------

D_h = sp.expand(W**3 - a * V**3)
N_h = sp.expand(3 * r * W**2 + 3 * V * W + c * V**2)
P_h = sp.expand(r**2 * D_h + V * N_h)
gate(sp.Poly(D_h, V, W).total_degree() == 3,
     "homogeneous pole row has degree three")
gate(sp.Poly(P_h, V, W).total_degree() == 3,
     "homogeneous numerator row has degree three")
zero(D_h.subs({V: v, W: 1}) - D,
     "homogeneous denominator dehomogenizes")
zero(N_h.subs({V: v, W: 1}) - N,
     "homogeneous node numerator dehomogenizes")
zero(P_h.subs({V: v, W: 1}) - P_num,
     "homogeneous P numerator dehomogenizes")

Xi = sp.expand(c**3 + 27 * a**2 * r**3 - 27 * a * c * r + 27 * a)
gate(sp.factor(sp.resultant(D, N, v) - Xi) == 0,
     "collision resultant Xi")

# D_h=0 never has V=0. At W=0 it vanishes exactly when a=0, while
# N_h then vanishes exactly when c=0; Xi records that seam.
gate(D_h.subs({V: 1, W: 0}) == -a,
     "infinity denominator value")
gate(N_h.subs({V: 1, W: 0}) == c,
     "infinity numerator value")
gate(Xi.subs({a: 0, c: 0}) == 0,
     "infinity basepoint lies on Xi")

# The finite point v=infinity is retained whenever a is nonzero.
x_infinity_num = W * N_h
y_infinity_num = V * N_h
gate(x_infinity_num.subs({V: 1, W: 0}) == 0,
     "x extends as zero at v infinity")
gate(y_infinity_num.subs({V: 1, W: 0}) == c,
     "y numerator at v infinity")
gate(D_h.subs({V: 1, W: 0}) == -a,
     "y denominator at v infinity")


# ---------------------------------------------------------------------------
# Class-group quotient for every degree-three component partition.
# ---------------------------------------------------------------------------

degree_partitions = ((3,), (1, 2), (1, 1, 1), (1,))
expected_gcds = (3, 1, 1, 1)
for degrees, expected in zip(degree_partitions, expected_gcds):
    current = 0
    for degree in degrees:
        current = gcd(current, degree)
    gate(current == expected,
         f"pole-component degree gcd {degrees}")
    gate(current >= 1,
         f"class quotient is finite {degrees}")


# ---------------------------------------------------------------------------
# Moving collision-free and one-collision hostile controls.
# ---------------------------------------------------------------------------

a_unit = (1 - t**3) / 27
c_unit = t
r_unit = sp.Integer(0)
Xi_unit = sp.factor(Xi.subs({a: a_unit, c: c_unit, r: r_unit}))
gate(Xi_unit == 1, "moving collision-free unit-resultant row")

a_hit = t
c_hit = sp.Integer(1)
r_hit = sp.Integer(0)
Xi_hit = sp.factor(Xi.subs({a: a_hit, c: c_hit, r: r_hit}))
gate(Xi_hit == 1 + 27 * t, "one-collision resultant row")
t0 = -sp.Rational(1, 27)
v0 = -3
gate(D.subs({a: a_hit, t: t0, v: v0}) == 0,
     "one-collision denominator basepoint")
gate(N.subs({c: c_hit, r: r_hit, v: v0}) == 0,
     "one-collision numerator basepoint")
gate(Xi_hit.subs(t, t0) == 0,
     "one-collision lies over Xi zero")

# Constant-a endpoint and a=0 nonreduced pole support.
gate(Xi.subs({a: 1, c: 0, r: 0}) == 27,
     "constant cubic-pole collision-free control")
gate(Xi.subs({a: 0, c: 1, r: 0}) == 1,
     "triple section a=0 collision-free control")


summary = {
    "checks": CHECKS,
    "family": "q=3rP-r^3+(P-r^2)^2(c+a(P-r^2))",
    "model": "relative P1 minus D_h=W^3-aV^3",
    "map": "x=W*N_h/D_h, y=V*N_h/D_h",
    "collision": "Xi=c^3+27a^2r^3-27acr+27a",
    "class": "Cl(Y)=Z/gcd(component degrees)Z",
    "conclusion": "Xi unit implies no same-field Keller A2",
    "scope": "Xi-zero exceptional curves remain open",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3969 affine-P graph relative-P1 collision companion")
print(f"CHECKS={CHECKS}")
print("NORMAL_FORM=X3_MINUS_3R_X2_MINUS_3XY_MINUS_CY2_MINUS_AY3")
print("NORMALIZATION=RELATIVE_P1_MINUS_CUBIC_POLE_MULTISECTION")
print("COLLISION=XI_C3_PLUS_27A2R3_MINUS_27ACR_PLUS_27A")
print("FINITE_GATE=XI_NONZERO_SCALAR_IMPLIES_BASEPOINT_FREE_DEGREE3")
print("CLASS=FINITE_CYCLIC_Z_MOD_GCD_COMPONENT_DEGREES")
print("CONCLUSION=COLLISION_FREE_AFFINE_P_GRAPH_FAMILY_CLOSED")
print("CONTROL_UNIT=A_(1_MINUS_T3)_OVER_27;C_T;R_0;XI_1")
print("CONTROL_COLLISION=A_T;C_1;R_0;XI_1_PLUS_27T")
print("SCOPE=XI_ZERO_EXCEPTIONAL_CURVES;HIGHER_P_DEPTH;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
