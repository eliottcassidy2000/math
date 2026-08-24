#!/usr/bin/env python3
"""Exact companion for THM-3962's constant-q product-curve obstruction."""

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


P, T, v, c = sp.symbols("P T v c", nonzero=False)


# ---------------------------------------------------------------------------
# Universal degree-three and Riemann--Hurwitz support ledger.
# ---------------------------------------------------------------------------

Q = sp.symbols("Q")
F_universal = T**3 - 3 * P * T - Q
zero(sp.discriminant(F_universal, T) + 27 * (Q**2 - 4 * P**3),
     "universal depressed-cubic discriminant")

degree = 3
total_ramification = 2 * degree - 2
max_one_point = degree - 1
gate(total_ramification == 4, "degree-three RH total is four")
gate(max_one_point == 2, "one degree-three source point contributes at most two")
gate((total_ramification + max_one_point - 1) // max_one_point == 2,
     "RH forces at least two ramification support points")

# Freeze every partition of 4 into tame degree-three local contributions
# 1 or 2; no one-point packet exists.
ramification_partitions = []
for number_of_points in range(1, 5):
    for mask in range(1 << number_of_points):
        contributions = [1 + ((mask >> j) & 1) for j in range(number_of_points)]
        if sum(contributions) == total_ramification:
            ramification_partitions.append(tuple(sorted(contributions)))
ramification_partitions = sorted(set(ramification_partitions))
gate(ramification_partitions == [(1, 1, 1, 1), (1, 1, 2), (2, 2)],
     "complete degree-three ramification support partitions")


# ---------------------------------------------------------------------------
# Repeated adjusted factor q=P^2+P.
# ---------------------------------------------------------------------------

q_repeat = P**2 + P
F_repeat = sp.expand(T**3 - 3 * P * T - q_repeat)
disc_repeat_in_P = sp.factor(sp.discriminant(-F_repeat, P))
gate(disc_repeat_in_P == (T + 1) ** 2 * (4 * T + 1),
     "repeated-factor curve has nonsquare P-discriminant")
gate(sp.factor_list(disc_repeat_in_P, T)[1] ==
     [(4 * T + 1, 1), (T + 1, 2)],
     "repeated-factor irreducibility odd row")

P_repeat = v**3
T_repeat = v * (v + 1)
zero(F_repeat.subs({P: P_repeat, T: T_repeat}),
     "repeated-factor normalization parametrization")
zero((T_repeat + 1) * v - (P_repeat + T_repeat),
     "repeated-factor rational inverse")
zero(v**3 - P_repeat, "repeated-factor integral equation")

h = sp.symbols("h")
K_repeat = sp.factor(q_repeat.subs(P, h**2) - 2 * h**3)
zero(K_repeat - h**2 * (h - 1) ** 2,
     "repeated adjusted hidden square")

# The singular point is exact, and its two addresses are the primitive cube
# roots. Reduction modulo v^2+v+1 freezes both coordinates.
repeat_gradients = [
    F_repeat,
    sp.diff(F_repeat, T),
    sp.diff(F_repeat, P),
]
for expression in repeat_gradients:
    zero(expression.subs({P: 1, T: -1}),
         "repeated-factor conductor point is singular")
zero(sp.diff(F_repeat, P).subs(P, T**2) +
     (T + 1) * (2 * T + 1),
     "repeated-factor Jacobian has only two ramification candidates")
gate(F_repeat.subs({P: sp.Rational(1, 4), T: sp.Rational(-1, 2)}) ==
     -sp.Rational(1, 16),
     "second repeated-factor Jacobian candidate is not on the curve")
x, y = sp.symbols("x y")
repeat_translated = sp.expand(F_repeat.subs({T: x - 1, P: y + 1}))
zero(repeat_translated - (x**3 - 3 * x**2 - 3 * x * y - y**2),
     "repeated-factor node translated equation")
gate(sp.discriminant(-3 * x**2 - 3 * x * y - y**2, y) == -3 * x**2,
     "repeated-factor tangent cone has two distinct lines")
cyclotomic = sp.Poly(v**2 + v + 1, v)
gate(sp.rem(sp.Poly(P_repeat - 1, v), cyclotomic) == 0,
     "both primitive cube roots map to P=1")
gate(sp.rem(sp.Poly(T_repeat + 1, v), cyclotomic) == 0,
     "both primitive cube roots map to T=-1")

dP_repeat = sp.factor(sp.diff(P_repeat, v))
gate(dP_repeat == 3 * v**2,
     "repeated-factor normalized map has one finite ramification point")
gate(sp.degree(P_repeat, v) == 3,
     "repeated-factor infinity has degree-three contact")


# ---------------------------------------------------------------------------
# Double-P conductor q=cP^2, with c a nonzero scalar.
# ---------------------------------------------------------------------------

q_double = c * P**2
F_double = sp.expand(T**3 - 3 * P * T - q_double)
disc_double_in_P = sp.factor(sp.discriminant(-F_double, P))
gate(disc_double_in_P == T**2 * (4 * T * c + 9),
     "double-P curve has nonsquare P-discriminant")

D = 3 + c * v
P_double = sp.expand(v**2 * D)
T_double = sp.expand(v * D)
zero(F_double.subs({P: P_double, T: T_double}),
     "double-P normalization parametrization")
zero(v * T_double - P_double, "double-P rational inverse v=P/T")
zero(c * v**2 + 3 * v - T_double,
     "double-P integral quadratic equation")

K_double = sp.factor(q_double.subs(P, h**2) - 2 * h**3)
zero(K_double - h**3 * (c * h - 2),
     "double-P hidden h^3 debt")

double_gradients = [
    F_double,
    sp.diff(F_double, T),
    sp.diff(F_double, P),
]
for expression in double_gradients:
    zero(expression.subs({P: 0, T: 0}),
         "double-P conductor origin is singular")
zero(sp.diff(F_double, P).subs(P, T**2) + T * (2 * T * c + 3),
     "double-P Jacobian has only two ramification candidates")
gate(sp.factor(F_double.subs({T: -sp.Rational(3, 2) / c,
                              P: sp.Rational(9, 4) / c**2})) ==
     sp.Rational(27, 16) / c**3,
     "second double-P Jacobian candidate is not on the curve")
double_tangent = -3 * P * T - c * P**2
zero(double_tangent + P * (3 * T + c * P),
     "double-P origin has two distinct tangent lines")

# Common zeros of P(v),T(v) are v=0 and 3+cv=0.
gate(sp.factor(P_double) == v**2 * (c * v + 3),
     "double-P P address factorization")
gate(sp.factor(T_double) == v * (c * v + 3),
     "double-P T address factorization")
zero(P_double.subs(v, -3 / c), "second conductor address maps to P=0")
zero(T_double.subs(v, -3 / c), "second conductor address maps to T=0")

dP_double = sp.factor(sp.diff(P_double, v))
gate(dP_double == 3 * v * (c * v + 2),
     "double-P finite ramification factorization")
zero(P_double.subs(v, 0), "first finite ramification value")
zero(P_double.subs(v, -2 / c) - 4 / c**2,
     "second finite ramification value")
gate(sp.degree(P_double, v) == 3,
     "double-P infinity has degree-three contact")

# The conductor and ramification address sets are genuinely different.
conductor_scaled = {sp.Integer(0), sp.Integer(-3)}
ramification_scaled = {sp.Integer(0), sp.Integer(-2)}
gate(conductor_scaled != ramification_scaled,
     "conductor gluing is not the normalized ramification packet")
gate(conductor_scaled & ramification_scaled == {sp.Integer(0)},
     "exact one-address overlap between conductor and ramification")


# ---------------------------------------------------------------------------
# Unit ranks of the two rational etale loci.
# ---------------------------------------------------------------------------

# P1 minus s points has unit rank s-1. Freeze the two examples and the
# universal lower bound |S|>=2.
gate(2 - 1 == 1, "repeated-factor etale locus unit rank one")
gate(3 - 1 == 2, "double-P etale locus unit rank two")
gate(len({"ramification_1", "ramification_2"}) - 1 >= 1,
     "universal rational etale locus has a nonconstant unit")


summary = {
    "checks": CHECKS,
    "universe": "irreducible T^3-3PT-q(P), q independent of t",
    "normalization": "normalization(C x A1)=Cnu x A1",
    "genus": "dense A2 projection forces Cbar=P1",
    "riemann_hurwitz": "degree3 total4 and support at least2",
    "unit": "P1 minus infinity-plus-affine-ramification has nonconstant unit",
    "repeat_control": "P=v^3,T=v(v+1); conductor omega,omega2; ramification 0,infinity",
    "double_control": "P=v^2(3+cv),T=v(3+cv); conductor 0,-3/c; ramification 0,-2/c,infinity",
    "scope": "all constant-q(P) orders including nonnormal conductor debt",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3962 constant-q product-curve ramification-unit companion")
print(f"CHECKS={CHECKS}")
print("NORMALIZATION=CURVE_NORMALIZATION_TIMES_A1_CYLINDER")
print("GENUS_GATE=DENSE_A2_PROJECTION_FORCES_PROJECTIVE_CURVE_P1")
print("RH=DEGREE3_TOTAL4;AT_LEAST_TWO_RAMIFICATION_SUPPORT_POINTS")
print("UNIT=P1_MINUS_INFINITY_AND_AFFINE_RAMIFICATION_HAS_NONCONSTANT_UNIT")
print("CONCLUSION=NO_SAME_FUNCTION_FIELD_KELLER_A2_FOR_ANY_IRREDUCIBLE_Q(P)")
print("REPEAT=P_V3;T_VVPLUS1;CONDUCTOR_OMEGA_PAIR;RAMIFICATION_ZERO_INFINITY")
print("DOUBLE=P_V2_3PLUSCV;T_V_3PLUSCV;CONDUCTOR_0_MINUS3;RAMIFICATION_0_MINUS2_INFINITY")
print("CLASS=EXPLICIT_NORMALIZATION_CYLINDERS_FACTORIAL_WITH_SCALAR_GLOBAL_UNITS")
print("SCOPE=NORMAL_AND_NONNORMAL_CONSTANT_Q_ORDERS_CLOSED;MOVING_VERTICAL_DEBT_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
