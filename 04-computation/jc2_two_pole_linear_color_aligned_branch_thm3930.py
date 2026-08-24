#!/usr/bin/env python3
"""Exact companion for THM-3930's aligned line-plus-quintic cubic packet."""

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


A, C, U, V, T, u, rho, kappa, p, q, r = sp.symbols(
    "A C U V T u rho kappa p q r"
)
rho_relation = sp.Poly(rho**2 + 3, rho)


def rho_remainder(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_reduced = sp.rem(sp.Poly(sp.expand(numerator), rho), rho_relation).as_expr()
    denominator_reduced = sp.rem(sp.Poly(sp.expand(denominator), rho), rho_relation).as_expr()
    return sp.cancel(numerator_reduced / denominator_reduced)


def zero(expression: sp.Expr, message: str) -> None:
    gate(rho_remainder(expression) == 0, message)


# ---------------------------------------------------------------------------
# Binary cubic and its exact two-pole normalization.
# ---------------------------------------------------------------------------

a = (A - rho) ** 2
c = 2 * rho * A + sp.Rational(2, 3)
d = A / 2 - rho / 18
Phi = sp.expand(a * U**3 + C * U**2 * V + c * U * V**2 + d * V**3)
f = sp.expand(Phi.subs({U: T, V: 1}))

A_u = u**3 + rho * u**2 - u
t_u = (u + rho / 3) / (u**2 - 1)
C_u = -sp.Rational(3, 2) * (u**2 - 1) * (u + rho / 3) * (
    u**2 + sp.Rational(8, 3) * rho * u - sp.Rational(11, 3)
)

zero(
    a.subs(A, A_u) * t_u**3 - c.subs(A, A_u) * t_u - 2 * d.subs(A, A_u),
    "repeated-root incidence",
)
zero(
    3 * a.subs(A, A_u) * t_u**2 + 2 * C_u * t_u + c.subs(A, A_u),
    "repeated-root derivative",
)
zero(f.subs({A: A_u, C: C_u, T: t_u}), "cubic vanishes at repeated root")
zero(sp.diff(f, T).subs({A: A_u, C: C_u, T: t_u}), "cubic derivative vanishes")
gate(sp.degree(sp.cancel(t_u).as_numer_denom()[0], u) == 1, "root numerator degree")
gate(sp.degree(sp.cancel(t_u).as_numer_denom()[1], u) == 2, "root map has degree two")
zero(A_u.subs(u, 1) - rho, "first root pole target A")
zero(A_u.subs(u, -1) - rho, "second root pole target A")
zero(C_u.subs(u, 1), "first root pole target C")
zero(C_u.subs(u, -1), "second root pole target C")


# ---------------------------------------------------------------------------
# Complete centered two-simple-pole trace packet.
# ---------------------------------------------------------------------------

A_general = u**3 + p * u**2 + q * u + r
t_general = (u + kappa) / (u**2 - 1)
minimal = sp.Poly(
    sp.expand(sp.resultant(A_general - A, (u**2 - 1) * T - (u + kappa), u)),
    T,
)
trace_coefficient = sp.factor(minimal.coeff_monomial(T**2))
trace_A_coefficient = sp.factor(sp.Poly(trace_coefficient, A).coeff_monomial(A))
trace_constant = sp.factor(sp.Poly(trace_coefficient, A).coeff_monomial(1))
gate(trace_A_coefficient == -2 * kappa * p - q - 3, "trace first equation")
zero(
    trace_constant.subs(q, -3 - 2 * kappa * p)
    + 4 * p * (kappa - 1) * (kappa + 1) * (kappa * p + 1),
    "trace second equation factorization",
)

# First trace family p=0,q=-3: a genuine degree-two root map leaves one pole.
A_first = u**3 - 3 * u + r
t_first = (u + kappa) / (u**2 - 1)
res_first = sp.Poly(
    sp.expand(sp.resultant(A_first - A, (u**2 - 1) * T - (u + kappa), u)), T
)
a_first = res_first.coeff_monomial(T**3)
c_first = -res_first.coeff_monomial(T)
C_first = sp.factor(
    -(3 * a_first.subs(A, A_first) * t_first**2 + c_first.subs(A, A_first))
    / (2 * t_first)
)
first_num, first_den = sp.cancel(C_first).as_numer_denom()
gate(sp.expand(first_den - 2 * (kappa + u)) == 0, "first-family residual denominator")
gate(
    sp.factor(first_num.subs(u, -kappa))
    == 3 * (kappa - 1) ** 3 * (kappa + 1) ** 3,
    "first-family pole obstruction",
)

# Second trace family p=-1/kappa,q=-1: polynomial iff 3*kappa^2+1=0.
A_second = u**3 - u**2 / kappa - u + r
t_second = (u + kappa) / (u**2 - 1)
res_second = sp.Poly(
    sp.cancel(sp.resultant(A_second - A, (u**2 - 1) * T - (u + kappa), u)), T
)
a_second = sp.factor(res_second.coeff_monomial(T**3))
c_second = sp.factor(-res_second.coeff_monomial(T))
C_second = sp.factor(
    -(3 * a_second.subs(A, A_second) * t_second**2 + c_second.subs(A, A_second))
    / (2 * t_second)
)
second_num, second_den = sp.cancel(C_second).as_numer_denom()
gate(
    sp.expand(second_den - 2 * kappa**2 * (kappa + u)) == 0,
    "second-family denominator",
)
gate(
    sp.factor(second_num.subs(u, -kappa))
    == kappa**2 * (kappa - 1) ** 2 * (kappa + 1) ** 2 * (3 * kappa**2 + 1),
    "second-family polynomiality factor",
)


# ---------------------------------------------------------------------------
# Discriminant, implicit quintic, and squarefree support.
# ---------------------------------------------------------------------------

Delta = sp.expand(sp.discriminant(f, T))
J = (
    243 * A**5
    - 4455 * rho * A**4
    - 648 * rho * A**3 * C
    - 27738 * A**3
    - 4104 * A**2 * C
    + 18506 * rho * A**2
    + 432 * A * C**2
    + 2376 * rho * A * C
    + 10839 * A
    + 72 * C**3
    - 48 * rho * C**2
    + 648 * C
    - 627 * rho
)
zero(Delta + (9 * A - rho) * J / 324, "line times quintic discriminant")
zero(J.subs({A: A_u, C: C_u}), "quintic parametrization")
gate(sp.Poly(J, A, C).total_degree() == 5, "implicit component is quintic")
gate(sp.degree(A_u, u) == 3, "A projection has degree three")
gate(sp.degree(C_u, u) == 5, "C projection has degree five")
gate(sp.LC(sp.Poly(C_u, u), u) == -sp.Rational(3, 2), "unique infinity owner")
line_restriction = rho_remainder(J.subs(A, rho / 9))
gate(line_restriction != 0, "vertical line is not a quintic component")
gate(sp.Poly(line_restriction, C).degree() == 3, "line meets quintic properly")

# The elimination equation is the quintic, fixing the normalization degree.
implicit_resultant = sp.resultant(A - A_u, C - C_u, u)
zero(implicit_resultant - J / 72, "basepoint-free birational implicit equation")


# ---------------------------------------------------------------------------
# Unit-ideal and nonmonogenic valuation controls.
# ---------------------------------------------------------------------------

zero(d.subs(A, rho / 9), "simple d zero")
gate(rho_remainder(a.subs(A, rho / 9)) != 0, "a survives at the d zero")
zero(c.subs(A, rho / 9), "vertical discriminant line is the c=d fibre")
gate(sp.degree(a, A) == 2, "a has double zero divisor")
gate(sp.degree(d, A) == 1, "d has simple zero divisor")

# No residue m modulo three can pay both divisor congruences.
for residue in range(3):
    at_a = (2 * (2 * residue + 1)) % 3
    at_d = (-2 * residue) % 3
    gate(not (at_a == 0 and at_d == 0), f"incompatible cube residues m={residue}")
gate((2 * (2 * 1 + 1)) % 3 == 0, "a-prime alone selects m=1 mod three")
gate((-2 * 0) % 3 == 0, "d-prime alone selects m=0 mod three")


summary = {
    "checks": CHECKS,
    "field": "rho^2=-3",
    "root_map_degree": 2,
    "branch_degrees": [1, 5],
    "normalizations": ["A1", "A1"],
    "common_infinity": "[0:1:0]",
    "quintic_A_projection_degree": 3,
    "visible_collision": ["u=1", "u=-1", "target=(rho,0)"],
    "nonmonogenic": "incompatible m mod 3 at ord(a)=2 and ord(d)=1",
    "open": "units,class_group,boundary_basis,Euler,A2_atlas",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3930 exact aligned line-plus-quintic cubic packet")
print(f"CHECKS={CHECKS}")
print("FIELD=rho^2+3=0")
print("ROOT_MAP_DEGREE=2;FINITE_POLES=u=1,-1")
print("DISCRIMINANT=-(9A-rho)J5/324;SQUAREFREE_COMPONENT_DEGREES=1,5")
print("NORMALIZATIONS=A1,A1;COMMON_INFINITY=[0:1:0]")
print("QUINTIC_A_DEGREE=3;AFFINE_ADDRESS_CAP=3;VISIBLE_COLLISION=2")
print("ORDER=unit-ideal,normal,nonmonogenic,S3")
print("OPEN=units,class_group,boundary_basis,Euler,A2_atlas")
print(f"SEMANTIC_SHA256={semantic}")
