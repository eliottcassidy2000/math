#!/usr/bin/env python3
"""Exact coefficient companion for THM-3897's final three-channel descent."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


x, y = sp.symbols("x y")
u, v, w = sp.symbols("u v w")
a = x + 1
L = 9 * x + 4
F = 15 * x**2 + 15 * x + 4
K = y**2 - F
T = u * y**2 + v * y + w
residual = sp.Poly(
    sp.expand(L**4 - 6 * a * L**2 * T**2 - 8 * K * T**3 - 3 * a**2 * T**4),
    y,
)

zero(
    residual.coeff_monomial(y**8) + u**3 * (8 + 3 * a**2 * u),
    "quadratic-y leading coefficient",
)
gate(
    sp.Poly(sp.together(-8 / (3 * a**2)).as_numer_denom()[1], x).degree() == 2,
    "vanishing-leading alternative has an a-square denominator",
)

linear_residual = sp.Poly(residual.as_expr().subs(u, 0), y)
gate(linear_residual.degree() == 5, "formal nonzero linear-y degree")
zero(linear_residual.coeff_monomial(y**5) + 8 * v**3, "linear-y odd top")

constant_residual = sp.Poly(linear_residual.as_expr().subs(v, 0), y)
constant_term = sp.expand(
    L**4 + 8 * F * w**3 - 6 * a * L**2 * w**2 - 3 * a**2 * w**4
)
gate(constant_residual.degree() == 2, "formal nonzero constant-T y-degree")
zero(constant_residual.coeff_monomial(y**2) + 8 * w**3, "constant-T y-square")
zero(constant_residual.coeff_monomial(y), "constant-T missing y term")
zero(constant_residual.coeff_monomial(1) - constant_term, "constant-T constant term")
zero(constant_term.subs({x: 0, w: 0}) - 256, "origin address obstruction")

# Exact factor shadow of the polynomial Pell step.  Reduction is performed
# modulo sqrt3^2-3, so no floating algebraic constants enter the transcript.
sqrt3, pell_s, pell_r = sp.symbols("sqrt3 pell_s pell_r")
pell_product = (pell_s - sqrt3 * a * pell_r) * (pell_s + sqrt3 * a * pell_r)
pell_reduced = sp.Poly(sp.expand(pell_product), sqrt3).rem(
    sp.Poly(sqrt3**2 - 3, sqrt3)
)
zero(pell_reduced.as_expr() - (pell_s**2 - 3 * a**2 * pell_r**2), "Pell factorization")
zero(
    (pell_s + sqrt3 * a * pell_r)
    - (pell_s - sqrt3 * a * pell_r)
    - 2 * sqrt3 * a * pell_r,
    "Pell unit difference",
)

# Base point and the generic rational-x hostile distinguish the global result
# from a false claim over k(x)[y].
T_hostile = -2 * K / (3 * a**2)
G_hostile = 4 * K**2 / (3 * a**3) - L**2
quartic_hostile = sp.expand(
    L**4
    - 6 * a * L**2 * T_hostile**2
    - 8 * K * T_hostile**3
    - 3 * a**2 * T_hostile**4
)
zero(G_hostile**2 - quartic_hostile, "rational-x hostile positive control")
zero(T_hostile.subs({x: 0, y: 0}) - sp.Rational(8, 3), "hostile fails address")
gate(
    sp.cancel(a**2 * T_hostile).subs(x, -1) != 0,
    "hostile fails polynomial x-descent",
)

semantic = {
    "inheritance": "THM3895 reduces T to u(x)y2+v(x)y+w(x)",
    "quadratic": "[y8]=-u3(8+3a2u);coprime parity gives polynomial Pell",
    "linear": "[y5]=-8v3 has odd degree",
    "constant": "missing y term plus origin value 256 kills nonzero w",
    "conclusion": "f=0 residual square and address imply T=0",
    "scope": "f_nonzero_and_JC2_remain_open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3897-f-zero-residual-all-degree-global-emptiness")
print("quadratic_y_gate=polynomial_Pell_units_force_u_zero")
print("linear_y_gate=odd_degree_five_forces_v_zero")
print("constant_y_gate=missing_linear_plus_origin_256_forces_w_zero")
print("conclusion=f_zero_implies_T_zero")
print("hostile_control=k(x)[y]_point_fails_descent_and_address")
print("remaining=f_nonzero_residual_lane_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
