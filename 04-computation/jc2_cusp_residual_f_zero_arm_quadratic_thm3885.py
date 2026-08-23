#!/usr/bin/env python3
"""Exact companion for THM-3885's f=0 arm and quadratic packet."""

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
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
T = sp.symbols("T")

# THM-3881's full residual equation at f=0.
residual = sp.expand(
    L**4 - 6 * a * L**2 * T**2 - 8 * K * T**3 - 3 * a**2 * T**4
)
zero(
    residual
    - (
        L**4
        + 2 * (3 * K * T + 3 * a * L**2 + a**2 * T**2) * 0
        + (8 * K * T + 6 * a * L**2 + 3 * a**2 * T**2) * (-T**2)
    ),
    "f-zero specialization of THM-3881 residual",
)

# The a=0 arm gives a polynomial Catalan equation.
arm = sp.expand(residual.subs(x, -1))
zero(arm - (625 - 8 * (y**2 - 4) * T**3), "a-zero arm identity")
c, alpha, beta = sp.symbols("c alpha beta")
arm_constant = sp.expand(arm.subs(T, c))
zero(
    arm_constant - (625 + 32 * c**3 - 8 * c**3 * y**2),
    "constant arm polynomial",
)
constant_square_groebner = sp.groebner(
    [
        alpha**2 + 8 * c**3,
        alpha * beta,
        beta**2 - 625 - 32 * c**3,
    ],
    alpha,
    beta,
    c,
    order="lex",
)
expected_constant_square_groebner = sp.groebner(
    [
        alpha**2 + 8 * c**3,
        alpha * beta,
        alpha * (32 * c**3 + 625),
        beta**2 - 32 * c**3 - 625,
        beta * c**3,
        c**3 * (32 * c**3 + 625),
    ],
    alpha,
    beta,
    c,
    order="lex",
)
gate(
    constant_square_groebner == expected_constant_square_groebner,
    "constant arm square classification",
)
zero(arm_constant.subs(c, 0) - 25**2, "zero arm positive control")
zero(
    arm_constant.subs(c**3, -sp.Rational(625, 32))
    - (sp.Rational(25, 2) * y) ** 2,
    "exceptional arm positive control",
)

# Every total-degree-at-most-two survivor has the following form after the arm
# theorem.  The x=0 odd coefficient kills q; y-degree then kills the rest.
p, q = sp.symbols("p q")
T_quadratic = sp.expand(c + a * (p * x + q * y - c))
zero(T_quadratic.subs({x: 0, y: 0}), "quadratic address")
zero(T_quadratic.subs(x, -1) - c, "quadratic arm constant")
quadratic_residual = sp.expand(residual.subs(T, T_quadratic))
quadratic_x_zero = sp.Poly(quadratic_residual.subs(x, 0), y)
gate(quadratic_x_zero.degree() == 5, "quadratic x-zero symbolic degree")
zero(
    quadratic_x_zero.coeff_monomial(y**5) + 8 * q**3,
    "quadratic x-zero odd coefficient",
)
T_one_color = sp.expand(T_quadratic.subs(q, 0))
gate(not T_one_color.has(y), "quadratic one-color reduction")
one_color_residual = sp.Poly(
    sp.expand(residual.subs(T, T_one_color)), y
)
gate(one_color_residual.degree() == 2, "quadratic one-color y degree")
zero(one_color_residual.coeff_monomial(y), "quadratic missing y-linear term")
zero(
    one_color_residual.coeff_monomial(y**2) + 8 * T_one_color**3,
    "quadratic unavoidable y-square coefficient",
)
zero(
    one_color_residual.coeff_monomial(1).subs(x, 0) - 256,
    "quadratic nonzero constant part",
)

# At L=0 the residual is a product of two coprime square-class buckets after
# removing d=gcd(tau,K0).  The last identity is the complete root-polarization
# parametrization of v^2-gamma^2*u^2=8e.
tau = sp.symbols("tau")
x_L = -sp.Rational(4, 9)
a_L = sp.Rational(5, 9)
K_L = y**2 - sp.Rational(8, 27)
L_arm = sp.expand(residual.subs({x: x_L, T: tau}))
zero(
    L_arm + tau**3 * (8 * K_L + sp.Rational(25, 27) * tau),
    "L-zero product identity",
)
d, e, sigma = sp.symbols("d e sigma")
zero(
    -(d * sigma) ** 3
    * (8 * d * e + sp.Rational(25, 27) * d * sigma)
    + d**4 * sigma**3 * (8 * e + sp.Rational(25, 27) * sigma),
    "L-zero gcd extraction",
)
lam, e_minus, e_plus = sp.symbols("lam e_minus e_plus", nonzero=True)
gamma = 5 / (3 * sp.sqrt(3))
u = sp.expand((8 * e_plus / lam - lam * e_minus) / (2 * gamma))
v = sp.expand((8 * e_plus / lam + lam * e_minus) / 2)
zero(
    v**2 - gamma**2 * u**2 - 8 * e_minus * e_plus,
    "L-zero root polarization",
)
polarized_residual = sp.expand(
    -d**4 * u**6 * (8 * e_minus * e_plus + gamma**2 * u**2)
)
zero(
    polarized_residual - (sp.sqrt(-1) * d**2 * u**3 * v) ** 2,
    "L-zero polarized square control",
)


semantic = {
    "lane": "f=0 residual L4-6aL2T2-8KT3-3a2T4",
    "a_arm": "T(-1,y)=c; c=0 or c3=-625/32",
    "quadratic": "address+arm+x0 odd coefficient+y-degree force T=0",
    "L_arm": "tau=d*u2;v2-(25/27)u2=8e;root partition of e",
    "scope": "nonlinear T=c+aU with deg U at least two remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure")
print("a_zero_arm_values=c=0_or_c^3=-625/32")
print("f_zero_total_degree_at_most_two_survivors=0")
print("L_zero_boundary=root_polarization_of_squarefree_divisor")
print("nonlinear_f_zero_lane=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
