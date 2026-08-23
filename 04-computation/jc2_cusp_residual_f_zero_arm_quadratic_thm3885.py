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

# Every total-degree-at-most-three survivor has the following complete form
# after the a=0 arm theorem.  Its x=0 restriction has degree at most two.
p, q = sp.symbols("p q")
alpha_c, beta_c, gamma_c = sp.symbols("alpha_c beta_c gamma_c")
U_cubic = (
    alpha_c * x**2
    + beta_c * x * y
    + q * y**2
    + gamma_c * x
    + p * y
    - c
)
T_cubic = sp.expand(c + a * U_cubic)
zero(T_cubic.subs({x: 0, y: 0}), "cubic address")
zero(T_cubic.subs(x, -1) - c, "cubic a-zero arm constant")
tau_x_zero = sp.expand(T_cubic.subs(x, 0))
zero(tau_x_zero - (p * y + q * y**2), "cubic x-zero restriction")
x_zero_residual = sp.expand(residual.subs({x: 0, T: tau_x_zero}))
zero(
    x_zero_residual
    - (256 - 96 * tau_x_zero**2 - 8 * (y**2 - 4) * tau_x_zero**3
       - 3 * tau_x_zero**4),
    "cubic x-zero residual identity",
)

# Fix the square-root sign by G(0)=16.  Coefficients y^2,y^3,y^4 determine
# the displayed recurrence; the next three coefficients give E5,E6,E7.
g2 = -3 * p**2
g3 = p * (p**2 - 6 * q)
g4 = -sp.Rational(3, 8) * (p**4 - 8 * p**2 * q + 8 * q**2)
G_x_zero = 16 + g2 * y**2 + g3 * y**3 + g4 * y**4
x_zero_difference = sp.Poly(
    sp.expand(G_x_zero**2 - x_zero_residual), y
)
for degree in range(5):
    zero(
        x_zero_difference.coeff_monomial(y**degree),
        f"cubic square-root recurrence through degree {degree}",
    )
z = sp.symbols("z")
E5 = 3 * z**2 - 24 * z * q - 4 * z + 48 * q**2
E6 = (
    13 * z**3
    - 120 * z**2 * q
    + 288 * z * q**2
    + 96 * z * q
    - 128 * q**3
)
E7 = z**3 - 14 * z**2 * q + 56 * z * q**2 - 64 * q**3 - 32 * q**2
zero(
    x_zero_difference.coeff_monomial(y**5)
    + 2 * p * E5.subs(z, p**2),
    "cubic y5 remainder",
)
zero(
    x_zero_difference.coeff_monomial(y**6)
    - sp.Rational(1, 4) * E6.subs(z, p**2),
    "cubic y6 remainder",
)
zero(
    x_zero_difference.coeff_monomial(y**7)
    + sp.Rational(3, 4) * p * E7.subs(z, p**2),
    "cubic y7 remainder",
)
zero(
    x_zero_difference.coeff_monomial(y**6).subs(p, 0)
    + 32 * q**3,
    "cubic p-zero seam forces q",
)
cubic_remainder_groebner = sp.groebner([E5, E6, E7], z, q, order="lex")
expected_cubic_remainder_groebner = sp.groebner(
    [z, q**2], z, q, order="lex"
)
gate(
    cubic_remainder_groebner == expected_cubic_remainder_groebner,
    "cubic nonzero-p remainder ideal",
)

# Once p=q=0, the L=0 arm is affine in y.  Its odd leading coefficient kills
# beta_c.  The remaining global polynomial depends only on x and is killed by
# the same y-degree-two/no-y-linear-term argument as the quadratic subcell.
T_after_x_zero = sp.expand(T_cubic.subs({p: 0, q: 0}))
A_L_cubic = (
    sp.Rational(4, 9) * c
    + sp.Rational(80, 729) * alpha_c
    - sp.Rational(20, 81) * gamma_c
)
B_L_cubic = -sp.Rational(20, 81) * beta_c
zero(
    T_after_x_zero.subs(x, -sp.Rational(4, 9))
    - (A_L_cubic + B_L_cubic * y),
    "cubic L-zero affine restriction",
)
affine_L_residual = sp.Poly(
    sp.expand(residual.subs({x: -sp.Rational(4, 9), T: A_L_cubic + B_L_cubic * y})),
    y,
)
gate(affine_L_residual.degree() == 5, "cubic L-zero affine symbolic degree")
zero(
    affine_L_residual.coeff_monomial(y**5) + 8 * B_L_cubic**3,
    "cubic L-zero odd leading coefficient",
)
A_L_symbol = sp.symbols("A_L_symbol")
constant_L_residual = sp.expand(
    residual.subs({x: -sp.Rational(4, 9), T: A_L_symbol})
)
zero(
    constant_L_residual
    + A_L_symbol**3
    * (8 * y**2 + (25 * A_L_symbol - 64) / 27),
    "cubic L-zero constant classification identity",
)
zero(
    constant_L_residual.subs(A_L_symbol, sp.Rational(64, 25))
    + 8 * sp.Rational(64, 25) ** 3 * y**2,
    "cubic L-zero exceptional positive control",
)
T_one_color = sp.expand(T_after_x_zero.subs(beta_c, 0))
gate(not T_one_color.has(y), "cubic one-color reduction")
one_color_residual = sp.Poly(sp.expand(residual.subs(T, T_one_color)), y)
gate(one_color_residual.degree() == 2, "cubic one-color y degree")
zero(one_color_residual.coeff_monomial(y), "cubic missing y-linear term")
zero(
    one_color_residual.coeff_monomial(y**2) + 8 * T_one_color**3,
    "cubic unavoidable y-square coefficient",
)
zero(
    one_color_residual.coeff_monomial(1).subs(x, 0) - 256,
    "cubic nonzero constant part",
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
    "cubic": "x0 recurrence+L0 odd degree+global y-degree force T=0",
    "L_arm": "tau=d*u2;v2-(25/27)u2=8e;root partition of e",
    "scope": "nonlinear T=c+aU with deg U at least three remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure")
print("a_zero_arm_values=c=0_or_c^3=-625/32")
print("f_zero_total_degree_at_most_three_nonzero_survivors=0")
print("L_zero_boundary=root_polarization_of_squarefree_divisor")
print("nonlinear_f_zero_lane=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
