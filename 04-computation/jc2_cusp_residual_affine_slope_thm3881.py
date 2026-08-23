#!/usr/bin/env python3
"""Exact audit for the affine-slope closure in THM-3881's rank-two norm."""

from __future__ import annotations

import hashlib
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.expand(expression) == 0, label)


def homogeneous_part(expression: sp.Expr, degree: int, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    scale = sp.symbols("scale")
    scaled = sp.Poly(sp.expand(expression.subs({x: scale * x, y: scale * y})), scale)
    return sp.expand(scaled.coeff_monomial(scale**degree))


x, y, F = sp.symbols("x y F")
alpha, beta, gamma = sp.symbols("alpha beta gamma")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
P = a * L**2
h = alpha * x + beta * y + gamma
T = h * F
r = a * T + K * F
A = K * T + a * P * F
B = P * F**2 - T**2
S = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * F
    + (8 * A + 6 * P + 3 * r**2) * B
)

R = a * h + K
M = K * h + a * P
expected_coefficients = {
    0: L**4,
    1: 6 * P * L**2,
    2: 6 * M * L**2 + 6 * P * (P - h**2),
    3: 2 * R**2 * L**2 + 8 * M * (P - h**2),
    4: 3 * R**2 * (P - h**2),
}
poly_F = sp.Poly(S, F)
gate(poly_F.degree() == 4, "quartic in the module coordinate")
for degree, expected in expected_coefficients.items():
    zero(poly_F.coeff_monomial(F**degree) - expected, f"F coefficient {degree}")

h1 = alpha * x + beta * y
K2 = y**2 - 15 * x**2
P3 = 81 * x**3
R2 = K2 + x * h1
M4 = 81 * x**4
leading_expected = {
    4: 3 * R2**2 * P3,
    3: 8 * M4 * P3,
}
zero(
    homogeneous_part(expected_coefficients[4], 7, x, y) - leading_expected[4],
    "degree-seven F4 leading form",
)
zero(
    homogeneous_part(expected_coefficients[3], 7, x, y) - leading_expected[3],
    "degree-seven F3 leading form",
)
for degree, cap in ((4, 7), (3, 7), (2, 6), (1, 5), (0, 4)):
    coefficient = expected_coefficients[degree]
    gate(sp.Poly(coefficient, x, y).total_degree() <= cap, f"coefficient degree cap F{degree}")

# For deg(f)=n>=1, the F4 contribution has degree 4n+7, while every
# lower-F contribution is strictly smaller.  Freeze the complete integer
# ledger over a hostile range of n rather than relying on prose arithmetic.
degree_rows = []
for n in range(1, 101):
    degrees = (4 * n + 7, 3 * n + 7, 2 * n + 6, n + 5, 4)
    gate(degrees[0] % 2 == 1, f"odd top degree n={n}")
    gate(degrees[0] > max(degrees[1:]), f"unique top degree n={n}")
    degree_rows.append((n, degrees[0], max(degrees[1:])))

# The constant nonzero f=c edge has two degree-seven contributions.  Their
# sum cannot cancel: R2 has y^2 coefficient one, so the x^3*y^4 coefficient
# below is 243*c^4.
c = sp.symbols("c", nonzero=True)
constant_leading = sp.expand(c**4 * leading_expected[4] + c**3 * leading_expected[3])
constant_formula = sp.expand(P3 * c**3 * (3 * c * R2**2 + 648 * x**4))
zero(constant_leading - constant_formula, "constant-f degree-seven form")
gate(
    sp.Poly(constant_leading, x, y).coeff_monomial(x**3 * y**4) == 243 * c**4,
    "constant-f noncancellation",
)

# Direct hostile controls with nonhomogeneous affine slopes and several f.
controls = (
    (0, x),
    (3, x + y + 1),
    (2 * x - 5 * y + 7, x**2 + y + 1),
    (-x + 4 * y - 2, x**3 - 2 * x * y + y**2 + 3),
)
control_rows = []
for index, (h_value, f_value) in enumerate(controls):
    direct = sp.expand(S.subs({alpha: sp.expand(h_value).coeff(x), beta: sp.expand(h_value).coeff(y), gamma: sp.expand(h_value).subs({x: 0, y: 0}), F: f_value}))
    n = sp.Poly(f_value, x, y).total_degree()
    expected_degree = 4 * n + 7 if n >= 1 else 7
    gate(sp.Poly(direct, x, y).total_degree() == expected_degree, f"direct degree control {index}")
    gate(expected_degree % 2 == 1, f"direct odd control {index}")
    control_rows.append((index, n, expected_degree))

semantic = (
    "THM-3881 rank-two residual restricted to T=h f",
    "h arbitrary affine-linear polynomial",
    "nonconstant f has unique top degree 4deg(f)+7",
    "constant nonzero f has noncancelling x^3 y^4 term",
    "polynomial squares have even total degree",
)
semantic_sha256 = hashlib.sha256(repr(semantic).encode()).hexdigest()
print("JC2_RANK_TWO_AFFINE_SLOPE_PRIMARY_20260823")
print("status=PROVED_DEGREE_OBSTRUCTION;JC2_OPEN")
print("coefficient_degree_caps=F4:7,F3:7,F2:6,F1:5,F0:4")
print("nonconstant_degree_rule=degS=4degf+7")
print("constant_edge=x3y4_coefficient=243*c^4")
print("degree_rows=" + repr((degree_rows[0], degree_rows[-1])).replace(" ", ""))
print("control_rows=" + repr(tuple(control_rows)).replace(" ", ""))
print(f"semantic_sha256={semantic_sha256}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
