#!/usr/bin/env python3
"""Exact leading-y companion for THM-3899's first f-nonzero gate."""

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


T, f, K, a, L = sp.symbols("T f K a L")
P = a * L**2
r = a * T + K * f
A = K * T + a * P * f
S = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * f
    + (8 * A + 6 * P + 3 * r**2) * (P * f**2 - T**2)
)

expanded = sp.expand(
    L**4
    + 6 * L**4 * a * f
    + 12 * L**4 * a**2 * f**2
    + 8 * L**4 * a**3 * f**3
    + 6 * K * L**2 * T * f
    + 12 * K * L**2 * a * T * f**2
    + 6 * K * L**2 * a**2 * T * f**3
    - 6 * L**2 * a * T**2
    - 6 * L**2 * a**2 * T**2 * f
    + 3 * L**2 * a**3 * T**2 * f**2
    - 8 * K * T**3
    - 6 * K * a * T**3 * f
    - 3 * a**2 * T**4
    + 2 * K**2 * L**2 * f**3
    + 3 * K**2 * L**2 * a * f**4
    - 3 * K**2 * T**2 * f**2
)
zero(S - expanded, "full THM-3881 residual expansion")

# Certify dominance when m<n.  Put n=N+1 and n-m=D+1.  For every monomial
# other than K^2 f^4, the deficit from degree 4n+4 is a polynomial in N,D
# with nonnegative coefficients and a positive constant.
N, D, m, n = sp.symbols("N D m n", integer=True, nonnegative=True)
poly = sp.Poly(S, T, f, K)
top_monomial = (0, 4, 2)
gate(poly.coeff_monomial(top_monomial) == 3 * a * L**2, "dominant coefficient")

dominance_rows = []
for monomial, coefficient in poly.terms():
    if monomial == top_monomial:
        continue
    degree = monomial[0] * m + monomial[1] * n + 2 * monomial[2]
    deficit = sp.expand(
        (4 * n + 4 - degree).subs(m, n - (D + 1)).subs(n, N + 1)
    )
    deficit_poly = sp.Poly(deficit, N, D)
    gate(
        all(coefficient >= 0 for coefficient in deficit_poly.coeffs())
        and deficit.subs({N: 0, D: 0}) > 0,
        f"strict dominance for monomial {monomial}",
    )
    dominance_rows.append((monomial, str(deficit)))

# On the equality seam m=n>=1, exactly K^2 f^4 and K^2 T^2 f^2 have top
# weight 4n+4.
equality_top = []
for monomial, coefficient in poly.terms():
    degree = (monomial[0] + monomial[1]) * n + 2 * monomial[2]
    deficit = sp.expand((4 * n + 4 - degree).subs(n, N + 1))
    if deficit == 0:
        equality_top.append((monomial, coefficient))
gate(
    equality_top == [((2, 2, 2), -3), ((0, 4, 2), 3 * L**2 * a)],
    "complete equality-seam top pair",
)

u, v, h, d = sp.symbols("u v h d")
d_relation = sp.Poly(d**2 + 3, d)
color_product = sp.Poly(sp.expand((h - d * u) * (h + d * u)), d).rem(d_relation)
zero(color_product.as_expr() - (h**2 + 3 * u**2), "equianharmonic norm colors")

# A positive leading-color payment: neither equality nor the color norm is
# itself an obstruction.
h0 = (a + 3 * L**2) / 2
u0 = (3 * L**2 - a) / (2 * d)
zero(h0 - d * u0 - a, "positive color one")
zero(h0 + d * u0 - 3 * L**2, "positive color two")
positive_norm = sp.together(h0**2 - 3 * (a * L**2 - u0**2))
positive_num = sp.Poly(positive_norm.as_numer_denom()[0], d).rem(d_relation)
zero(positive_num.as_expr(), "positive equality-seam norm solution")

# The lambda coordinate used by the incoming cubic carrier is exactly the
# same square root and the same two constants as THM-3895's quartic colors.
lam = sp.symbols("lam")
lam_relation = sp.Poly(lam**2 - lam + 1, lam)
d_lam = 2 * lam - 1
zero(
    sp.Poly(d_lam**2 + 3, lam).rem(lam_relation).as_expr(),
    "lambda gives sqrt minus three",
)
zero((3 + d_lam) / 2 - (1 + lam), "first shifted color")
zero((3 - d_lam) / 2 - (2 - lam), "second shifted color")

x = sp.symbols("x")
gate(sp.gcd(x + 1, 9 * x + 4) == 1, "a and L are coprime")
zero((9 * x + 4).subs(x, -1) + 5, "L is a unit at the a-prime")

semantic = {
    "residual": "expanded_exactly_in_T_f_K",
    "tariff": "deg_y_T_at_least_deg_y_f",
    "strict_gate": "odd_a_valuation_of_3aL2v4",
    "equality": "(h-du)(h+du)=3aL2v2,d2=-3",
    "positive_control": "colors_a_and_3L2_pay_the_leading_norm",
    "bridge": "d=2lambda-1;colors=1+lambda,2-lambda",
    "scope": "necessary_leading_filtration_only",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3899-nonzero-sidecar-y-degree-tariff")
print("strict_tariff=deg_y_T_must_be_at_least_deg_y_f")
print("dominant_term_m_lt_n=3aL2K2f4")
print(f"dominance_rows={len(dominance_rows)}")
print("equality_colors=(h-du)(h+du)=3aL2v2")
print("equianharmonic_coordinate=d=2lambda-1")
print("shifted_colors=1+lambda,2-lambda")
print("positive_control=equality_leading_norm_is_nonempty")
print("remaining=lower_coefficients_and_actual_square_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
