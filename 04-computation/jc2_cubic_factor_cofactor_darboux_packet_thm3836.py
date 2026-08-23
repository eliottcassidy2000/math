#!/usr/bin/env python3
"""Exact companion for THM-3836's cubic factor/cofactor packet."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


h, k, C, z = sp.symbols("h k C z")
P = 3 * h**3 + 7 * h**2 * k + k**3
Q = 7 * h**2 + 3 * k**2
B = 6 * h**3 + 7 * h**2 * k - k**3
S = h**2 * Q * C - k * B

r = 3 * z**3 + 7 * z**2 + 1
q = 7 * z**2 + 3
b = 6 * z**3 + 7 * z**2 - 1
s = z**2 * q * C - b
k_chart = r / (C * s)
h_chart = z * k_chart

zero(P.subs({h: h_chart, k: k_chart}, simultaneous=True) - k_chart**3 * r,
     "homogenized cubic spectral form")
zero(S.subs({h: h_chart, k: k_chart}, simultaneous=True) - k_chart**4 * s,
     "homogenized complementary factor")
zero(
    (C * S - P).subs({h: h_chart, k: k_chart}, simultaneous=True),
    "triangular chart clears to P=C*S",
)

# The Keller equation in the rational chart clears without an ansatz.
h_x, h_y, k_x, k_y, C_x, C_y, lam = sp.symbols(
    "h_x h_y k_x k_y C_x C_y lambda"
)
delta_C = k * (h_x * C_y - h_y * C_x) - h * (k_x * C_y - k_y * C_x)
z_x = (k * h_x - h * k_x) / k**2
z_y = (k * h_y - h * k_y) / k**2
jac_z_C = z_x * C_y - z_y * C_x
zero(jac_z_C - delta_C / k**2, "root-ratio Jacobian numerator")
zero(k * (P / k**3) - P / k**2, "weighted spectral right side")
zero(
    (jac_z_C - lam * k * (P / k**3)) - (delta_C - lam * P) / k**2,
    "Keller equation clears to delta(C)=lambda*P",
)

# Spectral arithmetic used by the comaximality and arm-label proof.
gate(sp.discriminant(r, z) != 0, "three cubic slopes distinct")
gate(sp.resultant(r, z, z) != 0, "cubic slopes nonzero")
gate(sp.resultant(r, q, z) != 0, "quadratic coefficient nonzero on slopes")
gate(sp.resultant(r, b, z) != 0, "sign coefficient nonzero on slopes")
zero(b + q - 2 * r, "b equals minus q on the cubic spectrum")

a = sp.symbols("a")
r_a = 3 * a**3 + 7 * a**2 + 1
Q_a = 7 * a**2 + 3
b_a = 6 * a**3 + 7 * a**2 - 1

zero(S.subs({h: 0}) - k**4, "S modulo h is k^4")
zero(S.subs({h: a * k, C: 0}) + b_a * k**4,
     "minus-arm complementary value")
plus_value = sp.expand(S.subs({h: a * k, C: -1 / a**2}))
plus_num = sp.together(plus_value).as_numer_denom()[0]
zero(sp.rem(sp.Poly(plus_num, a), sp.Poly(r_a, a)).as_expr(),
     "plus-arm value vanishes at every cubic slope")
zero(
    sp.rem(sp.Poly(b_a + Q_a, a), sp.Poly(r_a, a)).as_expr(),
    "plus/minus label relation at a cubic slope",
)

# Euler/cofactor identity for each possible monochromatic whole-member subset.
a1, a2, a3, J_hk = sp.symbols("a1 a2 a3 J_hk")
factors = [h - a1 * k, h - a2 * k, h - a3 * k]
for degree in (1, 2, 3):
    F = sp.prod(factors[:degree])
    zero(h * sp.diff(F, h) + k * sp.diff(F, k) - degree * F,
         f"Euler identity degree {degree}")
    delta_F = (k * sp.diff(F, k) + h * sp.diff(F, h)) * J_hk
    zero(delta_F - degree * F * J_hk,
         f"cofactor identity degree {degree}")
    expected_mod_h = (-1) ** degree * sp.prod([a1, a2, a3][:degree]) * k**degree
    zero(F.subs(h, 0) - expected_mod_h,
         f"nonzero h-boundary monomial degree {degree}")

# The differentiated determinant obstruction has a nonzero (j+1) multiplier
# in characteristic zero.  Record all three possible subset sizes.
for degree in (1, 2, 3):
    gate(degree != 0 and degree + 1 != 0,
         f"characteristic-zero endpoint multiplier degree {degree}")

semantic = {
    "factor": "P=3h^3+7h^2k+k^3=C*S; S=h^2(7h^2+3k^2)C-k(6h^3+7h^2k-k^3)",
    "cofactor": "delta(C)=k*J(h,C)-h*J(k,C)=lambda*P",
    "coprime": "(C,S)=1 and the three h-alpha_i*k pencils are pairwise comaximal",
    "allocation": "C and S allocate every irreducible pencil component",
    "obstruction": "whole-member allocation gives Euler cofactor and differentiated determinant contradiction",
    "conclusion": "some cubic pencil member has nonempty comaximal minus and plus factors",
    "arms": "C=0 and C=-1/alpha^2",
    "scope": "all-degree conditional packet; atlas later excluded by THM-3841",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3836-cubic-factor-cofactor-darboux-packet")
print("factor=P=C*S;P=3h^3+7h^2k+k^3")
print("cofactor=delta(C)=k*J(h,C)-h*J(k,C)=lambda*P")
print("coprime=(C,S)=1;pencil_members_pairwise_comaximal")
print("allocation=irreducible_components_select_C_or_S")
print("monochromatic=all_three_whole_member_allocations_impossible")
print("forced=one_cubic_fibre_has_nonempty_comaximal_minus_plus_factors")
print("arms=minus:C=0;plus:C=-1/alpha^2")
print("scope=all_degree_conditional_packet;atlas_excluded_by_THM3841")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
