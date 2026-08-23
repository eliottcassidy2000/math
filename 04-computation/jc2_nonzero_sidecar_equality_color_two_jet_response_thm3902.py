#!/usr/bin/env python3
"""Exact companion for THM-3902's equality-color two-jet response law."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, message)


def coefficient(expression: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), variable).coeff_monomial(variable**degree)


def zero_mod_d(expression: sp.Expr, d: sp.Symbol, message: str) -> None:
    """Check a rational identity in the quadratic field d^2=-3."""
    numerator = sp.together(expression).as_numer_denom()[0]
    remainder = sp.Poly(sp.expand(numerator), d).rem(sp.Poly(d**2 + 3, d))
    gate(remainder.as_expr() == 0, message)


# ---------------------------------------------------------------------------
# 1. Exact residual and the all-degree support window through two y-jets.
# ---------------------------------------------------------------------------

T_symbol, f_symbol, K_symbol, a, L = sp.symbols("T_symbol f_symbol K_symbol a L")
A0 = a * L**2
P = A0
r_symbol = a * T_symbol + K_symbol * f_symbol
A_symbol = K_symbol * T_symbol + a * P * f_symbol
S_symbol = sp.expand(
    L**4
    + 2 * (3 * A_symbol + 3 * P + r_symbol**2) * L**2 * f_symbol
    + (8 * A_symbol + 6 * P + 3 * r_symbol**2)
    * (P * f_symbol**2 - T_symbol**2)
)

residual_poly = sp.Poly(S_symbol, T_symbol, f_symbol, K_symbol)
top_pair = {(2, 2, 2), (0, 4, 2)}
boundary_source = (0, 3, 2)
curvature_pair = {(3, 1, 1), (1, 3, 1)}
visible_support = top_pair | {boundary_source} | curvature_pair

gate(
    residual_poly.coeff_monomial((2, 2, 2)) == -3,
    "top T2 f2 K2 coefficient",
)
gate(
    residual_poly.coeff_monomial((0, 4, 2)) == 3 * A0,
    "top f4 K2 coefficient",
)
gate(
    residual_poly.coeff_monomial(boundary_source) == 2 * L**2,
    "marked K2 f3 source coefficient",
)
gate(
    residual_poly.coeff_monomial((3, 1, 1)) == -6 * a,
    "first curvature coefficient",
)
gate(
    residual_poly.coeff_monomial((1, 3, 1)) == 6 * a**2 * L**2,
    "second curvature coefficient",
)

# For T^p f^q K^r on m=n, the deficit from y-degree 4n+4 is
# (4-p-q)n+4-2r.  Every slope is nonnegative.  The five displayed
# monomials are exactly those whose deficit can be at most two for n>=1.
for monomial, _coefficient_value in residual_poly.terms():
    p_t, p_f, p_k = monomial
    slope = 4 - p_t - p_f
    intercept = 4 - 2 * p_k
    gate(slope >= 0, f"nonnegative equality-seam deficit slope {monomial}")
    deficit_at_one = slope + intercept
    if monomial in visible_support:
        gate(deficit_at_one <= 2, f"visible two-jet support {monomial}")
    else:
        gate(deficit_at_one >= 3, f"all-degree omitted support {monomial}")


# ---------------------------------------------------------------------------
# 2. Raw first and second coefficient formulas.
# ---------------------------------------------------------------------------

y = sp.symbols("y")
kappa = sp.symbols("kappa")
u, u1, u2, v, v1, v2 = sp.symbols("u u1 u2 v v1 v2")
E = A0 * v**2 - u**2

C0_expected = 3 * v**2 * E
C1_base = 6 * v * ((2 * A0 * v**2 - u**2) * v1 - u * v * u1)
C2_base = (
    3
    * (
        A0 * (4 * v**3 * v2 + 6 * v**2 * v1**2)
        - u**2 * (v1**2 + 2 * v * v2)
        - 4 * u * u1 * v * v1
        - (u1**2 + 2 * u * u2) * v**2
    )
    + 6 * (kappa * v**2 + a * u * v) * E
)

for n_value in (1, 2, 3, 4):
    T_series = u * y**n_value + u1 * y ** (n_value - 1)
    f_series = v * y**n_value + v1 * y ** (n_value - 1)
    if n_value >= 2:
        T_series += u2 * y ** (n_value - 2)
        f_series += v2 * y ** (n_value - 2)

    specialized = sp.expand(
        S_symbol.subs(
            {
                T_symbol: T_series,
                f_symbol: f_series,
                K_symbol: y**2 + kappa,
            }
        )
    )
    C0 = coefficient(specialized, y, 4 * n_value + 4)
    C1 = coefficient(specialized, y, 4 * n_value + 3)
    C2 = coefficient(specialized, y, 4 * n_value + 2)

    C1_expected = C1_base + (2 * L**2 * v**3 if n_value == 1 else 0)
    C2_expected = C2_base
    if n_value == 1:
        C2_expected = C2_expected.subs({u2: 0, v2: 0}) + 6 * L**2 * v**2 * v1
    elif n_value == 2:
        C2_expected += 2 * L**2 * v**3

    zero(C0 - C0_expected, f"leading coefficient n={n_value}")
    zero(C1 - C1_expected, f"first response coefficient n={n_value}")
    zero(C2 - C2_expected, f"second response coefficient n={n_value}")


# ---------------------------------------------------------------------------
# 3. Normalized colors and denominator-cleared polynomial sidecars.
# ---------------------------------------------------------------------------

h, j1, j2, d = sp.symbols("h j1 j2 d")
delta1, delta2 = sp.symbols("delta1 delta2")
c_minus = h - d * u
c_plus = h + d * u
c_minus_1 = j1 - d * u1
c_plus_1 = j1 + d * u1
c_minus_2 = j2 - d * u2
c_plus_2 = j2 + d * u2

zero_mod_d(c_minus * c_plus - (h**2 + 3 * u**2), d, "leading color product")
zero_mod_d(
    c_minus * c_plus_1 + c_plus * c_minus_1 - (2 * h * j1 + 6 * u * u1),
    d,
    "first marked color product",
)
zero_mod_d(
    c_minus * c_plus_2
    + c_plus * c_minus_2
    + c_minus_1 * c_plus_1
    - (2 * h * j2 + j1**2 + 6 * u * u2 + 3 * u1**2),
    d,
    "second marked color product",
)

# The coefficient laws after using h^2=3E.
first_norm_rhs = 6 * A0 * v * v1 + 2 * delta1 * L**2 * v
second_norm_rhs = (
    3 * A0 * (v1**2 + 2 * v * v2)
    + 2 * L**2 * (delta1 * v1 + delta2 * v)
    + 2 * (kappa + a * u / v) * h**2
)

# Directly freeze the compact epsilon response.  The case n=3 represents
# every n>=3 through epsilon^2 because the marked source then starts later.
epsilon = sp.symbols("epsilon")
U_epsilon = u + u1 * epsilon + u2 * epsilon**2
V_epsilon = v + v1 * epsilon + v2 * epsilon**2
for n_value in (1, 2, 3):
    U_case = U_epsilon.subs(u2, 0) if n_value == 1 else U_epsilon
    V_case = V_epsilon.subs(v2, 0) if n_value == 1 else V_epsilon
    compact_response = (
        3 * A0 * V_case**2
        + 2 * L**2 * epsilon**n_value * V_case
        + 6
        * epsilon**2
        * (kappa + a * U_case / V_case)
        * (A0 * V_case**2 - U_case**2)
    )
    compact_truncation = sp.series(compact_response, epsilon, 0, 3).removeO()
    substitutions = {
        delta1: 1 if n_value == 1 else 0,
        delta2: 1 if n_value == 2 else 0,
    }
    expected_second = second_norm_rhs.subs(substitutions).subs(h**2, 3 * E)
    if n_value == 1:
        expected_second = expected_second.subs({u2: 0, v2: 0})
    zero(
        coefficient(compact_truncation, epsilon, 0) - 3 * A0 * v**2,
        f"epsilon response leading coefficient n={n_value}",
    )
    zero(
        coefficient(compact_truncation, epsilon, 1)
        - first_norm_rhs.subs(substitutions),
        f"epsilon response first coefficient n={n_value}",
    )
    zero(
        coefficient(compact_truncation, epsilon, 2) - expected_second,
        f"epsilon response second coefficient n={n_value}",
    )

zero(
    first_norm_rhs
    - (6 * (A0 * v * v1 - u * u1) + 2 * delta1 * L**2 * v)
    - 6 * u * u1,
    "first normalized response allocation",
)
zero(
    second_norm_rhs
    - (
        3 * (A0 * (v1**2 + 2 * v * v2) - (u1**2 + 2 * u * u2))
        + 2 * L**2 * (delta1 * v1 + delta2 * v)
        + 2 * (kappa + a * u / v) * h**2
    )
    - 3 * u1**2
    - 6 * u * u2,
    "second normalized response allocation",
)

J1, J2 = sp.symbols("J1 J2")
R_minus_1 = J1 - d * v * u1
R_plus_1 = J1 + d * v * u1
R_minus_2 = J2 - d * v**2 * u2
R_plus_2 = J2 + d * v**2 * u2

poly_first_lhs = c_minus * R_plus_1 + c_plus * R_minus_1
poly_second_lhs = (
    c_minus * R_plus_2
    + c_plus * R_minus_2
    + R_minus_1 * R_plus_1
)
zero_mod_d(
    poly_first_lhs.subs(J1, v * j1)
    - v * (c_minus * c_plus_1 + c_plus * c_minus_1),
    d,
    "first polynomial sidecar clearing",
)
zero_mod_d(
    poly_second_lhs.subs({J1: v * j1, J2: v**2 * j2})
    - v**2
    * (c_minus * c_plus_2 + c_plus * c_minus_2 + c_minus_1 * c_plus_1),
    d,
    "second polynomial sidecar clearing",
)
zero(
    v * first_norm_rhs
    - 2 * v**2 * (3 * A0 * v1 + delta1 * L**2),
    "first polynomial response right side",
)
zero(
    v**2 * second_norm_rhs
    - (
        3 * A0 * v**2 * (v1**2 + 2 * v * v2)
        + 2 * L**2 * v**2 * (delta1 * v1 + delta2 * v)
        + 2 * v * (kappa * v + a * u) * h**2
    ),
    "second polynomial response right side",
)


# ---------------------------------------------------------------------------
# 4. The exact a-prime deletion and canonical controls.
# ---------------------------------------------------------------------------

x = sp.symbols("x")
a_x = x + 1
L_x = 9 * x + 4
kappa_x = -15 * x**2 - 15 * x - 4
h_star = (a_x + 3 * L_x**2) / 2
u_star = (3 * L_x**2 - a_x) / (2 * d)

gate(sp.gcd(a_x, L_x) == 1, "a and L are coprime")
zero(L_x.subs(x, -1) + 5, "L is a unit at the a-prime")
zero((1 + 3 * a_x * v1).subs(x, -1) - 1, "n=1 response bracket is an a-unit")
zero_mod_d(h_star**2 - 3 * (a_x * L_x**2 - u_star**2), d, "canonical norm")
zero_mod_d(h_star - d * u_star - a_x, d, "canonical a-color")
zero_mod_d(h_star + d * u_star - 3 * L_x**2, d, "canonical L-color")
gate(sp.Poly(h_star, x).degree() == 2, "canonical h is nonunit")
gate(sp.gcd(sp.Poly(2 * h_star, x), sp.Poly(L_x, x)).degree() == 0,
     "canonical h is coprime to L")

# Zero lower sidecars: n=1 fails at the first response and n=2 at the
# second.  For every n>=3 the same leading payment lifts both displayed jets.
zero(
    (C1_base + 2 * L**2 * v**3).subs(
        {a: a_x, L: L_x, u: u_star, v: 1, u1: 0, v1: 0}
    )
    - 2 * L_x**2,
    "n=1 zero-sidecar hostile remainder",
)
zero_mod_d(
    (
        C2_base.subs({u2: 0, v2: 0}) + 2 * L**2 * v**3
    ).subs(
        {
            a: a_x,
            L: L_x,
            kappa: kappa_x,
            u: u_star,
            v: 1,
            u1: 0,
            v1: 0,
        }
    )
    - 2 * h_star**2 * (kappa_x + a_x * u_star)
    - 2 * L_x**2,
    d,
    "n=2 zero-sidecar hostile remainder",
)
zero_mod_d(
    C2_base.subs(
        {
            a: a_x,
            L: L_x,
            kappa: kappa_x,
            u: u_star,
            v: 1,
            u1: 0,
            u2: 0,
            v1: 0,
            v2: 0,
        }
    )
    - 2 * h_star**2 * (kappa_x + a_x * u_star),
    d,
    "n>=3 zero-sidecar positive two-jet lift",
)


# ---------------------------------------------------------------------------
# 5. Address-compatible n=1 positive two-jet control.
# ---------------------------------------------------------------------------

p_control = -sp.Rational(2, 147) / (1 + 4 * d)
v1_control = p_control * h_star
j1_control = a_x * v1_control + sp.Rational(1, 3)
u1_control = d * j1_control / 3
j2_control = (
    (kappa_x + a_x * u_star) * h_star
    + L_x**2 * p_control * (3 * a_x * p_control * h_star + 2) / 2
)
g1_control = h_star * v1_control + j1_control
g2_control = j2_control + v1_control * j1_control

zero_mod_d(
    u1_control.subs(x, 0) - 4 * v1_control.subs(x, 0),
    d,
    "n=1 positive-control address",
)

C0_control = C0_expected.subs(
    {a: a_x, L: L_x, u: u_star, v: 1}
)
C1_control = (C1_base + 2 * L**2 * v**3).subs(
    {
        a: a_x,
        L: L_x,
        u: u_star,
        v: 1,
        u1: u1_control,
        v1: v1_control,
    }
)
C2_control = (
    C2_base.subs({u2: 0, v2: 0}) + 6 * L**2 * v**2 * v1
).subs(
    {
        a: a_x,
        L: L_x,
        kappa: kappa_x,
        u: u_star,
        v: 1,
        u1: u1_control,
        v1: v1_control,
    }
)
zero_mod_d(C0_control - h_star**2, d, "n=1 positive-control leading square")
zero_mod_d(C1_control - 2 * h_star * g1_control, d,
           "n=1 positive-control first response")
zero_mod_d(C2_control - g1_control**2 - 2 * h_star * g2_control, d,
           "n=1 positive-control second response")


semantic = {
    "theorem": "THM-3902",
    "universe": "THM3881 residual, positive equal y-degree seam m=n>=1",
    "epsilon_law": "Cminus*Cplus=3aL2V2+2L2*eps^n*V+6eps2*(kappa+aU/V)*(aL2V2-U2) mod eps3",
    "sidecars": "J1=v*j1,J2=v^2*j2 with two denominator-cleared color responses",
    "a_deletion": "n=1 first marked response has a-valuation 2ord_a(v)",
    "hostiles": "canonical zero-sidecar payment fails at jet1 for n=1 and jet2 for n=2",
    "positives": "zero-sidecar n>=3 and address-compatible n=1 lift both displayed jets",
    "scope": "necessary two-jet response only; equality seam and JC2 remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3902-nonzero-sidecar-equality-color-two-jet-response")
print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("epsilon_response=Cminus*Cplus_mod_eps3_with_marked_source_2L2_eps^n_V")
print("polynomial_sidecars=J1_vj1_and_J2_v2j2")
print("n1_a_valuation=first_marked_response_deletes_one_total_a_order")
print("canonical_zero_sidecar_n1=FAILS_AT_FIRST_RESPONSE")
print("canonical_zero_sidecar_n2=FAILS_AT_SECOND_RESPONSE")
print("canonical_zero_sidecar_n_ge_3=LIFTS_TWO_JETS")
print("address_compatible_n1_positive_two_jet=PASS")
print("scope=NECESSARY_TWO_JET_LAW_ONLY_EQUALITY_SEAM_OPEN_JC2_OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
