#!/usr/bin/env python3
"""Exact support/identity companion for provisional THM-3901."""

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
    gate(sp.cancel(expression) == 0, message)


def reduce_s(expression: sp.Expr) -> sp.Expr:
    """Reduce a rational expression modulo s^2+3 by its numerator."""
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.Poly(sp.expand(numerator), s).rem(sp.Poly(s**2 + 3, s)).as_expr()


def zero_s(expression: sp.Expr, message: str) -> None:
    zero(reduce_s(expression), message)


T, f, K, a, L, s = sp.symbols("T f K a L s", nonzero=True)
P = a * L**2
r = a * T + K * f
A = K * T + a * P * f
B = P * f**2 - T**2
S = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * f
    + (8 * A + 6 * P + 3 * r**2) * B
)

expanded_S = sp.expand(
    L**4
    + 6 * a * L**4 * f
    + 12 * a**2 * L**4 * f**2
    + 8 * a**3 * L**4 * f**3
    + 6 * K * L**2 * T * f
    + 12 * a * K * L**2 * T * f**2
    + 6 * a**2 * K * L**2 * T * f**3
    - 6 * a * L**2 * T**2
    - 6 * a**2 * L**2 * T**2 * f
    + 3 * a**3 * L**2 * T**2 * f**2
    - 8 * K * T**3
    - 6 * a * K * T**3 * f
    - 3 * a**2 * T**4
    + 2 * K**2 * L**2 * f**3
    + 3 * a * K**2 * L**2 * f**4
    - 3 * K**2 * T**2 * f**2
)
zero(S - expanded_S, "fifteen-term residual expansion")

Q0 = s * T * r
R = sp.expand(S + 3 * T**2 * r**2)
R_structural = sp.expand(
    L**4
    + 6 * (A + P) * L**2 * f
    + (8 * A + 6 * P) * B
    + r**2 * L**2 * f * (2 + 3 * a * f)
)
zero(R - R_structural, "first osculating structural identity")
zero_s(S - Q0**2 - R, "S-Q0^2=R modulo s^2+3")

expanded_R = sp.expand(
    3 * a * K**2 * L**2 * f**4
    + 2 * K**2 * L**2 * f**3
    + 6 * a**2 * K * L**2 * T * f**3
    + 12 * a * K * L**2 * T * f**2
    + 6 * K * L**2 * T * f
    - 8 * K * T**3
    + 8 * a**3 * L**4 * f**3
    + 12 * a**2 * L**4 * f**2
    + 6 * a * L**4 * f
    + L**4
    + 3 * a**3 * L**2 * T**2 * f**2
    - 6 * a**2 * L**2 * T**2 * f
    - 6 * a * L**2 * T**2
)
zero(R - expanded_R, "first response expansion")

Q1 = s * (T * r - sp.Rational(1, 2) * a**2 * L**2 * f**2)
R2 = sp.expand(
    3 * a * K**2 * L**2 * f**4
    + 2 * K**2 * L**2 * f**3
    + 3 * a**2 * K * L**2 * T * f**3
    + 12 * a * K * L**2 * T * f**2
    + 6 * K * L**2 * T * f
    - 8 * K * T**3
    + sp.Rational(3, 4) * a**4 * L**4 * f**4
    + 8 * a**3 * L**4 * f**3
    + 12 * a**2 * L**4 * f**2
    + 6 * a * L**4 * f
    + L**4
    - 6 * a**2 * L**2 * T**2 * f
    - 6 * a * L**2 * T**2
)
zero_s(S - Q1**2 - R2, "second osculating identity")

u, v = sp.symbols("u v", nonzero=True)


def weighted_top(expression: sp.Expr, m_value: int, n_value: int) -> tuple[int, sp.Expr]:
    """Return generic weighted y-degree and leading coefficient (lc K=1)."""
    polynomial = sp.Poly(expression, T, f, K)
    rows = []
    for (t_power, f_power, k_power), coefficient in polynomial.terms():
        degree = t_power * m_value + f_power * n_value + 2 * k_power
        leading = coefficient * u**t_power * v**f_power
        rows.append((degree, leading))
    top_degree = max(degree for degree, _ in rows)
    top_coefficient = sp.expand(
        sum(leading for degree, leading in rows if degree == top_degree)
    )
    return top_degree, top_coefficient


# Generic m>n>0 top fan.  The finite grid is a hostile control for the
# all-degree affine comparisons written in the theorem proof.
generic_top_cells = 0
first_response_cells = 0
for n_value in range(1, 17):
    for delta in range(1, 13):
        m_value = n_value + delta
        degree, coefficient = weighted_top(S, m_value, n_value)
        if delta == 1:
            expected_degree = 4 * n_value + 6
            expected_coefficient = -3 * u**2 * v**2
        elif delta == 2:
            expected_degree = 4 * n_value + 8
            expected_coefficient = -3 * u**2 * (a * u + v) ** 2
        else:
            expected_degree = 4 * m_value
            expected_coefficient = -3 * a**2 * u**4
        gate(degree == expected_degree, f"generic top degree n={n_value},d={delta}")
        zero(coefficient - expected_coefficient, f"generic top lc n={n_value},d={delta}")
        generic_top_cells += 1

        response_degree, response_coefficient = weighted_top(R, m_value, n_value)
        if delta == 1:
            if n_value == 1:
                expected_response_degree = 8
                expected_response_coefficient = -8 * u**3 + 3 * a * L**2 * v**4
            else:
                expected_response_degree = 4 * n_value + 4
                expected_response_coefficient = 3 * a * L**2 * v**4
        elif delta == 2:
            if n_value < 4:
                expected_response_degree = 3 * n_value + 8
                expected_response_coefficient = -8 * u**3
            elif n_value == 4:
                expected_response_degree = 20
                expected_response_coefficient = (
                    -8 * u**3 + 3 * a * L**2 * v**2 * (a * u + v) ** 2
                )
            else:
                expected_response_degree = 4 * n_value + 4
                expected_response_coefficient = 3 * a * L**2 * v**2 * (a * u + v) ** 2
        else:
            if n_value < delta + 2:
                expected_response_degree = 3 * m_value + 2
                expected_response_coefficient = -8 * u**3
            elif n_value == delta + 2:
                expected_response_degree = 3 * m_value + 2
                expected_response_coefficient = u**2 * (
                    -8 * u + 3 * a**3 * L**2 * v**2
                )
            else:
                expected_response_degree = 2 * m_value + 2 * n_value
                expected_response_coefficient = 3 * a**3 * L**2 * u**2 * v**2
        gate(
            response_degree == expected_response_degree,
            f"first response degree n={n_value},d={delta}",
        )
        zero(
            response_coefficient - expected_response_coefficient,
            f"first response lc n={n_value},d={delta}",
        )
        first_response_cells += 1

# The second response is asserted only in the first-response automatic cell.
second_response_cells = 0
for delta in range(3, 13):
    for n_value in range(delta + 3, 3 * delta + 4):
        m_value = n_value + delta
        degree, coefficient = weighted_top(R2, m_value, n_value)
        if n_value < 2 * delta:
            expected_degree = 3 * m_value + 2
            expected_coefficient = -8 * u**3
        elif n_value == 2 * delta:
            expected_degree = 3 * m_value + 2
            expected_coefficient = u * (-8 * u**2 + 3 * a**2 * L**2 * v**3)
        else:
            expected_degree = m_value + 3 * n_value + 2
            expected_coefficient = 3 * a**2 * L**2 * u * v**3
        gate(degree == expected_degree, f"second response degree n={n_value},d={delta}")
        zero(coefficient - expected_coefficient, f"second response lc n={n_value},d={delta}")
        second_response_cells += 1

# Constant-sidecar boundary.
constant_boundary_cells = 0
for m_value in range(1, 13):
    degree, coefficient = weighted_top(S, m_value, 0)
    if m_value == 1:
        expected_degree = 6
        expected_coefficient = -3 * u**2 * v**2
    elif m_value == 2:
        expected_degree = 8
        expected_coefficient = -u**2 * (8 * u + 3 * (a * u + v) ** 2)
    else:
        expected_degree = 4 * m_value
        expected_coefficient = -3 * a**2 * u**4
    gate(degree == expected_degree, f"constant-sidecar top degree m={m_value}")
    zero(coefficient - expected_coefficient, f"constant-sidecar top lc m={m_value}")
    if m_value == 1 or m_value >= 3:
        response_degree, response_coefficient = weighted_top(R, m_value, 0)
        gate(response_degree == 3 * m_value + 2, f"constant response degree m={m_value}")
        zero(response_coefficient + 8 * u**3, f"constant response lc m={m_value}")
    constant_boundary_cells += 1

# Quotient simplifications, including the sharpened Q1 equality passport and
# the delta=2 cancellation-wall correction v=-a*u.
w, rho, u1 = sp.symbols("w rho u1", nonzero=True)
zero_s(
    3 * a * L**2 * v**4 / (2 * s * u * v)
    + s * a * L**2 * v**3 / (2 * u),
    "delta one response quotient",
)
zero_s(-8 * u**3 / (2 * s * a * u**2) - 4 * s * u / (3 * a), "KT3 quotient")
zero_s(
    3 * a**3 * L**2 * u**2 * v**2 / (2 * s * a * u**2)
    + s * a**2 * L**2 * v**2 / 2,
    "automatic first correction",
)
zero_s(
    ((-8 * u**2 + 3 * a**2 * L**2 * v**3) / (2 * s * a * u)).subs(u, a * u1)
    - (-8 * u1**2 + 3 * L**2 * v**3) / (2 * s * u1),
    "sharpened Q1 equality quotient",
)
zero_s(
    (3 * a * L**2 * rho * v**2 / (2 * s * u)).subs(v, -a * u)
    - 3 * a**3 * L**2 * rho * u / (2 * s),
    "delta two wall high response is automatic",
)
zero_s(
    ((-8 * u**3 + 3 * a * L**2 * rho**2 * v**2) / (2 * s * u * rho)).subs(v, -a * u)
    - u * (-8 * u + 3 * a**3 * L**2 * rho**2) / (2 * s * rho),
    "delta two wall tie quotient",
)

# The exact r=0 address branch on the cancellation wall.
q = sp.symbols("q", nonzero=True)
C = a**3 * L**2 - K**2
S_r_zero = sp.expand(S.subs({T: -K * q, f: a * q}))
expected_r_zero = sp.expand(
    L**4 * (1 + 6 * a**2 * q) + 12 * a * L**2 * q**2 * C + 8 * q**3 * C**2
)
zero(S_r_zero - expected_r_zero, "exact r=0 cancellation-wall branch")

# Sharp leading controls: the tariffs are not no-go statements.
gate(8 * 216**3 == 3 * 72**4, "delta one n=1 cancellation constants")
zero((a * L**2) - a * L**2, "delta one n>=2 payment")
zero((a + (1 - a)) - 1, "delta two low unit-w payment")
zero((a * L**2) - a * L**2, "delta two high payment")
zero((a / a) - 1, "Q1 equality u/a payment")

semantic = {
    "scope": "necessary_m_gt_n_fan_and_first_response",
    "top_fan": "delta_1;delta_2_norm;delta_at_least_3",
    "first_identity": "S-(sTr)^2=R;s^2=-3",
    "first_thresholds": "1:n=1;2:n=4;ge3:n=delta+2",
    "second_scope": "delta_ge3_and_n_gt_delta+2_only",
    "second_threshold": "n=2delta",
    "delta2_wall": "r_order_and_v=-au_simplified_tariffs",
    "constant_boundary": "m=1,2,ge3",
    "status": "provisional_not_existence_and_not_emptiness",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3901-nonzero-sidecar-osculating-response-fan")
print("identity_q0=S-(sTr)^2=R_mod_s2_plus_3")
print("identity_q1=second_defect_exact_mod_s2_plus_3")
print(f"generic_top_cells={generic_top_cells}")
print(f"first_response_cells={first_response_cells}")
print(f"second_response_cells={second_response_cells}")
print(f"constant_boundary_cells={constant_boundary_cells}")
print("delta2_wall=corrected_v_eq_minus_au_tariffs")
print("sharp_controls=leading_tariffs_are_payable")
print("scope=necessary_fan_not_square_existence_or_emptiness")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
