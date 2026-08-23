#!/usr/bin/env python3
"""Exact companion for THM-3887's binary-cubic pencil contact gate."""

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
    gate(sp.expand(expression) == 0, message)


def disc(a: sp.Expr, b: sp.Expr, c: sp.Expr, d: sp.Expr) -> sp.Expr:
    return sp.expand(
        b**2 * c**2
        - 4 * a * c**3
        - 4 * b**3 * d
        - 27 * a**2 * d**2
        + 18 * a * b * c * d
    )


t = sp.symbols("t")
alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")

# Double-plus-simple root: F=X^2Y.
double_contact = sp.Poly(
    disc(t * alpha, 1 + t * beta, t * gamma, t * delta), t
)
zero(double_contact.coeff_monomial(t) + 4 * delta, "double-root t1")
zero(
    double_contact.coeff_monomial(t**2) - (gamma**2 - 12 * beta * delta),
    "double-root t2",
)
double_delta_zero = sp.Poly(
    double_contact.as_expr().subs({delta: 0, gamma: 1}), t
)
gate(
    double_delta_zero.coeff_monomial(1) == 0
    and double_delta_zero.coeff_monomial(t) == 0
    and double_delta_zero.coeff_monomial(t**2) == 1,
    "double-root delta-zero contact two",
)
zero(
    double_contact.as_expr().subs({delta: 0, gamma: 0}),
    "double-root contained-line boundary",
)

# Triple root: F=X^3.
triple_contact = sp.Poly(
    disc(1 + t * alpha, t * beta, t * gamma, t * delta), t
)
zero(triple_contact.coeff_monomial(t), "triple-root missing t1")
zero(triple_contact.coeff_monomial(t**2) + 27 * delta**2, "triple-root t2")
zero(
    triple_contact.coeff_monomial(t**3)
    - (-54 * alpha * delta**2 + 18 * beta * gamma * delta - 4 * gamma**3),
    "triple-root t3",
)
triple_delta_zero = sp.Poly(
    triple_contact.as_expr().subs({delta: 0, gamma: 1}), t
)
gate(
    triple_delta_zero.coeff_monomial(1) == 0
    and triple_delta_zero.coeff_monomial(t) == 0
    and triple_delta_zero.coeff_monomial(t**2) == 0
    and triple_delta_zero.coeff_monomial(t**3) == -4,
    "triple-root delta-zero contact three",
)
zero(
    triple_contact.as_expr().subs({delta: 0, gamma: 0}),
    "triple-root contained-line boundary",
)

# Rank-one first jet.  Scale every base variable by t: ell*H has order one,
# while Q has order two.  The degree-five coefficient must carry ell^3.
ell = sp.symbols("ell")
h0, h1, h2, h3 = sp.symbols("h0 h1 h2 h3")
q0, q1, q2, q3 = sp.symbols("q0 q1 q2 q3")
rank_one = sp.Poly(
    disc(
        t * ell * h0 + t**2 * q0,
        t * ell * h1 + t**2 * q1,
        t * ell * h2 + t**2 * q2,
        t * ell * h3 + t**2 * q3,
    ),
    t,
)
H_disc = disc(h0, h1, h2, h3)
zero(rank_one.coeff_monomial(t**4) - ell**4 * H_disc, "quartic rank-one jet")
degree_five = sp.factor(rank_one.coeff_monomial(t**5))
gate(sp.rem(degree_five, ell**3, ell) == 0, "quintic jet divisible by ell^3")
polar = sp.expand(degree_five / ell**3)
gate(not polar.has(ell), "quintic polar independent of ell")

# Direct differential check of the polar term.
s = sp.symbols("s")
directional = sp.diff(
    disc(h0 + s * q0, h1 + s * q1, h2 + s * q2, h3 + s * q3), s
).subs(s, 0)
zero(polar - directional, "quintic jet is the discriminant differential")

# Hostile and positive controls.
zero(disc(1, 0, 0, 0), "triple-root singular control")
zero(disc(1, 0, 0, 1) + 27, "nonsingular binary cubic control")
zero(disc(1, 3 * t, 3 * t**2, t**3), "triple-root family control")

semantic = {
    "object": "binary cubic discriminant at a common coefficient zero",
    "pencil_contact": "finite contact at a singular member is at most three",
    "fourth_power": "forces rank-one first jet",
    "quintic_jet": "divisible by the cube of the unique tangent",
    "consequence": "reduced irreducible locally unibranch discriminant has degree at least six",
    "scope": "multibranch quintics and unit-ideal nonmonogenic index forms remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("theorem=THM-3887-binary-cubic-common-zero-quintic-unibranch-obstruction")
print("finite_pencil_discriminant_contact_at_most=3")
print("fourth_power_tangent_cone=rank_one_coefficient_jet")
print("quintic_jet_divisor=ell^3")
print("common_zero_reduced_irreducible_unibranch_degree_floor=6")
print("multibranch_quintic=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
