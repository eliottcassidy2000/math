#!/usr/bin/env python3
"""Exact companion for THM-3831's intrinsic spectral-fibre atlas."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


def numerator(expr: object) -> sp.Expr:
    return sp.expand(sp.cancel(expr).as_numer_denom()[0])  # type: ignore[arg-type]


a, k, z = sp.symbols("a k z")
Q = 7 * a**2 + 3
P = 3 * a**3 + 7 * a**2 + 1
b = (a + 1) * (2 * a + 1) * (3 * a - 1)

check(sp.gcd(sp.Poly(Q, a), sp.Poly(P, a)).degree() == 0,
      "quadratic and cubic spectral slopes are disjoint")
check(sp.gcd(sp.Poly(P, a), sp.Poly(a * Q, a)).degree() == 0,
      "a and Q are units on every cubic spectral slope")
check(sp.gcd(sp.Poly(Q, a), sp.Poly(a * (9 * a + 14), a)).degree() == 0,
      "quadratic-slope Laurent denominators are units")
check(sp.gcd(sp.Poly(Q, a), sp.Poly(b, a)).degree() == 0,
      "quadratic slopes avoid every B3 zero")
check(sp.discriminant(P, a) != 0 and sp.discriminant(Q, a) != 0,
      "all five spectral slopes are distinct")
check(sp.rem(b + Q, P, a) == 0,
      "on a cubic slope the B3 coefficient is minus Q")

# Use the actual THM-3811 root chart, not merely the cleared lift laws.  Its
# J-localization is the saturated affine-modification chart and contains every
# nonzero spectral fibre, including its D=0 points.
AA, CC, u = sp.symbols("AA CC u")
G_root = u**3 + 7 * u + 3
F_root = AA * G_root - CC * (CC + u**2)
J_root = AA * (3 * u**2 + 7) - 2 * CC * u
same(a**3 * G_root.subs(u, 1 / a), P,
     "marked-root chart specializes to the cubic slope polynomial")
same(a**3 * F_root.subs(u, 1 / a),
     AA * P - a * CC * (a**2 * CC + 1),
     "exact specialized root-chart equation")
same(a**2 * J_root.subs(u, 1 / a), AA * Q - 2 * a * CC,
     "exact specialized root-chart Jacobian")
check(sp.gcd(sp.Poly(a, a), sp.Poly(P * Q, a)).degree() == 0,
      "all spectral slopes avoid zero, so the missing P1 divisor is disjoint")
same((1 + a**2 * CC) - a**2 * CC, 1,
     "cubic root-chart component ideals are comaximal")


def reconstruction_relations(
    A: sp.Expr,
    C: sp.Expr,
    omega: sp.Expr,
    theta: sp.Expr,
    D: sp.Expr,
    h: sp.Expr,
    kk: sp.Expr,
    m: sp.Expr,
) -> list[sp.Expr]:
    """SL2/lift laws plus the original cubic algebra and different."""
    return [
        C * kk - m * h - 1,
        D * (7 * h**2 + 3 * kk**2) - 1 - 2 * C * kk,
        h * D * (9 * h + 14 * kk) - kk * m - 3 * h * C**2,
        A - h * D,
        omega - kk * D,
        theta - (m - 14 * A) / 3,
        omega**2 + 7 * A**2 - C * omega + A * theta,
        omega * theta - 3 * A**2 + A * C**2,
        theta**2 - 3 * A * C + C**3
        - (C**2 - 3 * A) * omega + 7 * A * theta,
        D - (C * omega - 3 * A * theta - 14 * A**2),
    ]


# Cubic spectral slopes.  Put z=Ck.  The compatibility numerator factors
# exactly as z(k+a^2 z) after imposing P(a)=0.
N0 = sp.expand(a**2 * (9 * a + 14) + Q)
N1 = sp.expand(2 * a**2 * (9 * a + 14) - Q)
same(N0, 3 * P, "constant compatibility coefficient")
same(N1, 3 * b, "linear compatibility coefficient")
check(sp.rem(N1 + 3 * Q, P, a) == 0,
      "cubic compatibility is z(k+a^2 z)")

h_c = a * k
C_c = z / k
m_c = (z - 1) / (a * k)
D_c = (1 + 2 * z) / (Q * k**2)
A_c = h_c * D_c
omega_c = k * D_c
theta_c = (m_c - 14 * A_c) / 3
J_c = 1 / (a * k)

G_cubic = sp.groebner(
    [P, z * (k + a**2 * z)], z, k, a, order="grevlex", domain=sp.QQ
)
for index, relation in enumerate(
    reconstruction_relations(A_c, C_c, omega_c, theta_c, D_c, h_c, k, m_c)
):
    check(G_cubic.reduce(numerator(relation))[1] == 0,
          f"cubic universal reconstruction relation {index}")

for expression, label in (
    (F_root.subs({u: 1 / a, AA: A_c, CC: C_c}),
     "cubic family lies in the actual root chart"),
    (J_root.subs({u: 1 / a, AA: A_c, CC: C_c}) - J_c,
     "cubic family has the correct localized root-chart J"),
    (D_c - A_c * J_c,
     "cubic family reconstructs D=AJ on the saturated chart"),
):
    check(G_cubic.reduce(numerator(expression))[1] == 0, label)

same((k + a**2 * z) - a**2 * z, k,
     "cubic component ideals are comaximal after k inversion")

# The normalized discriminant root is regular on every spectral fibre because
# k is a unit there.  It labels the two cubic components with opposite signs.
q_actual_c = 7 * h_c**2 + 3 * k**2
B0_c = (
    k**5 - 7 * h_c**2 * k**3 - 3 * h_c**2 * k**2
    - 6 * h_c**3 * k**2 - 7 * h_c**4
)
W_c = sp.cancel((h_c**2 * q_actual_c**2 * D_c + B0_c) / k)
B3_c = (h_c + k) * (2 * h_c + k) * (3 * h_c - k)
G_minus = sp.groebner([P, z], z, k, a, order="grevlex", domain=sp.QQ)
G_plus = sp.groebner(
    [P, k + a**2 * z], z, k, a, order="grevlex", domain=sp.QQ
)
check(G_minus.reduce(numerator(W_c + k * B3_c))[1] == 0,
      "C=0 cubic component has the minus square-root sign")
check(G_plus.reduce(numerator(W_c - k * B3_c))[1] == 0,
      "C=-1/a^2 cubic component has the plus square-root sign")

# Quadratic spectral slopes.  Q(a)=0 forces z=Ck=-1/2 and leaves one Laurent
# component.  All original laws are checked, not just the three lift laws.
h_q = a * k
C_q = -1 / (2 * k)
m_q = -3 / (2 * a * k)
D_q = (14 * k + 3) / (4 * (9 * a + 14) * k**3)
A_q = h_q * D_q
omega_q = k * D_q
theta_q = (m_q - 14 * A_q) / 3
J_q = 1 / (a * k)
G_quadratic = sp.groebner([Q], a, order="grevlex", domain=sp.QQ)


def zero_mod_quadratic(expr: object, label: str) -> None:
    poly = sp.Poly(numerator(expr), k)
    check(all(G_quadratic.reduce(coefficient)[1] == 0
              for coefficient in poly.all_coeffs()), label)


for index, relation in enumerate(
    reconstruction_relations(A_q, C_q, omega_q, theta_q, D_q, h_q, k, m_q)
):
    zero_mod_quadratic(relation, f"quadratic reconstruction relation {index}")

zero_mod_quadratic(
    F_root.subs({u: 1 / a, AA: A_q, CC: C_q}),
    "quadratic family lies in the actual root chart",
)
zero_mod_quadratic(
    J_root.subs({u: 1 / a, AA: A_q, CC: C_q}) - J_q,
    "quadratic family has the correct localized root-chart J",
)
zero_mod_quadratic(D_q - A_q * J_q,
                   "quadratic family reconstructs D=AJ")
zero_mod_quadratic(A_q.subs(k, a**2 / 2),
                   "quadratic Laurent chart retains its D-zero point")
zero_mod_quadratic(D_q.subs(k, a**2 / 2),
                   "quadratic exceptional point really has D zero")

q_actual_q = 7 * h_q**2 + 3 * k**2
B0_q = (
    k**5 - 7 * h_q**2 * k**3 - 3 * h_q**2 * k**2
    - 6 * h_q**3 * k**2 - 7 * h_q**4
)
W_q = sp.cancel((h_q**2 * q_actual_q**2 * D_q + B0_q) / k)
B3_q = (h_q + k) * (2 * h_q + k) * (3 * h_q - k)
zero_mod_quadratic(W_q + k * B3_q,
                   "quadratic component has only the minus sign")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "quadratic_slopes": "Q(a)=7a^2+3=0: one Gm, Ck=-1/2, W=-kB3",
    "cubic_slopes": "P(a)=3a^3+7a^2+1=0: Gm_minus disjoint_union Gm_plus",
    "cubic_equation": "z=Ck; z(k+a^2z)=0",
    "components": "minus C=0; plus C=-1/a^2; exact saturated THM-3811 root-chart quotients including D=0",
    "atlas": "irreducible h forces one cubic slope and nonempty pullback of both components",
    "scope": "necessary atlas passport only; no atlas and no Keller pair constructed",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit")
print("slopes=two_quadratic_plus_three_cubic")
print("quadratic=one_Gm;Ck=-1/2;sign=minus")
print("cubic=Gm_minus_disjoint_union_Gm_plus;equation=z*(k+a^2*z)")
print("labels=minus:C=0;plus:C=-1/a^2")
print("reconstruction=actual_root_chart_saturation_plus_all_original_laws;D_zero_retained")
print("atlas=irreducible_h_forces_cubic_slope_and_hits_both_components")
print("scope=no_atlas_or_Keller_pair_constructed")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
