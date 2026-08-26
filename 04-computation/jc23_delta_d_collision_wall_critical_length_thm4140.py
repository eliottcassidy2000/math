#!/usr/bin/env python3
"""Primary exact audit for THM-4140's Delta_D collision wall.

The calculation works in the normalized (X,T) chart.  It parameterizes the
live repeated quadratic edge, constructs the complete source polynomial,
factors the critical resultant, and audits the T=0 and T=-1/6 strata.

No Python ``assert`` is used, so optimized mode retains every truth gate.
"""

from __future__ import annotations

import hashlib

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


X, T, r = sp.symbols("X T r")

# On the normalized delta-only wall, kappa=2848/45-(7/6)delta and
# phi^2=4*kappa*delta.  Writing the repeated edge as
# delta*(Z+r)^2 gives the following rational parameterization.  The live
# gate kappa*delta!=0 is exactly r!=0; 6*r^2+7 cannot vanish on this chart.
delta = sp.Rational(5696, 15) / (6 * r**2 + 7)
kappa = sp.cancel(delta * r**2)
phi = sp.cancel(2 * delta * r)

require(
    sp.factor(kappa + sp.Rational(7, 6) * delta - sp.Rational(2848, 45)) == 0,
    "forced delta/kappa row changed",
)
require(sp.factor(phi**2 - 4 * kappa * delta) == 0, "Delta_D wall changed")
require(sp.factor(kappa / delta - r**2) == 0, "repeated-root parameter changed")

P = T + X**2 * T**2
Y = X * T * P
G = sp.expand(
    -sp.Rational(1, 2) * X**2 * T
    - 3 * P
    + sp.Rational(8, 3) * P**2
    - sp.Rational(1376, 135) * P**3
    + kappa * Y**2
    + phi * P**2 * Y
    + delta * P**4
)

GX = sp.diff(G, X)
GT = sp.diff(G, T)
f = sp.cancel(GX / T)
h = GT
require(not sp.denom(f).has(X, T), "G_X/T is not polynomial over Q(r)")
require(sp.factor(GX - T * f) == 0, "lost the exact T factor")
require(sp.Poly(f, X).degree() == 7, "wrong generic X-degree for f")
require(sp.Poly(h, X).degree() == 8, "wrong generic X-degree for h")
expected_x_lead = 8 * delta * T**7
require(sp.factor(sp.Poly(f, X).LC() - expected_x_lead) == 0, "wrong f leading row")
require(sp.factor(sp.Poly(h, X).LC() - expected_x_lead) == 0, "wrong h leading row")


def remainder_at(poly: sp.Expr, t_value: sp.Rational, modulus: sp.Expr) -> sp.Expr:
    specialized = sp.cancel(poly.subs(T, t_value))
    require(
        not sp.denom(specialized).has(X),
        "specialization is not polynomial over Q(r)",
    )
    return sp.factor(sp.rem(sp.Poly(specialized, X), sp.Poly(modulus, X)).as_expr())


hessian_det = sp.factor(sp.det(sp.hessian(G, (X, T))))

# T=0 belongs to the actual critical ideal (T*f,h), despite being a
# degree-drop artefact for Res_X(f,h).  It contributes two Morse points.
require(sp.factor(f.subs(T, 0)) == -X, "wrong f at T=0")
require(
    sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0,
    "wrong h at T=0",
)
require(remainder_at(G, sp.Rational(0), X**2 + 6) == 0, "wrong T=0 value")
require(
    remainder_at(hessian_det, sp.Rational(0), X**2 + 6) == 6,
    "T=0 pair is not Morse",
)

# P=0 at T=-1/6, X^2=6, so every maximum-eight response term vanishes.
universal_modulus = X**2 - 6
require(remainder_at(f, -sp.Rational(1, 6), universal_modulus) == 0, "universal f failed")
require(remainder_at(h, -sp.Rational(1, 6), universal_modulus) == 0, "universal h failed")
require(
    remainder_at(G - sp.Rational(1, 2), -sp.Rational(1, 6), universal_modulus) == 0,
    "wrong universal critical value",
)
require(
    remainder_at(hessian_det, -sp.Rational(1, 6), universal_modulus) == -6,
    "universal pair is not Morse",
)

# Exact critical resultant.  T^42 is the same generic-X-degree drop seen in
# THM-4130.  The universal pair gives (6T+1)^2; the strict residual now has
# degree fourteen throughout the live r!=0 wall.
resultant_xt = sp.factor(sp.resultant(f, h, X))
Q14_expr = sp.cancel(resultant_xt / (T**42 * (6 * T + 1) ** 2))
Q14 = sp.Poly(Q14_expr, T)
require(Q14.degree() == 14, "strict residual is not degree fourteen")

expected_q14_lead = (
    sp.Integer(48947063421231612219454924509816611957374976)
    * r**10
    / (sp.Integer(8649755859375) * (6 * r**2 + 7) ** 9)
)
expected_q14_constant = -(
    sp.Integer(139887797298879228245180416)
    / (sp.Integer(3796875) * (6 * r**2 + 7) ** 6)
)
require(sp.factor(Q14.LC() - expected_q14_lead) == 0, "Q14 leading row changed")
require(sp.factor(Q14.TC() - expected_q14_constant) == 0, "Q14 constant row changed")
require(
    sp.factor(resultant_xt - T**42 * (6 * T + 1) ** 2 * Q14.as_expr()) == 0,
    "resultant reconstruction failed",
)

monic_coefficients = [sp.factor(coefficient / Q14.LC()) for coefficient in Q14.all_coeffs()]
coefficient_payload = "\n".join(map(str, monic_coefficients)) + "\n"
coefficient_sha256 = hashlib.sha256(coefficient_payload.encode()).hexdigest()
require(
    coefficient_sha256 == "1ef62ff920010fcc5eb64c2a8e09e90a03b6c1b4b0b49beb17f4c63ec509ec74",
    "Q14 coefficient ledger changed",
)

# A fully rational live point supplies a squarefree positive control.
Q14_r1 = sp.Poly(Q14.as_expr().subs(r, 1), T)
clear_denominator = sp.lcm([sp.denom(value) for value in Q14_r1.all_coeffs()])
Q14_mod_101 = sp.Poly(sp.expand(clear_denominator * Q14_r1.as_expr()), T, modulus=101)
require(sp.gcd(Q14_mod_101, Q14_mod_101.diff()).degree() == 0, "r=1 control not squarefree mod 101")
require(kappa.subs(r, 0) == 0, "r=0 did not expose the excluded kappa=0 wall")

semantic_lines = (
    "scope=delta_only_exact_M8_DeltaD_wall",
    "parameter=delta=5696/(15(6r^2+7));kappa=delta*r^2;phi=2delta*r;r!=0;6r^2+7!=0",
    "critical_projection_residual_degree=14",
    "universal_t=-1/6:length=2",
    "omitted_t=0:length=2",
    "affine_critical_length=18",
    "critical_values=target_nodes_only_under_Keller",
    "verdict=FINITE_CRITICAL_LEDGER_ONLY;DeltaD_wall=OPEN;JC2=OPEN",
)
semantic_sha256 = hashlib.sha256(("\n".join(semantic_lines) + "\n").encode()).hexdigest()

print("THM4140_DELTA_D_CRITICAL_PRIMARY")
print("wall=delta=5696/(15(6r^2+7));kappa=delta*r^2;phi=2delta*r;r!=0;6r^2+7!=0")
print("critical_resultant=T^42*(6T+1)^2*Q14;degree_Q14=14")
print("Q14_monic_coefficient_sha256=" + coefficient_sha256)
print("universal=T=-1/6,X^2=6,length=2,value=1/2")
print("omitted=T=0,X^2=-6,length=2,value=0")
print("positive_control=r=1;Q14_squarefree_mod_101=YES")
print("hostile_gate=r=0=>kappa=0;EXCLUDED_FROM_LIVE_WALL")
print("semantic_sha256=" + semantic_sha256)
print("verdict=affine_critical_length=18;DeltaD_wall=OPEN;JC2=OPEN")
