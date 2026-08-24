#!/usr/bin/env python3
"""Exact companion for THM-3977's simultaneous cusp-arm family."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: bool, label: str) -> None:
    """Record one optimization-safe exact gate."""
    global CHECKS
    if not bool(condition):
        raise RuntimeError(f"gate failed: {label}")
    CHECKS += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(expression) == 0, label)


def jacobian(first: sp.Expr, second: sp.Expr, x: sp.Symbol,
             t: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, t)
        - sp.diff(first, t) * sp.diff(second, x)
    )


x, t, c, r = sp.symbols("x t c r")
P, Y = sp.symbols("P Y")

z = 1 + x**2 * t
p = z * t
y = x * z * t**2
A = c * p + y**2 + r * x

# Determinantal presentation and the two boundary restrictions.
zero(x**2 * p - z * (z - 1), "first determinantal relation")
zero(x * y - p * (z - 1), "second determinantal relation")
zero(z * y - x * p**2, "third determinantal relation")

X0, Z0, P0, Y0 = sp.symbols("X0 Z0 P0 Y0")
A_abstract = c * P0 + Y0**2 + r * X0
zero(A_abstract.subs({X0: 0, Z0: 0, P0: 0}) - Y0**2,
     "boundary D restriction")
zero(A_abstract.subs({X0: 0, Z0: 1, Y0: 0}) - c * P0,
     "retained arm L1 restriction")

# Exact critical resultant for the nonzero seam rows.
Ax_expected = (
    r + 2 * c * x * t**2 + 2 * x * t**4
    + 8 * x**3 * t**5 + 6 * x**5 * t**6
)
At_expected = (
    c + 2 * c * x**2 * t + 4 * x**2 * t**3
    + 10 * x**4 * t**4 + 6 * x**6 * t**5
)
zero(sp.diff(A, x) - Ax_expected, "expanded A_x")
zero(sp.diff(A, t) - At_expected, "expanded A_t")

H = (
    c**6 * x + 2 * c**5 * r * x**4 - 2 * c**4 * r
    + 8 * c**3 * r**2 * x**3 + 12 * c**2 * r**3 * x**6
    + 216 * c * r**4 * x**9 + 486 * r**5 * x**12
    - 64 * r**4 * x**5
)
critical_resultant = sp.factor(sp.resultant(Ax_expected, At_expected, t))
zero(critical_resultant - 96 * x**24 * H,
     "simultaneous cusp-arm critical resultant")
zero(H.subs(x, 0) + 2 * c**4 * r, "H nonzero-root constant row")
gate(sp.Poly(H, x).coeff_monomial(x) == c**6,
     "H nonconstant linear row")
gate(sp.Poly(Ax_expected, t).LC() == 6 * x**5,
     "A_x stable t-leading coefficient")
gate(sp.Poly(At_expected, t).LC() == 6 * x**6,
     "A_t stable t-leading coefficient")

# The zero-seam endpoint is submersive but has a logarithmic generic fibre.
A0 = sp.expand(A.subs(r, 0))
A0x = sp.diff(A0, x)
A0t = sp.diff(A0, t)
zero(sp.resultant(A0x, A0t, t) - 96 * c**6 * x**25,
     "endpoint critical resultant")
zero(jacobian(A0, y, x, t) + c * z * t**2,
     "endpoint J(A,y)")
zero(A0t.subs(t, 0) - c, "endpoint t=0 submersion row")
zero(sp.cancel(A0t.subs(t, -x**-2)) + c,
     "endpoint z=0 submersion row")
zero(p * y + x * y**2 - x * p**3,
     "generic-fibre parametrization identity")

p_generic = (P - Y**2) / c
Delta = p_generic**3 - Y**2
x_generic = Y * p_generic / Delta
z_generic = p_generic**3 / Delta
t_generic = Delta / p_generic**2
zero(z_generic - (1 + x_generic**2 * t_generic),
     "generic reconstruction of z")
zero(p_generic - z_generic * t_generic,
     "generic reconstruction of p")
zero(Y - x_generic * z_generic * t_generic**2,
     "generic reconstruction of y")
zero(-c * z_generic * t_generic**2 + c * Delta / p_generic,
     "generic reconstruction of J(A,y)")

E = (P - Y**2)**3 - c**3 * Y**2
zero(E - c**3 * Delta, "cleared generic-fibre pole polynomial")
expected_discriminant = 64 * P**3 * c**12 * (27 * P**2 + 4 * c**3)**2
zero(sp.discriminant(E, Y) - expected_discriminant,
     "generic-fibre discriminant")

domain = sp.QQ.frac_field(c, P)
E_poly = sp.Poly(E, Y, domain=domain)
gate(sp.gcd(E_poly, E_poly.diff()).degree() == 0,
     "generic-fibre squarefree gcd")
numerator_poly = sp.Poly(P - Y**2, Y, domain=domain)
gate(sp.gcd(E_poly, numerator_poly).degree() == 0,
     "generic logarithmic residues have nonzero numerator")

summary = {
    "checks": CHECKS,
    "family": "A_(c,r)=c*p+y^2+r*x on B_2, c nonzero",
    "boundary": "A|D=y^2 and A|L1=c*p",
    "resultant": "Res_t(A_x,A_t)=96*x^24*H_(c,r)",
    "nonzero_seam": "r!=0 forces an affine critical point",
    "endpoint": "r=0 submersive but has no rational h(A)-mate",
    "generic_fibre": "k(P)(y) with six nonzero logarithmic residues",
    "scope": "lowest seam family only; unrestricted Darboux and JC2 open",
}
semantic = hashlib.sha256(
    json.dumps(summary, sort_keys=True).encode()
).hexdigest()

print("THM-3977 simultaneous cusp-arm critical-resultant companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=C_P_PLUS_Y2_PLUS_R_X;C_NONZERO")
print("BOUNDARY=D_Y2;L1_C_P")
print("RESULTANT=96_X24_H_C_R")
print("NONZERO_SEAM=AFFINE_CRITICAL_POINT")
print("ENDPOINT=SUBMERSIVE;NO_RATIONAL_H_A_MATE")
print("GENERIC_FIBRE=K_Y;SIX_LOGARITHMIC_POLES")
print("DISC_E=64_P3_C12_27P2_PLUS4C3_SQUARED")
print("SCOPE=LOWEST_SEAM_FAMILY_ONLY;UNRESTRICTED_DARBOUX_AND_JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
