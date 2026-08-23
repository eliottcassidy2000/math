#!/usr/bin/env python3
"""Exact companion for THM-3792's pure first-normal nodal obstruction."""

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


r, z, e = sp.symbols("r z e")
c = sp.symbols("c", nonzero=True)
variables = (r, z, e)
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 * c**3 + 6 * r * e],
        [-9 * z**2, -3 * c**3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


# Generic Hamiltonian identities.
f = sp.Function("f")(e)
A = e**2 + z * f
H = f - 3 * e * sp.diff(f, e)
K = c**3 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(critical[0] - (-3 * f * r**2 - 9 * (2 * e + z * sp.diff(f, e)) * z**2),
     "generic critical r component")
zero(critical[1] + 3 * (2 * e + z * sp.diff(f, e)) * K,
     "generic critical z component")
zero(critical[2] - 3 * f * K, "generic critical e component")


# Every polynomial first-normal Bezout completion has the required nonzero
# constant, and the displayed h-parametrization is the full syzygy family.
t = sp.symbols("t")
h = sp.Function("h")(t)
f_jet = -1 / (3 * c**3) + 2 * t * h
g_jet = -t / (2 * c**3) + (3 * t**2 - 1) * h
bezout = sp.factor(3 * c**3 * (f_jet * (3 * t**2 - 1) - 2 * t * g_jet))
gate(bezout == 1, "all-polynomial Bezout family")
zero(f_jet.subs(t, 0) + 1 / (3 * c**3), "forced Bezout constant")


# Exact reconstruction at a root of G_f.  Treat F=f(eta) and HH=H_f(eta)
# as independent symbols modulo c^6 HH^3+864 eta^7.
eta, F, HH = sp.symbols("eta F HH", nonzero=True)
relation = c**6 * HH**3 + 864 * eta**7
r0 = -c**3 / (2 * eta)
z0 = 6 * eta**2 / HH
fp0 = (F - HH) / (3 * eta)


def reduce_relation(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.factor(sp.rem(numerator, relation, HH))


surface = r**2 * e - z**3 + c**3 * r
surface_at = surface.subs({r: r0, z: z0, e: eta})
gate(reduce_relation(surface_at) == 0, "reconstructed point lies on surface")
zero(K.subs({r: r0, e: eta}), "reconstructed K sheet")
z_cube = z0**3 + c**6 / (4 * eta)
gate(reduce_relation(z_cube) == 0, "reconstructed z cube")
remaining = c**6 * HH + 24 * eta**3 * z0**2
gate(reduce_relation(remaining) == 0, "remaining critical equation")

critical_r_reconstructed = -3 * F * r0**2 - 9 * (2 * eta + z0 * fp0) * z0**2
gate(reduce_relation(critical_r_reconstructed) == 0,
     "reconstructed critical r component")


# Degree and no-exception controls.  A generic monic degree-m polynomial is
# enough to verify the exact leading coefficient and the never-equal degrees;
# the proof handles arbitrary lower coefficients coefficientwise.
degree_rows: list[str] = []
for m in range(257):
    if m == 0:
        polynomial = sp.Integer(5)
    else:
        polynomial = e**m + 2 * e + 5 if m != 1 else 3 * e + 5
    transformed = sp.expand(polynomial - 3 * e * sp.diff(polynomial, e))
    gate(sp.degree(transformed, e) == m, f"H degree m={m}")
    expected_lc = (1 - 3 * m) * sp.LC(sp.Poly(polynomial, e))
    gate(sp.LC(sp.Poly(transformed, e)) == expected_lc,
         f"H leading coefficient m={m}")
    G = sp.expand(c**6 * transformed**3 + 864 * e**7)
    gate(sp.degree(G, e) == max(3 * m, 7), f"G degree m={m}")
    gate(sp.expand(G.subs(e, 0)) == c**6 * polynomial.subs(e, 0)**3,
         f"G constant m={m}")
    gate(3 * m != 7, f"degree noncollision m={m}")
    gate(1 - 3 * m != 0, f"H nonzero leading factor m={m}")
    degree_rows.append(f"{m}:{sp.degree(G,e)}:{expected_lc}")


# Canonical THM-3790 boundary: the e-equation and recovered z-equation agree
# exactly with the seven-point polynomial there.
f0 = -1 / (3 * c**3)
G0 = sp.factor(c**6 * f0**3 + 864 * e**7)
zero(G0 - (864 * e**7 - 1 / (27 * c**3)), "canonical G polynomial")
z_canonical = sp.factor(6 * e**2 / f0)
zero(z_canonical + 18 * c**3 * e**2, "canonical z recovery")
canonical_z_poly = 8 * z_canonical**7 + 9 * c**15
numerator = sp.together(canonical_z_poly).as_numer_denom()[0]
canonical_relation = sp.together(G0).as_numer_denom()[0]
gate(sp.rem(numerator, canonical_relation, e) == 0,
     "canonical seven-point boundary")


semantic = {
    "carrier": "A_f=e^2+z*f(e), f(0)!=0",
    "critical_polynomial": "G_f=c^6*(f-3e*f')^3+864e^7",
    "quantifiers": "algebraically_closed_characteristic_zero;c_nonzero;all_polynomial_f",
    "recovery": "r=-c^3/(2eta);z=6eta^2/(f(eta)-3eta*f'(eta))",
    "scope": "all pure first-normal nodal Bezout lifts fail; I^2 corrections open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
degree_blob = "\n".join(degree_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3792-pure-first-normal-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+c^3r;field=algebraically_closed_characteristic_zero;c!=0")
print("carrier=A_f=e^2+z*f(e);f(0)!=0")
print("critical_sheet=K=c^3+2re=0")
print("critical_polynomial=G_f=c^6*(f-3e*fprime)^3+864e^7")
print("recovery=e=eta;r=-c^3/(2eta);z=6eta^2/(f-3eta*fprime)")
print("degree_control=m[0,256];deg_G=max(3m,7);3m_ne_7")
print("canonical_control=f=-1/(3c^3);8z^7+9c^15=0")
print("closed=all_polynomial_first_normal_Bezout_parameters")
print("open=I^2_corrections;different_arm_immersion;higher_boundary_bidegree")
print(f"degree_sha256={hashlib.sha256(degree_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
