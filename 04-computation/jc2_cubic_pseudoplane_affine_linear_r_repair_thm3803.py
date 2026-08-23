#!/usr/bin/env python3
"""Exact companion for THM-3803's affine-linear r-repair obstruction."""

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


r, z, e, u = sp.symbols("r z e u")
b0, b1 = sp.symbols("b0 b1", nonzero=True)
variables = (r, z, e)
surface = r**2 * e - z**3 + r
poisson = sp.Matrix(
    [
        [0, 3 * r**2, 9 * z**2],
        [-3 * r**2, 0, 3 + 6 * r * e],
        [-9 * z**2, -3 - 6 * r * e, 0],
    ]
)


def bracket(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    dl = sp.Matrix([sp.diff(left, q) for q in variables])
    dr = sp.Matrix([sp.diff(right, q) for q in variables])
    return sp.expand((dl.T * poisson * dr)[0])


# Hamiltonian reduction for the affine-linear repair.
g = b0 + b1 * e
gp = b1
A = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
     "r Hamiltonian")
zero(critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
     "z Hamiltonian")
zero(critical[2] - (9 * g * z**2 - K_source),
     "e Hamiltonian")
zero(bracket(A, surface), "surface Casimir")

K = 1 + 2 * u
P = sp.expand(g * u**2 - K * (2 * e**3 + u * e * gp))
Q = sp.expand(e**2 * K**3 - 729 * g**3 * u**2 * (1 + u)**2)
zero(sp.LC(sp.Poly(P, u)) - (b0 - b1 * e), "P leading coefficient")
zero(sp.LC(sp.Poly(Q, u)) + 729 * g**3, "Q leading coefficient")


# Exact degree-eleven residual resultant.
resultant = sp.factor(sp.resultant(P, Q, u))
H = sp.factor(resultant / (e**4 * g**3))
coefficients = (
    b0 * (1 - 729 * b0 * b1**2),
    -2916 * b0**3 - b1,
    2916 * b0**2 * (729 * b0**3 + b1),
    2916 * b0 * b1 * (2187 * b0**3 - b1),
    8748 * b0**2 * (729 * b0 * b1**2 - 2),
    8748 * b0 * (972 * b0**3 + 243 * b0 * b1**3 + 2 * b1),
    25509168 * b0**3 * b1,
    11664 * b0 * (2187 * b0 * b1**2 - 2),
    11664 * (729 * b0**3 + 729 * b0 * b1**3 + 2 * b1),
    25509168 * b0**2 * b1,
    25509168 * b0 * b1**2,
    8503056 * b1**3,
)
H_table = sp.expand(sum(coefficient * e**j
                        for j, coefficient in enumerate(coefficients)))
zero(H - H_table, "residual coefficient table")
zero(resultant - e**4 * g**3 * H_table, "resultant factorization")
gate(sp.degree(H_table, e) == 11, "residual degree eleven")
for j, coefficient in enumerate(coefficients):
    zero(sp.expand(H_table).coeff(e, j) - coefficient,
         f"residual coefficient e^{j}")

# If all residual roots lay on e*g=0, factorization over the algebraic
# closure would be mu*e^a*g^(11-a).  The top ratio forces a=8, whereas the
# first three coefficients force ord_0(H)<=2.
h0, h1, h2 = coefficients[:3]
zero(
    h2.subs(b1, -2916 * b0**3) + 2916 * 2187 * b0**5,
    "h0/h1 seam leaves h2 nonzero",
)
top_ratio_cross = sp.expand(coefficients[10] * b1
                            - 3 * b0 * coefficients[11])
zero(top_ratio_cross, "actual next-to-top ratio is 3*b0/b1")
for order in range(3):
    gate(11 - order != 3, f"boundary factor ratio contradiction order={order}")


# Common-root boundary and exact source reconstruction.  At the off-boundary
# residual root, Q keeps degree four, so the homogeneous resultant has no
# root at infinity even if P happens to drop from degree two.
zero(Q.subs(u, 0) - e**2, "u=0 boundary")
zero(Q.subs(u, -1) + e**2, "u=-1 boundary")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * g**3 / 16,
     "K=0 boundary")
gate(sp.factor(sp.LC(sp.Poly(Q, u))) == -729 * g**3,
     "Q infinity value nonzero off g=0")

r_rec = u / e
z_rec = 9 * g * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * g) + Q / (9 * g * e**2 * K**2),
     "Q reconstructs z square")
zero(
    surface.subs({r: r_rec, z: z_rec})
    - u * (1 + u) * Q / (e**3 * K**3),
    "Q reconstructs surface",
)
zero(
    critical[2].subs({r: r_rec, z: z_rec})
    + Q / (e**2 * K**2),
    "Q kills e Hamiltonian",
)
zero(
    critical[1].subs({r: r_rec, z: z_rec}) - 3 * P / e**2,
    "P kills z Hamiltonian",
)
casimir_rec = sp.factor(
    K * critical[0].subs({r: r_rec, z: z_rec})
    - 3 * z_rec**2 * critical[1].subs({r: r_rec, z: z_rec})
    + r_rec**2 * critical[2].subs({r: r_rec, z: z_rec})
)
zero(casimir_rec, "reconstructed Casimir relation")


# Exact specialized controls exercise both a generic row and the seam where
# P loses its quadratic term at e=b0/b1; Q still excludes infinity.
control_rows: list[str] = []
for c0 in range(1, 10):
    for c1 in range(1, 10):
        Hc = sp.Poly(H_table.subs({b0: c0, b1: c1}), e)
        gate(Hc.degree() == 11, f"control degree {c0},{c1}")
        gate(Hc.LC() != 0, f"control leading {c0},{c1}")
        gate(sp.rem(Hc, sp.Poly(e * (c0 + c1 * e), e)) != 0,
             f"control is not supported as a simple boundary factor {c0},{c1}")
        control_rows.append(f"{c0}:{c1}:{Hc.TC()}:{Hc.LC()}")

semantic = {
    "carrier": "A=e^2-z/3+r*(b0+b1*e); c=1; (b0,b1)!=(0,0)",
    "resultant": "Res_u(P,Q)=e^4*(b0+b1e)^3*H; deg_e(H)=11",
    "boundary": "if roots(H) subset V(e*g), top ratio forces ord_0(H)=8 but ord_0(H)<=2",
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "quadratic-or-higher nonmonomial g; mixed r/z^2 corrections; different arm profile",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
control_blob = "\n".join(control_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3803-affine-linear-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(b0+b1*e);(b0,b1)!=(0,0)")
print("resultant=e^4*(b0+b1*e)^3*H(e);degree_H=11")
print("top_ratio=[e10]H/[e11]H=3*b0/b1")
print("origin_order=ord_0(H)<=2;boundary_only_would_force_ord_0(H)=8")
print("root=off_e*(b0+b1e)=0;common_u_finite;u_notin_0,-1,-1/2")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=degree_at_least_2_nonmonomial_g;mixed_r_z2;different_arm_profile")
print(f"control_sha256={hashlib.sha256(control_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
