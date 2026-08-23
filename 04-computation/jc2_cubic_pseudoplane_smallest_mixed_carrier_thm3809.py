#!/usr/bin/env python3
"""Exact companion for THM-3809's smallest mixed-carrier obstruction."""

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
b, h = sp.symbols("b h")
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


# Hamiltonian packet and its Casimir relation.
A = e**2 - z / 3 + b * r + h * z**2
L = 1 - 6 * h * z
K = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
expected_critical = (
    L * r**2 - 18 * e * z**2,
    3 * b * r**2 - 6 * e * K,
    9 * b * z**2 - L * K,
)
for index, (actual, expected) in enumerate(zip(critical, expected_critical)):
    zero(actual - expected, f"Hamiltonian component {index}")
zero(bracket(A, surface), "surface is a Casimir")
zero(
    K * critical[0] - 3 * z**2 * critical[1] + r**2 * critical[2],
    "Hamiltonian Casimir syzygy",
)


# At a critical point z cannot vanish.  With u=r/z, the first Hamiltonian
# gives e=L*u^2/18.  Surface membership and the e-Hamiltonian become p1,p2.
p1 = sp.expand(18 * z**2 - u**4 * L * z - 18 * u)
p2 = sp.expand(81 * b * z**2 - u**3 * L**2 * z - 9 * L)
resultant = sp.factor(sp.resultant(p1, p2, u))
R = sp.expand(
    59049 * b**4 * z**8
    + 13122 * b**3 * z**6 * L
    - 162 * b * z**2 * L**3
    - 8 * z**7 * L**5
    - 9 * L**4
)
zero(
    resultant - 729 * z**3 * (6 * h * z - 1) ** 3 * R,
    "slope resultant factorization",
)


# Generic mixed cell b*h!=0: R itself has a root away from z*L=0.
R_poly = sp.Poly(R, z)
gate(R_poly.degree() == 12, "generic residual degree")
zero(R_poly.LC() - 62208 * h**5, "generic residual leading coefficient")
gate(R_poly.TC() == -9, "generic residual constant coefficient")
zero(
    R.subs(z, 1 / (6 * h)) - 9 * b**4 / (256 * h**8),
    "generic residual excludes L=0",
)


# Parameter axes are recomputed in their specialized coefficient rings.
S = sp.expand(48 * h * z**8 - 8 * z**7 - 9)
zero(R.subs(b, 0) - L**4 * S, "b=0 residual factorization")
gate(sp.Poly(S, z).degree() == 8, "b=0 residual degree")
zero(sp.Poly(S, z).LC() - 48 * h, "b=0 residual leading coefficient")
gate(sp.Poly(S, z).TC() == -9, "b=0 residual constant coefficient")
gate(sp.cancel(S.subs(z, 1 / (6 * h))) == -9, "b=0 excludes L=0")
zero(
    sp.resultant(p1.subs(b, 0), p2.subs(b, 0), u)
    - resultant.subs(b, 0),
    "b=0 specialized resultant",
)

R_h0 = sp.Poly(sp.expand(R.subs(h, 0)), z)
gate(R_h0.degree() == 8, "h=0,b!=0 residual degree")
zero(R_h0.LC() - 59049 * b**4, "h=0,b!=0 leading coefficient")
gate(R_h0.TC() == -9, "h=0,b!=0 constant coefficient")
zero(
    sp.resultant(p1.subs(h, 0), p2.subs(h, 0), u)
    - resultant.subs(h, 0),
    "h=0 specialized resultant",
)

R_origin = sp.expand(R.subs({b: 0, h: 0}))
gate(R_origin == -8 * z**7 - 9, "parameter-origin residual")
gate(sp.Poly(R_origin, z).degree() == 7, "parameter-origin degree")
gate(sp.Poly(R_origin, z).TC() == -9, "parameter-origin constant")
zero(
    sp.resultant(p1.subs({b: 0, h: 0}), p2.subs({b: 0, h: 0}), u)
    - resultant.subs({b: 0, h: 0}),
    "parameter-origin specialized resultant",
)


# Every safe residual root gives a finite common u: both u-degrees remain
# fixed away from z*L=0.  Reconstruction uses no division at all.
zero(sp.LC(sp.Poly(p1, u)) + L * z, "p1 leading coefficient")
zero(sp.LC(sp.Poly(p2, u)) + L**2 * z, "p2 leading coefficient")
gate(sp.Poly(p1, u).degree() == 4, "p1 u-degree")
gate(sp.Poly(p2, u).degree() == 3, "p2 u-degree")
gate(sp.Poly(p1.subs(u, 0), z) == sp.Poly(18 * z**2, z), "exclude u=0")

r_rec = u * z
e_rec = L * u**2 / 18
reconstruction = {r: r_rec, e: e_rec}
zero(
    surface.subs(reconstruction) + z * p1 / 18,
    "p1 reconstructs the source surface",
)
zero(critical[0].subs(reconstruction), "reconstruction kills r-Hamiltonian")
zero(
    critical[2].subs(reconstruction) - p2 / 9,
    "p2 kills e-Hamiltonian",
)
zero(
    critical[1].subs(reconstruction) - u**2 * p2 / 27,
    "p2 kills z-Hamiltonian",
)


semantic = {
    "carrier": "A=e^2-z/3+b*r+h*z^2; all b,h",
    "compression": "u=r/z; L=1-6hz; e=L*u^2/18",
    "equations": "p1=18z^2-u^4Lz-18u; p2=81bz^2-u^3L^2z-9L",
    "resultant": "729z^3(6hz-1)^3R",
    "residual": "R=59049b^4z^8+13122b^3z^6L-162bz^2L^3-8z^7L^5-9L^4",
    "generic": "degR=12; LC=62208h^5; R(0)=-9; R(L=0)=9b^4/(256h^8)",
    "axes": "b=0:S=48hz^8-8z^7-9; h=0:degR=8; origin:R=-8z^7-9",
    "recovery": "r=uz; e=L*u^2/18; F=-zp1/18; Cr=0; Ce=p2/9; Cz=u^2p2/27",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3809-smallest-mixed-r-z2-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+b*r+h*z^2;all_b_h")
print("compression=u=r/z;L=1-6hz;e=L*u^2/18")
print("equations=p1=18z^2-u^4Lz-18u;p2=81bz^2-u^3L^2z-9L")
print("resultant=729z^3(6hz-1)^3R")
print("R=59049b^4z^8+13122b^3z^6L-162bz^2L^3-8z^7L^5-9L^4")
print("generic=deg12;LC=62208h^5;R0=-9;R_L0=9b^4/(256h^8)")
print("axes=b0:S=48hz^8-8z^7-9;h0:deg8;origin:R=-8z^7-9")
print("recovery=r=uz;e=L*u^2/18;F=-zp1/18;Cr=0;Ce=p2/9;Cz=u^2p2/27")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
