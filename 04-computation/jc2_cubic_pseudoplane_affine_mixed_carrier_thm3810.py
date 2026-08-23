#!/usr/bin/env python3
"""Exact companion for THM-3810's affine-profile mixed-carrier obstruction."""

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
a, b, h = sp.symbols("a b h")
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


A = e**2 - z / 3 + r * (a * e + b) + h * z**2
L = 1 - 6 * h * z
K = 1 + 2 * r * e
g = a * e + b
critical = [sp.factor(bracket(A, q)) for q in variables]
expected_critical = (
    L * r**2 - 9 * z**2 * (2 * e + a * r),
    3 * g * r**2 - 3 * K * (2 * e + a * r),
    9 * g * z**2 - L * K,
)
for index, (actual, expected) in enumerate(zip(critical, expected_critical)):
    zero(actual - expected, f"Hamiltonian component {index}")
zero(bracket(A, surface), "surface is a Casimir")
zero(
    K * critical[0] - 3 * z**2 * critical[1] + r**2 * critical[2],
    "Hamiltonian Casimir syzygy",
)


# With u=r/z, C_r=0 solves linearly for e even after the affine r-profile.
M = sp.expand(L * u**2 - 9 * a * u * z)
p1 = sp.expand(18 * z**2 - u**2 * z * M - 18 * u)
p2 = sp.expand(162 * b * z**2 + 9 * a * M * z**2 - 18 * L - 2 * L * u * z * M)
resultant = sp.factor(sp.resultant(p1, p2, u))
T = sp.expand(
    32 * L**5 * z**7
    + 36 * L**4
    + 324 * L**3 * a**2 * z**5
    - 3888 * L**3 * a * b * z**8
    + 648 * L**3 * b * z**2
    + 729 * L**2 * a**4 * z**10
    - 2916 * L**2 * a**2 * b * z**7
    - 13122 * L * a**3 * b * z**6
    + 104976 * L * a**2 * b**2 * z**9
    - 52488 * L * b**3 * z**6
    - 59049 * a**5 * b * z**11
    + 118098 * a**3 * b**2 * z**8
    - 236196 * b**4 * z**8
)
zero(
    resultant + 2916 * z**3 * (6 * h * z - 1) ** 3 * T,
    "affine-profile slope resultant",
)
gate(sp.Poly(T, z).TC() == 36, "residual constant coefficient")


# If h!=0 and all residual roots lay on L=0, algebraic closedness and
# T(0)=36 would force T=36*L^d.  Coefficients successively force d=4,
# b=0, a=0, and then leave a nonzero z^7*L^5 term.
zero(sp.expand(T).coeff(z, 1) + 864 * h, "boundary-only degree witness")
zero(
    sp.expand(T - 36 * L**4).coeff(z, 2) - 648 * b,
    "boundary-only b witness",
)
zero(
    sp.expand((T - 36 * L**4).subs(b, 0)).coeff(z, 5) - 324 * a**2,
    "boundary-only a witness",
)
zero(
    T.subs({a: 0, b: 0}) - 36 * L**4 - 32 * L**5 * z**7,
    "terminal pure-power contradiction",
)
gate(sp.Poly(32 * L**5 * z**7, z).degree() == 12,
     "terminal contradiction is nonzero for h!=0")


# At h=0 there is no L-boundary.  Recompute the resultant, and prove the
# residual is nonconstant from its z^2 coefficient if b!=0 and its z^7
# coefficient if b=0.
T_h0 = sp.expand(T.subs(h, 0))
gate(sp.Poly(T_h0, z).TC() == 36, "h=0 residual constant")
zero(sp.expand(T_h0).coeff(z, 2) - 648 * b, "h=0,b!=0 degree witness")
gate(sp.expand(T_h0.subs(b, 0)).coeff(z, 7) == 32,
     "h=0,b=0 degree witness")
zero(
    sp.resultant(p1.subs(h, 0), p2.subs(h, 0), u)
    - resultant.subs(h, 0),
    "h=0 specialized resultant",
)


# Safe residual roots retain the quartic/cubic u-degrees and reconstruct a
# source critical point without a parameter denominator.
zero(sp.LC(sp.Poly(p1, u)) + L * z, "p1 leading coefficient")
zero(sp.LC(sp.Poly(p2, u)) + 2 * L**2 * z, "p2 leading coefficient")
gate(sp.Poly(p1, u).degree() == 4, "p1 u-degree")
gate(sp.Poly(p2, u).degree() == 3, "p2 u-degree")
gate(sp.Poly(p1.subs(u, 0), z) == sp.Poly(18 * z**2, z), "exclude u=0")

r_rec = u * z
e_rec = M / 18
reconstruction = {r: r_rec, e: e_rec}
zero(
    surface.subs(reconstruction) + z * p1 / 18,
    "p1 reconstructs the source surface",
)
zero(critical[0].subs(reconstruction), "reconstruction kills r-Hamiltonian")
zero(
    critical[2].subs(reconstruction) - p2 / 18,
    "p2 kills e-Hamiltonian",
)
zero(
    critical[1].subs(reconstruction) - u**2 * p2 / 54,
    "p2 kills z-Hamiltonian",
)


semantic = {
    "carrier": "A=e^2-z/3+r*(a*e+b)+h*z^2; all a,b,h",
    "compression": "u=r/z; L=1-6hz; M=L*u^2-9auz; e=M/18",
    "equations": "p1=18z^2-u^2zM-18u; p2=162bz^2+9aMz^2-18L-2LuzM",
    "resultant": "-2916z^3(6hz-1)^3T",
    "boundary": "T0=36; z1=-864h forces d=4; z2 excess=648b; z5 excess=324a^2; terminal=32L^5z^7",
    "h0": "z2=648b; after b=0, z7=32",
    "recovery": "r=uz; e=M/18; F=-zp1/18; Cr=0; Ce=p2/18; Cz=u^2p2/54",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3810-affine-r-profile-constant-z2-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(a*e+b)+h*z^2;all_a_b_h")
print("compression=u=r/z;L=1-6hz;M=L*u^2-9auz;e=M/18")
print("equations=p1=18z^2-u^2zM-18u;p2=162bz^2+9aMz^2-18L-2LuzM")
print("resultant=-2916z^3(6hz-1)^3T")
print("boundary=T0=36;z1=-864h=>d4;z2_excess=648b;z5_excess=324a^2;terminal=32L^5z^7")
print("h0=z2=648b;b0_z7=32")
print("recovery=r=uz;e=M/18;F=-zp1/18;Cr=0;Ce=p2/18;Cz=u^2p2/54")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
