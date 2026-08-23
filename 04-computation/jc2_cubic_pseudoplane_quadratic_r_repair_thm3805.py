#!/usr/bin/env python3
"""Exact companion for THM-3805's uniform quadratic r-repair obstruction."""

from __future__ import annotations

import ast
import hashlib
import itertools
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
b0, b1, b2 = sp.symbols("b0 b1 b2")
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


g = b0 + b1 * e + b2 * e**2
gp = sp.diff(g, e)
A = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A, q)) for q in variables]
zero(critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
     "quadratic r Hamiltonian")
zero(critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
     "quadratic z Hamiltonian")
zero(critical[2] - (9 * g * z**2 - K_source),
     "quadratic e Hamiltonian")
zero(bracket(A, surface), "surface Casimir")


# Universal two-equation compression and residual resultant.
G, D = sp.symbols("G D")
K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u)**2)
resultant_generic = sp.factor(sp.resultant(P, Q, u))
H_generic = sp.factor(resultant_generic / (G**3 * e**4))
zero(resultant_generic - G**3 * e**4 * H_generic,
     "universal resultant factor")

H = sp.Poly(sp.expand(H_generic.subs({G: g, D: gp})), e)
gate(H.degree() == 14, "quadratic residual degree fourteen")
gate(H.LC() == 8503056 * b2**3, "quadratic residual leading coefficient")
gate(sp.expand(H.as_expr()).subs({b0: 2, b1: -3, b2: 5}) != 0,
     "positive residual control")


# If every H-root lies in V(e*g), factorization over an algebraically closed
# field makes H divide e*g*H'.  Compute the entire remainder and freeze the
# four top equations that already contradict b2!=0.
quotient, remainder = sp.div(
    sp.Poly(sp.expand(e * g * sp.diff(H.as_expr(), e)), e), H
)
expected_quotient = (
    16 * b0 + 2 * b1 * b2 + b2**3
    + (22 * b1 + 2 * b2**2) * e + 28 * b2 * e**2
) / 2
zero(quotient.as_expr() - expected_quotient, "logarithmic quotient")
gate(remainder.degree() == 13, "logarithmic remainder degree")

E13 = 8 * b0 - 2 * b1 * b2 - b2**3
E12 = (
    72 * b0 * b1 - 4 * b0 * b2**2 - 12 * b1**2 * b2
    - 4 * b1 * b2**3 + b2**5
)
E11 = (
    17496 * b0**2 * b2 + 29160 * b0 * b1**2
    - 7776 * b0 * b1 * b2**2 - 1458 * b0 * b2**4
    - 2916 * b1**3 * b2 + 729 * b1 * b2**5 + 80
)
E10 = (
    58320 * b0**2 * b1 * b2 - 3888 * b0**2 * b2**3
    + 21384 * b0 * b1**3 - 14580 * b0 * b1**2 * b2**2
    - 2430 * b0 * b1 * b2**4 + 243 * b0 * b2**6
    - 972 * b1**4 * b2 + 972 * b1**3 * b2**3
    + 729 * b1**2 * b2**5 + 64 * b1 - 56 * b2**2
)
zero(remainder.coeff_monomial(e**13) + 2125764 * b2**4 * E13,
     "top remainder E13")
zero(remainder.coeff_monomial(e**12) + 1062882 * b2**3 * E12,
     "top remainder E12")
zero(remainder.coeff_monomial(e**11) + 4374 * b2**2 * E11,
     "top remainder E11")
zero(remainder.coeff_monomial(e**10) + 4374 * b2 * E10,
     "top remainder E10")


# Exact two-branch contradiction after E13=E12=0.
b0_from_E13 = b2 * (2 * b1 + b2**2) / 8
zero(E13.subs(b0, b0_from_E13), "solve E13")
E12_reduced = sp.factor(E12.subs(b0, b0_from_E13))
zero(
    E12_reduced - b2 * (2 * b1 + b2**2) * (6 * b1 + b2**2) / 2,
    "E12 two-branch factorization",
)

branch_one = {b1: -b2**2 / 2, b0: sp.Integer(0)}
zero(E13.subs(branch_one), "branch one E13")
zero(E12.subs(branch_one), "branch one E12")
gate(E11.subs(branch_one) == 80, "branch one E11 contradiction")

branch_two = {b1: -b2**2 / 6, b0: b2**3 / 12}
zero(E13.subs(branch_two), "branch two E13")
zero(E12.subs(branch_two), "branch two E12")
E11_branch_two = sp.factor(E11.subs(branch_two))
E10_branch_two = sp.factor(E10.subs(branch_two))
zero(E11_branch_two - 5 * (27 * b2**7 + 32) / 2,
     "branch two E11 law")
zero(E10_branch_two + 5 * b2**2 * (81 * b2**7 + 80) / 6,
     "branch two E10 law")
zero(E10_branch_two.subs(b2**7, -sp.Rational(32, 27))
     - sp.Rational(40, 3) * b2**2,
     "branch two final contradiction")


# Exact reconstruction from any residual root eta with eta*g(eta)!=0.
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * G) + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
zero(surface.subs({r: r_rec, z: z_rec})
     - u * (1 + u) * Q / (e**3 * K**3),
     "Q reconstructs surface")
A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")
C_r_rec = sp.factor(r_rec**2 - 9 * z_rec**2 * A_e_symbol)
zero(
    C_r_rec - P / (G * e**2)
    - Q * (u**2 - P / G) / (e**4 * K**3),
    "P and Q kill r Hamiltonian",
)
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16,
     "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3,
     "Q leading coefficient excludes infinity")


# Hostile controls only: exhaustive finite-field scans of the full
# logarithmic remainder system.  They are not used in the characteristic-
# zero proof above.
remainder_numerators = [
    sp.together(remainder.coeff_monomial(e**j)).as_numer_denom()[0]
    for j in range(14)
]
remainder_evaluator = sp.lambdify((b0, b1, b2), remainder_numerators, "math")
scan_rows: list[str] = []
for prime in (5, 7, 11):
    survivors = 0
    for aa, bb, cc in itertools.product(range(prime), repeat=3):
        if cc == 0:
            continue
        values = remainder_evaluator(aa, bb, cc)
        if all(int(value) % prime == 0 for value in values):
            survivors += 1
    gate(survivors == 0, f"finite-field hostile p={prime}")
    scan_rows.append(f"{prime}:{survivors}")


semantic = {
    "carrier": "A=e^2-z/3+r*(b0+b1e+b2e^2); c=1",
    "new_case": "b2!=0; deg<=1 inherited from THM-3803",
    "resultant": "Res_u(P,Q)=e^4*g^3*H; deg(H)=14",
    "mechanism": "boundary-only roots imply H divides e*g*Hprime; top remainder has two impossible branches",
    "branches": "2b1+b2^2=0 gives E11=80; 6b1+b2^2=0 gives incompatible E11/E10",
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "degree>=3 g; mixed z^2/r corrections; different arm profile",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
scan_blob = "\n".join(scan_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3805-quadratic-r-repairs-of-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(b0+b1e+b2e^2)")
print("new_case=b2!=0;degree_le_1_inherited=THM-3803")
print("resultant=Res_u(P,Q)=e^4*g^3*H;deg_H=14;LC_H=8503056*b2^3")
print("boundary_criterion=all_H_roots_in_V(e*g)_implies_H_divides_e*g*Hprime")
print("branch_one=2b1+b2^2=0;E11=80")
print("branch_two=6b1+b2^2=0;b0=b2^3/12;27b2^7+32=0;E10=40b2^2/3")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=degree_ge_3_g;mixed_z2_r;different_arm_profile;Darboux_pair")
print(f"finite_field_scan_sha256={hashlib.sha256(scan_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
