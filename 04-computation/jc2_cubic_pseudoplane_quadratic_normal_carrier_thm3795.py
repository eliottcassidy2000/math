#!/usr/bin/env python3
"""Exact companion for THM-3795's quadratic-normal carrier obstruction."""

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


# Generic Hamiltonian identities for the full r-independent normal form.
f = sp.Function("f")(e)
h = sp.Function("h")(e)
A = e**2 + z * f + z**2 * h
A_e = sp.diff(A, e)
A_z = sp.diff(A, z)
K = c**3 + 2 * r * e
zero(bracket(A, r) - (-3 * r**2 * A_z - 9 * z**2 * A_e),
     "generic r Hamiltonian component")
zero(bracket(A, z) + 3 * K * A_e, "generic z Hamiltonian component")
zero(bracket(A, e) - 3 * K * A_z, "generic e Hamiltonian component")


# The dangerous curve K=0 is exactly a torus in the z-coordinate.
a = -c**6 / 4
e_torus = a / z**3
r_torus = 2 * z**3 / c**3
surface = r**2 * e - z**3 + c**3 * r
zero(surface.subs({r: r_torus, e: e_torus}), "torus lies on surface")
zero(K.subs({r: r_torus, e: e_torus}), "torus lies on K=0")
e_prime = sp.diff(e_torus, z)
zero(e_prime - 3 * c**6 / (4 * z**4), "torus e derivative")
zero(r_torus**2 * e_prime - 3 * z**2, "Hamiltonian-chain factor")

# Use independent placeholders for the two partial derivatives; this avoids
# relying on SymPy's representation of derivatives of composed Functions.
Ae_symbol, Az_symbol = sp.symbols("Ae_symbol Az_symbol")
phi_prime_chain = Az_symbol + e_prime * Ae_symbol
remaining = -3 * r_torus**2 * Az_symbol - 9 * z**2 * Ae_symbol
zero(remaining + 3 * r_torus**2 * phi_prime_chain,
     "remaining Hamiltonian equals -3r^2 Phi prime")


# Exact coefficient/exponent compiler for the Laurent derivative.  The rows
# are mathematical data, not random samples: every possible monomial index
# through 256 is checked, and the proof handles arbitrary finite degrees by
# the same closed formulas.
support_rows: list[str] = []
e2_exponent = -7
e2_coefficient = -3 * c**12 / 8
gate(e2_exponent % 3 == 2, "e2 residue class")
gate(e2_coefficient != 0, "e2 mandatory coefficient")

f_exponents: set[int] = set()
h_exponents: set[int] = set()
for j in range(257):
    f_exponent = -3 * j
    h_exponent = 1 - 3 * j
    f_factor = 1 - 3 * j
    h_factor = 2 - 3 * j
    gate(f_exponent % 3 == 0, f"f residue j={j}")
    gate(h_exponent % 3 == 1, f"h residue j={j}")
    gate(f_factor != 0, f"f derivative factor j={j}")
    gate(h_factor != 0, f"h derivative factor j={j}")
    gate(f_exponent != e2_exponent, f"f/e2 noncollision j={j}")
    gate(h_exponent != e2_exponent, f"h/e2 noncollision j={j}")
    gate(f_exponent not in h_exponents, f"f/h prior noncollision j={j}")
    gate(h_exponent not in f_exponents, f"h/f prior noncollision j={j}")
    f_exponents.add(f_exponent)
    h_exponents.add(h_exponent)
    support_rows.append(
        f"{j}:{f_exponent}:{f_factor}:{h_exponent}:{h_factor}"
    )

gate(len(f_exponents) == 257, "f exponents are unique")
gate(len(h_exponents) == 257, "h exponents are unique")
gate(f_exponents.isdisjoint(h_exponents), "f/h support classes disjoint")
gate(e2_exponent not in f_exponents | h_exponents,
     "e2 support class disjoint")
gate(0 in f_exponents, "f constant contributes exponent zero")


# Polynomial controls with every coefficient active.  They verify the exact
# Laurent formula, the two mandatory terms, and clearing to a nonconstant
# polynomial having nonzero constant term for a grid of unequal degrees.
degree_rows: list[str] = []
for m in range(17):
    for n in range(17):
        f_poly = sp.Integer(5) + sum(
            ((j + 2) * e**j for j in range(1, m + 1)), sp.Integer(0)
        )
        h_poly = sum(
            ((2 * j + 3) * e**j for j in range(n + 1)), sp.Integer(0)
        )
        phi = sp.expand(e_torus**2 + z * f_poly.subs(e, e_torus)
                        + z**2 * h_poly.subs(e, e_torus))
        derivative = sp.expand(sp.diff(phi, z))
        expected = -3 * c**12 / (8 * z**7)
        expected += sum(
            (1 - 3 * j) * (5 if j == 0 else j + 2) * a**j * z**(-3 * j)
            for j in range(m + 1)
        )
        expected += sum(
            (2 - 3 * j) * (2 * j + 3) * a**j * z**(1 - 3 * j)
            for j in range(n + 1)
        )
        zero(derivative - expected, f"Laurent formula m={m},n={n}")

        exponents = [e2_exponent]
        exponents.extend(-3 * j for j in range(m + 1))
        exponents.extend(1 - 3 * j for j in range(n + 1))
        low = min(exponents)
        cleared = sp.cancel(z**(-low) * derivative)
        cleared_poly = sp.Poly(cleared, z)
        gate(cleared_poly.as_expr().has(z),
             f"cleared derivative nonconstant m={m},n={n}")
        gate(cleared_poly.TC() != 0,
             f"cleared derivative nonzero constant m={m},n={n}")
        gate(sp.expand(derivative).coeff(z, 0) == 5,
             f"mandatory f0 coefficient m={m},n={n}")
        gate(sp.expand(derivative).coeff(z, e2_exponent) == e2_coefficient,
             f"mandatory e2 coefficient m={m},n={n}")
        degree_rows.append(f"{m}:{n}:{low}:{cleared_poly.degree()}")


# THM-3792 boundary: when h=0, the Laurent root condition reconstructs its
# critical polynomial.  We check the canonical constant first-normal row.
f0 = -1 / (3 * c**3)
phi0_prime = sp.factor(sp.diff(e_torus**2 + z * f0, z))
zero(phi0_prime - (-3 * c**12 / (8 * z**7) + f0),
     "canonical first-normal Laurent derivative")
canonical_equation = sp.factor(sp.together(phi0_prime * z**7))
gate(sp.degree(canonical_equation, z) == 7,
     "canonical seven-point boundary")
gate(canonical_equation.subs(z, 0) != 0,
     "canonical roots stay on torus")


# Sharp boundary: z^3 e^3 is the first new monomial whose derivative can
# collide with the e^2 term.  It cancels e^2 on K=0, but reduction by the
# surface relation makes its genuine r-dependence explicit.
sharp = e**2 - z / (3 * c**3) + 4 * z**3 * e**3 / c**6
sharp_torus = sp.factor(sharp.subs(e, e_torus))
zero(sharp_torus + z / (3 * c**3), "sharp torus cancellation")
zero(sp.diff(sharp_torus, z) + 1 / (3 * c**3),
     "sharp torus derivative is nonzero constant")
normal_reduction = 4 * (r**2 * e + c**3 * r) * e**3 / c**6
zero(normal_reduction - (4 * r**2 * e**4 / c**6 + 4 * r * e**3 / c**3),
     "sharp term has r-coupled normal form")
zero((4 * z**3 * e**3 / c**6 - normal_reduction).subs(
    z**3, r**2 * e + c**3 * r), "surface reduction of sharp term")
gate((3 - 1 - 3 * 3) == e2_exponent,
     "first cubic collision exponent")


# Cheapest genuinely r-coupled hostile control.  At c=1 the constant row
# A=e^2-z/3+b*r still has a critical point for every b!=0.  After putting
# u=r/z and eliminating e with the surface equation, two critical equations
# are p1=p2=0.  Their exact resultant has nonzero constant and positive
# degree; a common root has u,z nonzero, and the third Hamiltonian follows
# from the Casimir relation because F_z=-3z^2.
u, b = sp.symbols("u b", nonzero=True)
p1 = 18 * z**2 - u**4 * z - 18 * u
p2 = 81 * b * z**2 - u**3 * z - 9
gate(sp.Poly(p1, z).LC() == 18, "p1 has fixed quadratic degree")
gate(sp.Poly(p2, z).LC() == 81 * b, "p2 has fixed quadratic degree")
resultant = sp.factor(sp.resultant(p1, p2, z))
expected_resultant = 81 * (
    26244 * b**2 * u**2 + 9 * b * u**8 - 5832 * b * u
    - 2 * u**7 + 324
)
zero(resultant - expected_resultant, "constant r correction resultant")
gate(sp.Poly(resultant / 81, u).degree() == 8,
     "constant r resultant is nonconstant")
gate(sp.Poly(resultant / 81, u).TC() == 324,
     "constant r resultant roots have u nonzero")
gate(p2.subs(z, 0) == -9, "constant r common root has z nonzero")

e_from_surface = (z**2 - u) / (u**2 * z)
r_from_u = u * z
A_constant_r = e**2 - z / 3 + b * r
critical_constant_r = [
    sp.factor(bracket(A_constant_r, q).subs(c, 1).subs({r: r_from_u, e: e_from_surface}))
    for q in variables
]
zero(sp.factor(critical_constant_r[0] * u**2 / z) + p1,
     "p1 is the r Hamiltonian equation")
second_reduction = sp.factor(critical_constant_r[2] * u)
zero(9 * second_reduction - u * p2 + p1,
     "p2 is equivalent to the e Hamiltonian modulo p1")
surface_sub = surface.subs(c, 1).subs({r: r_from_u, e: e_from_surface})
zero(surface_sub, "constant r reconstruction lies on surface")
zero(bracket(A_constant_r, surface).subs(c, 1),
     "third Hamiltonian follows from the surface Casimir")


semantic = {
    "carrier": "A=e^2+z*f(e)+z^2*h(e), f(0)!=0",
    "critical_curve": "K=c^3+2re=0; e=-c^6/(4z^3); r=2z^3/c^3",
    "mechanism": "Phi_prime has disjoint exponent classes 2,0,1 modulo 3",
    "mandatory_terms": "(-3c^12/8)z^-7 and f(0)",
    "sharp_boundary": "z^3e^3 cancels the witness and is r-coupled in canonical normal form",
    "open": "nonconstant r-coupled coefficients; different arm profile",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
support_blob = "\n".join(support_rows).encode()
degree_blob = "\n".join(degree_rows).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3795-r-independent-quadratic-normal-nodal-carriers-have-critical-points")
print("surface=r^2e-z^3+c^3r;field=algebraically_closed_characteristic_zero;c!=0")
print("carrier=A=e^2+z*f(e)+z^2*h(e);f(0)!=0")
print("canonical_normal_form=k[r,e]+z*k[r,e]+z^2*k[r,e]")
print("critical_torus=K=c^3+2re=0;e=-c^6/(4z^3);r=2z^3/c^3")
print("laurent_residues=e2:2;zf:0;z2h:1_mod_3")
print("mandatory_terms=(-3c^12/8)z^-7;f(0)")
print("sharp_boundary=z^3e^3_collision;canonical_normal_form_is_r_coupled")
print("constant_r_hostile=c=1;all_b_nonzero;resultant_degree=8;constant=324")
print("open=nonconstant_r_coupled_coefficients;different_arm_profile;Darboux_pair")
print(f"support_sha256={hashlib.sha256(support_blob).hexdigest()}")
print(f"degree_sha256={hashlib.sha256(degree_blob).hexdigest()}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
