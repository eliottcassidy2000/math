#!/usr/bin/env python3
"""Exact companion for THM-3807's >=3-support cubic r-repair obstruction."""

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
b0, b1, b2, t = sp.symbols("b0 b1 b2 t")
G, D = sp.symbols("G D")
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


# Universal Hamiltonian compression.
g = b0 + b1 * e + b2 * e**2 + t * e**3
gp = sp.diff(g, e)
A_carrier = e**2 - z / 3 + r * g
K_source = 1 + 2 * r * e
critical = [sp.factor(bracket(A_carrier, q)) for q in variables]
zero(critical[0] - (r**2 - 9 * z**2 * (2 * e + r * gp)),
     "cubic r Hamiltonian")
zero(critical[1] - (3 * g * r**2 - 3 * K_source * (2 * e + r * gp)),
     "cubic z Hamiltonian")
zero(critical[2] - (9 * g * z**2 - K_source),
     "cubic e Hamiltonian")
zero(bracket(A_carrier, surface), "surface Casimir")

K = 1 + 2 * u
P = sp.expand(G * u**2 - K * (2 * e**3 + u * e * D))
Q = sp.expand(e**2 * K**3 - 729 * G**3 * u**2 * (1 + u)**2)
resultant_generic = sp.factor(sp.resultant(P, Q, u))
# Independent fixed-degree Sylvester determinant, with four shifted P rows
# and two shifted Q rows.
p_coefficients = sp.Poly(P, u).all_coeffs()
q_coefficients = sp.Poly(Q, u).all_coeffs()
sylvester_rows = []
for shift in range(4):
    sylvester_rows.append(
        [0] * shift + p_coefficients + [0] * (3 - shift)
    )
for shift in range(2):
    sylvester_rows.append(
        [0] * shift + q_coefficients + [0] * (1 - shift)
    )
zero(sp.factor(sp.Matrix(sylvester_rows).det()) - resultant_generic,
     "independent Sylvester determinant")
H_universal = sp.factor(resultant_generic / (G**3 * e**4))
zero(resultant_generic - G**3 * e**4 * H_universal,
     "universal resultant factor")
gate(len(sp.Poly(H_universal, G, D, e).terms()) == 20,
     "universal residual term count")


# Generic cubic stratum t != 0,1.  Polynomial division is over the localized
# coefficient ring Q[b0,b1,b2,t,t^-1,(t-1)^-1].
H_generic = sp.Poly(sp.expand(H_universal.subs({G: g, D: gp})), e)
gate(H_generic.degree() == 17, "generic cubic residual degree")
zero(H_generic.LC() - 8503056 * t**3 * (t - 1)**2,
     "generic cubic residual leading coefficient")
quotient, remainder = sp.div(
    sp.Poly(sp.expand(e * g * sp.diff(H_generic.as_expr(), e)), e),
    H_generic,
)
zero(
    e * g * sp.diff(H_generic.as_expr(), e)
    - quotient.as_expr() * H_generic.as_expr() - remainder.as_expr(),
    "generic Euclidean division identity",
)
gate(remainder.degree() == 16, "generic logarithmic remainder degree")

# If b2=0, >=3 support means b0*b1*t != 0, and the e^15 coefficient is
# already impossible.
zero(
    remainder.nth(15).subs(b2, 0)
    - 25509168 * b0 * b1 * t**3 * (t - 1),
    "generic b2=0 top contradiction",
)

# Assume b2 != 0.  The top coefficient solves b0, and the next coefficient
# factors into two exact branches A_factor=0 or B_factor=0.
E16 = (
    16 * b0 * (t - 1)**2 - 4 * b1 * b2 * (t - 1)**2
    + b2**3 * (t - 2)
)
zero(
    remainder.nth(16)
    + 1062882 * b2 * t**3 * E16 / (t - 1)**2,
    "generic top E16",
)
A_factor = 4 * b1 * (t - 1)**2 - b2**2 * (t - 2)
B_factor = 12 * b1 * (t - 1)**2 - b2**2 * (t - 2)
b0_solution = b2 * A_factor / (16 * (t - 1)**2)
zero(E16.subs(b0, b0_solution), "solve generic E16")
zero(
    remainder.nth(15).subs(b0, b0_solution)
    - 531441 * b2 * t**3 * A_factor * B_factor
    / (4 * (t - 1)**3),
    "generic second coefficient two-branch factorization",
)

# Branch A_factor=0.  Here b0=0 and b1 is forced.  The t=2 seam has only
# the two monomials e^2,e^3; outside it two coprime t-polynomials are forced.
branch_a = {
    b0: 0,
    b1: b2**2 * (t - 2) / (4 * (t - 1)**2),
}
q_a = 9 * t**2 - 20 * t + 20
p_a = 195 * t**4 - 838 * t**3 + 1380 * t**2 - 1080 * t + 352
zero(remainder.nth(14).subs(branch_a), "branch A e14 vanishes")
zero(
    remainder.nth(13).subs(branch_a)
    - 20412 * t**2 * (t - 2) * q_a,
    "branch A e13 law",
)
zero(
    remainder.nth(12).subs(branch_a)
    - 2916 * b2 * t * p_a / (t - 1),
    "branch A e12 law",
)
gate(sp.resultant(q_a, p_a, t) == 210898944,
     "branch A coprime resultant")
zero(sp.sympify(branch_a[b0]), "branch A forces b0=0")
zero(branch_a[b1].subs(t, 2), "branch A t=2 forces b1=0")

# Branch B_factor=0.  The e14 row gives t=2 or 10/9.  The former is the same
# binomial seam; the latter has a constant nonzero e13 remainder.
branch_b = {
    b0: -b2**3 * (t - 2) / (24 * (t - 1)**2),
    b1: b2**2 * (t - 2) / (12 * (t - 1)**2),
}
zero(
    remainder.nth(14).subs(branch_b)
    - 59049 * b2**6 * t**3 * (t - 2)**2 * (9 * t - 10)
    / (8 * (t - 1)**4),
    "branch B e14 law",
)
zero(branch_b[b0].subs(t, 2), "branch B t=2 forces b0=0")
zero(branch_b[b1].subs(t, 2), "branch B t=2 forces b1=0")
zero(
    remainder.nth(13).subs(branch_b).subs(t, sp.Rational(10, 9))
    + sp.Rational(1792000, 9),
    "branch B t=10/9 contradiction",
)


# Degree-drop stratum t=1 with b2 != 0.  Specialize H before division.
g_t1 = b0 + b1 * e + b2 * e**2 + e**3
H_t1 = sp.Poly(
    sp.expand(H_universal.subs({G: g_t1, D: sp.diff(g_t1, e)})), e
)
gate(H_t1.degree() == 15, "t=1,b2!=0 residual degree")
zero(H_t1.LC() - 2125764 * b2**2, "t=1,b2!=0 leading coefficient")
quotient_t1, remainder_t1 = sp.div(
    sp.Poly(sp.expand(e * g_t1 * sp.diff(H_t1.as_expr(), e)), e), H_t1
)
zero(
    e * g_t1 * sp.diff(H_t1.as_expr(), e)
    - quotient_t1.as_expr() * H_t1.as_expr() - remainder_t1.as_expr(),
    "t=1,b2!=0 Euclidean division",
)
zero(
    remainder_t1.nth(14) - 8503056 * b0 * (b0 + b1 * b2),
    "t=1 top branch law",
)
gate(remainder_t1.nth(13).subs(b0, 0) == -131220,
     "t=1,b0=0 contradiction")
t1_branch = {b0: -b1 * b2}
zero(
    remainder_t1.nth(13).subs(t1_branch)
    - 26244 * (648 * b1**2 * b2**3 - 5),
    "t=1 nonzero branch e13 law",
)
zero(
    remainder_t1.nth(12).subs(t1_branch)
    - 17496 * b2 * (2916 * b1**2 * b2**3 - 5),
    "t=1 nonzero branch e12 law",
)
zero(
    remainder_t1.nth(12).subs(t1_branch).subs(
        b1**2, sp.Rational(5, 648) / b2**3
    ) - 306180 * b2,
    "t=1 nonzero branch final contradiction",
)


# Deeper degree drop t=1,b2=0.  >=3 support now means b0*b1 != 0.
g_t1_b20 = b0 + b1 * e + e**3
H_t1_b20 = sp.Poly(
    sp.expand(H_universal.subs({G: g_t1_b20, D: sp.diff(g_t1_b20, e)})), e
)
gate(H_t1_b20.degree() == 11, "t=1,b2=0 residual degree")
zero(H_t1_b20.LC() - 2125764 * b0**2,
     "t=1,b2=0 leading coefficient")
quotient_t1_b20, remainder_t1_b20 = sp.div(
    sp.Poly(
        sp.expand(e * g_t1_b20 * sp.diff(H_t1_b20.as_expr(), e)), e
    ),
    H_t1_b20,
)
zero(
    e * g_t1_b20 * sp.diff(H_t1_b20.as_expr(), e)
    - quotient_t1_b20.as_expr() * H_t1_b20.as_expr()
    - remainder_t1_b20.as_expr(),
    "t=1,b2=0 Euclidean division",
)
E10 = 3483891 * b0**7 - 54675 * b0**4 * b1 + 1
zero(
    remainder_t1_b20.nth(10) - 4 * E10 / (81 * b0**6),
    "t=1,b2=0 top E10",
)
b1_solution = (3483891 * b0**7 + 1) / (54675 * b0**4)
X = sp.symbols("X")
poly_a = 1480774572985482 * X**2 + 867784104 * X - 53
poly_b = 579994041785274 * X**2 - 28717497 * X + 4
zero(E10.subs(b1, b1_solution), "solve t=1,b2=0 E10")
zero(
    remainder_t1_b20.nth(9).subs(b1, b1_solution)
    + 4 * poly_a.subs(X, b0**7) / (4100625 * b0**8),
    "t=1,b2=0 e9 polynomial",
)
zero(
    remainder_t1_b20.nth(8).subs(b1, b1_solution)
    + 4 * poly_b.subs(X, b0**7) / (110716875 * b0**10),
    "t=1,b2=0 e8 polynomial",
)
gate(
    sp.resultant(poly_a, poly_b, X)
    == 2408049135195443414554999563281250,
    "t=1,b2=0 coprime resultant",
)


# Repeated-root control for the one-way boundary divisibility criterion.
alpha, beta = sp.symbols("alpha beta")
boundary_g = (e - alpha) * (e - beta)
boundary_h = e**3 * (e - alpha)**2 * (e - beta)
gate(
    sp.rem(
        sp.Poly(sp.expand(e * boundary_g * sp.diff(boundary_h, e)), e),
        sp.Poly(boundary_h, e),
    ).is_zero,
    "repeated-root boundary divisibility control",
)


# Exact reconstruction from every off-boundary residual root.
r_rec = u / e
z_rec = 9 * G * u * (1 + u) / (e * K)
zero(z_rec**2 - K / (9 * G) + Q / (9 * G * e**2 * K**2),
     "Q reconstructs z square")
zero(
    surface.subs({r: r_rec, z: z_rec})
    - u * (1 + u) * Q / (e**3 * K**3),
    "Q reconstructs surface",
)
A_e_symbol = 2 * e + r_rec * D
C_e_rec = sp.factor(9 * G * z_rec**2 - K)
C_z_rec = sp.factor(3 * G * r_rec**2 - 3 * K * A_e_symbol)
zero(C_e_rec + Q / (e**2 * K**2), "Q kills e Hamiltonian")
zero(C_z_rec - 3 * P / e**2, "P kills z Hamiltonian")
zero(Q.subs(u, 0) - e**2, "exclude u=0")
zero(Q.subs(u, -1) + e**2, "exclude u=-1")
zero(Q.subs(u, -sp.Rational(1, 2)) + 729 * G**3 / 16,
     "exclude K=0")
zero(sp.LC(sp.Poly(Q, u)) + 729 * G**3,
     "Q leading coefficient excludes infinity")


semantic = {
    "carrier": "A=e^2-z/3+r*sum_(i=0)^3 b_i e^i; c=1; support>=3",
    "inheritance": "b3=0 is THM-3805; b3!=0 handled here without THM-3806",
    "generic": "t!=0,1; b2=0 killed at e15; b2!=0 splits A/B; only t=2 binomial seam",
    "drop": "t=1; b2!=0 degrees 15 and two-row contradiction; b2=0 degree 11 and coprime X=b0^7 quadratics",
    "criterion": "boundary-only roots imply H divides e*g*Hprime with multiplicity",
    "recovery": "r=u/e; z=9g*u*(u+1)/(e*(1+2u))",
    "open": "support<=2 assigned separately; degree>=4; mixed corrections; rational-pole mates",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3807-trinomial-and-full-cubic-r-repairs-have-critical-points")
print("surface=r^2e-z^3+r;c=1;field=algebraically_closed_characteristic_zero")
print("carrier=A=e^2-z/3+r*(b0+b1e+b2e^2+b3e^3);support_at_least_3")
print("generic=b3_notin_0,1;degree_H=17;branches=A_or_B;only_t2_boundary_is_binomial")
print("drop=b3=1;b2_nonzero_degree15;b2_zero_degree11;both_empty")
print("boundary_criterion=all_H_roots_in_V(e*g)_implies_H_divides_e*g*Hprime")
print("recovery=r=u/e;z=9g*u*(u+1)/(e*(1+2u))")
print("open=support_at_most_2_separate;degree_ge_4;mixed_corrections;rational_poles")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
