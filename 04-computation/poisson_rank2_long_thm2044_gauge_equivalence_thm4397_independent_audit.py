#!/usr/bin/env python3
"""Independent exact audit for the THM-4397 gauge equivalence.

This script is self-contained.  It reconstructs both polynomial Poisson maps
from their displayed formulas and does not import the primary THM-4397
certificate or any repository mathematics.
"""

from __future__ import annotations

import sympy as sp


def check(condition: bool, label: str) -> None:
    """Optimization-safe replacement for ``assert``."""
    if not bool(condition):
        raise RuntimeError(f"FAIL: {label}")


x, q, p, z, e = sp.symbols("x q p z e")
variables = (x, q, p, z)
s = x * q


def poisson(f: sp.Expr, h: sp.Expr) -> sp.Expr:
    """Bracket with {p,x}={z,q}=1."""
    return sp.expand(
        sp.diff(f, p) * sp.diff(h, x)
        - sp.diff(f, x) * sp.diff(h, p)
        + sp.diff(f, z) * sp.diff(h, q)
        - sp.diff(f, q) * sp.diff(h, z)
    )


def terms(f: sp.Expr) -> int:
    return len(sp.Poly(sp.expand(f), *variables, domain=sp.QQ).terms())


def total_degree(f: sp.Expr) -> int:
    return sp.Poly(sp.expand(f), *variables, domain=sp.QQ).total_degree()


def value(polynomials: tuple[sp.Expr, ...], point: tuple[sp.Rational, ...]):
    substitution = dict(zip(variables, point))
    return tuple(sp.factor(f.subs(substitution)) for f in polynomials)


# Shared determinant-one three-dimensional core, with e independent.
y_e = q - x * e / 3
v_e = x * y_e
w_e = 1 + v_e
sigma_e = sp.expand(
    y_e + 3 * x * w_e**2 * e + 3 * x * y_e**2 * (4 + 3 * v_e)
)
pi_half_e = sp.expand(
    (w_e**3 * e + y_e**2 * w_e * (4 + 3 * v_e)) / 2
)

# Common elementary source data.
a = 1 - 3 * s
c = (1 + 3 * s) / 2
L = 3 * x**2 * p + 2 * a * z
D0 = c * p - 3 * q**2 * z
R = sp.expand(x * (2 - 3 * s))

# The earlier THM-2044 gauge.
G = 252 * s**3 + 1008 * s**2 + 1379 * s + 659
g = -q**2 * G / 140
ell = sp.expand(L + g)
repo_connection_B_e = sp.expand(
    2 * e**4 * x**6 * (3 * s - 2)
    + e**3 * x**4 * (-90 * s**2 - 30 * s + 55)
    + e**2 * x**2 * (540 * s**3 + 720 * s**2 - 120 * s - 270)
    + e * (-1620 * s**4 - 3780 * s**3 - 1890 * s**2 + 810 * s + 540)
    + q**2 * (2430 * s**3 + 8100 * s**2 + 8640 * s + 2430)
)
H_repo_e = sp.expand(-e * repo_connection_B_e / 1620)
T_repo = sp.expand(sigma_e.subs(e, ell))
S_repo = sp.expand(pi_half_e.subs(e, ell))
D_repo = sp.expand(D0 + H_repo_e.subs(e, ell))
phi_repo = (R, T_repo, D_repo, S_repo)

# Long's arXiv:2608.23777v1 gauge.
beta = sp.expand(L - 9 * q**2)
H_paper_e = sp.expand(
    y_e**4 * (18 * v_e**2 + 78 * v_e + 125) / 20
    + 3 * e * y_e**2 * (v_e**3 + 5 * v_e**2 + 10 * v_e - 5) / 10
    - e**2 * (9 * v_e + 2) / 6
    - x**2 * e**3 / 6
)
S_paper = sp.expand(sigma_e.subs(e, beta))
T_paper = sp.expand(-pi_half_e.subs(e, beta))
D_paper = sp.expand(D0 + H_paper_e.subs(e, beta))
phi_paper = (R, T_paper, D_paper, S_paper)

# Negative generating-function orientation.  A_U carries paper-source points
# to THM-2044-source points, and Phi_paper = K o Phi_repo o A_U.
U = sp.expand(-q**3 * (63 * x**2 * q**2 + 318 * x * q + 601) / 840)
U_x = sp.expand(sp.diff(U, x))
U_q = sp.expand(sp.diff(U, q))
P_A = sp.expand(p + U_x)
Z_A = sp.expand(z + U_q)

print("THM-4397 INDEPENDENT GAUGE-EQUIVALENCE AUDIT")
print("convention: {p,x}={z,q}=1; A_U sends (p,z) to (p+U_x,z+U_q)")
print("U =", sp.factor(U))
print("U_x =", sp.factor(U_x))
print("U_q =", sp.factor(U_q))

# Classical and exact-Weyl checks for the momentum translation.
canonical_checks = (
    poisson(P_A, x) == 1,
    poisson(Z_A, q) == 1,
    poisson(x, q) == 0,
    poisson(P_A, q) == 0,
    poisson(Z_A, x) == 0,
    poisson(P_A, Z_A) == 0,
)
check(all(canonical_checks), "A_U classical canonical relations")
check(sp.expand(sp.diff(U_q, x) - sp.diff(U_x, q)) == 0, "closed gradient")
print("A_U classical canonical relations: PASS")

# In the Weyl convention [p,x]=[z,q]=1, multiplication by U_x,U_q gives
# [p+U_x,z+U_q]=d_x U_q-d_q U_x; all other relations are immediate.
weyl_relation_residuals = (
    sp.Integer(1) - 1,
    sp.Integer(1) - 1,
    sp.Integer(0),
    sp.Integer(0),
    sp.diff(U_q, x) - sp.diff(U_x, q),
)
check(all(sp.expand(residual) == 0 for residual in weyl_relation_residuals),
      "exact Weyl lift")
print("exact Weyl lift [p+U_x,z+U_q]=0 and mixed canonical relations: PASS")

# Two independent adapted-coordinate identities.  They avoid a blind giant
# expansion by comparing the L and D0 increments before substitution.
L_increment = sp.expand(3 * x**2 * U_x + 2 * a * U_q)
D0_increment = sp.expand(c * U_x - 3 * q**2 * U_q)
H_difference_e = sp.expand(H_paper_e - H_repo_e)
expected_H_difference = sp.expand(q**4 * (18 * s**2 + 78 * s + 125) / 20)
check(sp.expand(H_difference_e - expected_H_difference) == 0,
      "connection difference independent of e")
check(sp.expand(D0_increment - expected_H_difference) == 0,
      "D0 translation equals connection difference")
ell_after_A = sp.expand(ell + L_increment)
check(sp.expand(ell_after_A - beta) == 0, "ell after A_U equals beta")
D_repo_after_A = sp.expand(
    D0 + D0_increment + H_repo_e.subs(e, ell_after_A)
)
check(sp.expand(D_repo_after_A - D_paper) == 0,
      "D_repo after A_U equals D_paper")
print("ell o A_U = beta: PASS")
print("D_repo o A_U = D_paper: PASS")
print("H_paper(e)-H_repo(e) =", sp.factor(H_difference_e))

# The common core and the target quarter-turn K(r,t,d,s)=(r,-s,d,t)
# now prove the complete right-left identity without coefficient guessing.
T_repo_after_A = sp.expand(sigma_e.subs(e, ell_after_A))
S_repo_after_A = sp.expand(pi_half_e.subs(e, ell_after_A))
K_phi_repo_after_A = (
    R,
    sp.expand(-S_repo_after_A),
    D_repo_after_A,
    T_repo_after_A,
)
check(all(sp.expand(left - right) == 0
          for left, right in zip(K_phi_repo_after_A, phi_paper)),
      "Phi_paper = K o Phi_repo o A_U")
omega = sp.Matrix([
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [-1, 0, 0, 0],
    [0, -1, 0, 0],
])
jacobian_K = sp.Matrix([
    [1, 0, 0, 0],
    [0, 0, 0, -1],
    [0, 0, 1, 0],
    [0, 1, 0, 0],
])
check(jacobian_K.T * omega * jacobian_K == omega,
      "target quarter-turn K is symplectic")
print("K(r,t,d,s)=(r,-s,d,t) is symplectic: PASS")
print("Phi_paper = K o Phi_repo o A_U: PASS")

# Direct source-independent checks on the paper representative.
paper_relations = (
    (D_paper, R, 1),
    (S_paper, T_paper, 1),
    (R, S_paper, 0),
    (R, T_paper, 0),
    (D_paper, S_paper, 0),
    (D_paper, T_paper, 0),
)
check(all(poisson(left, right) == expected
          for left, right, expected in paper_relations),
      "six paper Poisson identities")
print("paper six Poisson identities by direct expansion: PASS")
print("paper term counts (R,T,D,S):", tuple(terms(f) for f in phi_paper))
print("paper total degrees (R,T,D,S):", tuple(total_degree(f) for f in phi_paper))
print("repo term counts (R,T,D,S):", tuple(terms(f) for f in phi_repo))

# Independent determinant and exact reduced fibre calculation in the three
# algebraically independent core variables.
X, Y, E = sp.symbols("X Y E")
V = X * Y
W = 1 + V
R_core = sp.expand(2 * X - 3 * X**2 * Y - X**3 * E)
S_core = sp.expand(Y + 3 * X * W**2 * E + 3 * X * Y**2 * (4 + 3 * V))
T_core = sp.expand(-(W**3 * E + Y**2 * W * (4 + 3 * V)) / 2)
core_det = sp.factor(sp.Matrix([R_core, T_core, S_core]).jacobian([X, Y, E]).det())
check(core_det == 1, "core determinant")
groebner = sp.groebner(
    [R_core, S_core, T_core - sp.Rational(1, 8)],
    E, Y, X,
    order="lex",
    domain=sp.QQ,
)
actual_basis = [sp.expand(poly.as_expr()) for poly in groebner.polys]
expected_basis = [
    E - sp.Rational(27, 4) * X**2 + sp.Rational(1, 4),
    Y + sp.Rational(3, 2) * X,
    X**3 - X,
]
check(actual_basis == expected_basis, "reduced three-point Groebner basis")
check(sp.gcd(X**3 - X, sp.diff(X**3 - X, X)) == 1,
      "reduced fibre eliminant")
print("core determinant: 1 PASS")
print("core reduced Groebner basis:")
for polynomial in actual_basis:
    print(" ", sp.factor(polynomial))

paper_points = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), -sp.Rational(1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(247, 96), -sp.Rational(89, 64)),
    (-sp.Rational(1), -sp.Rational(2, 3), sp.Rational(247, 96), -sp.Rational(89, 64)),
)
repo_points_expected = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), -sp.Rational(1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(224839, 90720), -sp.Rational(173417, 60480)),
    (-sp.Rational(1), -sp.Rational(2, 3), sp.Rational(224839, 90720), -sp.Rational(173417, 60480)),
)

transported_points = []
for point in paper_points:
    xx, qq, pp, zz = point
    base = {x: xx, q: qq}
    transported_points.append((
        xx,
        qq,
        sp.factor(pp + U_x.subs(base)),
        sp.factor(zz + U_q.subs(base)),
    ))
transported_points = tuple(transported_points)
check(transported_points == repo_points_expected, "three-point source transport")

paper_target = (sp.Rational(0), sp.Rational(1, 8), sp.Rational(0), sp.Rational(0))
repo_target = (sp.Rational(0), sp.Rational(0), sp.Rational(0), -sp.Rational(1, 8))
check(all(value(phi_paper, point) == paper_target for point in paper_points),
      "paper common fibre")
check(all(value(phi_repo, point) == repo_target for point in transported_points),
      "repo common fibre")
K_repo_target = (repo_target[0], -repo_target[3], repo_target[2], repo_target[1])
check(K_repo_target == paper_target, "target transport")
print("three paper-source points transported by A_U:")
for source, target in zip(paper_points, transported_points):
    print(" ", source, "->", target)
print("common target transport:", repo_target, "->", paper_target, "PASS")
print("exact fibre is reduced of cardinality three and transports bijectively: PASS")
print("THM-4397 INDEPENDENT AUDIT: PASS")
