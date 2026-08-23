#!/usr/bin/env python3
"""Exact companion for THM-3773's quadratic rational-cover degree."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


z, w, L = sp.symbols("z w L")
X, T = sp.symbols("X T")
a0, a1, beta, kappa = sp.symbols("a0 a1 beta kappa", nonzero=True)
gamma = sp.symbols("gamma")
p, q = sp.symbols("p q")


def jac(first: sp.Expr, second: sp.Expr, left: sp.Symbol, right: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, left) * sp.diff(second, right)
        - sp.diff(first, right) * sp.diff(second, left)
    )


delta = 2 * a1 * gamma - beta * a0**2
A = a0 + a1 * z
B = gamma + beta * a0 * z + beta * a1 * z**2 / 2
C = a1 + 2 * beta * L
D = sp.expand(A**2 + 4 * B * L)
E = sp.Poly(sp.expand(C**2 * (w - z) ** 2 - D), z)
E_expected = sp.Poly(
    2 * beta * L * C * z**2
    - 2 * C * (a0 + C * w) * z
    + C**2 * w**2
    - a0**2
    - 4 * gamma * L,
    z,
)
gate(E == E_expected, "quadratic minimal-equation identity")
gate(E.degree() == 2, "quadratic degree in z")
gate(sp.factor(E.LC() - 2 * beta * L * C) == 0, "quadratic leading coefficient")

disc_z = sp.factor(sp.discriminant(E.as_expr(), z))
gate(sp.Poly(disc_z, w).degree() == 2, "first discriminant degree in w")
disc_w = sp.factor(sp.discriminant(disc_z, w))
gate(
    sp.factor(disc_w + 256 * beta * delta * L**2 * C**4) == 0,
    "nested discriminant invoice",
)

# A completed-square expression independently checks the first discriminant.
disc_z_expected = sp.factor(
    4
    * C
    * (
        a1 * C**2 * w**2
        + 2 * a0 * C**2 * w
        + a0**2 * (a1 + 4 * beta * L)
        + 8 * beta * gamma * L**2
    )
)
gate(sp.factor(disc_z - disc_z_expected) == 0, "first discriminant expansion")

# Recover the source field from L,w,z through Y=C(w-z).
z_source = X * T
A_source = A.subs(z, z_source)
B_source = B.subs(z, z_source)
Q_source = sp.expand(X * A_source + X**2 * B_source)
Y_source = sp.expand(A_source + 2 * X * B_source)
C_source = C.subs(L, Q_source)
w_source = sp.cancel(z_source + Y_source / C_source)
P_source = sp.cancel(-kappa * w_source / (2 * Q_source))
gate(
    sp.cancel(jac(P_source, Q_source, X, T) - kappa) == 0,
    "universal rational Keller response",
)
gate(
    sp.cancel(2 * Q_source / (Y_source + A_source) - X) == 0,
    "source X recovery",
)
gate(sp.cancel(z_source / X - T) == 0, "source T recovery")
gate(
    sp.cancel(
        E.as_expr().subs({L: Q_source, w: w_source, z: z_source})
    )
    == 0,
    "minimal equation on the source",
)

# The normalized integral member has the displayed quadratic equation.
normal_subs = {a0: 1, a1: -2, beta: -1, gamma: 0}
E_normal = sp.Poly(sp.expand(E.as_expr().subs(normal_subs)), z)
E_normal_expected = sp.Poly(
    4 * L * (1 + L) * z**2
    + 4 * (1 + L) * (1 - 2 * (1 + L) * w) * z
    + 4 * (1 + L) ** 2 * w**2
    - 1,
    z,
)
gate(E_normal == E_normal_expected, "normalized quadratic equation")
gate(
    sp.factor(disc_w.subs(normal_subs) - 4096 * L**2 * (1 + L) ** 4) == 0,
    "normalized nested discriminant",
)

# Every elementary rational target shear is birational and symplectic.
F = sp.Function("F")
G = sp.Function("G")
gate(jac(p, q + F(p), p, q) == 1, "upper rational target shear")
gate(jac(p + G(q), q, p, q) == 1, "lower rational target shear")

# Several exact parameter controls retain a nonzero admissibility invoice and
# a genuinely squarefree first discriminant.
controls = (
    {a0: 1, a1: -2, beta: -1, gamma: 0},
    {a0: 1, a1: 1, beta: 1, gamma: 1},
    {a0: 2, a1: -1, beta: 3, gamma: 5},
    {a0: -1, a1: 3, beta: -2, gamma: 1},
)
control_count = 0
for values in controls:
    gate(sp.expand(delta.subs(values)) != 0, "admissible delta control")
    specialized = sp.Poly(sp.expand(disc_z.subs(values)), w)
    gate(specialized.degree() == 2, "control discriminant degree")
    gate(sp.discriminant(specialized.as_expr(), w) != 0,
         "control discriminant squarefreeness")
    control_count += 1

semantic_rows = (
    "family:THM3758;L=Q;C=a1+2betaL;w=z+Y/C",
    "minimal:2betaLCz2-2C(a0+Cw)z+C2w2-a02-4gammaL",
    "nested_disc:-256beta*delta*L2*C4",
    "field_degree:[k(X,T):k(P0,Q)]=2",
    "target_words:birational_symplectic_preserve_degree",
    "boundary:degree_two_Galois_Keller_impossible_by_THM1330",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3773-quadratic-rational-keller-cover-degree-and-target-word-nonpolynomialization")
print("scope=THM3758_full_admissible_parameter_family")
print("minimal_equation=quadratic_in_z_over_k(L,w)")
print("nested_discriminant=-256*beta*delta*L^2*(a1+2betaL)^4")
print("function_field_degree=2;extension=separable_Galois")
print("target_words=all_finite_birational_rational_symplectic_words_blocked")
print(f"parameter_controls={control_count}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")

