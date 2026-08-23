#!/usr/bin/env python3
"""Exact companion for THM-3758's rational-exact split-fibre family."""

from __future__ import annotations

import ast
import hashlib
from itertools import product
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


X, T, z, L = sp.symbols("X T z L")
a0, a1, beta, gamma, kappa = sp.symbols(
    "a0 a1 beta gamma kappa", nonzero=True
)


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def total_degree_basis(bound: int) -> list[sp.Expr]:
    return [
        X**i * T**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]


# Universal parameter identities in the independent (X,z) chart.
A = a0 + a1 * z
B = gamma + beta * a0 * z + beta * a1 * z**2 / 2
delta = 2 * a1 * gamma - beta * a0**2
gate(sp.expand(sp.diff(B, z) - beta * A) == 0, "B'=beta A")
gate(
    sp.expand(2 * B * sp.diff(A, z) - A * sp.diff(B, z) - delta) == 0,
    "smoothness Bezout constant",
)
gate(
    sp.factor(sp.resultant(A, B, z) - a1 * delta / 2) == 0,
    "A,B coprimality resultant",
)

Q_chart = sp.expand(X * A + X**2 * B)
Y_chart = sp.expand(A + 2 * X * B)
D = sp.expand(A**2 + 4 * B * L)
C = a1 + 2 * beta * L
gate(
    sp.expand(Y_chart**2 - (A**2 + 4 * B * Q_chart)) == 0,
    "quadratic generic-fibre identity",
)
gate(sp.expand(sp.diff(D, z) - 2 * A * C) == 0,
     "quadratic exactness derivative")
gate(
    sp.factor(sp.discriminant(D, z) + 8 * delta * L * C) == 0,
    "generic-fibre nonsquare discriminant",
)

# The exact primitive is easiest to verify after returning to z=XT.
z_source = X * T
A_source = A.subs(z, z_source)
B_source = B.subs(z, z_source)
Q = sp.expand(X * A_source + X**2 * B_source)
Y = sp.expand(A_source + 2 * X * B_source)
P0 = sp.cancel(
    -kappa / (2 * Q)
    * (z_source + Y / (a1 + 2 * beta * Q))
)
gate(sp.cancel(jac(P0, Q) - kappa) == 0,
     "universal rational Jacobian mate")

# The two zero-fibre components see opposite principal coefficients.
S = z + Y_chart / C
axis_value = sp.cancel(S.subs({X: 0, z: 0, L: 0}))
residual_value = sp.cancel(
    S.subs({X: -A / B, L: 0})
)
gate(sp.cancel(axis_value - a0 / a1) == 0,
     "axis principal coefficient")
gate(sp.cancel(residual_value + a0 / a1) == 0,
     "residual principal coefficient")
gate(sp.cancel(axis_value + residual_value) == 0,
     "opposite split-fibre principal parts")


# Frozen parameter controls include the normalized integral hostile and
# several rational members of the same all-degree family.
parameter_controls = (
    (1, -2, -1, 0),
    (1, 1, 1, 0),
    (1, 2, 1, 1),
    (2, -1, 1, 0),
    (-1, 2, -1, 1),
)
smooth_controls = 0
rational_controls = 0
for av0, av1, bv, gv in parameter_controls:
    deltav = 2 * av1 * gv - bv * av0**2
    gate(av0 != 0 and av1 != 0 and bv != 0 and deltav != 0,
         "admissible parameter control")
    qv = sp.expand(Q.subs({a0: av0, a1: av1, beta: bv, gamma: gv}))
    critical = sp.groebner(
        (sp.diff(qv, X), sp.diff(qv, T)), X, T, order="lex"
    )
    gate(critical.contains(sp.Integer(1)), "smooth parameter control")
    pv = P0.subs({a0: av0, a1: av1, beta: bv, gamma: gv, kappa: 1})
    gate(sp.cancel(jac(pv, qv) - 1) == 0,
         "rational-mate parameter control")
    gate(sp.factor(qv / X).has(X), "nontrivial residual zero-fibre component")
    smooth_controls += 1
    rational_controls += 1


# Sharp failures: delta=0 creates a torus critical point, while a0=0 loses
# the protected X-axis derivative.
boundary_delta = {
    a0: 1,
    a1: -2,
    beta: -1,
    gamma: sp.Rational(1, 4),
}
q_delta = sp.expand(Q.subs(boundary_delta))
gate(
    sp.diff(q_delta, X).subs({X: 2, T: sp.Rational(1, 2)}) == 0
    and sp.diff(q_delta, T).subs({X: 2, T: sp.Rational(1, 2)}) == 0,
    "delta-zero torus critical boundary",
)
q_axis = sp.expand(Q.subs({a0: 0, a1: 1, beta: 1, gamma: 1}))
gate(
    sp.diff(q_axis, X).subs({X: 0, T: 0}) == 0
    and sp.diff(q_axis, T).subs({X: 0, T: 0}) == 0,
    "a0-zero axis critical boundary",
)


# Full polynomial-mate systems for the normalized hostile and three nearby
# parameter points.  The proof is all-degree; this is an independent finite
# guard on the primitive-torsor conclusion.
mate_obstructions = 0
rank_gap_min = None
for profile_index, parameters in enumerate(parameter_controls[:4]):
    av0, av1, bv, gv = parameters
    qv = sp.expand(Q.subs({a0: av0, a1: av1, beta: bv, gamma: gv}))
    for bound in range(0, 13):
        basis = total_degree_basis(bound)
        coefficients = sp.symbols(
            f"p_{profile_index}_{bound}_0:{len(basis)}"
        )
        prospective = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        equation = sp.Poly(jac(prospective, qv) - 1, X, T)
        rows = [coefficient for _, coefficient in equation.terms()]
        matrix, rhs = sp.linear_eq_to_matrix(rows, coefficients)
        gap = matrix.row_join(rhs).rank() - matrix.rank()
        gate(gap == 1, "bounded polynomial-mate obstruction")
        rank_gap_min = gap if rank_gap_min is None else min(rank_gap_min, gap)
        mate_obstructions += 1


# The normalized integral witness and its exact source support.
normalized = sp.expand(Q.subs({a0: 1, a1: -2, beta: -1, gamma: 0}))
expected_normalized = X - 2 * X**2 * T - X**3 * T + X**4 * T**2
gate(sp.expand(normalized - expected_normalized) == 0,
     "normalized integral hostile")
support = tuple(exponent for exponent, _ in sp.Poly(normalized, X, T).terms())
gate(support == ((4, 2), (3, 1), (2, 1), (1, 0)),
     "normalized support")


semantic_rows = (
    "family:z=XT;A=a0+a1*z;B=gamma+beta*a0*z+(beta*a1/2)*z^2",
    "parameter:a0*a1*beta*delta_nonzero;delta=2*a1*gamma-beta*a0^2",
    "smooth:Bprime=beta*A;2*B*A'-A*B'=delta",
    "fibre:Y=A+2*X*B;Y^2=A^2+4*B*L;dY=(a1+2*beta*L)*A*dz/Y",
    "mate:P0=-c/(2L)*(z+Y/(a1+2*beta*L))",
    "nonregular:opposite_principal_parts_on_X=0_and_A+X*B=0",
    "normalized:X-2*X^2*T-X^3*T+X^4*T^2",
)
semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("theorem=THM-3758-quadratic-radial-carrier-rational-exact-split-fibre")
print("scope=algebraically_closed_characteristic_zero;four_parameter_family")
print("geometry=smooth_reducible_noncoordinate;delta_nonzero")
print("rational_mate=explicit_complete_primitive_torsor")
print("polynomial_mate=none_by_opposite_zero_fibre_principal_parts")
print("mechanism=quadratic_generic_fibre_exactness_then_vertical_split")
print(
    f"smooth_controls={smooth_controls};rational_controls={rational_controls};"
    f"mate_obstructions={mate_obstructions};rank_gap_min={rank_gap_min}"
)
print("boundaries=delta_zero_torus_critical;a0_zero_axis_critical")
print("normalized=X-2*X^2*T-X^3*T+X^4*T^2;total_degree=6")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
