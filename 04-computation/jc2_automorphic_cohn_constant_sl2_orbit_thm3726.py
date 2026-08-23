#!/usr/bin/env python3
"""Exact companion for THM-3726's constant SL2 orbit classification."""

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


X, T, L, S = sp.symbols("X T L S")
a, b, c, d, h = sp.symbols("a b c d h")


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(
    first: sp.Expr,
    second: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, x)
    )


M0 = sp.Matrix(((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2)))
R = sp.Matrix(((a, b), (c, d)))
N = sp.expand(M0 * R)
alpha = (sp.expand(N[0, 0]), sp.expand(N[0, 1]))
beta = (sp.expand(N[1, 0]), sp.expand(N[1, 1]))
lower = (
    sp.expand(beta[0] + h * alpha[0]),
    sp.expand(beta[1] + h * alpha[1]),
)
top = (
    sp.expand(alpha[0] + h * beta[0]),
    sp.expand(alpha[1] + h * beta[1]),
)
lower_formula = 2 * (
    (a + c * h - d) * X + ((4 * a - d) * h - b) * T
)
top_formula = 2 * (
    ((a - d) * h + c) * X + (4 * a - d - b * h) * T
)
gate(sp.expand(N.det() - (a * d - b * c)) == 0, "constant orbit determinant")
gate(sp.expand(curl(lower) - lower_formula) == 0, "lower constant closure")
gate(sp.expand(curl(top) - top_formula) == 0, "top constant closure")

# The compatibility resultant is exactly the square discriminant minus the
# determinant.  On SL2 it forces epsilon=2a-d to be +1 or -1.
lower_compatibility = sp.expand(
    c * b - (4 * a - d) * (d - a)
)
top_compatibility = sp.expand(
    (a - d) * (4 * a - d) + b * c
)
discriminant_error = sp.expand((2 * a - d) ** 2 - (a * d - b * c))
gate(sp.expand(lower_compatibility - discriminant_error) == 0,
     "lower discriminant resultant")
gate(sp.expand(top_compatibility - discriminant_error) == 0,
     "top discriminant resultant")

semantic_rows: list[str] = []

# Universal lower parametrization.  The closure equations are built in and
# determinant one is equivalent to epsilon^2=1.
epsilon, C, H = sp.symbols("epsilon C H")
A_lower = epsilon + C * H
D_lower = 2 * A_lower - epsilon
B_lower = (2 * A_lower + epsilon) * H
R_lower = sp.Matrix(((A_lower, B_lower), (C, D_lower)))
N_lower = sp.expand(M0 * R_lower)
alpha_lower = (N_lower[0, 0], N_lower[0, 1])
beta_lower = (N_lower[1, 0], N_lower[1, 1])
L_lower = epsilon * (X + 2 * H * T)
S_lower = (C * X + (2 * A_lower + epsilon) * T) / 3
Q_lower = sp.expand(L_lower + L_lower**2 * S_lower)
gate(sp.expand(R_lower.det() - epsilon**2) == 0,
     "lower parametrized determinant")
gate(curl((beta_lower[0] + H * alpha_lower[0],
           beta_lower[1] + H * alpha_lower[1])) == 0,
     "lower parametrized closure")
for sign in (-1, 1):
    substitution = {epsilon: sign}
    gate(jac(L_lower, S_lower, X, T).subs(substitution) == 1,
         "lower Broughton source coordinates")
    gate(
        sp.expand(
            sp.diff(Q_lower, X)
            - (beta_lower[0] + H * alpha_lower[0])
        ).subs(substitution) == 0,
        "lower Broughton X derivative",
    )
    gate(
        sp.expand(
            sp.diff(Q_lower, T)
            - (beta_lower[1] + H * alpha_lower[1])
        ).subs(substitution) == 0,
        "lower Broughton T derivative",
    )
    inverse = sp.solve(
        (
            sp.Eq(L, L_lower.subs(substitution)),
            sp.Eq(S, S_lower.subs(substitution)),
        ),
        (X, T),
        dict=True,
    )[0]
    debt = sp.factor(curl(alpha_lower).subs(substitution).subs(inverse))
    gate(debt == 6 * S, "lower universal debt")
    semantic_rows.append(
        f"lower-sign={sign}:"
        + hashlib.sha256(
            (sp.srepr(Q_lower.subs(substitution)) + sp.srepr(debt)).encode()
        ).hexdigest()
    )

# Universal top parametrization.
B, H_top = sp.symbols("B H_top")
A_top = (B * H_top - epsilon) / 2
D_top = 2 * A_top - epsilon
C_top = (A_top - epsilon) * H_top
R_top = sp.Matrix(((A_top, B), (C_top, D_top)))
N_top = sp.expand(M0 * R_top)
alpha_top = (N_top[0, 0], N_top[0, 1])
beta_top = (N_top[1, 0], N_top[1, 1])
L_top = epsilon * (H_top * X + 2 * T)
S_top = ((A_top - epsilon) * X + B * T) / 3
Q_top = sp.expand(L_top + L_top**2 * S_top)
gate(sp.expand(R_top.det() - epsilon**2) == 0,
     "top parametrized determinant")
gate(curl((alpha_top[0] + H_top * beta_top[0],
           alpha_top[1] + H_top * beta_top[1])) == 0,
     "top parametrized closure")
for sign in (-1, 1):
    substitution = {epsilon: sign}
    gate(jac(L_top, S_top, X, T).subs(substitution) == 1,
         "top Broughton source coordinates")
    gate(
        sp.expand(
            sp.diff(Q_top, X)
            - (alpha_top[0] + H_top * beta_top[0])
        ).subs(substitution) == 0,
        "top Broughton X derivative",
    )
    gate(
        sp.expand(
            sp.diff(Q_top, T)
            - (alpha_top[1] + H_top * beta_top[1])
        ).subs(substitution) == 0,
        "top Broughton T derivative",
    )
    inverse = sp.solve(
        (
            sp.Eq(L, L_top.subs(substitution)),
            sp.Eq(S, S_top.subs(substitution)),
        ),
        (X, T),
        dict=True,
    )[0]
    debt = sp.factor(curl(beta_top).subs(substitution).subs(inverse))
    gate(debt == -6 * S, "top universal debt")
    semantic_rows.append(
        f"top-sign={sign}:"
        + hashlib.sha256(
            (sp.srepr(Q_top.subs(substitution)) + sp.srepr(debt)).encode()
        ).hexdigest()
    )

# Hostile integer SL2 grid checks the raw coefficient systems, including all
# triangular and zero-entry boundaries in the box.
integer_matrices = 0
closed_rows = 0
for a0 in range(-4, 5):
    for b0 in range(-4, 5):
        for c0 in range(-4, 5):
            for d0 in range(-4, 5):
                if a0 * d0 - b0 * c0 != 1:
                    continue
                integer_matrices += 1
                lower_equations = (
                    a0 + c0 * h - d0,
                    (4 * a0 - d0) * h - b0,
                )
                top_equations = (
                    (a0 - d0) * h + c0,
                    4 * a0 - d0 - b0 * h,
                )
                lower_solutions = sp.solve(lower_equations, h, dict=True)
                top_solutions = sp.solve(top_equations, h, dict=True)
                for solutions in (lower_solutions, top_solutions):
                    if solutions:
                        closed_rows += 1
                        gate((2 * a0 - d0) ** 2 == 1,
                             "hostile closure discriminant")
                if (2 * a0 - d0) ** 2 != 1:
                    gate(not lower_solutions and not top_solutions,
                         "hostile off-discriminant gate")

# In the normalized coordinates, every surviving constant combination has
# the inherited Broughton target.  Recheck its finite cokernel.
Broughton = L + L**2 * S
for bound in range(0, 13):
    basis = [
        L**i * S**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]
    coefficients = sp.symbols(f"f{bound}_0:{len(basis)}")
    candidate = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, basis)
    )
    equation = sp.Poly(jac(candidate, Broughton, L, S) - 6 * S, L, S)
    system = [coefficient for _, coefficient in equation.terms()]
    gate(sp.linsolve(system, coefficients) == sp.EmptySet,
         "finite Broughton cokernel")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3726-automorphic-Cohn-constant-SL2-orbit-classification")
print("scope=arbitrary_R_in_SL2(k);constant_exposed_left_parameter")
print("closure_locus=(2a-d)^2=1_with_typed_linear_compatibility")
print("normal_form=Q=L+L^2S;J(L,S)=1")
print("remaining_debts=lower:6S;top:-6S;Broughton_nonentry")
print(f"hostile_integer_matrices={integer_matrices};closed_rows={closed_rows}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
