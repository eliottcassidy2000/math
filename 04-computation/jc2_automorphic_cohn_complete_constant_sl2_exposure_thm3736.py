#!/usr/bin/env python3
"""Exact companion for THM-3736's complete constant-right Cohn orbit."""

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


X, T, Z = sp.symbols("X T Z")
# Keep the general SL2 entries assumption-free: a or d may vanish before the
# valuation reduction.  Reciprocals appear only in the triangular branch,
# where det(R)=1 itself forces a,d nonzero.
a, b, c, d = sp.symbols("a b c d")
h_function = sp.Function("h_function")(X, T)


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


M0 = sp.Matrix(((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2)))
R = sp.Matrix(((a, b), (c, d)))
N = sp.expand(M0 * R)
alpha = (N[0, 0], N[0, 1])
beta = (N[1, 0], N[1, 1])
gate(sp.expand(N.det() - (a * d - b * c)) == 0,
     "general constant-right determinant")

# The two full closure equations.  No degree or support assumption is made.
lower = (
    beta[0] + h_function * alpha[0],
    beta[1] + h_function * alpha[1],
)
upper = (
    alpha[0] + h_function * beta[0],
    alpha[1] + h_function * beta[1],
)
lower_pde = (
    (4 * a * T**2 + c * (2 * X * T - 1)) * sp.diff(h_function, T)
    - (4 * b * T**2 + d * (2 * X * T - 1)) * sp.diff(h_function, X)
    + 2 * (c * X + (4 * a - d) * T) * h_function
    + 2 * ((a - d) * X - b * T)
)
upper_pde = (
    (a * (1 + 2 * X * T) + c * X**2) * sp.diff(h_function, T)
    - (b * (1 + 2 * X * T) + d * X**2) * sp.diff(h_function, X)
    + 2 * ((a - d) * X - b * T) * h_function
    + 2 * (c * X + (4 * a - d) * T)
)
gate(sp.expand(curl(lower) - lower_pde) == 0,
     "universal lower closure PDE")
gate(sp.expand(curl(upper) - upper_pde) == 0,
     "universal upper closure PDE")

# At the highest total degree, both orientations reduce to the same linear
# vector field eigen-equation.  The divisors T and X expose c and b,
# respectively, by a valuation descent.
h_top = sp.Function("h_top")(X, T)
w_h = (
    (2 * a * T + c * X) * sp.diff(h_top, T)
    - (2 * b * T + d * X) * sp.diff(h_top, X)
)
lower_top = 2 * T * w_h + 2 * (c * X + (4 * a - d) * T) * h_top
upper_top = X * w_h + 2 * ((a - d) * X - b * T) * h_top
F_lower = 2 * T * h_top
F_upper = 2 * X * h_top


def w(expression: sp.Expr) -> sp.Expr:
    return sp.expand(
        (2 * a * T + c * X) * sp.diff(expression, T)
        - (2 * b * T + d * X) * sp.diff(expression, X)
    )


gate(sp.expand(w(F_lower) + (2 * a - d) * F_lower - lower_top) == 0,
     "lower top-degree eigenvector identity")
gate(sp.expand(w(F_upper) + (2 * a - d) * F_upper - 2 * upper_top) == 0,
     "upper top-degree eigenvector identity")

# Exact triangular gauge identities after the valuation descent.  For the
# lower row c=0; for the upper row b=0.  The exceptional denominator is
# handled by the Jordan argument in the proof and by finite hostile checks
# below.
D = 2 * a + 1 / a
lambda_lower = 2 * b / D
kappa_lower = b / D
gate(sp.simplify(
    4 * a * lambda_lower - 4 * b + 2 * (1 / a) * lambda_lower
) == 0, "lower source-shear derivative coefficient")
gate(sp.simplify(
    (4 * a - 1 / a) * kappa_lower
    - (a - 1 / a) * lambda_lower - b
) == 0, "lower source-shear forcing coefficient")

mu_upper = c / D
kappa_upper = 2 * c / D
gate(sp.simplify(c - (1 / a) * mu_upper - 2 * a * mu_upper) == 0,
     "upper source-shear derivative coefficient")
gate(sp.simplify(
    (a - 1 / a) * kappa_upper + c
    - (4 * a - 1 / a) * mu_upper
) == 0, "upper source-shear forcing coefficient")

semantic_rows: list[str] = []

# Constructive positive controls for both triangular-gauge towers.  These
# check the raw matrix rows rather than only the reduced profile ODE.
for r in range(1, 9):
    A_lower = sp.Rational(r + 1, 2 * r)
    a_lower = sp.sqrt(A_lower)
    d_lower = 1 / a_lower
    K_lower = sp.expand((1 + 2 * Z / r) ** r)
    H_lower = sp.cancel((K_lower - 1 - 2 * Z) / (4 * Z**2))
    phi_lower = sp.cancel(
        r * ((1 + 2 * Z / r) ** (r + 1) - 1)
        / (2 * (r + 1) * Z)
    )
    for b_value in (-3, -1, 0, 2):
        R_lower = sp.Matrix(((a_lower, b_value), (0, d_lower)))
        N_lower = sp.expand(M0 * R_lower)
        alpha_lower = (N_lower[0, 0], N_lower[0, 1])
        beta_lower = (N_lower[1, 0], N_lower[1, 1])
        lam = sp.simplify(2 * b_value / (2 * a_lower + d_lower))
        kap = sp.simplify(b_value / (2 * a_lower + d_lower))
        Y = sp.expand(X + lam * T)
        h_lower = sp.expand(kap + Y**2 * H_lower.subs(Z, Y * T))
        Q_lower = sp.expand(a_lower * Y * phi_lower.subs(Z, Y * T))
        row_lower = (
            sp.expand(beta_lower[0] + h_lower * alpha_lower[0]),
            sp.expand(beta_lower[1] + h_lower * alpha_lower[1]),
        )
        gate(sp.simplify(R_lower.det() - 1) == 0,
             f"lower triangular determinant r={r},b={b_value}")
        gate(jac(Y, T) == 1,
             f"lower source shear r={r},b={b_value}")
        gate(sp.simplify(curl(row_lower)) == 0,
             f"lower triangular closure r={r},b={b_value}")
        gate(sp.simplify(sp.diff(Q_lower, X) - row_lower[0]) == 0,
             f"lower triangular Q_X r={r},b={b_value}")
        gate(sp.simplify(sp.diff(Q_lower, T) - row_lower[1]) == 0,
             f"lower triangular Q_T r={r},b={b_value}")

    A_upper = sp.Rational(r, 2 * (r + 1))
    a_upper = sp.sqrt(A_upper)
    d_upper = 1 / a_upper
    K_upper = sp.expand(-(1 - 2 * Z / r) ** r)
    H_upper = sp.cancel((K_upper + 1 - 2 * Z) / Z**2)
    phi_upper = sp.cancel(
        -r * (1 - (1 - 2 * Z / r) ** (r + 1))
        / (2 * (r + 1) * Z)
    )
    for c_value in (-2, 0, 1, 3):
        R_upper = sp.Matrix(((a_upper, 0), (c_value, d_upper)))
        N_upper = sp.expand(M0 * R_upper)
        alpha_upper = (N_upper[0, 0], N_upper[0, 1])
        beta_upper = (N_upper[1, 0], N_upper[1, 1])
        mu = sp.simplify(c_value / (2 * a_upper + d_upper))
        kap = sp.simplify(2 * c_value / (2 * a_upper + d_upper))
        V = sp.expand(T + mu * X)
        h_upper = sp.expand(kap + V**2 * H_upper.subs(Z, X * V))
        Q_upper = sp.expand(d_upper * V * phi_upper.subs(Z, X * V))
        row_upper = (
            sp.expand(alpha_upper[0] + h_upper * beta_upper[0]),
            sp.expand(alpha_upper[1] + h_upper * beta_upper[1]),
        )
        gate(sp.simplify(R_upper.det() - 1) == 0,
             f"upper triangular determinant r={r},c={c_value}")
        gate(jac(X, V) == 1,
             f"upper source shear r={r},c={c_value}")
        gate(sp.simplify(curl(row_upper)) == 0,
             f"upper triangular closure r={r},c={c_value}")
        gate(sp.simplify(sp.diff(Q_upper, X) - row_upper[0]) == 0,
             f"upper triangular Q_X r={r},c={c_value}")
        gate(sp.simplify(sp.diff(Q_upper, T) - row_upper[1]) == 0,
             f"upper triangular Q_T r={r},c={c_value}")

    semantic_rows.append(
        f"r={r}:"
        + hashlib.sha256(
            "|".join(
                sp.srepr(item)
                for item in (
                    A_lower, H_lower, phi_lower,
                    A_upper, H_upper, phi_upper,
                )
            ).encode()
        ).hexdigest()
    )

# Hostile finite controls for the valuation descent.  In each exact integer
# SL2 matrix, the top homogeneous eigen-equation has no nonzero solution with
# the required divisor whenever the forbidden off-triangular entry is nonzero.
integer_matrices: list[tuple[int, int, int, int]] = []
for av in range(-3, 4):
    for bv in range(-3, 4):
        for cv in range(-3, 4):
            for dv in range(-3, 4):
                if av * dv - bv * cv == 1:
                    integer_matrices.append((av, bv, cv, dv))


def only_zero_linear_solution(expression: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    polynomial = sp.Poly(sp.expand(expression), X, T)
    equations = [coefficient for _, coefficient in polynomial.terms()]
    matrix, _ = sp.linear_eq_to_matrix(equations, variables)
    return matrix.rank() == len(variables)


valuation_cases = 0
for av, bv, cv, dv in integer_matrices:
    substitution = {a: av, b: bv, c: cv, d: dv}
    for m in range(1, 6):
        degree = m + 1
        if cv != 0:
            coefficients = sp.symbols(f"lf_{av}_{bv}_{cv}_{dv}_{m}_1:{degree + 1}")
            F = sum(
                coefficients[i - 1] * X ** (degree - i) * T**i
                for i in range(1, degree + 1)
            )
            equation = (w(F) + (2 * a - d) * F).subs(substitution)
            gate(only_zero_linear_solution(equation, coefficients),
                 "lower forbidden-c valuation hostile")
            valuation_cases += 1
        if bv != 0:
            coefficients = sp.symbols(f"uf_{av}_{bv}_{cv}_{dv}_{m}_1:{degree + 1}")
            F = sum(
                coefficients[i - 1] * X**i * T ** (degree - i)
                for i in range(1, degree + 1)
            )
            equation = (w(F) + (2 * a - d) * F).subs(substitution)
            gate(only_zero_linear_solution(equation, coefficients),
                 "upper forbidden-b valuation hostile")
            valuation_cases += 1

# The exceptional gauge denominator 2a+d=0 gives a scalar plus a nilpotent
# operator on every homogeneous degree.  Exact matrices at a^2=-1/2 guard
# both signs and both triangular orientations.
I = sp.I
for a_jordan in (I / sp.sqrt(2), -I / sp.sqrt(2)):
    d_jordan = -2 * a_jordan
    gate(sp.simplify(a_jordan * d_jordan - 1) == 0,
         "Jordan determinant boundary")
    for m in range(1, 8):
        degree = m + 1
        for b_jordan in (-2, 0, 3):
            coefficients = sp.symbols(f"jl_{m}_{b_jordan}_0:{degree + 1}")
            F = sum(
                coefficients[i] * X ** (degree - i) * T**i
                for i in range(degree + 1)
            )
            jordan_w = (
                2 * a_jordan * (T * sp.diff(F, T) + X * sp.diff(F, X))
                - 2 * b_jordan * T * sp.diff(F, X)
            )
            equation = sp.expand(jordan_w + 4 * a_jordan * F)
            gate(only_zero_linear_solution(equation, coefficients),
                 "lower Jordan hostile boundary")
        for c_jordan in (-3, 0, 1):
            coefficients = sp.symbols(f"ju_{m}_{c_jordan}_0:{degree + 1}")
            F = sum(
                coefficients[i] * X ** (degree - i) * T**i
                for i in range(degree + 1)
            )
            jordan_w = (
                2 * a_jordan * (T * sp.diff(F, T) + X * sp.diff(F, X))
                + c_jordan * X * sp.diff(F, T)
            )
            equation = sp.expand(jordan_w + 4 * a_jordan * F)
            gate(only_zero_linear_solution(equation, coefficients),
                 "upper Jordan hostile boundary")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3736-automorphic-Cohn-complete-constant-SL2-exposure")
print("scope=arbitrary_constant_R_in_SL2;arbitrary_polynomial_exposed_parameter")
print("reduction=constant_Broughton_locus_or_nonconstant_triangular_binomial_gauge")
print("nonconstant_lower=forces_c=0;nonconstant_upper=forces_b=0")
print("exceptional_2a+d=0=Jordan_nilpotent_nonclosure")
print("mate=none_in_every_classified_case")
print(f"hostile_integer_sl2_matrices={len(integer_matrices)}")
print(f"valuation_cases={valuation_cases}")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
