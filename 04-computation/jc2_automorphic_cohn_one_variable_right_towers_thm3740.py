#!/usr/bin/env python3
"""Exact companion for THM-3740's one-variable right-shear towers."""

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
A_symbol = sp.symbols("A_symbol", nonzero=True)
h_function = sp.Function("h_function")(X, T)
v_function = sp.Function("v_function")(X, T)
u_function = sp.Function("u_function")(X, T)


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


M0 = sp.Matrix(((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2)))

# Derive both full closure PDEs before specializing the right parameter to one
# variable.  Multiplication by a removes the diagonal square root.
a_symbol = sp.symbols("a_symbol", nonzero=True)
A_relation = {a_symbol**2: A_symbol}
D_symbol = sp.diag(a_symbol, 1 / a_symbol)
E_plus = sp.Matrix(((1, v_function), (0, 1)))
E_minus = sp.Matrix(((1, 0), (u_function, 1)))

N_plus = sp.expand(M0 * D_symbol * E_plus)
alpha_plus = (N_plus[0, 0], N_plus[0, 1])
beta_plus = (N_plus[1, 0], N_plus[1, 1])
lower_plus = (
    beta_plus[0] + h_function * alpha_plus[0],
    beta_plus[1] + h_function * alpha_plus[1],
)
L_A = (
    4 * A_symbol * T**2 * sp.diff(h_function, T)
    - (2 * X * T - 1) * sp.diff(h_function, X)
    + 2 * (4 * A_symbol - 1) * T * h_function
    + 2 * (A_symbol - 1) * X
)
lower_formula = L_A - A_symbol * sp.diff(
    v_function * (1 + 2 * X * T + 4 * T**2 * h_function), X
)
gate(sp.simplify(
    a_symbol * curl(lower_plus)
    - lower_formula.subs(A_symbol, a_symbol**2)
) == 0, "universal E-plus/lower closure PDE")

N_minus = sp.expand(M0 * D_symbol * E_minus)
alpha_minus = (N_minus[0, 0], N_minus[0, 1])
beta_minus = (N_minus[1, 0], N_minus[1, 1])
upper_minus = (
    alpha_minus[0] + h_function * beta_minus[0],
    alpha_minus[1] + h_function * beta_minus[1],
)
U_A = (
    A_symbol * (1 + 2 * X * T) * sp.diff(h_function, T)
    - X**2 * sp.diff(h_function, X)
    + 2 * (A_symbol - 1) * X * h_function
    + 2 * (4 * A_symbol - 1) * T
)
upper_formula = U_A + sp.diff(
    u_function * (2 * X * T - 1 + X**2 * h_function), T
)
gate(sp.simplify(
    a_symbol * curl(upper_minus)
    - upper_formula.subs(A_symbol, a_symbol**2)
) == 0, "universal E-minus/upper closure PDE")

# Crossed orientations on the one-variable right-data boundary.
v_of_T = sp.Function("v_of_T")(T)
u_of_X = sp.Function("u_of_X")(X)
N_plus_cross = sp.expand(
    M0 * D_symbol * sp.Matrix(((1, v_of_T), (0, 1)))
)
alpha_plus_cross = (N_plus_cross[0, 0], N_plus_cross[0, 1])
beta_plus_cross = (N_plus_cross[1, 0], N_plus_cross[1, 1])
upper_plus_cross = (
    alpha_plus_cross[0] + h_function * beta_plus_cross[0],
    alpha_plus_cross[1] + h_function * beta_plus_cross[1],
)
upper_cross_formula = U_A - A_symbol * v_of_T * (
    (1 + 2 * X * T) * sp.diff(h_function, X) + 2 * T * h_function
)
gate(sp.simplify(
    a_symbol * curl(upper_plus_cross)
    - upper_cross_formula.subs(A_symbol, a_symbol**2)
) == 0, "one-variable E-plus/upper crossed PDE")

N_minus_cross = sp.expand(
    M0 * D_symbol * sp.Matrix(((1, 0), (u_of_X, 1)))
)
alpha_minus_cross = (N_minus_cross[0, 0], N_minus_cross[0, 1])
beta_minus_cross = (N_minus_cross[1, 0], N_minus_cross[1, 1])
lower_minus_cross = (
    beta_minus_cross[0] + h_function * alpha_minus_cross[0],
    beta_minus_cross[1] + h_function * alpha_minus_cross[1],
)
lower_cross_formula = L_A + u_of_X * (
    (2 * X * T - 1) * sp.diff(h_function, T) + 2 * X * h_function
)
gate(sp.simplify(
    a_symbol * curl(lower_minus_cross)
    - lower_cross_formula.subs(A_symbol, a_symbol**2)
) == 0, "one-variable E-minus/lower crossed PDE")

semantic_rows: list[str] = []
right_profiles = (
    sp.Integer(0),
    sp.Integer(1),
    T + 2 * T**2,
    -3 + T**3 - 2 * T**5,
)
dual_profiles = tuple(profile.xreplace({T: X}) for profile in right_profiles)

for r in range(1, 9):
    # Lower tower after an arbitrary v(T).
    A = sp.Rational(r + 1, 2 * r)
    a = sp.sqrt(A)
    K = sp.expand((1 + 2 * Z / r) ** r)
    H = sp.cancel((K - 1 - 2 * Z) / (4 * Z**2))
    phi = sp.cancel(
        r * ((1 + 2 * Z / r) ** (r + 1) - 1)
        / (2 * (r + 1) * Z)
    )
    for index, v in enumerate(right_profiles):
        polynomial = sp.Poly(v, T)
        B = sp.Integer(0)
        for (degree,), coefficient in polynomial.terms():
            denominator = sp.Rational(degree + 1) + 1 / (2 * A)
            gate(denominator != 0,
                 f"lower transport denominator r={r},degree={degree}")
            B += coefficient * T ** (degree + 1) / denominator
        B = sp.expand(B)
        gate(sp.expand(sp.diff(B, T) + B / (2 * A * T) - v) == 0,
             f"lower transport ODE r={r},profile={index}")
        gate(sp.expand(B.subs(T, 0)) == 0,
             f"lower transport origin r={r},profile={index}")

        Y = sp.expand(X + B)
        h = sp.expand(B / (2 * T) + Y**2 * H.subs(Z, Y * T))
        Q = sp.expand(a * Y * phi.subs(Z, Y * T))
        N = sp.expand(M0 * sp.diag(a, 1 / a) * sp.Matrix(((1, v), (0, 1))))
        alpha = (N[0, 0], N[0, 1])
        beta = (N[1, 0], N[1, 1])
        row = (
            sp.expand(beta[0] + h * alpha[0]),
            sp.expand(beta[1] + h * alpha[1]),
        )
        gate(sp.simplify(N.det() - 1) == 0,
             f"lower right determinant r={r},profile={index}")
        gate(jac(Y, T) == 1,
             f"lower nonlinear source translation r={r},profile={index}")
        gate(sp.simplify(curl(row)) == 0,
             f"lower raw closure r={r},profile={index}")
        gate(sp.simplify(sp.diff(Q, X) - row[0]) == 0,
             f"lower Q_X r={r},profile={index}")
        gate(sp.simplify(sp.diff(Q, T) - row[1]) == 0,
             f"lower Q_T r={r},profile={index}")
        # The exact connection defect canceled by the transport ODE.
        diagonal_second = sp.expand(
            (Y**2 + Y**2 * H.subs(Z, Y * T) * (2 * Y * T - 1)) / a
        )
        unsheared_second = sp.expand(
            (X**2 + h * (2 * X * T - 1)) / a
        )
        gate(sp.simplify(
            unsheared_second - diagonal_second + B * K.subs(Z, Y * T) / (2 * a * T)
        ) == 0, f"lower connection defect r={r},profile={index}")

    # Upper dual tower after an arbitrary u(X).
    A_top = sp.Rational(r, 2 * (r + 1))
    a_top = sp.sqrt(A_top)
    K_top = sp.expand(-(1 - 2 * Z / r) ** r)
    H_top = sp.cancel((K_top + 1 - 2 * Z) / Z**2)
    phi_top = sp.cancel(
        -r * (1 - (1 - 2 * Z / r) ** (r + 1))
        / (2 * (r + 1) * Z)
    )
    for index, u in enumerate(dual_profiles):
        polynomial = sp.Poly(u, X)
        C_profile = sp.Integer(0)
        for (degree,), coefficient in polynomial.terms():
            denominator = sp.Rational(degree + 1) + 2 * A_top
            gate(denominator != 0,
                 f"upper transport denominator r={r},degree={degree}")
            C_profile += coefficient * X ** (degree + 1) / denominator
        C_profile = sp.expand(C_profile)
        gate(sp.expand(
            sp.diff(C_profile, X) + 2 * A_top * C_profile / X - u
        ) == 0, f"upper transport ODE r={r},profile={index}")
        gate(sp.expand(C_profile.subs(X, 0)) == 0,
             f"upper transport origin r={r},profile={index}")

        V = sp.expand(T + C_profile)
        h = sp.expand(2 * C_profile / X + V**2 * H_top.subs(Z, X * V))
        Q = sp.expand((1 / a_top) * V * phi_top.subs(Z, X * V))
        N = sp.expand(
            M0 * sp.diag(a_top, 1 / a_top)
            * sp.Matrix(((1, 0), (u, 1)))
        )
        alpha = (N[0, 0], N[0, 1])
        beta = (N[1, 0], N[1, 1])
        row = (
            sp.expand(alpha[0] + h * beta[0]),
            sp.expand(alpha[1] + h * beta[1]),
        )
        gate(sp.simplify(N.det() - 1) == 0,
             f"upper right determinant r={r},profile={index}")
        gate(jac(X, V) == 1,
             f"upper nonlinear source translation r={r},profile={index}")
        gate(sp.simplify(curl(row)) == 0,
             f"upper raw closure r={r},profile={index}")
        gate(sp.simplify(sp.diff(Q, X) - row[0]) == 0,
             f"upper Q_X r={r},profile={index}")
        gate(sp.simplify(sp.diff(Q, T) - row[1]) == 0,
             f"upper Q_T r={r},profile={index}")

    semantic_rows.append(
        f"r={r}:"
        + hashlib.sha256(
            "|".join(
                sp.srepr(item)
                for item in (A, H, phi, B, A_top, H_top, phi_top, C_profile)
            ).encode()
        ).hexdigest()
    )

# Hostile bounded controls for the sharp v_X and u_T gates.  The theorem's
# leading-degree proof is unbounded; these exact affine systems guard mixed
# low-degree and boundary cases.
mixed_v_profiles = (X, X + T, X**2 + 1, X * T + T**2)
mixed_u_profiles = (T, T + X, T**2 - 1, X * T + X**2)
hostile_systems = 0
for r in range(1, 4):
    for A, v in (
        (sp.Rational(r + 1, 2 * r), profile)
        for profile in mixed_v_profiles
    ):
        for bound in range(0, 5):
            basis = [
                X**i * T**j
                for i in range(bound + 1)
                for j in range(bound + 1 - i)
            ]
            coefficients = sp.symbols(f"lv_{r}_{hostile_systems}_0:{len(basis)}")
            h = sum(coefficient * monomial for coefficient, monomial in zip(coefficients, basis))
            equation = sp.Poly(
                (
                    4 * A * T**2 * sp.diff(h, T)
                    - (2 * X * T - 1) * sp.diff(h, X)
                    + 2 * (4 * A - 1) * T * h
                    + 2 * (A - 1) * X
                    - A * sp.diff(v * (1 + 2 * X * T + 4 * T**2 * h), X)
                ), X, T
            )
            system = [coefficient for _, coefficient in equation.terms()]
            gate(sp.linsolve(system, coefficients) == sp.EmptySet,
                 "hostile mixed-X right parameter")
            hostile_systems += 1

for r in range(1, 4):
    for A, u in (
        (sp.Rational(r, 2 * (r + 1)), profile)
        for profile in mixed_u_profiles
    ):
        for bound in range(0, 5):
            basis = [
                X**i * T**j
                for i in range(bound + 1)
                for j in range(bound + 1 - i)
            ]
            coefficients = sp.symbols(f"uu_{r}_{hostile_systems}_0:{len(basis)}")
            h = sum(coefficient * monomial for coefficient, monomial in zip(coefficients, basis))
            equation = sp.Poly(
                (
                    A * (1 + 2 * X * T) * sp.diff(h, T)
                    - X**2 * sp.diff(h, X)
                    + 2 * (A - 1) * X * h
                    + 2 * (4 * A - 1) * T
                    + sp.diff(u * (2 * X * T - 1 + X**2 * h), T)
                ), X, T
            )
            system = [coefficient for _, coefficient in equation.terms()]
            gate(sp.linsolve(system, coefficients) == sp.EmptySet,
                 "hostile mixed-T right parameter")
            hostile_systems += 1

# The transport maps have a one-dimensional cokernel at negative resonant A,
# but the closure equation itself has no survivor there.  These exact systems
# guard the exceptional branch, while the theorem uses the all-degree leading
# ODE to exclude it uniformly.
exceptional_systems = 0
for resonance in range(0, 4):
    A_lower_bad = -sp.Rational(1, 2 * (resonance + 1))
    v_bad = T**resonance
    A_upper_bad = -sp.Rational(resonance + 1, 2)
    u_bad = X**resonance
    for m in range(2, 8):
        for j in range(0, 8):
            gate(A_lower_bad != sp.Rational(m + 1, 2 * (j + 2)),
                 "lower exceptional leading-ODE separation")
            gate(A_upper_bad != sp.Rational(j + 2, 2 * (m + 1)),
                 "upper exceptional leading-ODE separation")
    for bound in range(0, 6):
        basis = [
            X**i * T**j
            for i in range(bound + 1)
            for j in range(bound + 1 - i)
        ]
        coefficients = sp.symbols(
            f"elb_{resonance}_{bound}_0:{len(basis)}"
        )
        h = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        lower_equation = sp.Poly(
            4 * A_lower_bad * T**2 * sp.diff(h, T)
            - (2 * X * T - 1) * sp.diff(h, X)
            + 2 * (4 * A_lower_bad - 1) * T * h
            + 2 * (A_lower_bad - 1) * X
            - A_lower_bad * sp.diff(
                v_bad * (1 + 2 * X * T + 4 * T**2 * h), X
            ), X, T
        )
        system = [coefficient for _, coefficient in lower_equation.terms()]
        gate(sp.linsolve(system, coefficients) == sp.EmptySet,
             "lower exceptional transport hostile")
        exceptional_systems += 1

        coefficients = sp.symbols(
            f"eub_{resonance}_{bound}_0:{len(basis)}"
        )
        h = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        upper_equation = sp.Poly(
            A_upper_bad * (1 + 2 * X * T) * sp.diff(h, T)
            - X**2 * sp.diff(h, X)
            + 2 * (A_upper_bad - 1) * X * h
            + 2 * (4 * A_upper_bad - 1) * T
            + sp.diff(
                u_bad * (2 * X * T - 1 + X**2 * h), T
            ), X, T
        )
        system = [coefficient for _, coefficient in upper_equation.terms()]
        gate(sp.linsolve(system, coefficients) == sp.EmptySet,
             "upper exceptional transport hostile")
        exceptional_systems += 1

# Exact crossed-orientation boundary controls.  Nonconstant one-variable
# right data leave only the depth-one h=0 sheets; constants add precisely the
# second THM-3726 sheet.
cross_controls = 0
for v in (T, 1 + T**2, -2 + T**3):
    A = sp.Rational(1, 4)
    a = sp.sqrt(A)
    N = sp.expand(M0 * sp.diag(a, 1 / a) * sp.Matrix(((1, v), (0, 1))))
    alpha = (N[0, 0], N[0, 1])
    gate(curl(alpha) == 0, "nonconstant E-plus upper depth-one cross")
    cross_controls += 1
for u in (X, 1 + X**2, -2 + X**3):
    A = sp.Integer(1)
    a = sp.sqrt(A)
    N = sp.expand(M0 * sp.diag(a, 1 / a) * sp.Matrix(((1, 0), (u, 1))))
    beta = (N[1, 0], N[1, 1])
    gate(curl(beta) == 0, "nonconstant E-minus lower depth-one cross")
    cross_controls += 1
for constant in (-3, -1, 1, 2):
    # E-plus upper: A=1/4,h=0 and A=1,h=3/v.
    for A, h in ((sp.Rational(1, 4), 0), (sp.Integer(1), sp.Rational(3, constant))):
        a = sp.sqrt(A)
        N = sp.expand(
            M0 * sp.diag(a, 1 / a)
            * sp.Matrix(((1, constant), (0, 1)))
        )
        alpha = (N[0, 0], N[0, 1])
        beta = (N[1, 0], N[1, 1])
        gate(curl((alpha[0] + h * beta[0], alpha[1] + h * beta[1])) == 0,
             "constant E-plus upper crossed sheets")
        cross_controls += 1
    # E-minus lower: A=1,h=0 and A=1/4,h=3/(4u).
    for A, h in ((sp.Integer(1), 0),
                 (sp.Rational(1, 4), sp.Rational(3, 4 * constant))):
        a = sp.sqrt(A)
        N = sp.expand(
            M0 * sp.diag(a, 1 / a)
            * sp.Matrix(((1, 0), (constant, 1)))
        )
        alpha = (N[0, 0], N[0, 1])
        beta = (N[1, 0], N[1, 1])
        gate(curl((beta[0] + h * alpha[0], beta[1] + h * alpha[1])) == 0,
             "constant E-minus lower crossed sheets")
        cross_controls += 1

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3740-automorphic-Cohn-one-variable-right-shear-towers")
print(
    "scope=arbitrary_right_parameters_in_compatible_exposures;"
    "all_one_variable_crossed_exposures"
)
print("sharp_gate=lower_forces_v_X=0;upper_forces_u_T=0")
print("transport=Y=X+B(T)_and_V=T+C(X)_with_unique_triangular_ODEs")
print("inheritance=source_automorphic_images_of_THM3734;no_polynomial_mates")
print(
    f"tested_depths=1..8;hostile_systems={hostile_systems};"
    f"exceptional_systems={exceptional_systems};cross_controls={cross_controls}"
)
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
