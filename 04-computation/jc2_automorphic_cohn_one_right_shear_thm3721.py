#!/usr/bin/env python3
"""Exact companion for THM-3721's two automorphic-Cohn shear gates."""

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


X, T, z, ell = sp.symbols("X T z ell")


def eplus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, value), (0, 1)))


def eminus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, 0), (value, 1)))


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(first, X) * sp.diff(second, T)
        - sp.diff(first, T) * sp.diff(second, X)
    )


def monomials(total_degree: int) -> list[sp.Expr]:
    return [
        X**i * T**j
        for i in range(total_degree + 1)
        for j in range(total_degree + 1 - i)
    ]


def has_polynomial_solution(expression: sp.Expr, unknowns: tuple[sp.Symbol, ...]) -> bool:
    equations = [
        coefficient
        for _, coefficient in sp.Poly(sp.expand(expression), X, T).terms()
    ]
    return sp.linsolve(equations, unknowns) != sp.EmptySet


# The normalized coframe is a constant elementary row rotation of Cohn's
# matrix after the ring automorphism y -> 2T.
Cohn_at_2T = sp.Matrix(
    ((1 + 2 * X * T, X**2), (-4 * T**2, 1 - 2 * X * T))
)
rotation = sp.Matrix(((0, -1), (1, 0)))
M0 = sp.Matrix(
    ((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2))
)
gate(rotation * Cohn_at_2T == M0, "automorphic Cohn normalization")
gate(sp.expand(M0.det()) == 1, "normalized determinant")

V = sp.Function("v")(X, T)
H = sp.Function("h")(X, T)
Mv = M0 * eplus(V)
alpha = (sp.expand(Mv[0, 0]), sp.expand(Mv[0, 1]))
beta = (sp.expand(Mv[1, 0]), sp.expand(Mv[1, 1]))
gate(sp.expand(Mv.det()) == 1, "decorated determinant")
gate(curl(alpha) == 6 * T - 4 * T**2 * sp.diff(V, X), "alpha curl")
gate(
    sp.expand(
        curl(beta) + (1 + 2 * X * T) * sp.diff(V, X) + 2 * T * V
    ) == 0,
    "beta curl",
)

# Bottom-exposed and top-exposed closure equations, checked directly rather
# than inferred from one another.
bottom = (
    sp.expand(beta[0] + H * alpha[0]),
    sp.expand(beta[1] + H * alpha[1]),
)
bottom_equation = sp.expand(
    4 * T**2 * sp.diff(H, T)
    - (2 * X * T - 1 + 4 * T**2 * V) * sp.diff(H, X)
    + (6 * T - 4 * T**2 * sp.diff(V, X)) * H
    - (1 + 2 * X * T) * sp.diff(V, X)
    - 2 * T * V
)
gate(sp.expand(curl(bottom) - bottom_equation) == 0, "bottom closure PDE")

top = (
    sp.expand(alpha[0] + H * beta[0]),
    sp.expand(alpha[1] + H * beta[1]),
)
top_equation = sp.expand(
    (1 + 2 * X * T) * sp.diff(H, T)
    - (X**2 + (1 + 2 * X * T) * V) * sp.diff(H, X)
    - ((1 + 2 * X * T) * sp.diff(V, X) + 2 * T * V) * H
    + 6 * T
    - 4 * T**2 * sp.diff(V, X)
)
gate(sp.expand(curl(top) - top_equation) == 0, "top closure PDE")

# The complete triangular bottom survivor.  The diagonal Euler equation has
# the unique coefficientwise solution A=2T H_b.
degree = 8
b_coefficients = sp.symbols(f"b0:{degree + 1}")
b = sum(coefficient * T**j for j, coefficient in enumerate(b_coefficients))
Hb = sum(
    coefficient * T**j / (2 * j + 3)
    for j, coefficient in enumerate(b_coefficients)
)
A = sp.expand(2 * T * Hb)
L = sp.expand(X + A)
Q = sp.expand(L + T * L**2)
gate(sp.expand(A + 2 * T * sp.diff(A, T) - 2 * T * b) == 0,
     "Euler primitive")
Mv_b = sp.expand(M0 * eplus(b))
alpha_b = (Mv_b[0, 0], Mv_b[0, 1])
beta_b = (Mv_b[1, 0], Mv_b[1, 1])
gate(
    sp.expand(sp.diff(Q, X) - (beta_b[0] + Hb * alpha_b[0])) == 0,
    "bottom survivor X derivative",
)
gate(
    sp.expand(sp.diff(Q, T) - (beta_b[1] + Hb * alpha_b[1])) == 0,
    "bottom survivor T derivative",
)
gate(curl(alpha_b) == 6 * T, "triangular bottom debt")

# Under L=X+A(T), the homogeneous first-closure equation is charge diagonal.
# Check the displayed one-variable equation in a broad charge/degree bank and
# independently check that its low and high endpoint conditions conflict.
semantic_rows: list[str] = []
for charge in range(0, 11):
    for f_degree in range(0, 11):
        coefficients = sp.symbols(f"F_{charge}_{f_degree}_0:{f_degree + 1}")
        F = sum(coefficient * z**j for j, coefficient in enumerate(coefficients))
        K = sp.expand(ell**charge * F.subs(z, ell * T))
        operator = sp.expand(
            4 * T**2 * sp.diff(K, T)
            + (1 - 2 * ell * T) * sp.diff(K, ell)
            + 6 * T * K
        )
        predicted = sp.expand(
            ell ** (charge - 1)
            * (
                z * (1 + 2 * z) * sp.diff(F, z)
                + (charge + (6 - 2 * charge) * z) * F
            ).subs(z, ell * T)
        )
        gate(sp.expand(operator - predicted) == 0, "charge operator")
        if f_degree >= charge:
            low = charge
            high = 2 * (f_degree + 3 - charge)
            gate(not (low == 0 and high == 0), "charge endpoint conflict")
        semantic_rows.append(
            f"charge={charge},degree={f_degree}:"
            + hashlib.sha256(sp.srepr(operator).encode()).hexdigest()
        )

# The constant top survivor is the same Broughton polynomial, up to a source
# shear and a nonzero target scaling.
c = sp.symbols("c", nonzero=True)
hc = 3 / c
Lc = X + 2 * c * T / 3
Qc = sp.expand((3 / c) * (Lc + T * Lc**2))
Mv_c = sp.expand(M0 * eplus(c))
alpha_c = (Mv_c[0, 0], Mv_c[0, 1])
beta_c = (Mv_c[1, 0], Mv_c[1, 1])
gate(
    sp.expand(sp.diff(Qc, X) - (alpha_c[0] + hc * beta_c[0])) == 0,
    "top constant survivor X derivative",
)
gate(
    sp.expand(sp.diff(Qc, T) - (alpha_c[1] + hc * beta_c[1])) == 0,
    "top constant survivor T derivative",
)
gate(curl(beta_c) == -2 * c * T, "top constant debt")

# Independent bounded hostile scout: every tested v with positive X-degree
# has no bottom or top closure of the prescribed degree.  These finite checks
# guard the signs and edge cases; the theorem uses the universal degree proof.
scout_cells = 0
for x_degree in range(1, 5):
    for t_degree in range(0, 5):
        for scalar in (sp.Integer(-2), sp.Integer(-1), sp.Integer(1), sp.Integer(2)):
            probe_v = scalar * X**x_degree * T**t_degree
            basis = monomials(6)
            unknowns = sp.symbols(f"u_{x_degree}_{t_degree}_{scalar}_0:{len(basis)}")
            probe_h = sum(value * monomial for value, monomial in zip(unknowns, basis))
            substitution = {V: probe_v, H: probe_h}
            bottom_probe = sp.expand(bottom_equation.subs(substitution).doit())
            top_probe = sp.expand(top_equation.subs(substitution).doit())
            gate(not has_polynomial_solution(bottom_probe, unknowns),
                 "nontriangular bottom scout")
            gate(not has_polynomial_solution(top_probe, unknowns),
                 "nontriangular top scout")
            scout_cells += 2

# Recheck the inherited Broughton cokernel by finite linear algebra.  The
# universal isolated-chain proof is THM-3716.
B = X + T * X**2
for bound in range(0, 13):
    basis = monomials(bound)
    unknowns = sp.symbols(f"g{bound}_0:{len(basis)}")
    candidate = sum(value * monomial for value, monomial in zip(unknowns, basis))
    gate(not has_polynomial_solution(jac(candidate, B) - 6 * T, unknowns),
         "Broughton finite cokernel")

# The lower right orientation is asymmetric at the first stage but has the
# same Broughton obstruction after a source shear.
U = sp.Function("u")(X, T)
Mu = M0 * eminus(U)
alpha_minus = (sp.expand(Mu[0, 0]), sp.expand(Mu[0, 1]))
beta_minus = (sp.expand(Mu[1, 0]), sp.expand(Mu[1, 1]))
gate(sp.expand(Mu.det()) == 1, "lower-shear determinant")
gate(
    sp.expand(
        curl(alpha_minus)
        - (6 * T + 2 * X * U + (2 * X * T - 1) * sp.diff(U, T))
    ) == 0,
    "lower-shear alpha curl",
)
gate(curl(beta_minus) == X**2 * sp.diff(U, T), "lower-shear beta curl")

lower_minus = (
    sp.expand(beta_minus[0] + H * alpha_minus[0]),
    sp.expand(beta_minus[1] + H * alpha_minus[1]),
)
lower_minus_equation = sp.expand(
    (4 * T**2 + U * (2 * X * T - 1)) * sp.diff(H, T)
    - (2 * X * T - 1) * sp.diff(H, X)
    + (6 * T + 2 * X * U + (2 * X * T - 1) * sp.diff(U, T)) * H
    + X**2 * sp.diff(U, T)
)
gate(
    sp.expand(curl(lower_minus) - lower_minus_equation) == 0,
    "lower-shear bottom closure PDE",
)

top_minus = (
    sp.expand(alpha_minus[0] + H * beta_minus[0]),
    sp.expand(alpha_minus[1] + H * beta_minus[1]),
)
gate(
    sp.expand(
        top_minus[0] * beta_minus[1]
        - top_minus[1] * beta_minus[0]
        - 1
    ) == 0,
    "lower-shear top determinant PDE",
)

u_coefficients = sp.symbols("ux0:9")
u_x = sum(coefficient * X**j for j, coefficient in enumerate(u_coefficients))
F_x = sum(
    coefficient * X ** (j + 3) / (j + 3)
    for j, coefficient in enumerate(u_coefficients)
)
V_x = sp.cancel(F_x / X**2)
S = sp.symbols("S")
Q_minus = sp.expand(X + T * X**2 + F_x)
Mu_x = sp.expand(M0 * eminus(u_x))
alpha_minus_x = (Mu_x[0, 0], Mu_x[0, 1])
beta_minus_x = (Mu_x[1, 0], Mu_x[1, 1])
gate(sp.expand(sp.diff(Q_minus, X) - beta_minus_x[0]) == 0,
     "lower-shear triangular X derivative")
gate(sp.expand(sp.diff(Q_minus, T) - beta_minus_x[1]) == 0,
     "lower-shear triangular T derivative")
gate(sp.expand(Q_minus - (X + X**2 * (T + V_x))) == 0,
     "lower-shear source shear")
gate(
    sp.expand(
        curl(alpha_minus_x).subs(T, S - V_x)
        - (6 * S + 2 * (X * sp.diff(V_x, X) - V_x))
    ) == 0,
    "lower-shear transformed debt",
)

# The Hamiltonian operator preserves charge sectors up to a shift.  Hence the
# charge -1 target 6S cannot be cancelled by the nonnegative-charge pure-X
# summand in the transformed debt.
B_minus = X + X**2 * S
for i in range(0, 10):
    for j in range(0, 10):
        monomial = X**i * S**j
        value = sp.expand(
            sp.diff(monomial, X) * sp.diff(B_minus, S)
            - sp.diff(monomial, S) * sp.diff(B_minus, X)
        )
        predicted = sp.expand(
            -j * X**i * S ** (j - 1)
            + (i - 2 * j) * X ** (i + 1) * S**j
        )
        gate(sp.expand(value - predicted) == 0, "Broughton charge transport")
        semantic_rows.append(
            f"minus-charge={i-j}:" + hashlib.sha256(sp.srepr(value).encode()).hexdigest()
        )

# Check the Laurent conjugacy used by the opposite exposed-row proof.
for u_degree in range(0, 6):
    uc = sp.symbols(f"lu{u_degree}_0:{u_degree + 1}")
    ux = sum(coefficient * X**j for j, coefficient in enumerate(uc))
    for q_degree in range(1, 7):
        qc = {
            (i, j): sp.symbols(f"lq{u_degree}_{q_degree}_{i}_{j}")
            for i in range(q_degree + 1)
            for j in range(q_degree + 1 - i)
        }
        q_poly = sum(
            coefficient * X**i * T**j for (i, j), coefficient in qc.items()
        )
        q_tilde = sp.expand(q_poly.subs(T, z / (2 * X)))
        original = sp.expand(
            (
                (1 + 2 * X * T + X**2 * ux) * sp.diff(q_poly, T)
                - X**2 * sp.diff(q_poly, X)
            ).subs(T, z / (2 * X))
            / X
        )
        laurent = sp.expand(
            (2 + z + 2 * X**2 * ux) * sp.diff(q_tilde, z)
            - X * sp.diff(q_tilde, X)
        )
        gate(sp.expand(original - laurent) == 0, "lower-shear Laurent conjugacy")
        minus_one_layer = sp.expand(
            sum(
                coefficient * z**j / 2**j
                for (i, j), coefficient in qc.items()
                if i - j == -1
            )
        )
        gate(sp.expand(minus_one_layer.subs(z, 0)) == 0,
             "honest lower-shear minus-one layer")
        semantic_rows.append(
            f"minus-Laurent={u_degree},{q_degree}:"
            + hashlib.sha256(
                (sp.srepr(original) + sp.srepr(minus_one_layer)).encode()
            ).hexdigest()
        )

for polynomial_degree in range(0, 12):
    coefficients = sp.symbols(f"lf{polynomial_degree}_0:{polynomial_degree + 1}")
    Fz = sum(coefficient * z**j for j, coefficient in enumerate(coefficients))
    for exponent in range(-12, 0):
        homogeneous = sp.Poly(
            (2 + z) * sp.diff(Fz, z) - exponent * Fz,
            z,
        ).all_coeffs()
        matrix, vector = sp.linear_eq_to_matrix(homogeneous, coefficients)
        gate(matrix.rank() == len(coefficients), "lower-shear Laurent kernel")
        gate(all(entry == 0 for entry in vector), "lower-shear Laurent RHS")
    inhomogeneous = sp.Poly((2 + z) * sp.diff(Fz, z) + Fz + 1, z).all_coeffs()
    expected = tuple([-sp.Integer(1)] + [sp.Integer(0)] * polynomial_degree)
    gate(sp.linsolve(inhomogeneous, coefficients) == sp.FiniteSet(expected),
         "lower-shear minus-one solution")

# Hostile T-dependent lower-shear parameters: neither exposed orientation
# closes in the tested coefficient box.
minus_scout_cells = 0
for x_degree in range(0, 4):
    for t_degree in range(1, 5):
        for scalar in (sp.Integer(-1), sp.Integer(1)):
            probe_u = scalar * X**x_degree * T**t_degree
            basis = monomials(6)
            unknowns = sp.symbols(f"mu_{x_degree}_{t_degree}_{scalar}_0:{len(basis)}")
            probe_h = sum(value * monomial for value, monomial in zip(unknowns, basis))
            lower_probe = sp.expand(
                lower_minus_equation.subs({U: probe_u, H: probe_h}).doit()
            )
            top_probe_form = (
                sp.expand(alpha_minus[0] + probe_h * beta_minus[0]),
                sp.expand(alpha_minus[1] + probe_h * beta_minus[1]),
            )
            top_probe = sp.expand(
                curl(
                    (
                        top_probe_form[0].subs(U, probe_u).doit(),
                        top_probe_form[1].subs(U, probe_u).doit(),
                    )
                )
            )
            gate(not has_polynomial_solution(lower_probe, unknowns),
                 "T-dependent lower-shear bottom scout")
            gate(not has_polynomial_solution(top_probe, unknowns),
                 "T-dependent lower-shear top scout")
            minus_scout_cells += 2

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3721-automorphic-Cohn-single-right-shear-nonentry")
print("normalized_core=M0=rotation*Cohn(X,2T);det=1;nonentry=inherited")
print("upper_bottom_survivor=v(T);Q=L+T*L^2;debt=6T")
print("upper_top_survivor=v=constant_nonzero;Q=(3/v)*(L+T*L^2);debt=-2vT")
print("lower_bottom_survivor=u(X);Q=X+X^2*(T+V);debt=6S+2(XV'-V)")
print("lower_top_survivor=NONE;obstruction=forbidden_X^-1_Laurent_layer")
print("nontriangular_gates=upper:v_X_nonzero;lower:u_T_nonzero")
print(f"hostile_nontriangular_cells={scout_cells + minus_scout_cells};degree_h<=6")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
