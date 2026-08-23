#!/usr/bin/env python3
"""Exact companion for THM-3723's constant C3 Cohn resonance."""

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
u, v, h = sp.symbols("u v h", nonzero=True)


def eplus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, value), (0, 1)))


def eminus(value: sp.Expr) -> sp.Matrix:
    return sp.Matrix(((1, 0), (value, 1)))


def curl(one_form: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(sp.diff(one_form[0], T) - sp.diff(one_form[1], X))


def jac(first: sp.Expr, second: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(first, x) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, x)
    )


def monomials(total_degree: int) -> list[sp.Expr]:
    return [
        X**i * T**j
        for i in range(total_degree + 1)
        for j in range(total_degree + 1 - i)
    ]


def linear_solution_set(
    expression: sp.Expr,
    unknowns: tuple[sp.Symbol, ...],
) -> sp.Set:
    equations = [
        coefficient
        for _, coefficient in sp.Poly(sp.expand(expression), X, T).terms()
    ]
    return sp.linsolve(equations, unknowns)


M0 = sp.Matrix(((4 * T**2, 2 * X * T - 1), (1 + 2 * X * T, X**2)))
R = eplus(v) * eminus(u)
N = sp.expand(M0 * R)
alpha = (sp.expand(N[0, 0]), sp.expand(N[0, 1]))
beta = (sp.expand(N[1, 0]), sp.expand(N[1, 1]))
gate(sp.expand(R.det()) == 1, "right word determinant")
gate(sp.expand(N.det()) == 1, "decorated determinant")

lower = (
    sp.expand(beta[0] + h * alpha[0]),
    sp.expand(beta[1] + h * alpha[1]),
)
top = (
    sp.expand(alpha[0] + h * beta[0]),
    sp.expand(alpha[1] + h * beta[1]),
)
lower_formula = 2 * (
    u * (h + v) * X + (4 * h * u * v + 3 * h - v) * T
)
top_formula = 2 * (
    u * (h * v + 1) * X + (-h * v + 4 * u * v + 3) * T
)
gate(sp.expand(curl(lower) - lower_formula) == 0, "lower exposed curl")
gate(sp.expand(curl(top) - top_formula) == 0, "top exposed curl")

# With arbitrary polynomial right parameters, retain the derivative terms and
# check both exact first-closure expressions directly.
U = sp.Function("U")(X, T)
V = sp.Function("V")(X, T)
H = sp.Function("H")(X, T)
N_polynomial = sp.expand(M0 * eplus(V) * eminus(U))
alpha_polynomial = (
    sp.expand(N_polynomial[0, 0]),
    sp.expand(N_polynomial[0, 1]),
)
beta_polynomial = (
    sp.expand(N_polynomial[1, 0]),
    sp.expand(N_polynomial[1, 1]),
)
lower_polynomial = curl(
    (
        sp.expand(beta_polynomial[0] + H * alpha_polynomial[0]),
        sp.expand(beta_polynomial[1] + H * alpha_polynomial[1]),
    )
)
top_polynomial = curl(
    (
        sp.expand(alpha_polynomial[0] + H * beta_polynomial[0]),
        sp.expand(alpha_polynomial[1] + H * beta_polynomial[1]),
    )
)
gate(
    sp.expand(lower_polynomial.subs({U: u, V: v, H: h}).doit() - lower_formula) == 0,
    "polynomial lower specializes to constant curl",
)
gate(
    sp.expand(top_polynomial.subs({U: u, V: v, H: h}).doit() - top_formula) == 0,
    "polynomial top specializes to constant curl",
)

# In the nondegenerate constant cell u,v != 0, the two coefficient systems
# have the unique common resonance uv=-1 with h=-v or h=u respectively.
gate(sp.expand(lower_formula.subs(h, -v).subs(u, -1 / v)) == 0,
     "lower C3 closure")
gate(sp.expand(top_formula.subs(h, u).subs(u, -1 / v)) == 0,
     "top C3 closure")
gate(
    sp.solve(
        (
            sp.expand(lower_formula.coeff(X) / 2),
            sp.expand(lower_formula.coeff(T) / 2),
        ),
        (h, u),
        dict=True,
    ) == [{h: -v, u: -1 / v}],
    "lower symbolic classification",
)
gate(
    sp.solve(
        (
            sp.expand(top_formula.coeff(X) / 2),
            sp.expand(top_formula.coeff(T) / 2),
        ),
        (h, u),
        dict=True,
    ) == [{h: -1 / v, u: -1 / v}],
    "top symbolic classification",
)

R_c3 = sp.simplify(R.subs(u, -1 / v))
gate(sp.simplify(R_c3**3 + sp.eye(2)) == sp.zeros(2), "projective order three")
gate(sp.simplify(R_c3**2 - R_c3 + sp.eye(2)) == sp.zeros(2),
     "C3 quadratic law")

# Both exposed orientations produce the same gradient line.  In the source
# coordinates below it is exactly the Broughton polynomial L+L^2S.
N_c3 = sp.simplify(N.subs(u, -1 / v))
alpha_c3 = (sp.expand(N_c3[0, 0]), sp.expand(N_c3[0, 1]))
beta_c3 = (sp.expand(N_c3[1, 0]), sp.expand(N_c3[1, 1]))
gamma = (
    sp.expand(beta_c3[0] - v * alpha_c3[0]),
    sp.expand(beta_c3[1] - v * alpha_c3[1]),
)
eta = (
    sp.expand(alpha_c3[0] - beta_c3[0] / v),
    sp.expand(alpha_c3[1] - beta_c3[1] / v),
)
L_source = 2 * v * T - X
M_source = X + v * T
S_source = -M_source / (3 * v)
Q = sp.expand(L_source + L_source**2 * S_source)
gate(jac(L_source, S_source, X, T) == 1, "resonance source coordinates")
gate(sp.expand(gamma[0] - sp.diff(Q, X)) == 0, "lower potential X")
gate(sp.expand(gamma[1] - sp.diff(Q, T)) == 0, "lower potential T")
gate(sp.expand(eta[0] + sp.diff(Q, X) / v) == 0, "top potential X")
gate(sp.expand(eta[1] + sp.diff(Q, T) / v) == 0, "top potential T")

inverse = sp.solve(
    (sp.Eq(L, L_source), sp.Eq(S, S_source)),
    (X, T),
    dict=True,
)[0]
alpha_debt = sp.factor(curl(alpha_c3).subs(inverse))
beta_debt = sp.factor(curl(beta_c3).subs(inverse))
gate(alpha_debt == 6 * S, "lower remaining Broughton debt")
gate(beta_debt == 6 * v * S, "top remaining Broughton debt")

# Finite Hamiltonian-cokernel guard in the normalized source coordinates.
B = L + L**2 * S
semantic_rows: list[str] = []
for bound in range(0, 13):
    broughton_basis = [
        L**i * S**j
        for i in range(bound + 1)
        for j in range(bound + 1 - i)
    ]
    coefficients = sp.symbols(f"f{bound}_0:{len(broughton_basis)}")
    candidate = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, broughton_basis)
    )
    equation = sp.Poly(jac(candidate, B, L, S) - 6 * S, L, S)
    system = [coefficient for _, coefficient in equation.terms()]
    gate(sp.linsolve(system, coefficients) == sp.EmptySet,
         "finite Broughton cokernel")
    semantic_rows.append(
        f"degree={bound}:" + hashlib.sha256(sp.srepr(equation.as_expr()).encode()).hexdigest()
    )

# Hostile rational constant grid: solve each exposed coefficient pair and
# compare with the exact uv=-1 classification.
values = (
    sp.Rational(-3), sp.Rational(-2), sp.Rational(-1),
    sp.Rational(-1, 2), sp.Rational(-1, 3),
    sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(1),
    sp.Rational(2), sp.Rational(3),
)
grid_cells = 0
for u0 in values:
    for v0 in values:
        lower_coefficients = (
            sp.expand(lower_formula.subs({u: u0, v: v0}) / 2).coeff(X),
            sp.expand(lower_formula.subs({u: u0, v: v0}) / 2).coeff(T),
        )
        top_coefficients = (
            sp.expand(top_formula.subs({u: u0, v: v0}) / 2).coeff(X),
            sp.expand(top_formula.subs({u: u0, v: v0}) / 2).coeff(T),
        )
        lower_solutions = sp.solve(lower_coefficients, h, dict=True)
        top_solutions = sp.solve(top_coefficients, h, dict=True)
        expected_lower = [{h: -v0}] if u0 * v0 == -1 else []
        expected_top = [{h: u0}] if u0 * v0 == -1 else []
        gate(lower_solutions == expected_lower, "lower hostile classification")
        gate(top_solutions == expected_top, "top hostile classification")
        grid_cells += 2

# Independent nonconstant hostile grid.  The universal proof uses leading
# forms and directional-derivative kernels; this bounded solve guards all
# three degree regimes and both exposed orientations.
probes = (
    sp.Integer(-1), sp.Integer(1), sp.Integer(2),
    X, T, X + T + 1, X**2 + T, X * T + 1,
)
nonconstant_cells = 0
basis = monomials(6)
for u_index, u_probe in enumerate(probes):
    for v_index, v_probe in enumerate(probes):
        if not ({X, T} & (u_probe.free_symbols | v_probe.free_symbols)):
            continue
        coefficients = sp.symbols(
            f"p_{u_index}_{v_index}_0:{len(basis)}"
        )
        h_probe = sum(
            coefficient * monomial
            for coefficient, monomial in zip(coefficients, basis)
        )
        substitutions = {U: u_probe, V: v_probe, H: h_probe}
        lower_probe = sp.expand(lower_polynomial.subs(substitutions).doit())
        top_probe = sp.expand(top_polynomial.subs(substitutions).doit())
        gate(linear_solution_set(lower_probe, coefficients) == sp.EmptySet,
             "nonconstant lower hostile gate")
        gate(linear_solution_set(top_probe, coefficients) == sp.EmptySet,
             "nonconstant top hostile gate")
        semantic_rows.append(
            f"nonconstant={u_index},{v_index}:"
            + hashlib.sha256(
                (sp.srepr(lower_probe) + sp.srepr(top_probe)).encode()
            ).hexdigest()
        )
        nonconstant_cells += 2

# At selected resonances, allowing a degree-six h introduces no homogeneous
# solutions: the exact solution remains the constant h=-v or h=u.
for v0 in (sp.Rational(-2), sp.Rational(-1), sp.Rational(1), sp.Rational(2)):
    u0 = -1 / v0
    coefficients = sp.symbols(f"r_{str(v0).replace('/', '_')}_0:{len(basis)}")
    h_probe = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, basis)
    )
    substitutions = {U: u0, V: v0, H: h_probe}
    lower_solution = linear_solution_set(
        sp.expand(lower_polynomial.subs(substitutions).doit()), coefficients
    )
    top_solution = linear_solution_set(
        sp.expand(top_polynomial.subs(substitutions).doit()), coefficients
    )
    expected_lower = tuple([-v0] + [sp.Integer(0)] * (len(basis) - 1))
    expected_top = tuple([u0] + [sp.Integer(0)] * (len(basis) - 1))
    gate(lower_solution == sp.FiniteSet(expected_lower),
         "resonant lower homogeneous kernel")
    gate(top_solution == sp.FiniteSet(expected_top),
         "resonant top homogeneous kernel")

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

semantic = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()
print("theorem=THM-3723-automorphic-Cohn-C3-complete-two-right-nonentry")
print("right_word=E_+(v)E_-(u);nonzero_polynomials")
print("nonconstant_gate=either_parameter_nonconstant_prevents_first_closure")
print("constant_resonance=uv=-1")
print("projective_law=R^3=-I;R^2-R+I=0")
print("closed_gradient=Q=L+L^2S;J(L,S)=1")
print("remaining_debts=6S,6vS;Broughton_cokernel_nonzero")
print(f"hostile_constant_cells={grid_cells}")
print(f"hostile_nonconstant_cells={nonconstant_cells};degree_h<=6")
print(f"semantic_sha256={semantic}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
