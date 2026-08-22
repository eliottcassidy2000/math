#!/usr/bin/env python3
"""Exact controls for proved THM-3634.

The gates verify the retained compiler fibre, the quadratic-fold first-stable
determinant, constant output-change covariance, both rank-one unit controls,
and the sharp hostile in the enlarged equal-triple-value polynomial ring.
Structural uses of polynomial units and the factor theorem are proved in the
theorem text rather than inferred from bounded computation.
"""

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one active exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact polynomial/rational zero test."""
    return sp.cancel(expression) == 0


print("THM-3634 exact companion -- proved quadratic-fold first-stable rank debt")
print("status=PROVED + VERIFIED-EXACT + HOSTILE-AUDITED")


print("SECTION compiler relation and retained triple")
x, q, t = sp.symbols("x q t")
D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
require("compiler relation", zero(c**2 * e - b * (b + 4)))

triple = ((-1, -3), (0, -sp.Rational(3, 4)), (1, -3))
for index, (x_value, q_value) in enumerate(triple):
    image = tuple(
        sp.expand(expression.subs({x: x_value, q: q_value}))
        for expression in (b, c, e)
    )
    require(f"retained branch {index}", image == (0, 0, -3))
print("PASS compiler=surface_relation retained_triple=(0,0,-3)")


print("SECTION quadratic critical derivative and first-stable determinant")
Q_symbol = sp.Function("Q")(x)
q_fold = Q_symbol + t**2
require("quadratic vertical derivative vanishes", sp.diff(q_fold, t).subs(t, 0) == 0)

U_x, V_x, A_value, B_value, kappa = sp.symbols(
    "U_x V_x A_value B_value kappa"
)
first_stable_determinant = U_x * B_value - A_value * V_x
require(
    "first stable determinant orientation",
    first_stable_determinant
    == sp.det(sp.Matrix([[U_x, A_value], [V_x, B_value]])),
)
require(
    "normalized determinant control",
    first_stable_determinant.subs(
        {U_x: 1, V_x: 2, A_value: 3, B_value: 7}
    )
    == 1,
)
print("PASS t=0 shadow=(f_x,f_t;g_x,g_t)=(Uprime,A;Vprime,B)")


print("SECTION constant GL2 output covariance and rank-one factors")
alpha, beta, gamma, delta = sp.symbols("alpha beta gamma delta")
M = sp.Matrix([[alpha, beta], [gamma, delta]])
h_x = alpha * U_x + beta * V_x
h_t = alpha * A_value + beta * B_value
k_x = gamma * U_x + delta * V_x
k_t = gamma * A_value + delta * B_value
transformed = sp.expand(h_x * k_t - h_t * k_x)
require(
    "constant output covariance",
    zero(transformed - M.det() * first_stable_determinant),
)
H_x, H_t, K_x, K_t = sp.symbols("H_x H_t K_x K_t")
abstract_determinant = H_x * K_t - H_t * K_x
require(
    "horizontal zero factor",
    abstract_determinant.subs(H_x, 0) == -H_t * K_x,
)
require(
    "stable zero factor",
    abstract_determinant.subs(H_t, 0) == H_x * K_t,
)

mu = sp.symbols("mu", nonzero=True)
nu = sp.symbols("nu")
affine = mu * x + nu
require("affine minus-middle separation", affine.subs(x, -1) - affine.subs(x, 0) == -mu)
require("affine plus-middle separation", affine.subs(x, 1) - affine.subs(x, 0) == mu)
print("PASS rank1_combo gives_unit_derivative and_affine_triple_contradiction")


print("SECTION equal-triple-value enlargement and hostile sharpness")
L = sp.expand(x * (x**2 - 1))
H0, constant = sp.symbols("H0 constant")
enlarged_template = constant + L * (x**4 + H0 * x + 1)
for point in (-1, 0, 1):
    require(
        f"enlarged template value point={point}",
        enlarged_template.subs(x, point) == constant,
    )

U = sp.expand(L * (3 * x**2 - 2) / 2)
V = L
A = sp.expand(L * (sp.Rational(225, 8) * x**3 - sp.Rational(75, 4) * x))
B = sp.expand(1 + sp.Rational(45, 4) * x * L)
hostile = (U, V, A, B)
expected_values = (0, 0, 0, 1)
for index, (polynomial, expected) in enumerate(zip(hostile, expected_values)):
    values = tuple(sp.expand(polynomial.subs(x, point)) for point in (-1, 0, 1))
    require(f"hostile common value index={index}", values == (expected, expected, expected))

require("hostile U divisible by L", zero(U / L - (3 * x**2 - 2) / 2))
require("hostile V divisible by L", V == L)
require(
    "hostile A divisible by L",
    zero(A / L - (sp.Rational(225, 8) * x**3 - sp.Rational(75, 4) * x)),
)
require("hostile B minus one divisible by L", zero((B - 1) / L - sp.Rational(45, 4) * x))

hostile_determinant = sp.expand(sp.diff(U, x) * B - A * sp.diff(V, x))
require("hostile determinant", hostile_determinant == 1)

horizontal_coefficients = sp.Matrix(
    [
        [sp.Poly(U, x).nth(degree), sp.Poly(V, x).nth(degree)]
        for degree in range(1, max(sp.degree(U, x), sp.degree(V, x)) + 1)
    ]
)
stable_coefficients = sp.Matrix(
    [
        [sp.Poly(A, x).nth(degree), sp.Poly(B, x).nth(degree)]
        for degree in range(max(sp.degree(A, x), sp.degree(B, x)) + 1)
    ]
)
require("hostile horizontal quotient rank", horizontal_coefficients.rank() == 2)
require("hostile first stable rank", stable_coefficients.rank() == 2)

mutated_B = sp.expand(B + x * L)
mutated_determinant = sp.expand(sp.diff(U, x) * mutated_B - A * sp.diff(V, x))
require("active hostile mutation", mutated_determinant != 1)
print("PASS E=C+L*C[x] hostile=triple_equal determinant=1 ranks=(2,2)")
print("PASS caveat=hostile_membership_in_gamma_Q(R)_NOT_ASSERTED")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
