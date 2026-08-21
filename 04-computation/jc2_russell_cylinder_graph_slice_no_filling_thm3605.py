#!/usr/bin/env python3
"""Exact controls for provisional THM-3605.

The universal graph obstructions are proved by the weight and divisibility
arguments in the theorem.  This companion verifies every displayed algebraic
identity, both quotient-ring compositions of the cylinder maps, finite hostile
ansatz rows, and positive punctured controls with exact SymPy arithmetic.
"""

from itertools import product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact gate and raise with a stable label on failure."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expr):
    """Exact rational-function zero test."""
    return sp.cancel(expr) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


def jacobian(first, second, first_var, second_var):
    """Ordinary two-variable Jacobian in the displayed variable order."""
    return sp.expand(
        sp.diff(first, first_var) * sp.diff(second, second_var)
        - sp.diff(first, second_var) * sp.diff(second, first_var)
    )


print("THM-3605 exact companion -- provisional Russell-cylinder graph-slice gate")
print("status=finite exact controls; universal graph/no-mate conclusions are proof-driven")


print("SECTION polynomial Russell specialization and inverse")
b, c, e, w = sp.symbols("b c e w")
B, C, Y, S = sp.symbols("B C Y S")
T = sp.symbols("T")


def P(value):
    return sp.expand(value * (value + 4))


def g(value):
    return sp.expand(value - P(value) * (value + 2) / 8)


def H_factor(value):
    return sp.expand((value - 2) * (value + 6) / 64)


U = (T + 2) / 8
V = -sp.Rational(1, 4)
require("Bezout identity", same(U * sp.diff(P(T), T) + V * P(T), 1))
require("g expanded", same(g(T), -T**2 * (T + 6) / 8))
require("P(g) factorization", same(P(g(T)), P(T) ** 2 * H_factor(T)))

# Forward point map Y_2 x A1_w -> Y_1 x A1_S.
C_forward = c
B_forward = b + c * w
Y_forward = c * e + (2 * b + 4) * w + c * w**2
S_forward = ((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8

# Polynomial inverse point map.
b_inverse = g(B) + C**2 * S
w_inverse = Y * (B + 2) / 8 - C * S
e_inverse = Y**2 * H_factor(B) + 2 * (g(B) + 2) * S + C**2 * S**2


def reduce_y2(expr):
    """Reduce a polynomial modulo b^2+4b-c^2e, monic in b."""
    dividend = sp.Poly(sp.expand(expr), b)
    relation = sp.Poly(b**2 + 4 * b - c**2 * e, b)
    return sp.cancel(sp.rem(dividend, relation).as_expr())


def reduce_y1(expr):
    """Reduce a polynomial modulo B^2+4B-CY, monic in B."""
    dividend = sp.Poly(sp.expand(expr), B)
    relation = sp.Poly(B**2 + 4 * B - C * Y, B)
    return sp.cancel(sp.rem(dividend, relation).as_expr())


require(
    "forward target relation",
    reduce_y2(C_forward * Y_forward - P(B_forward)) == 0,
)
require(
    "inverse source relation",
    reduce_y1(C**2 * e_inverse - P(b_inverse)) == 0,
)

forward_substitution = {
    C: C_forward,
    B: B_forward,
    Y: Y_forward,
    S: S_forward,
}
for label, inverse_coordinate, expected in (
    ("b", b_inverse, b),
    ("e", e_inverse, e),
    ("w", w_inverse, w),
):
    composed = inverse_coordinate.subs(forward_substitution, simultaneous=True)
    require(
        f"inverse-after-forward {label}",
        reduce_y2(composed - expected) == 0,
    )

inverse_substitution = {
    c: C,
    b: b_inverse,
    e: e_inverse,
    w: w_inverse,
}
for label, forward_coordinate, expected in (
    ("C", C_forward, C),
    ("B", B_forward, B),
    ("Y", Y_forward, Y),
    ("S", S_forward, S),
):
    composed = forward_coordinate.subs(inverse_substitution, simultaneous=True)
    require(
        f"forward-after-inverse {label}",
        reduce_y1(composed - expected) == 0,
    )

for label, coordinate, variables in (
    ("C forward polynomial", C_forward, (b, c, e, w)),
    ("B forward polynomial", B_forward, (b, c, e, w)),
    ("Y forward polynomial", Y_forward, (b, c, e, w)),
    ("S forward polynomial", S_forward, (b, c, e, w)),
    ("b inverse polynomial", b_inverse, (B, C, Y, S)),
    ("e inverse polynomial", e_inverse, (B, C, Y, S)),
    ("w inverse polynomial", w_inverse, (B, C, Y, S)),
):
    require(label, sp.denom(sp.cancel(coordinate)) == 1)
    require(f"{label} parsed", sp.Poly(sp.expand(coordinate), *variables) is not None)

e_open = P(b) / c**2
S_open = sp.cancel(S_forward.subs(e, e_open))
volume_matrix = sp.Matrix(
    [[sp.diff(value, variable) for variable in (c, b, w)]
     for value in (C_forward, B_forward, S_open)]
)
volume_determinant = sp.factor(volume_matrix.det())
require("open-chart cylinder determinant", same(volume_determinant, -1 / c))
require(
    "residue-volume coefficient",
    same(volume_determinant / C_forward, -1 / c**2),
)
print("PASS cylinder quotient compositions=7 polynomial-coordinate gates=14 volume gates=2")


print("SECTION THM-3561 stabilization and collision curve")
x, q = sp.symbols("x q")
D = 1 + x**2 * q
a = q / D**2
c_compiler = x * D * (D + 2)
b_compiler = (D - 1) * (D + 2) ** 2
e_compiler = q * (D + 3)

require(
    "compiler b=ac^2",
    same(b_compiler, a * c_compiler**2),
)
require(
    "compiler e=a(b+4)",
    same(e_compiler, a * (b_compiler + 4)),
)
require(
    "compiler Y2 relation",
    same(c_compiler**2 * e_compiler, P(b_compiler)),
)
require(
    "rational Keller determinant",
    same(jacobian(a, c_compiler, x, q), -3),
)

for label, coordinate in (
    ("b compiler polynomial", b_compiler),
    ("c compiler polynomial", c_compiler),
    ("e compiler polynomial", e_compiler),
):
    require(label, sp.denom(sp.cancel(coordinate)) == 1)
    require(f"{label} parsed", sp.Poly(sp.expand(coordinate), x, q) is not None)

# Pullback of dc wedge db/c^2; its negative is the transported stable volume.
omega2_coefficient = sp.cancel(
    jacobian(c_compiler, b_compiler, x, q) / c_compiler**2
)
require("compiler residue two-form", same(omega2_coefficient, 3))
require("stable volume after cylinder", same(-omega2_coefficient, -3))

compiler_forward = {
    b: b_compiler,
    c: c_compiler,
    e: e_compiler,
}
C_stable = sp.expand(C_forward.subs(compiler_forward))
B_stable = sp.expand(B_forward.subs(compiler_forward))
Y_stable = sp.expand(Y_forward.subs(compiler_forward))
S_stable = sp.expand(S_forward.subs(compiler_forward))
for label, coordinate in (
    ("C stable polynomial", C_stable),
    ("B stable polynomial", B_stable),
    ("Y stable polynomial", Y_stable),
    ("S stable polynomial", S_stable),
):
    require(label, sp.denom(sp.cancel(coordinate)) == 1)
    require(f"{label} parsed", sp.Poly(sp.expand(coordinate), x, q, w) is not None)

collision_points = (
    (sp.Integer(0), -sp.Rational(3, 4)),
    (sp.Integer(1), sp.Integer(-3)),
    (sp.Integer(-1), sp.Integer(-3)),
)
for index, (x_value, q_value) in enumerate(collision_points):
    substitution = {x: x_value, q: q_value}
    require(f"collision a row={index}", same(a.subs(substitution), -sp.Rational(3, 4)))
    require(f"collision c row={index}", c_compiler.subs(substitution) == 0)
    require(f"collision b row={index}", b_compiler.subs(substitution) == 0)
    require(f"collision e row={index}", e_compiler.subs(substitution) == -3)
    require(f"collision C curve row={index}", C_stable.subs(substitution) == 0)
    require(f"collision B curve row={index}", B_stable.subs(substitution) == 0)
    require(f"collision Y curve row={index}", same(Y_stable.subs(substitution), 4 * w))
    require(
        f"collision S curve row={index}",
        same(S_stable.subs(substitution), 3 * (w**2 - 1) / 4),
    )
print("PASS compiler identities=12 collision rows=3x8 stable polynomial gates=8")


print("SECTION collision-bearing arm and polynomial source-graph hostiles")
A_stable = sp.cancel(B_stable / C_stable)
arm_base = sp.cancel(a * c_compiler)
require("root-zero arm transport", same(A_stable, arm_base + w))
require("explicit arm pole formula", same(arm_base, x * q * (D + 2) / D))

arm_numerator, arm_denominator = sp.fraction(sp.cancel(arm_base))
require("arm pole denominator", same(arm_denominator, D))
require("arm numerator coprime", sp.gcd(arm_numerator, arm_denominator) == 1)

source_graphs = (
    sp.Integer(0),
    sp.Integer(1),
    x,
    q,
    x * q,
    x**2 + q,
    3 * x - 2 * q + 5,
    x**3 * q**2 + 2 * x * q - 7,
)
for index, h_source in enumerate(source_graphs):
    A_graph = sp.cancel(arm_base + h_source)
    numerator, denominator = sp.fraction(A_graph)
    require(f"source graph pole denominator row={index}", same(denominator, D))
    require(
        f"source graph pole coprime row={index}",
        sp.gcd(numerator, denominator) == 1,
    )
    arm_form = jacobian(c_compiler, A_graph, x, q)
    expected_form = 3 * c_compiler - jacobian(h_source, c_compiler, x, q)
    require(f"source graph arm-form identity row={index}", same(arm_form, expected_form))

alpha, beta, gamma = sp.symbols("alpha beta gamma")
h_affine = alpha * x + beta * q + gamma
affine_values = [
    sp.expand(h_affine.subs({x: x_value, q: q_value}))
    for x_value, q_value in collision_points
]
affine_solution = sp.solve(
    [affine_values[0] - affine_values[1], affine_values[0] - affine_values[2]],
    (alpha, beta),
    dict=True,
)
require("affine collision classification", affine_solution == [{alpha: 0, beta: 0}])
constant_arm_form = jacobian(c_compiler, arm_base + gamma, x, q)
require("constant graph form is 3c", same(constant_arm_form, 3 * c_compiler))
require("constant graph form is nonunit", sp.Poly(c_compiler, x, q).total_degree() > 0)
print("PASS root-zero arm identities=4 polynomial graph hostiles=8x3 affine gates=3")


print("SECTION weight identity and source-graph ODE obstruction")
u = sp.symbols("u")
p = sp.expand((u + 1) * (u + 3))
F_samples = (
    sp.Integer(1),
    u,
    1 + u**2,
    u * (u + 2),
    2 - 3 * u + u**3,
)
weight_rows = 0
for m, F_sample in product(range(-5, 6), F_samples):
    f_xq = x**m * F_sample.subs(u, x**2 * q)
    c_xq = x * p.subs(u, x**2 * q)
    expected = x ** (m + 2) * (
        m * F_sample * sp.diff(p, u) - sp.diff(F_sample, u) * p
    ).subs(u, x**2 * q)
    require(
        f"weight Hamiltonian formula m={m} Frow={F_samples.index(F_sample)}",
        same(jacobian(f_xq, c_xq, x, q), expected),
    )
    weight_rows += 1

K0 = sp.symbols("K0")
source_ode_numerator = -u * (u + 3) ** 2 + K0
require(
    "source ODE antiderivative",
    same(sp.diff(source_ode_numerator, u), -3 * p),
)
require("source ODE root minus3", source_ode_numerator.subs(u, -3) == K0)
require("source ODE forced K0", source_ode_numerator.subs({u: -1, K0: 0}) == 4)

# Positive one-arm boundary: p_1=u admits a lawful weight-minus-one solution.
p_one = u
F_one = -sp.Rational(3, 2) * u
h_one = sp.cancel(x ** (-1) * F_one.subs(u, x**2 * q))
c_one = x * p_one.subs(u, x**2 * q)
require("one-arm ODE positive", same(-sp.diff(F_one * p_one, u), 3 * p_one))
require("one-arm h polynomial", sp.Poly(h_one, x, q) is not None)
require("one-arm direct positive", same(jacobian(h_one, c_one, x, q), 3 * c_one))


def finite_ansatz_has_no_solution(total_degree, right_hand_side, prefix):
    """Exact bounded hostile for a polynomial Hamiltonian equation."""
    coefficients = []
    candidate = 0
    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            coefficient = sp.symbols(f"{prefix}_{total_degree}_{i}_{j}")
            coefficients.append(coefficient)
            candidate += coefficient * x**i * q**j
    equation = sp.Poly(
        sp.expand(jacobian(candidate, c_compiler, x, q) - right_hand_side),
        x,
        q,
    )
    solution = sp.linsolve(equation.coeffs(), coefficients)
    return solution is sp.EmptySet


for degree in range(0, 9):
    require(
        f"bounded source graph no-unit degree={degree}",
        finite_ansatz_has_no_solution(degree, 3 * c_compiler - 1, "h"),
    )
print(f"PASS weight formula rows={weight_rows} ODE gates=6 bounded no-unit degrees=0..8")


print("SECTION complete punctured-Darboux graph classification")
aa, cc = sp.symbols("aa cc")
mu_values = (-3, -1, 1, 2, sp.Rational(5, 2))
G_values = (
    sp.Integer(0),
    cc,
    cc**2 + 2,
    cc**3 - 2 * cc + 5,
)
punctured_rows = 0
for mu_value, G_value in product(mu_values, G_values):
    H_graph = (mu_value - cc) * aa + G_value
    A_output = sp.expand(cc * aa + H_graph)
    A_expected = sp.expand(mu_value * aa + G_value)
    b_chart = aa * cc**2
    e_chart = aa * (b_chart + 4)
    chart_substitution = {b: b_chart, c: cc, e: e_chart, w: H_graph}
    B_output = sp.expand(B_forward.subs(chart_substitution))
    Y_output = sp.expand(Y_forward.subs(chart_substitution))
    S_output = sp.expand(S_forward.subs(chart_substitution))
    S_expected = sp.expand(
        cc * A_expected**3 / 8
        + 3 * A_expected**2 / 4
        + (A_expected - G_value) / mu_value
    )

    require(
        f"punctured derivative classification mu={mu_value} Grow={G_values.index(G_value)}",
        same(cc + sp.diff(H_graph, aa), mu_value),
    )
    require(
        f"punctured A triangular mu={mu_value} Grow={G_values.index(G_value)}",
        same(A_output, A_expected),
    )
    require(
        f"punctured chart Jacobian mu={mu_value} Grow={G_values.index(G_value)}",
        jacobian(A_output, cc, aa, cc) == mu_value,
    )
    require(
        f"punctured B graph mu={mu_value} Grow={G_values.index(G_value)}",
        same(B_output, cc * A_expected),
    )
    require(
        f"punctured Y graph mu={mu_value} Grow={G_values.index(G_value)}",
        same(Y_output, A_expected * (cc * A_expected + 4)),
    )
    require(
        f"punctured S graph mu={mu_value} Grow={G_values.index(G_value)}",
        same(S_output, S_expected),
    )

    source_A = sp.cancel(mu_value * a + G_value.subs(cc, c_compiler))
    require(
        f"punctured source Jacobian mu={mu_value} Grow={G_values.index(G_value)}",
        same(jacobian(source_A, c_compiler, x, q), -3 * mu_value),
    )
    source_numerator, source_denominator = sp.fraction(source_A)
    denominator_ratio = sp.cancel(source_denominator / D**2)
    require(
        f"punctured double pole mu={mu_value} Grow={G_values.index(G_value)}",
        denominator_ratio != 0 and not denominator_ratio.has(x, q),
    )
    require(
        f"punctured double pole coprime mu={mu_value} Grow={G_values.index(G_value)}",
        sp.gcd(source_numerator, source_denominator) == 1,
    )

    collision_A = A_expected.subs({aa: -sp.Rational(3, 4), cc: 0})
    collision_w = H_graph.subs({aa: -sp.Rational(3, 4), cc: 0})
    collision_S = S_expected.subs({aa: -sp.Rational(3, 4), cc: 0})
    require(
        f"punctured collision w=A mu={mu_value} Grow={G_values.index(G_value)}",
        same(collision_w, collision_A),
    )
    require(
        f"punctured collision S curve mu={mu_value} Grow={G_values.index(G_value)}",
        same(collision_S, 3 * (collision_A**2 - 1) / 4),
    )
    punctured_rows += 1

affine_a, affine_c, affine_constant = sp.symbols(
    "affine_a affine_c affine_constant"
)
H_affine_punctured = affine_a * aa + affine_c * cc + affine_constant
affine_form_coefficient = cc + sp.diff(H_affine_punctured, aa)
require("punctured affine coefficient varies", sp.diff(affine_form_coefficient, cc) == 1)
require("punctured affine graph never unit", affine_form_coefficient.has(cc))

# Cheapest positive punctured member mu=1,G=0.
H_positive = (1 - cc) * aa
A_positive = sp.expand(cc * aa + H_positive)
S_positive = sp.expand(cc * A_positive**3 / 8 + 3 * A_positive**2 / 4 + A_positive)
require("positive punctured A=a", A_positive == aa)
require("positive punctured Jacobian", jacobian(A_positive, cc, aa, cc) == 1)
require(
    "positive punctured collision stable value",
    H_positive.subs({aa: -sp.Rational(3, 4), cc: 0}) == -sp.Rational(3, 4),
)
require(
    "positive punctured collision S",
    S_positive.subs({aa: -sp.Rational(3, 4), cc: 0}) == -sp.Rational(21, 64),
)
print(f"PASS punctured graph rows={punctured_rows} gates_per_row=11 affine=2 positive=4")


print("SECTION fixed-c no-polynomial-mate obstruction")
I = u**3 / 3 + 2 * u**2 + 3 * u
require("p antiderivative", same(sp.diff(I, u), p))
require("antiderivative at minus3", I.subs(u, -3) == 0)
require("antiderivative at minus1", I.subs(u, -1) == -sp.Rational(4, 3))

lambda_values = (-5, -1, 1, 2, sp.Rational(7, 3))
mate_rows = 0
for lambda_value in lambda_values:
    integrated = -lambda_value * I + K0
    require(
        f"mate integrated ODE lambda={lambda_value}",
        same(sp.diff(integrated, u), -lambda_value * p),
    )
    require(
        f"mate root minus3 forces K0 lambda={lambda_value}",
        integrated.subs(u, -3) == K0,
    )
    require(
        f"mate root minus1 hostile lambda={lambda_value}",
        integrated.subs({u: -1, K0: 0}) == sp.Rational(4, 3) * lambda_value,
    )
    require(
        f"mate root incompatibility lambda={lambda_value}",
        integrated.subs({u: -1, K0: 0}) != 0,
    )
    mate_rows += 1

for degree in range(0, 9):
    require(
        f"bounded fixed-c no-mate degree={degree}",
        finite_ansatz_has_no_solution(degree, sp.Integer(1), "f"),
    )

# Zero right-hand side is the sharp nonunit boundary.
require("zero-mate boundary", jacobian(c_compiler**2 + 1, c_compiler, x, q) == 0)
print(f"PASS fixed-c root rows={mate_rows}x4 bounded no-mate degrees=0..8 zero boundary=1")


print("SECTION status and scope sentinels")
require("source puncture is nonempty", D.subs({x: 1, q: -1}) == 0)
require("collision avoids puncture r0", D.subs({x: 0, q: -sp.Rational(3, 4)}) == 1)
require("collision avoids puncture rplus", D.subs({x: 1, q: -3}) == -2)
require("collision avoids puncture rminus", D.subs({x: -1, q: -3}) == -2)
print("PASS dimension=3 stable; graph scopes only; non-graph and S-mixing projections remain OPEN")
print(f"ALL EXACT CHECKS PASSED gates={CHECKS}")
