#!/usr/bin/env python3
"""Exact controls for proved THM-3612.

The all-regular-function theorem is local and proof-driven.  This companion
uses exact truncated collision germs rather than expanding the full high-degree
fold maps.  It verifies the first-jet branch table, the complete exceptional
quadratic-jet parametrization, cancellation of arbitrary cubic target jets,
the second-difference obstruction, and the sharp tuned controls.
"""

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact gate and fail with a stable label."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact rational-function zero test."""
    return sp.cancel(expression) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


def jacobian(first, second, first_var, second_var):
    """Ordinary two-variable Jacobian in the displayed order."""
    return sp.expand(
        sp.diff(first, first_var) * sp.diff(second, second_var)
        - sp.diff(first, second_var) * sp.diff(second, first_var)
    )


def truncate(expression, variables, degree):
    """Keep terms of total degree at most degree in the listed variables."""
    polynomial = sp.Poly(sp.expand(expression), *variables)
    kept = []
    for monomial, coefficient in polynomial.terms():
        if sum(monomial) <= degree:
            term = coefficient
            for variable, exponent in zip(variables, monomial):
                term *= variable**exponent
            kept.append(term)
    return sp.expand(sum(kept, sp.Integer(0)))


def homogeneous_part(expression, variables, degree):
    """Extract one exact total-degree layer."""
    polynomial = sp.Poly(sp.expand(expression), *variables)
    kept = []
    for monomial, coefficient in polynomial.terms():
        if sum(monomial) == degree:
            term = coefficient
            for variable, exponent in zip(variables, monomial):
                term *= variable**exponent
            kept.append(term)
    return sp.expand(sum(kept, sp.Integer(0)))


print("THM-3612 exact companion -- proved even-fold collision-jet gate")
print("status=verified exact controls plus independent hostile audit; all-regular proof is formal-jet driven")


print("SECTION compiler, Russell relation, and collision")
x, q, w = sp.symbols("x q w")
u = x**2 * q
D = 1 + u
b = sp.expand(u * (u + 3) ** 2)
c = sp.expand(x * (u + 1) * (u + 3))
e = sp.expand(q * (u + 4))
B = sp.expand(b + c * w)
C = c
Y = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
S = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)
Z = sp.expand(S + sp.Rational(3, 4))

require("compiler c2e relation", zero(c**2 * e - b * (b + 4)))
require("Russell target relation", zero(C * Y - B * (B + 4)))
for label, expression in (("B", B), ("C", C), ("Y", Y), ("S", S)):
    require(f"Russell {label} polynomial", sp.denom(sp.cancel(expression)) == 1)

collision_source = (
    (sp.Integer(-1), sp.Integer(-3)),
    (sp.Integer(0), -sp.Rational(3, 4)),
    (sp.Integer(1), sp.Integer(-3)),
)
collision_target = (sp.Integer(0), sp.Integer(0), sp.Integer(0), -sp.Rational(3, 4))
for index, (x_value, q_value) in enumerate(collision_source):
    values = tuple(
        sp.expand(expression.subs({x: x_value, q: q_value, w: 0}))
        for expression in (B, C, Y, S)
    )
    require(f"collision target row={index}", values == collision_target)

target_B, target_C, target_Y = sp.symbols("target_B target_C target_Y")
target_relation = target_C * target_Y - target_B * (target_B + 4)
require(
    "target smooth B derivative",
    sp.diff(target_relation, target_B).subs(
        {target_B: 0, target_C: 0, target_Y: 0}
    )
    == -4,
)
print("PASS compiler_and_collision_gates=10 smooth_derivative=-4")


print("SECTION closed fold, non-graph involution, and residue sign")
t = sp.symbols("t")
Q0 = -sp.Rational(3, 4) - sp.Rational(9, 4) * x**2
Qstar = (
    -sp.Rational(3, 4)
    - sp.Rational(27, 4) * x**2
    + sp.Rational(9, 2) * x**4
)
Qdag = (
    -sp.Rational(3, 4)
    - sp.Rational(27, 2) * x**2
    + 18 * x**4
    - sp.Rational(27, 4) * x**6
)

fold_controls = (
    ("Q0", Q0, -sp.Rational(9, 2), -sp.Rational(9, 2)),
    ("Qstar", Qstar, sp.Rational(9, 2), sp.Rational(81, 2)),
    ("Qdag", Qdag, sp.Rational(9, 2), -sp.Rational(27, 2)),
)
for label, polynomial, expected_first, expected_second in fold_controls:
    require(f"{label} even", zero(polynomial.subs(x, -x) - polynomial))
    require(f"{label} value zero", same(polynomial.subs(x, 0), -sp.Rational(3, 4)))
    require(f"{label} value plus", same(polynomial.subs(x, 1), -3))
    require(f"{label} value minus", same(polynomial.subs(x, -1), -3))
    require(f"{label} first jet", same(sp.diff(polynomial, x).subs(x, 1), expected_first))
    require(f"{label} second jet", same(sp.diff(polynomial, x, 2).subs(x, 1), expected_second))

q_fold = sp.expand(Q0 + t**2)
require("fold q invariant under sheet involution", same(q_fold.subs(t, -t), q_fold))
require("stable coordinate changes under sheet involution", same((-t) - t, -2 * t))
require("stable coordinate not sheet invariant", zero(t - (-t)) == False)

H = q - Q0 - w**2
normal_matrix = sp.Matrix(
    [[sp.diff(expression, variable) for variable in (x, q, w)] for expression in (H, x, w)]
)
require("fold principal normal determinant", normal_matrix.det() == -1)
require("fold residue coefficient", 3 * normal_matrix.det() == -3)

for label, polynomial, _, _ in fold_controls:
    fold_substitution = {q: polynomial + t**2, w: t}
    for index, x_value in enumerate((-1, 0, 1)):
        values = tuple(
            sp.expand(expression.subs(fold_substitution).subs({x: x_value, t: 0}))
            for expression in (B, C, Y, S)
        )
        require(f"{label} fold collision row={index}", values == collision_target)
print("PASS fold_geometry_gates=32 residue=3dx^dt generic_sheet_degree=2")


print("SECTION universal common-cotangent first-jet gate")
v = sp.symbols("v")


def source_gradient(expression, x_value, q_value):
    """Gradient in source coordinates (x,q,w) at a collision point."""
    return tuple(
        sp.expand(sp.diff(expression, variable).subs({x: x_value, q: q_value, w: 0}))
        for variable in (x, q, w)
    )


branch_data = (
    ("minus", -1, -3, -v),
    ("zero", 0, -sp.Rational(3, 4), 0),
    ("plus", 1, -3, v),
)
tangent_rows = {}
for label, x_value, q_value, slope in branch_data:
    rows = []
    for expression in (B, C, Y, Z):
        gradient = source_gradient(expression, x_value, q_value)
        dx_coefficient = sp.expand(gradient[0] + slope * gradient[1])
        dt_coefficient = gradient[2]
        rows.append((dx_coefficient, dt_coefficient))
    tangent_rows[label] = rows

expected_tangent_rows = {
    "minus": (
        (0, 0),
        (12 - 2 * v, 0),
        (-36 + 6 * v, 4),
        ((v - 9) / 2, 0),
    ),
    "zero": ((0, 0), (3, 0), (-9, 4), (0, 0)),
    "plus": (
        (0, 0),
        (12 - 2 * v, 0),
        (-36 + 6 * v, 4),
        ((9 - v) / 2, 0),
    ),
}
for label in ("minus", "zero", "plus"):
    for coordinate, (actual, expected) in enumerate(
        zip(tangent_rows[label], expected_tangent_rows[label])
    ):
        require(
            f"tangent row={label} coordinate={coordinate}",
            same(actual[0], expected[0]) and same(actual[1], expected[1]),
        )

Acoef, Kcoef, Rcoef = sp.symbols("Acoef Kcoef Rcoef")


def two_form_value(rows):
    """Evaluate A dC^dY+K dC^dZ+R dY^dZ on (dx,dt)."""
    c_row, y_row, z_row = rows[1], rows[2], rows[3]

    def determinant(left, right):
        return sp.expand(left[0] * right[1] - left[1] * right[0])

    return sp.expand(
        Acoef * determinant(c_row, y_row)
        + Kcoef * determinant(c_row, z_row)
        + Rcoef * determinant(y_row, z_row)
    )


j_minus = two_form_value(tangent_rows["minus"])
j_zero = two_form_value(tangent_rows["zero"])
j_plus = two_form_value(tangent_rows["plus"])
require("first jet middle value", same(j_zero, 12 * Acoef))
require(
    "first jet plus value",
    same(j_plus, (48 - 8 * v) * Acoef + 2 * (v - 9) * Rcoef),
)
require(
    "first jet minus value",
    same(j_minus, (48 - 8 * v) * Acoef + 2 * (9 - v) * Rcoef),
)
require("side difference factor", same(j_plus - j_minus, 4 * (v - 9) * Rcoef))
require(
    "side-middle after R zero",
    same((j_plus - j_zero).subs(Rcoef, 0), 4 * (9 - 2 * v) * Acoef),
)
require("v=9 hostile side", same(j_plus.subs(v, 9), -24 * Acoef))
require("v=9 hostile middle", same(j_zero.subs(v, 9), 12 * Acoef))
require(
    "v=9/2 tangent match",
    same(j_plus.subs({v: sp.Rational(9, 2), Rcoef: 0}), 12 * Acoef)
    and same(j_minus.subs({v: sp.Rational(9, 2), Rcoef: 0}), 12 * Acoef),
)

hostile_first_jet_rows = 0
for numerator in range(-18, 37):
    test_v = sp.Rational(numerator, 2)
    if test_v == sp.Rational(9, 2):
        continue
    coefficient_matrix = sp.Matrix(
        [
            [sp.diff(j_plus - j_zero, Acoef), sp.diff(j_plus - j_zero, Rcoef)],
            [sp.diff(j_minus - j_zero, Acoef), sp.diff(j_minus - j_zero, Rcoef)],
        ]
    ).subs(v, test_v)
    nullspace = coefficient_matrix.nullspace()
    common_nonzero = False
    for vector in nullspace:
        if sp.expand(j_zero.subs({v: test_v, Acoef: vector[0], Rcoef: vector[1]})) != 0:
            common_nonzero = True
    require(f"finite first-jet hostile v={test_v}", common_nonzero == False)
    hostile_first_jet_rows += 1
print(f"PASS first_jet_symbolic=20 hostile_rows={hostile_first_jet_rows} unique_escape=9/2")


print("SECTION exceptional Hessian solution and cubic second difference")
xi = sp.symbols("xi")
r0, r1 = sp.symbols("r0 r1")
beta, ell = sp.symbols("beta ell")
target_c, target_y, target_z = sp.symbols("target_c target_y target_z")
target_variables = (target_c, target_y, target_z)


def local_branch(x_value, q_value, slope, second_derivative, degree):
    """Truncated exceptional fold germ in target parameters (C,Y,Z)."""
    local_q = q_value + slope * xi + second_derivative * xi**2 / 2 + t**2
    substitution = {x: x_value + xi, q: local_q, w: t}
    return tuple(
        truncate(expression.subs(substitution), (xi, t), degree)
        for expression in (C, Y, Z)
    )


exceptional_branches = (
    local_branch(-1, -3, -sp.Rational(9, 2), r1, 3),
    local_branch(0, -sp.Rational(3, 4), 0, r0, 3),
    local_branch(1, -3, sp.Rational(9, 2), r1, 3),
)

g_CC, g_CY, g_CZ, g_YY, g_YZ, g_ZZ = sp.symbols(
    "g_CC g_CY g_CZ g_YY g_YZ g_ZZ"
)
f_CC = -g_CY - 3 * beta / 8
f_CY = -g_YY - 7 * beta / 32
f_CZ = 1 + 4 * r1 / 27 - 7 * beta**2 / 64 - beta * g_YY / 2 - g_YZ / 2
f_YY = ell
f_YZ = beta * ell
f_ZZ = beta**2 * ell

F_quadratic = sp.expand(
    target_c
    + f_CC * target_c**2 / 2
    + f_CY * target_c * target_y
    + f_CZ * target_c * target_z
    + f_YY * target_y**2 / 2
    + f_YZ * target_y * target_z
    + f_ZZ * target_z**2 / 2
)
G_quadratic = sp.expand(
    target_y
    + beta * target_z
    + g_CC * target_c**2 / 2
    + g_CY * target_c * target_y
    + g_CZ * target_c * target_z
    + g_YY * target_y**2 / 2
    + g_YZ * target_y * target_z
    + g_ZZ * target_z**2 / 2
)


def compose_target(expression, branch, degree):
    """Compose a target Taylor polynomial with one truncated source germ."""
    substitution = dict(zip(target_variables, branch))
    return truncate(expression.subs(substitution), (xi, t), degree)


raw_f = sp.symbols("raw_f_0:6")
raw_g = sp.symbols("raw_g_0:6")
quadratic_pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def quadratic_form(coefficients):
    """Target quadratic with Hessian-normalized diagonal coefficients."""
    terms = []
    for coefficient, (left, right) in zip(coefficients, quadratic_pairs):
        divisor = 2 if left == right else 1
        terms.append(
            coefficient * target_variables[left] * target_variables[right] / divisor
        )
    return sp.expand(sum(terms, sp.Integer(0)))


F_raw = target_c + quadratic_form(raw_f)
G_raw = target_y + beta * target_z + quadratic_form(raw_g)
raw_equations = []
for branch in exceptional_branches:
    raw_local_F = compose_target(F_raw, branch, 2)
    raw_local_G = compose_target(G_raw, branch, 2)
    raw_local_J = truncate(jacobian(raw_local_F, raw_local_G, xi, t), (xi, t), 1)
    raw_equations.extend(
        (
            sp.Poly(raw_local_J - 12, xi, t).coeff_monomial((1, 0)),
            sp.Poly(raw_local_J - 12, xi, t).coeff_monomial((0, 1)),
        )
    )
raw_matrix, raw_rhs = sp.linear_eq_to_matrix(raw_equations, raw_f + raw_g)
require("exceptional Hessian coefficient rank", raw_matrix.rank() == 5)
require(
    "exceptional Hessian augmented rank",
    raw_matrix.row_join(raw_rhs).rank() == 5,
)


for index, branch in enumerate(exceptional_branches):
    local_F = compose_target(F_quadratic, branch, 2)
    local_G = compose_target(G_quadratic, branch, 2)
    local_J = truncate(jacobian(local_F, local_G, xi, t), (xi, t), 1)
    require(f"exceptional constant row={index}", same(local_J.subs({xi: 0, t: 0}), 12))
    require(
        f"exceptional xi-linear row={index}",
        sp.Poly(local_J - 12, xi, t).coeff_monomial((1, 0)) == 0,
    )
    require(
        f"exceptional t-linear row={index}",
        sp.Poly(local_J - 12, xi, t).coeff_monomial((0, 1)) == 0,
    )

cubic_monomials = []
for c_degree in range(4):
    for y_degree in range(4 - c_degree):
        z_degree = 3 - c_degree - y_degree
        cubic_monomials.append(
            target_c**c_degree * target_y**y_degree * target_z**z_degree
        )
f_cubic_coefficients = sp.symbols(f"f_cubic_0:{len(cubic_monomials)}")
g_cubic_coefficients = sp.symbols(f"g_cubic_0:{len(cubic_monomials)}")
F_cubic = sp.expand(
    F_quadratic
    + sum(
        coefficient * monomial
        for coefficient, monomial in zip(f_cubic_coefficients, cubic_monomials)
    )
)
G_cubic = sp.expand(
    G_quadratic
    + sum(
        coefficient * monomial
        for coefficient, monomial in zip(g_cubic_coefficients, cubic_monomials)
    )
)

vertical_t2_coefficients = []
for index, branch in enumerate(exceptional_branches):
    local_F = compose_target(F_cubic, branch, 3)
    local_G = compose_target(G_cubic, branch, 3)
    local_J = truncate(jacobian(local_F, local_G, xi, t), (xi, t), 2)
    coefficient = sp.Poly(local_J - 12, xi, t).coeff_monomial((0, 2))
    vertical_t2_coefficients.append(coefficient)
    require(
        f"cubic leaves constant row={index}",
        sp.Poly(local_J, xi, t).coeff_monomial((0, 0)) == 12,
    )
    require(
        f"cubic leaves t-linear row={index}",
        sp.Poly(local_J - 12, xi, t).coeff_monomial((0, 1)) == 0,
    )

second_difference = sp.factor(
    vertical_t2_coefficients[0]
    - 2 * vertical_t2_coefficients[1]
    + vertical_t2_coefficients[2]
)
expected_second_difference = -sp.Rational(16, 3) * (2 * r1 + 27)
require("universal cubic second difference", same(second_difference, expected_second_difference))
cubic_symbols = set(f_cubic_coefficients) | set(g_cubic_coefficients)
require("all arbitrary cubic coefficients cancel", second_difference.free_symbols.isdisjoint(cubic_symbols))
require(
    "Qstar minus576 invoice",
    second_difference.subs(r1, sp.Rational(81, 2)) == -576,
)
require(
    "Qdag tuned invoice vanishes",
    second_difference.subs(r1, -sp.Rational(27, 2)) == 0,
)

hostile_second_jet_rows = 0
for numerator in range(-40, 41):
    test_r = sp.Rational(numerator, 2)
    value = sp.expand(expected_second_difference.subs(r1, test_r))
    expected_zero = test_r == -sp.Rational(27, 2)
    require(f"finite second-jet hostile r={test_r}", (value == 0) == expected_zero)
    hostile_second_jet_rows += 1
print(
    "PASS exceptional_hessian_rank=5 cubic_coefficients=20 "
    f"second_jet_hostiles={hostile_second_jet_rows} invariant=-16(2r+27)/3"
)


print("SECTION sharp Qstar first-order survivor")
Fstar = sp.expand(target_c + 7 * target_c * target_z)
Gstar = target_y
qstar_branches = (
    local_branch(-1, -3, -sp.Rational(9, 2), sp.Rational(81, 2), 3),
    local_branch(0, -sp.Rational(3, 4), 0, -sp.Rational(27, 2), 3),
    local_branch(1, -3, sp.Rational(9, 2), sp.Rational(81, 2), 3),
)
qstar_vertical = []
for index, branch in enumerate(qstar_branches):
    local_F = compose_target(Fstar, branch, 3)
    local_G = compose_target(Gstar, branch, 3)
    local_J = truncate(jacobian(local_F, local_G, xi, t), (xi, t), 2)
    require(f"Qstar survivor constant row={index}", local_J.subs({xi: 0, t: 0}) == 12)
    require(
        f"Qstar survivor first jet row={index}",
        homogeneous_part(local_J - 12, (xi, t), 1) == 0,
    )
    qstar_vertical.append(sp.Poly(local_J - 12, xi, t).coeff_monomial((0, 2)))
require("Qstar vertical coefficients", tuple(qstar_vertical) == (-141, 147, -141))
require("Qstar vertical second difference", qstar_vertical[0] - 2 * qstar_vertical[1] + qstar_vertical[2] == -576)
print("PASS Qstar_survivor_gates=8 vertical=(-141,147,-141)")


print("SECTION doubly tuned Qdag second-order survivor")
Fdag = sp.expand(
    target_c
    - target_c * target_z
    - sp.Rational(209, 144) * target_c**3
    + target_c**2 * target_y / 16
    + sp.Rational(7, 64) * target_c * target_y**2
    - sp.Rational(755, 108) * target_c * target_z**2
)
Gdag = target_y

qdag_fold = {q: Qdag + t**2, w: t}
qdag_branches = []
for x_value in (-1, 0, 1):
    substitution = {x: x_value + xi}
    qdag_branches.append(
        tuple(
            truncate(expression.subs(qdag_fold).subs(substitution), (xi, t), 4)
            for expression in (C, Y, Z)
        )
    )

qdag_third_layers = []
for index, branch in enumerate(qdag_branches):
    local_F = compose_target(Fdag, branch, 4)
    local_G = compose_target(Gdag, branch, 4)
    local_J = truncate(jacobian(local_F, local_G, xi, t), (xi, t), 3)
    require(f"Qdag second-order survivor row={index}", truncate(local_J - 12, (xi, t), 2) == 0)
    qdag_third_layers.append(homogeneous_part(local_J - 12, (xi, t), 3))
require("Qdag survivor has higher error minus", qdag_third_layers[0] != 0)
require("Qdag survivor has higher error plus", qdag_third_layers[2] != 0)
require("Qdag third layers antisymmetric", same(qdag_third_layers[0], -qdag_third_layers[2]))
require("Qdag middle third layer zero", qdag_third_layers[1] == 0)
print("PASS Qdag_second_order_survivor_gates=7 higher_jet_nonzero=TRUE")


print("SECTION affine closure and coordinate-plane positive controls")
affine_difference_matrix = sp.Matrix(
    [[1, -sp.Rational(9, 4)], [-1, -sp.Rational(9, 4)]]
)
require("affine collision points noncollinear", affine_difference_matrix.det() == -sp.Rational(9, 2))

x_zero_values = tuple(sp.expand(expression.subs(x, 0)) for expression in (B, C, Y, S))
require(
    "x=0 target coordinates",
    x_zero_values == (0, 0, 4 * w, q + sp.Rational(3, 4) * w**2),
)
require(
    "x=0 YS Jacobian",
    jacobian(x_zero_values[2], x_zero_values[3], q, w) == -4,
)
require(
    "x=0 triangular inverse q",
    same(
        x_zero_values[3] - 3 * x_zero_values[2] ** 2 / 64,
        q,
    ),
)

W_target = sp.expand(Y * (B + 2) / 8 - C * S)
require("inverse cylinder stable coordinate", same(W_target, w))
q_zero_C = sp.expand(C.subs(q, 0))
q_zero_W = sp.expand(W_target.subs(q, 0))
require("q=0 CW coordinates", (q_zero_C, q_zero_W) == (3 * x, w))
require("q=0 CW Jacobian", jacobian(q_zero_C, q_zero_W, x, w) == 3)
print("PASS affine_and_positive_controls=7 x0_jac=-4 q0_jac=3")


print("SECTION source AST gate")
source_path = Path(__file__)
source_tree = ast.parse(source_path.read_text(encoding="utf-8"))
assertion_nodes = sum(1 for node in ast.walk(source_tree) if isinstance(node, ast.Assert))
require("AST assertion nodes absent", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
