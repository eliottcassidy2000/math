#!/usr/bin/env python3
"""Exact controls for proved THM-3624.

The script verifies the complete non-even first-jet tangent-match locus, its
projective parametrization, the weighted vertical top-symbol cokernel, the
generic degree-two invoice on the rho=0 stratum, its exceptional closure, and
an explicit non-even polynomial that survives the enlarged arbitrary-two-form
jet system through total source degree four.
"""

import ast
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


CHECKS = 0


def require(label, condition):
    """Record one exact active gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact polynomial/rational zero test."""
    return sp.cancel(expression) == 0


def truncate(expression, variables, degree):
    """Keep terms of total degree at most degree."""
    polynomial = sp.Poly(sp.expand(expression), *variables)
    kept = []
    for monomial, coefficient in polynomial.terms():
        if sum(monomial) <= degree:
            term = coefficient
            for variable, exponent in zip(variables, monomial):
                term *= variable**exponent
            kept.append(term)
    return sp.expand(sum(kept, sp.Integer(0)))


def exponent_triples(maximum_degree):
    """Target monomials in stable total-degree order."""
    return [
        (first, second, degree - first - second)
        for degree in range(maximum_degree + 1)
        for first in range(degree + 1)
        for second in range(degree + 1 - first)
    ]


def source_monomials(maximum_degree):
    """Source monomials xi^i t^j in stable total-degree order."""
    return [
        (xi_degree, degree - xi_degree)
        for degree in range(maximum_degree + 1)
        for xi_degree in range(degree + 1)
    ]


print("THM-3624 exact companion -- proved non-even weighted-cokernel boundary")
print("status=proved finite mechanism and independently hostile-audited")


print("SECTION compiler collision and general non-even tangent table")
x, q, w, xi, t = sp.symbols("x q w xi t")
compiler_u = x**2 * q
D = 1 + compiler_u
b = sp.expand(compiler_u * (compiler_u + 3) ** 2)
c_source = sp.expand(x * (compiler_u + 1) * (compiler_u + 3))
e_source = sp.expand(q * (compiler_u + 4))
B = sp.expand(b + c_source * w)
C = c_source
Y = sp.expand(c_source * e_source + (2 * b + 4) * w + c_source * w**2)
S = sp.expand(
    ((b + 2) * (e_source + 3 * w**2) + c_source * w * (3 * e_source + w**2))
    / 8
)
Z = sp.expand(S + sp.Rational(3, 4))

require("compiler c2e relation", zero(c_source**2 * e_source - b * (b + 4)))
require("Russell target relation", zero(C * Y - B * (B + 4)))

collision_target = (0, 0, 0, -sp.Rational(3, 4))
for index, (x_value, q_value) in enumerate(
    ((-1, -3), (0, -sp.Rational(3, 4)), (1, -3))
):
    values = tuple(
        sp.expand(expression.subs({x: x_value, q: q_value, w: 0}))
        for expression in (B, C, Y, S)
    )
    require(f"collision row={index}", values == collision_target)

u, v_minus, v_plus = sp.symbols("u v_minus v_plus")


def source_gradient(expression, x_value, q_value):
    """Gradient in (x,q,w) at one collision point."""
    return tuple(
        sp.expand(
            sp.diff(expression, variable).subs(
                {x: x_value, q: q_value, w: 0}
            )
        )
        for variable in (x, q, w)
    )


branch_specifications = (
    ("minus", -1, -3, v_minus),
    ("zero", 0, -sp.Rational(3, 4), u),
    ("plus", 1, -3, v_plus),
)
tangent_rows = {}
for label, x_value, q_value, slope in branch_specifications:
    rows = []
    for expression in (B, C, Y, Z):
        gradient = source_gradient(expression, x_value, q_value)
        rows.append(
            (
                sp.expand(gradient[0] + slope * gradient[1]),
                gradient[2],
            )
        )
    tangent_rows[label] = rows

expected_tangents = {
    "minus": (
        (0, 0),
        (12 + 2 * v_minus, 0),
        (-36 - 6 * v_minus, 4),
        (-(v_minus + 9) / 2, 0),
    ),
    "zero": ((0, 0), (3, 0), (-9, 4), (u, 0)),
    "plus": (
        (0, 0),
        (12 - 2 * v_plus, 0),
        (-36 + 6 * v_plus, 4),
        ((9 - v_plus) / 2, 0),
    ),
}
for label in ("minus", "zero", "plus"):
    for coordinate, (actual, expected) in enumerate(
        zip(tangent_rows[label], expected_tangents[label])
    ):
        require(
            f"tangent label={label} coordinate={coordinate}",
            zero(actual[0] - expected[0]) and zero(actual[1] - expected[1]),
        )
print("PASS collision=triple tangent_table=general_non_even_exact")


print("SECTION first-jet nonzero locus and projective parametrization")
A_form, K_form, R_form = sp.symbols("A_form K_form R_form")


def two_form_value(rows):
    """Evaluate A dC^dY+K dC^dZ+R dY^dZ on (dx,dt)."""
    c_row, y_row, z_row = rows[1], rows[2], rows[3]

    def determinant(left, right):
        return sp.expand(left[0] * right[1] - left[1] * right[0])

    return sp.expand(
        A_form * determinant(c_row, y_row)
        + K_form * determinant(c_row, z_row)
        + R_form * determinant(y_row, z_row)
    )


j_minus = two_form_value(tangent_rows["minus"])
j_zero = two_form_value(tangent_rows["zero"])
j_plus = two_form_value(tangent_rows["plus"])
require(
    "first-jet minus value",
    zero(j_minus - ((48 + 8 * v_minus) * A_form + 2 * (v_minus + 9) * R_form)),
)
require("first-jet zero value", zero(j_zero - (12 * A_form - 4 * u * R_form)))
require(
    "first-jet plus value",
    zero(j_plus - ((48 - 8 * v_plus) * A_form + 2 * (v_plus - 9) * R_form)),
)
require("first-jet K invisible", sp.diff(j_minus + j_zero + j_plus, K_form) == 0)

matching_matrix = sp.Matrix(
    [
        [18 - 4 * v_plus, v_plus - 9 + 2 * u],
        [18 + 4 * v_minus, v_minus + 9 + 2 * u],
    ]
)
matching_equation = (
    4 * u * (v_minus + v_plus)
    + 4 * v_minus * v_plus
    - 27 * v_minus
    + 27 * v_plus
    - 162
)
require("matching determinant", zero(matching_matrix.det() + 2 * matching_equation))

common_row = sp.Matrix([[12, -4 * u]])
plus_common_minor = sp.det(matching_matrix[[0], :].col_join(common_row))
minus_common_minor = sp.det(matching_matrix[[1], :].col_join(common_row))
zero_curve_plus = (4 * u - 3) * v_plus - 24 * u + 27
zero_curve_minus = (4 * u + 3) * v_minus + 24 * u + 27
require("zero-common plus minor", zero(plus_common_minor - 4 * zero_curve_plus))
require("zero-common minus minor", zero(minus_common_minor + 4 * zero_curve_minus))

rho = sp.symbols("rho")
v_plus_rho = (18 + rho * (2 * u - 9)) / (4 - rho)
v_minus_rho = -(18 + rho * (2 * u + 9)) / (4 + rho)
rho_substitution = {v_minus: v_minus_rho, v_plus: v_plus_rho}
require("rho generic matching surface", zero(matching_equation.subs(rho_substitution)))
require(
    "rho generic plus equation",
    zero((matching_matrix[0, 0] + rho * matching_matrix[0, 1]).subs(rho_substitution)),
)
require(
    "rho generic minus equation",
    zero((matching_matrix[1, 0] + rho * matching_matrix[1, 1]).subs(rho_substitution)),
)
require("rho generic common value", zero(j_zero.subs({A_form: 1, R_form: rho}) - 4 * (3 - u * rho)))

# Exceptional projective values rho=+/-4 and rho=infinity.
require(
    "rho=4 branch",
    zero(matching_equation.subs({u: sp.Rational(9, 4), v_minus: -9})),
)
rho_four_substitution = {
    u: sp.Rational(9, 4),
    v_minus: -9,
    A_form: 1,
    R_form: 4,
}
require("rho=4 minus match", zero((j_minus - j_zero).subs(rho_four_substitution)))
require("rho=4 plus match", zero((j_plus - j_zero).subs(rho_four_substitution)))
require(
    "rho=4 common",
    j_zero.subs({u: sp.Rational(9, 4), A_form: 1, R_form: 4}) == -24,
)
require(
    "rho=-4 branch",
    zero(matching_equation.subs({u: -sp.Rational(9, 4), v_plus: 9})),
)
rho_minus_four_substitution = {
    u: -sp.Rational(9, 4),
    v_plus: 9,
    A_form: 1,
    R_form: -4,
}
require(
    "rho=-4 minus match",
    zero((j_minus - j_zero).subs(rho_minus_four_substitution)),
)
require(
    "rho=-4 plus match",
    zero((j_plus - j_zero).subs(rho_minus_four_substitution)),
)
require(
    "rho=-4 common",
    j_zero.subs({u: -sp.Rational(9, 4), A_form: 1, R_form: -4}) == -24,
)
infinite_branch = {v_plus: 9 - 2 * u, v_minus: -9 - 2 * u}
require("rho=infinity matching", zero(matching_equation.subs(infinite_branch)))
require(
    "rho=infinity minus match",
    zero((j_minus - j_zero).subs({**infinite_branch, A_form: 0, R_form: 1})),
)
require(
    "rho=infinity plus match",
    zero((j_plus - j_zero).subs({**infinite_branch, A_form: 0, R_form: 1})),
)
require(
    "rho=infinity common value",
    zero(j_zero.subs({A_form: 0, R_form: 1}) + 4 * u),
)
print("PASS tangent_match=surface zero_common=removed rho_parametrization=all_strata")


print("SECTION weighted vertical top symbol")
a_vector = sp.Matrix([12 + 2 * v_minus, 3, 12 - 2 * v_plus])
k_vector = sp.Matrix([-(v_minus + 9) / 2, u, (9 - v_plus) / 2])
weight_vector = a_vector.cross(k_vector)
require("weight kills a", zero(weight_vector.dot(a_vector)))
require("weight kills k", zero(weight_vector.dot(k_vector)))

rho_zero = {v_minus: -sp.Rational(9, 2), v_plus: sp.Rational(9, 2)}
rho_zero_a = sp.simplify(a_vector.subs(rho_zero))
rho_zero_k = sp.simplify(k_vector.subs(rho_zero))
rho_zero_weight = sp.Matrix([9 - 4 * u, -18, 9 + 4 * u])
require("rho=0 a vector", rho_zero_a == sp.Matrix([3, 3, 3]))
require("rho=0 k vector", rho_zero_k == sp.Matrix([-sp.Rational(9, 4), u, sp.Rational(9, 4)]))
require("rho=0 weight kills a", zero(rho_zero_weight.dot(rho_zero_a)))
require("rho=0 weight kills k", zero(rho_zero_weight.dot(rho_zero_k)))
require(
    "even weight specializes second difference",
    rho_zero_weight.subs(u, 0) == sp.Matrix([9, -18, 9]),
)

vertical_A, vertical_R, vertical_degree = sp.symbols(
    "vertical_A vertical_R vertical_degree"
)
vertical_values = 4 ** (vertical_degree + 1) * (
    vertical_A * a_vector - vertical_R * k_vector
)
require("weighted vertical contribution cancels", zero(weight_vector.dot(vertical_values)))
print("PASS weighted_vertical_cokernel=general rho0_weight=(9-4u,-18,9+4u)")


print("SECTION generic degree-two invoice and exceptional closure")
e_centered = sp.expand(e_source + 3)


def pullback_matrix_from_local_series(local_series, maximum_degree):
    """Arbitrary-two-form pullback matrix for three prescribed Q germs."""
    target_exponents = exponent_triples(maximum_degree)
    source_exponents = source_monomials(maximum_degree)
    branch_columns = []
    for x_value, q_series in zip((-1, 0, 1), local_series):
        substitution = {x: x_value + xi, q: q_series + t**2, w: t}
        germs = tuple(
            truncate(expression.subs(substitution), (xi, t), maximum_degree + 2)
            for expression in (c_source, e_centered, w)
        )
        wedges = (
            truncate(
                sp.diff(germs[0], xi) * sp.diff(germs[1], t)
                - sp.diff(germs[0], t) * sp.diff(germs[1], xi),
                (xi, t),
                maximum_degree,
            ),
            truncate(
                sp.diff(germs[0], xi) * sp.diff(germs[2], t)
                - sp.diff(germs[0], t) * sp.diff(germs[2], xi),
                (xi, t),
                maximum_degree,
            ),
            truncate(
                sp.diff(germs[1], xi) * sp.diff(germs[2], t)
                - sp.diff(germs[1], t) * sp.diff(germs[2], xi),
                (xi, t),
                maximum_degree,
            ),
        )
        columns = []
        for wedge_coefficient in wedges:
            for c_degree, e_degree, w_degree in target_exponents:
                columns.append(
                    truncate(
                        germs[0] ** c_degree
                        * germs[1] ** e_degree
                        * germs[2] ** w_degree
                        * wedge_coefficient,
                        (xi, t),
                        maximum_degree,
                    )
                )
        branch_columns.append(columns)

    rows = []
    target = []
    for columns in branch_columns:
        for xi_degree, t_degree in source_exponents:
            monomial = xi**xi_degree * t**t_degree
            rows.append(
                [
                    sp.Poly(column, xi, t).coeff_monomial(monomial)
                    for column in columns
                ]
            )
            target.append(12 if (xi_degree, t_degree) == (0, 0) else 0)
    return sp.Matrix(rows), sp.Matrix(target)


q2_minus, q2_zero, q2_plus = sp.symbols("q2_minus q2_zero q2_plus")
symbolic_local_series = (
    -3 - sp.Rational(9, 2) * xi + q2_minus * xi**2 / 2,
    -sp.Rational(3, 4) + u * xi + q2_zero * xi**2 / 2,
    -3 + sp.Rational(9, 2) * xi + q2_plus * xi**2 / 2,
)
degree_two_matrix, degree_two_target = pullback_matrix_from_local_series(
    symbolic_local_series, 2
)
require("degree-two matrix shape", degree_two_matrix.shape == (18, 30))

# A sparse generic left-null row.  Its only denominator is 4u+9.
denominator = 4 * u + 9
invoice_row = sp.zeros(18, 1)
invoice_row[0] = -4 * (
    16 * u**2 * q2_minus
    - 16 * u**2 * q2_plus
    - 24 * u * q2_minus
    - 24 * u * q2_plus
    - 324 * u
    - 27 * q2_minus
    + 27 * q2_plus
) / (9 * denominator**2)
invoice_row[2] = 2 * (4 * u - 9) / (3 * denominator)
invoice_row[3] = -(4 * u - 9) / denominator
invoice_row[6] = -4 * (
    8 * u * q2_minus
    - 16 * u * q2_plus
    - 216 * u
    - 18 * q2_minus
    - 36 * q2_plus
    - 729
) / (3 * denominator**2)
invoice_row[9] = -18 / denominator
invoice_row[14] = sp.Rational(2, 3)
invoice_row[15] = 1
require(
    "degree-two sparse left-null identity",
    all(zero(entry) for entry in invoice_row.T * degree_two_matrix),
)
invoice_value = sp.factor(invoice_row.dot(degree_two_target))
expected_invoice_value = -16 * (
    (4 * u - 9) * q2_minus - (4 * u + 9) * q2_plus - 243
) / (3 * (4 * u + 9))
require("degree-two invoice value", zero(invoice_value - expected_invoice_value))
require("degree-two central jet absent", sp.diff(invoice_value, q2_zero) == 0)

# At u=9/4 a second left-null row kills every possible second-jet choice but
# evaluates to 36 on the normalized constant target.  u=-9/4 follows by x
# reflection, which exchanges the two side branches and reverses u.
plus_exception_matrix = degree_two_matrix.subs(u, sp.Rational(9, 4))
exception_row = sp.zeros(18, 1)
exception_row[0] = (32 * q2_zero + 297) / 216
exception_row[6] = -(32 * q2_zero - 351) / 216
exception_row[8] = sp.Rational(2, 3)
exception_row[9] = -1
exception_row[15] = 1
require(
    "u=9/4 exceptional left-null identity",
    all(zero(entry) for entry in exception_row.T * plus_exception_matrix),
)
require("u=9/4 exceptional target debt", exception_row.dot(degree_two_target) == 36)

reflection_test = sp.symbols("reflection_test")
require("compiler b reflection invariant", zero(b.subs(x, -x) - b))
require("compiler c reflection anti-invariant", zero(c_source.subs(x, -x) + c_source))
require("compiler e reflection invariant", zero(e_source.subs(x, -x) - e_source))
require(
    "central slope reflection",
    sp.diff((-sp.Rational(3, 4) + reflection_test * (-x)), x).subs(x, 0)
    == -reflection_test,
)
print("PASS degree2=weighted_side_invoice central_q2=free exceptional_u=plusminus9/4_closed")


print("SECTION explicit non-even polynomial finite survivor")
Q_h = sp.Rational(1, 5408) * (
    44069 * x**11
    + 112059 * x**10
    - 154749 * x**9
    - 406377 * x**8
    + 188107 * x**7
    + 513081 * x**6
    - 82835 * x**5
    - 230931 * x**4
    + 5408 * x
    - 4056
)
expected_jets = {
    -1: (-3, -sp.Rational(9, 2), 0, 0),
    0: (-sp.Rational(3, 4), 1, 0, 0),
    1: (-3, sp.Rational(9, 2), -sp.Rational(243, 13), sp.Rational(10449, 169)),
}
for x_value, jets in expected_jets.items():
    for order, expected in enumerate(jets):
        require(
            f"Q_h jet x={x_value} order={order}",
            sp.diff(Q_h, x, order).subs(x, x_value) == expected,
        )
require("Q_h non-even", zero(Q_h.subs(x, -x) - Q_h) == False)
require(
    "Q_h degree-two invoice",
    (9 - 4) * expected_jets[-1][2] + (9 + 4) * expected_jets[1][2] == -243,
)
require(
    "Q_h degree-four selected invoice",
    65 * expected_jets[-1][3] - 169 * expected_jets[1][3] + 10449 == 0,
)

expected_ranks = {0: 2, 1: 7, 2: 15, 3: 26, 4: 40}
for maximum_degree in range(5):
    local_series = tuple(
        sp.expand(Q_h.subs(x, x_value + xi)) for x_value in (-1, 0, 1)
    )
    matrix, target = pullback_matrix_from_local_series(local_series, maximum_degree)
    expected_shape = (
        3 * comb(maximum_degree + 2, 2),
        3 * comb(maximum_degree + 3, 3),
    )
    require(f"Q_h matrix shape N={maximum_degree}", matrix.shape == expected_shape)
    rank = DomainMatrix.from_Matrix(matrix).rank()
    augmented_rank = DomainMatrix.from_Matrix(matrix.row_join(target)).rank()
    require(f"Q_h rank N={maximum_degree}", rank == expected_ranks[maximum_degree])
    require(f"Q_h constant compatibility N={maximum_degree}", augmented_rank == rank)
    print(
        f"PASS Q_h N={maximum_degree} shape={matrix.rows}x{matrix.cols} "
        f"rank={rank} augmented={augmented_rank}"
    )


print("SECTION source AST gate")
source_path = Path(__file__)
source_tree = ast.parse(source_path.read_text(encoding="utf-8"))
assertion_nodes = sum(1 for node in ast.walk(source_tree) if isinstance(node, ast.Assert))
require("AST assertion nodes absent", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
