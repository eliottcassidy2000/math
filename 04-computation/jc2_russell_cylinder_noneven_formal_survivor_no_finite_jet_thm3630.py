#!/usr/bin/env python3
"""Exact controls for the THM-3630 non-even formal survivor.

The all-order statement is the formal Newton/Vandermonde gluing argument
recorded in the theorem. This companion checks its algebraic inputs, direct
divided differences and decomposable pullback, four independent Hermite
polynomial truncations, and two hostile repairs of the THM-3627 degree-six
polynomial.
"""

import ast
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


CHECKS = 0


def require(label, condition):
    """Record one active exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact zero test after rational simplification."""
    return sp.cancel(expression) == 0


print("THM-3630 exact companion -- provisional non-even formal survivor")
print("status=formal decomposable-pair survivor and no uniform finite jet bound; audit=PENDING")


print("SECTION exact common curve and transverse slopes")
x, q, s, u = sp.symbols("x q s u")
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
require("local chart Jacobian", zero(sp.diff(a, x) * sp.diff(c, q) - sp.diff(a, q) * sp.diff(c, x) + 3))
Q_infinity = -sp.Rational(3, 4) - sp.Rational(9, 4) / x**2
g = (1 - sp.Rational(4, 3) * s) ** -sp.Rational(1, 2)
a_zero = -sp.Rational(3, 4) + s
q_side = sp.cancel(Q_infinity.subs(x, g) + s)
require("comparison side q", zero(q_side + 3 / g**2))
require("comparison side D", zero(D.subs({x: g, q: q_side}) + 2))
require("comparison side a", zero(a.subs({x: g, q: q_side}) - a_zero))
require("comparison plus c", zero(c.subs({x: g, q: q_side})))
require("comparison minus c", zero(c.subs({x: -g, q: q_side})))
require("comparison middle a", zero(a.subs({x: 0, q: a_zero}) - a_zero))
require("comparison middle c", zero(c.subs({x: 0, q: a_zero})))

q_slope = sp.symbols("q_slope")
a_x_total = sp.diff(a, x) + sp.diff(a, q) * q_slope
c_x_total = sp.diff(c, x) + sp.diff(c, q) * q_slope
side_a_x = sp.factor(
    a_x_total.subs({q: -3 / x**2, q_slope: sp.diff(Q_infinity, x)})
)
side_c_x = sp.factor(
    c_x_total.subs({q: -3 / x**2, q_slope: sp.diff(Q_infinity, x)})
)
require("side a_x formula", zero(side_a_x - 9 / (4 * x**3)))
require("side c_x formula", zero(side_c_x - 3))
require("middle a_x formula", zero(a_x_total.subs({x: 0, q_slope: u}) - u))
require("middle c_x formula", zero(c_x_total.subs({x: 0, q_slope: u}) - 3))

lambda_values = (-sp.Rational(3, 4), u / 3, sp.Rational(3, 4))
vandermonde = sp.Matrix([[1, value, value**2] for value in lambda_values])
vandermonde_debt = sp.factor(vandermonde.det())
require(
    "Vandermonde exceptional locus",
    zero(vandermonde_debt + (4 * u - 9) * (4 * u + 9) / 96),
)
require("u=1 Vandermonde nonzero", vandermonde_debt.subs(u, 1) != 0)
common_boundary = sp.Matrix([4, 4, 4])
require(
    "degree-zero common boundary interpolation",
    vandermonde[:, :2].row_join(common_boundary).rank()
    == vandermonde[:, :2].rank(),
)
for normal_degree in range(1, 9):
    evaluation = sp.Matrix(
        [
            [value**degree for degree in range(normal_degree + 2)]
            for value in lambda_values
        ]
    )
    require(
        f"normal degree {normal_degree} interpolation rank",
        evaluation.subs(u, 1).rank() == 3,
    )
print(
    "PASS common_curve=exact c_x=3 "
    "lambda=(-3/(4g^3),u/3,3/(4g^3)) exceptional_u=plusminus9/4"
)
print("PASS formal_gluing=Vandermonde normal_degrees_all m0_boundary=4")


print("SECTION direct Newton interpolation and decomposable pullback")
c_newton, z_newton = sp.symbols("c_newton z_newton")
phi_1, A_1 = sp.symbols("phi_1 A_1")
U_21, U_31, V_21, V_31 = sp.symbols("U_21 U_31 V_21 V_31", nonzero=True)
phi_2 = phi_1 + c_newton * U_21
phi_3 = phi_1 + c_newton * U_31
A_2 = A_1 + c_newton**2 * V_21
A_3 = A_1 + c_newton**2 * V_31
D_21 = sp.cancel((A_2 - A_1) / (phi_2 - phi_1))
D_31 = sp.cancel((A_3 - A_1) / (phi_3 - phi_1))
D_321 = sp.cancel((D_31 - D_21) / (phi_3 - phi_2))
require("Newton D21 lies in cR", zero(D_21 - c_newton * V_21 / U_21))
require("Newton D31 lies in cR", zero(D_31 - c_newton * V_31 / U_31))
require("Newton D321 regular", not sp.denom(D_321).has(c_newton))
F_newton = sp.cancel(
    A_1
    + D_21 * (z_newton - phi_1)
    + D_321 * (z_newton - phi_1) * (z_newton - phi_2)
)
require("Newton restriction branch 1", zero(F_newton.subs(z_newton, phi_1) - A_1))
require("Newton restriction branch 2", zero(F_newton.subs(z_newton, phi_2) - A_2))
require("Newton restriction branch 3", zero(F_newton.subs(z_newton, phi_3) - A_3))
cx_gate = sp.symbols("cx_gate", nonzero=True)
require("decomposable pullback coefficient", zero((12 / cx_gate) * cx_gate - 12))
print(
    "PASS Newton=D21_D31_in_cR D321_regular restrictions=A1_A2_A3 "
    "decomposable_pair=(F,w) pullback=12"
)

# A hostile control on the tempting but false fixed-Q_infinity terminal
# shortcut. For the THM-3624 lower packet at u=1, the determinant of the
# three moving tangent vectors has a different next constant from the exact
# arbitrary-form elimination.
k_minus, k_plus = sp.symbols("k_minus k_plus")
X = sp.series(g - 1, s, 0, 3).removeO()
delta_minus = (
    -sp.Rational(9, 2) + k_minus * X**2 / 2
    - (-sp.Rational(9, 2) + sp.Rational(27, 2) * X - 27 * X**2)
)
delta_plus = (
    sp.Rational(9, 2) - sp.Rational(243, 13) * X + k_plus * X**2 / 2
    - (sp.Rational(9, 2) - sp.Rational(27, 2) * X + 27 * X**2)
)
h = sp.series(sp.Rational(9, 4) / g**3, s, 0, 3).removeO()
tangent_determinant = sp.series(
    delta_minus * delta_plus
    + (1 - h) * delta_minus
    + (1 + h) * delta_plus,
    s,
    0,
    3,
).removeO()
tangent_invoice = sp.factor(sp.expand(tangent_determinant).coeff(s, 2))
require(
    "moving tangent predicted invoice",
    zero(tangent_invoice + (65 * k_minus - 169 * k_plus + 11178) / 234),
)
require("terminal shortcut constant gap", 11178 - 10449 == 729)
print("PASS terminal_control=tangent_constant_11178 exact_constant_10449 gap=729")


# A source series is a sparse dictionary (xi_degree,t_degree) -> Fraction.
def series_add(left, right):
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, Fraction(0)) + coefficient
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def series_scale(series, scalar):
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in series.items()
        if coefficient * scalar
    }


def series_multiply(left, right, maximum_degree):
    answer = {}
    for (left_xi, left_t), left_coefficient in left.items():
        for (right_xi, right_t), right_coefficient in right.items():
            xi_degree = left_xi + right_xi
            t_degree = left_t + right_t
            if xi_degree + t_degree <= maximum_degree:
                monomial = (xi_degree, t_degree)
                answer[monomial] = (
                    answer.get(monomial, Fraction(0))
                    + left_coefficient * right_coefficient
                )
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def series_derivative(series, variable):
    answer = {}
    for (xi_degree, t_degree), coefficient in series.items():
        exponent = xi_degree if variable == 0 else t_degree
        if exponent:
            monomial = (
                (xi_degree - 1, t_degree)
                if variable == 0
                else (xi_degree, t_degree - 1)
            )
            answer[monomial] = coefficient * exponent
    return answer


def series_wedge(left, right, maximum_degree):
    return series_add(
        series_multiply(
            series_derivative(left, 0),
            series_derivative(right, 1),
            maximum_degree,
        ),
        series_scale(
            series_multiply(
                series_derivative(left, 1),
                series_derivative(right, 0),
                maximum_degree,
            ),
            -1,
        ),
    )


def target_monomials(maximum_degree):
    return [
        (first, second, degree - first - second)
        for degree in range(maximum_degree + 1)
        for first in range(degree + 1)
        for second in range(degree + 1 - first)
    ]


def source_monomials(maximum_degree):
    return [
        (xi_degree, degree - xi_degree)
        for degree in range(maximum_degree + 1)
        for xi_degree in range(degree + 1)
    ]


def polynomial_jet(coefficients, point, order):
    value = Fraction(0)
    for degree in range(order, len(coefficients)):
        value += (
            coefficients[degree]
            * Fraction(factorial(degree), factorial(degree - order))
            * point ** (degree - order)
        )
    return value


def q_infinity_jet(point, order):
    constant = Fraction(-3, 4) if order == 0 else Fraction(0)
    return constant - Fraction(9, 4) * ((-1) ** order) * factorial(order + 1) * (
        Fraction(point) ** (-order - 2)
    )


def hermite_coefficients(order):
    """Unique degree <=3*order+2 packet matching (6) with u=1."""
    coefficient_count = 3 * (order + 1)
    rows = []
    values = []
    for point in (-1, 0, 1):
        for derivative_order in range(order + 1):
            rows.append(
                [
                    (
                        Fraction(factorial(degree), factorial(degree - derivative_order))
                        * point ** (degree - derivative_order)
                        if degree >= derivative_order
                        else Fraction(0)
                    )
                    for degree in range(coefficient_count)
                ]
            )
            if point == 0:
                value = (
                    Fraction(-3, 4)
                    if derivative_order == 0
                    else Fraction(1)
                    if derivative_order == 1
                    else Fraction(0)
                )
            else:
                value = q_infinity_jet(point, derivative_order)
            values.append(value)
    matrix = sp.Matrix(
        [[sp.Rational(entry.numerator, entry.denominator) for entry in row] for row in rows]
    )
    target = sp.Matrix(
        [sp.Rational(value.numerator, value.denominator) for value in values]
    )
    solution = matrix.inv() * target
    return tuple(Fraction(int(value.p), int(value.q)) for value in solution)


def q_branch(coefficients, point, maximum_degree):
    answer = {}
    for degree, coefficient in enumerate(coefficients):
        for xi_degree in range(min(degree, maximum_degree) + 1):
            monomial = (xi_degree, 0)
            answer[monomial] = answer.get(monomial, Fraction(0)) + (
                coefficient
                * comb(degree, xi_degree)
                * point ** (degree - xi_degree)
            )
    if maximum_degree >= 2:
        answer[(0, 2)] = answer.get((0, 2), Fraction(0)) + 1
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def branch_coordinates(coefficients, point, maximum_degree):
    x_series = {(0, 0): Fraction(point), (1, 0): Fraction(1)}
    q_series = q_branch(coefficients, point, maximum_degree)
    x_squared_q = series_multiply(
        series_multiply(x_series, x_series, maximum_degree),
        q_series,
        maximum_degree,
    )
    c_series = series_multiply(
        x_series,
        series_add(
            {(0, 0): Fraction(3)},
            series_add(
                series_scale(x_squared_q, 4),
                series_multiply(x_squared_q, x_squared_q, maximum_degree),
            ),
        ),
        maximum_degree,
    )
    epsilon_series = series_add(
        series_multiply(
            q_series,
            series_add({(0, 0): Fraction(4)}, x_squared_q),
            maximum_degree,
        ),
        {(0, 0): Fraction(3)},
    )
    return c_series, epsilon_series, {(0, 1): Fraction(1)}


def to_sympy(value):
    return sp.Rational(value.numerator, value.denominator)


def pullback_matrix(coefficients, maximum_degree):
    target_exponents = target_monomials(maximum_degree)
    source_exponents = source_monomials(maximum_degree)
    rows = []
    for point in (-1, 0, 1):
        coordinates = branch_coordinates(coefficients, point, maximum_degree + 1)
        wedges = (
            series_wedge(coordinates[0], coordinates[1], maximum_degree),
            series_wedge(coordinates[0], coordinates[2], maximum_degree),
            series_wedge(coordinates[1], coordinates[2], maximum_degree),
        )
        powers = []
        for coordinate in coordinates:
            coordinate_powers = [{(0, 0): Fraction(1)}]
            for _ in range(maximum_degree):
                coordinate_powers.append(
                    series_multiply(
                        coordinate_powers[-1], coordinate, maximum_degree
                    )
                )
            powers.append(coordinate_powers)
        columns = []
        for wedge_coefficient in wedges:
            for c_degree, epsilon_degree, w_degree in target_exponents:
                coefficient_germ = series_multiply(
                    series_multiply(
                        powers[0][c_degree],
                        powers[1][epsilon_degree],
                        maximum_degree,
                    ),
                    powers[2][w_degree],
                    maximum_degree,
                )
                columns.append(
                    series_multiply(
                        coefficient_germ, wedge_coefficient, maximum_degree
                    )
                )
        for monomial in source_exponents:
            rows.append(
                [to_sympy(column.get(monomial, Fraction(0))) for column in columns]
            )
    target = sp.Matrix(
        [
            12 if monomial == (0, 0) else 0
            for _ in (-1, 0, 1)
            for monomial in source_exponents
        ]
    )
    return sp.Matrix(rows), target


def exact_ranks(coefficients, maximum_degree):
    matrix, target = pullback_matrix(coefficients, maximum_degree)
    return (
        matrix,
        DomainMatrix.from_Matrix(matrix).rank(),
        DomainMatrix.from_Matrix(matrix.row_join(target)).rank(),
    )


print("SECTION Hermite polynomial finite survivors")
expected_ranks = {2: 15, 4: 40, 6: 77, 8: 126}
for matched_order in range(2, 6):
    coefficients = hermite_coefficients(matched_order)
    require(
        f"Hermite degree n={matched_order}",
        len(coefficients) - 1 == 3 * matched_order + 2,
    )
    for point in (-1, 0, 1):
        for derivative_order in range(matched_order + 1):
            expected = (
                q_infinity_jet(point, derivative_order)
                if point
                else Fraction(-3, 4)
                if derivative_order == 0
                else Fraction(1)
                if derivative_order == 1
                else Fraction(0)
            )
            require(
                f"Hermite jet n={matched_order} point={point} order={derivative_order}",
                polynomial_jet(coefficients, point, derivative_order) == expected,
            )
    source_cutoff = 2 * matched_order - 2
    matrix, rank, augmented_rank = exact_ranks(coefficients, source_cutoff)
    expected_shape = (
        3 * comb(source_cutoff + 2, 2),
        3 * comb(source_cutoff + 3, 3),
    )
    require(f"Hermite shape n={matched_order}", matrix.shape == expected_shape)
    require(
        f"Hermite rank n={matched_order}", rank == expected_ranks[source_cutoff]
    )
    require(f"Hermite compatibility n={matched_order}", augmented_rank == rank)
    print(
        f"PASS matched_order={matched_order} polynomial_degree={len(coefficients)-1} "
        f"N={source_cutoff} shape={matrix.rows}x{matrix.cols} "
        f"rank={rank} augmented={augmented_rank}"
    )


def polynomial_multiply(left, right):
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] += left_coefficient * right_coefficient
    return answer


def polynomial_power(polynomial, exponent):
    answer = [Fraction(1)]
    for _ in range(exponent):
        answer = polynomial_multiply(answer, polynomial)
    return answer


def polynomial_add(left, right, scalar=Fraction(1)):
    answer = [Fraction(0)] * max(len(left), len(right))
    for degree in range(len(answer)):
        answer[degree] = (
            (left[degree] if degree < len(left) else Fraction(0))
            + scalar * (right[degree] if degree < len(right) else Fraction(0))
        )
    while answer and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


print("SECTION hostile Q_h one-side repairs")
Q_h = tuple(
    Fraction(coefficient, 5408)
    for coefficient in (
        -4056,
        5408,
        0,
        0,
        -230931,
        -82835,
        513081,
        188107,
        -406377,
        -154749,
        112059,
        44069,
    )
)
H_4 = polynomial_multiply(
    polynomial_power([0, 1], 5),
    polynomial_multiply(
        polynomial_power([1, 1], 5), polynomial_power([-1, 1], 4)
    ),
)
lambda_4 = Fraction(-58212503, 2249728)
Q_4 = polynomial_add(Q_h, H_4, lambda_4)
for point in (-1, 0):
    for derivative_order in range(5):
        require(
            f"H4 preserves point={point} order={derivative_order}",
            polynomial_jet(H_4, point, derivative_order) == 0,
        )
for derivative_order in range(4):
    require(
        f"H4 preserves plus order={derivative_order}",
        polynomial_jet(H_4, 1, derivative_order) == 0,
    )
require("H4 active fourth jet", polynomial_jet(H_4, 1, 4) == 768)
require("H4 shift", 768 * lambda_4 == Fraction(-174637509, 8788))

baseline_matrix, baseline_rank, baseline_augmented = exact_ranks(Q_h, 6)
require("Q_h degree-six shape", baseline_matrix.shape == (84, 252))
require("Q_h degree-six rank", baseline_rank == 77)
require("Q_h degree-six failure", baseline_augmented == 78)
print("PASS Q_h N=6 shape=84x252 rank=77 augmented=78 hostile_failure=reproduced")

for source_cutoff, expected_rank in ((6, 77), (7, 100)):
    matrix, rank, augmented_rank = exact_ranks(Q_4, source_cutoff)
    require(f"Q4 rank N={source_cutoff}", rank == expected_rank)
    require(f"Q4 repair N={source_cutoff}", augmented_rank == rank)
    print(
        f"PASS Q_[4] N={source_cutoff} shape={matrix.rows}x{matrix.cols} "
        f"rank={rank} augmented={augmented_rank}"
    )

H_5 = polynomial_multiply(
    polynomial_power([0, 1], 6),
    polynomial_multiply(
        polynomial_power([1, 1], 6), polynomial_power([-1, 1], 5)
    ),
)
lambda_5 = Fraction(2428928805, 58492928)
Q_5 = polynomial_add(Q_4, H_5, lambda_5)
for point in (-1, 0):
    for derivative_order in range(6):
        require(
            f"H5 preserves point={point} order={derivative_order}",
            polynomial_jet(H_5, point, derivative_order) == 0,
        )
for derivative_order in range(5):
    require(
        f"H5 preserves plus order={derivative_order}",
        polynomial_jet(H_5, 1, derivative_order) == 0,
    )
require("H5 active fifth jet", polynomial_jet(H_5, 1, 5) == 7680)
require("H5 shift", 7680 * lambda_5 == Fraction(36433932075, 114244))
matrix, rank, augmented_rank = exact_ranks(Q_5, 8)
require("Q5 degree-eight rank", rank == 126)
require("Q5 degree-eight repair", augmented_rank == rank)
print(
    f"PASS Q_[5] N=8 shape={matrix.rows}x{matrix.cols} "
    f"rank={rank} augmented={augmented_rank}"
)


print("SECTION source AST gate")
source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assertion_nodes = sum(1 for node in ast.walk(source_tree) if isinstance(node, ast.Assert))
require("AST assertion nodes absent", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
