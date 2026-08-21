#!/usr/bin/env python3
"""Exact finite closure of the THM-3624 non-even hostile polynomial.

The companion builds the enlarged arbitrary-target-two-form pullback matrices
over the rationals through total source degree six.  It verifies compatibility
through degree five and an exact rank jump at degree six, then checks a sparse
left-cokernel certificate directly against every degree-six column.
"""

import ast
from fractions import Fraction
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


# Ascending coefficients of the THM-3624 polynomial Q_h.
Q_COEFFICIENTS = tuple(
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


def polynomial_jet(point, order):
    """Exact derivative Q_h^(order)(point)."""
    value = Fraction(0)
    for degree in range(order, len(Q_COEFFICIENTS)):
        falling = 1
        for factor in range(degree - order + 1, degree + 1):
            falling *= factor
        value += Q_COEFFICIENTS[degree] * falling * point ** (degree - order)
    return value


# A truncated source series is a dict (xi_degree,t_degree) -> Fraction.
def series_add(left, right):
    """Add two sparse source series."""
    answer = dict(left)
    for monomial, coefficient in right.items():
        answer[monomial] = answer.get(monomial, Fraction(0)) + coefficient
    return {monomial: coefficient for monomial, coefficient in answer.items() if coefficient}


def series_scale(series, scalar):
    """Scale a sparse source series."""
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in series.items()
        if coefficient * scalar
    }


def series_multiply(left, right, maximum_degree):
    """Multiply and truncate by total source degree."""
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
    """Differentiate in variable 0=xi or 1=t."""
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
    """Coefficient of dleft wedge dright on dxi wedge dt."""
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
    """Exponents of c^a epsilon^b w^d in stable degree order."""
    return [
        (first, second, degree - first - second)
        for degree in range(maximum_degree + 1)
        for first in range(degree + 1)
        for second in range(degree + 1 - first)
    ]


def source_monomials(maximum_degree):
    """Exponents of xi^a t^b in stable degree order."""
    return [
        (xi_degree, degree - xi_degree)
        for degree in range(maximum_degree + 1)
        for xi_degree in range(degree + 1)
    ]


def q_branch(point, maximum_degree):
    """Truncated Q_h(point+xi)+t^2."""
    answer = {}
    for degree, coefficient in enumerate(Q_COEFFICIENTS):
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


def branch_coordinates(point, maximum_degree):
    """Local germs (c,e+3,w) for x=point+xi and q=Q_h(x)+t^2."""
    x_series = {(0, 0): Fraction(point), (1, 0): Fraction(1)}
    q_series = q_branch(point, maximum_degree)
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
    """Convert a Fraction to an exact SymPy rational."""
    return sp.Rational(value.numerator, value.denominator)


def pullback_matrix(maximum_degree):
    """Build P_N and the normalized common constant column tau_N."""
    target_exponents = target_monomials(maximum_degree)
    source_exponents = source_monomials(maximum_degree)
    rows = []

    for point in (-1, 0, 1):
        coordinates = branch_coordinates(point, maximum_degree + 1)
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


print("THM-3627 exact companion -- provisional non-even hostile degree-six closure")
print("status=finite exact arbitrary-two-form obstruction; audit=PENDING")


print("SECTION hostile polynomial and collision jets")
expected_jets = {
    -1: (
        Fraction(-3),
        Fraction(-9, 2),
        Fraction(0),
        Fraction(0),
        Fraction(2005977, 338),
        Fraction(-28734945, 338),
        Fraction(89757090, 169),
        Fraction(-326387565, 338),
    ),
    0: (
        Fraction(-3, 4),
        Fraction(1),
        Fraction(0),
        Fraction(0),
        Fraction(-692793, 676),
        Fraction(-1242525, 676),
        Fraction(23088645, 338),
        Fraction(59253705, 338),
    ),
    1: (
        Fraction(-3),
        Fraction(9, 2),
        Fraction(-243, 13),
        Fraction(10449, 169),
        Fraction(5647875, 338),
        Fraction(99668685, 338),
        Fraction(480254085, 169),
        Fraction(6097132755, 338),
    ),
}
for point, jets in expected_jets.items():
    for order, expected in enumerate(jets):
        require(
            f"Q_h jet point={point} order={order}",
            polynomial_jet(point, order) == expected,
        )
require("Q_h is non-even", Q_COEFFICIENTS[1] != 0)
print("PASS Q_h=non_even collision_jets=orders_0_through_7_exact")


print("SECTION arbitrary-target-two-form rank staircase")
expected_ranks = {0: 2, 1: 7, 2: 15, 3: 26, 4: 40, 5: 57, 6: 77}
expected_augmented_ranks = {0: 2, 1: 7, 2: 15, 3: 26, 4: 40, 5: 57, 6: 78}
degree_six_matrix = None
degree_six_target = None
for maximum_degree in range(7):
    matrix, target = pullback_matrix(maximum_degree)
    expected_shape = (
        3 * comb(maximum_degree + 2, 2),
        3 * comb(maximum_degree + 3, 3),
    )
    rank = DomainMatrix.from_Matrix(matrix).rank()
    augmented_rank = DomainMatrix.from_Matrix(matrix.row_join(target)).rank()
    require(f"matrix shape N={maximum_degree}", matrix.shape == expected_shape)
    require(f"matrix rank N={maximum_degree}", rank == expected_ranks[maximum_degree])
    require(
        f"augmented rank N={maximum_degree}",
        augmented_rank == expected_augmented_ranks[maximum_degree],
    )
    print(
        f"PASS N={maximum_degree} shape={matrix.rows}x{matrix.cols} "
        f"rank={rank} augmented={augmented_rank}"
    )
    if maximum_degree == 6:
        degree_six_matrix = matrix
        degree_six_target = target

require(
    "positive survival control through N=5",
    expected_ranks[5] == expected_augmented_ranks[5],
)
require(
    "hostile degree-six rank jump",
    expected_augmented_ranks[6] == expected_ranks[6] + 1,
)


print("SECTION sparse degree-six left-cokernel certificate")
certificate_entries = {
    0: Fraction(242002972, 2313441),
    2: Fraction(-319114, 59319),
    3: Fraction(6189, 2197),
    5: Fraction(1504, 1521),
    7: Fraction(-412, 507),
    9: Fraction(-40, 351),
    10: Fraction(36, 169),
    12: Fraction(20, 117),
    16: Fraction(-10, 39),
    21: Fraction(5, 13),
    28: Fraction(141280234, 257049),
    30: Fraction(-16512, 2197),
    31: Fraction(-6189, 2197),
    35: Fraction(-216, 169),
    38: Fraction(-36, 169),
    49: Fraction(-18, 13),
    58: Fraction(13016, 4563),
    61: Fraction(128, 117),
    63: Fraction(32, 39),
    65: Fraction(8, 27),
    68: Fraction(4, 9),
    72: Fraction(2, 3),
    77: Fraction(1),
}
certificate = sp.zeros(degree_six_matrix.rows, 1)
for row, coefficient in certificate_entries.items():
    certificate[row] = to_sympy(coefficient)
require("certificate support size", len(certificate_entries) == 23)
require(
    "certificate annihilates every degree-six column",
    all(entry == 0 for entry in certificate.T * degree_six_matrix),
)
certificate_debt = sp.factor(certificate.dot(degree_six_target))
require(
    "certificate target debt",
    certificate_debt == sp.Rational(465700024, 59319),
)

# An active hostile mutation: deleting the final plus-branch t^6 weight breaks
# the annihilator, so the certificate is not passing through a vacuous row.
mutated_certificate = certificate.copy()
mutated_certificate[77] = 0
require(
    "mutated certificate fails annihilation",
    any(entry != 0 for entry in mutated_certificate.T * degree_six_matrix),
)
print(
    "PASS left_cokernel=support_23 all_252_columns_killed "
    f"target_debt={certificate_debt} mutation=detected"
)


print("SECTION source AST gate")
source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assertion_nodes = sum(1 for node in ast.walk(source_tree) if isinstance(node, ast.Assert))
require("AST assertion nodes absent", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
