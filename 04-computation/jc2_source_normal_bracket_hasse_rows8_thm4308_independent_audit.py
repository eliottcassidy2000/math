#!/usr/bin/env python3
"""Standard-library independent audit for the proposed THM-4308.

This script imports no primary certificate and no computer-algebra package.
It works over a sparse polynomial ring over ``fractions.Fraction`` and checks
the fixed source-normal bracket/residual rows through row eight, the finite
depth-two/depth-three projections, and an exact gate-interior numerical jet.

The scope is deliberately finite: projected membership through t^8 is not an
all-row B_2 lift, a Keller pair, seam entry, or a consequence for JC(2).
"""

from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
from math import comb, gcd, lcm


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


PARAMETERS = (
    "Delta",
    "Phi",
    "Theta",
    "eta",
    "zeta3",
    "upsilon5",
    "xi10",
    "alpha11",
    "beta11",
    "U",
    "W",
    "Z",
)
PARAMETER_INDEX = {name: index for index, name in enumerate(PARAMETERS)}
ZERO_EXPONENT = (0,) * len(PARAMETERS)


class ParameterPolynomial:
    """A sparse element of Q[PARAMETERS]."""

    def __init__(self, terms: dict[tuple[int, ...], F] | None = None) -> None:
        self.terms = {
            monomial: F(coefficient)
            for monomial, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def coerce(value: object) -> "ParameterPolynomial":
        if isinstance(value, ParameterPolynomial):
            return value
        return ParameterPolynomial({ZERO_EXPONENT: F(value)})

    def __add__(self, other: object) -> "ParameterPolynomial":
        other_poly = self.coerce(other)
        result = self.terms.copy()
        for monomial, coefficient in other_poly.terms.items():
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return ParameterPolynomial(result)

    __radd__ = __add__

    def __neg__(self) -> "ParameterPolynomial":
        return ParameterPolynomial(
            {monomial: -coefficient for monomial, coefficient in self.terms.items()}
        )

    def __sub__(self, other: object) -> "ParameterPolynomial":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "ParameterPolynomial":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "ParameterPolynomial":
        other_poly = self.coerce(other)
        result: dict[tuple[int, ...], F] = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other_poly.terms.items():
                monomial = tuple(
                    left + right
                    for left, right in zip(left_monomial, right_monomial)
                )
                result[monomial] = (
                    result.get(monomial, F(0))
                    + left_coefficient * right_coefficient
                )
        return ParameterPolynomial(result)

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | F) -> "ParameterPolynomial":
        scalar_fraction = F(scalar)
        require(scalar_fraction != 0, "nonzero scalar divisor")
        return ParameterPolynomial(
            {
                monomial: coefficient / scalar_fraction
                for monomial, coefficient in self.terms.items()
            }
        )

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __eq__(self, other: object) -> bool:
        return self.terms == self.coerce(other).terms

    def substitute(
        self, parameter: str | int, replacement: object
    ) -> "ParameterPolynomial":
        index = (
            PARAMETER_INDEX[parameter]
            if isinstance(parameter, str)
            else parameter
        )
        replacement_poly = self.coerce(replacement)
        result = ParameterPolynomial()
        for monomial, coefficient in self.terms.items():
            exponent = monomial[index]
            base_monomial = list(monomial)
            base_monomial[index] = 0
            term = ParameterPolynomial({tuple(base_monomial): coefficient})
            for _ in range(exponent):
                term = term * replacement_poly
            result = result + term
        return result

    def substitute_many(
        self, replacements: dict[str, object]
    ) -> "ParameterPolynomial":
        result = self
        for parameter in PARAMETERS:
            if parameter in replacements:
                result = result.substitute(parameter, replacements[parameter])
        return result

    def as_fraction(self) -> F:
        require(
            set(self.terms).issubset({ZERO_EXPONENT}),
            "specialization is rational",
        )
        return self.terms.get(ZERO_EXPONENT, F(0))

    def uses(self, parameter: str) -> bool:
        index = PARAMETER_INDEX[parameter]
        return any(monomial[index] for monomial in self.terms)


def parameter(name: str) -> ParameterPolynomial:
    exponent = [0] * len(PARAMETERS)
    exponent[PARAMETER_INDEX[name]] = 1
    return ParameterPolynomial({tuple(exponent): F(1)})


(
    Delta,
    Phi,
    Theta,
    eta,
    zeta3,
    upsilon5,
    xi10,
    alpha11,
    beta11,
    U,
    W,
    Z,
) = (parameter(name) for name in PARAMETERS)


XPolynomial = dict[int, ParameterPolynomial]
Rows = dict[int, XPolynomial]


def x_add(left: XPolynomial, right: XPolynomial) -> XPolynomial:
    result = left.copy()
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, ParameterPolynomial()) + coefficient
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def x_negate(poly: XPolynomial) -> XPolynomial:
    return {degree: -coefficient for degree, coefficient in poly.items()}


def x_subtract(left: XPolynomial, right: XPolynomial) -> XPolynomial:
    return x_add(left, x_negate(right))


def x_scale(poly: XPolynomial, scalar: object) -> XPolynomial:
    return {
        degree: coefficient * scalar
        for degree, coefficient in poly.items()
        if coefficient * scalar
    }


def x_multiply(left: XPolynomial, right: XPolynomial) -> XPolynomial:
    result: XPolynomial = {}
    for left_degree, left_coefficient in left.items():
        for right_degree, right_coefficient in right.items():
            degree = left_degree + right_degree
            result[degree] = (
                result.get(degree, ParameterPolynomial())
                + left_coefficient * right_coefficient
            )
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def x_derivative(poly: XPolynomial) -> XPolynomial:
    return {
        degree - 1: coefficient * degree
        for degree, coefficient in poly.items()
        if degree
    }


def x_shift(poly: XPolynomial, shift: int) -> XPolynomial:
    require(all(degree + shift >= 0 for degree in poly), "nonnegative x shift")
    return {degree + shift: coefficient for degree, coefficient in poly.items()}


def x_coefficient(poly: XPolynomial, degree: int) -> ParameterPolynomial:
    return poly.get(degree, ParameterPolynomial())


def x_constant(value: object) -> XPolynomial:
    coefficient = ParameterPolynomial.coerce(value)
    return {0: coefficient} if coefficient else {}


def x_substitute(
    poly: XPolynomial, parameter_name: str, replacement: object
) -> XPolynomial:
    result: XPolynomial = {}
    for degree, coefficient in poly.items():
        new_coefficient = coefficient.substitute(parameter_name, replacement)
        if new_coefficient:
            result[degree] = new_coefficient
    return result


def rows_substitute(
    rows: Rows, replacements: dict[str, object]
) -> Rows:
    result: Rows = {}
    for row, poly in rows.items():
        new_poly = poly
        for parameter_name in PARAMETERS:
            if parameter_name in replacements:
                new_poly = x_substitute(
                    new_poly, parameter_name, replacements[parameter_name]
                )
        result[row] = new_poly
    return result


def copy_rows(rows: Rows) -> Rows:
    return {row: poly.copy() for row, poly in rows.items()}


def rational_rref(matrix: list[list[F]]) -> tuple[list[list[F]], list[int]]:
    result = [list(map(F, row)) for row in matrix]
    row_count = len(result)
    column_count = len(result[0]) if row_count else 0
    pivots: list[int] = []
    pivot_row = 0
    for column in range(column_count):
        chosen = next(
            (row for row in range(pivot_row, row_count) if result[row][column]),
            None,
        )
        if chosen is None:
            continue
        result[pivot_row], result[chosen] = result[chosen], result[pivot_row]
        pivot = result[pivot_row][column]
        result[pivot_row] = [entry / pivot for entry in result[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not result[row][column]:
                continue
            multiplier = result[row][column]
            result[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(result[row], result[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break
    return result, pivots


def rational_rank(matrix: list[list[F]]) -> int:
    return len(rational_rref(matrix)[1])


def rational_nullspace(matrix: list[list[F]]) -> list[list[F]]:
    reduced, pivots = rational_rref(matrix)
    column_count = len(matrix[0])
    free_columns = [column for column in range(column_count) if column not in pivots]
    result: list[list[F]] = []
    for free_column in free_columns:
        vector = [F(0)] * column_count
        vector[free_column] = F(1)
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -reduced[row][free_column]
        result.append(vector)
    return result


def transpose(matrix: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*matrix)]


def primitive_integer_vector(vector: list[F]) -> list[int]:
    denominator = 1
    for entry in vector:
        denominator = lcm(denominator, entry.denominator)
    integers = [int(entry * denominator) for entry in vector]
    divisor = reduce(gcd, (abs(entry) for entry in integers if entry))
    integers = [entry // divisor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    return [-entry for entry in integers] if first_nonzero < 0 else integers


def solve_rational_matrix_with_parameter_rhs(
    matrix: list[list[F]], rhs: list[ParameterPolynomial]
) -> tuple[list[ParameterPolynomial], list[ParameterPolynomial], int]:
    row_count = len(matrix)
    column_count = len(matrix[0])
    augmented: list[list[object]] = [
        list(matrix[row]) + [rhs[row]] for row in range(row_count)
    ]
    pivots: list[tuple[int, int]] = []
    pivot_row = 0
    for column in range(column_count):
        chosen = next(
            (
                row
                for row in range(pivot_row, row_count)
                if augmented[row][column]
            ),
            None,
        )
        if chosen is None:
            continue
        augmented[pivot_row], augmented[chosen] = (
            augmented[chosen],
            augmented[pivot_row],
        )
        pivot = F(augmented[pivot_row][column])
        for index in range(column_count):
            augmented[pivot_row][index] = F(augmented[pivot_row][index]) / pivot
        augmented[pivot_row][-1] = augmented[pivot_row][-1] / pivot
        for row in range(row_count):
            if row == pivot_row or not augmented[row][column]:
                continue
            multiplier = F(augmented[row][column])
            for index in range(column_count):
                augmented[row][index] = (
                    F(augmented[row][index])
                    - multiplier * F(augmented[pivot_row][index])
                )
            augmented[row][-1] = (
                augmented[row][-1]
                - multiplier * augmented[pivot_row][-1]
            )
        pivots.append((pivot_row, column))
        pivot_row += 1
    solution = [ParameterPolynomial() for _ in range(column_count)]
    for row, column in pivots:
        solution[column] = augmented[row][-1]
    residuals = [augmented[row][-1] for row in range(pivot_row, row_count)]
    return solution, residuals, pivot_row


def proportional_ratio(
    actual: ParameterPolynomial, expected: ParameterPolynomial
) -> F | None:
    require(bool(expected), "nonzero expected equation")
    monomial = next(iter(expected.terms))
    ratio = actual.terms.get(monomial, F(0)) / expected.terms[monomial]
    return ratio if actual == expected * ratio else None


def solve_for_parameter(
    equation: ParameterPolynomial, parameter_name: str
) -> ParameterPolynomial:
    index = PARAMETER_INDEX[parameter_name]
    coefficient = F(0)
    remainder: dict[tuple[int, ...], F] = {}
    for monomial, value in equation.terms.items():
        if monomial[index]:
            require(
                monomial[index] == 1 and sum(monomial) == 1,
                f"{parameter_name} occurs with constant linear coefficient",
            )
            coefficient += value
        else:
            remainder[monomial] = value
    require(coefficient != 0, f"solve {parameter_name}")
    return -ParameterPolynomial(remainder) / coefficient


def source_rows() -> tuple[Rows, ParameterPolynomial]:
    """Expand G=-u/2+H directly from p,y monomials."""

    K = F(2848, 45) - F(7, 6) * Delta
    source_terms = (
        (-3, 1, 0),
        (F(8, 3), 2, 0),
        (F(-1376, 135), 3, 0),
        (K, 0, 2),
        (Phi, 2, 1),
        (Delta, 4, 0),
        (Theta, 1, 2),
        (eta, 3, 1),
        (zeta3, 0, 3),
        (upsilon5, 5, 0),
        (xi10, 2, 2),
        (alpha11, 4, 1),
        (beta11, 1, 3),
        (U, 6, 0),
        (W, 3, 2),
        (Z, 0, 4),
    )
    result: Rows = {1: {2: ParameterPolynomial.coerce(F(-1, 2))}}
    for coefficient, p_power, y_power in source_terms:
        coefficient_poly = ParameterPolynomial.coerce(coefficient)
        for binomial_index in range(p_power + y_power + 1):
            row = p_power + 2 * y_power + binomial_index
            degree = y_power + 2 * binomial_index
            result.setdefault(row, {})[degree] = (
                result.setdefault(row, {}).get(degree, ParameterPolynomial())
                + coefficient_poly
                * comb(p_power + y_power, binomial_index)
            )
    result = {
        row: {
            degree: coefficient
            for degree, coefficient in poly.items()
            if coefficient
        }
        for row, poly in result.items()
    }
    return result, K


def inherited_rows(K: ParameterPolynomial) -> tuple[Rows, Rows]:
    A: Rows = {
        0: {0: ParameterPolynomial.coerce(1), 2: ParameterPolynomial.coerce(F(1, 4))},
        1: {0: ParameterPolynomial.coerce(F(4, 3)), 2: ParameterPolynomial.coerce(2)},
        2: {
            0: ParameterPolynomial.coerce(F(-32, 9)),
            2: ParameterPolynomial.coerce(F(-4, 5)),
        },
        3: {
            0: ParameterPolynomial.coerce(F(2176, 135)),
            1: -Phi / 2,
            2: F(1088, 315) - F(4, 7) * K,
            4: ParameterPolynomial.coerce(F(-32, 15)),
        },
    }
    C: Rows = {
        0: {
            1: ParameterPolynomial.coerce(F(-3, 4)),
            3: ParameterPolynomial.coerce(F(-1, 8)),
        },
        1: {
            1: ParameterPolynomial.coerce(-4),
            3: ParameterPolynomial.coerce(F(-3, 2)),
        },
        2: {
            1: ParameterPolynomial.coerce(F(88, 15)),
            3: ParameterPolynomial.coerce(F(-12, 5)),
        },
        3: {
            0: F(3, 4) * Phi,
            1: F(-8128, 315) + F(6, 7) * K,
            2: F(3, 8) * Phi,
            3: F(736, 105) + F(3, 7) * K,
            5: ParameterPolynomial.coerce(F(8, 5)),
        },
    }
    return A, C


Q_POLY: XPolynomial = {
    0: ParameterPolynomial.coerce(-3),
    2: ParameterPolynomial.coerce(F(-1, 2)),
}


def bracket_lower(A: Rows, C: Rows, row: int) -> XPolynomial:
    result: XPolynomial = {}
    for index in range(1, row):
        result = x_add(
            result,
            x_add(
                x_scale(
                    x_multiply(x_derivative(A[index]), C[row - index]),
                    row - index,
                ),
                x_scale(
                    x_multiply(A[index], x_derivative(C[row - index])),
                    -index,
                ),
            ),
        )
    return result


def residual_lower(A: Rows, C: Rows, row: int) -> XPolynomial:
    result: XPolynomial = {}
    for index in range(1, row):
        result = x_add(result, x_multiply(C[index], C[row - index]))
    for first in range(row):
        for second in range(row):
            third = row - first - second
            if 0 <= third < row:
                result = x_subtract(
                    result,
                    x_multiply(x_multiply(A[first], A[second]), A[third]),
                )
    return result


def invariant_residual(A: Rows, C: Rows, row: int) -> XPolynomial:
    return x_subtract(
        residual_lower(A, C, row),
        x_scale(x_multiply(Q_POLY, bracket_lower(A, C, row)), F(1, row)),
    )


def bracket_particular(A: Rows, C: Rows, row: int) -> tuple[XPolynomial, XPolynomial]:
    determinant_target = x_scale(bracket_lower(A, C, row), F(-1, row))
    A_constant = F(4, 3) * x_coefficient(determinant_target, 0)
    A_row = x_constant(A_constant)
    remainder = x_subtract(
        determinant_target,
        x_scale(
            x_multiply(
                {
                    0: ParameterPolynomial.coerce(2),
                    2: ParameterPolynomial.coerce(1),
                },
                A_row,
            ),
            F(3, 8),
        ),
    )
    require(x_coefficient(remainder, 0) == 0, f"row {row} division by x")
    C_row = x_scale(x_shift(remainder, -1), 2)
    return A_row, C_row


def tangent_matrix(next_row: int) -> list[list[F]]:
    """Matrix of theta -> q' theta - (q/next_row) theta'."""

    result = [[F(0) for _ in range(next_row)] for _ in range(next_row + 1)]
    for degree in range(next_row):
        if degree:
            result[degree - 1][degree] += F(3 * degree, next_row)
        result[degree + 1][degree] += F(degree, 2 * next_row) - 1
    return result


def tangent_poly(coefficients: list[ParameterPolynomial]) -> XPolynomial:
    return {
        degree: coefficient
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    }


def bracket_operator_matrix(row: int) -> list[list[F]]:
    a_columns = row + 2
    c_columns = row + 3
    result = [
        [F(0) for _ in range(a_columns + c_columns)]
        for _ in range(row + 4)
    ]
    for degree in range(a_columns):
        result[degree][degree] += F(3 * row, 4)
        result[degree + 2][degree] += F(3 * row, 8)
    for degree in range(c_columns):
        result[degree + 1][a_columns + degree] += F(row, 2)
    return result


def jacobian_coefficient(A: Rows, C: Rows, power: int) -> XPolynomial:
    result: XPolynomial = {}
    for first in range(power + 2):
        second = power + 1 - first
        if first not in A or second not in C:
            continue
        result = x_add(
            result,
            x_add(
                x_scale(
                    x_multiply(x_derivative(A[first]), C[second]), second
                ),
                x_scale(
                    x_multiply(A[first], x_derivative(C[second])), -first
                ),
            ),
        )
    return result


def P_coefficient(A: Rows, C: Rows, row: int) -> XPolynomial:
    result: XPolynomial = {}
    for first in range(row + 1):
        second = row - first
        result = x_add(result, x_multiply(C[first], C[second]))
    for first in range(row + 1):
        for second in range(row - first + 1):
            third = row - first - second
            result = x_subtract(
                result,
                x_multiply(x_multiply(A[first], A[second]), A[third]),
            )
    result = x_add(result, x_scale(A[row], F(3, 4)))
    if row == 0:
        result = x_add(result, x_constant(F(1, 4)))
    return result


def check_direct_replay(A: Rows, C: Rows, G: Rows, label: str) -> None:
    for power in range(8):
        expected = x_constant(1) if power == 0 else {}
        require(
            jacobian_coefficient(A, C, power) == expected,
            f"{label} Jacobian t^{power}",
        )
    for row in range(9):
        require(
            P_coefficient(A, C, row) == G.get(row, {}),
            f"{label} residual t^{row}",
        )


def depth_projection_matrix(
    depth: int,
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]], list[list[F]]]:
    coordinates = [
        (row, degree)
        for row in range(9)
        for degree in range(row + depth + 1)
    ]
    monomials: list[tuple[int, int, int, int]] = []
    for u_power in range(depth + 1):
        for x_power in range(depth - u_power + 1):
            for y_power in range((8 - u_power) // 2 + 1):
                for p_power in range(8 - u_power - 2 * y_power + 1):
                    monomials.append((x_power, u_power, p_power, y_power))
    coordinate_index = {
        coordinate: index for index, coordinate in enumerate(coordinates)
    }
    matrix = [
        [F(0) for _ in monomials]
        for _ in coordinates
    ]
    for column, (x_power, u_power, p_power, y_power) in enumerate(monomials):
        for binomial_index in range(p_power + y_power + 1):
            row = u_power + p_power + 2 * y_power + binomial_index
            degree = (
                x_power
                + 2 * u_power
                + y_power
                + 2 * binomial_index
            )
            if row <= 8:
                matrix[coordinate_index[(row, degree)]][column] += F(
                    comb(p_power + y_power, binomial_index)
                )
    return coordinates, monomials, matrix


def row_coordinate_vector(rows: Rows, coordinates: list[tuple[int, int]]) -> list[ParameterPolynomial]:
    return [x_coefficient(rows.get(row, {}), degree) for row, degree in coordinates]


def dot_parameter(
    rational_vector: list[F], parameter_vector: list[ParameterPolynomial]
) -> ParameterPolynomial:
    return sum(
        (
            coefficient * value
            for coefficient, value in zip(rational_vector, parameter_vector)
        ),
        ParameterPolynomial(),
    )


Affine = tuple[ParameterPolynomial, tuple[F, ...]]


def affine_add(left: Affine, right: Affine) -> Affine:
    return (
        left[0] + right[0],
        tuple(a + b for a, b in zip(left[1], right[1])),
    )


def affine_scale(value: Affine, scalar: int | F) -> Affine:
    return (value[0] * scalar, tuple(scalar * entry for entry in value[1]))


def affine_equal(left: Affine, right: Affine) -> bool:
    return left[0] == right[0] and left[1] == right[1]


def terminal_affine_coefficient(
    base_row: XPolynomial, derivative: XPolynomial, degree: int
) -> Affine:
    tangent = tuple(
        x_coefficient(derivative, degree - theta_degree).as_fraction()
        if degree - theta_degree >= 0
        else F(0)
        for theta_degree in range(9)
    )
    return x_coefficient(base_row, degree), tangent


def functional_affine_value(
    functional: dict[tuple[int, int], F],
    rows: Rows,
    derivative: XPolynomial,
) -> Affine:
    result: Affine = (ParameterPolynomial(), (F(0),) * 9)
    for (row, degree), coefficient in functional.items():
        value = (
            terminal_affine_coefficient(rows[8], derivative, degree)
            if row == 8
            else (x_coefficient(rows[row], degree), (F(0),) * 9)
        )
        result = affine_add(result, affine_scale(value, coefficient))
    return result


def constant_affine(value: object) -> Affine:
    return ParameterPolynomial.coerce(value), (F(0),) * 9


def functional_annihilates_matrix(
    functional: dict[tuple[int, int], F],
    coordinates: list[tuple[int, int]],
    matrix: list[list[F]],
) -> bool:
    vector = [functional.get(coordinate, F(0)) for coordinate in coordinates]
    return all(
        sum(vector[row] * matrix[row][column] for row in range(len(vector))) == 0
        for column in range(len(matrix[0]))
    )


def specialize_rows(rows: Rows, values: dict[str, object]) -> Rows:
    return rows_substitute(rows, values)


def rational_coordinate_vector(
    rows: Rows, coordinates: list[tuple[int, int]]
) -> list[F]:
    return [
        x_coefficient(rows.get(row, {}), degree).as_fraction()
        for row, degree in coordinates
    ]


def augmented_rank(matrix: list[list[F]], vector: list[F]) -> int:
    return rational_rank(
        [row + [vector[index]] for index, row in enumerate(matrix)]
    )


def main() -> None:
    G, K = source_rows()
    A, C = inherited_rows(K)

    expected_G = {
        4: {
            0: Delta,
            1: Phi,
            2: K - F(1376, 45),
            4: ParameterPolynomial.coerce(F(8, 3)),
        },
        5: {
            0: upsilon5,
            1: eta,
            2: 4 * Delta + Theta,
            3: 3 * Phi,
            4: 2 * K - F(1376, 45),
        },
        6: {
            0: U,
            1: alpha11,
            2: 5 * upsilon5 + xi10,
            3: 4 * eta + zeta3,
            4: 6 * Delta + 3 * Theta,
            5: 3 * Phi,
            6: K - F(1376, 135),
        },
        7: {
            2: 6 * U + W,
            3: 5 * alpha11 + beta11,
            4: 10 * upsilon5 + 4 * xi10,
            5: 6 * eta + 3 * zeta3,
            6: 4 * Delta + 3 * Theta,
            7: Phi,
        },
        8: {
            4: 15 * U + 5 * W + Z,
            5: 10 * alpha11 + 4 * beta11,
            6: 10 * upsilon5 + 6 * xi10,
            7: 4 * eta + 3 * zeta3,
            8: Delta + Theta,
        },
    }
    for row in range(4, 9):
        require(G[row] == expected_G[row], f"direct source expansion G{row}")

    gradient_left = x_add(
        x_scale(x_multiply(A[0], A[0]), -3),
        x_constant(F(3, 4)),
    )
    gradient_right = x_scale(C[0], 2)
    require(
        gradient_left == x_scale(x_multiply(Q_POLY, x_derivative(C[0])), -1),
        "gradient A component",
    )
    require(
        gradient_right == x_multiply(Q_POLY, x_derivative(A[0])),
        "gradient C component",
    )

    bracket_rank_rows = []
    for row in range(4, 9):
        matrix = bracket_operator_matrix(row)
        rank = rational_rank(matrix)
        nullity = len(matrix[0]) - rank
        require(rank == row + 4, f"bracket operator row {row} surjective")
        require(nullity == row + 1, f"bracket operator row {row} kernel")
        bracket_rank_rows.append((row, len(matrix), len(matrix[0]), rank, nullity))

    E = {
        5: 2025 * upsilon5 + 9000 * Delta + 1350 * Theta + 184832,
        6: (
            200475 * U
            + 109350 * xi10
            - 5593860 * Delta
            - 529200 * Theta
            - 137763328
        ),
        7: (
            801900 * W
            + 1782000 * Delta * Delta
            + 156163200 * Delta
            + 868725 * Phi * Phi
            + 27390480 * Theta
            - 3434400 * xi10
            + 12891824128
        ),
        8: (
            21651300 * Z
            - 225022050 * Delta * Delta
            - 59073300 * Delta * Theta
            - 9512522400 * Delta
            + 34749000 * Phi * Phi
            + 39092625 * Phi * eta
            + 940522560 * Theta
            + 185376600 * xi10
            - 1112446017536
        ),
    }
    raw_factors = {
        5: F(7, 24300),
        6: F(7, 6561000),
        7: F(1, 4374000),
        8: F(1, 275562000),
    }
    expected_cokernels = {
        5: [21, 0, 14, 0, 36, 0],
        6: [77, 0, 42, 0, 84, 0, 360],
        7: [143, 0, 66, 0, 108, 0, 360, 0],
        8: [715, 0, 286, 0, 396, 0, 1080, 0, 5040],
    }
    new_parameters = {5: "upsilon5", 6: "U", 7: "W", 8: "Z"}
    response_substitutions: dict[str, ParameterPolynomial] = {}

    require(invariant_residual(A, C, 4) == G[4], "inherited G4 compatibility")
    for row in range(4, 8):
        A[row], C[row] = bracket_particular(A, C, row)
        next_row = row + 1
        base_invariant = invariant_residual(A, C, next_row)

        matrix = tangent_matrix(next_row)
        require(rational_rank(matrix) == next_row, f"M{next_row} rank")
        left_kernel = rational_nullspace(transpose(matrix))
        require(len(left_kernel) == 1, f"M{next_row} cokernel dimension")
        require(
            primitive_integer_vector(left_kernel[0])
            == expected_cokernels[next_row],
            f"M{next_row} primitive cokernel",
        )

        for theta_degree in range(next_row):
            theta_basis = {theta_degree: ParameterPolynomial.coerce(1)}
            A[row] = x_add(A[row], x_multiply(theta_basis, x_derivative(A[0])))
            C[row] = x_add(C[row], x_multiply(theta_basis, x_derivative(C[0])))
            direct_change = x_subtract(
                invariant_residual(A, C, next_row), base_invariant
            )
            matrix_change = {
                degree: ParameterPolynomial.coerce(matrix[degree][theta_degree])
                for degree in range(next_row + 1)
                if matrix[degree][theta_degree]
            }
            require(
                direct_change == matrix_change,
                f"M{next_row} direct tangent column {theta_degree}",
            )
            A[row] = x_subtract(
                A[row], x_multiply(theta_basis, x_derivative(A[0]))
            )
            C[row] = x_subtract(
                C[row], x_multiply(theta_basis, x_derivative(C[0]))
            )

        target = x_subtract(G[next_row], base_invariant)
        theta_solution, residuals, rank = solve_rational_matrix_with_parameter_rhs(
            matrix,
            [x_coefficient(target, degree) for degree in range(next_row + 1)],
        )
        require(rank == next_row, f"M{next_row} solve rank")
        require(len(residuals) == 1, f"M{next_row} one compatibility")
        ratio = proportional_ratio(residuals[0], E[next_row])
        require(ratio == raw_factors[next_row], f"E{next_row} raw factor")

        new_parameter = new_parameters[next_row]
        replacement = solve_for_parameter(residuals[0], new_parameter)
        response_substitutions[new_parameter] = replacement
        A = rows_substitute(A, {new_parameter: replacement})
        C = rows_substitute(C, {new_parameter: replacement})
        G = rows_substitute(G, {new_parameter: replacement})
        theta_solution = [
            coefficient.substitute(new_parameter, replacement)
            for coefficient in theta_solution
        ]
        theta = tangent_poly(theta_solution)
        A[row] = x_add(A[row], x_multiply(theta, x_derivative(A[0])))
        C[row] = x_add(C[row], x_multiply(theta, x_derivative(C[0])))
        require(
            invariant_residual(A, C, next_row) == G[next_row],
            f"G{next_row} compatibility after E{next_row}",
        )
        require(max(A[row], default=-1) <= row + 1, f"A{row} sharp cap")
        require(max(C[row], default=-1) <= row + 2, f"C{row} sharp cap")

    A[8], C[8] = bracket_particular(A, C, 8)
    require(max(A[8], default=-1) <= 9, "A8 sharp cap")
    require(max(C[8], default=-1) <= 10, "C8 sharp cap")
    check_direct_replay(A, C, G, "symbolic pre-Hasse")

    pre_hasse_A = copy_rows(A)
    pre_hasse_C = copy_rows(C)
    pre_hasse_G = copy_rows(G)

    coordinates_2, monomials_2, matrix_2 = depth_projection_matrix(2)
    coordinates_3, monomials_3, matrix_3 = depth_projection_matrix(3)
    rank_2 = rational_rank(matrix_2)
    rank_3 = rational_rank(matrix_3)
    left_kernel_2 = rational_nullspace(transpose(matrix_2))
    left_kernel_3 = rational_nullspace(transpose(matrix_3))
    require(
        (len(coordinates_2), len(monomials_2), rank_2, len(left_kernel_2))
        == (63, 131, 51, 12),
        "depth-two projection dimensions",
    )
    require(
        (len(coordinates_3), len(monomials_3), rank_3, len(left_kernel_3))
        == (72, 204, 63, 9),
        "depth-three projection dimensions",
    )

    A_delta = {(2, 0): F(-4), (3, 2): F(3), (4, 4): F(-2), (5, 6): F(1)}
    A_zeta = {(4, 1): F(-10), (5, 3): F(6), (6, 5): F(-3), (7, 7): F(1)}
    A_c89 = {(4, 0): F(15), (5, 2): F(-10), (6, 4): F(6), (7, 6): F(-3), (8, 8): F(1)}
    C_c89 = {(4, 1): F(15), (5, 3): F(-10), (6, 5): F(6), (7, 7): F(-3), (8, 9): F(1)}
    A_c810 = {(4, 1): F(-15), (5, 3): F(8), (6, 5): F(-3), (8, 9): F(1)}
    C_c810 = {(3, 0): F(-6), (4, 2): F(5), (5, 4): F(-4), (6, 6): F(3), (7, 8): F(-2), (8, 10): F(1)}
    for functional in (A_delta, A_zeta, A_c89, A_c810):
        require(
            functional_annihilates_matrix(functional, coordinates_2, matrix_2),
            "displayed P2 functional",
        )
    for functional in (C_c89, C_c810):
        require(
            functional_annihilates_matrix(functional, coordinates_3, matrix_3),
            "displayed P3 functional",
        )

    A_derivative = x_derivative(A[0])
    C_derivative = x_derivative(C[0])
    zero_tangent = (F(0),) * 9
    c89 = terminal_affine_coefficient(C[8], C_derivative, 9)
    c810 = terminal_affine_coefficient(C[8], C_derivative, 10)
    L_delta = functional_affine_value(A_delta, A, A_derivative)
    L_zeta = functional_affine_value(A_zeta, A, A_derivative)
    L_A_c89 = functional_affine_value(A_c89, A, A_derivative)
    L_C_c89 = functional_affine_value(C_c89, C, C_derivative)
    L_A_c810 = functional_affine_value(A_c810, A, A_derivative)
    L_C_c810 = functional_affine_value(C_c810, C, C_derivative)
    require(
        affine_equal(
            L_delta,
            constant_affine(-(15 * Delta - 896) / 45),
        ),
        "Delta Hasse identity",
    )
    require(
        affine_equal(
            L_zeta,
            constant_affine(F(3, 20) * (3 * Phi + 2 * zeta3)),
        ),
        "zeta3 Hasse identity",
    )
    require(
        affine_equal(
            L_A_c89,
            affine_scale(
                affine_add(
                    constant_affine(
                        13215 * Delta
                        + 7950 * Theta
                        + 6075 * xi10
                        - 2583808
                    ),
                    affine_scale(c89, -2475),
                ),
                F(4, 7425),
            ),
        ),
        "A c89 Hasse identity",
    )
    require(
        affine_equal(
            L_C_c89,
            affine_scale(
                affine_add(
                    constant_affine(
                        6900 * Delta
                        - 13425 * Theta
                        - 12150 * xi10
                        + 3159808
                    ),
                    affine_scale(c89, 4950),
                ),
                F(1, 4950),
            ),
        ),
        "C c89 Hasse identity",
    )
    require(
        affine_equal(
            L_A_c810,
            affine_scale(
                affine_add(
                    constant_affine(27 * Phi + 27 * eta + 45 * zeta3),
                    affine_scale(c810, -40),
                ),
                F(1, 30),
            ),
        ),
        "A c810 Hasse identity",
    )
    require(
        affine_equal(
            L_C_c810,
            affine_scale(
                affine_add(
                    constant_affine(-27 * eta - 27 * zeta3),
                    affine_scale(c810, 40),
                ),
                F(1, 40),
            ),
        ),
        "C c810 Hasse identity",
    )

    Hasse_matrix: list[list[F]] = []
    Hasse_base: list[ParameterPolynomial] = []
    A_vector = row_coordinate_vector(A, coordinates_2)
    C_vector = row_coordinate_vector(C, coordinates_3)
    for left_vector, coordinate_vector, coordinates, derivative in (
        *(
            (left_vector, A_vector, coordinates_2, A_derivative)
            for left_vector in left_kernel_2
        ),
        *(
            (left_vector, C_vector, coordinates_3, C_derivative)
            for left_vector in left_kernel_3
        ),
    ):
        Hasse_base.append(dot_parameter(left_vector, coordinate_vector))
        tangent_row: list[F] = []
        for theta_degree in range(9):
            theta_basis = {theta_degree: ParameterPolynomial.coerce(1)}
            terminal_change = x_multiply(theta_basis, derivative)
            value = sum(
                (
                    left_vector[index]
                    * x_coefficient(terminal_change, degree).as_fraction()
                    for index, (row, degree) in enumerate(coordinates)
                    if row == 8
                ),
                F(0),
            )
            tangent_row.append(value)
        Hasse_matrix.append(tangent_row)
    terminal_hasse_rank = rational_rank(Hasse_matrix)
    require(terminal_hasse_rank == 2, "terminal Hasse tangent rank two")
    require(9 - terminal_hasse_rank == 7, "terminal affine dimension seven")

    hasse_substitutions = {
        "Delta": F(896, 15),
        "Theta": F(512, 75),
        "zeta3": F(-3, 2) * Phi,
    }
    hasse_rhs = [
        -value.substitute_many(hasse_substitutions) for value in Hasse_base
    ]
    terminal_theta, terminal_residuals, terminal_rank = (
        solve_rational_matrix_with_parameter_rhs(Hasse_matrix, hasse_rhs)
    )
    require(terminal_rank == 2, "terminal Hasse solve rank")
    require(
        all(not residual for residual in terminal_residuals),
        "three Hasse conditions sufficient",
    )

    A = rows_substitute(A, hasse_substitutions)
    C = rows_substitute(C, hasse_substitutions)
    G = rows_substitute(G, hasse_substitutions)
    terminal_theta_poly = tangent_poly(terminal_theta)
    A[8] = x_add(A[8], x_multiply(terminal_theta_poly, x_derivative(A[0])))
    C[8] = x_add(C[8], x_multiply(terminal_theta_poly, x_derivative(C[0])))
    A_vector = row_coordinate_vector(A, coordinates_2)
    C_vector = row_coordinate_vector(C, coordinates_3)
    require(
        all(not dot_parameter(left, A_vector) for left in left_kernel_2),
        "all P2 projection equations",
    )
    require(
        all(not dot_parameter(left, C_vector) for left in left_kernel_3),
        "all P3 projection equations",
    )
    require(
        x_coefficient(C[8], 9)
        == (1215 * xi10 - 348032) / 495,
        "terminal c89",
    )
    require(
        x_coefficient(C[8], 10)
        == F(27, 40) * (eta - F(3, 2) * Phi),
        "terminal c810",
    )
    check_direct_replay(A, C, G, "symbolic Hasse-compatible")

    final_response = {
        name: value.substitute_many(hasse_substitutions)
        for name, value in response_substitutions.items()
    }
    expected_response = {
        "upsilon5": ParameterPolynomial.coerce(F(-731648, 2025)),
        "U": (475515904 - 109350 * xi10) / 200475,
        "W": -(
            4343625 * Phi * Phi
            - 17172000 * xi10
            + 143826305024
        )
        / 4009500,
        "Z": (
            12506118074368
            - 173745000 * Phi * Phi
            - 195463125 * Phi * eta
            - 926883000 * xi10
        )
        / 108256500,
    }
    for name, expected in expected_response.items():
        require(final_response[name] == expected, f"final response {name}")
    require(
        not any(
            final_response[name].uses(parameter_name)
            for name in ("U", "W", "Z")
            for parameter_name in ("alpha11", "beta11")
        ),
        "top response forgets alpha11 beta11",
    )

    zero_free_parameters = {
        "Phi": 0,
        "eta": 0,
        "xi10": 0,
        "alpha11": 0,
        "beta11": 0,
    }
    A_numeric = specialize_rows(A, zero_free_parameters)
    C_numeric = specialize_rows(C, zero_free_parameters)
    G_numeric = specialize_rows(G, zero_free_parameters)
    check_direct_replay(A_numeric, C_numeric, G_numeric, "numeric interior")
    A_numeric_vector = rational_coordinate_vector(A_numeric, coordinates_2)
    C_numeric_vector = rational_coordinate_vector(C_numeric, coordinates_3)
    require(
        augmented_rank(matrix_2, A_numeric_vector) == rank_2,
        "numeric P2 projected membership",
    )
    require(
        augmented_rank(matrix_3, C_numeric_vector) == rank_3,
        "numeric P3 projected membership",
    )
    require(x_coefficient(C_numeric[8], 9).as_fraction() == F(-348032, 495), "numeric c89")
    require(x_coefficient(C_numeric[8], 10).as_fraction() == 0, "numeric c810")

    U_value = expected_response["U"].substitute_many(zero_free_parameters).as_fraction()
    W_value = expected_response["W"].substitute_many(zero_free_parameters).as_fraction()
    Z_value = expected_response["Z"].substitute_many(zero_free_parameters).as_fraction()
    Lambda_value = U_value + W_value + Z_value
    D_top_value = W_value * W_value - 4 * U_value * Z_value
    require(U_value == F(475515904, 200475), "numeric U")
    require(W_value == F(-35956576256, 1002375), "numeric W")
    require(Z_value == F(3126529518592, 27064125), "numeric Z")
    require(Lambda_value == F(443979321344, 5412825), "numeric Lambda")
    require(
        D_top_value == F(5173344945126466650112, 27128402296875),
        "numeric Dtop",
    )
    require(U_value * Z_value * Lambda_value * D_top_value != 0, "numeric gate interior")

    bracket_hostile_ranks = []
    hostile_channels = {5: 0, 6: 0, 7: 2, 8: 4}
    for row in range(5, 9):
        matrix = tangent_matrix(row)
        hostile = [F(0)] * (row + 1)
        hostile[hostile_channels[row]] = F(1)
        rank = rational_rank(
            [matrix[index] + [hostile[index]] for index in range(row + 1)]
        )
        require(rank == row + 1, f"E{row} one-coordinate hostile")
        bracket_hostile_ranks.append((row, rank))

    hostile_base = {
        "Phi": 0,
        "eta": 0,
        "xi10": 0,
        "alpha11": 0,
        "beta11": 0,
    }
    Hasse_hostiles = (
        (
            "Delta",
            {
                **hostile_base,
                "Delta": F(896, 15) + 1,
                "Theta": F(512, 75),
                "zeta3": 0,
            },
            L_delta,
            F(-1, 3),
        ),
        (
            "zeta3",
            {
                **hostile_base,
                "Delta": F(896, 15),
                "Theta": F(512, 75),
                "zeta3": 1,
            },
            L_zeta,
            F(3, 10),
        ),
    )
    for name, values, affine_value, expected in Hasse_hostiles:
        hostile_A = specialize_rows(pre_hasse_A, values)
        hostile_C = specialize_rows(pre_hasse_C, values)
        hostile_G = specialize_rows(pre_hasse_G, values)
        check_direct_replay(hostile_A, hostile_C, hostile_G, f"{name} Hasse hostile")
        residual = affine_value[0].substitute_many(values).as_fraction()
        require(residual == expected, f"{name} Hasse hostile residual")

    theta_hostile_values = {
        **hostile_base,
        "Delta": F(896, 15),
        "Theta": F(512, 75) + 1,
        "zeta3": 0,
    }
    theta_hostile_A = specialize_rows(pre_hasse_A, theta_hostile_values)
    theta_hostile_C = specialize_rows(pre_hasse_C, theta_hostile_values)
    theta_hostile_G = specialize_rows(pre_hasse_G, theta_hostile_values)
    check_direct_replay(
        theta_hostile_A,
        theta_hostile_C,
        theta_hostile_G,
        "Theta Hasse hostile",
    )
    theta_obstruction = (
        33330 * Delta + 2475 * Theta - 2007808
    ).substitute_many(theta_hostile_values).as_fraction()
    require(theta_obstruction == 2475, "Theta Hasse hostile elimination")

    print("THM4308_SOURCE_NORMAL_BRACKET_HASSE_ROWS8_INDEPENDENT_V1")
    print("GAUGE a=1 gamma=-1/2 I=3/4; K=2848/45-(7/6)Delta")
    print("G4=Delta+Phi*x+(K-1376/45)*x^2+(8/3)*x^4")
    print("G5=upsilon5+eta*x+(4Delta+Theta)*x^2+3Phi*x^3+(2K-1376/45)*x^4")
    print("G6=U+alpha11*x+(5upsilon5+xi10)*x^2+(4eta+zeta3)*x^3+(6Delta+3Theta)*x^4+3Phi*x^5+(K-1376/135)*x^6")
    print("G7=(6U+W)*x^2+(5alpha11+beta11)*x^3+(10upsilon5+4xi10)*x^4+(6eta+3zeta3)*x^5+(4Delta+3Theta)*x^6+Phi*x^7")
    print("G8=(15U+5W+Z)*x^4+(10alpha11+4beta11)*x^5+(10upsilon5+6xi10)*x^6+(4eta+3zeta3)*x^7+(Delta+Theta)*x^8")
    print(f"BRACKET_OPERATOR_ROWS={bracket_rank_rows}")
    print(f"TANGENT_COKERNELS={expected_cokernels}")
    print("E5=2025upsilon5+9000Delta+1350Theta+184832=0; raw=(7/24300)E5")
    print("E6=200475U+109350xi10-5593860Delta-529200Theta-137763328=0; raw=(7/6561000)E6")
    print("E7=801900W+1782000Delta^2+156163200Delta+868725Phi^2+27390480Theta-3434400xi10+12891824128=0; raw=(1/4374000)E7")
    print("E8=21651300Z-225022050Delta^2-59073300DeltaTheta-9512522400Delta+34749000Phi^2+39092625Phi*eta+940522560Theta+185376600xi10-1112446017536=0; raw=(1/275562000)E8")
    print("P2_PROJECTION=63x131 rank=51 left_nullity=12")
    print("P3_PROJECTION=72x204 rank=63 left_nullity=9")
    print("HASSE_IFF Delta=896/15; Theta=512/75; zeta3=-3Phi/2; terminal_rank=2; terminal_affine_dim=7")
    print("TERMINAL c89=(1215xi10-348032)/495; c810=27(eta+zeta3)/40")
    print("RESPONSE K=-32/5; upsilon5=-731648/2025")
    print("RESPONSE U=(475515904-109350xi10)/200475")
    print("RESPONSE W=-(4343625Phi^2-17172000xi10+143826305024)/4009500")
    print("RESPONSE Z=(12506118074368-173745000Phi^2-195463125Phi*eta-926883000xi10)/108256500")
    print(f"NUMERIC_INTERIOR U={U_value}; W={W_value}; Z={Z_value}")
    print(f"NUMERIC_INTERIOR Lambda={Lambda_value}; Dtop={D_top_value}")
    print("NUMERIC_REPLAY J[t^0..t^7],P-G[t^0..t^8],P2,P3=PASS; c89=-348032/495; c810=0")
    print(f"BRACKET_HOSTILES augmented_ranks={bracket_hostile_ranks}; one-coordinate perturbations=REJECT")
    print("HASSE_HOSTILES Delta+1:-1/3; Theta+1:2475; zeta3+1:3/10; bracket/P replay=PASS")
    print(f"CHECKS={CHECKS}")
    print("SCOPE finite weight<=12 fixed-gauge jet and pi_<=8(P2/P3) only; no all-row B2 lift, Keller pair, wall forcing, seam entry, JC2, or DC2")


if __name__ == "__main__":
    main()
