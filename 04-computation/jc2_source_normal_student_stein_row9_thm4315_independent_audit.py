#!/usr/bin/env python3
"""Independent exact audit for THM-4315.

This certificate uses only the Python standard library.  It imports neither
the THM-4308/THM-4315 primary scripts nor a computer-algebra package.  Sparse
polynomials over ``Fraction`` reconstruct the source-normal rows, the
Student--Stein cokernel, the row-nine obstruction, the cubic-corner
specialization, and the row-nine depth projections.

The result is a finite row-nine theorem.  In particular, the Markov language
below is an exact generator/coboundary interpretation of the bracket
operator, not a probabilistic assertion about arbitrary Keller pairs.
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


# ---------------------------------------------------------------------------
# Sparse parameter polynomials and Laurent specializations
# ---------------------------------------------------------------------------

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
        require(scalar_fraction != 0, "parameter polynomial nonzero divisor")
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

    def substitute(self, parameter_name: str, replacement: object) -> "ParameterPolynomial":
        index = PARAMETER_INDEX[parameter_name]
        replacement_poly = self.coerce(replacement)
        result = ParameterPolynomial()
        for monomial, coefficient in self.terms.items():
            exponent = monomial[index]
            base = list(monomial)
            base[index] = 0
            term = ParameterPolynomial({tuple(base): coefficient})
            for _ in range(exponent):
                term = term * replacement_poly
            result = result + term
        return result


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


class Laurent:
    """Sparse Laurent polynomial in one variable phi over Q."""

    def __init__(self, terms: dict[int, F] | None = None) -> None:
        self.terms = {
            int(exponent): F(coefficient)
            for exponent, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def coerce(value: object) -> "Laurent":
        if isinstance(value, Laurent):
            return value
        return Laurent({0: F(value)})

    def __add__(self, other: object) -> "Laurent":
        other_poly = self.coerce(other)
        result = self.terms.copy()
        for exponent, coefficient in other_poly.terms.items():
            result[exponent] = result.get(exponent, F(0)) + coefficient
        return Laurent(result)

    __radd__ = __add__

    def __neg__(self) -> "Laurent":
        return Laurent({exponent: -coefficient for exponent, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Laurent":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "Laurent":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "Laurent":
        other_poly = self.coerce(other)
        result: dict[int, F] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in other_poly.terms.items():
                exponent = left_exponent + right_exponent
                result[exponent] = (
                    result.get(exponent, F(0))
                    + left_coefficient * right_coefficient
                )
        return Laurent(result)

    __rmul__ = __mul__

    def __truediv__(self, scalar: int | F) -> "Laurent":
        scalar_fraction = F(scalar)
        require(scalar_fraction != 0, "Laurent nonzero scalar divisor")
        return Laurent(
            {
                exponent: coefficient / scalar_fraction
                for exponent, coefficient in self.terms.items()
            }
        )

    def __pow__(self, exponent: int) -> "Laurent":
        require(exponent >= 0, "Laurent nonnegative power")
        result = Laurent.coerce(1)
        for _ in range(exponent):
            result = result * self
        return result

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __eq__(self, other: object) -> bool:
        return self.terms == self.coerce(other).terms


PHI_LAURENT = Laurent({1: F(1)})


def evaluate_parameter_polynomial(
    polynomial: ParameterPolynomial, values: dict[str, Laurent]
) -> Laurent:
    result = Laurent()
    for monomial, coefficient in polynomial.terms.items():
        term = Laurent.coerce(coefficient)
        for name, exponent in zip(PARAMETERS, monomial):
            if exponent:
                term = term * (values[name] ** exponent)
        result = result + term
    return result


# ---------------------------------------------------------------------------
# Rational linear algebra
# ---------------------------------------------------------------------------


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


def transpose(matrix: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*matrix)]


def rational_matrix_multiply(
    left: list[list[F]], right: list[list[F]]
) -> list[list[F]]:
    require(bool(left) and bool(right), "nonempty rational matrix product")
    require(len(left[0]) == len(right), "rational matrix product dimensions")
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def rational_nullspace(matrix: list[list[F]]) -> list[list[F]]:
    reduced, pivots = rational_rref(matrix)
    column_count = len(matrix[0]) if matrix else 0
    free_columns = [column for column in range(column_count) if column not in pivots]
    result: list[list[F]] = []
    for free_column in free_columns:
        vector = [F(0)] * column_count
        vector[free_column] = F(1)
        for row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -reduced[row][free_column]
        result.append(vector)
    return result


def primitive_integer_vector(vector: list[F]) -> list[int]:
    denominator = 1
    for entry in vector:
        denominator = lcm(denominator, entry.denominator)
    integers = [int(entry * denominator) for entry in vector]
    divisor = reduce(gcd, (abs(entry) for entry in integers if entry))
    integers = [entry // divisor for entry in integers]
    first_nonzero = next(entry for entry in integers if entry)
    return [-entry for entry in integers] if first_nonzero < 0 else integers


def solve_rational_matrix_with_ring_rhs(
    matrix: list[list[F]], rhs: list[object]
) -> tuple[list[object], list[object], int]:
    row_count = len(matrix)
    column_count = len(matrix[0])
    augmented: list[list[object]] = [
        list(map(F, matrix[row])) + [rhs[row]] for row in range(row_count)
    ]
    pivots: list[tuple[int, int]] = []
    pivot_row = 0
    for column in range(column_count):
        chosen = next(
            (row for row in range(pivot_row, row_count) if augmented[row][column]),
            None,
        )
        if chosen is None:
            continue
        augmented[pivot_row], augmented[chosen] = augmented[chosen], augmented[pivot_row]
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
                augmented[row][-1] - multiplier * augmented[pivot_row][-1]
            )
        pivots.append((pivot_row, column))
        pivot_row += 1
    zero = rhs[0] - rhs[0]
    solution = [zero for _ in range(column_count)]
    for row, column in pivots:
        solution[column] = augmented[row][-1]
    residuals = [augmented[row][-1] for row in range(pivot_row, row_count)]
    return solution, residuals, pivot_row


def affine_combination(
    base: list[object], basis: list[list[F]], coordinates: list[object]
) -> list[object]:
    require(len(basis) == len(coordinates), "affine basis coordinate count")
    result = list(base)
    for vector, coordinate in zip(basis, coordinates):
        require(len(vector) == len(result), "affine basis vector dimension")
        result = [
            value + coefficient * coordinate
            for value, coefficient in zip(result, vector)
        ]
    return result


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
                f"{parameter_name} has constant linear coefficient",
            )
            coefficient += value
        else:
            remainder[monomial] = value
    require(coefficient != 0, f"solve {parameter_name}")
    return -ParameterPolynomial(remainder) / coefficient


# ---------------------------------------------------------------------------
# Sparse x-polynomials and source-normal bracket rows
# ---------------------------------------------------------------------------

XPolynomial = dict[int, object]
Rows = dict[int, XPolynomial]


def x_add(left: XPolynomial, right: XPolynomial) -> XPolynomial:
    result = left.copy()
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, 0) + coefficient
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
                result.get(degree, 0) + left_coefficient * right_coefficient
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


def x_coefficient(poly: XPolynomial, degree: int) -> object:
    return poly.get(degree, 0)


def x_constant(value: object) -> XPolynomial:
    return {0: value} if value else {}


def rows_substitute(
    rows: Rows, parameter_name: str, replacement: ParameterPolynomial
) -> Rows:
    return {
        row: {
            degree: coefficient.substitute(parameter_name, replacement)
            for degree, coefficient in poly.items()
            if coefficient.substitute(parameter_name, replacement)
        }
        for row, poly in rows.items()
    }


def rows_evaluate(rows: Rows, values: dict[str, Laurent]) -> Rows:
    return {
        row: {
            degree: evaluate_parameter_polynomial(coefficient, values)
            for degree, coefficient in poly.items()
            if evaluate_parameter_polynomial(coefficient, values)
        }
        for row, poly in rows.items()
    }


def source_rows() -> tuple[Rows, ParameterPolynomial]:
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
                + coefficient_poly * comb(p_power + y_power, binomial_index)
            )
    return {
        row: {degree: coefficient for degree, coefficient in poly.items() if coefficient}
        for row, poly in result.items()
    }, K


def inherited_rows(K: ParameterPolynomial) -> tuple[Rows, Rows]:
    A: Rows = {
        0: {0: ParameterPolynomial.coerce(1), 2: ParameterPolynomial.coerce(F(1, 4))},
        1: {0: ParameterPolynomial.coerce(F(4, 3)), 2: ParameterPolynomial.coerce(2)},
        2: {0: ParameterPolynomial.coerce(F(-32, 9)), 2: ParameterPolynomial.coerce(F(-4, 5))},
        3: {
            0: ParameterPolynomial.coerce(F(2176, 135)),
            1: -Phi / 2,
            2: F(1088, 315) - F(4, 7) * K,
            4: ParameterPolynomial.coerce(F(-32, 15)),
        },
    }
    C: Rows = {
        0: {1: ParameterPolynomial.coerce(F(-3, 4)), 3: ParameterPolynomial.coerce(F(-1, 8))},
        1: {1: ParameterPolynomial.coerce(-4), 3: ParameterPolynomial.coerce(F(-3, 2))},
        2: {1: ParameterPolynomial.coerce(F(88, 15)), 3: ParameterPolynomial.coerce(F(-12, 5))},
        3: {
            0: F(3, 4) * Phi,
            1: F(-8128, 315) + F(6, 7) * K,
            2: F(3, 8) * Phi,
            3: F(736, 105) + F(3, 7) * K,
            5: ParameterPolynomial.coerce(F(8, 5)),
        },
    }
    return A, C


Q_BOUNDARY: XPolynomial = {
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
    sample = next(iter(A[0].values()))
    boundary = (
        {0: Laurent.coerce(-3), 2: Laurent.coerce(F(-1, 2))}
        if isinstance(sample, Laurent)
        else Q_BOUNDARY
    )
    return x_subtract(
        residual_lower(A, C, row),
        x_scale(x_multiply(boundary, bracket_lower(A, C, row)), F(1, row)),
    )


def bracket_particular(A: Rows, C: Rows, row: int) -> tuple[XPolynomial, XPolynomial]:
    determinant_target = x_scale(bracket_lower(A, C, row), F(-1, row))
    A_constant = F(4, 3) * x_coefficient(determinant_target, 0)
    A_row = x_constant(A_constant)
    remainder = x_subtract(
        determinant_target,
        x_scale(x_multiply({0: 2, 2: 1}, A_row), F(3, 8)),
    )
    require(x_coefficient(remainder, 0) == 0, f"row {row} exact division by x")
    C_row = x_scale(x_shift(remainder, -1), 2)
    return A_row, C_row


def tangent_matrix(next_row: int) -> list[list[F]]:
    """Matrix of theta -> -x theta +(x^2+6)theta'/(2m)."""
    result = [[F(0) for _ in range(next_row)] for _ in range(next_row + 1)]
    for degree in range(next_row):
        if degree:
            result[degree - 1][degree] += F(3 * degree, next_row)
        result[degree + 1][degree] += F(degree, 2 * next_row) - 1
    return result


def tangent_poly(coefficients: list[object]) -> XPolynomial:
    return {
        degree: coefficient
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    }


def dot_ring(vector: list[F], values: list[object]) -> object:
    result = values[0] - values[0]
    for coefficient, value in zip(vector, values):
        result = result + coefficient * value
    return result


# ---------------------------------------------------------------------------
# Student--Stein moment functional
# ---------------------------------------------------------------------------


def student_moments(m: int) -> list[F]:
    """Moments through degree m of density proportional to (x^2+6)^(-m-1)."""
    moments = [F(0)] * (m + 1)
    moments[0] = F(1)
    for degree in range(2, m + 1, 2):
        r = degree // 2
        moments[degree] = (
            F(6 * (2 * r - 1), 2 * m - 2 * r + 1) * moments[degree - 2]
        )
    return moments


def student_primitive_cokernel(m: int) -> list[int]:
    return primitive_integer_vector(student_moments(m))


def student_audit() -> dict[int, list[int]]:
    displayed: dict[int, list[int]] = {}
    for m in range(2, 25):
        matrix = tangent_matrix(m)
        moments = student_moments(m)
        require(rational_rank(matrix) == m, f"Student tangent rank m={m}")
        require(len(rational_nullspace(transpose(matrix))) == 1, f"Student cokernel m={m}")
        require(
            all(
                sum(moments[row] * matrix[row][column] for row in range(m + 1)) == 0
                for column in range(m)
            ),
            f"Student Stein annihilation m={m}",
        )
        # The Pearson diffusion generator is
        # L_m f=(x^2+6)f''-2mx f'.  Its monomial expectation vanishes.
        for degree in range(1, m + 1):
            stationary_generator_mean = (
                -degree * (2 * m - degree + 1) * moments[degree]
                + (6 * degree * (degree - 1) * moments[degree - 2] if degree >= 2 else 0)
            )
            require(stationary_generator_mean == 0, f"Pearson generator mean m={m} k={degree}")
        displayed[m] = student_primitive_cokernel(m)
    require(displayed[5] == [21, 0, 14, 0, 36, 0], "row-five Student weights")
    require(displayed[6] == [77, 0, 42, 0, 84, 0, 360], "row-six Student weights")
    require(displayed[7] == [143, 0, 66, 0, 108, 0, 360, 0], "row-seven Student weights")
    require(displayed[8] == [715, 0, 286, 0, 396, 0, 1080, 0, 5040], "row-eight Student weights")
    require(
        displayed[9] == [12155, 0, 4290, 0, 5148, 0, 11880, 0, 45360, 0],
        "row-nine Student weights",
    )
    return displayed


# ---------------------------------------------------------------------------
# Exact integer-polynomial arithmetic for Q, R, and modular certificates
# ---------------------------------------------------------------------------


def trim(poly: list[int | F]) -> list[int | F]:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def poly_add(left: list[int], right: list[int]) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    return list(map(int, trim(result)))


def poly_scale(poly: list[int], scalar: int) -> list[int]:
    return list(map(int, trim([scalar * coefficient for coefficient in poly])))


def poly_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return list(map(int, trim(result)))


def poly_shift(poly: list[int], amount: int) -> list[int]:
    return [0] * amount + poly


def primitive_integer_polynomial(poly: list[int]) -> tuple[int, list[int]]:
    content = reduce(gcd, (abs(coefficient) for coefficient in poly if coefficient))
    primitive = [coefficient // content for coefficient in poly]
    if primitive[-1] < 0:
        primitive = [-coefficient for coefficient in primitive]
    return content, primitive


def fraction_poly_divmod(left: list[F], right: list[F]) -> tuple[list[F], list[F]]:
    left = list(map(F, trim(left)))
    right = list(map(F, trim(right)))
    require(right != [F(0)], "polynomial nonzero divisor")
    if len(left) < len(right):
        return [F(0)], left
    quotient = [F(0)] * (len(left) - len(right) + 1)
    remainder = left
    while remainder != [F(0)] and len(remainder) >= len(right):
        degree = len(remainder) - len(right)
        coefficient = remainder[-1] / right[-1]
        quotient[degree] = coefficient
        for index, value in enumerate(right):
            remainder[index + degree] -= coefficient * value
        remainder = list(map(F, trim(remainder)))
    return list(map(F, trim(quotient))), remainder


def fraction_poly_gcd(left: list[int], right: list[int]) -> list[F]:
    a = list(map(F, trim(left)))
    b = list(map(F, trim(right)))
    while b != [F(0)]:
        _, remainder = fraction_poly_divmod(a, b)
        a, b = b, remainder
    leading = a[-1]
    return [coefficient / leading for coefficient in a]


def derivative_integer(poly: list[int]) -> list[int]:
    return list(map(int, trim([index * coefficient for index, coefficient in enumerate(poly)][1:])))


def mod_poly(poly: list[int], prime: int) -> list[int]:
    return list(map(int, trim([coefficient % prime for coefficient in poly])))


def mod_poly_add(left: list[int], right: list[int], prime: int) -> list[int]:
    return mod_poly(poly_add(left, right), prime)


def mod_poly_multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    return mod_poly(poly_multiply(left, right), prime)


def fraction_mod(value: int | F, prime: int) -> int:
    rational = F(value)
    denominator = rational.denominator % prime
    require(denominator != 0, "finite-field rational denominator")
    return (rational.numerator % prime) * pow(denominator, -1, prime) % prime


def laurent_mod(polynomial: Laurent, phi_value: int, prime: int) -> int:
    require(phi_value % prime != 0, "finite-field Laurent parameter nonzero")
    result = 0
    for exponent, coefficient in polynomial.terms.items():
        phi_power = (
            pow(phi_value, exponent, prime)
            if exponent >= 0
            else pow(pow(phi_value, -1, prime), -exponent, prime)
        )
        result = (result + fraction_mod(coefficient, prime) * phi_power) % prime
    return result


def scalar_mod(value: object, prime: int, phi_value: int) -> int:
    if isinstance(value, Laurent):
        return laurent_mod(value, phi_value, prime)
    return fraction_mod(F(value), prime)


def modular_rank(
    matrix: list[list[object]], prime: int, phi_value: int
) -> int:
    reduced = [
        [scalar_mod(value, prime, phi_value) for value in row] for row in matrix
    ]
    row_count = len(reduced)
    column_count = len(reduced[0]) if row_count else 0
    pivot_row = 0
    for column in range(column_count):
        chosen = next(
            (row for row in range(pivot_row, row_count) if reduced[row][column]),
            None,
        )
        if chosen is None:
            continue
        reduced[pivot_row], reduced[chosen] = reduced[chosen], reduced[pivot_row]
        inverse = pow(reduced[pivot_row][column], -1, prime)
        reduced[pivot_row] = [value * inverse % prime for value in reduced[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not reduced[row][column]:
                continue
            multiplier = reduced[row][column]
            reduced[row] = [
                (value - multiplier * pivot) % prime
                for value, pivot in zip(reduced[row], reduced[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def polynomial_value_mod(poly: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(poly):
        result = (result * value + coefficient) % prime
    return result


def parameter_polynomial_mod(
    polynomial: ParameterPolynomial, point: dict[str, int], prime: int
) -> int:
    result = 0
    for monomial, coefficient in polynomial.terms.items():
        term = fraction_mod(coefficient, prime)
        for name, exponent in zip(PARAMETERS, monomial):
            term = term * pow(point[name], exponent, prime) % prime
        result = (result + term) % prime
    return result


def algebraic_gate_audit() -> tuple[list[int], list[int]]:
    # alpha=N(X)/(D Phi^3), and alpha^2=-143 A(X)B(X)/C.
    numerator = [
        -74378924775425263164981248,
        -14082869793796263936000,
        6971519208442078125,
    ]
    alpha_denominator = 396452079682031250
    alpha_square_denominator = 6824533500000
    factor_a = [404652032, 14625]
    factor_b = [13056802816, 820125]
    raw_q = poly_add(
        poly_scale(poly_multiply(numerator, numerator), alpha_square_denominator),
        poly_scale(
            poly_shift(poly_multiply(factor_a, factor_b), 3),
            143 * alpha_denominator * alpha_denominator,
        ),
    )
    content, q_poly = primitive_integer_polynomial(raw_q)
    expected_q = [
        44257795605986960276324945338517826145242081533100032,
        16759499408238096044037088463875607198378754048000,
        -6709927871370175861935855782936495259648000000,
        137633556412285429978153325875719168000000000,
        14163685391496771581808513584548828125000,
        316016952601619726458584136962890625,
    ]
    require(content == 853066687500, "Q primitive content")
    require(q_poly == expected_q, "Q exact primitive coefficients")
    require(
        fraction_poly_gcd(q_poly, derivative_integer(q_poly)) == [F(1)],
        "Q squarefree over Q",
    )
    forbidden = poly_multiply([0, 1], poly_multiply(factor_a, factor_b))
    require(
        fraction_poly_gcd(q_poly, forbidden) == [F(1)],
        "Q avoids Phi/U/c2 forbidden factors",
    )

    high_contact = [
        2970579390109346274816679296272171008,
        2284603892441775363795663716352000,
        164114458618573873612800000000,
        7547170421607067494140625,
    ]
    require(
        fraction_poly_gcd(high_contact, forbidden) == [F(1)],
        "high-contact cubic avoids forbidden factors",
    )
    prime = 19
    s_poly = [1, 9, -1]
    t_poly = [3, 9, -1, 4, -8]
    bezout = mod_poly_add(
        mod_poly_multiply(s_poly, q_poly, prime),
        mod_poly_multiply(t_poly, high_contact, prime),
        prime,
    )
    require(bezout == [1], "mod-19 Q/high-contact Bezout identity")
    require(len(mod_poly(q_poly, prime)) == 6, "Q degree retained mod 19")
    require(len(mod_poly(high_contact, prime)) == 4, "R degree retained mod 19")
    return q_poly, high_contact


# ---------------------------------------------------------------------------
# Depth projections through row nine
# ---------------------------------------------------------------------------


def depth_projection_matrix(
    depth: int, max_row: int
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]], list[list[F]]]:
    coordinates = [
        (row, degree)
        for row in range(max_row + 1)
        for degree in range(row + depth + 1)
    ]
    monomials: list[tuple[int, int, int, int]] = []
    for u_power in range(depth + 1):
        for x_power in range(depth - u_power + 1):
            for y_power in range((max_row - u_power) // 2 + 1):
                for p_power in range(max_row - u_power - 2 * y_power + 1):
                    monomials.append((x_power, u_power, p_power, y_power))
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    matrix = [[F(0) for _ in monomials] for _ in coordinates]
    for column, (x_power, u_power, p_power, y_power) in enumerate(monomials):
        for binomial_index in range(p_power + y_power + 1):
            row = u_power + p_power + 2 * y_power + binomial_index
            degree = x_power + 2 * u_power + y_power + 2 * binomial_index
            if row <= max_row:
                matrix[coordinate_index[(row, degree)]][column] += F(
                    comb(p_power + y_power, binomial_index)
                )
    return coordinates, monomials, matrix


def row_coordinate_vector(rows: Rows, coordinates: list[tuple[int, int]]) -> list[object]:
    return [x_coefficient(rows.get(row, {}), degree) for row, degree in coordinates]


def depth_fibre_data(A: Rows, C: Rows, max_row: int) -> dict[str, object]:
    """Return the exact affine depth system on the terminal tangent row."""

    coordinates_2, monomials_2, matrix_2 = depth_projection_matrix(2, max_row)
    coordinates_3, monomials_3, matrix_3 = depth_projection_matrix(3, max_row)
    rank_2 = rational_rank(matrix_2)
    rank_3 = rational_rank(matrix_3)
    left_2 = rational_nullspace(transpose(matrix_2))
    left_3 = rational_nullspace(transpose(matrix_3))

    A_vector = row_coordinate_vector(A, coordinates_2)
    C_vector = row_coordinate_vector(C, coordinates_3)
    A0_derivative = x_derivative(A[0])
    C0_derivative = x_derivative(C[0])
    system: list[list[F]] = []
    base: list[Laurent] = []
    for left_vector, coordinates, vector, derivative in (
        *((left, coordinates_2, A_vector, A0_derivative) for left in left_2),
        *((left, coordinates_3, C_vector, C0_derivative) for left in left_3),
    ):
        base.append(dot_ring(left_vector, vector))
        row_coefficients: list[F] = []
        for theta_degree in range(max_row + 1):
            change = x_multiply({theta_degree: Laurent.coerce(1)}, derivative)
            value = F(0)
            for index, (row, degree) in enumerate(coordinates):
                if row == max_row:
                    coefficient = x_coefficient(change, degree)
                    if coefficient:
                        require(
                            isinstance(coefficient, Laurent)
                            and set(coefficient.terms).issubset({0}),
                            "terminal tangent coefficient rational",
                        )
                        value += left_vector[index] * coefficient.terms.get(0, F(0))
            row_coefficients.append(value)
        system.append(row_coefficients)

    return {
        "max_row": max_row,
        "coordinates_2": coordinates_2,
        "coordinates_3": coordinates_3,
        "monomials_2": monomials_2,
        "monomials_3": monomials_3,
        "projection_2": matrix_2,
        "projection_3": matrix_3,
        "rank_2": rank_2,
        "rank_3": rank_3,
        "left_2": left_2,
        "left_3": left_3,
        "tangent": system,
        "rhs": [-value for value in base],
    }


def depth_residual_system(
    A: Rows,
    C: Rows,
) -> tuple[
    tuple[int, int, int, int],
    tuple[int, int, int, int],
    int,
    int,
    dict[str, object],
]:
    data = depth_fibre_data(A, C, 9)
    coordinates_2 = data["coordinates_2"]
    coordinates_3 = data["coordinates_3"]
    monomials_2 = data["monomials_2"]
    monomials_3 = data["monomials_3"]
    matrix_2 = data["projection_2"]
    matrix_3 = data["projection_3"]
    rank_2 = data["rank_2"]
    rank_3 = data["rank_3"]
    left_2 = data["left_2"]
    left_3 = data["left_3"]
    system = data["tangent"]
    rhs = data["rhs"]
    require(
        (len(coordinates_2), len(monomials_2), rank_2, len(left_2))
        == (75, 160, 59, 16),
        "row-nine P2 projection dimensions",
    )
    require(
        (len(coordinates_3), len(monomials_3), rank_3, len(left_3))
        == (85, 251, 73, 12),
        "row-nine P3 projection dimensions",
    )

    terminal_rank = rational_rank(system)
    solution, residuals, augmented_rank = solve_rational_matrix_with_ring_rhs(
        system, rhs
    )
    require(terminal_rank == 3, "row-nine terminal depth rank three")
    require(augmented_rank == 3, "row-nine augmented depth rank three")
    require(all(not residual for residual in residuals), "row-nine no extra depth consistency")
    require(10 - terminal_rank == 7, "row-nine depth affine dimension seven")

    theta = tangent_poly(solution)
    A[9] = x_add(A[9], x_multiply(theta, x_derivative(A[0])))
    C[9] = x_add(C[9], x_multiply(theta, x_derivative(C[0])))
    final_A = row_coordinate_vector(A, coordinates_2)
    final_C = row_coordinate_vector(C, coordinates_3)
    require(all(not dot_ring(left, final_A) for left in left_2), "all P2 residuals vanish")
    require(all(not dot_ring(left, final_C) for left in left_3), "all P3 residuals vanish")
    return (
        (len(coordinates_2), len(monomials_2), rank_2, len(left_2)),
        (len(coordinates_3), len(monomials_3), rank_3, len(left_3)),
        terminal_rank,
        10 - terminal_rank,
        data,
    )


def finite_field_proposal_sidecar(
    values: dict[str, Laurent],
    q_poly: list[int],
    gate_equations: dict[int, ParameterPolynomial],
    e9_numerator: ParameterPolynomial,
    depth_8: dict[str, object],
    g9_restricted: list[list[F]],
    g9_rhs: list[Laurent],
    depth_9: dict[str, object],
) -> tuple[tuple[int, ...], tuple[int, ...], int, int, int]:
    """Audit one explicitly specified uniform affine proposal over F_19.

    The returned exponents belong to the sequential proposal whose coordinates
    are sampled independently and uniformly at each displayed affine-linear
    gate.  They are not intrinsic probabilities of the source, fibre tower, or
    any putative Keller pair.
    """

    prime = 19
    phi_value = 8
    point = tuple(
        laurent_mod(values[name], phi_value, prime)
        for name in (
            "Delta",
            "Theta",
            "Phi",
            "eta",
            "zeta3",
            "upsilon5",
            "xi10",
            "U",
            "W",
            "Z",
            "alpha11",
            "beta11",
        )
    )
    expected_without_k = (4, 1, 8, 6, 7, 9, 4, 12, 14, 12, 15, 4)
    require(point == expected_without_k, "F19 corner point coordinates")
    k_value = fraction_mod(F(-32, 5), prime)
    require(k_value == 5, "F19 K coordinate")
    displayed_point = point[:2] + (k_value,) + point[2:]
    require(
        displayed_point == (4, 1, 5, 8, 6, 7, 9, 4, 12, 14, 12, 15, 4),
        "F19 displayed positive point",
    )

    x_value = phi_value * phi_value % prime
    require(x_value == 7, "F19 X=Phi^2")
    require(polynomial_value_mod(q_poly, x_value, prime) == 0, "F19 Q survivor")
    factor_u = (820125 * x_value + 13056802816) % prime
    factor_c2 = (14625 * x_value + 404652032) % prime
    require((x_value, factor_u, factor_c2) == (7, 10, 10), "F19 forbidden controls")

    coordinate_values = dict(
        zip(
            (
                "Delta",
                "Theta",
                "Phi",
                "eta",
                "zeta3",
                "upsilon5",
                "xi10",
                "U",
                "W",
                "Z",
                "alpha11",
                "beta11",
            ),
            point,
        )
    )
    for row, equation in gate_equations.items():
        require(
            parameter_polynomial_mod(equation, coordinate_values, prime) == 0,
            f"F19 E{row} gate",
        )
    require(
        parameter_polynomial_mod(e9_numerator, coordinate_values, prime) == 0,
        "F19 full E9 gate",
    )
    delta_value = coordinate_values["Delta"]
    theta_value = coordinate_values["Theta"]
    eta_value = coordinate_values["eta"]
    zeta_value = coordinate_values["zeta3"]
    xi_value = coordinate_values["xi10"]
    u_value = coordinate_values["U"]
    alpha_value = coordinate_values["alpha11"]
    phi_value_mod = coordinate_values["Phi"]
    cleared_corner_equations = (
        15 * delta_value - 896,
        75 * theta_value - 512,
        2 * zeta_value + 3 * phi_value_mod,
        12798000 * xi_value - 4343625 * phi_value_mod**2 - 124805668864,
        107983125 * phi_value_mod * eta_value
        - 2091705253888
        + 258703875 * phi_value_mod**2,
        57591000 * u_value
        + 13 * (820125 * phi_value_mod**2 + 13056802816),
        396452079682031250 * phi_value_mod**3 * alpha_value
        - (
            6971519208442078125 * phi_value_mod**4
            - 14082869793796263936000 * phi_value_mod**2
            - 74378924775425263164981248
        ),
    )
    require(
        all(value % prime == 0 for value in cleared_corner_equations),
        "F19 cleared gate/corner formulas",
    )
    require((coordinate_values["W"] + 2 * coordinate_values["U"]) % prime == 0, "F19 W=-2U")
    require((coordinate_values["Z"] - coordinate_values["U"]) % prime == 0, "F19 Z=U")
    require(
        (coordinate_values["alpha11"] + coordinate_values["beta11"]) % prime == 0,
        "F19 beta=-alpha",
    )
    c2_value = (coordinate_values["upsilon5"] + coordinate_values["xi10"]) % prime
    require(c2_value == 13, "F19 c2 coordinate")
    require(
        (
            coordinate_values["alpha11"] ** 2
            - 4 * coordinate_values["U"] * c2_value
        )
        % prime
        == 0,
        "F19 k1 square relation",
    )
    rho_value = (
        -coordinate_values["alpha11"]
        * pow(2 * coordinate_values["U"] % prime, -1, prime)
    ) % prime
    l1_value = (
        coordinate_values["eta"]
        + coordinate_values["zeta3"]
        + rho_value * coordinate_values["upsilon5"]
    ) % prime
    require(l1_value == 5, "F19 L1 nonzero positive control")

    # Projection and affine-system ranks are preserved at the positive point.
    require(modular_rank(depth_8["projection_2"], prime, phi_value) == 51, "F19 row-eight P2 rank")
    require(modular_rank(depth_8["projection_3"], prime, phi_value) == 63, "F19 row-eight P3 rank")
    require(modular_rank(depth_9["projection_2"], prime, phi_value) == 59, "F19 row-nine P2 rank")
    require(modular_rank(depth_9["projection_3"], prime, phi_value) == 73, "F19 row-nine P3 rank")

    depth_8_rank = modular_rank(depth_8["tangent"], prime, phi_value)
    depth_8_augmented = modular_rank(
        [row + [rhs] for row, rhs in zip(depth_8["tangent"], depth_8["rhs"])],
        prime,
        phi_value,
    )
    g9_rank = modular_rank(g9_restricted, prime, phi_value)
    g9_augmented = modular_rank(
        [row + [rhs] for row, rhs in zip(g9_restricted, g9_rhs)],
        prime,
        phi_value,
    )
    depth_9_rank = modular_rank(depth_9["tangent"], prime, phi_value)
    depth_9_augmented = modular_rank(
        [row + [rhs] for row, rhs in zip(depth_9["tangent"], depth_9["rhs"])],
        prime,
        phi_value,
    )
    require((depth_8_rank, depth_8_augmented) == (2, 2), "F19 row-eight depth consistency")
    require((g9_rank, g9_augmented) == (7, 7), "F19 restricted G9 consistency")
    require((depth_9_rank, depth_9_augmented) == (3, 3), "F19 row-nine depth consistency")

    bracket_increments = tuple(
        modular_rank(tangent_matrix(m), prime, phi_value) for m in range(5, 9)
    )
    increments = bracket_increments + (depth_8_rank, g9_rank, depth_9_rank)
    require(increments == (5, 6, 7, 8, 2, 7, 3), "F19 proposal acceptance increments")
    row_8_exponent = sum(increments[:5])
    row_9_exponent = sum(increments)
    require((row_8_exponent, row_9_exponent) == (28, 38), "F19 proposal survival exponents")
    return displayed_point, increments, row_8_exponent, row_9_exponent, l1_value


# ---------------------------------------------------------------------------
# Full row-nine reconstruction
# ---------------------------------------------------------------------------


def main() -> None:
    student_weights = student_audit()
    G, K = source_rows()
    A, C = inherited_rows(K)

    expected_source_g9: XPolynomial = {
        6: 20 * U + 10 * W + 4 * Z,
        7: 10 * alpha11 + 6 * beta11,
        8: 5 * upsilon5 + 4 * xi10,
        9: eta + zeta3,
    }
    require(G[9] == expected_source_g9, "direct sparse source row G9")

    equations = {
        5: 2025 * upsilon5 + 9000 * Delta + 1350 * Theta + 184832,
        6: 200475 * U + 109350 * xi10 - 5593860 * Delta - 529200 * Theta - 137763328,
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
    new_parameters = {5: "upsilon5", 6: "U", 7: "W", 8: "Z"}

    require(invariant_residual(A, C, 4) == G[4], "inherited G4 bracket compatibility")
    for row in range(4, 8):
        A[row], C[row] = bracket_particular(A, C, row)
        next_row = row + 1
        target = x_subtract(G[next_row], invariant_residual(A, C, next_row))
        matrix = tangent_matrix(next_row)
        cokernel = [F(value) for value in student_weights[next_row]]
        obstruction = dot_ring(
            cokernel,
            [x_coefficient(target, degree) for degree in range(next_row + 1)],
        )
        expected = equations[next_row]
        # Only proportionality is needed here; exact row-nine normalization is
        # separately fixed below.
        shared_monomial = next(iter(expected.terms))
        ratio = obstruction.terms.get(shared_monomial, F(0)) / expected.terms[shared_monomial]
        require(ratio != 0 and obstruction == expected * ratio, f"recover E{next_row}")
        replacement = solve_for_parameter(obstruction, new_parameters[next_row])
        A = rows_substitute(A, new_parameters[next_row], replacement)
        C = rows_substitute(C, new_parameters[next_row], replacement)
        G = rows_substitute(G, new_parameters[next_row], replacement)
        target = x_subtract(G[next_row], invariant_residual(A, C, next_row))
        theta_solution, residuals, rank = solve_rational_matrix_with_ring_rhs(
            matrix,
            [x_coefficient(target, degree) for degree in range(next_row + 1)],
        )
        require(rank == next_row, f"row {next_row} tangent solve rank")
        require(all(not residual for residual in residuals), f"row {next_row} tangent solve")
        theta = tangent_poly(theta_solution)
        A[row] = x_add(A[row], x_multiply(theta, x_derivative(A[0])))
        C[row] = x_add(C[row], x_multiply(theta, x_derivative(C[0])))
        require(invariant_residual(A, C, next_row) == G[next_row], f"match G{next_row}")

    A[8], C[8] = bracket_particular(A, C, 8)
    target_9 = x_subtract(G[9], invariant_residual(A, C, 9))
    row_9_obstruction = dot_ring(
        [F(value) for value in student_weights[9]],
        [x_coefficient(target_9, degree) for degree in range(10)],
    )
    e9_numerator = (
        68294026800 * Delta * Delta
        + 3653910000 * Delta * Theta
        - 5288166000 * Delta * xi10
        + 176911616000 * Delta
        + 1547488800 * Phi * Phi
        + 3987447750 * Phi * alpha11
        + 24602292000 * Phi * eta
        + 4222003500 * Phi * zeta3
        + 2258685000 * Theta * Theta
        + 225494104640 * Theta
        + 1993723875 * eta * eta
        + 263331993600 * xi10
        + 105193437167616
    )
    require(row_9_obstruction == e9_numerator / 328050, "exact E9 normalization")

    xi_corner = F(4343625, 12798000) * PHI_LAURENT**2 + F(124805668864, 12798000)
    eta_corner = Laurent({-1: F(2091705253888, 107983125), 1: F(-2839, 1185)})
    u_corner = Laurent({2: F(-13 * 820125, 57591000), 0: F(-13 * 13056802816, 57591000)})
    alpha_corner = Laurent(
        {
            1: F(6971519208442078125, 396452079682031250),
            -1: F(-14082869793796263936000, 396452079682031250),
            -3: F(-74378924775425263164981248, 396452079682031250),
        }
    )
    values = {
        "Delta": Laurent.coerce(F(896, 15)),
        "Phi": PHI_LAURENT,
        "Theta": Laurent.coerce(F(512, 75)),
        "eta": eta_corner,
        "zeta3": F(-3, 2) * PHI_LAURENT,
        "upsilon5": Laurent.coerce(F(-731648, 2025)),
        "xi10": xi_corner,
        "alpha11": alpha_corner,
        "beta11": -alpha_corner,
        "U": u_corner,
        "W": -2 * u_corner,
        "Z": u_corner,
    }
    require(
        not evaluate_parameter_polynomial(row_9_obstruction, values),
        "corner alpha solves E9",
    )
    alpha_coefficient = evaluate_parameter_polynomial(
        row_9_obstruction.substitute("alpha11", 1)
        - row_9_obstruction.substitute("alpha11", 0),
        values,
    )
    require(
        alpha_coefficient and set(alpha_coefficient.terms) == {1},
        "E9 genuinely solves alpha",
    )

    q_poly, high_contact = algebraic_gate_audit()

    A_corner = rows_evaluate(A, values)
    C_corner = rows_evaluate(C, values)
    G_corner = rows_evaluate(G, values)
    depth_8 = depth_fibre_data(A_corner, C_corner, 8)
    p2_dims_8 = (
        len(depth_8["coordinates_2"]),
        len(depth_8["monomials_2"]),
        depth_8["rank_2"],
        len(depth_8["left_2"]),
    )
    p3_dims_8 = (
        len(depth_8["coordinates_3"]),
        len(depth_8["monomials_3"]),
        depth_8["rank_3"],
        len(depth_8["left_3"]),
    )
    require(p2_dims_8 == (63, 131, 51, 12), "row-eight P2 projection dimensions")
    require(p3_dims_8 == (72, 204, 63, 9), "row-eight P3 projection dimensions")
    depth_8_rank = rational_rank(depth_8["tangent"])
    theta_8_base, depth_8_residuals, depth_8_solve_rank = (
        solve_rational_matrix_with_ring_rhs(depth_8["tangent"], depth_8["rhs"])
    )
    require(depth_8_rank == 2, "row-eight terminal depth rank two")
    require(depth_8_solve_rank == 2, "row-eight depth solve rank two")
    require(all(not value for value in depth_8_residuals), "row-eight depth consistency")
    depth_8_kernel = rational_nullspace(depth_8["tangent"])
    require(len(depth_8_kernel) == 7, "fixed-source J8 affine dimension seven")
    theta_8_base_poly = tangent_poly(theta_8_base)
    A_corner[8] = x_add(
        A_corner[8], x_multiply(theta_8_base_poly, x_derivative(A_corner[0]))
    )
    C_corner[8] = x_add(
        C_corner[8], x_multiply(theta_8_base_poly, x_derivative(C_corner[0]))
    )

    # Restrict the G9 tangent map to the seven-dimensional row-eight depth
    # fibre.  Full rank fixes one and only one J8 point that can underlie J9.
    depth_8_basis_matrix = transpose(depth_8_kernel)
    g9_restricted = rational_matrix_multiply(
        tangent_matrix(9), depth_8_basis_matrix
    )
    require(rational_rank(g9_restricted) == 7, "G9 rank seven on J8")
    target_9_corner = x_subtract(
        G_corner[9], invariant_residual(A_corner, C_corner, 9)
    )
    g9_rhs = [x_coefficient(target_9_corner, degree) for degree in range(10)]
    g9_coordinates, g9_residuals, g9_rank = solve_rational_matrix_with_ring_rhs(
        g9_restricted, g9_rhs
    )
    require(g9_rank == 7, "restricted G9 solve rank")
    require(all(not residual for residual in g9_residuals), "restricted G9 consistency")
    theta_8_kernel = affine_combination(
        [Laurent() for _ in range(9)], depth_8_kernel, g9_coordinates
    )
    theta_8_kernel_poly = tangent_poly(theta_8_kernel)
    A_corner[8] = x_add(
        A_corner[8], x_multiply(theta_8_kernel_poly, x_derivative(A_corner[0]))
    )
    C_corner[8] = x_add(
        C_corner[8], x_multiply(theta_8_kernel_poly, x_derivative(C_corner[0]))
    )
    require(
        invariant_residual(A_corner, C_corner, 9) == G_corner[9],
        "corner row-nine residual matched",
    )
    truncation_image_dimension = len(depth_8_kernel) - rational_rank(g9_restricted)
    require(truncation_image_dimension == 0, "J9 to J8 truncation image dimension zero")
    # The row-eight projected-depth terminal conditions are automatically
    # respected by the unique G9 tangent on this corner.
    expected_c89 = F(1215, 495) * xi_corner - F(348032, 495)
    expected_c810 = F(27, 40) * (eta_corner + F(-3, 2) * PHI_LAURENT)
    require(x_coefficient(C_corner[8], 9) == expected_c89, "corner row-eight c89")
    require(x_coefficient(C_corner[8], 10) == expected_c810, "corner row-eight c810")

    A_corner[9], C_corner[9] = bracket_particular(A_corner, C_corner, 9)
    p2_dims, p3_dims, terminal_rank, terminal_dimension, depth_9 = (
        depth_residual_system(A_corner, C_corner)
    )
    positive_point, proposal_increments, row8_exponent, row9_exponent, l1_value = (
        finite_field_proposal_sidecar(
            values,
            q_poly,
            equations,
            e9_numerator,
            depth_8,
            g9_restricted,
            g9_rhs,
            depth_9,
        )
    )

    print("THM-4315 STUDENT--STEIN ROW-NINE INDEPENDENT AUDIT: PASS")
    print("implementation=stdlib Fraction + sparse parameter/Laurent polynomials; no primary import; no CAS")
    print("Student_law=mu_m(dx) proportional (x^2+6)^(-(m+1)) dx; audited m=2..24")
    print("Stein_operator=A_m(theta)=(x^2+6)theta'-2m*x*theta; D_m=A_m/(2m)")
    print("Pearson_generator=L_m(f)=(x^2+6)f''-2m*x*f'; stationary monomial means vanish")
    print(f"row9_cokernel={student_weights[9]}")
    print("G9=x^6*((20U+10W+4Z)+(10alpha11+6beta11)x+(5upsilon5+4xi10)x^2+(eta+zeta3)x^3)")
    print("E9_normalization=Student_cokernel(target)=E9/328050")
    print("corner_E9_alpha=verified_exact")
    print(f"Q_degree={len(q_poly)-1}; Q_content=853066687500; Q_squarefree=YES; Q_forbidden_gcd=1")
    print(f"high_contact_degree={len(high_contact)-1}; mod19_Bezout=s*Q+t*R=1")
    print(f"P2_row9=ambient:{p2_dims[0]} columns:{p2_dims[1]} rank:{p2_dims[2]} left_nullity:{p2_dims[3]}")
    print(f"P3_row9=ambient:{p3_dims[0]} columns:{p3_dims[1]} rank:{p3_dims[2]} left_nullity:{p3_dims[3]}")
    print(f"terminal_depth_rank={terminal_rank}; augmented_rank={terminal_rank}; affine_dimension={terminal_dimension}")
    print("fixed_source_fibres=J8:A7 J9:A7; G9_rank_on_J8=7; truncation_J9_to_J8_image_dimension=0")
    print("stochastic_scope=sequential independent uniform affine-coordinate proposal at the fixed compatible F19 source")
    print(f"F19_proposal_increments={proposal_increments}; survival_through_row8=19^-{row8_exponent}; survival_through_row9=19^-{row9_exponent}")
    print(f"F19_positive_point={positive_point}; L1={l1_value}; projection_and_affine_ranks_preserved=YES")
    print("full_fibre_uniform_marginals_projectively_consistent=NO; pushforward_from_J9_is_delta_on_J8")
    print("proposal_caveat=these exponents are proposal-dependent acceptance fractions, not intrinsic source probabilities or an extinction theorem")
    print("consequence=L1=0 high-contact cubic is incompatible with row-nine bracket lift")
    print("scope=finite row-nine bracket/depth gate only; no all-row lift, termination, seam entry, JC2, or DC2")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
