#!/usr/bin/env python3
"""Independent exact audit of the source-normal row-ten obstruction.

This certificate uses only the Python standard library.  It imports neither a
THM-4308/THM-4315 computation nor a computer-algebra package.  Sparse Laurent
polynomials in ``Phi`` reconstruct the source rows directly from

    p=t(1+x^2 t),       y=x t p,

transport the inherited bracket data through G9, and evaluate the row-ten
Student--Stein cokernel.  Integer-polynomial arithmetic then verifies the
claimed cubic normalization and explicit modular Bezout certificates.

The consequence is deliberately finite and local: the audited cubic corner
cannot extend through the row-ten bracket equation.  No row-ten depth
projection, all-row lifting statement, seam entry, or Jacobian-conjecture
conclusion is asserted.
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
# Sparse Laurent and x-polynomial arithmetic
# ---------------------------------------------------------------------------


class Laurent:
    """Sparse Laurent polynomial in Phi over Q."""

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
        return Laurent(
            {exponent: -coefficient for exponent, coefficient in self.terms.items()}
        )

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
        divisor = F(scalar)
        require(divisor != 0, "Laurent nonzero rational divisor")
        return Laurent(
            {exponent: coefficient / divisor for exponent, coefficient in self.terms.items()}
        )

    def __pow__(self, exponent: int) -> "Laurent":
        require(exponent >= 0, "Laurent nonnegative power")
        result = Laurent.coerce(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __eq__(self, other: object) -> bool:
        return self.terms == self.coerce(other).terms


PHI = Laurent({1: F(1)})
XPolynomial = dict[int, Laurent]
Rows = dict[int, XPolynomial]


def x_add(left: XPolynomial, right: XPolynomial) -> XPolynomial:
    result = left.copy()
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, Laurent()) + coefficient
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
                result.get(degree, Laurent())
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


def x_coefficient(poly: XPolynomial, degree: int) -> Laurent:
    return poly.get(degree, Laurent())


def x_constant(value: object) -> XPolynomial:
    coefficient = Laurent.coerce(value)
    return {0: coefficient} if coefficient else {}


# ---------------------------------------------------------------------------
# Rational linear algebra and Student--Stein cokernels
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


def tangent_matrix(m: int) -> list[list[F]]:
    """Matrix of D_m theta=-x theta +(x^2+6)theta'/(2m)."""
    result = [[F(0) for _ in range(m)] for _ in range(m + 1)]
    for degree in range(m):
        if degree:
            result[degree - 1][degree] += F(3 * degree, m)
        result[degree + 1][degree] += F(degree, 2 * m) - 1
    return result


def student_moments(m: int) -> list[F]:
    """Moments through m for density proportional to (x^2+6)^(-m-1)."""
    moments = [F(0)] * (m + 1)
    moments[0] = F(1)
    for degree in range(2, m + 1, 2):
        r = degree // 2
        moments[degree] = (
            F(6 * (2 * r - 1), 2 * m - 2 * r + 1) * moments[degree - 2]
        )
    return moments


def primitive_integer_vector(vector: list[F]) -> list[int]:
    denominator = 1
    for entry in vector:
        denominator = lcm(denominator, entry.denominator)
    integers = [int(entry * denominator) for entry in vector]
    divisor = reduce(gcd, (abs(entry) for entry in integers if entry))
    integers = [entry // divisor for entry in integers]
    first = next(entry for entry in integers if entry)
    return [-entry for entry in integers] if first < 0 else integers


def student_audit() -> dict[int, list[int]]:
    displayed: dict[int, list[int]] = {}
    for m in range(2, 21):
        matrix = tangent_matrix(m)
        moments = student_moments(m)
        require(rational_rank(matrix) == m, f"tangent rank m={m}")
        require(
            all(
                sum(moments[row] * matrix[row][column] for row in range(m + 1))
                == 0
                for column in range(m)
            ),
            f"Student--Stein annihilation m={m}",
        )
        for degree in range(1, m + 1):
            # L_m f=(x^2+6)f''-2m*x*f'.
            generator_mean = (
                -degree * (2 * m - degree + 1) * moments[degree]
                + (
                    6 * degree * (degree - 1) * moments[degree - 2]
                    if degree >= 2
                    else 0
                )
            )
            require(generator_mean == 0, f"Pearson generator mean m={m} k={degree}")
        displayed[m] = primitive_integer_vector(moments)
    require(
        displayed[9] == [12155, 0, 4290, 0, 5148, 0, 11880, 0, 45360, 0],
        "row-nine Student weights",
    )
    require(
        displayed[10]
        == [46189, 0, 14586, 0, 15444, 0, 30888, 0, 99792, 0, 489888],
        "row-ten Student weights",
    )
    return displayed


def solve_rational_matrix_with_laurent_rhs(
    matrix: list[list[F]], rhs: list[Laurent]
) -> tuple[list[Laurent], list[Laurent], int]:
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
        augmented[pivot_row][-1] = Laurent.coerce(augmented[pivot_row][-1]) / pivot
        for row in range(row_count):
            if row == pivot_row or not augmented[row][column]:
                continue
            multiplier = F(augmented[row][column])
            for index in range(column_count):
                augmented[row][index] = (
                    F(augmented[row][index])
                    - multiplier * F(augmented[pivot_row][index])
                )
            augmented[row][-1] = Laurent.coerce(augmented[row][-1]) - (
                Laurent.coerce(augmented[pivot_row][-1]) * multiplier
            )
        pivots.append((pivot_row, column))
        pivot_row += 1
    solution = [Laurent() for _ in range(column_count)]
    for row, column in pivots:
        solution[column] = Laurent.coerce(augmented[row][-1])
    residuals = [
        Laurent.coerce(augmented[row][-1]) for row in range(pivot_row, row_count)
    ]
    return solution, residuals, pivot_row


def tangent_poly(coefficients: list[Laurent]) -> XPolynomial:
    return {
        degree: coefficient
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    }


def dot_laurent(vector: list[int], values: list[Laurent]) -> Laurent:
    result = Laurent()
    for coefficient, value in zip(vector, values):
        result = result + coefficient * value
    return result


# ---------------------------------------------------------------------------
# Source-normal rows and bracket transport
# ---------------------------------------------------------------------------


def corner_values() -> dict[str, Laurent]:
    xi = F(4343625, 12798000) * PHI**2 + F(124805668864, 12798000)
    eta = Laurent(
        {-1: F(2091705253888, 107983125), 1: F(-2839, 1185)}
    )
    u_value = Laurent(
        {
            2: F(-13 * 820125, 57591000),
            0: F(-13 * 13056802816, 57591000),
        }
    )
    alpha = Laurent(
        {
            1: F(6971519208442078125, 396452079682031250),
            -1: F(-14082869793796263936000, 396452079682031250),
            -3: F(-74378924775425263164981248, 396452079682031250),
        }
    )
    return {
        "Delta": Laurent.coerce(F(896, 15)),
        "Phi": PHI,
        "Theta": Laurent.coerce(F(512, 75)),
        "eta": eta,
        "zeta3": F(-3, 2) * PHI,
        "upsilon5": Laurent.coerce(F(-731648, 2025)),
        "xi10": xi,
        "alpha11": alpha,
        "beta11": -alpha,
        "U": u_value,
        "W": -2 * u_value,
        "Z": u_value,
    }


def source_rows(values: dict[str, Laurent]) -> Rows:
    delta = values["Delta"]
    kappa = F(2848, 45) - F(7, 6) * delta
    require(kappa == F(-32, 5), "corner K=-32/5")
    source_terms = (
        (Laurent.coerce(-3), 1, 0),
        (Laurent.coerce(F(8, 3)), 2, 0),
        (Laurent.coerce(F(-1376, 135)), 3, 0),
        (kappa, 0, 2),
        (values["Phi"], 2, 1),
        (delta, 4, 0),
        (values["Theta"], 1, 2),
        (values["eta"], 3, 1),
        (values["zeta3"], 0, 3),
        (values["upsilon5"], 5, 0),
        (values["xi10"], 2, 2),
        (values["alpha11"], 4, 1),
        (values["beta11"], 1, 3),
        (values["U"], 6, 0),
        (values["W"], 3, 2),
        (values["Z"], 0, 4),
    )
    rows: Rows = {1: {2: Laurent.coerce(F(-1, 2))}}
    # Since p=t(1+x^2t) and y=xtp, a monomial p^a y^b is
    # x^b t^(a+2b) (1+x^2t)^(a+b).
    for coefficient, p_power, y_power in source_terms:
        for binomial_index in range(p_power + y_power + 1):
            row = p_power + 2 * y_power + binomial_index
            degree = y_power + 2 * binomial_index
            rows.setdefault(row, {})[degree] = (
                rows.setdefault(row, {}).get(degree, Laurent())
                + coefficient * comb(p_power + y_power, binomial_index)
            )
    return {
        row: {degree: coefficient for degree, coefficient in poly.items() if coefficient}
        for row, poly in rows.items()
    }


def inherited_rows(values: dict[str, Laurent]) -> tuple[Rows, Rows]:
    kappa = F(2848, 45) - F(7, 6) * values["Delta"]
    a_rows: Rows = {
        0: {0: Laurent.coerce(1), 2: Laurent.coerce(F(1, 4))},
        1: {0: Laurent.coerce(F(4, 3)), 2: Laurent.coerce(2)},
        2: {0: Laurent.coerce(F(-32, 9)), 2: Laurent.coerce(F(-4, 5))},
        3: {
            0: Laurent.coerce(F(2176, 135)),
            1: -values["Phi"] / 2,
            2: Laurent.coerce(F(1088, 315) - F(4, 7) * kappa),
            4: Laurent.coerce(F(-32, 15)),
        },
    }
    c_rows: Rows = {
        0: {1: Laurent.coerce(F(-3, 4)), 3: Laurent.coerce(F(-1, 8))},
        1: {1: Laurent.coerce(-4), 3: Laurent.coerce(F(-3, 2))},
        2: {1: Laurent.coerce(F(88, 15)), 3: Laurent.coerce(F(-12, 5))},
        3: {
            0: F(3, 4) * values["Phi"],
            1: Laurent.coerce(F(-8128, 315) + F(6, 7) * kappa),
            2: F(3, 8) * values["Phi"],
            3: Laurent.coerce(F(736, 105) + F(3, 7) * kappa),
            5: Laurent.coerce(F(8, 5)),
        },
    }
    return a_rows, c_rows


BOUNDARY: XPolynomial = {
    0: Laurent.coerce(-3),
    2: Laurent.coerce(F(-1, 2)),
}


def bracket_lower(a_rows: Rows, c_rows: Rows, row: int) -> XPolynomial:
    result: XPolynomial = {}
    for index in range(1, row):
        result = x_add(
            result,
            x_add(
                x_scale(
                    x_multiply(x_derivative(a_rows[index]), c_rows[row - index]),
                    row - index,
                ),
                x_scale(
                    x_multiply(a_rows[index], x_derivative(c_rows[row - index])),
                    -index,
                ),
            ),
        )
    return result


def residual_lower(a_rows: Rows, c_rows: Rows, row: int) -> XPolynomial:
    result: XPolynomial = {}
    for index in range(1, row):
        result = x_add(result, x_multiply(c_rows[index], c_rows[row - index]))
    for first in range(row):
        for second in range(row):
            third = row - first - second
            if 0 <= third < row:
                result = x_subtract(
                    result,
                    x_multiply(
                        x_multiply(a_rows[first], a_rows[second]), a_rows[third]
                    ),
                )
    return result


def invariant_residual(a_rows: Rows, c_rows: Rows, row: int) -> XPolynomial:
    return x_subtract(
        residual_lower(a_rows, c_rows, row),
        x_scale(x_multiply(BOUNDARY, bracket_lower(a_rows, c_rows, row)), F(1, row)),
    )


def bracket_particular(
    a_rows: Rows, c_rows: Rows, row: int
) -> tuple[XPolynomial, XPolynomial]:
    determinant_target = x_scale(bracket_lower(a_rows, c_rows, row), F(-1, row))
    a_constant = F(4, 3) * x_coefficient(determinant_target, 0)
    a_row = x_constant(a_constant)
    remainder = x_subtract(
        determinant_target,
        x_scale(
            x_multiply(
                {0: Laurent.coerce(2), 2: Laurent.coerce(1)}, a_row
            ),
            F(3, 8),
        ),
    )
    require(x_coefficient(remainder, 0) == 0, f"row {row} exact division by x")
    c_row = x_scale(x_shift(remainder, -1), 2)
    return a_row, c_row


def transport_to_next_source_row(
    a_rows: Rows,
    c_rows: Rows,
    source: Rows,
    row: int,
    student_weights: list[int],
) -> Laurent:
    """Choose row ``row`` and match source row ``row+1`` when compatible."""
    a_rows[row], c_rows[row] = bracket_particular(a_rows, c_rows, row)
    next_row = row + 1
    target = x_subtract(source[next_row], invariant_residual(a_rows, c_rows, next_row))
    obstruction = dot_laurent(
        student_weights,
        [x_coefficient(target, degree) for degree in range(next_row + 1)],
    )
    solution, residuals, rank = solve_rational_matrix_with_laurent_rhs(
        tangent_matrix(next_row),
        [x_coefficient(target, degree) for degree in range(next_row + 1)],
    )
    require(rank == next_row, f"G{next_row} tangent rank")
    if not obstruction:
        require(all(not residual for residual in residuals), f"G{next_row} compatibility")
        theta = tangent_poly(solution)
        a_rows[row] = x_add(
            a_rows[row], x_multiply(theta, x_derivative(a_rows[0]))
        )
        c_rows[row] = x_add(
            c_rows[row], x_multiply(theta, x_derivative(c_rows[0]))
        )
        require(
            invariant_residual(a_rows, c_rows, next_row) == source[next_row],
            f"exact source match G{next_row}",
        )
    else:
        require(any(residuals), f"G{next_row} incompatible obstruction detected")
    return obstruction


# ---------------------------------------------------------------------------
# Integer-polynomial normalization and modular Bezout certificates
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
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] += left_coefficient * right_coefficient
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
    dividend = list(map(F, trim(left)))
    divisor = list(map(F, trim(right)))
    require(divisor != [F(0)], "polynomial nonzero divisor")
    if len(dividend) < len(divisor):
        return [F(0)], dividend
    quotient = [F(0)] * (len(dividend) - len(divisor) + 1)
    remainder = dividend
    while remainder != [F(0)] and len(remainder) >= len(divisor):
        degree = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[degree] = coefficient
        for index, value in enumerate(divisor):
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
    return list(
        map(int, trim([index * coefficient for index, coefficient in enumerate(poly)][1:]))
    )


def mod_poly(poly: list[int], prime: int) -> list[int]:
    return list(map(int, trim([coefficient % prime for coefficient in poly])))


def mod_poly_add(left: list[int], right: list[int], prime: int) -> list[int]:
    return mod_poly(poly_add(left, right), prime)


def mod_poly_multiply(left: list[int], right: list[int], prime: int) -> list[int]:
    return mod_poly(poly_multiply(left, right), prime)


P10 = [
    -2955996966894649715961849999793382752256,
    -49516799750570385919467992383488000,
    -114199708075156991490870528000000,
    68842335386673891964107421875,
]


def algebraic_audit() -> tuple[list[int], list[int], list[int]]:
    # Rebuild Q independently by clearing alpha^2=-143*A(X)*B(X)/C,
    # where alpha=N(X)/(D*Phi^3) and X=Phi^2.
    alpha_numerator = [
        -74378924775425263164981248,
        -14082869793796263936000,
        6971519208442078125,
    ]
    alpha_denominator = 396452079682031250
    alpha_square_denominator = 6824533500000
    factor_a = [404652032, 14625]
    factor_b = [13056802816, 820125]
    raw_q = poly_add(
        poly_scale(
            poly_multiply(alpha_numerator, alpha_numerator),
            alpha_square_denominator,
        ),
        poly_scale(
            poly_shift(poly_multiply(factor_a, factor_b), 3),
            143 * alpha_denominator * alpha_denominator,
        ),
    )
    content_q, q_poly = primitive_integer_polynomial(raw_q)
    expected_q = [
        44257795605986960276324945338517826145242081533100032,
        16759499408238096044037088463875607198378754048000,
        -6709927871370175861935855782936495259648000000,
        137633556412285429978153325875719168000000000,
        14163685391496771581808513584548828125000,
        316016952601619726458584136962890625,
    ]
    require(content_q == 853066687500, "Q clearing content")
    require(q_poly == expected_q, "Q exact primitive coefficients")

    content_p, primitive_p = primitive_integer_polynomial(P10)
    require(content_p == 1 and primitive_p == P10, "P10 primitive normalization")
    require(
        fraction_poly_gcd(P10, derivative_integer(P10)) == [F(1)],
        "P10 squarefree over Q",
    )
    require(fraction_poly_gcd(P10, q_poly) == [F(1)], "P10/Q gcd over Q")

    # Every nonconstant forbidden denominator is represented here:
    # X=Phi^2, A(X)=14625X+404652032 (c2), and
    # B(X)=820125X+13056802816 (U).  Rational integer denominators are nonzero.
    forbidden = poly_multiply([0, 1], poly_multiply(factor_a, factor_b))
    require(
        fraction_poly_gcd(P10, forbidden) == [F(1)],
        "P10 avoids Phi/c2/U forbidden factors",
    )

    # Explicit simultaneous certificates over F_23.  Coefficients are stored
    # in ascending order and may be represented by centered residues.
    prime = 23
    s_p = [2, 10, -2, -10, -1]
    t_q = [2, -10, 5]
    u_p = [7, -7, -4]
    v_f = [-10, -3, -11]
    require(
        mod_poly_add(
            mod_poly_multiply(s_p, P10, prime),
            mod_poly_multiply(t_q, q_poly, prime),
            prime,
        )
        == [1],
        "mod-23 P10/Q Bezout certificate",
    )
    require(
        mod_poly_add(
            mod_poly_multiply(u_p, P10, prime),
            mod_poly_multiply(v_f, forbidden, prime),
            prime,
        )
        == [1],
        "mod-23 P10/forbidden Bezout certificate",
    )
    require(len(mod_poly(P10, prime)) == 4, "P10 degree retained mod 23")
    require(len(mod_poly(q_poly, prime)) == 6, "Q degree retained mod 23")
    require(len(mod_poly(forbidden, prime)) == 4, "forbidden degree retained mod 23")
    return q_poly, forbidden, primitive_p


# ---------------------------------------------------------------------------
# Full source-to-row-ten audit
# ---------------------------------------------------------------------------


def main() -> None:
    student_weights = student_audit()
    values = corner_values()
    source = source_rows(values)
    a_rows, c_rows = inherited_rows(values)

    # These identities are read directly from the sparse substitution, not
    # entered as source rows.
    require(
        source[9]
        == {
            6: 4 * values["U"],
            7: 4 * values["alpha11"],
            8: 5 * values["upsilon5"] + 4 * values["xi10"],
            9: values["eta"] + values["zeta3"],
        },
        "direct sparse source row G9",
    )
    require(
        source[10]
        == {
            8: values["U"],
            9: values["alpha11"],
            10: values["upsilon5"] + values["xi10"],
        },
        "direct sparse source row G10",
    )
    require(
        invariant_residual(a_rows, c_rows, 4) == source[4],
        "inherited G4 compatibility",
    )

    # Gate/corner values make the E5,...,E9 Student obstructions vanish.
    for row in range(4, 9):
        obstruction = transport_to_next_source_row(
            a_rows, c_rows, source, row, student_weights[row + 1]
        )
        require(not obstruction, f"corner E{row + 1}=0")

    expected_c89 = F(1215, 495) * values["xi10"] - F(348032, 495)
    expected_c810 = F(27, 40) * (values["eta"] + values["zeta3"])
    require(x_coefficient(c_rows[8], 9) == expected_c89, "transported row-eight c89")
    require(x_coefficient(c_rows[8], 10) == expected_c810, "transported row-eight c810")

    # The next transport is incompatible.  Its cokernel is invariant under
    # every possible theta_9, so no choice of the new row-nine tangent can
    # alter the obstruction.
    obstruction_10 = transport_to_next_source_row(
        a_rows, c_rows, source, 9, student_weights[10]
    )
    normalization_denominator = 2518243204518514160156250
    expected_obstruction = Laurent(
        {
            2 * index - 4: F(143 * coefficient, normalization_denominator)
            for index, coefficient in enumerate(P10)
        }
    )
    require(
        obstruction_10 == expected_obstruction,
        "symbolic E10=143*P10(Phi^2)/(D*Phi^4)",
    )
    require(obstruction_10, "row-ten obstruction nonzero Laurent polynomial")

    q_poly, forbidden, p10 = algebraic_audit()

    print("JC2 SOURCE-NORMAL STUDENT--STEIN ROW-TEN INDEPENDENT AUDIT: PASS")
    print("implementation=stdlib Fraction + sparse Laurent/x polynomials; no primary import; no CAS")
    print("source=p=t(1+x^2*t), y=x*t*p; G4..G10 reconstructed directly")
    print("transport=bracket particular + exact Student tangent solve through G9")
    print("Student_law=mu_m(dx) proportional (x^2+6)^(-(m+1)) dx; audited m=2..20")
    print("Stein_operator=A_m(theta)=(x^2+6)theta'-2m*x*theta; D_m=A_m/(2m)")
    print("Pearson_generator=L_m(f)=(x^2+6)f''-2m*x*f'; stationary monomial means vanish")
    print(f"row10_cokernel={student_weights[10]}")
    print("G10=U*x^8+alpha11*x^9+(upsilon5+xi10)*x^10 on the audited corner")
    print(
        "E10=143*P10(Phi^2)/(2518243204518514160156250*Phi^4) verified as a Laurent identity"
    )
    print(f"P10_ascending={p10}")
    print(f"P10_degree={len(p10)-1}; P10_content=1; P10_squarefree=YES")
    print(f"Q_degree={len(q_poly)-1}; gcd(P10,Q)=1")
    print(f"forbidden_degree={len(forbidden)-1}; forbidden=X*(14625X+404652032)*(820125X+13056802816)")
    print("mod23_P10_Q_Bezout=(-X^4-10X^3-2X^2+10X+2)*P10+(5X^2-10X+2)*Q=1")
    print("mod23_P10_F_Bezout=(-4X^2-7X+7)*P10+(-11X^2-3X-10)*F=1")
    print("consequence=the audited exact cubic corner has no row-ten bracket lift")
    print("scope=finite row-ten bracket necessity only; depth-10, all-row lift, seam entry, JC2, and DC2 remain unclaimed")
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
