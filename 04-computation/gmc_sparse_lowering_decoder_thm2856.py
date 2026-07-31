#!/usr/bin/env python3
"""Exact companion for THM-2856.

All truth-bearing arithmetic is over integers or fractions.  The script has
no optional dependencies and uses explicit exception gates, so ``python``
and ``python -O`` execute the same checks.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb, factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add_poly(left, right):
    size = max(len(left), len(right))
    out = [Fraction(0) for _ in range(size)]
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return trim(out)


def scale_poly(poly, scalar):
    return trim([scalar * value for value in poly])


def mul_poly(left, right):
    out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return trim(out)


def pow_poly(poly, exponent):
    out = [Fraction(1)]
    base = list(poly)
    power = exponent
    while power:
        if power & 1:
            out = mul_poly(out, base)
        base = mul_poly(base, base)
        power //= 2
    return out


def derivative(poly):
    if len(poly) <= 1:
        return [Fraction(0)]
    return trim([index * poly[index] for index in range(1, len(poly))])


def factorial_readout(poly):
    return sum(value * factorial(index) for index, value in enumerate(poly))


def rising(value, degree):
    out = 1
    for offset in range(degree):
        out *= value + offset
    return out


def determinant_bareiss(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "determinant matrix not square")
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap_index = None
            for row_index in range(pivot_index + 1, size):
                if work[row_index][pivot_index] != 0:
                    swap_index = row_index
                    break
            require(swap_index is not None, "unexpected singular Bareiss pivot")
            work[pivot_index], work[swap_index] = (
                work[swap_index],
                work[pivot_index],
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row_index in range(pivot_index + 1, size):
            for column_index in range(pivot_index + 1, size):
                numerator = (
                    work[row_index][column_index] * pivot
                    - work[row_index][pivot_index]
                    * work[pivot_index][column_index]
                )
                require(
                    numerator % previous == 0,
                    "Bareiss exact-division failure",
                )
                work[row_index][column_index] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def rank_fraction(matrix):
    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    require(
        all(len(row) == column_count for row in work),
        "rank matrix is ragged",
    )
    pivot_row = 0
    for column in range(column_count):
        found = None
        for row in range(pivot_row, row_count):
            if work[row][column] != 0:
                found = row
                break
        if found is None:
            continue
        work[pivot_row], work[found] = work[found], work[pivot_row]
        pivot = work[pivot_row][column]
        work[pivot_row] = [entry / pivot for entry in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][j] - factor * work[pivot_row][j]
                    for j in range(column_count)
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def invert_fraction(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "inverse matrix not square")
    work = [
        [Fraction(value) for value in row]
        + [Fraction(int(row_index == column_index)) for column_index in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot_row = None
        for row in range(column, size):
            if work[row][column] != 0:
                pivot_row = row
                break
        require(pivot_row is not None, "matrix unexpectedly singular")
        work[column], work[pivot_row] = work[pivot_row], work[column]
        pivot = work[column][column]
        work[column] = [entry / pivot for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            if factor:
                work[row] = [
                    work[row][j] - factor * work[column][j]
                    for j in range(2 * size)
                ]
    return [row[size:] for row in work]


def monic_laguerre(degree):
    return [
        Fraction(
            (-1) ** (degree + power)
            * factorial(degree)
            * comb(degree, power),
            factorial(power),
        )
        for power in range(degree + 1)
    ]


def divmod_poly(dividend, divisor):
    numerator = trim(dividend)
    denominator = trim(divisor)
    require(denominator != [0], "polynomial division by zero")
    if len(numerator) < len(denominator):
        return [Fraction(0)], numerator
    quotient = [Fraction(0) for _ in range(len(numerator) - len(denominator) + 1)]
    lead = denominator[-1]
    while numerator != [0] and len(numerator) >= len(denominator):
        shift = len(numerator) - len(denominator)
        scalar = numerator[-1] / lead
        quotient[shift] += scalar
        for index, value in enumerate(denominator):
            numerator[shift + index] -= scalar * value
        numerator = trim(numerator)
    return trim(quotient), trim(numerator)


def gcd_poly(left, right):
    a = trim(left)
    b = trim(right)
    while b != [0]:
        _, remainder = divmod_poly(a, b)
        a, b = b, remainder
    require(a != [0], "zero polynomial gcd")
    return scale_poly(a, Fraction(1, 1) / a[-1])


def reduce_monic(poly, modulus):
    quotient, remainder = divmod_poly(poly, modulus)
    del quotient
    dimension = len(modulus) - 1
    return remainder + [Fraction(0) for _ in range(dimension - len(remainder))]


def multiply_by_u(vector, modulus):
    return reduce_monic([Fraction(0)] + list(vector), modulus)


def representative_derivative(vector, dimension):
    raw = derivative(vector)
    return raw + [Fraction(0) for _ in range(dimension - len(raw))]


def vandermonde(exponents):
    out = 1
    for left_index in range(len(exponents)):
        for right_index in range(left_index + 1, len(exponents)):
            out *= exponents[right_index] - exponents[left_index]
    return out


def sparse_polynomial(exponents, coefficients):
    degree = max(exponents)
    out = [Fraction(0) for _ in range(degree + 1)]
    for exponent, coefficient in zip(exponents, coefficients):
        out[exponent] += Fraction(coefficient, factorial(exponent))
    return trim(out)


print("THM-2856 SPARSE LOWERING / LAGUERRE DEFECT AUDIT")

support_digest = sha256()
support_cells = 0
four_slot_cells = 0
for slot_count in range(1, 7):
    for support in combinations(range(13), slot_count):
        matrix = [
            [rising(exponent, degree) for exponent in support]
            for degree in range(slot_count)
        ]
        determinant = determinant_bareiss(matrix)
        expected = vandermonde(support)
        require(
            determinant == expected,
            f"rising-factorial determinant mismatch on {support}",
        )
        support_digest.update(f"{support}:{determinant}\n".encode("ascii"))
        support_cells += 1

for support in combinations(range(16), 4):
    derivative_restriction = [
        [
            rising(support[column], degree) - rising(support[-1], degree)
            for column in range(3)
        ]
        for degree in range(1, 4)
    ]
    require(
        determinant_bareiss(derivative_restriction) != 0,
        f"mean-zero lowering bank failed on {support}",
    )
    four_slot_cells += 1

print(f"rising_vandermonde_cells={support_cells}")
print(f"rising_vandermonde_digest={support_digest.hexdigest()}")
print(f"four_slot_mean_zero_rank3_cells={four_slot_cells}")

sample_support = (0, 3, 8, 15)
sample_coefficients = (Fraction(5, 7), Fraction(-11, 13), Fraction(17, 19), Fraction(-23, 29))
sample = sparse_polynomial(sample_support, sample_coefficients)
sample_derivative = derivative(sample)
sample_observations = [factorial_readout(sample)]
sample_multipliers = [factorial_readout(sample)]
for degree in range(1, 4):
    lowering_poly = [Fraction(0) for _ in range(degree)] + sample_derivative
    multiplier_poly = [Fraction(0) for _ in range(degree)] + sample
    previous_multiplier_poly = [Fraction(0) for _ in range(degree - 1)] + sample
    lowering_value = factorial_readout(lowering_poly)
    multiplier_value = factorial_readout(multiplier_poly)
    previous_multiplier = factorial_readout(previous_multiplier_poly)
    require(
        lowering_value == multiplier_value - degree * previous_multiplier,
        f"lowering/multiplier triangular identity failed at degree {degree}",
    )
    formula_value = sum(
        coefficient * rising(exponent, degree)
        for exponent, coefficient in zip(sample_support, sample_coefficients)
    )
    require(
        lowering_value == formula_value,
        f"direct sparse lowering formula failed at degree {degree}",
    )
    sample_observations.append(lowering_value)
    sample_multipliers.append(multiplier_value)

sample_matrix = [
    [rising(exponent, degree) for exponent in sample_support]
    for degree in range(4)
]
sample_inverse = invert_fraction(sample_matrix)
sample_decoded = [
    sum(sample_inverse[row][column] * sample_observations[column] for column in range(4))
    for row in range(4)
]
require(sample_decoded == list(sample_coefficients), "sample sparse decoder failed")
print(f"sample_support={sample_support}")
print(f"sample_observations={tuple(sample_observations)}")
print(f"sample_decoded={tuple(sample_decoded)}")

chain_controls = [
    [Fraction(0), Fraction(6), Fraction(0), Fraction(-1)],
    [Fraction(3, 5), Fraction(-7, 11), Fraction(2, 3)],
    [Fraction(-2), Fraction(0), Fraction(5, 7), Fraction(1, 13)],
]
chain_cases = 0
for control_index, polynomial in enumerate(chain_controls):
    for power in range(6):
        left = factorial_readout(
            mul_poly(pow_poly(polynomial, power), derivative(polynomial))
        )
        moment = factorial_readout(pow_poly(polynomial, power + 1))
        boundary = polynomial[0] ** (power + 1)
        right = Fraction(moment - boundary, power + 1)
        require(
            left == right,
            f"chain-lowering identity failed on control {control_index}, power {power}",
        )
        chain_cases += 1

first_moment_hostile = chain_controls[0]
require(factorial_readout(first_moment_hostile) == 0, "mean-zero hostile lost")
require(first_moment_hostile[0] == 0, "boundary-zero hostile lost")
require(
    factorial_readout(derivative(first_moment_hostile)) == 0,
    "bare first lowering should be blind",
)
require(first_moment_hostile != [Fraction(0)], "hostile became zero")
print(f"chain_identity_cases={chain_cases}")
print("bare_chain_hostile=6s-s^3:mean0_boundary0_lowering0_nonzero")

commutator_cases = 0
selector_cases = 0
quotient_digest = sha256()
for dimension in range(1, 15):
    relation = monic_laguerre(dimension)
    relation_derivative = derivative(relation)
    reduced_derivative = reduce_monic(relation_derivative, relation)
    require(
        gcd_poly(relation, relation_derivative) == [Fraction(1)],
        f"Laguerre relation not squarefree at D={dimension}",
    )

    defect_columns = []
    for basis_index in range(dimension):
        basis = [Fraction(0) for _ in range(dimension)]
        basis[basis_index] = Fraction(1)
        derivative_after_multiply = representative_derivative(
            multiply_by_u(basis, relation),
            dimension,
        )
        multiply_after_derivative = multiply_by_u(
            representative_derivative(basis, dimension),
            relation,
        )
        commutator = [
            derivative_after_multiply[index] - multiply_after_derivative[index]
            for index in range(dimension)
        ]
        expected = list(basis)
        if basis_index == dimension - 1:
            expected = [
                expected[index] - reduced_derivative[index]
                for index in range(dimension)
            ]
        require(
            commutator == expected,
            f"rank-one commutator sign failure at D={dimension}, basis={basis_index}",
        )
        defect_columns.append(
            [
                Fraction(int(row == basis_index)) - commutator[row]
                for row in range(dimension)
            ]
        )
        commutator_cases += 1

    defect_matrix = [
        [defect_columns[column][row] for column in range(dimension)]
        for row in range(dimension)
    ]
    require(
        rank_fraction(defect_matrix) == 1,
        f"commutator defect rank is not one at D={dimension}",
    )
    require(
        defect_matrix[-1][-1] == dimension,
        f"commutator trace cancellation failed at D={dimension}",
    )

    hankel = [
        [Fraction(factorial(row + column)) for column in range(dimension)]
        for row in range(dimension)
    ]
    hankel_inverse = invert_fraction(hankel)
    top_selector = [
        hankel_inverse[row][dimension - 1] for row in range(dimension)
    ]
    for monomial_degree in range(dimension):
        selector_value = sum(
            top_selector[index] * factorial(monomial_degree + index)
            for index in range(dimension)
        )
        require(
            selector_value == int(monomial_degree == dimension - 1),
            f"top coefficient selector failed at D={dimension}, monomial={monomial_degree}",
        )
        selector_cases += 1

    # The zero quotient class has representatives 0 and ell_D, whose
    # ordinary derivatives disagree by the nonzero unit ell_D' modulo ell_D.
    require(
        reduced_derivative != [Fraction(0) for _ in range(dimension)],
        f"representative-derivative hostile vanished at D={dimension}",
    )
    quotient_digest.update(
        (
            f"{dimension}:{relation}:{reduced_derivative}:"
            f"{top_selector}\n"
        ).encode("ascii")
    )

print(f"laguerre_commutator_basis_cases={commutator_cases}")
print(f"laguerre_top_selector_cases={selector_cases}")
print(f"laguerre_quotient_digest={quotient_digest.hexdigest()}")
print("representative_hostile=0~ell_D_but_0_prime!=ell_D_prime_mod_ell_D")
print("all_checks=PASS")
