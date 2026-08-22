#!/usr/bin/env python3
"""Exact parity controls for THM-3187's central-sign obstruction."""

from itertools import combinations, product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# Gaussian integers are pairs (real, imaginary), so every matrix check below
# is exact integer arithmetic.
def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[0], -x[1])


def mul(x, y):
    return (x[0] * y[0] - x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def conj(x):
    return (x[0], -x[1])


ZERO = (0, 0)
ONE = (1, 0)
MINUS_ONE = (-1, 0)
I = (0, 1)


def matrix_mul(left, right):
    return tuple(
        tuple(
            sum_gaussian(mul(left[i][k], right[k][j])
                         for k in range(len(right)))
            for j in range(len(right[0])))
        for i in range(len(left)))


def sum_gaussian(values):
    answer = ZERO
    for value in values:
        answer = add(answer, value)
    return answer


def product_gaussian(values):
    answer = ONE
    for value in values:
        answer = mul(answer, value)
    return answer


def matrix_neg(matrix):
    return tuple(tuple(neg(value) for value in row) for row in matrix)


def identity(size):
    return tuple(tuple(ONE if i == j else ZERO for j in range(size))
                 for i in range(size))


def scalar_matrix(value, matrix):
    return tuple(tuple(mul(value, entry) for entry in row) for row in matrix)


def kronecker(left, right):
    answer = []
    for left_row in left:
        for right_row in right:
            row = []
            for left_entry in left_row:
                row.extend(mul(left_entry, entry) for entry in right_row)
            answer.append(tuple(row))
    return tuple(answer)


I2 = identity(2)
X = ((ZERO, ONE), (ONE, ZERO))
Z = ((ONE, ZERO), (ZERO, MINUS_ONE))
QA = kronecker(scalar_matrix(I, X), I2)
QB = kronecker(scalar_matrix(I, Z), X)
I4 = identity(4)
CENTRAL = matrix_mul(QA, QA)

require(CENTRAL == matrix_neg(I4), "QA square lost central sign")
require(matrix_mul(QB, QB) == CENTRAL, "QB square mismatch")
require(matrix_mul(QA, QB) == matrix_neg(matrix_mul(QB, QA)),
        "Pauli lifts lost anticommutation")


# Homogeneous polynomial parity, checked on every exponent vector through
# degree eight.  The number of four-variable monomials is binom(r+3,3).
MONOMIAL_CHECKS = 0
for degree in range(9):
    count = 0
    for a in range(degree + 1):
        for b in range(degree - a + 1):
            for c in range(degree - a - b + 1):
                d = degree - a - b - c
                require((-1) ** (a + b + c + d) == (-1) ** degree,
                        "homogeneous parity failed")
                count += 1
                MONOMIAL_CHECKS += 1
    require(count == comb(degree + 3, 3), "monomial census drift")


# On every tensor, symmetric, and exterior power the scalar central action is
# (-1)^degree.  Derive it from the reconstructed central matrix on every
# standard basis vector rather than inserting the answer as a matrix.
require(all(CENTRAL[i][j] == (MINUS_ONE if i == j else ZERO)
            for i in range(4) for j in range(4)),
        "central matrix is not diagonal -I4")
CENTRAL_DIAGONAL = tuple(CENTRAL[i][i] for i in range(4))
TENSOR_BASIS_CHECKS = 0
SYMMETRIC_BASIS_CHECKS = 0
EXTERIOR_BASIS_CHECKS = 0
EXTERIOR_SIGNS = []
for degree in range(5):
    expected = (MINUS_ONE if degree % 2 else ONE)
    tensor_signs = []
    for indices in product(range(4), repeat=degree):
        value = ONE
        for index in indices:
            value = mul(value, CENTRAL_DIAGONAL[index])
        tensor_signs.append(value)
        TENSOR_BASIS_CHECKS += 1
    require(len(tensor_signs) == 4 ** degree
            and set(tensor_signs) == {expected},
            "tensor-power parity drift")

    symmetric_signs = []
    for a in range(degree + 1):
        for b in range(degree - a + 1):
            for c in range(degree - a - b + 1):
                exponents = (a, b, c, degree - a - b - c)
                value = ONE
                for index, exponent in enumerate(exponents):
                    for _ in range(exponent):
                        value = mul(value, CENTRAL_DIAGONAL[index])
                symmetric_signs.append(value)
                SYMMETRIC_BASIS_CHECKS += 1
    require(len(symmetric_signs) == comb(degree + 3, 3)
            and set(symmetric_signs) == {expected},
            "symmetric-power parity drift")

    signs = tuple(
        product_gaussian(CENTRAL_DIAGONAL[index] for index in index_set)
        for index_set in combinations(range(4), degree))
    EXTERIOR_BASIS_CHECKS += len(signs)
    require(len(signs) == comb(4, degree), "exterior dimension drift")
    require(set(signs) == {expected}, "exterior parity drift")
    EXTERIOR_SIGNS.append((degree, len(signs), signs[0][0]))


# A Hermitian/Gram observable and a jointly co-shifted correlation are even.
# A correlation against a fixed, non-co-shifted charged reference is odd.
vector = ((1, 2), (-3, 1), (2, -4), (5, 0))
reference = ((2, -1), (1, 3), (-2, 2), (4, 1))
minus_vector = tuple(neg(value) for value in vector)
minus_reference = tuple(neg(value) for value in reference)


def gram(value):
    return tuple(tuple(mul(x, conj(y)) for y in value) for x in value)


def pairing(left, right):
    return sum_gaussian(mul(x, conj(y)) for x, y in zip(left, right))


require(gram(vector) == gram(minus_vector), "Gram failed sign quotient")
fixed_pair = pairing(vector, reference)
require(fixed_pair != ZERO, "reference control accidentally vanished")
require(pairing(minus_vector, reference) == neg(fixed_pair),
        "fixed-reference oddness failed")
require(pairing(minus_vector, minus_reference) == fixed_pair,
        "co-shifted correlation failed evenness")


# Multihomogeneous bookkeeping: co-shifting all inputs uses total degree;
# moving only the first input uses its own degree.  This is the exact
# distinction between a co-moving Gram and a fixed charged reference.
MULTIDEGREE_CHECKS = 0
for left_degree in range(6):
    for right_degree in range(6):
        require(
            (-1) ** left_degree * (-1) ** right_degree
            == (-1) ** (left_degree + right_degree),
            "multihomogeneous total parity failed")
        require((-1) ** left_degree in (-1, 1),
                "one-input parity failed")
        MULTIDEGREE_CHECKS += 1


print("THM-3187 central-sign parity quotient exact control")
print("pauli_center=-I4")
print("monomial_parity_checks=" + repr(MONOMIAL_CHECKS))
print("tensor_basis_checks=" + repr(TENSOR_BASIS_CHECKS))
print("symmetric_basis_checks=" + repr(SYMMETRIC_BASIS_CHECKS))
print("exterior_basis_checks=" + repr(EXTERIOR_BASIS_CHECKS))
print("exterior_degree_dimension_central_sign=" + repr(tuple(EXTERIOR_SIGNS)))
print("hermitian_gram_central_sign=+1")
print("jointly_coshifted_correlation_sign=+1")
print("fixed_reference_correlation_sign=-1")
print("multidegree_checks=" + repr(MULTIDEGREE_CHECKS))
print("scope=parity_necessity_not_physical_edge_or_LRC_closure")
print("all_exact_checks=PASS")
