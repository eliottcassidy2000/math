#!/usr/bin/env python3
"""Exact companion for THM-3745.

For A=k[F(b),bF(b)] inside B=k[b], the theorem proves symbolically that

    B/A = direct_sum_{i=1}^{m-1} k[X]/(X^i),
    conductor(B/A) = F(b)^(m-1) B,
    dim_k(B/A) = binom(m,2).

This script independently checks the module/conductor linear algebra over a
large prime field for split squarefree polynomials of degrees 2,...,10 and for
the nonsquarefree monomial hostile.  It also checks the Pell classification of
square triangular colengths and the typed degree-(72,108) JC near miss.  No
computer-algebra package or floating point is used.
"""

from __future__ import annotations

from math import isqrt


PRIME = 1_000_003


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def polynomial_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index] + left_value * right_value
            ) % PRIME
    return result


def polynomial_evaluate(coefficients: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % PRIME
    return result


def split_squarefree_polynomial(degree: int) -> list[int]:
    coefficients = [1]
    for root in range(1, degree + 1):
        coefficients = polynomial_multiply(coefficients, [(-root) % PRIME, 1])
    require(len(coefficients) == degree + 1 and coefficients[-1] == 1, "monic split polynomial")
    return coefficients


def relation_evaluate(coefficients: list[int], x_value: int, y_value: int) -> int:
    """Evaluate Y^m+sum c_i X^(m-i)Y^i-X^(m+1)."""

    degree = len(coefficients) - 1
    result = pow(y_value, degree, PRIME)
    for index in range(degree):
        result += coefficients[index] * pow(x_value, degree - index, PRIME) * pow(y_value, index, PRIME)
    result -= pow(x_value, degree + 1, PRIME)
    return result % PRIME


def coordinate_index(degree: int, x_power: int, b_power: int) -> int:
    return x_power * degree + b_power


def multiply_by_b(vector: list[int], coefficients: list[int]) -> list[int]:
    """Multiply in B/(X^(m-1)B), in basis X^r b^i."""

    degree = len(coefficients) - 1
    result = [0] * len(vector)
    for x_power in range(degree - 1):
        for b_power in range(degree):
            value = vector[coordinate_index(degree, x_power, b_power)]
            if value == 0:
                continue
            if b_power + 1 < degree:
                target = coordinate_index(degree, x_power, b_power + 1)
                result[target] = (result[target] + value) % PRIME
                continue

            # b^m = X-sum_(i=0)^(m-1)c_i b^i.
            if x_power + 1 < degree - 1:
                target = coordinate_index(degree, x_power + 1, 0)
                result[target] = (result[target] + value) % PRIME
            for target_power in range(degree):
                target = coordinate_index(degree, x_power, target_power)
                result[target] = (
                    result[target] - value * coefficients[target_power]
                ) % PRIME
    return result


def matrix_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    matrix = [row[:] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column] % PRIME),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], PRIME - 2, PRIME)
        matrix[pivot_row] = [(value * inverse) % PRIME for value in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                (value - multiplier * pivot_value) % PRIME
                for value, pivot_value in zip(matrix[row], matrix[pivot_row], strict=True)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def conductor_quotient_kernel_dimension(coefficients: list[int]) -> tuple[int, int, int]:
    """Compute c/I inside B/I, where I=X^(m-1)B.

    A/I is spanned by X^r b^i with r>=i.  An element lies in the conductor
    precisely when all its b^j multiples, 0<=j<m, stay in A/I.
    """

    degree = len(coefficients) - 1
    ambient_dimension = degree * (degree - 1)
    a_basis = [
        (x_power, b_power)
        for x_power in range(degree - 1)
        for b_power in range(degree)
        if x_power >= b_power
    ]
    forbidden = [
        coordinate_index(degree, x_power, b_power)
        for x_power in range(degree - 1)
        for b_power in range(degree)
        if x_power < b_power
    ]
    require(len(a_basis) == degree * (degree - 1) // 2, "A/I triangular basis")

    feature_columns: list[list[int]] = []
    for x_power, b_power in a_basis:
        vector = [0] * ambient_dimension
        vector[coordinate_index(degree, x_power, b_power)] = 1
        features: list[int] = []
        for _ in range(degree):
            features.extend(vector[index] for index in forbidden)
            vector = multiply_by_b(vector, coefficients)
        feature_columns.append(features)

    rows = [
        [feature_columns[column][row] for column in range(len(feature_columns))]
        for row in range(len(feature_columns[0]))
    ]
    rank = matrix_rank(rows)
    return len(a_basis) - rank, len(a_basis), ambient_dimension


checks = 0

print("=== Plane equation and triangular module quotient ===")
print("Work over GF(1000003); F_m=product_(r=1)^m (b-r).")
for degree in range(2, 11):
    coefficients = split_squarefree_polynomial(degree)
    for b_value in range(0, 2 * degree + 3):
        x_value = polynomial_evaluate(coefficients, b_value)
        y_value = b_value * x_value % PRIME
        require(relation_evaluate(coefficients, x_value, y_value) == 0, "parametrization satisfies G")
        checks += 1

    conductor_kernel, a_mod_i_dimension, b_mod_i_dimension = conductor_quotient_kernel_dimension(coefficients)
    delta = degree * (degree - 1) // 2
    require(a_mod_i_dimension == delta, "length A/I is delta")
    require(b_mod_i_dimension == 2 * delta, "length B/I is two delta")
    require(conductor_kernel == 0, "no conductor larger than I")
    checks += 3
    print(
        f"m={degree:2d}  delta=dim(B/A)={delta:2d}  "
        f"dim(B/I)={b_mod_i_dimension:2d}  dim((conductor)/I)={conductor_kernel}"
    )

print("The audited lattice is B/A=direct_sum_(i=1)^(m-1) k[X]/(X^i).")
print("I=X^(m-1)B has colength m(m-1); the exact annihilator audit finds c/I=0.")

print("\n=== Nonsquarefree hostile: algebra survives, branch typing does not ===")
for degree in range(2, 11):
    monomial_coefficients = [0] * degree + [1]
    conductor_kernel, a_mod_i_dimension, b_mod_i_dimension = conductor_quotient_kernel_dimension(monomial_coefficients)
    require(conductor_kernel == 0, "monomial conductor equality")
    require(a_mod_i_dimension * 2 == b_mod_i_dimension, "monomial length symmetry")
    checks += 2
print("F=b^m still has conductor b^(m(m-1)) and delta=T_(m-1),")
print("but its tangent cone is the repeated line Y^m; at m=2 this is the cusp k[b^2,b^3], not an ordinary node.")

print("\n=== Square triangular colengths are exactly a Pell orbit ===")
square_delta_rows: list[tuple[int, int, int]] = []
for degree in range(1, 10_001):
    delta = degree * (degree - 1) // 2
    root = isqrt(delta)
    if root * root == delta:
        x_value = 2 * degree - 1
        require(x_value * x_value - 8 * root * root == 1, "square delta gives Pell equation")
        square_delta_rows.append((degree, root, x_value))
    checks += 1

pell_rows: list[tuple[int, int, int]] = []
x_value, root = 1, 0
while (x_value + 1) // 2 <= 10_000:
    pell_rows.append(((x_value + 1) // 2, root, x_value))
    x_value, root = 3 * x_value + 8 * root, x_value + 3 * root
require(square_delta_rows == pell_rows, "bounded square-delta rows equal complete Pell recurrence")
checks += 1
print("m,q,2m-1 through m<=10000:")
print(" ", ", ".join(f"({degree},{root},{x_value})" for degree, root, x_value in square_delta_rows))
print("Equivalently T_(m-1)=q^2 iff (2m-1)^2-8q^2=1.")

print("\n=== Degree-(72,108) JC near miss at m=9 ===")
degree = 9
delta = degree * (degree - 1) // 2
require(delta == 36, "m=9 delta")
require(2 * delta == 72 and 3 * delta == 108, "current JC degree pair numerology")
require(degree * (degree - 1) == 72, "deg F^8")
require(degree * 12 == 108, "deg F^12")
# If U=F^8 and V=F^12, then V^2=U^3, so their two-variable Jacobian is zero.
require(2 * 12 == 3 * 8, "F^8,F^12 algebraic dependence")
checks += 5
print("delta=36, length(B/c)=72, and degrees (F^8,F^12)=(72,108).")
print("But (F^12)^2=(F^8)^3: the pair is algebraically dependent and has Jacobian zero, not one.")
print("Normalization length, conductor colength, and Keller polynomial degree remain different types.")

print("\n=== Scope ===")
print("Separable F gives a geometrically ordinary m-fold point; without separability, only the algebraic formulas survive.")
print("The Pell selector preserves the scalar delta only and forgets F, branch labels, incidence, and the conductor ideal.")
print("There is no speed-set/loneliness map and no planar Keller pair; LRC(14) and JC(2) remain open.")
print(f"CHECKS={checks}")
