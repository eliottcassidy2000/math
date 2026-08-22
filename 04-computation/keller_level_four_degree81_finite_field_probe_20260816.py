#!/usr/bin/env python3
"""Finite-field degree-81 separability gate for the fixed Keller fourth iterate.

At one target over F_101, build three nested cubic inverse algebras of
dimensions 3, 9, and 27.  The norm of the fourth inverse cubic is recovered
from 82 exact modular determinants and checked off-grid.  A squarefree
degree-81 result certifies that the generic fourth x-eliminant has full
degree and that its 27 cubic blocks are separable and pairwise coprime.
"""

from __future__ import annotations

import hashlib


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


MODULUS = 101
TARGET = (1, 1, 1)
EXPECTED_LEDGER_SHA256 = "1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7"


def mod_inv(value: int) -> int:
    value %= MODULUS
    require(value != 0, "attempted inversion of zero")
    return pow(value, MODULUS - 2, MODULUS)


def matrix_det(matrix: list[list[int]]) -> int:
    work = [[value % MODULUS for value in row] for row in matrix]
    size = len(work)
    determinant = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % MODULUS
        inverse_pivot = mod_inv(pivot_value)
        for row in range(column + 1, size):
            if work[row][column] == 0:
                continue
            scale = work[row][column] * inverse_pivot % MODULUS
            for index in range(column, size):
                work[row][index] = (work[row][index] - scale * work[column][index]) % MODULUS
    return determinant % MODULUS


def matrix_solve(matrix: list[list[int]], target: list[int]) -> list[int]:
    size = len(matrix)
    work = [
        [value % MODULUS for value in matrix[row]] + [target[row] % MODULUS]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        require(pivot is not None, "attempted inversion of a zero divisor")
        work[column], work[pivot] = work[pivot], work[column]
        inverse_pivot = mod_inv(work[column][column])
        work[column] = [(value * inverse_pivot) % MODULUS for value in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            scale = work[row][column]
            work[row] = [
                (work[row][index] - scale * work[column][index]) % MODULUS
                for index in range(size + 1)
            ]
    return [work[row][-1] for row in range(size)]


class PrimeField:
    dim = 1

    @staticmethod
    def scalar(value: int) -> int:
        return value % MODULUS

    @staticmethod
    def add(left: int, right: int) -> int:
        return (left + right) % MODULUS

    @staticmethod
    def neg(value: int) -> int:
        return (-value) % MODULUS

    @staticmethod
    def mul(left: int, right: int) -> int:
        return left * right % MODULUS

    @staticmethod
    def power(value: int, exponent: int) -> int:
        return pow(value, exponent, MODULUS)

    @staticmethod
    def inv(value: int) -> int:
        return mod_inv(value)

    @staticmethod
    def flatten(value: int) -> list[int]:
        return [value % MODULUS]

    @staticmethod
    def unflatten(values: list[int]) -> int:
        require(len(values) == 1, "prime-field unflatten dimension changed")
        return values[0] % MODULUS

    @staticmethod
    def basis() -> list[int]:
        return [1]


class Cubic:
    """A finite base algebra extended by theta^3=p*theta+q."""

    def __init__(self, base, p, q, label: str):
        self.base = base
        self.p = p
        self.q = q
        self.label = label
        self.dim = 3 * base.dim
        zero = base.scalar(0)
        one = base.scalar(1)
        self.zero = (zero, zero, zero)
        self.one = (one, zero, zero)
        self.theta = (zero, one, zero)

    def scalar(self, value: int):
        return (self.base.scalar(value), self.base.scalar(0), self.base.scalar(0))

    def embed(self, value):
        return (value, self.base.scalar(0), self.base.scalar(0))

    def add(self, left, right):
        return tuple(self.base.add(left[index], right[index]) for index in range(3))

    def neg(self, value):
        return tuple(self.base.neg(coefficient) for coefficient in value)

    def mul(self, left, right):
        convolution = [self.base.scalar(0) for _ in range(5)]
        for i in range(3):
            for j in range(3):
                convolution[i + j] = self.base.add(
                    convolution[i + j], self.base.mul(left[i], right[j])
                )
        out0 = self.base.add(convolution[0], self.base.mul(convolution[3], self.q))
        out1 = self.base.add(
            convolution[1],
            self.base.add(
                self.base.mul(convolution[3], self.p),
                self.base.mul(convolution[4], self.q),
            ),
        )
        out2 = self.base.add(convolution[2], self.base.mul(convolution[4], self.p))
        return (out0, out1, out2)

    def power(self, value, exponent: int):
        result, factor = self.one, value
        while exponent:
            if exponent & 1:
                result = self.mul(result, factor)
            exponent //= 2
            if exponent:
                factor = self.mul(factor, factor)
        return result

    def flatten(self, value) -> list[int]:
        return sum((self.base.flatten(coefficient) for coefficient in value), [])

    def unflatten(self, values: list[int]):
        require(len(values) == self.dim, f"{self.label} unflatten dimension changed")
        width = self.base.dim
        return tuple(
            self.base.unflatten(values[index * width : (index + 1) * width])
            for index in range(3)
        )

    def basis(self):
        result = []
        for outer_degree in range(3):
            for base_element in self.base.basis():
                coefficients = [self.base.scalar(0) for _ in range(3)]
                coefficients[outer_degree] = base_element
                result.append(tuple(coefficients))
        return result

    def multiplication_matrix(self, value) -> list[list[int]]:
        columns = [self.flatten(self.mul(value, element)) for element in self.basis()]
        return [
            [columns[column][row] for column in range(self.dim)]
            for row in range(self.dim)
        ]

    def inv(self, value):
        solution = matrix_solve(self.multiplication_matrix(value), self.flatten(self.one))
        return self.unflatten(solution)

    def norm(self, value) -> int:
        return matrix_det(self.multiplication_matrix(value))


def sub(algebra, left, right):
    return algebra.add(left, algebra.neg(right))


def divide(algebra, numerator, denominator):
    return algebra.mul(numerator, algebra.inv(denominator))


def fmap(algebra, x, y, z):
    one, two, three, four = map(algebra.scalar, (1, 2, 3, 4))
    xy = algebra.mul(x, y)
    unit = algebra.add(one, xy)
    four_plus = algebra.add(four, algebra.mul(three, xy))
    first = algebra.add(
        algebra.mul(algebra.power(unit, 3), z),
        algebra.mul(algebra.mul(algebra.power(y, 2), unit), four_plus),
    )
    second = algebra.add(
        y,
        algebra.add(
            algebra.mul(three, algebra.mul(algebra.mul(x, algebra.power(unit, 2)), z)),
            algebra.mul(three, algebra.mul(algebra.mul(x, algebra.power(y, 2)), four_plus)),
        ),
    )
    third = sub(
        algebra,
        sub(algebra, algebra.mul(two, x), algebra.mul(three, algebra.mul(algebra.power(x, 2), y))),
        algebra.mul(algebra.power(x, 3), z),
    )
    return first, second, third


def l_value(algebra, a, b, c):
    return algebra.add(
        algebra.add(
            sub(
                algebra,
                algebra.mul(algebra.scalar(27), algebra.mul(algebra.power(a, 2), algebra.power(c, 2))),
                algebra.mul(algebra.scalar(18), algebra.mul(algebra.mul(a, b), c)),
            ),
            algebra.mul(algebra.scalar(16), a),
        ),
        sub(algebra, algebra.mul(algebra.power(b, 3), c), algebra.power(b, 2)),
    )


def inverse_coordinates(algebra, a, b, c, x):
    two, three, nine, twelve = map(algebra.scalar, (2, 3, 9, 12))
    denominator = algebra.add(
        algebra.add(
            algebra.mul(sub(algebra, algebra.mul(twelve, a), algebra.power(b, 2)), algebra.power(x, 2)),
            algebra.mul(b, x),
        ),
        two,
    )
    y_numerator = algebra.mul(
        algebra.mul(three, algebra.mul(a, x)),
        algebra.add(algebra.mul(sub(algebra, algebra.mul(nine, algebra.mul(a, c)), b), x), two),
    )
    y = sub(algebra, b, divide(algebra, y_numerator, denominator))
    z = divide(
        algebra,
        sub(
            algebra,
            sub(algebra, algebra.mul(two, x), c),
            algebra.mul(three, algebra.mul(algebra.power(x, 2), y)),
        ),
        algebra.power(x, 3),
    )
    return x, y, z


def make_extension(base, a, b, c, label: str):
    L = l_value(base, a, b, c)
    T = sub(base, base.scalar(4), base.mul(base.scalar(3), base.mul(b, c)))
    inverse_L = base.inv(L)
    p = base.neg(base.mul(T, inverse_L))
    q = base.mul(base.mul(base.scalar(2), c), inverse_L)
    extension = Cubic(base, p, q, label)
    derivative = sub(
        extension,
        extension.mul(extension.scalar(3), extension.power(extension.theta, 2)),
        extension.embed(p),
    )
    extension.inv(derivative)
    return extension


def trim(polynomial: list[int]) -> list[int]:
    result = [value % MODULUS for value in polynomial]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add_scaled(target: list[int], source: list[int], scale: int) -> None:
    for index, coefficient in enumerate(source):
        target[index] = (target[index] + scale * coefficient) % MODULUS


def poly_multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i + j] = (result[i + j] + x * y) % MODULUS
    return trim(result)


def poly_evaluate(polynomial: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * value + coefficient) % MODULUS
    return result


def poly_derivative(polynomial: list[int]) -> list[int]:
    return trim([(index * polynomial[index]) % MODULUS for index in range(1, len(polynomial))])


def poly_remainder(numerator: list[int], denominator: list[int]) -> list[int]:
    out, divisor = trim(numerator), trim(denominator)
    require(divisor != [0], "polynomial division by zero")
    inverse_lead = mod_inv(divisor[-1])
    while out != [0] and len(out) >= len(divisor):
        shift = len(out) - len(divisor)
        scale = out[-1] * inverse_lead % MODULUS
        for index, coefficient in enumerate(divisor):
            out[index + shift] = (out[index + shift] - scale * coefficient) % MODULUS
        out = trim(out)
    return out


def poly_gcd(left: list[int], right: list[int]) -> list[int]:
    x, y = trim(left), trim(right)
    while y != [0]:
        x, y = y, poly_remainder(x, y)
    inverse_lead = mod_inv(x[-1])
    return trim([(inverse_lead * coefficient) % MODULUS for coefficient in x])


def interpolate_consecutive(values: list[int]) -> list[int]:
    """Newton interpolation at 0,1,...,d over F_p, ascending coefficients."""

    degree_bound = len(values) - 1
    require(degree_bound < MODULUS, "interpolation grid wraps modulo p")
    polynomial = [0] * (degree_bound + 1)
    basis = [1]
    differences = [value % MODULUS for value in values]
    for order in range(degree_bound + 1):
        poly_add_scaled(polynomial, basis, differences[0])
        differences = [
            (differences[index + 1] - differences[index]) % MODULUS
            for index in range(len(differences) - 1)
        ]
        if order < degree_bound:
            basis = poly_multiply(basis, [(-order) % MODULUS, 1])
            inverse_factorial_step = mod_inv(order + 1)
            basis = [(coefficient * inverse_factorial_step) % MODULUS for coefficient in basis]
    return trim(polynomial)


base = PrimeField()
target0 = tuple(base.scalar(value) for value in TARGET)
K1 = make_extension(base, *target0, "K1")
q1 = inverse_coordinates(K1, *(K1.scalar(value) for value in TARGET), K1.theta)
require(fmap(K1, *q1) == tuple(K1.scalar(value) for value in TARGET), "K1 inverse graph failed")

K2 = make_extension(K1, *q1, "K2")
q1_in_K2 = tuple(K2.embed(value) for value in q1)
q2 = inverse_coordinates(K2, *q1_in_K2, K2.theta)
require(fmap(K2, *q2) == q1_in_K2, "K2 inverse graph failed")

K3 = make_extension(K2, *q2, "K3")
q2_in_K3 = tuple(K3.embed(value) for value in q2)
q3 = inverse_coordinates(K3, *q2_in_K3, K3.theta)
require(fmap(K3, *q3) == q2_in_K3, "K3 inverse graph failed")

L3 = l_value(K3, *q3)
T3 = sub(K3, K3.scalar(4), K3.mul(K3.scalar(3), K3.mul(q3[1], q3[2])))
C3 = K3.neg(K3.mul(K3.scalar(2), q3[2]))


def norm_value(value: int) -> int:
    element = K3.add(
        K3.mul(L3, K3.scalar(value**3)),
        K3.add(K3.mul(T3, K3.scalar(value)), C3),
    )
    return K3.norm(element)


grid_values = [norm_value(value) for value in range(82)]
P4 = interpolate_consecutive(grid_values)
require(len(P4) - 1 == 81, "the fourth eliminant lost degree")
for value, expected in enumerate(grid_values):
    require(poly_evaluate(P4, value) == expected, "interpolation grid mismatch")
require(poly_evaluate(P4, 82) == norm_value(82), "off-grid determinant hostile failed")
require(poly_gcd(P4, poly_derivative(P4)) == [1], "degree-81 eliminant is not squarefree")

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(P4))
digest = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(digest == EXPECTED_LEDGER_SHA256, "degree-81 coefficient ledger changed")

print("== F_101 level-four degree-81 separability gate ==")
print(f"target={TARGET}; tower dimensions=(3,9,27)")
print("three inverse graphs: PASS; all cubic derivatives and chart denominators are units")
print("fourth norm-product: degree=81; squarefree; off-grid determinant PASS")
print(f"ascending-coefficient sha256={digest}")
print("therefore generic fourth x-eliminant degree/separability and block coprimality are nonempty")
print("scope: finite-field genericity gate only; no Norm(J) numerator or discriminant expansion")
print("all exact checks passed")
