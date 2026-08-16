#!/usr/bin/env python3
"""Independent F_101 audit of the level-four degree-81 Keller gate.

This implementation deliberately shares neither the submitted coefficient-
tuple cubic algebra, its Python Gaussian elimination, nor its consecutive-
point Newton interpolation.  Every algebra element is represented by its
FLINT regular-representation matrix.  The degree-81 norm polynomial is then
recovered from all 100 nonzero field values by multiplicative Fourier
inversion on F_101^*, with X=0 retained as an unused control.
"""

from __future__ import annotations

import hashlib

from flint import nmod_mat, nmod_poly


P = 101
TARGET = (1, 1, 1)
EXPECTED_HASH = "1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def block_matrix(blocks: list[list[nmod_mat]], modulus: int) -> nmod_mat:
    block_rows = len(blocks)
    block_columns = len(blocks[0])
    size = blocks[0][0].nrows()
    require(all(len(row) == block_columns for row in blocks), "ragged block matrix")
    require(
        all(block.nrows() == size and block.ncols() == size for row in blocks for block in row),
        "incompatible block size",
    )
    rows = []
    for outer_row in range(block_rows):
        for inner_row in range(size):
            row = []
            for outer_column in range(block_columns):
                row.extend(
                    int(blocks[outer_row][outer_column][inner_row, inner_column])
                    for inner_column in range(size)
                )
            rows.append(row)
    return nmod_mat(rows, modulus)


class RegularAlgebra:
    """A nested cubic algebra represented only by multiplication matrices."""

    def __init__(self, base=None, p_element=None, q_element=None, label="F101"):
        self.base = base
        self.label = label
        if base is None:
            self.dim = 1
            self.one = nmod_mat([[1]], P)
            self.zero = nmod_mat([[0]], P)
            self.theta = None
            return

        n = base.dim
        self.dim = 3 * n
        zero = nmod_mat(n, n, P)
        identity = base.one
        # theta*(v0+v1*theta+v2*theta^2)
        #   = q*v2 + (v0+p*v2)*theta + v1*theta^2.
        self.theta = block_matrix(
            [
                [zero, zero, q_element],
                [identity, zero, p_element],
                [zero, identity, zero],
            ],
            P,
        )
        self.one = nmod_mat(self.dim, self.dim, P)
        for index in range(self.dim):
            self.one[index, index] = 1
        self.zero = nmod_mat(self.dim, self.dim, P)

    def scalar(self, value: int) -> nmod_mat:
        return self.one * (value % P)

    def embed(self, value: nmod_mat) -> nmod_mat:
        require(self.base is not None, "prime field has no lower algebra to embed")
        require(value.nrows() == self.base.dim and value.ncols() == self.base.dim, "bad embedding")
        zero = nmod_mat(self.base.dim, self.base.dim, P)
        return block_matrix(
            [[value, zero, zero], [zero, value, zero], [zero, zero, value]],
            P,
        )

    @staticmethod
    def add(left: nmod_mat, right: nmod_mat) -> nmod_mat:
        return left + right

    @staticmethod
    def neg(value: nmod_mat) -> nmod_mat:
        return -value

    @staticmethod
    def mul(left: nmod_mat, right: nmod_mat) -> nmod_mat:
        return left * right

    def power(self, value: nmod_mat, exponent: int) -> nmod_mat:
        result, factor = self.one, value
        while exponent:
            if exponent & 1:
                result = result * factor
            exponent >>= 1
            if exponent:
                factor = factor * factor
        return result

    def inverse(self, value: nmod_mat, label: str) -> nmod_mat:
        determinant = int(value.det()) % P
        require(determinant != 0, f"nonunit encountered: {label}")
        return value.inv()

    @staticmethod
    def norm(value: nmod_mat) -> int:
        return int(value.det()) % P


def sub(algebra: RegularAlgebra, left: nmod_mat, right: nmod_mat) -> nmod_mat:
    return algebra.add(left, algebra.neg(right))


def divide(algebra: RegularAlgebra, numerator: nmod_mat, denominator: nmod_mat, label: str) -> nmod_mat:
    return algebra.mul(numerator, algebra.inverse(denominator, label))


def fmap(algebra: RegularAlgebra, x: nmod_mat, y: nmod_mat, z: nmod_mat):
    one, two, three, four = (algebra.scalar(value) for value in (1, 2, 3, 4))
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


def l_value(algebra: RegularAlgebra, a: nmod_mat, b: nmod_mat, c: nmod_mat) -> nmod_mat:
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


def inverse_coordinates(
    algebra: RegularAlgebra,
    a: nmod_mat,
    b: nmod_mat,
    c: nmod_mat,
    x: nmod_mat,
    label: str,
):
    two, three, nine, twelve = (algebra.scalar(value) for value in (2, 3, 9, 12))
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
    y = sub(algebra, b, divide(algebra, y_numerator, denominator, f"{label} y denominator"))
    z = divide(
        algebra,
        sub(
            algebra,
            sub(algebra, algebra.mul(two, x), c),
            algebra.mul(three, algebra.mul(algebra.power(x, 2), y)),
        ),
        algebra.power(x, 3),
        f"{label} x^3 denominator",
    )
    return x, y, z


def make_extension(base: RegularAlgebra, a: nmod_mat, b: nmod_mat, c: nmod_mat, label: str):
    L = l_value(base, a, b, c)
    T = sub(base, base.scalar(4), base.mul(base.scalar(3), base.mul(b, c)))
    inverse_L = base.inverse(L, f"{label} leading L")
    p_element = base.neg(base.mul(T, inverse_L))
    q_element = base.mul(base.mul(base.scalar(2), c), inverse_L)
    extension = RegularAlgebra(base, p_element, q_element, label)
    derivative = sub(
        extension,
        extension.mul(extension.scalar(3), extension.power(extension.theta, 2)),
        extension.embed(p_element),
    )
    extension.inverse(derivative, f"{label} cubic derivative")
    return extension


K0 = RegularAlgebra()
t0 = tuple(K0.scalar(value) for value in TARGET)
K1 = make_extension(K0, *t0, "K1")
q1 = inverse_coordinates(K1, *(K1.scalar(value) for value in TARGET), K1.theta, "K1")
require(fmap(K1, *q1) == tuple(K1.scalar(value) for value in TARGET), "K1 inverse graph")

K2 = make_extension(K1, *q1, "K2")
q1_up = tuple(K2.embed(value) for value in q1)
q2 = inverse_coordinates(K2, *q1_up, K2.theta, "K2")
require(fmap(K2, *q2) == q1_up, "K2 inverse graph")

K3 = make_extension(K2, *q2, "K3")
q2_up = tuple(K3.embed(value) for value in q2)
q3 = inverse_coordinates(K3, *q2_up, K3.theta, "K3")
require(fmap(K3, *q3) == q2_up, "K3 inverse graph")

L3 = l_value(K3, *q3)
T3 = sub(K3, K3.scalar(4), K3.mul(K3.scalar(3), K3.mul(q3[1], q3[2])))
C3 = K3.neg(K3.mul(K3.scalar(2), q3[2]))


def norm_at(value: int) -> int:
    element = K3.add(
        K3.mul(L3, K3.scalar(pow(value, 3, P))),
        K3.add(K3.mul(T3, K3.scalar(value)), C3),
    )
    return K3.norm(element)


# Multiplicative Fourier inversion on F_101^*: for 0<=k<100,
# c_k=(1/100) sum_x P(x)x^{-k}=-sum_x P(x)x^{-k} mod 101.
nonzero_values = {value: norm_at(value) for value in range(1, P)}
coefficients_0_to_99 = []
for exponent in range(P - 1):
    character_sum = sum(
        polynomial_value * pow(value, (-exponent) % (P - 1), P)
        for value, polynomial_value in nonzero_values.items()
    ) % P
    coefficients_0_to_99.append((-character_sum) % P)

require(all(value == 0 for value in coefficients_0_to_99[82:]), "Fourier tail above degree 81 is nonzero")
coefficients = coefficients_0_to_99[:82]
require(coefficients[-1] != 0, "degree-81 leading coefficient vanished")
require(coefficients[0] == norm_at(0), "unused X=0 determinant control failed")

polynomial = nmod_poly(coefficients, P)
require(polynomial.degree() == 81, "FLINT polynomial degree changed")
require(polynomial.gcd(polynomial.derivative()) == nmod_poly([1], P), "degree-81 polynomial is not squarefree")

# A direct Horner pass is independent of the Fourier recovery formula and
# checks every determinant value used as spectral input.
for value, expected in nonzero_values.items():
    require(int(polynomial(value)) % P == expected, f"Horner mismatch at X={value}")

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(coefficients))
digest = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(digest == EXPECTED_HASH, "independent degree-81 ledger disagrees with candidate")

print("== independent F_101 level-four Fourier/FLINT audit ==")
print(f"target={TARGET}; regular-matrix dimensions=(3,9,27)")
print("inverse graphs and all leading/derivative/chart units: PASS")
print("100-point multiplicative Fourier inversion: degrees 82..99 vanish; degree=81")
print("unused X=0 determinant and all-point Horner controls: PASS")
print("FLINT derivative gcd=1; the 81-root norm-product is squarefree")
print(f"ascending-coefficient sha256={digest}")
print("scope: a good-reduction genericity witness, not a factorization or image theorem for G")
print("all independent exact checks passed")
