#!/usr/bin/env python3
"""Exact degree-27 squarefree specialization without a discriminant.

At target (1,1,1), this builds the first two cubic inverse algebras (dimensions
3 and 9), then interpolates the norm of the third inverse cubic from exact
multiplication determinants.  A degree-plus-one off-grid determinant is an
independent interpolation hostile.  ``require`` remains active under
``python -O``.
"""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


class QQ:
    dim = 1

    @staticmethod
    def scalar(value):
        return sp.Rational(value)

    @staticmethod
    def add(x, y):
        return x + y

    @staticmethod
    def neg(x):
        return -x

    @staticmethod
    def mul(x, y):
        return x * y

    @staticmethod
    def power(x, exponent):
        return x**exponent

    @staticmethod
    def inv(x):
        require(x != 0, "zero inverse in QQ")
        return 1 / x

    @staticmethod
    def flatten(x):
        return [sp.Rational(x)]

    @staticmethod
    def unflatten(values):
        require(len(values) == 1, "QQ unflatten dimension changed")
        return sp.Rational(values[0])

    @staticmethod
    def basis():
        return [sp.Integer(1)]


class Cubic:
    """A base algebra extended by theta^3=p*theta+q."""

    def __init__(self, base, p, q, label):
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

    def scalar(self, value):
        return (self.base.scalar(value), self.base.scalar(0), self.base.scalar(0))

    def embed(self, value):
        return (value, self.base.scalar(0), self.base.scalar(0))

    def add(self, x, y):
        return tuple(self.base.add(x[index], y[index]) for index in range(3))

    def neg(self, x):
        return tuple(self.base.neg(value) for value in x)

    def mul(self, x, y):
        convolution = [self.base.scalar(0) for _ in range(5)]
        for i in range(3):
            for j in range(3):
                convolution[i + j] = self.base.add(
                    convolution[i + j], self.base.mul(x[i], y[j])
                )
        # theta^3=p*theta+q and theta^4=p*theta^2+q*theta.
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

    def power(self, x, exponent):
        result, factor = self.one, x
        while exponent:
            if exponent & 1:
                result = self.mul(result, factor)
            exponent //= 2
            if exponent:
                factor = self.mul(factor, factor)
        return result

    def flatten(self, x):
        return sum((self.base.flatten(value) for value in x), [])

    def unflatten(self, values):
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

    def multiplication_matrix(self, x):
        columns = [self.flatten(self.mul(x, element)) for element in self.basis()]
        return sp.Matrix(
            self.dim, self.dim, lambda row, column: columns[column][row]
        )

    def inv(self, x):
        matrix = self.multiplication_matrix(x)
        target = sp.Matrix(self.flatten(self.one))
        return self.unflatten(list(matrix.inv() * target))

    def norm(self, x):
        return sp.factor(self.multiplication_matrix(x).det(method="domain-ge"))


def sub(algebra, x, y):
    return algebra.add(x, algebra.neg(y))


def divide(algebra, x, y):
    return algebra.mul(x, algebra.inv(y))


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


def L_value(algebra, a, b, c):
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


def make_extension(base, a, b, c, label):
    L = L_value(base, a, b, c)
    T = sub(base, base.scalar(4), base.mul(base.scalar(3), base.mul(b, c)))
    inverse_L = base.inv(L)
    p = base.neg(base.mul(T, inverse_L))
    q = base.mul(base.mul(base.scalar(2), c), inverse_L)
    return Cubic(base, p, q, label), L, T


def polynomial_norm_by_values(algebra, coefficients, degree):
    """Interpolate one exact norm and test it at the first unused value."""

    points = []
    for value in range(degree + 1):
        element = algebra.scalar(0)
        for power, coefficient in coefficients:
            element = algebra.add(
                element, algebra.mul(coefficient, algebra.scalar(value**power))
            )
        points.append((sp.Integer(value), algebra.norm(element)))
    X = sp.symbols("X")
    polynomial = sp.Poly(sp.interpolate(points, X), X, domain=sp.QQ)
    require(polynomial.degree() <= degree, "interpolation exceeded the norm degree bound")
    for value, expected in points:
        require(polynomial.eval(value) == expected, "interpolation grid mismatch")
    hostile_value = degree + 1
    hostile_element = algebra.scalar(0)
    for power, coefficient in coefficients:
        hostile_element = algebra.add(
            hostile_element,
            algebra.mul(coefficient, algebra.scalar(hostile_value**power)),
        )
    require(
        polynomial.eval(hostile_value) == algebra.norm(hostile_element),
        "off-grid determinant hostile failed",
    )
    return polynomial


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_RAW_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
H_LEDGER_SHA256 = "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"
P3_LEDGER_SHA256 = "fa8ba9f1cb850116c347f6e31100d1902dea3ee1c11c2b6548f7280aa9f01d50"
P3_DENOMINATOR = sp.Integer(32687327622020737269760000000000000)

raw = H_PATH.read_bytes()
require(hashlib.sha256(raw).hexdigest() == H_RAW_SHA256, "transported H pickle changed")
H = pickle.loads(raw)
a, b, c = sp.symbols("a b c")
H_poly = sp.Poly(H, a, b, c)
H_ledger = "\n".join(
    f"{monomial}:{coefficient}" for monomial, coefficient in H_poly.terms()
)
require(
    hashlib.sha256(H_ledger.encode("ascii")).hexdigest() == H_LEDGER_SHA256,
    "transported H coefficient ledger changed",
)

target = tuple(map(sp.Integer, (1, 1, 1)))
X = sp.symbols("X")
first_cubic = sp.Poly(25 * X**3 + X - 2, X)
require(first_cubic.gcd(first_cubic.diff()).degree() == 0, "first cubic is not squarefree")
K1, L0, T0 = make_extension(QQ, *target, "K1")
q1 = inverse_coordinates(K1, *(K1.scalar(value) for value in target), K1.theta)
require(fmap(K1, *q1) == tuple(K1.scalar(value) for value in target), "level-one inverse graph failed")

L1 = L_value(K1, *q1)
T1 = sub(K1, K1.scalar(4), K1.mul(K1.scalar(3), K1.mul(q1[1], q1[2])))
P2 = polynomial_norm_by_values(
    K1,
    ((3, L1), (1, T1), (0, K1.neg(K1.mul(K1.scalar(2), q1[2])))),
    9,
)
_denominator2, P2_integer = P2.clear_denoms(convert=True)
_content2, primitive2 = P2_integer.primitive()
require(
    primitive2.degree() == 9 and primitive2.gcd(primitive2.diff()).degree() == 0,
    "degree-nine core is not squarefree",
)
H111 = sp.Integer(H.subs({a: 1, b: 1, c: 1}))
require(primitive2.LC() == H111, "degree-nine leading coefficient changed")

K2, _L1_again, _T1_again = make_extension(K1, *q1, "K2")
q1_in_K2 = tuple(K2.embed(value) for value in q1)
q2 = inverse_coordinates(K2, *q1_in_K2, K2.theta)
require(fmap(K2, *q2) == q1_in_K2, "level-two inverse graph failed")

L2 = L_value(K2, *q2)
T2 = sub(K2, K2.scalar(4), K2.mul(K2.scalar(3), K2.mul(q2[1], q2[2])))
C2 = K2.neg(K2.mul(K2.scalar(2), q2[2]))
P3 = polynomial_norm_by_values(K2, ((3, L2), (1, T2), (0, C2)), 27)
denominator3, P3_integer = P3.clear_denoms(convert=True)
content3, primitive3 = P3_integer.primitive()
require(denominator3 == P3_DENOMINATOR and content3 == 1, "P3 normalization changed")
require(primitive3.degree() == 27, "P3 lost degree")
require(primitive3.gcd(primitive3.diff()).degree() == 0, "P3 is not squarefree")
ledger3 = "\n".join(
    f"{position}:{coefficient}"
    for position, coefficient in enumerate(primitive3.all_coeffs())
)
digest3 = hashlib.sha256(ledger3.encode("ascii")).hexdigest()
require(digest3 == P3_LEDGER_SHA256, "P3 coefficient ledger changed")

print("== finite-exact level-three squarefree tower ==")
print(f"transported H raw/ledger sha256={H_RAW_SHA256}/{H_LEDGER_SHA256}")
print(f"K1: E=({L0})w^3+({T0})w-2; squarefree; inverse graph PASS")
print(f"P2: degree=9 squarefree; primitive lead={primitive2.LC()}; off-grid norm PASS")
print("K2: dimension=9; inverse graph PASS")
print(
    f"P3: degree=27 squarefree; content={content3}; denominator={denominator3}; off-grid norm PASS"
)
print(f"P3 primitive ordered-coefficient sha256={digest3}")
print(f"P3 leading coefficient digits={len(str(abs(int(primitive3.LC()))))}")
print("one target certifies generic separability and pairwise block coprimality")
print("scope: no degree-27 discriminant was formed; no general JC implication")
print("all exact checks passed")
