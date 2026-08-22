#!/usr/bin/env python3
"""Exact audit for the constant (b,c,k) extension of THM-3212.

All characteristic-zero arithmetic is performed directly in Q[u]/(q)[x].
Independent finite-field reductions certify the squarefree controls.  No
repository executable is imported.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import ClassVar
import hashlib
import math

import sympy as sp


Q = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass(frozen=True)
class Alg:
    c: tuple[Q, Q, Q]
    relation: ClassVar[tuple[Q, Q, Q]] = (Q(0), Q(0), Q(0))

    @staticmethod
    def make(a=0, b=0, c=0) -> "Alg":
        return Alg((Q(a), Q(b), Q(c)))

    def __add__(self, other) -> "Alg":
        other = as_alg(other)
        return Alg(tuple(self.c[i] + other.c[i] for i in range(3)))

    __radd__ = __add__

    def __neg__(self) -> "Alg":
        return Alg(tuple(-entry for entry in self.c))

    def __sub__(self, other) -> "Alg":
        return self + (-as_alg(other))

    def __rsub__(self, other) -> "Alg":
        return as_alg(other) - self

    def __mul__(self, other) -> "Alg":
        other = as_alg(other)
        raw = [Q(0)] * 5
        for i, left in enumerate(self.c):
            for j, right in enumerate(other.c):
                raw[i + j] += left * right
        r0, r1, r2 = Alg.relation
        for degree in (4, 3):
            value = raw[degree]
            raw[degree] = Q(0)
            raw[degree - 3] += value * r0
            raw[degree - 2] += value * r1
            raw[degree - 1] += value * r2
        return Alg(tuple(raw[:3]))

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Alg":
        require(exponent >= 0, "negative algebra exponent")
        answer, base, power = ONE_A, self, exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def inverse(self) -> "Alg":
        require(self != ZERO_A, "division by zero")
        basis = (ONE_A, U_A, U_A * U_A)
        columns = [(self * item).c for item in basis]
        matrix = [
            [columns[column][row] for column in range(3)]
            + [Q(1 if row == 0 else 0)]
            for row in range(3)
        ]
        for pivot in range(3):
            row = next((i for i in range(pivot, 3) if matrix[i][pivot]), None)
            require(row is not None, "nonunit accessory element")
            matrix[pivot], matrix[row] = matrix[row], matrix[pivot]
            scale = matrix[pivot][pivot]
            matrix[pivot] = [entry / scale for entry in matrix[pivot]]
            for i in range(3):
                if i == pivot:
                    continue
                scale = matrix[i][pivot]
                if scale:
                    matrix[i] = [
                        matrix[i][j] - scale * matrix[pivot][j]
                        for j in range(4)
                    ]
        return Alg(tuple(matrix[i][3] for i in range(3)))

    def __truediv__(self, other) -> "Alg":
        return self * as_alg(other).inverse()

    def __rtruediv__(self, other) -> "Alg":
        return as_alg(other) / self

    def is_zero(self) -> bool:
        return self == ZERO_A


def as_alg(value) -> Alg:
    return value if isinstance(value, Alg) else Alg.make(value)


ZERO_A = Alg.make(0)
ONE_A = Alg.make(1)
U_A = Alg.make(0, 1)


@dataclass(frozen=True)
class Poly:
    c: tuple[Alg, ...]

    @staticmethod
    def make(coefficients) -> "Poly":
        entries = [as_alg(entry) for entry in coefficients]
        while entries and entries[-1].is_zero():
            entries.pop()
        return Poly(tuple(entries))

    @property
    def degree(self) -> int:
        return len(self.c) - 1

    @property
    def lc(self) -> Alg:
        require(self.degree >= 0, "zero polynomial has no leading coefficient")
        return self.c[-1]

    def is_zero(self) -> bool:
        return self.degree < 0

    def __add__(self, other) -> "Poly":
        other = as_poly(other)
        size = max(len(self.c), len(other.c))
        return Poly.make(
            (self.c[i] if i < len(self.c) else ZERO_A)
            + (other.c[i] if i < len(other.c) else ZERO_A)
            for i in range(size)
        )

    __radd__ = __add__

    def __neg__(self) -> "Poly":
        return Poly.make(-entry for entry in self.c)

    def __sub__(self, other) -> "Poly":
        return self + (-as_poly(other))

    def __rsub__(self, other) -> "Poly":
        return as_poly(other) - self

    def __mul__(self, other) -> "Poly":
        other = as_poly(other)
        if self.is_zero() or other.is_zero():
            return ZERO_P
        entries = [ZERO_A] * (len(self.c) + len(other.c) - 1)
        for i, left in enumerate(self.c):
            for j, right in enumerate(other.c):
                entries[i + j] = entries[i + j] + left * right
        return Poly.make(entries)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "Poly":
        require(exponent >= 0, "negative polynomial exponent")
        answer, base, power = ONE_P, self, exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer

    def scalar(self, value) -> "Poly":
        factor = as_alg(value)
        return Poly.make(entry * factor for entry in self.c)

    def derivative(self) -> "Poly":
        return Poly.make(self.c[i] * i for i in range(1, len(self.c)))

    def divmod(self, divisor: "Poly") -> tuple["Poly", "Poly"]:
        require(not divisor.is_zero(), "division by zero polynomial")
        remainder = list(self.c)
        quotient = [ZERO_A] * max(0, self.degree - divisor.degree + 1)
        inverse_lc = divisor.lc.inverse()
        while remainder and len(remainder) - 1 >= divisor.degree:
            shift = len(remainder) - 1 - divisor.degree
            factor = remainder[-1] * inverse_lc
            quotient[shift] = factor
            for index, entry in enumerate(divisor.c):
                remainder[index + shift] = remainder[index + shift] - factor * entry
            while remainder and remainder[-1].is_zero():
                remainder.pop()
        return Poly.make(quotient), Poly.make(remainder)

    def exact_quotient(self, divisor: "Poly") -> "Poly":
        quotient, remainder = self.divmod(divisor)
        require(remainder.is_zero(), "expected exact polynomial quotient")
        return quotient

    def monic(self) -> "Poly":
        return self if self.is_zero() else self.scalar(self.lc.inverse())


def as_poly(value) -> Poly:
    return value if isinstance(value, Poly) else Poly.make([value])


ZERO_P = Poly.make([])
ONE_P = Poly.make([1])
X_P = Poly.make([0, 1])


def gcd(left: Poly, right: Poly) -> Poly:
    while not right.is_zero():
        _, remainder = left.divmod(right)
        left, right = right, remainder
    return left.monic()


def evaluate(poly: Poly, value: Alg) -> Alg:
    answer = ZERO_A
    for entry in reversed(poly.c):
        answer = answer * value + entry
    return answer


def serialize(poly: Poly) -> bytes:
    rows = []
    for entry in poly.c:
        rows.append(",".join(f"{q.numerator}/{q.denominator}" for q in entry.c))
    return (";".join(rows) + "\n").encode()


def configure(c3: int, c2: int, c1: int, c0: int) -> None:
    Alg.relation = (Q(-c0, c3), Q(-c1, c3), Q(-c2, c3))


def Kpoly(V: Poly, A: Poly, dV: Poly, b, c, k) -> Poly:
    b, c, k = as_alg(b), as_alg(c), as_alg(k)
    return (
        (V**6).scalar(64 * k**3)
        + (V**5).scalar(-96 * b * k**2 + 64 * c * k)
        + (V**4 * dV).scalar(32 * b**2 * c * k**2)
        + (V**4).scalar(32 * b**2 * k - 32 * b * c)
        + (V**4 * A * dV).scalar(48 * b * k**2)
        - 16 * V**4 * A
        + (V**3 * dV).scalar(-8 * b**4 * k**2 - 32 * b**3 * c * k + 32 * b**2 * c**2)
        + (V**3 * A * dV).scalar(-24 * b**2 * k + 16 * b * c)
        + (V**2 * dV).scalar(8 * b**5 * k - 8 * b**4 * c)
        + (V**2 * A * dV**2).scalar(16 * b**3 * c * k)
        - (V**2 * A * dV).scalar(4 * b**3)
        + (V**2 * A**2 * dV**2).scalar(12 * b**2 * k)
        - (V * A * dV**2).scalar(4 * b**5 * k)
        + (A**2 * dV**3).scalar(2 * b**4 * c)
        + (A**3 * dV**3).scalar(b**3)
    )


def universal_parameter_coefficients(V: Poly, A: Poly, dV: Poly):
    entries: dict[tuple[int, int, int], Poly] = {}

    def add(key, value):
        entries[key] = entries.get(key, ZERO_P) + value

    add((0, 0, 3), 64 * V**6)
    add((1, 0, 2), -96 * V**5 + 48 * V**4 * A * dV)
    add((0, 1, 1), 64 * V**5)
    add((2, 1, 2), 32 * V**4 * dV)
    add((2, 0, 1), 32 * V**4 - 24 * V**3 * A * dV + 12 * V**2 * A**2 * dV**2)
    add((1, 1, 0), -32 * V**4 + 16 * V**3 * A * dV)
    add((0, 0, 0), -16 * V**4 * A)
    add((4, 0, 2), -8 * V**3 * dV)
    add((3, 1, 1), -32 * V**3 * dV + 16 * V**2 * A * dV**2)
    add((2, 2, 0), 32 * V**3 * dV)
    add((5, 0, 1), 8 * V**2 * dV - 4 * V * A * dV**2)
    add((4, 1, 0), -8 * V**2 * dV + 2 * A**2 * dV**3)
    add((3, 0, 0), -4 * V**2 * A * dV + A**3 * dV**3)
    return entries


def universal_resultant_check() -> None:
    y, V, W, dV, b, c, k = sp.symbols("y V W dV b c k")
    L = y**2 + b * y + c * V
    R1 = 2 * L * (2 * y + b) + W
    R2 = k * V**3 + V**2 * y - b * dV * y * L
    derived = sp.expand(sp.resultant(R1, R2, y))
    A = sp.symbols("A")
    quotient = sp.cancel(derived.subs(W, V * A) / V**3)
    expected = (
        64 * k**3 * V**6
        + (-96 * b * k**2 + 64 * c * k) * V**5
        + 32 * b**2 * c * k**2 * V**4 * dV
        + (32 * b**2 * k - 32 * b * c) * V**4
        + 48 * b * k**2 * V**4 * A * dV
        - 16 * V**4 * A
        + (-8 * b**4 * k**2 - 32 * b**3 * c * k + 32 * b**2 * c**2) * V**3 * dV
        + (-24 * b**2 * k + 16 * b * c) * V**3 * A * dV
        + (8 * b**5 * k - 8 * b**4 * c) * V**2 * dV
        + 16 * b**3 * c * k * V**2 * A * dV**2
        - 4 * b**3 * V**2 * A * dV
        + 12 * b**2 * k * V**2 * A**2 * dV**2
        - 4 * b**5 * k * V * A * dV**2
        + 2 * b**4 * c * A**2 * dV**3
        + b**3 * A**3 * dV**3
    )
    require(sp.expand(quotient - expected) == 0, "universal resultant differs")
    print("universal_resultant=PASS")


def local_jet_symbolic_check() -> None:
    """Derive the S-wall coefficient and the cubic/quintic escape pair."""
    t, b, c, k, r = sp.symbols("t b c k r")
    v1, v2, v3, v4, v5, v6 = sp.symbols("v1:7")
    a2, a3, a4, a5, a6 = sp.symbols("a2:7")
    V = v1 * t + v2 * t**2 + v3 * t**3 + v4 * t**4 + v5 * t**5 + v6 * t**6
    A = 2 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5 + a6 * t**6
    dV = sp.diff(V, t)
    response = sp.expand(2 * V * sp.diff(A, t) - A * dV - 2 * V)
    solved = {}
    for degree, variable in zip(range(2, 7), (a2, a3, a4, a5, a6)):
        equation = sp.expand(response.subs(solved)).coeff(t, degree)
        solved[variable] = sp.solve(equation, variable)[0]
    A = sp.expand(A.subs(solved))
    K = (
        64 * k**3 * V**6
        + (-96 * b * k**2 + 64 * c * k) * V**5
        + 32 * b**2 * c * k**2 * V**4 * dV
        + (32 * b**2 * k - 32 * b * c) * V**4
        + 48 * b * k**2 * V**4 * A * dV
        - 16 * V**4 * A
        + (-8 * b**4 * k**2 - 32 * b**3 * c * k + 32 * b**2 * c**2) * V**3 * dV
        + (-24 * b**2 * k + 16 * b * c) * V**3 * A * dV
        + (8 * b**5 * k - 8 * b**4 * c) * V**2 * dV
        + 16 * b**3 * c * k * V**2 * A * dV**2
        - 4 * b**3 * V**2 * A * dV
        + 12 * b**2 * k * V**2 * A**2 * dV**2
        - 4 * b**5 * k * V * A * dV**2
        + 2 * b**4 * c * A**2 * dV**3
        + b**3 * A**3 * dV**3
    )
    coefficient3 = sp.expand(K).coeff(t, 3)
    expected3 = (
        -sp.Rational(8, 3)
        * b**2
        * v1**2
        * (b * k - 2 * c)
        * (4 * b**2 * v2 + 3 * (b * k + 2 * c) * v1**2)
    )
    require(sp.factor(coefficient3 - expected3) == 0, "S-wall t^3 coefficient")

    eta = -2 * v2 / (3 * v1**2)
    tau_K = sp.expand(K.subs({c: r * b**2, k: b * (2 * eta - 2 * r)}))
    L0 = 3 * r * v1**2 + v2
    A2 = 30 * r**2 * v1**4 + 10 * r * v1**2 * v2 + 18 * v1 * v3 - 16 * v2**2
    coefficient4 = sp.cancel(tau_K.coeff(t, 4))
    expected4 = sp.Rational(64, 45) * b**3 / v1 * L0 * (b**3 * A2 - 15 * v1**3)
    require(sp.factor(coefficient4 - expected4) == 0, "S-wall t^4 coefficient")

    J3 = (
        315 * r**3 * v1**6
        + 105 * r**2 * v1**4 * v2
        + r * (315 * v1**3 * v3 - 280 * v1**2 * v2**2)
        - 90 * v1**2 * v4
        + 219 * v1 * v2 * v3
        - 128 * v2**3
    )
    J5 = (
        113400 * r**5 * v1**10
        + 354375 * r**4 * v1**8 * v2
        + r**3 * (59535 * v1**7 * v3 + 148680 * v1**6 * v2**2)
        + r**2 * (12150 * v1**6 * v4 + 256770 * v1**5 * v2 * v3 - 204165 * v1**4 * v2**3)
        + r
        * (
            -15750 * v1**5 * v5
            - 29700 * v1**4 * v2 * v4
            + 11529 * v1**4 * v3**2
            + 202554 * v1**3 * v2**2 * v3
            - 157474 * v1**2 * v2**4
        )
        - 5250 * v1**3 * v2 * v5
        - 7290 * v1**3 * v3 * v4
        - 4770 * v1**2 * v2**2 * v4
        + 13077 * v1**2 * v2 * v3**2
        + 30177 * v1 * v2**3 * v3
        - 25712 * v2**5
    )
    h_value = 15 * v1**3 / A2

    def replace_b_cubes(expression):
        numerator, denominator = sp.fraction(sp.cancel(expression))
        polynomial = sp.Poly(numerator, b, r)
        minimum = min(monomial[0] for monomial, _ in polynomial.terms())
        require(
            all((monomial[0] - minimum) % 3 == 0 for monomial, _ in polynomial.terms()),
            "unexpected b grading",
        )
        replaced = sum(
            coefficient * h_value ** ((power_b - minimum) // 3) * r**power_r
            for (power_b, power_r), coefficient in polynomial.terms()
        )
        return sp.cancel(b**minimum * replaced / denominator)

    coefficient5 = replace_b_cubes(tau_K.coeff(t, 5))
    coefficient6 = replace_b_cubes(tau_K.coeff(t, 6))
    expected5 = -sp.Rational(640, 7) * v1**4 * L0 * J3 / A2**2
    expected6 = -sp.Rational(128, 63) * v1**3 * J5 / A2**2
    require(sp.factor(coefficient5 - expected5) == 0, "S-wall cubic coefficient")
    require(sp.factor(coefficient6 - expected6) == 0, "S-wall quintic coefficient")
    print("local_jet_elimination=PASS")


def jet_polynomials(v1, v2, v3, v4, v5) -> tuple[Poly, Poly]:
    J3 = (
        (X_P**3).scalar(315 * v1**6)
        + (X_P**2).scalar(105 * v1**4 * v2)
        + X_P.scalar(315 * v1**3 * v3 - 280 * v1**2 * v2**2)
        + ONE_P.scalar(-90 * v1**2 * v4 + 219 * v1 * v2 * v3 - 128 * v2**3)
    )
    J5 = (
        (X_P**5).scalar(113400 * v1**10)
        + (X_P**4).scalar(354375 * v1**8 * v2)
        + (X_P**3).scalar(59535 * v1**7 * v3 + 148680 * v1**6 * v2**2)
        + (X_P**2).scalar(
            12150 * v1**6 * v4
            + 256770 * v1**5 * v2 * v3
            - 204165 * v1**4 * v2**3
        )
        + X_P.scalar(
            -15750 * v1**5 * v5
            - 29700 * v1**4 * v2 * v4
            + 11529 * v1**4 * v3**2
            + 202554 * v1**3 * v2**2 * v3
            - 157474 * v1**2 * v2**4
        )
        + ONE_P.scalar(
            -5250 * v1**3 * v2 * v5
            - 7290 * v1**3 * v3 * v4
            - 4770 * v1**2 * v2**2 * v4
            + 13077 * v1**2 * v2 * v3**2
            + 30177 * v1 * v2**3 * v3
            - 25712 * v2**5
        )
    )
    return J3, J5


def prs_degrees(left: Poly, right: Poly):
    degrees = []
    while not right.is_zero():
        degrees.append((left.degree, right.degree))
        _, remainder = left.divmod(right)
        left, right = right, remainder
    return tuple(degrees), left.monic()


def build_case(name: str):
    if name == "4111":
        configure(100, 244, 237, 44)
        u = U_A
        a, b = 4, 1
        v = (8 * u**2 + 9 * u + 8) / 7
        s0 = 5 * (u + 1) / 7
        A0 = 80 * v**2 * (u + 1) / 343
        extras = (9, 0)
        spectrum = ((Q(2), 1), (Q(-1, 2), 1), (Q(-2), 3))
        eta_expected = Alg(
            (Q(-4786369, 131072), Q(-2131157, 32768), Q(-873425, 32768))
        )
    else:
        configure(75, -89, -31, 61)
        u = U_A
        a, b = 3, 2
        v = (24 * u**2 - 16 * u - 16) / 21
        s0 = (5 * u - 4) / 7
        A0 = 9 * v**2 * (5 * u - 4) / 343
        extras = (6, 3)
        spectrum = ((Q(2), 1), (Q(-2, 3), 1), (Q(-1), 1), (Q(-2), 2))
        eta_expected = Alg(
            (Q(-139515593, 1679616), Q(340202737, 3359232), Q(50530025, 3359232))
        )
    C = -7 * A0
    q2 = X_P**2 - X_P.scalar(u) + v
    D = X_P**a * (X_P - 1) ** b * q2
    T = X_P * (X_P - 1) * q2
    E = (
        (X_P - 1) * q2 * a
        + X_P * q2 * b
        + X_P * (X_P - 1) * (2 * X_P - u)
    ).scalar(Q(1, 7))
    S = X_P + s0
    V = (S * D * T**2).scalar(4 / C**2)
    A = (S * E * T).scalar(2 / C)
    dV = V.derivative()
    g = S * T
    boundary = S**3 * T**8 * X_P**extras[0] * (X_P - 1) ** extras[1]
    require(2 * V * A.derivative() - A * dV == 2 * V, "response derivative identity")
    require(V.degree == 16 and A.degree == 8 and boundary.degree == 44, "degree data")
    for value, multiplicity in spectrum:
        require(gcd(A.derivative() - value, g).degree == multiplicity, "A' spectrum")
    jets = []
    derivative = V
    for order in range(1, 6):
        derivative = derivative.derivative()
        jets.append(evaluate(derivative, -s0) / math.factorial(order))
    v1, v2, v3, v4, v5 = jets
    eta = -2 * v2 / (3 * v1**2)
    require(eta == eta_expected, "eta formula")
    coefficients = universal_parameter_coefficients(V, A, dV)
    quotient_degrees = []
    for key in sorted(coefficients):
        quotient_degrees.append((key, coefficients[key].exact_quotient(boundary).degree))
    require(len(quotient_degrees) == 13, "universal coefficient count")
    J3, J5 = jet_polynomials(v1, v2, v3, v4, v5)
    degrees, common = prs_degrees(J5, J3)
    require(degrees == ((5, 3), (3, 2), (2, 1), (1, 0)), "jet PRS profile")
    require(common == ONE_P, "jet cubic/quintic gcd")
    controls = []
    for label, bb, cc, kk, expected_degree, expected_gcdg in (
        ("generic_52", 1, 0, 1, 52, 0),
        ("k0_28", 1, 1, 0, 28, 0),
        ("k0_S_escape", 1, eta, 0, 28, 1),
        ("k1_S_escape", 1, 0, 2 * eta, 52, 1),
    ):
        H = Kpoly(V, A, dV, bb, cc, kk).exact_quotient(boundary)
        require(H.degree == expected_degree, f"{label} degree")
        require(gcd(H, g).degree == expected_gcdg, f"{label} boundary gcd")
        if expected_degree == 28:
            require(gcd(H, H.derivative()).degree == 0, f"{label} squarefree")
        controls.append((label, H.degree, expected_gcdg, hashlib.sha256(serialize(H)).hexdigest()))
    print(f"case={name}")
    print(f"A_prime_spectrum={tuple((str(value), mult) for value, mult in spectrum)}")
    print(f"eta_coefficients={eta.c}")
    print(f"universal_parameter_coefficient_quotient_degrees={tuple(quotient_degrees)}")
    print(f"jet_PRS_degrees={degrees} jet_gcd=1")
    print(f"J3_sha256={hashlib.sha256(serialize(J3)).hexdigest()}")
    print(f"J5_sha256={hashlib.sha256(serialize(J5)).hexdigest()}")
    for row in controls:
        print(f"control={row}")


def finite_case(name: str) -> None:
    x = sp.symbols("x")
    if name == "4111":
        p, u, a, b, eta = 113, 85, 4, 1, 101
        v = (8 * u**2 + 9 * u + 8) * pow(7, -1, p) % p
        s0 = 5 * (u + 1) * pow(7, -1, p) % p
        A0 = 80 * v**2 * (u + 1) * pow(343, -1, p) % p
        extras = (9, 0)
    else:
        p, u, a, b, eta = 101, 64, 3, 2, 43
        v = (24 * u**2 - 16 * u - 16) * pow(21, -1, p) % p
        s0 = (5 * u - 4) * pow(7, -1, p) % p
        A0 = 9 * v**2 * (5 * u - 4) * pow(343, -1, p) % p
        extras = (6, 3)
    C = -7 * A0 % p
    P = lambda expression: sp.Poly(expression, x, modulus=p)
    q2 = P(x**2 - u * x + v)
    D = P(x**a * (x - 1) ** b) * q2
    T = P(x * (x - 1)) * q2
    E = P(
        (
            a * (x - 1) * (x**2 - u * x + v)
            + b * x * (x**2 - u * x + v)
            + x * (x - 1) * (2 * x - u)
        )
        * pow(7, -1, p)
    )
    S = P(x + s0)
    V = (S * D * T**2).mul_ground(4 * pow(C**2, -1, p) % p)
    A = (S * E * T).mul_ground(2 * pow(C, -1, p) % p)
    dV = V.diff()
    g = S * T
    boundary = S**3 * T**8 * P(x) ** extras[0] * P(x - 1) ** extras[1]

    def finite_K(bb, cc, kk):
        return (
            (V**6).mul_ground(64 * kk**3)
            + (V**5).mul_ground(-96 * bb * kk**2 + 64 * cc * kk)
            + (V**4 * dV).mul_ground(32 * bb**2 * cc * kk**2)
            + (V**4).mul_ground(32 * bb**2 * kk - 32 * bb * cc)
            + (V**4 * A * dV).mul_ground(48 * bb * kk**2)
            - (V**4 * A).mul_ground(16)
            + (V**3 * dV).mul_ground(-8 * bb**4 * kk**2 - 32 * bb**3 * cc * kk + 32 * bb**2 * cc**2)
            + (V**3 * A * dV).mul_ground(-24 * bb**2 * kk + 16 * bb * cc)
            + (V**2 * dV).mul_ground(8 * bb**5 * kk - 8 * bb**4 * cc)
            + (V**2 * A * dV**2).mul_ground(16 * bb**3 * cc * kk)
            - (V**2 * A * dV).mul_ground(4 * bb**3)
            + (V**2 * A**2 * dV**2).mul_ground(12 * bb**2 * kk)
            - (V * A * dV**2).mul_ground(4 * bb**5 * kk)
            + (A**2 * dV**3).mul_ground(2 * bb**4 * cc)
            + (A**3 * dV**3).mul_ground(bb**3)
        )

    rows = []
    for label, cc, kk, degree, boundary_degree in (
        ("generic_52", 0, 1, 52, 0),
        ("k0_28", 1, 0, 28, 0),
        ("k0_S_escape", eta, 0, 28, 1),
        ("k1_S_escape", 0, 2 * eta, 52, 1),
    ):
        H = sp.exquo(finite_K(1, cc, kk), boundary)
        require(H.degree() == degree, "finite degree")
        require(sp.gcd(H, g).degree() == boundary_degree, "finite boundary gcd")
        require(sp.gcd(H, H.diff()).degree() == 0, "finite squarefree")
        rows.append((label, degree, boundary_degree, 0))
    print(f"finite_case={name} good_reduction=(p={p},u={u},eta={eta}) rows={tuple(rows)}")


def main() -> None:
    universal_resultant_check()
    local_jet_symbolic_check()
    build_case("4111")
    build_case("3211")
    finite_case("4111")
    finite_case("3211")
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
