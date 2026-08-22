#!/usr/bin/env python3
"""Exact companion for the THM-3237 degree-nine infinity-wall candidate.

The family is

    P_t(x,z) = (V(x) z^2 + (1+t*x^9) z)^2 + A(x) z + x.

Characteristic-zero calculations take place directly in each cubic
accessory algebra Q[u]/(q).  The two finite-field calculations are good
reductions used as squarefreeness and boundary-disjointness certificates.
No repository executable is imported and no assertion statement is used.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd as integer_gcd
from typing import ClassVar

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
        columns = [(self * item).c for item in (ONE_A, U_A, U_A * U_A)]
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


def as_poly(value) -> Poly:
    return value if isinstance(value, Poly) else Poly.make([value])


ZERO_P = Poly.make([])
ONE_P = Poly.make([1])
X_P = Poly.make([0, 1])


@dataclass(frozen=True)
class ParamPoly:
    """Polynomial in t whose coefficients are accessory polynomials in x."""

    c: tuple[Poly, ...]

    @staticmethod
    def make(coefficients) -> "ParamPoly":
        entries = [as_poly(entry) for entry in coefficients]
        while entries and entries[-1].is_zero():
            entries.pop()
        return ParamPoly(tuple(entries))

    def is_zero(self) -> bool:
        return not self.c

    def __add__(self, other) -> "ParamPoly":
        other = as_param(other)
        size = max(len(self.c), len(other.c))
        return ParamPoly.make(
            (self.c[i] if i < len(self.c) else ZERO_P)
            + (other.c[i] if i < len(other.c) else ZERO_P)
            for i in range(size)
        )

    __radd__ = __add__

    def __neg__(self) -> "ParamPoly":
        return ParamPoly.make(-entry for entry in self.c)

    def __sub__(self, other) -> "ParamPoly":
        return self + (-as_param(other))

    def __rsub__(self, other) -> "ParamPoly":
        return as_param(other) - self

    def __mul__(self, other) -> "ParamPoly":
        other = as_param(other)
        if self.is_zero() or other.is_zero():
            return ZERO_T
        entries = [ZERO_P] * (len(self.c) + len(other.c) - 1)
        for i, left in enumerate(self.c):
            for j, right in enumerate(other.c):
                entries[i + j] = entries[i + j] + left * right
        return ParamPoly.make(entries)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "ParamPoly":
        require(exponent >= 0, "negative parameter exponent")
        answer, base, power = ONE_T, self, exponent
        while power:
            if power & 1:
                answer = answer * base
            base = base * base
            power //= 2
        return answer


def as_param(value) -> ParamPoly:
    return value if isinstance(value, ParamPoly) else ParamPoly.make([value])


ZERO_T = ParamPoly.make([])
ONE_T = ParamPoly.make([ONE_P])
T_T = ParamPoly.make([ZERO_P, ONE_P])


def configure(c3: int, c2: int, c1: int, c0: int) -> None:
    Alg.relation = (Q(-c0, c3), Q(-c1, c3), Q(-c2, c3))


def multiplication_matrix(value: Alg) -> tuple[tuple[Q, ...], ...]:
    columns = [(value * item).c for item in (ONE_A, U_A, U_A * U_A)]
    return tuple(tuple(columns[column][row] for column in range(3)) for row in range(3))


def determinant3(matrix) -> Q:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def norm(value: Alg) -> Q:
    return determinant3(multiplication_matrix(value))


def primitive_wall_norm(gamma: Alg) -> tuple[int, ...]:
    t = sp.symbols("t")
    matrix = multiplication_matrix(gamma)
    symbolic = sp.Matrix([
        [
            sp.Rational(matrix[i][j].numerator, matrix[i][j].denominator) * t
            - (2 if i == j else 0)
            for j in range(3)
        ]
        for i in range(3)
    ])
    coefficients = sp.Poly(sp.expand(symbolic.det()), t).all_coeffs()
    denominators = [int(entry.q) for entry in coefficients]
    common_denominator = 1
    for denominator in denominators:
        common_denominator = common_denominator * denominator // integer_gcd(common_denominator, denominator)
    integers = [int(entry * common_denominator) for entry in coefficients]
    content = 0
    for entry in integers:
        content = integer_gcd(content, abs(entry))
    integers = [entry // content for entry in integers]
    if integers[0] < 0:
        integers = [-entry for entry in integers]
    return tuple(integers)


def serialize(poly: Poly) -> bytes:
    rows = []
    for entry in poly.c:
        rows.append(",".join(f"{q.numerator}/{q.denominator}" for q in entry.c))
    return (";".join(rows) + "\n").encode()


def K_general(V, A, B, J):
    """The V^(-3)-saturated y-resultant for C_0=0 and E_0=x."""
    return (
        -A**3 * J**3
        + 12 * A**2 * J**2 * V**2
        - 4 * A * B**3 * J**2 * V
        + 4 * A * B**2 * J * V**2
        + 24 * A * B * J * V**3
        - 48 * A * J * V**4
        - 16 * A * V**4
        - 8 * B**4 * J * V**2
        + 8 * B**3 * J * V**3
        + 32 * B**2 * V**4
        - 96 * B * V**5
        + 64 * V**6
    )


def universal_checks() -> None:
    y, V, A, B, J = sp.symbols("y V A B J")
    L = y**2 + B * y
    R1 = 2 * L * (2 * y + B) + V * A
    R2 = V**3 + V**2 * y + J * y * L
    require(
        sp.expand(sp.resultant(R1, R2, y) - V**3 * K_general(V, A, B, J)) == 0,
        "universal y-resultant",
    )

    x, t = sp.symbols("x t")
    v, v15, v14, a = sp.symbols("v v15 v14 a", nonzero=True)
    a7 = a * v15 / (2 * v)
    a6 = a * (4 * v * v14 - v15**2) / (8 * v**2)
    Vjet = v * x**16 + v15 * x**15 + v14 * x**14
    Ajet = a * x**8 + a7 * x**7 + a6 * x**6
    B9 = 1 + t * x**9
    J9 = 2 * Vjet * sp.diff(B9, x) - B9 * sp.diff(Vjet, x)
    K9 = sp.Poly(sp.expand(K_general(Vjet, Ajet, B9, J9)), x)
    c99 = sp.factor(K9.coeff_monomial(x**99))
    c98 = sp.factor(K9.coeff_monomial(x**98))
    wall = v / a
    kappa = 4 * v * v14 - 3 * v15**2
    c97_wall = sp.factor(K9.coeff_monomial(x**97).subs(t, wall))
    require(c99 == -16 * t**4 * v**3 * (a * t - v), "x^99 infinity coefficient")
    require(c98 == -72 * t**4 * v**2 * v15 * (a * t - v), "x^98 infinity coefficient")
    require(c97_wall == -2 * v**6 * kappa / a**4, "x^97 wall coefficient")

    ledger = []
    for degree in range(10):
        j_degree = 15 if degree == 0 else (degree + 15 if degree < 8 else (22 if degree == 8 else 24))
        bounds = (
            24 + 3 * j_degree,
            48 + 2 * j_degree,
            24 + 3 * degree + 2 * j_degree,
            40 + 2 * degree + j_degree,
            56 + degree + j_degree,
            72 + j_degree,
            72,
            32 + 4 * degree + j_degree,
            48 + 3 * degree + j_degree,
            64 + 2 * degree,
            80 + degree,
            96,
        )
        ledger.append(max(bounds))
    require(tuple(ledger) == (96,) * 9 + (99,), "degree-nine minimality ledger")
    print("universal_resultant=PASS")
    print("infinity_coefficients=x99,x98_share_(a*t-v);x97_wall=-2*v^6*kappa/a^4")
    print("degree_ledger_d=0..9=" + repr(tuple(ledger)))


def family_in_t(V: Poly, A: Poly, boundary: Poly, degree: int) -> tuple[Poly, ...]:
    dV = V.derivative()
    x_power = X_P**degree
    B = ONE_T + T_T * x_power
    J_linear = (V * (X_P ** (degree - 1))).scalar(2 * degree) - x_power * dV
    J = as_param(-dV) + T_T * J_linear
    K = K_general(as_param(V), as_param(A), B, J)
    return tuple(entry.exact_quotient(boundary) for entry in K.c)


def x_coefficient(parameter_coefficients: tuple[Poly, ...], degree: int) -> tuple[Alg, ...]:
    entries = []
    for coefficient in parameter_coefficients:
        entries.append(coefficient.c[degree] if degree < len(coefficient.c) else ZERO_A)
    while entries and entries[-1] == ZERO_A:
        entries.pop()
    return tuple(entries)


def parameter_evaluate(coefficients: tuple[Alg, ...], value: Alg) -> Alg:
    answer = ZERO_A
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def parameter_derivative_evaluate(coefficients: tuple[Alg, ...], value: Alg) -> Alg:
    answer = ZERO_A
    for degree in range(len(coefficients) - 1, 0, -1):
        answer = answer * value + coefficients[degree] * degree
    return answer


def specialized_residual(V: Poly, A: Poly, boundary: Poly, value: Alg) -> tuple[Poly, Poly]:
    B = ONE_P + (X_P**9).scalar(value)
    J = V.scalar(18 * value) * X_P**8 - B * V.derivative()
    return B, K_general(V, A, B, J).exact_quotient(boundary)


def build_case(name: str) -> None:
    if name == "4111":
        configure(100, 244, 237, 44)
        u = U_A
        exponent_a, exponent_b = 4, 1
        accessory_v = (8 * u**2 + 9 * u + 8) / 7
        s0 = 5 * (u + 1) / 7
        A0 = 80 * accessory_v**2 * (u + 1) / 343
        extras = (9, 0)
        expected_norm = Q(494424620106921, 3276800000000)
        expected_wall = (1250000000, 5131065625, 7626867738, 4689453125)
    else:
        configure(75, -89, -31, 61)
        u = U_A
        exponent_a, exponent_b = 3, 2
        accessory_v = (24 * u**2 - 16 * u - 16) / 21
        s0 = (5 * u - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * u - 4) / 343
        extras = (6, 3)
        expected_norm = Q(
            -102215864475620375014841556549072265625,
            12116574790945106558976,
        )
        expected_wall = (23887872, 293981184, -339674272, -14068359375)

    gamma = -7 * A0
    q2 = X_P**2 - X_P.scalar(u) + accessory_v
    D = X_P**exponent_a * (X_P - 1) ** exponent_b * q2
    T = X_P * (X_P - 1) * q2
    E = (
        (X_P - 1) * q2 * exponent_a
        + X_P * q2 * exponent_b
        + X_P * (X_P - 1) * (2 * X_P - u)
    ).scalar(Q(1, 7))
    S = X_P + s0
    V = (S * D * T**2).scalar(4 / gamma**2)
    A = (S * E * T).scalar(2 / gamma)
    boundary = S**3 * T**8 * X_P**extras[0] * (X_P - 1) ** extras[1]
    require(V.degree == 16 and A.degree == 8 and boundary.degree == 44, "degree data")
    require(V.lc == 4 / gamma**2 and A.lc == 2 / gamma, "leading response coefficients")
    require(2 * V * A.derivative() - A * V.derivative() == 2 * V, "response identity")

    H8 = family_in_t(V, A, boundary, 8)
    H9 = family_in_t(V, A, boundary, 9)
    degrees8 = tuple(coefficient.degree for coefficient in H8)
    degrees9 = tuple(coefficient.degree for coefficient in H9)
    require(max(degrees8) == 52, "degree-eight residual stays at 52")
    require(max(degrees9) == 55, "degree-nine residual reaches 55")

    v = V.c[16]
    v15 = V.c[15]
    v14 = V.c[14]
    a = A.c[8]
    kappa = 4 * v * v14 - 3 * v15**2
    require(v == 4 / gamma**2 and a == 2 / gamma, "gamma normalization")
    require(norm(kappa) == expected_norm and expected_norm != 0, "nonzero kappa norm")
    wall_norm = primitive_wall_norm(gamma)
    require(wall_norm == expected_wall, "primitive wall norm")

    h55 = x_coefficient(H9, 55)
    h54 = x_coefficient(H9, 54)
    h53 = x_coefficient(H9, 53)
    expected_h55 = (ZERO_A,) * 4 + (4096 / gamma**8, -2048 / gamma**7)
    require(h55 == expected_h55, "exact h55 parameter polynomial")
    t_wall = 2 / gamma
    require(parameter_evaluate(h55, t_wall) == ZERO_A, "h55 wall cancellation")
    require(parameter_evaluate(h54, t_wall) == ZERO_A, "h54 wall cancellation")
    require(
        parameter_evaluate(h53, t_wall) == -512 * kappa / gamma**8,
        "h53 wall carry",
    )
    require(
        parameter_derivative_evaluate(h55, t_wall) == -32768 / gamma**11,
        "transverse wall derivative",
    )

    B_zero, H_zero = specialized_residual(V, A, boundary, ZERO_A)
    B_generic, H_generic = specialized_residual(V, A, boundary, as_alg(2))
    B_wall, H_wall = specialized_residual(V, A, boundary, t_wall)
    require(B_zero == ONE_P and H_zero.degree == 52, "constant support control")
    require(B_generic.degree == 9 and H_generic.degree == 55, "generic nonlinear control")
    require(B_wall.degree == 9 and H_wall.degree == 53, "wall degree drop by two")

    reciprocal_delta_coefficient = -64 / (gamma**3 * kappa)
    reciprocal_epsilon_coefficient = -64 / (gamma**4 * kappa)
    require(
        reciprocal_epsilon_coefficient * gamma == reciprocal_delta_coefficient,
        "epsilon and delta reciprocal laws",
    )

    print(f"case={name}")
    print(f"degree8_parameter_quotients={degrees8}")
    print(f"degree9_parameter_quotients={degrees9}")
    print(f"wall_norm_coefficients={wall_norm}")
    print(f"kappa_norm={expected_norm.numerator}/{expected_norm.denominator}")
    print("degrees=(t=0:52,t=2:55,t=2/Gamma:53)")
    print(f"generic_H_sha256={sha256(serialize(H_generic)).hexdigest()}")
    print(f"wall_H_sha256={sha256(serialize(H_wall)).hexdigest()}")
    print("reciprocal_law=w^2~-64*delta/(Gamma^3*kappa)=-64*epsilon/(Gamma^4*kappa)")


def finite_case(name: str) -> None:
    x = sp.symbols("x")
    if name == "4111":
        p, u, exponent_a, exponent_b = 113, 85, 4, 1
        q_coefficients = (100, 244, 237, 44)
        accessory_v = (8 * u**2 + 9 * u + 8) * pow(7, -1, p) % p
        s0 = 5 * (u + 1) * pow(7, -1, p) % p
        A0 = 80 * accessory_v**2 * (u + 1) * pow(343, -1, p) % p
        extras = (9, 0)
        expected_wall = 85
    else:
        p, u, exponent_a, exponent_b = 101, 64, 3, 2
        q_coefficients = (75, -89, -31, 61)
        accessory_v = (24 * u**2 - 16 * u - 16) * pow(21, -1, p) % p
        s0 = (5 * u - 4) * pow(7, -1, p) % p
        A0 = 9 * accessory_v**2 * (5 * u - 4) * pow(343, -1, p) % p
        extras = (6, 3)
        expected_wall = 89

    c3, c2, c1, c0 = q_coefficients
    require((c3 * u**3 + c2 * u**2 + c1 * u + c0) % p == 0, "accessory reduction root")
    require((3 * c3 * u**2 + 2 * c2 * u + c1) % p != 0, "simple accessory reduction")
    gamma = -7 * A0 % p
    require(gamma != 0, "gamma reduction is a unit")
    P = lambda expression: sp.Poly(expression, x, modulus=p)
    q2 = P(x**2 - u * x + accessory_v)
    D = P(x**exponent_a * (x - 1) ** exponent_b) * q2
    T = P(x * (x - 1)) * q2
    E = P(
        (
            exponent_a * (x - 1) * (x**2 - u * x + accessory_v)
            + exponent_b * x * (x**2 - u * x + accessory_v)
            + x * (x - 1) * (2 * x - u)
        )
        * pow(7, -1, p)
    )
    S = P(x + s0)
    V = (S * D * T**2).mul_ground(4 * pow(gamma**2, -1, p) % p)
    A = (S * E * T).mul_ground(2 * pow(gamma, -1, p) % p)
    boundary = S**3 * T**8 * P(x) ** extras[0] * P(x - 1) ** extras[1]
    g = S * T

    def residual(parameter: int):
        B = P(1 + parameter * x**9)
        Bprime = P(9 * parameter * x**8)
        J = V.mul_ground(2) * Bprime - B * V.diff()
        K = K_general(V, A, B, J)
        return B, sp.exquo(K, boundary)

    wall = 2 * pow(gamma, -1, p) % p
    require(wall == expected_wall, "wall reduction")
    B_zero, H_zero = residual(0)
    B_generic, H_generic = residual(2)
    B_wall, H_wall = residual(wall)
    require(H_zero.degree() == 52, "finite t=0 degree")
    require(H_generic.degree() == 55, "finite generic degree")
    require(H_wall.degree() == 53, "finite wall degree")

    controls = []
    for label, B, H in (
        ("generic", B_generic, H_generic),
        ("wall", B_wall, H_wall),
    ):
        gcd_B_g = sp.gcd(B, g).degree()
        gcd_H_g = sp.gcd(H, g).degree()
        gcd_H_derivative = sp.gcd(H, H.diff()).degree()
        require((gcd_B_g, gcd_H_g, gcd_H_derivative) == (0, 0, 0), f"finite {label} gcd gates")
        controls.append((label, H.degree(), gcd_B_g, gcd_H_g, gcd_H_derivative))

    v = int(V.nth(16)) % p
    v15 = int(V.nth(15)) % p
    v14 = int(V.nth(14)) % p
    kappa = (4 * v * v14 - 3 * v15 * v15) % p
    require(kappa != 0, "finite kappa carry")
    print(
        f"finite_case={name} good_reduction=(p={p},u={u}) "
        f"wall_t={wall} kappa={kappa} controls={tuple(controls)}"
    )


def main() -> None:
    universal_checks()
    build_case("4111")
    build_case("3211")
    finite_case("4111")
    finite_case("3211")
    print("SCOPE=PROVED_CANDIDATE_SYMBOLIC_PLUS_FINITE_EXACT_CONTROLS;JC2=NOT_CLAIMED")
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
