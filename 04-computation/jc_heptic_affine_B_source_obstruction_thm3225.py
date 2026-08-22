#!/usr/bin/env python3
"""Exact candidate audit for the affine-B extension of THM-3212.

The family is

    P_B(x,z) = (V(x) z^2 + B(x) z)^2 + A(x) z + x,
    deg(B) <= 1.

All characteristic-zero accessory arithmetic is performed in Q[u]/(q)[x].
Finite-field reductions are used only as independent squarefreeness controls.
No repository executable is imported and no assert statement is used.
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


@dataclass(frozen=True)
class ParamPoly:
    """Polynomial in the affine parameter t with Poly coefficients."""

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


def K_general(V, A, B, J):
    """The V^(-3)-saturated y-resultant for c=0,k=1."""
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


def universal_symbolic_check() -> None:
    y, V, A, B, J = sp.symbols("y V A B J")
    L = y**2 + B * y
    R1 = 2 * L * (2 * y + B) + V * A
    R2 = V**3 + V**2 * y + J * y * L
    expected = K_general(V, A, B, J)
    derived = sp.factor(sp.resultant(R1, R2, y))
    require(sp.expand(derived - V**3 * expected) == 0, "universal resultant")

    Vp, Bp, Ap = sp.symbols("Vp Bp Ap")
    raw = 2 * L * (Vp * y**2 + Bp * V * y) + Ap * y * V**2 + V**3
    covariant = 2 * V * Bp - B * Vp
    response_Ap = (A * Vp + 2 * V) / (2 * V)
    reduction = sp.cancel(raw.subs(Ap, response_Ap) - R2.subs(J, covariant) - Vp * y * R1 / 2)
    require(sp.expand(reduction) == 0, "gradient-ideal reduction")

    w, b, e, v = sp.symbols("w b e v", nonzero=True)
    for m in (3, 4, 5, 6):
        local_V = v * w**m
        local_A = sp.Rational(2, 2 - m) * w
        local_B = b + e * w
        local_J = 2 * local_V * e - local_B * sp.diff(local_V, w)
        local_K = sp.expand(K_general(local_V, local_A, local_B, local_J))
        coefficient = local_K.coeff(w, 3 * m - 1)
        target = sp.Rational(16 * m * (m - 1), m - 2) * b**5 * v**3
        require(sp.factor(coefficient - target) == 0, f"T-root leading coefficient m={m}")

    v1, v2, v3, v4, v5 = sp.symbols("v1:6", nonzero=True)
    a2, a3, a4, a5 = sp.symbols("a2:6")
    local_V = v1 * w + v2 * w**2 + v3 * w**3 + v4 * w**4 + v5 * w**5
    local_A = 2 * w + a2 * w**2 + a3 * w**3 + a4 * w**4 + a5 * w**5
    response = sp.expand(2 * local_V * sp.diff(local_A, w) - local_A * sp.diff(local_V, w) - 2 * local_V)
    solved = {}
    for degree, variable in zip(range(2, 6), (a2, a3, a4, a5)):
        equation = sp.expand(response.subs(solved)).coeff(w, degree)
        solved[variable] = sp.solve(equation, variable)[0]
    local_A = sp.expand(local_A.subs(solved))
    local_B = b + e * w
    local_J = 2 * local_V * e - local_B * sp.diff(local_V, w)
    local_K = sp.expand(K_general(local_V, local_A, local_B, local_J))
    require(all(local_K.coeff(w, degree) == 0 for degree in range(3)), "S boundary order at least three")
    wall = 4 * b * v2 - 6 * e * v1 + 3 * v1**2
    target3 = -sp.Rational(8, 3) * b**4 * v1**2 * wall
    require(sp.factor(local_K.coeff(w, 3) - target3) == 0, "S-wall coefficient")

    escape_e = (4 * b * v2 + 3 * v1**2) / (6 * v1)
    Q4 = b**3 * (-54 * v1 * v3 + 8 * v2**2) - 30 * b**2 * v1**2 * v2 + 45 * v1**3
    Q5 = (
        b**5 * (3240 * v1**2 * v4 + 11016 * v1 * v2 * v3 - 1552 * v2**3)
        + b**4 * (3969 * v1**3 * v3 + 6132 * v1**2 * v2**2)
        + b**3 * (1260 * v1**4 * v2 - 4536 * v1**2 * v3 + 672 * v1 * v2**2)
        - 12600 * b**2 * v1**3 * v2
        + 945 * b * v1**5
        + 3780 * v1**4
    )
    coefficient4 = sp.factor(sp.cancel(local_K.coeff(w, 4).subs(e, escape_e)))
    coefficient5 = sp.factor(sp.cancel(local_K.coeff(w, 5).subs(e, escape_e)))
    require(sp.factor(coefficient4 - sp.Rational(16, 45) * b**2 * v1 * Q4) == 0, "S escape Q4")
    require(sp.factor(coefficient5 + sp.Rational(8, 945) * Q5) == 0, "S escape Q5")
    print("universal_resultant_and_local_jets=PASS")


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
    else:
        configure(75, -89, -31, 61)
        u = U_A
        a, b = 3, 2
        v = (24 * u**2 - 16 * u - 16) / 21
        s0 = (5 * u - 4) / 7
        A0 = 9 * v**2 * (5 * u - 4) / 343
        extras = (6, 3)
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

    # Universal B=1+t*x coefficient-by-coefficient boundary division.
    Vt, At = as_param(V), as_param(A)
    Bt = ONE_T + T_T * X_P
    Jt = as_param(-dV) + T_T * (2 * V - X_P * dV)
    Kt = K_general(Vt, At, Bt, Jt)
    H_coefficients = tuple(entry.exact_quotient(boundary) for entry in Kt.c)
    quotient_degrees = tuple(entry.degree for entry in H_coefficients)
    require(len(H_coefficients) == 6, "affine parameter degree five")
    require(quotient_degrees[0] == 52, "constant coefficient residual degree")
    require(max(quotient_degrees[1:]) < 52, "uniform degree-52 leading term")

    derivatives = []
    derivative = V
    for order in range(1, 5):
        derivative = derivative.derivative()
        derivatives.append(evaluate(derivative, -s0) / math.factorial(order))
    v1, v2, v3, v4 = derivatives
    require(v1 != ZERO_A and v2 != ZERO_A, "nonzero S jets")

    Q4 = (
        (X_P**3).scalar(-54 * v1 * v3 + 8 * v2**2)
        + (X_P**2).scalar(-30 * v1**2 * v2)
        + ONE_P.scalar(45 * v1**3)
    )
    Q5 = (
        (X_P**5).scalar(3240 * v1**2 * v4 + 11016 * v1 * v2 * v3 - 1552 * v2**3)
        + (X_P**4).scalar(3969 * v1**3 * v3 + 6132 * v1**2 * v2**2)
        + (X_P**3).scalar(1260 * v1**4 * v2 - 4536 * v1**2 * v3 + 672 * v1 * v2**2)
        + (X_P**2).scalar(-12600 * v1**3 * v2)
        + X_P.scalar(945 * v1**5)
        + ONE_P.scalar(3780 * v1**4)
    )
    degrees, common = prs_degrees(Q5, Q4)
    require(degrees == ((5, 3), (3, 2), (2, 1), (1, 0)), "S-escape PRS profile")
    require(common == ONE_P, "S-escape Q4/Q5 gcd")

    # Structured S-escape hostile inside B=1+t*x.
    denominator = 4 * s0 * v2 + 6 * v1
    require(denominator != ZERO_A, "escape denominator")
    t_escape = (4 * v2 + 3 * v1**2) / denominator
    B_escape = X_P.scalar(t_escape) + 1
    J_escape = V.scalar(2 * t_escape) - B_escape * dV
    H_escape = K_general(V, A, B_escape, J_escape).exact_quotient(boundary)
    require(H_escape.degree == 52, "escape residual degree")
    require(evaluate(B_escape, -s0) != ZERO_A, "escape B nonzero at S")
    require(evaluate(g, -1 / t_escape) != ZERO_A, "escape B unit modulo g")
    quotient, remainder = H_escape.divmod(S)
    require(remainder.is_zero(), "escape root at S")
    _, second_remainder = quotient.divmod(S)
    require(not second_remainder.is_zero(), "escape multiplicity exactly one")

    # Exact nonconstant control B=1+x; finite reductions certify Morse-ness.
    B_control = 1 + X_P
    J_control = 2 * V - B_control * dV
    H_control = K_general(V, A, B_control, J_control).exact_quotient(boundary)
    require(H_control.degree == 52, "nonconstant control degree")
    require(evaluate(g, -1) != ZERO_A, "nonconstant control B unit modulo g")

    print(f"case={name}")
    print(f"parameter_quotient_degrees={quotient_degrees}")
    print(f"S_escape_PRS_degrees={degrees} gcd=1")
    print(f"Q4_sha256={hashlib.sha256(serialize(Q4)).hexdigest()}")
    print(f"Q5_sha256={hashlib.sha256(serialize(Q5)).hexdigest()}")
    print(f"t_escape_coefficients={t_escape.c}")
    print(f"escape_H_degree=52 escape_S_multiplicity=1 escape_H_sha256={hashlib.sha256(serialize(H_escape)).hexdigest()}")
    print(f"control_B=1+x H_degree=52 H_sha256={hashlib.sha256(serialize(H_control)).hexdigest()}")


def finite_case(name: str) -> None:
    x = sp.symbols("x")
    if name == "4111":
        p, u, a, b = 113, 85, 4, 1
        v = (8 * u**2 + 9 * u + 8) * pow(7, -1, p) % p
        s0 = 5 * (u + 1) * pow(7, -1, p) % p
        A0 = 80 * v**2 * (u + 1) * pow(343, -1, p) % p
        extras = (9, 0)
    else:
        p, u, a, b = 101, 64, 3, 2
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

    def finite_K(B, J):
        return (
            -A**3 * J**3
            + (A**2 * J**2 * V**2).mul_ground(12)
            - (A * B**3 * J**2 * V).mul_ground(4)
            + (A * B**2 * J * V**2).mul_ground(4)
            + (A * B * J * V**3).mul_ground(24)
            - (A * J * V**4).mul_ground(48)
            - (A * V**4).mul_ground(16)
            - (B**4 * J * V**2).mul_ground(8)
            + (B**3 * J * V**3).mul_ground(8)
            + (B**2 * V**4).mul_ground(32)
            - (B * V**5).mul_ground(96)
            + (V**6).mul_ground(64)
        )

    B_control = P(1 + x)
    J_control = V.mul_ground(2) - B_control * dV
    H_control = sp.exquo(finite_K(B_control, J_control), boundary)
    require(H_control.degree() == 52, "finite control degree")
    require(sp.gcd(B_control, g).degree() == 0, "finite control B unit")
    require(sp.gcd(H_control, g).degree() == 0, "finite control boundary disjoint")
    require(sp.gcd(H_control, H_control.diff()).degree() == 0, "finite control squarefree")

    v1 = int(dV.eval(-s0)) % p
    v2 = int(dV.diff().eval(-s0)) * pow(2, -1, p) % p
    t_escape = (4 * v2 + 3 * v1 * v1) * pow((4 * s0 * v2 + 6 * v1) % p, -1, p) % p
    B_escape = P(1 + t_escape * x)
    J_escape = V.mul_ground(2 * t_escape) - B_escape * dV
    H_escape = sp.exquo(finite_K(B_escape, J_escape), boundary)
    require(H_escape.degree() == 52, "finite escape degree")
    require(sp.gcd(B_escape, g).degree() == 0, "finite escape B unit")
    require(sp.gcd(H_escape, g).degree() == 1, "finite escape one S root")
    require(sp.gcd(H_escape, H_escape.diff()).degree() == 0, "finite escape squarefree")

    B_collision = P(1 - x)
    require(sp.gcd(B_collision, g).degree() == 1, "finite collision hostile")
    print(
        f"finite_case={name} good_reduction=(p={p},u={u}) "
        f"control=(degree=52,gcd_g=0,gcd_derivative=0) "
        f"escape=(t={t_escape},degree=52,gcd_g=1,gcd_derivative=0) "
        f"collision_B=1-x_gcd_g=1"
    )


def main() -> None:
    universal_symbolic_check()
    build_case("4111")
    build_case("3211")
    finite_case("4111")
    finite_case("3211")
    print("FAILED_CHECKS=NONE")


if __name__ == "__main__":
    main()
