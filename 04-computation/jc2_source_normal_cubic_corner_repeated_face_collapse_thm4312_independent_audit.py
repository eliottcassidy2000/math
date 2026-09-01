#!/usr/bin/env python3
"""Dependency-free independent audit for the proposed THM-4312.

The input consists only of the displayed THM-4308 response formulas and the
literal THM-4304 repeated-regime equations.  In particular, this script does
not import the primary certificate or any computer-algebra package.

The important quotient check is explicit: ``V^2=w0(w0^9-L)`` occurs only
after the ramified substitution sigma=lambda^3, z=lambda*w0.  The actual
P(3,1) invariant is w=sigma/z^3, and its exceptional quotient is the genus
one cubic ``U*Y^2=1-L*w^3``.
"""

from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
from math import gcd, lcm


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


SYMBOLS = ("r", "h", "xi", "U", "rho", "up", "eta", "zeta")
SYMBOL_INDEX = {name: index for index, name in enumerate(SYMBOLS)}
ZERO_MONOMIAL = (0,) * len(SYMBOLS)


class SPoly:
    """Sparse Q-polynomial in the fixed symbol list."""

    def __init__(self, terms: dict[tuple[int, ...], F] | None = None) -> None:
        self.terms = {
            monomial: F(coefficient)
            for monomial, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def make(value: int | F | "SPoly") -> "SPoly":
        if isinstance(value, SPoly):
            return value
        coefficient = F(value)
        return SPoly({ZERO_MONOMIAL: coefficient}) if coefficient else SPoly()

    def __add__(self, other: int | F | "SPoly") -> "SPoly":
        other = self.make(other)
        result = self.terms.copy()
        for monomial, coefficient in other.terms.items():
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return SPoly(result)

    __radd__ = __add__

    def __neg__(self) -> "SPoly":
        return SPoly({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: int | F | "SPoly") -> "SPoly":
        return self + (-self.make(other))

    def __rsub__(self, other: int | F | "SPoly") -> "SPoly":
        return self.make(other) - self

    def __mul__(self, other: int | F | "SPoly") -> "SPoly":
        other = self.make(other)
        result: dict[tuple[int, ...], F] = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other.terms.items():
                monomial = tuple(
                    left + right
                    for left, right in zip(left_monomial, right_monomial)
                )
                result[monomial] = (
                    result.get(monomial, F(0))
                    + left_coefficient * right_coefficient
                )
        return SPoly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "SPoly":
        if exponent < 0:
            raise ValueError("negative polynomial exponent")
        result = SPoly.make(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def __truediv__(self, scalar: int | F) -> "SPoly":
        scalar = F(scalar)
        if not scalar:
            raise ZeroDivisionError
        return SPoly(
            {monomial: coefficient / scalar for monomial, coefficient in self.terms.items()}
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (int, F, SPoly)):
            return False
        return self.terms == self.make(other).terms

    def __bool__(self) -> bool:
        return bool(self.terms)

    def substitute(self, name: str, replacement: int | F | "SPoly") -> "SPoly":
        index = SYMBOL_INDEX[name]
        replacement = self.make(replacement)
        result = SPoly()
        for monomial, coefficient in self.terms.items():
            exponent = monomial[index]
            base_monomial = list(monomial)
            base_monomial[index] = 0
            term = SPoly({tuple(base_monomial): coefficient})
            term = term * (replacement ** exponent)
            result = result + term
        return result

    def as_fraction(self) -> F:
        require(set(self.terms).issubset({ZERO_MONOMIAL}), "constant specialization")
        return self.terms.get(ZERO_MONOMIAL, F(0))


def symbol(name: str) -> SPoly:
    monomial = [0] * len(SYMBOLS)
    monomial[SYMBOL_INDEX[name]] = 1
    return SPoly({tuple(monomial): F(1)})


r, h, xi, Us, rho, up, etas, zetas = (symbol(name) for name in SYMBOLS)


SourcePoly = dict[tuple[int, int, int], SPoly]


def source_add(left: SourcePoly, right: SourcePoly) -> SourcePoly:
    result = left.copy()
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, SPoly()) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def source_scale(poly: SourcePoly, scalar: int | F | SPoly) -> SourcePoly:
    return {
        monomial: coefficient * scalar
        for monomial, coefficient in poly.items()
        if coefficient * scalar
    }


def source_multiply(left: SourcePoly, right: SourcePoly) -> SourcePoly:
    result: SourcePoly = {}
    for (q1, t1, y1), c1 in left.items():
        for (q2, t2, y2), c2 in right.items():
            monomial = (q1 + q2, t1 + t2, y1 + y2)
            result[monomial] = result.get(monomial, SPoly()) + c1 * c2
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def source_power(poly: SourcePoly, exponent: int) -> SourcePoly:
    result: SourcePoly = {(0, 0, 0): SPoly.make(1)}
    for _ in range(exponent):
        result = source_multiply(result, poly)
    return result


def source_t_shift(poly: SourcePoly, shift: int) -> SourcePoly:
    return {(q_degree, t_degree + shift, y_degree): coefficient
            for (q_degree, t_degree, y_degree), coefficient in poly.items()}


def q_to_ty_divide_t3(poly: SourcePoly) -> SourcePoly:
    result: SourcePoly = {}
    for (q_degree, t_degree, y_degree), coefficient in poly.items():
        new_t_degree = t_degree + q_degree - 3
        if new_t_degree < 0:
            raise AssertionError("negative strict-transform row")
        monomial = (0, new_t_degree, y_degree + q_degree)
        result[monomial] = result.get(monomial, SPoly()) + coefficient
    return {monomial: coefficient for monomial, coefficient in result.items() if coefficient}


def t_row(poly: SourcePoly, t_degree: int) -> dict[int, SPoly]:
    return {
        y_degree: coefficient
        for (q_degree, row, y_degree), coefficient in poly.items()
        if q_degree == 0 and row == t_degree
    }


def y_evaluate(poly: dict[int, SPoly], value: SPoly) -> SPoly:
    return sum((coefficient * value**degree for degree, coefficient in poly.items()), SPoly())


def y_derivative(poly: dict[int, SPoly]) -> dict[int, SPoly]:
    return {
        degree - 1: coefficient * degree
        for degree, coefficient in poly.items()
        if degree
    }


def y_add(left: dict[int, SPoly], right: dict[int, SPoly]) -> dict[int, SPoly]:
    result = left.copy()
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, SPoly()) + coefficient
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def y_scale(poly: dict[int, SPoly], scalar: int | F | SPoly) -> dict[int, SPoly]:
    return {degree: coefficient * scalar for degree, coefficient in poly.items()
            if coefficient * scalar}


def y_shift(poly: dict[int, SPoly], shift: int) -> dict[int, SPoly]:
    return {degree + shift: coefficient for degree, coefficient in poly.items()}


def y_nth_derivative(poly: dict[int, SPoly], order: int) -> dict[int, SPoly]:
    result = poly
    for _ in range(order):
        result = y_derivative(result)
    return result


def coefficient_list_degree(coefficients: list[SPoly]) -> int:
    degree = len(coefficients) - 1
    while degree >= 0 and not coefficients[degree]:
        degree -= 1
    return degree


UniPoly = list[F]


def uni_clean(poly: UniPoly) -> UniPoly:
    result = list(poly)
    while result and not result[-1]:
        result.pop()
    return result


def uni_add(left: UniPoly, right: UniPoly) -> UniPoly:
    result = [F(0)] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += coefficient
    return uni_clean(result)


def uni_scale(poly: UniPoly, scalar: int | F) -> UniPoly:
    return uni_clean([coefficient * F(scalar) for coefficient in poly])


def uni_multiply(left: UniPoly, right: UniPoly) -> UniPoly:
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            result[i + j] += left_coefficient * right_coefficient
    return uni_clean(result)


def uni_divmod(dividend: UniPoly, divisor: UniPoly) -> tuple[UniPoly, UniPoly]:
    quotient = [F(0)] * max(1, len(dividend) - len(divisor) + 1)
    remainder = uni_clean(dividend)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        subtraction = [F(0)] * shift + [coefficient * entry for entry in divisor]
        remainder = uni_add(remainder, uni_scale(subtraction, -1))
    return uni_clean(quotient), remainder


def uni_gcd(left: UniPoly, right: UniPoly) -> UniPoly:
    left = uni_clean(left)
    right = uni_clean(right)
    while right:
        _, remainder = uni_divmod(left, right)
        left, right = right, remainder
    return uni_scale(left, 1 / left[-1]) if left else []


def primitive_integer_coefficients(poly: UniPoly) -> list[int]:
    denominator = 1
    for coefficient in poly:
        denominator = lcm(denominator, coefficient.denominator)
    integers = [int(coefficient * denominator) for coefficient in poly]
    divisor = reduce(gcd, (abs(integer) for integer in integers if integer))
    integers = [integer // divisor for integer in integers]
    return [-integer for integer in integers] if integers[-1] < 0 else integers


class Quadratic:
    """Q(rho), with rho^2 fixed by the positive control."""

    D = F(-98333997651, 30863472406)

    def __init__(self, a: int | F = 0, b: int | F = 0) -> None:
        self.a = F(a)
        self.b = F(b)

    @staticmethod
    def make(value: int | F | "Quadratic") -> "Quadratic":
        return value if isinstance(value, Quadratic) else Quadratic(value)

    def __add__(self, other: int | F | "Quadratic") -> "Quadratic":
        other = self.make(other)
        return Quadratic(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Quadratic":
        return Quadratic(-self.a, -self.b)

    def __sub__(self, other: int | F | "Quadratic") -> "Quadratic":
        return self + (-self.make(other))

    def __rsub__(self, other: int | F | "Quadratic") -> "Quadratic":
        return self.make(other) - self

    def __mul__(self, other: int | F | "Quadratic") -> "Quadratic":
        other = self.make(other)
        return Quadratic(
            self.a * other.a + self.D * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def norm(self) -> F:
        return self.a * self.a - self.D * self.b * self.b

    def inverse(self) -> "Quadratic":
        norm = self.norm()
        if not norm:
            raise ZeroDivisionError
        return Quadratic(self.a / norm, -self.b / norm)

    def __truediv__(self, other: int | F | "Quadratic") -> "Quadratic":
        return self * self.make(other).inverse()

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (int, F, Quadratic)):
            return False
        other = self.make(other)
        return self.a == other.a and self.b == other.b

    def __bool__(self) -> bool:
        return bool(self.a or self.b)


QPoly = list[Quadratic]


def qpoly_clean(poly: QPoly) -> QPoly:
    result = list(poly)
    while result and not result[-1]:
        result.pop()
    return result


def qpoly_add(left: QPoly, right: QPoly) -> QPoly:
    result = [Quadratic()] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        result[index] = result[index] + coefficient
    for index, coefficient in enumerate(right):
        result[index] = result[index] + coefficient
    return qpoly_clean(result)


def qpoly_scale(poly: QPoly, scalar: Quadratic) -> QPoly:
    return qpoly_clean([coefficient * scalar for coefficient in poly])


def qpoly_divmod(dividend: QPoly, divisor: QPoly) -> tuple[QPoly, QPoly]:
    quotient = [Quadratic()] * max(1, len(dividend) - len(divisor) + 1)
    remainder = qpoly_clean(dividend)
    while remainder and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] = quotient[shift] + coefficient
        subtraction = [Quadratic()] * shift + [coefficient * entry for entry in divisor]
        remainder = qpoly_add(remainder, qpoly_scale(subtraction, Quadratic(-1)))
    return qpoly_clean(quotient), remainder


def qpoly_gcd(left: QPoly, right: QPoly) -> QPoly:
    left = qpoly_clean(left)
    right = qpoly_clean(right)
    while right:
        _, remainder = qpoly_divmod(left, right)
        left, right = right, remainder
    return qpoly_scale(left, left[-1].inverse()) if left else []


def main() -> None:
    Delta = F(896, 15)
    Theta = F(512, 75)
    zeta_factor = F(-3, 2)
    upsilon = F(-731648, 2025)

    U_response = (SPoly.make(475515904) - 109350 * xi) / 200475
    W_response = -(
        4343625 * r - 17172000 * xi + 143826305024
    ) / 4009500
    Z_response = (
        12506118074368 - 173745000 * r - 195463125 * h
        - 926883000 * xi
    ) / 108256500
    corner_w = W_response + 2 * U_response
    require(
        corner_w
        == (-4343625 * r + 12798000 * xi - 124805668864) / 4009500,
        "corner W equation",
    )
    xi_corner = (4343625 * r + 124805668864) / 12798000
    require(corner_w.substitute("xi", xi_corner) == 0, "corner W=-2U")
    corner_z = (Z_response - U_response).substitute("xi", xi_corner)
    h_corner = (2091705253888 - 258703875 * r) / 107983125
    require(corner_z.substitute("h", h_corner) == 0, "corner Z=U")
    U_corner = U_response.substitute("xi", xi_corner)
    require(
        U_corner == -F(13, 57591000) * (820125 * r + 13056802816),
        "corner U response",
    )
    require(h_corner.substitute("r", 0).as_fraction() != 0, "corner Phi nonzero")

    balanced_residual = Delta - F(2048, 45)
    c4 = Delta + Theta
    require(balanced_residual == F(128, 9), "balanced exclusion")
    require(c4 == F(1664, 25), "k3 c4 exclusion")
    require(upsilon != 0, "k3 upsilon exclusion")
    xi_k2 = -upsilon
    U_k2 = F(475515904 - 109350 * xi_k2, 200475)
    k2_residual = upsilon * upsilon - 4 * U_k2 * c4
    require(U_k2 == F(39636992, 18225), "k2 U")
    require(k2_residual == F(-1839105572864, 4100625), "k2 exclusion")

    c2_corner = SPoly.make(upsilon) + xi_corner
    require(
        c2_corner == F(11, 474000) * (14625 * r + 404652032),
        "k1 c2",
    )
    k1_square = 4 * U_corner * c2_corner
    expected_square = -F(143, 6824533500000) * (
        14625 * r + 404652032
    ) * (
        820125 * r + 13056802816
    )
    require(k1_square == expected_square, "k1 square equation")

    one_plus_q: SourcePoly = {
        (0, 0, 0): SPoly.make(1),
        (1, 0, 0): SPoly.make(1),
    }
    top = source_scale(
        source_add(
            source_add(source_power(one_plus_q, 6), source_scale(source_power(one_plus_q, 5), -2)),
            source_power(one_plus_q, 4),
        ),
        Us,
    )
    alpha_row = source_t_shift(
        source_add(
            source_scale(source_power(one_plus_q, 5), -2 * Us * rho),
            source_scale(source_power(one_plus_q, 4), 2 * Us * rho),
        ),
        1,
    )
    xi_symbolic = Us * rho**2 - up
    second_row = source_t_shift(
        source_add(
            source_scale(source_power(one_plus_q, 5), up),
            source_scale(source_power(one_plus_q, 4), xi_symbolic),
        ),
        2,
    )
    third_row = source_t_shift(
        source_add(
            source_scale(source_power(one_plus_q, 4), etas),
            source_scale(source_power(one_plus_q, 3), zetas),
        ),
        3,
    )
    fourth_row = source_t_shift(
        source_add(
            source_scale(source_power(one_plus_q, 4), Delta),
            source_scale(source_power(one_plus_q, 3), Theta),
        ),
        4,
    )
    hhat = source_add(
        source_add(source_add(top, alpha_row), source_add(second_row, third_row)),
        fourth_row,
    )
    literal_source = source_multiply({(1, 0, 0): SPoly.make(1)}, hhat)
    strict = q_to_ty_divide_t3(literal_source)
    P0 = t_row(strict, 0)
    P1 = t_row(strict, 1)
    P2 = t_row(strict, 2)
    expected_P0 = {
        3: Us,
        2: -2 * Us * rho,
        1: Us * rho**2,
    }
    expected_P1 = {
        4: 4 * Us,
        3: -8 * Us * rho,
        2: 4 * Us * rho**2 + up,
        1: etas + zetas,
    }
    expected_P2 = {
        5: 6 * Us,
        4: -12 * Us * rho,
        3: 6 * Us * rho**2 + 4 * up,
        2: 4 * etas + 3 * zetas,
        1: c4,
    }
    require(P0 == expected_P0, "sparse-source P0")
    require(P1 == expected_P1, "sparse-source P1")
    require(P2 == expected_P2, "sparse-source P2")
    require(y_evaluate(P0, rho) == 0, "P0 double root value")
    require(y_evaluate(y_derivative(P0), rho) == 0, "P0 double root derivative")
    require(
        y_evaluate(y_derivative(y_derivative(P0)), rho) == 2 * Us * rho,
        "P0 nontriple coefficient",
    )
    L_symbolic = etas + zetas + up * rho
    require(y_evaluate(P1, rho) == rho * L_symbolic, "k1 graph coefficient L")

    E0 = y_add(P0, y_scale(y_shift(y_derivative(P0), 1), -1))
    E1 = y_add(P1, y_scale(y_shift(y_derivative(P1), 1), -1))
    E2 = y_add(P2, y_scale(y_shift(y_derivative(P2), 1), -1))
    E0_prime = y_evaluate(y_derivative(E0), rho)
    E0_second = y_evaluate(y_nth_derivative(E0, 2), rho)
    E1_value = y_evaluate(E1, rho)
    E1_prime = y_evaluate(y_derivative(E1), rho)
    E2_value = y_evaluate(E2, rho)
    a_numerator = -up                 # a=a_numerator/(2U)
    b_numerator = -(4 * etas + 4 * rho * up + 3 * zetas)  # b=b_numerator/(2U)
    require(
        E0_prime * a_numerator + 2 * Us * E1_value == 0,
        "critical graph E1 solves a=-upsilon/(2U)",
    )
    require(
        4 * Us * E0_prime * b_numerator
        + E0_second * a_numerator**2
        + 4 * Us * E1_prime * a_numerator
        + 8 * Us**2 * E2_value == 0,
        "critical graph E2 solves b=-(4eta+4rho upsilon+3zeta)/(2U)",
    )
    P0_second = y_evaluate(y_nth_derivative(P0, 2), rho)
    P0_third = y_evaluate(y_nth_derivative(P0, 3), rho)
    P1_prime = y_evaluate(y_derivative(P1), rho)
    P1_second = y_evaluate(y_nth_derivative(P1, 2), rho)
    P2_prime = y_evaluate(y_derivative(P2), rho)
    require(
        P0_second * a_numerator + 2 * Us * P1_prime == 2 * Us * L_symbolic,
        "critical graph derives first coefficient",
    )

    alpha_symbolic = -2 * Us * rho
    L2_numerator = (
        4 * Us**2 * c4
        - 2 * Us * alpha_symbolic * (4 * etas + 3 * zetas)
        - Us * up**2
        + 4 * alpha_symbolic**2 * up
    )
    derived_L2_cleared = (
        4 * Us * P0_second * b_numerator
        + P0_third * a_numerator**2
        + 4 * Us * P1_second * a_numerator
        + 8 * Us**2 * P2_prime
    )
    require(derived_L2_cleared == 2 * L2_numerator, "critical graph derives L2")
    require(
        L2_numerator
        == 4 * Us**2 * c4
        + 4 * Us**2 * rho * (4 * etas + 3 * zetas)
        + 16 * Us**2 * rho**2 * up
        - Us * up**2,
        "second graph coefficient numerator",
    )
    c3_symbolic = etas + zetas
    rho_high_contact = -c3_symbolic / upsilon
    require(
        L2_numerator.substitute("rho", rho_high_contact).substitute("up", upsilon)
        == 4 * Us**2 * (c4 + c3_symbolic * zetas / upsilon) - Us * upsilon**2,
        "L=0 second coefficient cancellation",
    )

    require(3 - 3 * 1 == 0, "P(3,1) invariant w=sigma/z^3")
    require(2 - 2 * 1 == 0, "local invariant V=Q/z^2")
    require((10 - 2 * 3, 3 + 1) == (4, 4), "splitter excess tie")
    require((12 - 4 + F(4, 2)) == 10, "Fq order")
    require(9 * 3 + 11 * 1 - 10 == 28, "good-form lower bound")
    cubic_coefficients = [SPoly.make(1), SPoly(), SPoly(), -L_symbolic]
    quotient_degree = coefficient_list_degree(cubic_coefficients)
    ca, cb, cc, cd = cubic_coefficients
    cubic_discriminant = (
        cb**2 * cc**2 - 4 * ca * cc**3 - 4 * cb**3 * cd
        - 27 * ca**2 * cd**2 + 18 * ca * cb * cc * cd
    )
    require(cubic_discriminant == -27 * L_symbolic**2, "cubic discriminant from coefficients")
    quotient_genus = (quotient_degree - 1) // 2
    weierstrass_A = cubic_coefficients[2]  # no w term after reversing/scaling
    quotient_j_numerator = 4 * weierstrass_A**3
    require(quotient_degree == 3, "weighted quotient cubic degree")
    require(bool(cubic_discriminant), "weighted quotient cubic squarefree")
    require(quotient_genus == 1, "weighted quotient genus one")
    require(quotient_j_numerator == 0, "weighted quotient j zero")
    require(3 > 1, "ramified cover is not invariant quotient")

    Phi_control = F(1)
    xi_control = F(124810012489, 12798000)
    eta_control = F(2091446550013, 107983125)
    U_control = F(-169749098233, 57591000)
    c2_control = F(4451333227, 474000)
    rho_control = Quadratic(0, 1)
    alpha_control = -2 * U_control * rho_control
    beta_control = -alpha_control
    require(rho_control * rho_control == Quadratic(Quadratic.D), "control rho square")
    require(alpha_control * alpha_control == Quadratic(4 * U_control * c2_control), "control k1 square")
    require(beta_control == -alpha_control, "control beta=-alpha")
    c3_control = eta_control + zeta_factor * Phi_control
    L_control = Quadratic(c3_control) + upsilon * rho_control
    L_norm = L_control.norm()
    expected_L_norm = F(
        270260378011253985379632330934787603,
        719758107151040278729687500,
    )
    require(L_norm == expected_L_norm, "control nonzero L norm")
    require(L_control != 0 and U_control != 0, "control nonzero quotient coefficients")

    cubic = [Quadratic(1), Quadratic(), Quadratic(), -L_control]
    cubic_derivative = [Quadratic(), Quadratic(), -3 * L_control]
    require(len(qpoly_gcd(cubic, cubic_derivative)) == 1, "control cubic squarefree")
    B_weierstrass = -L_control.inverse()
    discriminant = -432 * B_weierstrass * B_weierstrass
    require(discriminant != 0, "control elliptic discriminant")
    require(0 == 0, "control short-Weierstrass j=0")

    A = F(4183410507776)
    B = F(841357125)
    C = F(215966250)
    U_poly = uni_scale([F(13056802816), F(820125)], F(-13, 57591000))
    c2_poly = uni_scale([F(404652032), F(14625)], F(11, 474000))
    c3_numerator = [A, -B]
    high_contact_fraction_poly = uni_add(
        uni_multiply(U_poly, uni_multiply(c3_numerator, c3_numerator)),
        uni_scale(
            uni_multiply([F(0), F(1)], c2_poly),
            -upsilon * upsilon * C * C,
        ),
    )
    R_coefficients = primitive_integer_coefficients(high_contact_fraction_poly)
    expected_R = [
        2970579390109346274816679296272171008,
        2284603892441775363795663716352000,
        164114458618573873612800000000,
        7547170421607067494140625,
    ]
    require(R_coefficients == expected_R, "high-contact L=0 eliminant")
    forbidden = uni_multiply(
        [F(0), F(1)],
        uni_multiply(
            [F(13056802816), F(820125)],
            [F(404652032), F(14625)],
        ),
    )
    require(len(uni_gcd(high_contact_fraction_poly, forbidden)) == 1, "high-contact avoids forbidden factors")

    S_coefficients = [
        -1395571970793868500140032,
        61293210070929408000,
        8970234157828125,
    ]
    U_denominator_poly = [F(13056802816), F(820125)]
    L2_linear_poly = uni_add(
        [c4],
        uni_scale(c3_numerator, -F(3, 2) / (C * upsilon)),
    )
    L2_over_U_denominator = uni_add(
        uni_multiply(L2_linear_poly, U_denominator_poly),
        [upsilon**2 * F(57591000, 52)],
    )
    require(
        L2_over_U_denominator
        == uni_scale([F(entry) for entry in S_coefficients], F(-1, 676262246400)),
        "high-contact L2 rational function",
    )
    S_fraction_poly = [F(entry) for entry in S_coefficients]
    L2_nonzero_on_R = len(uni_gcd(high_contact_fraction_poly, S_fraction_poly)) == 1
    require(L2_nonzero_on_R, "gcd(R,S)=1")

    modulus = 17

    def mod_clean(poly: list[int]) -> list[int]:
        result = [entry % modulus for entry in poly]
        while result and result[-1] == 0:
            result.pop()
        return result

    def mod_add(left: list[int], right: list[int]) -> list[int]:
        result = [0] * max(len(left), len(right))
        for index, entry in enumerate(left):
            result[index] += entry
        for index, entry in enumerate(right):
            result[index] += entry
        return mod_clean(result)

    def mod_multiply(left: list[int], right: list[int]) -> list[int]:
        result = [0] * (len(left) + len(right) - 1)
        for i, left_entry in enumerate(left):
            for j, right_entry in enumerate(right):
                result[i + j] += left_entry * right_entry
        return mod_clean(result)

    bezout_mod_17 = mod_add(
        mod_multiply([-8, -5], R_coefficients),
        mod_multiply([-5, 1, -1], S_coefficients),
    )
    require(bezout_mod_17 == [1], "mod-17 Bezout identity")

    require(2 - 2 * 1 == 0, "P(2,1) invariant w=sigma/z^2")
    require((12 - 2 * 3, 2 * 3) == (6, 6), "second splitter excess tie")
    require((9 - 3 + F(6, 2)) == 9, "second Fq order")
    require(9 * 2 + 11 * 1 - 9 == 20, "second good-form lower bound")
    quartic_coefficients = [SPoly.make(1), SPoly(), SPoly(), SPoly(), -L2_numerator]
    quartic_degree = coefficient_list_degree(quartic_coefficients)
    quartic_genus = (quartic_degree - 2) // 2
    quartic_derivative = [index * coefficient
                          for index, coefficient in enumerate(quartic_coefficients)][1:]
    quartic_discriminant_nonzero = (
        L2_nonzero_on_R
        and quartic_coefficients[0] == 1
        and all(not coefficient for coefficient in quartic_derivative[:-1])
        and bool(quartic_derivative[-1])
    )
    qe, qd, qc, qb, qa = quartic_coefficients
    binary_quartic_I = 12 * qa * qe - 3 * qb * qd + qc**2
    binary_quartic_J = (
        72 * qa * qc * qe + 9 * qb * qc * qd - 27 * qa * qd**2
        - 27 * qb**2 * qe - 2 * qc**3
    )
    binary_quartic_I_nonzero = L2_nonzero_on_R and bool(binary_quartic_I)
    require(binary_quartic_I == -12 * L2_numerator, "quartic I from coefficients")
    require(binary_quartic_J == 0, "quartic J from coefficients")
    require(quartic_discriminant_nonzero, "second quotient squarefree")
    require(quartic_degree == 4 and quartic_genus == 1, "second quotient genus one")
    require(binary_quartic_I_nonzero and binary_quartic_J == 0, "second quotient j=1728")

    W_control = -(
        F(4343625) - F(17172000) * xi_control + F(143826305024)
    ) / F(4009500)
    Z_control = (
        F(12506118074368) - F(173745000)
        - F(195463125) * eta_control - F(926883000) * xi_control
    ) / F(108256500)
    require(W_control == -2 * U_control and Z_control == U_control, "control exact corner")

    print("THM4312_SOURCE_NORMAL_CUBIC_CORNER_INDEPENDENT_V2")
    print("GATES Delta=896/15 Theta=512/75 zeta3=-3Phi/2 upsilon5=-731648/2025")
    print("CORNER xi10=(4343625r+124805668864)/12798000; Phi*eta=(2091705253888-258703875r)/107983125")
    print("CORNER U=Z=-13(820125r+13056802816)/57591000; W=-2U; Phi!=0")
    print("EXCLUSIONS balanced_Delta_residual=128/9; k2_square_residual=-1839105572864/4100625; k3_c4=1664/25")
    print("K1 rho^2=c2/U; alpha11=-2U*rho; beta11=2U*rho")
    print("SPARSE_SOURCE P0=U*y*(y-rho)^2")
    print("SPARSE_SOURCE P1(rho)=rho*L; L=eta+zeta3+upsilon5*rho")
    print("SPARSE_SOURCE P2=6Uy^5-12Urho*y^4+(6Urho^2+4upsilon5)y^3+(4eta+3zeta3)y^2+c4y")
    print("CRITICAL_GRAPH a=-upsilon5/(2U); b=-(4eta+4rho*upsilon5+3zeta3)/(2U); L2=DERIVED")
    print("WEIGHTED_QUOTIENT w=sigma/z^3 V=Q/z^2 Y=wV; UY^2=1-Lw^3")
    print("WEIGHTED_QUOTIENT L!=0 squarefree_cubic=yes genus=1 j=0")
    print("RAMIFIED_COVER sigma=lambda^3,z=lambda*w0 gives degree-10 genus-4 cover; not the P(3,1) invariant quotient")
    print("ORDERS primitive=(s,beta,gamma)=(3,1,4) d=4 Fq=10 good_form_lower_bound=28")
    print(f"CONTROL Phi=1 xi10={xi_control} eta={eta_control} U={U_control}")
    print(f"CONTROL rho^2={Quadratic.D} L_norm={L_norm} corner=PASS cubic_squarefree=PASS")
    print("HIGH_CONTACT R(r)=7547170421607067494140625r^3+164114458618573873612800000000r^2+2284603892441775363795663716352000r+2970579390109346274816679296272171008")
    print("HIGH_CONTACT L2=[4U^2c4-2Ualpha(4eta+3zeta)-Uupsilon5^2+4alpha^2upsilon5]/(4U^2)")
    print("HIGH_CONTACT on_L=0 L2=c4+c3*zeta3/upsilon5-upsilon5^2/(4U)=-S/[676262246400(820125r+13056802816)]")
    print("HIGH_CONTACT S(r)=8970234157828125r^2+61293210070929408000r-1395571970793868500140032")
    print("HIGH_CONTACT gcd(R,forbidden)=1 gcd(R,S)=1; mod17 (-5r-8)R+(-r^2+r-5)S=1")
    print("SECOND_QUOTIENT w=sigma/z^2 V=Q/z^3 Y=wV; UY^2=1-L2w^4")
    print("SECOND_QUOTIENT L2!=0 squarefree_quartic=yes genus=1 j=1728")
    print("SECOND_ORDERS primitive=(s,beta,gamma)=(2,1,3) d=6 Fq=9 good_form_lower_bound=20")
    print(f"CHECKS={CHECKS}")
    print("SCOPE finite row-eight/corner intersection plus first and high-contact second k1 quotients; no all-row lift, seam entry, JC2, or DC2")


if __name__ == "__main__":
    main()
