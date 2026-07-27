#!/usr/bin/env python3
"""Exact probe for the degree-18 central-wall R=0 closure.

This companion is deliberately theorem-neutral: it reserves no theorem ID and
uses only the Python standard library.  It verifies the algebraic data needed
for the audited reduction

    D = 25 B^2 / 126,  R = 20 B C + 21 W = 0,

including the two singular weighted ratios, nodal Hessians, S != 0, the
normalization by lines through the node, the first-flux square function, and
the finite-field gcd certificate giving twelve branch points.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import reduce
from math import gcd


@dataclass(frozen=True)
class K:
    """Element a+b*s of Q(s), s^2=3."""

    a: Q = Q(0)
    b: Q = Q(0)

    @staticmethod
    def coerce(value: object) -> "K":
        if isinstance(value, K):
            return value
        return K(Q(value), Q(0))

    def __add__(self, other: object) -> "K":
        other = K.coerce(other)
        return K(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "K":
        return K(-self.a, -self.b)

    def __sub__(self, other: object) -> "K":
        return self + (-K.coerce(other))

    def __rsub__(self, other: object) -> "K":
        return K.coerce(other) - self

    def __mul__(self, other: object) -> "K":
        other = K.coerce(other)
        return K(
            self.a * other.a + 3 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    __rmul__ = __mul__

    def inverse(self) -> "K":
        norm = self.a * self.a - 3 * self.b * self.b
        if norm == 0:
            raise ZeroDivisionError("zero divisor in Q(sqrt(3))")
        return K(self.a / norm, -self.b / norm)

    def __truediv__(self, other: object) -> "K":
        return self * K.coerce(other).inverse()

    def __rtruediv__(self, other: object) -> "K":
        return K.coerce(other) / self

    def __pow__(self, exponent: int) -> "K":
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        out = K(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                out *= base
            base *= base
            n >>= 1
        return out

    def conjugate(self) -> "K":
        return K(self.a, -self.b)

    def norm(self) -> Q:
        return self.a * self.a - 3 * self.b * self.b

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def __str__(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*sqrt(3)"
        sign = "+" if self.b > 0 else "-"
        return f"{self.a}{sign}{abs(self.b)}*sqrt(3)"


ZERO = K(0)
ONE = K(1)


def trim(poly: list[K]) -> list[K]:
    poly = list(poly)
    while poly and poly[-1].is_zero():
        poly.pop()
    return poly


def pconst(value: object) -> list[K]:
    value = K.coerce(value)
    return [] if value.is_zero() else [value]


def padd(left: list[K], right: list[K]) -> list[K]:
    out = [ZERO] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] = out[index] + value
    for index, value in enumerate(right):
        out[index] = out[index] + value
    return trim(out)


def pneg(poly: list[K]) -> list[K]:
    return [-value for value in poly]


def psub(left: list[K], right: list[K]) -> list[K]:
    return padd(left, pneg(right))


def pscale(poly: list[K], scalar: object) -> list[K]:
    scalar = K.coerce(scalar)
    return trim([scalar * value for value in poly])


def pmul(left: list[K], right: list[K]) -> list[K]:
    if not left or not right:
        return []
    out = [ZERO] * (len(left) + len(right) - 1)
    for i, x_value in enumerate(left):
        for j, y_value in enumerate(right):
            out[i + j] = out[i + j] + x_value * y_value
    return trim(out)


def ppow(poly: list[K], exponent: int) -> list[K]:
    if exponent < 0:
        raise ValueError("polynomial exponent must be nonnegative")
    out = pconst(1)
    base = poly
    n = exponent
    while n:
        if n & 1:
            out = pmul(out, base)
        base = pmul(base, base)
        n >>= 1
    return out


def pdivmod(dividend: list[K], divisor: list[K]) -> tuple[list[K], list[K]]:
    dividend = trim(dividend)
    divisor = trim(divisor)
    if not divisor:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [ZERO] * max(0, len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] = (
                remainder[index + shift] - coefficient * value
            )
        remainder = trim(remainder)
    return trim(quotient), remainder


def pmonic(poly: list[K]) -> list[K]:
    poly = trim(poly)
    return pscale(poly, poly[-1].inverse()) if poly else []


def pgcd(left: list[K], right: list[K]) -> list[K]:
    left = trim(left)
    right = trim(right)
    while right:
        _, remainder = pdivmod(left, right)
        left, right = right, remainder
    return pmonic(left)


def pderivative(poly: list[K]) -> list[K]:
    return trim([K(index) * poly[index] for index in range(1, len(poly))])


def pequal(left: list[K], right: list[K]) -> bool:
    return not trim(psub(left, right))


@dataclass(frozen=True)
class RF:
    """Reduced rational function over Q(sqrt(3))[lambda]."""

    numerator: tuple[K, ...]
    denominator: tuple[K, ...]

    @staticmethod
    def make(
        numerator: list[K], denominator: list[K] | None = None
    ) -> "RF":
        if denominator is None:
            denominator = pconst(1)
        numerator = trim(numerator)
        denominator = trim(denominator)
        if not denominator:
            raise ZeroDivisionError("zero rational-function denominator")
        if not numerator:
            return RF(tuple(), tuple(pconst(1)))
        common = pgcd(numerator, denominator)
        numerator, rem_num = pdivmod(numerator, common)
        denominator, rem_den = pdivmod(denominator, common)
        assert not rem_num and not rem_den
        scale = denominator[-1].inverse()
        numerator = pscale(numerator, scale)
        denominator = pscale(denominator, scale)
        return RF(tuple(numerator), tuple(denominator))

    @staticmethod
    def constant(value: object) -> "RF":
        return RF.make(pconst(value))

    def __add__(self, other: object) -> "RF":
        other = rf_coerce(other)
        numerator = padd(
            pmul(list(self.numerator), list(other.denominator)),
            pmul(list(other.numerator), list(self.denominator)),
        )
        denominator = pmul(list(self.denominator), list(other.denominator))
        return RF.make(numerator, denominator)

    __radd__ = __add__

    def __neg__(self) -> "RF":
        return RF.make(pneg(list(self.numerator)), list(self.denominator))

    def __sub__(self, other: object) -> "RF":
        return self + (-rf_coerce(other))

    def __rsub__(self, other: object) -> "RF":
        return rf_coerce(other) - self

    def __mul__(self, other: object) -> "RF":
        other = rf_coerce(other)
        return RF.make(
            pmul(list(self.numerator), list(other.numerator)),
            pmul(list(self.denominator), list(other.denominator)),
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "RF":
        other = rf_coerce(other)
        return RF.make(
            pmul(list(self.numerator), list(other.denominator)),
            pmul(list(self.denominator), list(other.numerator)),
        )

    def __rtruediv__(self, other: object) -> "RF":
        return rf_coerce(other) / self

    def __pow__(self, exponent: int) -> "RF":
        if exponent < 0:
            return RF.make(
                ppow(list(self.denominator), -exponent),
                ppow(list(self.numerator), -exponent),
            )
        return RF.make(
            ppow(list(self.numerator), exponent),
            ppow(list(self.denominator), exponent),
        )

    def is_zero(self) -> bool:
        return not self.numerator


def rf_coerce(value: object) -> RF:
    if isinstance(value, RF):
        return value
    return RF.constant(value)


def format_poly(poly: list[K], variable: str = "lambda") -> str:
    if not poly:
        return "0"
    terms: list[str] = []
    for degree, coefficient in enumerate(poly):
        if coefficient.is_zero():
            continue
        if degree == 0:
            terms.append(f"({coefficient})")
        elif degree == 1:
            terms.append(f"({coefficient})*{variable}")
        else:
            terms.append(f"({coefficient})*{variable}^{degree}")
    return " + ".join(terms)


def primitive_integer_pair_polynomial(poly: list[K]) -> list[tuple[int, int]]:
    denominator_lcm = 1
    for coefficient in poly:
        for value in (coefficient.a, coefficient.b):
            denominator_lcm = (
                denominator_lcm
                * value.denominator
                // gcd(denominator_lcm, value.denominator)
            )
    integer_pairs = [
        (
            int(coefficient.a * denominator_lcm),
            int(coefficient.b * denominator_lcm),
        )
        for coefficient in poly
    ]
    content = reduce(
        gcd,
        (abs(value) for pair in integer_pairs for value in pair),
    )
    integer_pairs = [
        (left // content, right // content)
        for left, right in integer_pairs
    ]
    return integer_pairs


def verify_blowup_identity() -> None:
    # Sparse multivariate polynomials in (V,y,B,C), exponent tuples length four.
    def m_add(
        left: dict[tuple[int, ...], int],
        right: dict[tuple[int, ...], int],
    ) -> dict[tuple[int, ...], int]:
        out = dict(left)
        for monomial, coefficient in right.items():
            out[monomial] = out.get(monomial, 0) + coefficient
            if out[monomial] == 0:
                del out[monomial]
        return out

    def m_mul(
        left: dict[tuple[int, ...], int],
        right: dict[tuple[int, ...], int],
    ) -> dict[tuple[int, ...], int]:
        out: dict[tuple[int, ...], int] = {}
        for monomial_left, coefficient_left in left.items():
            for monomial_right, coefficient_right in right.items():
                monomial = tuple(
                    x + y for x, y in zip(monomial_left, monomial_right)
                )
                out[monomial] = (
                    out.get(monomial, 0)
                    + coefficient_left * coefficient_right
                )
        return {key: value for key, value in out.items() if value}

    def m_scale(
        poly: dict[tuple[int, ...], int], scalar: int
    ) -> dict[tuple[int, ...], int]:
        return {key: scalar * value for key, value in poly.items()}

    V = {(1, 0, 0, 0): 1}
    y = {(0, 1, 0, 0): 1}
    B = {(0, 0, 1, 0): 1}
    C = {(0, 0, 0, 1): 1}
    y2 = m_mul(y, y)
    y3 = m_mul(y2, y)
    U = m_mul(y, V)
    p = m_mul(
        y2,
        m_add(m_scale(y2, 245), m_scale(B, 1890)),
    )
    q = m_mul(
        y2,
        m_add(
            m_add(m_scale(y3, 539), m_scale(m_mul(B, y), 11340)),
            m_scale(C, 183708),
        ),
    )
    substituted_F = m_add(
        m_add(m_mul(m_mul(U, U), U), m_scale(m_mul(p, U), 3)),
        m_scale(m_mul(y, q), -7),
    )
    strict_transform = {
        tuple(
            exponent - (3 if index == 1 else 0)
            for index, exponent in enumerate(monomial)
        ): coefficient
        for monomial, coefficient in substituted_F.items()
    }
    expected_H = m_add(
        m_add(
            m_add(m_mul(m_mul(V, V), V), m_scale(m_mul(y2, V), 735)),
            m_scale(m_mul(B, V), 5670),
        ),
        m_add(
            m_add(m_scale(y3, -3773), m_scale(m_mul(B, y), -79380)),
            m_scale(C, -1285956),
        ),
    )
    assert strict_transform == expected_H


def singular_data(sign: int) -> dict[str, K]:
    root = K(Q(-35, 2), Q(sign * 21, 2))
    z_value = K(Q(-7, 15), Q(sign * 7, 36))
    c_value = K(Q(119, 1944), Q(-sign * 7, 216))
    rho = c_value**2 / z_value**3
    expected_rho = K(
        Q(-250, 13041),
        Q(sign * 500, 117369),
    )
    assert 2 * root**2 + 70 * root - 49 == ZERO
    assert z_value == (10 * root - 77) / 540
    assert c_value == (7 - 3 * root) / 972
    assert rho == expected_rho
    # H/y^3 and both affine derivative equations at y=1.
    assert root**2 + 245 + 1890 * z_value == ZERO
    assert 1470 * root - 11319 - 79380 * z_value == ZERO
    assert (
        root**3
        + 735 * root
        + 5670 * z_value * root
        - 3773
        - 79380 * z_value
        - 1285956 * c_value
        == ZERO
    )
    hessian = (
        6 * root * (1470 * root - 22638) - 1470**2
    )
    expected_hessian = 1166886 * K(5, -sign * 4)
    assert hessian == expected_hessian
    assert hessian.norm() == Q(-23) * 1166886**2
    s_over_b5 = 2888 + 244944 * rho
    expected_s = K(Q(-41576, 23), Q(sign * 24000, 23))
    assert s_over_b5 == expected_s
    assert s_over_b5.norm() == Q(24512, 23)
    return {
        "r": root,
        "B": z_value,
        "C": c_value,
        "rho": rho,
        "hessian": hessian,
        "S_over_B5": s_over_b5,
    }


def parameterized_first_flux(data: dict[str, K]) -> tuple[list[K], list[K]]:
    lam = [ZERO, ONE]
    root = data["r"]
    b_value = data["B"]
    c_value = data["C"]
    d_value = Q(25, 126) * b_value**2

    q2 = padd(
        padd(
            pscale(ppow(lam, 2), 3 * root),
            pscale(lam, 1470),
        ),
        pconst((1470 * root - 22638) / 2),
    )
    q3 = padd(
        padd(ppow(lam, 3), pscale(lam, 735)),
        pconst(-3773),
    )
    y_num = psub(q3, q2)
    v_num = psub(pscale(q3, root), pmul(lam, q2))
    y_value = RF.make(y_num, q3)
    v_value = RF.make(v_num, q3)

    # Verify the line-through-node normalization exactly.
    h_value = (
        v_value**3
        + 735 * y_value**2 * v_value
        + v_value * (5670 * b_value)
        - 3773 * y_value**3
        - y_value * (79380 * b_value)
        - 1285956 * c_value
    )
    assert h_value.is_zero()

    u_value = (
        35 * y_value**2
        + 1080 * b_value
        - 4 * y_value * v_value
    ) / 1701
    n2_value = (
        45927 * u_value**2
        - u_value * (58320 * b_value)
        - 5670 * u_value * y_value**2
        + y_value**2 * (2160 * b_value)
        + 93312 * d_value
        - y_value * (15552 * c_value)
        + 35 * y_value**4
    )
    z_square = n2_value * -Q(2, 5103) / y_value

    expected_denominator = ppow(q3, 3)
    assert pequal(list(z_square.denominator), expected_denominator)
    primitive_numerator = primitive_integer_pair_polynomial(
        list(z_square.numerator)
    )
    n_poly = [K(left, right) for left, right in primitive_numerator]
    expected_z = RF.make(
        pscale(n_poly, Q(8, 729)),
        expected_denominator,
    )
    assert z_square == expected_z

    # The parameterization formulas themselves.
    expected_q2 = padd(
        padd(
            pscale(ppow(lam, 2), 3 * root),
            pscale(lam, 1470),
        ),
        pconst((1470 * root - 22638) / 2),
    )
    assert pequal(q2, expected_q2)
    assert pequal(q3, [K(-3773), K(735), ZERO, ONE])
    return n_poly, q3


def finite_field_certificate(
    integer_pairs: list[tuple[int, int]],
) -> tuple[list[tuple[int, int]], tuple[int, int, int]]:
    prime = 5

    def fadd(
        left: tuple[int, int], right: tuple[int, int]
    ) -> tuple[int, int]:
        return (
            (left[0] + right[0]) % prime,
            (left[1] + right[1]) % prime,
        )

    def fneg(value: tuple[int, int]) -> tuple[int, int]:
        return (-value[0] % prime, -value[1] % prime)

    def fmul(
        left: tuple[int, int], right: tuple[int, int]
    ) -> tuple[int, int]:
        return (
            (left[0] * right[0] + 3 * left[1] * right[1]) % prime,
            (left[0] * right[1] + left[1] * right[0]) % prime,
        )

    def finverse(value: tuple[int, int]) -> tuple[int, int]:
        norm = (value[0] * value[0] - 3 * value[1] * value[1]) % prime
        if norm == 0:
            raise ZeroDivisionError("zero in F_25")
        inverse_norm = pow(norm, -1, prime)
        return (
            value[0] * inverse_norm % prime,
            -value[1] * inverse_norm % prime,
        )

    def ftrim(poly: list[tuple[int, int]]) -> list[tuple[int, int]]:
        poly = list(poly)
        while poly and poly[-1] == (0, 0):
            poly.pop()
        return poly

    def fdivmod(
        dividend: list[tuple[int, int]],
        divisor: list[tuple[int, int]],
    ) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        dividend = ftrim(dividend)
        divisor = ftrim(divisor)
        quotient = [(0, 0)] * max(0, len(dividend) - len(divisor) + 1)
        remainder = list(dividend)
        inverse_lead = finverse(divisor[-1])
        while len(remainder) >= len(divisor):
            shift = len(remainder) - len(divisor)
            coefficient = fmul(remainder[-1], inverse_lead)
            quotient[shift] = coefficient
            for index, value in enumerate(divisor):
                remainder[index + shift] = fadd(
                    remainder[index + shift],
                    fneg(fmul(coefficient, value)),
                )
            remainder = ftrim(remainder)
        return ftrim(quotient), remainder

    def fgcd(
        left: list[tuple[int, int]],
        right: list[tuple[int, int]],
    ) -> list[tuple[int, int]]:
        left = ftrim(left)
        right = ftrim(right)
        while right:
            _, remainder = fdivmod(left, right)
            left, right = right, remainder
        inverse_lead = finverse(left[-1])
        return [fmul(value, inverse_lead) for value in left]

    reduced_n = [
        (left % prime, right % prime)
        for left, right in integer_pairs
    ]
    derivative_n = [
        (
            degree * integer_pairs[degree][0] % prime,
            degree * integer_pairs[degree][1] % prime,
        )
        for degree in range(1, len(integer_pairs))
    ]
    reduced_q3 = [(2, 0), (0, 0), (0, 0), (1, 0)]
    derivative_q3 = [(0, 0), (0, 0), (3, 0)]
    gcd_degrees = (
        len(fgcd(reduced_n, derivative_n)) - 1,
        len(fgcd(reduced_n, reduced_q3)) - 1,
        len(fgcd(reduced_q3, derivative_q3)) - 1,
    )
    return reduced_n, gcd_degrees


def main() -> None:
    verify_blowup_identity()
    print("central-wall substitution: F(yV,y)=y^3 H(V,y) [PASS]")

    infinity_discriminant = -4 * 735**3 - 27 * 3773**2
    assert infinity_discriminant == -1972620783
    print(f"infinity cubic discriminant = {infinity_discriminant} [PASS]")
    print("H_y(V,0) = -79380 B, so B!=0 makes every exceptional point smooth")

    plus = singular_data(+1)
    minus = singular_data(-1)
    assert plus["rho"].conjugate() == minus["rho"]
    assert plus["rho"] != minus["rho"]
    for label, data in (("+", plus), ("-", minus)):
        print(f"rho_{label} = {data['rho']}")
        print(f"  r_{label} = {data['r']}")
        print(f"  B/y^2 = {data['B']}")
        print(f"  C/y^3 = {data['C']}")
        print(f"  Hessian/y^2 = {data['hessian']}")
        print(f"  S/B^5 = {data['S_over_B5']}")
    print("two conjugate affine singular ratios and nonzero nodal Hessians [PASS]")
    print("S is nonzero at both singular ratios [PASS]")

    tangent_ratio = Q(-250, 15309)
    tangent_s_over_b5 = 2888 + 244944 * tangent_ratio
    assert tangent_s_over_b5 == -1112
    print(
        "central tangent-discriminant ratio = "
        f"{tangent_ratio}; S/B^5 = {tangent_s_over_b5} [CORRECTED WARNING]"
    )

    plus_n, q3 = parameterized_first_flux(plus)
    minus_n, minus_q3 = parameterized_first_flux(minus)
    assert pequal(q3, minus_q3)
    assert minus_n == [value.conjugate() for value in plus_n]

    expected_pairs = [
        (3967243811384, -2293376193024),
        (916257470940, -532736790012),
        (10954768986, -2475570258),
        (-11405129358, 6621991614),
        (185297175, -319618719),
        (18237996, 10847718),
        (-3634428, 1963332),
        (114660, -82908),
        (-5355, 3339),
        (-22, 12),
    ]
    assert [(int(value.a), int(value.b)) for value in plus_n] == expected_pairs
    print("line-through-node parameterization H=0 [PASS]")
    print("Q3(lambda) = lambda^3 + 735 lambda - 3773")
    print("Z_+(lambda) = (8/729) N_+(lambda) / Q3(lambda)^3")
    print("N_+(lambda) coefficients, low degree first:")
    for degree, value in enumerate(plus_n):
        print(f"  {degree}: {value}")
    print("N_- is the sqrt(3)-conjugate of N_+ [PASS]")

    reduced_n, gcd_degrees = finite_field_certificate(expected_pairs)
    expected_reduction = [
        (4, 1),
        (0, 3),
        (1, 2),
        (2, 4),
        (0, 1),
        (1, 3),
        (2, 2),
        (0, 2),
        (0, 4),
        (3, 2),
    ]
    assert reduced_n == expected_reduction
    assert gcd_degrees == (0, 0, 0)
    print("F_25 reduction of N_+, low degree first:")
    print(" ", reduced_n)
    print("gcd degrees (N,N'), (N,Q3), (Q3,Q3') =", gcd_degrees)
    print("finite-field squarefree/coprime certificate [PASS]")

    numerator_branch_points = len(plus_n) - 1
    denominator_branch_points = len(q3) - 1
    infinity_branch_points = 0
    branch_points = (
        numerator_branch_points
        + denominator_branch_points
        + infinity_branch_points
    )
    assert (numerator_branch_points, denominator_branch_points) == (9, 3)
    assert branch_points == 12
    genus = (branch_points - 2) // 2
    assert genus == 5
    print(
        "square cover branch count = "
        f"{numerator_branch_points}+{denominator_branch_points}"
        f"+{infinity_branch_points} = {branch_points}"
    )
    print(f"square-cover genus = {genus} [PASS]")
    print("ALL EXACT R=0 CENTRAL-WALL PROBE CHECKS PASSED")


if __name__ == "__main__":
    main()
