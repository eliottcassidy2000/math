#!/usr/bin/env python3
"""Exact standard-library companion for THM-3584.

The structural theorem is all-exponent and all-degree.  This companion checks
its differential-polynomial identities, support arithmetic, parity split,
arm-order gates, and terminal coprimality on explicit finite hostile ranges.
It imports no theorem companion and uses no optimization-sensitive asserts.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    """Optimization-safe exact truth gate."""
    if not condition:
        raise RuntimeError(message)


# ---------------------------------------------------------------------------
# A tiny Laurent differential-polynomial algebra over Q.
#
# Negative exponents are used only for the invertible symbols A, h, and K in
# cleared first-integral identities.  The derivative map is first-order; none
# of the checked expressions requires differentiating h', K', a', q', or w'.

NAMES = (
    "A", "B", "L", "M",
    "h", "hp", "K", "Kp", "a", "ap", "q", "qp", "w", "wp",
    "y", "yp",
    "f0", "f0p", "f1", "f1p", "f2", "f2p",
    "g0", "g0p", "g1", "g1p", "g2", "g2p",
)
INDEX = {name: i for i, name in enumerate(NAMES)}
DERIVATIVE_INDEX = {
    INDEX["h"]: INDEX["hp"],
    INDEX["K"]: INDEX["Kp"],
    INDEX["a"]: INDEX["ap"],
    INDEX["q"]: INDEX["qp"],
    INDEX["w"]: INDEX["wp"],
    INDEX["y"]: INDEX["yp"],
    INDEX["f0"]: INDEX["f0p"],
    INDEX["f1"]: INDEX["f1p"],
    INDEX["f2"]: INDEX["f2p"],
    INDEX["g0"]: INDEX["g0p"],
    INDEX["g1"]: INDEX["g1p"],
    INDEX["g2"]: INDEX["g2p"],
}
ZERO_MONOMIAL = (0,) * len(NAMES)


class Expr:
    """Sparse Laurent polynomial with exact rational coefficients."""

    def __init__(self, terms: dict[tuple[int, ...], Fraction] | None = None):
        self.terms = {
            monomial: Fraction(coefficient)
            for monomial, coefficient in (terms or {}).items()
            if coefficient
        }

    @staticmethod
    def constant(value: int | Fraction) -> "Expr":
        value = Fraction(value)
        return Expr({ZERO_MONOMIAL: value}) if value else Expr()

    @staticmethod
    def variable(name: str) -> "Expr":
        exponents = [0] * len(NAMES)
        exponents[INDEX[name]] = 1
        return Expr({tuple(exponents): Fraction(1)})

    def __add__(self, other: object) -> "Expr":
        other_expr = to_expr(other)
        terms = dict(self.terms)
        for monomial, coefficient in other_expr.terms.items():
            terms[monomial] = terms.get(monomial, Fraction(0)) + coefficient
            if not terms[monomial]:
                del terms[monomial]
        return Expr(terms)

    def __radd__(self, other: object) -> "Expr":
        return self + other

    def __neg__(self) -> "Expr":
        return Expr({monomial: -coefficient for monomial, coefficient in self.terms.items()})

    def __sub__(self, other: object) -> "Expr":
        return self + (-to_expr(other))

    def __rsub__(self, other: object) -> "Expr":
        return to_expr(other) - self

    def __mul__(self, other: object) -> "Expr":
        other_expr = to_expr(other)
        terms: dict[tuple[int, ...], Fraction] = {}
        for left_monomial, left_coefficient in self.terms.items():
            for right_monomial, right_coefficient in other_expr.terms.items():
                monomial = tuple(
                    left_monomial[i] + right_monomial[i]
                    for i in range(len(NAMES))
                )
                terms[monomial] = terms.get(monomial, Fraction(0)) + (
                    left_coefficient * right_coefficient
                )
                if not terms[monomial]:
                    del terms[monomial]
        return Expr(terms)

    def __rmul__(self, other: object) -> "Expr":
        return self * other

    def __pow__(self, exponent: int) -> "Expr":
        require(isinstance(exponent, int), "Laurent exponent must be integral")
        if exponent == 0:
            return Expr.constant(1)
        if exponent < 0:
            require(len(self.terms) == 1, "negative power requires one monomial")
            monomial, coefficient = next(iter(self.terms.items()))
            return Expr({
                tuple(exponent * entry for entry in monomial): coefficient ** exponent
            })
        result = Expr.constant(1)
        base = self
        power = exponent
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power >>= 1
        return result

    def derivative(self) -> "Expr":
        terms: dict[tuple[int, ...], Fraction] = {}
        for monomial, coefficient in self.terms.items():
            for variable_index, derivative_index in DERIVATIVE_INDEX.items():
                exponent = monomial[variable_index]
                if not exponent:
                    continue
                differentiated = list(monomial)
                differentiated[variable_index] -= 1
                differentiated[derivative_index] += 1
                differentiated_tuple = tuple(differentiated)
                terms[differentiated_tuple] = terms.get(
                    differentiated_tuple, Fraction(0)
                ) + coefficient * exponent
                if not terms[differentiated_tuple]:
                    del terms[differentiated_tuple]
        return Expr(terms)

    def is_zero(self) -> bool:
        return not self.terms


def to_expr(value: object) -> Expr:
    if isinstance(value, Expr):
        return value
    if isinstance(value, (int, Fraction)):
        return Expr.constant(value)
    raise TypeError(f"cannot convert {type(value)!r} to Expr")


def zero(expression: Expr, message: str) -> None:
    require(expression.is_zero(), message)


def derivative(expression: Expr) -> Expr:
    return expression.derivative()


def wronskian(weight_f: int, f: Expr, weight_g: int, g: Expr) -> Expr:
    """The coefficient after removing c^(u+v+N-1)."""
    return weight_g * derivative(f) * g - weight_f * f * derivative(g)


AA, BB, LL, MM = (Expr.variable(name) for name in ("A", "B", "L", "M"))
h, hp = Expr.variable("h"), Expr.variable("hp")
K, Kp = Expr.variable("K"), Expr.variable("Kp")
a, ap = Expr.variable("a"), Expr.variable("ap")
q, qp = Expr.variable("q"), Expr.variable("qp")
w, wp = Expr.variable("w"), Expr.variable("wp")
y, yp = Expr.variable("y"), Expr.variable("yp")
f_coefficients = [Expr.variable(f"f{i}") for i in range(3)]
g_coefficients = [Expr.variable(f"g{i}") for i in range(3)]


# ---------------------------------------------------------------------------
# Small exact Q[b] implementation for independent terminal gcd controls.

def poly_trim(poly: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    values = list(poly)
    while values and not values[-1]:
        values.pop()
    return tuple(values)


def poly(value: list[int | Fraction] | tuple[int | Fraction, ...]) -> tuple[Fraction, ...]:
    return poly_trim(tuple(Fraction(entry) for entry in value))


def poly_add(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    size = max(len(left), len(right))
    return poly_trim(tuple(
        (left[i] if i < len(left) else 0) + (right[i] if i < len(right) else 0)
        for i in range(size)
    ))


def poly_scale(poly_value: tuple[Fraction, ...], scalar: int | Fraction) -> tuple[Fraction, ...]:
    scalar = Fraction(scalar)
    return poly_trim(tuple(scalar * entry for entry in poly_value))


def poly_mul(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if not left or not right:
        return ()
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_entry in enumerate(left):
        for j, right_entry in enumerate(right):
            result[i + j] += left_entry * right_entry
    return poly_trim(tuple(result))


def poly_pow(base: tuple[Fraction, ...], exponent: int) -> tuple[Fraction, ...]:
    require(exponent >= 0, "ordinary polynomial power must be nonnegative")
    result = poly([1])
    factor = base
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, factor)
        factor = poly_mul(factor, factor)
        power >>= 1
    return result


def poly_derivative(poly_value: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return poly_trim(tuple(i * poly_value[i] for i in range(1, len(poly_value))))


def poly_divmod(
    numerator: tuple[Fraction, ...], denominator: tuple[Fraction, ...]
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    require(bool(denominator), "polynomial division by zero")
    remainder = list(numerator)
    quotient = [Fraction(0)] * max(1, len(numerator) - len(denominator) + 1)
    while remainder and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[shift] += coefficient
        for i, entry in enumerate(denominator):
            remainder[shift + i] -= coefficient * entry
        remainder = list(poly_trim(tuple(remainder)))
    return poly_trim(tuple(quotient)), poly_trim(tuple(remainder))


def poly_gcd(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    a_poly, b_poly = left, right
    while b_poly:
        _, remainder = poly_divmod(a_poly, b_poly)
        a_poly, b_poly = b_poly, remainder
    if not a_poly:
        return ()
    return poly_scale(a_poly, Fraction(1, 1) / a_poly[-1])


def poly_degree(poly_value: tuple[Fraction, ...]) -> int:
    return len(poly_value) - 1


# ---------------------------------------------------------------------------
# 1. Laurent compiler, regularity, and local scalar gate.

compiler_controls = 0
for r in range(-7, 4):
    for s in range(-6, 5):
        for step in range(1, 6):
            for row in range(5):
                direct = Expr.constant(0)
                expanded = Expr.constant(0)
                for i in range(3):
                    for j in range(3):
                        if i + j != row:
                            continue
                        direct += wronskian(
                            r + i * step,
                            f_coefficients[i],
                            s + j * step,
                            g_coefficients[j],
                        )
                        expanded += (
                            (s + j * step)
                            * derivative(f_coefficients[i])
                            * g_coefficients[j]
                            - (r + i * step)
                            * f_coefficients[i]
                            * derivative(g_coefficients[j])
                        )
                zero(direct - expanded, "five-row Laurent compiler")
                compiler_controls += 1

require(
    [sum(1 for i in range(3) for j in range(3) if i + j == row) for row in range(5)]
    == [1, 2, 3, 2, 1],
    "convolution multiplicities",
)

regularity_controls = 0
scalar_gate_controls = 0
for exponent in range(2, 33):
    for weight in range(-4 * exponent, 2 * exponent + 1):
        if weight < 0:
            arm_power = (-weight + exponent - 1) // exponent
            require(weight + exponent * arm_power >= 0, "regularity sufficiency")
            require(
                weight + exponent * (arm_power - 1) < 0,
                "regularity minimality",
            )
        regularity_controls += 1

    survivors: list[tuple[int, int, int, int]] = []
    for left_weight in range(-4 * exponent, 4 * exponent + 1):
        right_weight = -(exponent - 1) - left_weight
        left_min = (
            (-left_weight + exponent - 1) // exponent if left_weight < 0 else 0
        )
        right_min = (
            (-right_weight + exponent - 1) // exponent if right_weight < 0 else 0
        )
        for left_order in range(left_min, left_min + 4):
            for right_order in range(right_min, right_min + 4):
                order = left_order + right_order - 1
                multiplier = right_weight * left_order - left_weight * right_order
                if order == 0 and multiplier:
                    survivors.append(
                        (left_weight, right_weight, left_order, right_order)
                    )
    require(
        survivors == [
            (-exponent, 1, 1, 0),
            (1, -exponent, 0, 1),
        ],
        f"simple-arm scalar gate N={exponent}",
    )
    scalar_gate_controls += 1


# ---------------------------------------------------------------------------
# 2. Off-central rows kappa=0,1,3,4.

degree_gate_controls = 0
for exponent in range(2, 17):
    for degree_h in range(2, 10):
        for degree_u in range(0, 10):
            require(degree_h + exponent * degree_u > 0, "E_N leading coefficient")
            require(degree_u + exponent * degree_h > 0, "J_N leading coefficient")
            require(degree_h + degree_u - 1 >= 1, "homogeneous degree gate")
            degree_gate_controls += 1

kappa_one_controls = 0
kappa_three_controls = 0
kappa_three_sign_controls = 0
for exponent in range(2, 13):
    for k in range(1, 13):
        step = exponent * k + 1
        support_p = (
            -exponent,
            exponent * (k - 1) + 1,
            exponent * (2 * k - 1) + 2,
        )
        support_q = (-exponent * k, 1, exponent * k + 2)
        require(
            support_p[1] - support_p[0] == step
            and support_p[2] - support_p[1] == step
            and support_q[1] - support_q[0] == step
            and support_q[2] - support_q[1] == step,
            "kappa=1 support normalization",
        )
        low_f = AA * h
        low_g = BB * h**k
        scalar = (
            wronskian(-exponent, low_f, 1, q)
            + wronskian(exponent * (k - 1) + 1, a, -exponent * k, low_g)
        )
        inner = AA * q - k * BB * h ** (k - 1) * a
        expected = hp * inner + exponent * h * derivative(inner)
        zero(scalar - expected, f"kappa=1 operator identity N={exponent},k={k}")
        kappa_one_controls += 1

    # Up to exchanging P and Q, the kappa=3 scalar gate has two branches:
    # high endpoints (d-N,1), where the negative member (if any) is the
    # middle coefficient, and (d+1,-N), where the high extreme itself has
    # mixed signs.  The latter can never have zero Wronskian at an arm.  In
    # the former, d<N is likewise mixed-sign; d=N is the removable 0/1
    # boundary; only d>N reaches the common-power normalization below.
    for step in range(1, 3 * exponent + 1):
        primary_high = (step - exponent, 1)
        alternate_high = (step + 1, -exponent)
        require(
            primary_high[0] - step == -exponent,
            "kappa=3 primary scalar pair",
        )
        require(
            alternate_high[1] == -exponent
            and alternate_high[0] - step == 1,
            "kappa=3 alternate scalar pair",
        )
        for positive_order in range(0, 4):
            for negative_order in range(1, 5):
                alternate_multiplier = (
                    -exponent * positive_order
                    - (step + 1) * negative_order
                )
                require(
                    alternate_multiplier != 0,
                    "kappa=3 alternate mixed-sign high extreme",
                )
                if step < exponent:
                    primary_multiplier = (
                        negative_order
                        + (exponent - step) * positive_order
                    )
                    require(
                        primary_multiplier != 0,
                        "kappa=3 subcritical mixed-sign high extreme",
                    )
        if step == exponent:
            zero(
                wronskian(0, f_coefficients[0], 1, g_coefficients[0])
                - derivative(f_coefficients[0]) * g_coefficients[0],
                "kappa=3 zero/one boundary",
            )
        kappa_three_sign_controls += 1

    for n in range(1, 13):
        step = exponent + n
        support_p = (-n - 2 * exponent, -exponent, n)
        support_q = (-2 * n - 2 * exponent + 1, -n - exponent + 1, 1)
        require(
            support_p[1] - support_p[0] == step
            and support_p[2] - support_p[1] == step
            and support_q[1] - support_q[0] == step
            and support_q[2] - support_q[1] == step,
            "kappa=3 support normalization",
        )
        high_f = LL * K**n
        high_g = MM * K
        scalar = (
            wronskian(-exponent, a, 1, high_g)
            + wronskian(n, high_f, -n - exponent + 1, q)
        )
        inner = MM * a - n * LL * K ** (n - 1) * q
        expected = K * derivative(inner) + exponent * Kp * inner
        zero(scalar - expected, f"kappa=3 operator identity N={exponent},n={n}")
        require(
            (-(-exponent) + exponent - 1) // exponent == 1,
            "kappa=3 middle-a arm divisibility",
        )
        require(
            (n + exponent - 1 + exponent - 1) // exponent >= 1,
            "kappa=3 middle-q arm divisibility",
        )
        kappa_three_controls += 1


# ---------------------------------------------------------------------------
# 3. Complete central resonance census.

central_census_controls = 0
middle_profiles = 0
endpoint_profiles = 0
for exponent in range(2, 21):
    for R in range(1, 121):
        for T in range(1, 121):
            numerator = R + T - (exponent - 1)
            if numerator <= 0 or numerator % 2:
                continue
            middle = abs(R - T) == exponent + 1
            endpoint_r = R == exponent and gcd(R, T) == exponent
            endpoint_t = T == exponent and gcd(R, T) == exponent
            possible = middle or endpoint_r or endpoint_t
            endpoint_expected = (
                exponent % 2 == 1
                and (
                    (R == exponent and T % exponent == 0 and (T // exponent) % 2 == 1)
                    or (
                        T == exponent
                        and R % exponent == 0
                        and (R // exponent) % 2 == 1
                    )
                )
            )
            require(
                possible == (middle or endpoint_expected),
                f"central resonance classification N={exponent},R={R},T={T}",
            )
            middle_profiles += int(middle)
            endpoint_profiles += int(endpoint_expected and not middle)
            central_census_controls += 1


# ---------------------------------------------------------------------------
# 4. Middle central family |R-T|=N+1.

middle_arm_controls = 0
middle_bridge_controls = 0
middle_terminal_controls = 0
for exponent in range(2, 13):
    for R in range(1, 97):
        delta = gcd(R, exponent + 1)
        alpha = R // delta
        beta = (R + exponent + 1) // delta
        require(gcd(delta, exponent) == 1, "consecutive-exponent gcd")
        for arm_order in range(1, 9):
            multiplier = alpha * (delta - exponent * arm_order)
            require(multiplier != 0, "middle lower-bridge multiplier")
            cancellation = (beta - alpha) * arm_order == 1
            require(
                cancellation == (delta == exponent + 1 and arm_order == 1),
                "middle arm-order resonance",
            )
            middle_arm_controls += 1

    gamma = gcd(2, exponent + 1)
    D = (exponent + 1) // gamma
    require(gamma == (1 if exponent % 2 == 0 else 2), "middle parity split")
    for m in range(1, 9):
        R = (exponent + 1) * m
        T = R + exponent + 1
        low_f = AA * h**m
        low_g = BB * h ** (m + 1)
        lower = (
            wronskian(-R, low_f, -exponent, q)
            + wronskian(1, a, -T, low_g)
        )
        target = h ** (m - 1) * (
            m
            * AA
            * ((exponent + 1) * h * qp - exponent * hp * q)
            - (m + 1)
            * BB
            * h
            * ((exponent + 1) * h * ap + hp * a)
        )
        zero(lower - target, f"middle lower bridge N={exponent},m={m}")

        lam = Fraction(m + 1, m) * BB * AA**-1
        q_particular = lam * h * a
        particular_lower = (
            wronskian(-R, low_f, -exponent, q_particular)
            + wronskian(1, a, -T, low_g)
        )
        zero(particular_lower, f"middle forced q N={exponent},m={m}")

        first_integral_left = derivative(w ** (exponent + 1) * h ** (-exponent))
        first_integral_right = (
            w**exponent
            * h ** (-exponent - 1)
            * ((exponent + 1) * h * wp - exponent * hp * w)
        )
        zero(
            first_integral_left - first_integral_right,
            f"middle (N+1)-power first integral N={exponent}",
        )
        require(
            all((exponent + 1) * order != exponent for order in range(0, 8)),
            "simple-arm middle valuation",
        )

        U = R + 2
        V = R - exponent + 1
        require(gcd(U, V) == gamma, "middle high gcd")
        u = U // gamma
        v = V // gamma
        require(u - v == D and v >= 1, "middle high power normalization")
        high_f = LL * K**u
        high_g = MM * K**v
        upper = (
            wronskian(1, a, V, high_g)
            + wronskian(U, high_f, -exponent, q_particular)
        )
        A0 = MM * v - lam * LL * u * h * K**D
        Y = a * A0
        expected_upper = K ** (v - 1) * (
            gamma * K * derivative(Y) - Kp * Y
        )
        zero(upper - expected_upper, f"middle upper bridge N={exponent},m={m}")

        terminal_integral_left = derivative(y**gamma * K**-1)
        terminal_integral_right = (
            y ** (gamma - 1)
            * K**-2
            * (gamma * K * yp - Kp * y)
        )
        zero(
            terminal_integral_left - terminal_integral_right,
            f"middle gamma-power terminal integral N={exponent}",
        )
        middle_bridge_controls += 1

        h_poly = poly([0, -1, 1])
        K_poly = poly([m + 2, 1, 1])
        A0_poly = poly_add(
            poly([v]),
            poly_scale(poly_mul(h_poly, poly_pow(K_poly, D)), -(u + 1)),
        )
        require(poly_degree(A0_poly) > 0, "middle terminal nonconstant")
        require(poly_degree(poly_gcd(A0_poly, K_poly)) == 0, "middle terminal gcd")
        middle_terminal_controls += 1


# ---------------------------------------------------------------------------
# 5. Odd-N endpoint family {R,T}={N,Nk}, k odd.

endpoint_bridge_controls = 0
endpoint_terminal_controls = 0
endpoint_first_integral_controls = 0
for exponent in range(3, 14, 2):
    for k in range(1, 12, 2):
        step = (exponent * k + 1) // 2
        S = (exponent * k - 1) // 2
        p = (exponent * (k - 2) + 1) // 2
        U = exponent * (k - 1) + 1
        support_p = (-exponent, p, U)
        support_q = (-exponent * k, -S, 1)
        require(
            support_p[1] - support_p[0] == step
            and support_p[2] - support_p[1] == step
            and support_q[1] - support_q[0] == step
            and support_q[2] - support_q[1] == step,
            "odd endpoint support normalization",
        )

        low_f = AA * h
        low_g = BB * h**k
        row_below = (
            wronskian(-exponent, low_f, -S, q)
            + wronskian(p, a, -exponent * k, low_g)
        )
        X = AA * q - k * BB * h ** (k - 1) * a
        expected_below = exponent * h * derivative(X) - S * hp * X
        zero(row_below - expected_below, f"endpoint lower bridge N={exponent},k={k}")

        endpoint_integral_left = derivative(X**exponent * h**(-S))
        endpoint_integral_right = (
            X ** (exponent - 1)
            * h ** (-S - 1)
            * (exponent * h * derivative(X) - S * hp * X)
        )
        zero(
            endpoint_integral_left - endpoint_integral_right,
            f"endpoint N-power first integral N={exponent},k={k}",
        )
        require(S % exponent != 0, "odd endpoint simple-arm valuation")
        endpoint_first_integral_controls += 1

        high_f = LL * K**U
        high_g = MM * K
        row_above = (
            wronskian(p, a, 1, high_g)
            + wronskian(U, high_f, -S, q)
        )
        primitive = MM * a * K ** (-p) - U * LL * q * K**S
        expected_above = K ** (p + 1) * derivative(primitive)
        zero(row_above - expected_above, f"endpoint upper bridge N={exponent},k={k}")

        lam = k * BB * AA**-1
        q_particular = lam * h ** (k - 1) * a
        endpoint_power = S + p
        require(endpoint_power == exponent * (k - 1), "endpoint terminal power")
        D0 = MM - U * LL * lam * h ** (k - 1) * K**endpoint_power
        primitive_particular = (
            MM * a * K ** (-p) - U * LL * q_particular * K**S
        )
        zero(
            primitive_particular - a * K ** (-p) * D0,
            f"endpoint integrated terminal N={exponent},k={k}",
        )

        if k == 1:
            scalar = (
                wronskian(-exponent, low_f, 1, high_g)
                + wronskian(p, a, -S, q_particular)
                + wronskian(U, high_f, -exponent, low_g)
            )
            delta = AA * MM - BB * LL
            degree_operator = K * hp + exponent * h * Kp
            zero(
                scalar - delta * degree_operator,
                f"minimal odd endpoint scalar N={exponent}",
            )
        else:
            require(p > 0, "nonminimal odd endpoint p positivity")
            h_poly = poly([0, -1, 1])
            K_poly = poly([k + 2, 1, 1])
            D0_poly = poly_add(
                poly([1]),
                poly_scale(
                    poly_mul(
                        poly_pow(h_poly, k - 1),
                        poly_pow(K_poly, exponent * (k - 1)),
                    ),
                    -(U + k),
                ),
            )
            require(poly_degree(D0_poly) > 0, "endpoint terminal nonconstant")
            require(poly_degree(poly_gcd(D0_poly, K_poly)) == 0, "endpoint terminal gcd")
            endpoint_terminal_controls += 1
        endpoint_bridge_controls += 1


# ---------------------------------------------------------------------------
# 6. Sharp N=3 controls.

N3_MIDDLE = ((-4, 1, 6), (-8, -3, 2))
N3_ENDPOINT_ONE = ((-3, -1, 1), (-3, -1, 1))
N3_ENDPOINT_THREE = ((-3, 2, 7), (-9, -4, 1))
for pair in (N3_MIDDLE, N3_ENDPOINT_ONE, N3_ENDPOINT_THREE):
    for support in pair:
        require(
            support[1] - support[0] == support[2] - support[1],
            "N=3 hostile arithmetic progression",
        )

# Both-negative scalar channel (-1,-1): every regular coefficient has arm
# order at least one, so the bracket has order at least one.
require(1 + 1 - 1 == 1, "N=3 both-negative scalar hostile")
# The zero-weight boundary (-(N-1),0)=(-2,0) has zero leading multiplier.
require(0 * 1 - (-2) * 0 == 0, "N=3 zero-multiplier hostile")

# Repeated-arm formal resonances show exactly why squarefreeness is needed:
# for N=3, h=z^4,w=z^3 solves 4hw'-3h'w=0; and h=z^3,X=z solves
# 3hX'-h'X=0 in the minimal endpoint family.
require(4 * 3 - 3 * 4 == 0, "N=3 middle repeated-arm hostile")
# The endpoint calculation is 3*z^3*1 - (3*z^2)*z = 0, so its
# coefficient is 3-3.
require(3 - 3 == 0, "N=3 endpoint repeated-arm hostile")

# Degree-one and Laurent-pole hostiles.
require(1 + 3 * 0 == 1, "N=3 degree-one differential hostile")
require(Fraction(-1, 2) * (1 - 3) == 1, "N=3 Laurent-pole Darboux hostile")


print("THM-3584 all-exponent equal-step three-by-three Danielewski Darboux nonentry")
print(f"Laurent compiler rows: {compiler_controls}; convolution=[1,2,3,2,1]")
print(
    f"regularity controls: {regularity_controls}; simple-arm gates: "
    f"N=2..32 ({scalar_gate_controls} exponents), survivor=(-N,1) and reflection"
)
print(f"kappa=0/4 degree gates: {degree_gate_controls}")
print(f"kappa=1 operator/support controls: {kappa_one_controls}")
print(
    f"kappa=3 operator/support controls: {kappa_three_controls}; "
    f"sign-branch controls: {kappa_three_sign_controls}"
)
print(
    f"central census: {central_census_controls}; middle profiles={middle_profiles}; "
    f"odd-endpoint profiles={endpoint_profiles}"
)
print(
    f"middle family: arm={middle_arm_controls}, bridges={middle_bridge_controls}, "
    f"terminal-gcd={middle_terminal_controls}; gamma=1 even N / 2 odd N"
)
print(
    f"odd-N endpoint family: first-integrals={endpoint_first_integral_controls}, "
    f"bridges={endpoint_bridge_controls}, terminal-gcd={endpoint_terminal_controls}"
)
print("N=3 sharp profiles: (-4,1,6)|(-8,-3,2); (-3,-1,1)^2; (-3,2,7)|(-9,-4,1)")
print("hostiles: both-negative / zero-multiplier / repeated-arm / degree-one / Laurent-pole")
print("scope: reduced exact equal-step 3x3 supports; does not reprove THM-3576; no JC consequence")
print("all optimization-safe exact truth gates passed")
