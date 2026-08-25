#!/usr/bin/env python3
"""Rigorous margin-function frontier certificate for THM-4089.

This script certifies two facts about the *displayed elementary margin
function* in Christopher D. Long's ``p-adic-zeta-irrationality`` manuscript:

1. the manuscript's 22 fixed rational witnesses have positive margins;
2. after exact minimization in xi, the globally concave Y-objective has a
   negative maximum in each immediate next case
   (2,31), (3,13), (5,7), (7,5).

The second assertion uses derivative brackets and concave tangent upper
bounds, not a grid search.  It is a theorem about the formula, not a
verification of the manuscript's geometric or adelic proof.

The outward-rounded fixed-point interval primitives are adapted from
``hybrid_rational_interval_certificate.py`` at upstream commit b46a177,
Copyright (c) 2026 Christopher D. Long, under the MIT License.  This file adds
the derivative enclosure, global concavity certificate, next-case hostiles,
and independent output format.  No binary floating point or third-party
package is used.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt


PRECISION_DIGITS = 80
PRINT_DIGITS = 30
SCALE = 10 ** PRECISION_DIGITS
SERIES_STOP_UNITS = 100
Interval = tuple[int, int]
GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def floor_div(a: int, b: int) -> int:
    return a // b


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def i_fraction(value: Fraction) -> Interval:
    return (
        floor_div(SCALE * value.numerator, value.denominator),
        ceil_div(SCALE * value.numerator, value.denominator),
    )


def i_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def i_sub(left: Interval, right: Interval) -> Interval:
    return left[0] - right[1], left[1] - right[0]


def i_mul(left: Interval, right: Interval) -> Interval:
    products = (
        left[0] * right[0], left[0] * right[1],
        left[1] * right[0], left[1] * right[1],
    )
    return floor_div(min(products), SCALE), ceil_div(max(products), SCALE)


def i_mul_int(value: Interval, scalar: int) -> Interval:
    if scalar >= 0:
        return value[0] * scalar, value[1] * scalar
    return value[1] * scalar, value[0] * scalar


def i_mul_fraction(value: Interval, scalar: Fraction) -> Interval:
    n, d = scalar.numerator, scalar.denominator
    if n >= 0:
        return floor_div(value[0] * n, d), ceil_div(value[1] * n, d)
    return floor_div(value[1] * n, d), ceil_div(value[0] * n, d)


def positive_mul_rational(value: Interval, numerator: int, denominator: int) -> Interval:
    require(0 <= value[0] <= value[1], "nonnegative interval multiplier")
    return (
        floor_div(value[0] * numerator, denominator),
        ceil_div(value[1] * numerator, denominator),
    )


def sqrt_fraction_bounds(value: Fraction) -> Interval:
    require(value >= 0, "nonnegative square-root input")
    root = isqrt((value.numerator * SCALE * SCALE) // value.denominator)
    while (root + 1) ** 2 * value.denominator <= value.numerator * SCALE * SCALE:
        root += 1
    while root ** 2 * value.denominator > value.numerator * SCALE * SCALE:
        root -= 1
    if root ** 2 * value.denominator == value.numerator * SCALE * SCALE:
        return root, root
    return root, root + 1


def atan_rational_bounds(numerator: int, denominator: int) -> Interval:
    require(0 <= numerator < denominator, "atan argument in [0,1)")
    if numerator == 0:
        return 0, 0
    term = (
        floor_div(SCALE * numerator, denominator),
        ceil_div(SCALE * numerator, denominator),
    )
    partial = (0, 0)
    k = 0
    plus = True
    while True:
        partial = i_add(partial, term) if plus else i_sub(partial, term)
        next_term = positive_mul_rational(
            term,
            numerator * numerator * (2 * k + 1),
            denominator * denominator * (2 * k + 3),
        )
        if next_term[1] <= SERIES_STOP_UNITS:
            if plus:
                return partial[0] - next_term[1], partial[1]
            return partial[0], partial[1] + next_term[1]
        term = next_term
        plus = not plus
        k += 1
        require(k < 1_000_000, "atan convergence")


def log_integer_bounds(p: int) -> Interval:
    require(p >= 2, "log integer domain")
    numerator, denominator = p - 1, p + 1
    term = (
        floor_div(SCALE * numerator, denominator),
        ceil_div(SCALE * numerator, denominator),
    )
    partial = (0, 0)
    k = 0
    while True:
        partial = i_add(partial, term)
        next_term = positive_mul_rational(
            term,
            numerator * numerator * (2 * k + 1),
            denominator * denominator * (2 * k + 3),
        )
        tail = ceil_div(
            next_term[1] * denominator * denominator,
            denominator * denominator - numerator * numerator,
        )
        if tail <= SERIES_STOP_UNITS:
            return 2 * partial[0], 2 * (partial[1] + tail)
        term = next_term
        k += 1
        require(k < 1_000_000, "atanh convergence")


def pi_bounds() -> Interval:
    return i_sub(i_mul_int(atan_rational_bounds(1, 5), 16),
                 i_mul_int(atan_rational_bounds(1, 239), 4))


def acos_fraction_bounds(value: Fraction) -> Interval:
    require(0 < value < 1, "acos argument in (0,1)")
    tangent = sqrt_fraction_bounds((1 - value) / (1 + value))
    lower = atan_rational_bounds(tangent[0], SCALE)
    upper = atan_rational_bounds(tangent[1], SCALE)
    return 2 * lower[0], 2 * upper[1]


def euler_phi(n: int) -> int:
    result, residual, prime = n, n, 2
    while prime * prime <= residual:
        if residual % prime == 0:
            while residual % prime == 0:
                residual //= prime
            result -= result // prime
        prime += 1
    if residual > 1:
        result -= result // residual
    return result


def harmonic(k: int) -> Fraction:
    return sum((Fraction(1, j) for j in range(1, k + 1)), Fraction(0))


def xi_star(p: int, s: int) -> Fraction:
    return Fraction(12 * (s * (s + 1) - 1), 12 * s * s + (s - 1) * (p + 1))


def prime_window(s: int, xi: Fraction) -> Fraction:
    m = s + 1
    k = Fraction(s, 1) // xi
    residual = max(Fraction(0), Fraction(m) - Fraction(k + 1) * xi)
    return (Fraction(2 * m - 1, 2) * harmonic(k) - k * xi
            + residual * residual / Fraction(2 * (k + 1)))


def small_cost(p: int, s: int, xi: Fraction) -> Fraction:
    return s * s * xi - (s - 1) * (xi - Fraction(p + 1, 24) * xi * xi)


def tau_star(p: int, s: int) -> tuple[Fraction, Fraction]:
    xi = xi_star(p, s)
    require(1 < xi < Fraction(s + 1, s), f"xi in positive K=s-1 cell {(p,s)}")
    require(xi < Fraction(s, s - 1), f"xi floor cell {(p,s)}")
    require(Fraction(p + 1, 12) * xi < 1, f"xi Hasse branch {(p,s)}")
    # Exact derivative of S+sI on this cell.
    derivative = (s * s - (s - 1) * (1 - Fraction(p + 1, 12) * xi)
                  + s * s * (xi - 2))
    require(derivative == 0, f"exact xi stationary point {(p,s)}")
    tau = Fraction(2, (s + 1) ** 2) * (
        small_cost(p, s, xi) + s * prime_window(s, xi)
    )
    return tau, xi


def collision_interval(p: int, y: Fraction, pi_i: Interval) -> Interval:
    collision = (0, 0)
    c = p
    while c * y < 1:
        x = c * y
        bracket = i_sub(
            acos_fraction_bounds(x),
            i_mul_fraction(sqrt_fraction_bounds(1 - x * x), x),
        )
        term = i_mul_fraction(
            i_mul(pi_i, bracket), Fraction(4 * euler_phi(c), c * c)
        )
        collision = i_add(collision, term)
        c += p
    return collision


def margin_interval(
    p: int,
    s: int,
    y: Fraction,
    pi_i: Interval,
    logs: dict[int, Interval],
) -> Interval:
    tau, _xi = tau_star(p, s)
    require(0 < y < Fraction(1, p), f"Y domain {(p,s)}")
    width = i_sub(
        i_mul_fraction(logs[p], Fraction(12, p - 1)),
        i_mul_fraction(pi_i, 2 * y),
    )
    return i_sub(
        i_sub(i_mul_int(width, s), i_fraction(Fraction(s + 1) * tau)),
        collision_interval(p, y, pi_i),
    )


def derivative_interval(p: int, s: int, y: Fraction, pi_i: Interval) -> Interval:
    # d/dY [s Lambda_p(Y)-C_p(Y)]
    radical_sum = (0, 0)
    c = p
    while c * y < 1:
        radical_sum = i_add(
            radical_sum,
            i_mul_fraction(sqrt_fraction_bounds(1 - (c * y) ** 2),
                           Fraction(euler_phi(c), c)),
        )
        c += p
    bracket = i_add(i_fraction(Fraction(-2 * s)), i_mul_int(radical_sum, 8))
    return i_mul(pi_i, bracket)


def tangent_upper(
    p: int,
    s: int,
    left: Fraction,
    right: Fraction,
    pi_i: Interval,
    logs: dict[int, Interval],
) -> int:
    require(left < right, "nonempty derivative bracket")
    m_left = margin_interval(p, s, left, pi_i, logs)
    m_right = margin_interval(p, s, right, pi_i, logs)
    d_left = derivative_interval(p, s, left, pi_i)
    d_right = derivative_interval(p, s, right, pi_i)
    require(d_left[0] > 0, f"positive left derivative {(p,s)}")
    require(d_right[1] < 0, f"negative right derivative {(p,s)}")
    width = right - left
    upper_from_left = m_left[1] + ceil_div(d_left[1] * width.numerator, width.denominator)
    # d_right<0 and left-right<0, so the tangent correction is positive.
    upper_from_right = m_right[1] + ceil_div(
        d_right[0] * (left - right).numerator,
        (left - right).denominator,
    )
    upper = min(upper_from_left, upper_from_right)
    require(upper < 0, f"negative global tangent upper {(p,s)}")
    return upper


def decimal_lower(value: int, digits: int = PRINT_DIGITS) -> str:
    factor = 10 ** (PRECISION_DIGITS - digits)
    scaled = floor_div(value, factor)
    sign = "-" if scaled < 0 else ""
    scaled = abs(scaled)
    return f"{sign}{scaled // 10**digits}.{scaled % 10**digits:0{digits}d}"


def decimal_upper(value: int, digits: int = PRINT_DIGITS) -> str:
    factor = 10 ** (PRECISION_DIGITS - digits)
    scaled = ceil_div(value, factor)
    sign = "-" if scaled < 0 else ""
    scaled = abs(scaled)
    return f"{sign}{scaled // 10**digits}.{scaled % 10**digits:0{digits}d}"


POSITIVE_WITNESSES = {
    (2, 3): Fraction(16, 79), (2, 5): Fraction(6, 49),
    (2, 7): Fraction(1, 11), (2, 9): Fraction(3, 43),
    (2, 11): Fraction(3, 52), (2, 13): Fraction(2, 41),
    (2, 15): Fraction(3, 71), (2, 17): Fraction(3, 80),
    (2, 19): Fraction(1, 30), (2, 21): Fraction(1, 33),
    (2, 23): Fraction(1, 36), (2, 25): Fraction(1, 39),
    (2, 27): Fraction(1, 42), (2, 29): Fraction(1, 46),
    (3, 3): Fraction(4, 27), (3, 5): Fraction(7, 73),
    (3, 7): Fraction(1, 15), (3, 9): Fraction(1, 19),
    (3, 11): Fraction(3, 70), (5, 3): Fraction(7, 71),
    (5, 5): Fraction(1, 16), (7, 3): Fraction(1, 14),
}


NEXT_BRACKETS = {
    (2, 31): (Fraction(205498, 10**7), Fraction(205499, 10**7)),
    (3, 13): (Fraction(365378, 10**7), Fraction(365379, 10**7)),
    (5, 7): (Fraction(431772, 10**7), Fraction(431773, 10**7)),
    (7, 5): (Fraction(466407, 10**7), Fraction(466408, 10**7)),
}


def main() -> None:
    pi_i = pi_bounds()
    logs = {p: log_integer_bounds(p) for p in (2, 3, 5, 7)}
    positives: list[tuple[int, tuple[int, int]]] = []
    for case, y in POSITIVE_WITNESSES.items():
        interval = margin_interval(*case, y, pi_i, logs)
        require(interval[0] > 0, f"positive displayed witness {case}")
        positives.append((interval[0], case))
    smallest = min(positives)

    negative_uppers = {}
    for case, bracket in NEXT_BRACKETS.items():
        negative_uppers[case] = tangent_upper(*case, *bracket, pi_i, logs)

    print("THM-4089 HYBRID P-ADIC ZETA MARGIN FRONTIER CERTIFICATE")
    print(f"scale: 10^{PRECISION_DIGITS}")
    print(f"positive fixed witnesses: {len(positives)}/22")
    print(f"smallest positive case: {smallest[1]}")
    print(f"smallest positive lower: {decimal_lower(smallest[0])}")
    print("next cases: global maxima bounded by concave tangents")
    for case in sorted(negative_uppers):
        left, right = NEXT_BRACKETS[case]
        print(
            f"  {case}: Y in [{left},{right}], "
            f"global margin upper <= {decimal_upper(negative_uppers[case])}"
        )
    print("external 22-value irrationality theorem: NOT VERIFIED BY THIS CERTIFICATE")
    print(f"GATES: {GATES}")
    print("VERIFIED-EXACT: True")


if __name__ == "__main__":
    main()
