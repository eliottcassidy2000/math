#!/usr/bin/env python3
"""Exact controls for THM-3146's full linear Euclidean cancellation."""

from fractions import Fraction
from math import comb, factorial, isqrt


PRIME_LIMIT_EXCLUSIVE = 10_000
FULL_POLYGON_PRIMES = (5, 7, 11, 29, 37, 41, 101, 211)


def primes_below(limit: int):
    sieve = bytearray(b"\x01") * limit
    if limit:
        sieve[0] = 0
    if limit > 1:
        sieve[1] = 0
    for p in range(2, isqrt(limit - 1) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit:p] = b"\x00" * (((limit - 1 - start) // p) + 1)
    return [p for p in range(3, limit, 2) if sieve[p]]


def vp(n: int, p: int):
    if n == 0:
        return None
    value = 0
    while n % p == 0:
        value += 1
        n //= p
    return value


def moment_coefficient(n: int, d: int, j: int) -> int:
    return comb(n, j) * sum(
        comb(n - j, ell)
        * d ** (n - j - ell)
        * (-1) ** ell
        * factorial(2 * j + ell)
        for ell in range(n - j + 1)
    )


def moment_coefficients(n: int, d: int):
    return [moment_coefficient(n, d, j) for j in range(n + 1)]


def direct_bivariate_moment(n: int, d: int):
    poly = {(0, 0): 1}
    base = {(0, 0): d, (1, 0): -1, (2, 1): 1}
    for _ in range(n):
        product = {}
        for (t0, v0), a in poly.items():
            for (t1, v1), b in base.items():
                key = (t0 + t1, v0 + v1)
                product[key] = product.get(key, 0) + a * b
        poly = product
    answer = [0] * (n + 1)
    for (t_degree, v_degree), coefficient in poly.items():
        answer[v_degree] += coefficient * factorial(t_degree)
    return answer


def lower_hull(coefficients, p: int):
    points = [(j, vp(a, p)) for j, a in enumerate(coefficients) if a]
    hull = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if Fraction(y1 - y0, x1 - x0) >= Fraction(y2 - y1, x2 - x1):
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def slopes(hull):
    return [
        Fraction(y1 - y0, x1 - x0)
        for (x0, y0), (x1, y1) in zip(hull, hull[1:])
    ]


def fmt_slopes(values):
    return "[" + ",".join(str(value) for value in values) + "]"


def cleared_remainder(a, b, p: int):
    """Clear the p-unit denominator in B-(c*v+q)A."""
    denominator = 2 * p + 1
    leading_ratio = (2 * p + 4) * (2 * p + 3)
    constant_numerator = 2 * (p + 2) * (p + 3)
    full = []
    for j in range(p + 3):
        value = denominator * b[j]
        if j:
            value -= denominator * leading_ratio * a[j - 1]
        if j < len(a):
            value += constant_numerator * a[j]
        full.append(value)
    assert full[p + 2] == 0 and full[p + 1] == 0
    return full[: p + 1]


def symbolic_residue_controls():
    primes = primes_below(PRIME_LIMIT_EXCLUSIVE)
    endpoint_exceptions = []
    leading_exceptions = []
    midpoint_units = True
    for p in primes:
        if p < 5:
            continue
        m = (p - 1) // 2
        endpoint = 87 % p
        leading = 37 % p
        if endpoint == 0:
            endpoint_exceptions.append(p)
        if leading == 0:
            leading_exceptions.append(p)
        midpoint = (
            (-1) ** m
            * 12
            * pow(3, m + 2, p)
            * pow(m * (m - 1), -1, p)
        ) % p
        midpoint_units &= midpoint != 0
    return primes, endpoint_exceptions, leading_exceptions, midpoint_units


def full_polygon_controls():
    direct = True
    quotient = True
    leading = True
    penultimate = True
    hulls = True
    records = []
    for p in FULL_POLYGON_PRIMES:
        d = p + 3
        k = (p + 1) // 2
        a = moment_coefficients(p + 1, d)
        b = moment_coefficients(p + 2, d)
        if p <= 7:
            direct &= a == direct_bivariate_moment(p + 1, d)
            direct &= b == direct_bivariate_moment(p + 2, d)
        remainder = cleared_remainder(a, b, p)
        expected_leading = (
            -2
            * p
            * (p + 1)
            * (p + 2)
            * (24 * p + 37)
            * factorial(2 * p - 2)
        )
        expected_penultimate = (
            2
            * p
            * (p - 1)
            * (p + 1)
            * (p + 2)
            * (24 * p * p - 115 * p - 246)
            * factorial(2 * p - 4)
        )
        leading &= remainder[p] == expected_leading
        penultimate &= remainder[p - 1] == expected_penultimate
        denominator = 2 * p + 1
        leading_ratio = (2 * p + 4) * (2 * p + 3)
        constant_numerator = 2 * (p + 2) * (p + 3)
        quotient &= (
            Fraction(-constant_numerator, denominator)
            == Fraction(
                b[p + 1] - leading_ratio * a[p], a[p + 1]
            )
        )
        ha = lower_hull(a, p)
        hr = lower_hull(remainder, p)
        expected_a = [(0, 0), (1, 0), (p + 1, 2)]
        if p == 29:
            expected_r = [(0, 1), (k, 1), (p, 2)]
        elif p == 37:
            expected_r = [(0, 0), (k, 1), (p - 1, 2), (p, 3)]
        else:
            expected_r = [(0, 0), (k, 1), (p, 2)]
        hulls &= ha == expected_a and hr == expected_r
        records.append((p, ha, hr, slopes(ha), slopes(hr)))
    return direct, quotient, leading, penultimate, hulls, records


def main():
    primes, endpoint, leading_ex, midpoint = symbolic_residue_controls()
    print(
        "prime_universe=odd_primes_below_10000 "
        f"count={len(primes)} first={primes[0]} last={primes[-1]}"
    )
    print(f"constant_residue_87_exceptions={endpoint}")
    print(f"leading_residue_37_exceptions={leading_ex}")
    print(f"symbolic_midpoint_residue_nonzero={midpoint}")
    direct, quotient, leading, penultimate, hulls, records = full_polygon_controls()
    print(f"independent_bivariate_controls_p_le_7={direct}")
    print(f"full_linear_quotient_identity={quotient}")
    print(f"exact_leading_remainder_identity={leading}")
    print(f"exact_penultimate_remainder_identity={penultimate}")
    print(f"expected_regular_and_exceptional_hulls={hulls}")
    for p, ha, hr, sa, sr in records:
        print(
            f"p={p} A_hull={ha} R_hull={hr} "
            f"A_slopes={fmt_slopes(sa)} R_slopes={fmt_slopes(sr)} "
            f"overlap={fmt_slopes(sorted(set(sa).intersection(sr)))}"
        )


if __name__ == "__main__":
    main()
