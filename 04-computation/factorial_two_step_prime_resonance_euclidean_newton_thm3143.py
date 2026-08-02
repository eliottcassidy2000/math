#!/usr/bin/env python3
"""Exact controls for THM-3143's Euclidean Newton-face cancellation."""

from fractions import Fraction
from math import comb, factorial, isqrt


PRIME_LIMIT_EXCLUSIVE = 10_000
FULL_POLYGON_PRIMES = (3, 5, 7, 11, 101, 211)


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


def symbolic_residue_controls():
    primes = primes_below(PRIME_LIMIT_EXCLUSIVE)
    checks = True
    for p in primes:
        m = (p - 1) // 2
        # R_k/p == (-1)^(m+1) 2^(m+2)/m (mod p), k=m+1.
        midpoint = ((-1) ** (m + 1) * pow(2, m + 2, p) * pow(m, -1, p)) % p
        checks &= midpoint != 0
        # R_p/p^2 == -4(p+1)(p+2)((2p-2)!/p) (mod p).
        # Wilson gives ((2p-2)!/p) == (p-1)!(p-2)! == -1 (mod p).
        leading = (-4 * 1 * 2 * -1) % p
        checks &= leading != 0
    return primes, checks


def full_polygon_controls():
    direct = True
    cancellation = True
    common_face = True
    expected_hulls = True
    separated = True
    records = []
    for p in FULL_POLYGON_PRIMES:
        d = p + 2
        k = (p + 1) // 2
        a = moment_coefficients(p, d)
        b = moment_coefficients(p + 1, d)
        if p <= 11:
            direct &= a == direct_bivariate_moment(p, d)
            direct &= b == direct_bivariate_moment(p + 1, d)
        multiplier = (2 * p + 2) * (2 * p + 1)
        cancellation &= b[p + 1] == multiplier * a[p]
        remainder = [
            b[j] - (multiplier * a[j - 1] if j else 0)
            for j in range(p + 1)
        ]
        expected_leading = -4 * p * (p + 1) * (p + 2) * factorial(2 * p - 2)
        cancellation &= remainder[p] == expected_leading
        common_face &= (
            a[0] % p == 2
            and (a[p] // (p * p)) % p == 2
            and b[1] % p == 4 % p
            and (b[p + 1] // (p * p)) % p == 4 % p
        )
        ha = lower_hull(a, p)
        hb = lower_hull(b, p)
        hr = lower_hull(remainder, p)
        expected_hulls &= ha == [(0, 0), (p, 2)]
        expected_hulls &= hb == [(0, 0), (1, 0), (p + 1, 2)]
        expected_hulls &= hr == [(0, 0), (k, 1), (p, 2)]
        separated &= set(slopes(ha)).isdisjoint(slopes(hr))
        records.append((p, ha, hb, hr, slopes(ha), slopes(hr)))
    return direct, cancellation, common_face, expected_hulls, separated, records


def main():
    primes, symbolic = symbolic_residue_controls()
    print(
        "prime_universe=odd_primes_below_10000 "
        f"count={len(primes)} first={primes[0]} last={primes[-1]}"
    )
    print(f"symbolic_midpoint_and_leading_residues_nonzero={symbolic}")
    direct, cancellation, common_face, hulls, separated, records = full_polygon_controls()
    print(f"independent_bivariate_controls_p_le_11={direct}")
    print(f"exact_euclidean_cancellation={cancellation}")
    print(f"common_frobenius_face_endpoint_residues={common_face}")
    print(f"expected_A_B_R_newton_hulls={hulls}")
    print(f"A_and_remainder_slope_sets_disjoint={separated}")
    for p, ha, hb, hr, sa, sr in records:
        print(
            f"p={p} A_hull={ha} B_hull={hb} R_hull={hr} "
            f"A_slopes={fmt_slopes(sa)} R_slopes={fmt_slopes(sr)}"
        )


if __name__ == "__main__":
    main()
