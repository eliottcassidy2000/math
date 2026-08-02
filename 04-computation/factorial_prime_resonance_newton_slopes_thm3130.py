#!/usr/bin/env python3
"""Exact controls for THM-3130's prime-resonance Newton polygons."""

from fractions import Fraction
from math import comb, factorial, isqrt


PRIME_LIMIT = 223


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for d in range(3, isqrt(n) + 1, 2):
        if n % d == 0:
            return False
    return True


def vp(n: int, p: int):
    if n == 0:
        return None
    answer = 0
    while n % p == 0:
        answer += 1
        n //= p
    return answer


def moment_coefficient(n: int, d: int, j: int) -> int:
    """[v^j] L((d-t+v*t^2)^n), where L(t^k)=k!."""
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
    """Independent repeated-product construction for small controls."""
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
    points = [(j, vp(a, p)) for j, a in enumerate(coefficients) if a != 0]
    hull = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            old_slope = Fraction(y1 - y0, x1 - x0)
            new_slope = Fraction(y2 - y1, x2 - x1)
            if old_slope >= new_slope:
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


def prime_controls():
    primes = [p for p in range(3, PRIME_LIMIT + 1) if is_prime(p)]
    formula_direct = True
    reduction = True
    hulls = True
    disjoint = True
    for p in primes:
        a = moment_coefficients(p - 2, p)
        b = moment_coefficients(p - 1, p)
        if p <= 13:
            formula_direct &= a == direct_bivariate_moment(p - 2, p)
            formula_direct &= b == direct_bivariate_moment(p - 1, p)
        for n, coefficients in ((p - 2, a), (p - 1, b)):
            for j, coefficient in enumerate(coefficients):
                terminal = comb(n, j) * (-1) ** (n - j) * factorial(n + j)
                reduction &= (coefficient - terminal) % p == 0
        ha = lower_hull(a, p)
        hb = lower_hull(b, p)
        expected_a = [(0, 0), (1, 0)] if p == 3 else [(0, 0), (1, 0), (p - 2, 1)]
        expected_b = [(0, 0), (p - 1, 1)]
        hulls &= ha == expected_a and hb == expected_b
        disjoint &= set(slopes(ha)).isdisjoint(slopes(hb))
    return primes, formula_direct, reduction, hulls, disjoint


def hostile_composite(d: int, p: int):
    a = moment_coefficients(d - 2, d)
    b = moment_coefficients(d - 1, d)
    ha = lower_hull(a, p)
    hb = lower_hull(b, p)
    sa = slopes(ha)
    sb = slopes(hb)
    return ha, hb, sa, sb, sorted(set(sa).intersection(sb))


def main():
    primes, direct, reduction, hulls, disjoint = prime_controls()
    print(
        f"prime_universe=odd_primes_3_through_{PRIME_LIMIT} count={len(primes)} "
        f"first={primes[0]} last={primes[-1]}"
    )
    print(f"independent_bivariate_controls_p_le_13={direct}")
    print(f"terminal_term_reduction_mod_p={reduction}")
    print(f"expected_newton_hulls={hulls}")
    print(f"prime_slope_sets_disjoint={disjoint}")
    for d, p in ((4, 2), (6, 3)):
        ha, hb, sa, sb, overlap = hostile_composite(d, p)
        print(
            f"hostile_composite_d={d} p={p} A_hull={ha} B_hull={hb} "
            f"A_slopes={fmt_slopes(sa)} B_slopes={fmt_slopes(sb)} "
            f"overlap={fmt_slopes(overlap)}"
        )


if __name__ == "__main__":
    main()
