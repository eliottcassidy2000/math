#!/usr/bin/env python3
"""Exact controls for THM-3138's left-factorial resonance polygons."""

from fractions import Fraction
from hashlib import sha256
from itertools import product
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


def left_factorial_mod(p: int) -> int:
    """!p=sum_(k=0)^(p-1) k! modulo p."""
    running_factorial = 1
    answer = 1
    for k in range(1, p):
        running_factorial = running_factorial * k % p
        answer = (answer + running_factorial) % p
    return answer


def derangement_mod(n: int, modulus: int) -> int:
    """D_n modulo modulus, using D_n=n D_(n-1)+(-1)^n."""
    answer = 1
    for k in range(1, n + 1):
        answer = (k * answer + (-1 if k % 2 else 1)) % modulus
    return answer


def moment_coefficient(n: int, d: int, j: int) -> int:
    """[v^j] L((d-t+v*t^2)^n), with L(t^k)=k!."""
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
    """Independent repeated-product constructor for small controls."""
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


def lower_hull_from_valuations(values):
    """Lower hull of an explicit finite valuation word."""
    hull = []
    for point in enumerate(values):
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


def fmt_slopes(values):
    return "[" + ",".join(str(value) for value in values) + "]"


def residue_scan():
    primes = primes_below(PRIME_LIMIT_EXCLUSIVE)
    residues = []
    congruence = True
    for p in primes:
        residue = left_factorial_mod(p)
        residues.append((p, residue))
        congruence &= derangement_mod(p - 1, p) == residue
    semantic_bytes = "".join(f"{p}:{residue}\n" for p, residue in residues).encode()
    return primes, residues, congruence, sha256(semantic_bytes).hexdigest()


def full_polygon_controls():
    direct = True
    midpoint_units = True
    expected_hulls = True
    disjoint = True
    records = []
    for p in FULL_POLYGON_PRIMES:
        d = p + 1
        m = (p - 1) // 2
        first = moment_coefficients(p - 1, d)
        second = moment_coefficients(p, d)
        if p <= 11:
            direct &= first == direct_bivariate_moment(p - 1, d)
            direct &= second == direct_bivariate_moment(p, d)
        midpoint_units &= first[m] % p == (-1) ** (m + 1) % p
        first_hull = lower_hull(first, p)
        second_hull = lower_hull(second, p)
        expected_first = [(0, 0), (m, 0), (p - 1, 1)]
        expected_second = [(0, 0), (p, 2)]
        expected_hulls &= first_hull == expected_first
        expected_hulls &= second_hull == expected_second
        first_slopes = slopes(first_hull)
        second_slopes = slopes(second_hull)
        disjoint &= set(first_slopes).isdisjoint(second_slopes)
        records.append((p, first_hull, second_hull, first_slopes, second_slopes))
    return direct, midpoint_units, expected_hulls, disjoint, records


def synthetic_raised_endpoint_controls():
    """Exhaust the abstract valuation chamber used by the unconditional proof."""
    checked = 0
    valid = True
    for m in range(1, 6):
        target_slope = Fraction(2, 2 * m + 1)
        for left in product(range(4), repeat=m):
            for right in product(range(1, 4), repeat=m - 1):
                values = list(left) + [0] + list(right) + [1]
                current_slopes = slopes(lower_hull_from_valuations(values))
                valid &= current_slopes[-1] == Fraction(1, m)
                valid &= all(value <= 0 for value in current_slopes[:-1])
                valid &= target_slope not in current_slopes
                checked += 1
    return checked, valid


def main():
    primes, residues, congruence, digest = residue_scan()
    zeros = [p for p, residue in residues if residue == 0]
    new_starts = sum(p - 1 > 200 for p in primes)
    print(
        "prime_universe=odd_primes_below_10000 "
        f"count={len(primes)} first={primes[0]} last={primes[-1]}"
    )
    print(f"left_factorial_zero_residues={zeros}")
    print(f"derangement_left_factorial_congruence={congruence}")
    print(f"new_window_starts_r_gt_200={new_starts}")
    print(f"residue_semantic_sha256={digest}")
    print(
        "first_residues="
        + ",".join(f"{p}:{residue}" for p, residue in residues[:10])
    )
    print(
        "last_residues="
        + ",".join(f"{p}:{residue}" for p, residue in residues[-5:])
    )
    direct, midpoint, hulls, disjoint, records = full_polygon_controls()
    print(f"independent_bivariate_controls_p_le_11={direct}")
    print(f"midpoint_unit_formula={midpoint}")
    print(f"expected_newton_hulls={hulls}")
    print(f"slope_sets_disjoint={disjoint}")
    synthetic_count, synthetic_valid = synthetic_raised_endpoint_controls()
    print(
        f"synthetic_nonnegative_prefix_words={synthetic_count} "
        f"unconditional_slope_separation={synthetic_valid}"
    )
    for p, first_hull, second_hull, first_slopes, second_slopes in records:
        print(
            f"p={p} A_(p-1)_hull={first_hull} A_p_hull={second_hull} "
            f"slopes={fmt_slopes(first_slopes)}|{fmt_slopes(second_slopes)}"
        )


if __name__ == "__main__":
    main()
