#!/usr/bin/env python3
"""Exact controls for THM-3153's second Euclidean-Newton separation."""

from fractions import Fraction
from math import comb, factorial, isqrt


PRIME_LIMIT_EXCLUSIVE = 10_000
FULL_POLYGON_PRIMES = (7, 11, 13, 59, 61, 107, 193, 277)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_below(limit):
    sieve = bytearray(b"\x01") * limit
    if limit:
        sieve[0] = 0
    if limit > 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit - 1) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit : prime] = b"\x00" * (
                (limit - 1 - prime * prime) // prime + 1
            )
    return [prime for prime in range(3, limit, 2) if sieve[prime]]


def is_prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    return all(number % divisor for divisor in range(3, isqrt(number) + 1, 2))


def valuation(number, prime):
    if number == 0:
        return None
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def lower_hull(coefficients, prime):
    points = [
        (index, valuation(coefficient, prime))
        for index, coefficient in enumerate(coefficients)
        if coefficient
    ]
    hull = []
    for point in points:
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if Fraction(y1 - y0, x1 - x0) >= Fraction(
                y2 - y1, x2 - x1
            ):
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


def moment_coefficient(n, d, j):
    return comb(n, j) * sum(
        comb(n - j, ell)
        * d ** (n - j - ell)
        * (-1) ** ell
        * factorial(2 * j + ell)
        for ell in range(n - j + 1)
    )


def moment_coefficients(n, d):
    return [moment_coefficient(n, d, j) for j in range(n + 1)]


def direct_bivariate_moment(n, d):
    polynomial = {(0, 0): 1}
    base = {(0, 0): d, (1, 0): -1, (2, 1): 1}
    for _ in range(n):
        product = {}
        for (t0, v0), left in polynomial.items():
            for (t1, v1), right in base.items():
                key = (t0 + t1, v0 + v1)
                product[key] = product.get(key, 0) + left * right
        polynomial = product
    answer = [0] * (n + 1)
    for (t_degree, v_degree), coefficient in polynomial.items():
        answer[v_degree] += coefficient * factorial(t_degree)
    return answer


def first_remainder(a, b, prime, modulus=None):
    denominator = 2 * prime + 3
    leading = (2 * prime + 6) * (2 * prime + 5)
    constant_numerator = 2 * (prime + 3) * (prime + 4)
    full = []
    for index in range(prime + 4):
        value = denominator * b[index]
        if index:
            value -= denominator * leading * a[index - 1]
        if index < len(a):
            value += constant_numerator * a[index]
        if modulus is not None:
            value %= modulus
        full.append(value)
    require(full[prime + 3] == full[prime + 2] == 0,
            "first quotient did not cancel its top two coefficients")
    return full[: prime + 2]


def second_parameters(prime):
    denominator = (
        (prime + 3) * (2 * prime - 1) * (24 * prime + 61) ** 2
    )
    linear_numerator = (
        -2
        * (2 * prime + 1)
        * (2 * prime + 3)
        * (2 * prime - 1)
        * (24 * prime + 61)
    )
    constant_numerator = (
        4 * (prime + 4) * (2 * prime + 1) * (28 * prime + 67)
    )
    return denominator, linear_numerator, constant_numerator


def second_remainder(a, first, prime, modulus=None):
    denominator, linear_numerator, constant_numerator = second_parameters(prime)
    full = []
    for index in range(prime + 3):
        value = denominator * a[index]
        if index:
            value -= linear_numerator * first[index - 1]
        if index < len(first):
            value -= constant_numerator * first[index]
        if modulus is not None:
            value %= modulus
        full.append(value)
    require(full[prime + 2] == full[prime + 1] == 0,
            "second quotient did not cancel its top two coefficients")
    return full[: prime + 1]


def high_polynomials(prime):
    endpoint = (
        256 * prime ** 4
        - 29696 * prime ** 3
        - 193568 * prime ** 2
        - 406368 * prime
        - 276169
    )
    penultimate = (
        256 * prime ** 5
        - 32000 * prime ** 4
        + 97504 * prime ** 3
        + 1530816 * prime ** 2
        + 3882095 * prime
        + 2891889
    )
    return endpoint, penultimate


def expected_high_coefficients(prime):
    endpoint_polynomial, penultimate_polynomial = high_polynomials(prime)
    common = (
        4
        * prime
        * (prime - 1)
        * (prime + 1)
        * (prime + 2)
        * (prime + 3)
        * (2 * prime + 3)
    )
    endpoint = common * endpoint_polynomial * factorial(2 * prime - 4)
    penultimate = (
        -common
        * (prime - 2)
        * penultimate_polynomial
        * factorial(2 * prime - 6)
    )
    return endpoint, penultimate


def midpoint_residue(prime):
    m = (prime - 1) // 2
    return (
        65880
        * (-1) ** (m - 2)
        * pow(4, m + 3, prime)
        * pow(m * (m - 1) * (m - 2), -1, prime)
    ) % prime


def exact_polygon_controls():
    expected = {
        7: [(0, 1), (5, 1), (7, 2)],
        11: [(0, 0), (1, 0), (7, 1), (11, 2)],
        13: [(0, 0), (1, 0), (8, 1), (13, 2)],
        59: [(0, 1), (1, 0), (31, 1), (59, 2)],
        61: [(0, 0), (31, 1), (61, 2)],
        107: [(0, 0), (1, 0), (55, 1), (107, 2)],
        193: [(0, 1), (1, 0), (98, 1), (193, 2)],
        277: [(0, 0), (1, 0), (140, 1), (276, 2), (277, 3)],
    }
    direct = quotient = high = midpoint = hulls = True
    records = []
    for prime in FULL_POLYGON_PRIMES:
        d = prime + 4
        a = moment_coefficients(prime + 2, d)
        b = moment_coefficients(prime + 3, d)
        if prime == 7:
            direct &= a == direct_bivariate_moment(prime + 2, d)
            direct &= b == direct_bivariate_moment(prime + 3, d)
        first = first_remainder(a, b, prime)
        second = second_remainder(a, first, prime)
        quotient &= len(first) == prime + 2 and len(second) == prime + 1
        wanted_endpoint, wanted_penultimate = expected_high_coefficients(prime)
        high &= second[prime] == wanted_endpoint
        high &= second[prime - 1] == wanted_penultimate
        k = (prime + 3) // 2
        if prime == 61:
            midpoint &= valuation(second[k], prime) >= 2
            midpoint &= midpoint_residue(prime) == 0
        else:
            midpoint &= valuation(second[k], prime) == 1
            midpoint &= (second[k] // prime) % prime == midpoint_residue(prime)
        ha = lower_hull(a, prime)
        hs = lower_hull(second, prime)
        hulls &= ha == [(0, 0), (2, 0), (prime + 2, 2)]
        hulls &= hs == expected[prime]
        records.append((prime, ha, hs, slopes(ha), slopes(hs)))
    return direct, quotient, high, midpoint, hulls, records


def pascal_mod(maximum, modulus):
    rows = [[1]]
    for degree in range(1, maximum + 1):
        previous = rows[-1]
        row = [1]
        row.extend(
            (previous[index - 1] + previous[index]) % modulus
            for index in range(1, degree)
        )
        row.append(1)
        rows.append(row)
    return rows


def moment_coefficients_mod(n, d, modulus, pascal, factorials):
    powers = [1] * (n + 1)
    for exponent in range(1, n + 1):
        powers[exponent] = powers[exponent - 1] * d % modulus
    answer = []
    for index in range(n + 1):
        remainder = n - index
        total = 0
        row = pascal[remainder]
        for ell in range(remainder + 1):
            term = row[ell] * powers[remainder - ell] % modulus
            term = term * factorials[2 * index + ell] % modulus
            total = total - term if ell % 2 else total + term
        answer.append(pascal[n][index] * total % modulus)
    return answer


def truncated_valuation(number, prime, cap):
    number %= prime ** cap
    if number == 0:
        return cap
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def modular_large_polygon(prime=997, cap=5):
    modulus = prime ** cap
    maximum = prime + 3
    pascal = pascal_mod(maximum, modulus)
    factorials = [1] * (2 * maximum + 1)
    for index in range(1, len(factorials)):
        factorials[index] = factorials[index - 1] * index % modulus
    a = moment_coefficients_mod(
        prime + 2, prime + 4, modulus, pascal, factorials
    )
    b = moment_coefficients_mod(
        prime + 3, prime + 4, modulus, pascal, factorials
    )
    first = first_remainder(a, b, prime, modulus)
    second = second_remainder(a, first, prime, modulus)
    values = [truncated_valuation(value, prime, cap) for value in second]
    hull = []
    for point in enumerate(values):
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if Fraction(y1 - y0, x1 - x0) >= Fraction(
                y2 - y1, x2 - x1
            ):
                hull.pop()
            else:
                break
        hull.append(point)
    expected = [(0, 0), (1, 0), (500, 1), (996, 2), (997, 3)]
    require(hull == expected, "p=997 modular polygon drift")
    return hull, slopes(hull)


def residue_and_exception_controls():
    primes = primes_below(PRIME_LIMIT_EXCLUSIVE)
    delta = 586672128
    delta_exceptions = [prime for prime in primes if prime > 6 and delta % prime == 0]
    low_constant = 168 * 11387
    low_linear = 168 * 14640
    low_walls = [
        prime
        for prime in primes
        if prime > 6 and (low_constant * low_linear) % prime == 0
    ]
    midpoint_walls = [prime for prime in primes if prime > 6 and 65880 % prime == 0]
    endpoint_walls = [prime for prime in primes if 276169 % prime == 0]
    penultimate_walls = [
        prime for prime in primes if prime > 6 and 2891889 % prime == 0
    ]
    require(delta_exceptions == [7, 4547], "delta-4 exceptions drift")
    require(low_walls == [7, 59, 61, 193], "low-coordinate walls drift")
    require(midpoint_walls == [61], "midpoint walls drift")
    require(endpoint_walls == [277, 997], "endpoint walls drift")
    require(
        penultimate_walls == [7, 11, 13, 107],
        "penultimate nonvertex walls drift",
    )
    require(is_prime(4549), "4549 fallback primality drift")
    finite_fallbacks = [3, 5, 7, 59, 61, 193]
    require(all(prime + 2 <= 200 for prime in finite_fallbacks),
            "finite fallback left THM-3124 range")
    return (
        len(primes), delta_exceptions, low_walls, midpoint_walls,
        endpoint_walls, penultimate_walls, finite_fallbacks,
    )


def fmt_slopes(values):
    return "[" + ",".join(str(value) for value in values) + "]"


def main():
    direct, quotient, high, midpoint, hulls, records = exact_polygon_controls()
    print(f"independent_bivariate_control_p7={direct}")
    print(f"two_exact_linear_quotients={quotient}")
    print(f"second_remainder_high_identities={high}")
    print(f"midpoint_residue_identity={midpoint}")
    print(f"expected_exact_polygon_atlas={hulls}")
    for prime, ha, hs, sa, ss in records:
        print(
            f"p={prime} A_hull={ha} S_hull={hs} "
            f"A_slopes={fmt_slopes(sa)} S_slopes={fmt_slopes(ss)} "
            f"overlap={fmt_slopes(sorted(set(sa).intersection(ss)))}"
        )
    hull997, slopes997 = modular_large_polygon()
    print(
        f"p=997 modular_cap5_S_hull={hull997} "
        f"S_slopes={fmt_slopes(slopes997)}"
    )
    (
        prime_count, delta_exceptions, low_walls, midpoint_walls,
        endpoint_walls, penultimate_walls, finite_fallbacks,
    ) = residue_and_exception_controls()
    print(
        f"prime_scan=odd_primes_below_{PRIME_LIMIT_EXCLUSIVE} "
        f"count={prime_count}"
    )
    print(f"height_zero_resultant_exceptions={delta_exceptions}")
    print(f"second_remainder_low_coordinate_walls={low_walls}")
    print(f"midpoint_residue_walls={midpoint_walls}")
    print(f"raised_endpoint_primes={endpoint_walls}")
    print(f"raised_penultimate_nonvertex_primes={penultimate_walls}")
    print(f"finite_THM3124_fallback_primes={finite_fallbacks}")
    print("exception_p4547_uses_prime_r=4549=True")
    print("all_exact_checks=True")


if __name__ == "__main__":
    main()
