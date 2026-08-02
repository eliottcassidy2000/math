#!/usr/bin/env python3
"""Exact controls for THM-3148's fixed-offset Frobenius endpoint descent."""

from fractions import Fraction
from math import comb, factorial, isqrt


PRIME_LIMIT = 97
OFFSET_LIMIT = 8
TAIL_LIMIT = 6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def primes_up_to(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [p for p in range(2, limit + 1) if sieve[p]]


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


def fixed_endpoint(a, s):
    return moment_coefficients(a, s)


def normalized_endpoint(a, s):
    """Coefficients of L((1-t/s+u*t^2)^a) in ascending u degree."""
    return [
        comb(a, j)
        * sum(
            Fraction(
                comb(a - j, ell) * (-1) ** ell * factorial(2 * j + ell),
                s ** ell,
            )
            for ell in range(a - j + 1)
        )
        for j in range(a + 1)
    ]


def bareiss_determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for entry in range(column + 1, size):
                numerator = (
                    work[row][entry] * pivot
                    - work[row][column] * work[column][entry]
                )
                require(numerator % previous == 0, "Bareiss division lost exactness")
                work[row][entry] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def resultant(left, right):
    m = len(left) - 1
    n = len(right) - 1
    if m < 0 or n < 0:
        raise ValueError("zero polynomial has no resultant here")
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    matrix = []
    for shift in range(n):
        matrix.append(
            [0] * shift + left_descending + [0] * (n - 1 - shift)
        )
    for shift in range(m):
        matrix.append(
            [0] * shift + right_descending + [0] * (m - 1 - shift)
        )
    return bareiss_determinant(matrix)


def factor_integer(number):
    number = abs(number)
    factors = {}
    divisor = 2
    while divisor * divisor <= number:
        while number % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            number //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if number > 1:
        factors[number] = factors.get(number, 0) + 1
    return factors


def format_factorization(factors):
    if not factors:
        return "1"
    return "*".join(
        str(prime) if exponent == 1 else f"{prime}^{exponent}"
        for prime, exponent in sorted(factors.items())
    )


def trim_mod(polynomial, prime):
    answer = [coefficient % prime for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def remainder_mod(dividend, divisor, prime):
    dividend = trim_mod(dividend, prime)
    divisor = trim_mod(divisor, prime)
    require(divisor != [0], "division by the zero polynomial")
    while dividend != [0] and len(dividend) >= len(divisor):
        quotient = dividend[-1] * pow(divisor[-1], -1, prime) % prime
        shift = len(dividend) - len(divisor)
        for index, coefficient in enumerate(divisor):
            dividend[index + shift] = (
                dividend[index + shift] - quotient * coefficient
            ) % prime
        dividend = trim_mod(dividend, prime)
    return dividend


def gcd_mod(left, right, prime):
    left = trim_mod(left, prime)
    right = trim_mod(right, prime)
    while right != [0]:
        left, right = right, remainder_mod(left, right, prime)
    inverse = pow(left[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in left]


def congruence_controls():
    primes = primes_up_to(PRIME_LIMIT)
    cases = 0
    direct = True
    untruncated = True
    for prime in primes:
        for s in range(1, OFFSET_LIMIT + 1):
            for a in range(TAIL_LIMIT + 1):
                large = moment_coefficients(prime + a, prime + s)
                small = fixed_endpoint(a, s)
                expected = [
                    (s * coefficient) % prime for coefficient in small
                ] + [0] * prime
                require(
                    [coefficient % prime for coefficient in large] == expected,
                    "fixed-offset Frobenius congruence drift",
                )
                if prime > 2 * a and s % prime:
                    require(
                        small[-1] % prime == factorial(2 * a) % prime != 0,
                        "untruncated endpoint degree drift",
                    )
                cases += 1
                if prime <= 7 and s <= 4 and a <= 3:
                    direct &= large == direct_bivariate_moment(
                        prime + a, prime + s
                    )
    return len(primes), cases, direct, untruncated


def scaling_and_constant_controls():
    scaling = True
    recurrence = True
    checks = 0
    for s in range(1, 21):
        endpoints = [fixed_endpoint(a, s) for a in range(22)]
        for a in range(21):
            fixed = endpoints[a]
            normalized = normalized_endpoint(a, s)
            scaling &= all(
                fixed[j] == s ** (a - j) * normalized[j]
                for j in range(a + 1)
            )
            recurrence &= (
                endpoints[a + 1][0]
                == s ** (a + 1) - (a + 1) * fixed[0]
            )
            checks += 1
    return scaling, recurrence, checks


def resultant_controls():
    expected = {
        2: {},
        3: {2: 2, 29: 1},
        4: {2: 11, 3: 2, 7: 1, 4547: 1},
        5: {2: 16, 3: 7, 5: 2, 11: 2, 20747: 1, 249721: 1},
        6: {
            2: 47,
            3: 15,
            5: 4,
            7: 1,
            139: 1,
            3767: 1,
            12041: 1,
            807241: 1,
        },
    }
    records = []
    for s, wanted in expected.items():
        left = fixed_endpoint(s - 2, s)
        right = fixed_endpoint(s - 1, s)
        delta = resultant(left, right)
        factors = factor_integer(delta)
        require(factors == wanted, f"offset-{s} resultant factorization drift")
        endpoint_delta = s ** (2 * s - 3) * delta
        require(
            endpoint_delta
            == resultant(
                [s * coefficient for coefficient in left],
                [s * coefficient for coefficient in right],
            ),
            "outer-s resultant scaling drift",
        )
        exceptional = tuple(
            prime for prime in factors if prime > 2 * (s - 1)
        )
        records.append((s, delta, factors, exceptional))
    return records


def modular_exception_controls():
    expected = {
        (3, 29): [1, 1],
        (4, 7): [2, 1],
        (4, 4547): [1304, 1],
        (5, 11): [9, 1],
        (5, 20747): [15102, 1],
        (5, 249721): [26953, 1],
        (6, 139): [26, 1],
        (6, 3767): [2126, 1],
        (6, 12041): [3037, 1],
        (6, 807241): [784752, 1],
    }
    records = []
    for (s, prime), wanted in expected.items():
        left = fixed_endpoint(s - 2, s)
        right = fixed_endpoint(s - 1, s)
        common = gcd_mod(left, right, prime)
        require(common == wanted, "exceptional residual gcd drift")
        root = (-common[0]) % prime
        require(root != 0, "exceptional residual root became zero")
        records.append((s, prime, common, root))
    wall = gcd_mod(fixed_endpoint(2, 4), fixed_endpoint(3, 4), 61)
    require(wall == [1], "p=61 quotient wall became a residual exception")
    return records, wall


def main():
    prime_count, cases, direct, untruncated = congruence_controls()
    print(
        f"frobenius_universe=primes_le_{PRIME_LIMIT},"
        f"s_1_through_{OFFSET_LIMIT},a_0_through_{TAIL_LIMIT} "
        f"prime_count={prime_count} cases={cases}"
    )
    print(f"coefficientwise_frobenius_congruence=True")
    print(f"independent_small_bivariate_expansion={direct}")
    print(f"untruncated_leading_degree_controls={untruncated}")
    scaling, recurrence, checks = scaling_and_constant_controls()
    print(f"base_window_scaling_identity={scaling} checks={checks}")
    print(f"zero_root_exclusion_recurrence={recurrence} checks={checks}")
    records = resultant_controls()
    for s, delta, factors, exceptional in records:
        print(
            f"s={s} base_r={s-2} delta={delta} "
            f"factorization={format_factorization(factors)} "
            f"eligible_exceptional_primes={list(exceptional)}"
        )
    exception_records, wall = modular_exception_controls()
    for s, prime, common, root in exception_records:
        print(
            f"exception_s={s} p={prime} monic_gcd={common} "
            f"nonzero_root={root}"
        )
    print(f"quotient_denominator_wall_s=4_p=61_monic_gcd={wall}")
    print("all_exact_checks=True")


if __name__ == "__main__":
    main()
