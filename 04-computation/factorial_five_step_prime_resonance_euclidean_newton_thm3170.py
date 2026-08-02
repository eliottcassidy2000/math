#!/usr/bin/env python3
"""Exact controls for THM-3170's offset-five Euclidean--Newton flag.

The companion fixes the normalization

  d=p+5, A=A_(p+3), B=A_(p+4),
  R=(2p+5)B-[(2p+5)(2p+8)(2p+7)v-2(p+4)(p+5)]A,
  S=D A-(N1 v+N0)R,

where D,N1,N0 are returned by second_parameters().
"""

from fractions import Fraction
from math import comb, factorial, isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    return all(number % divisor for divisor in range(3, isqrt(number) + 1, 2))


def primes_up_to(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [number for number in range(2, limit + 1) if sieve[number]]


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


def bareiss_determinant(matrix):
    size = len(matrix)
    if not size:
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
                require(numerator % previous == 0,
                        "Bareiss division lost exactness")
                work[row][entry] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def resultant(left, right):
    m = len(left) - 1
    n = len(right) - 1
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    matrix = []
    for shift in range(n):
        matrix.append([0] * shift + left_descending + [0] * (n - 1 - shift))
    for shift in range(m):
        matrix.append([0] * shift + right_descending + [0] * (m - 1 - shift))
    return bareiss_determinant(matrix)


def trim_mod(polynomial, prime):
    answer = [coefficient % prime for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def remainder_mod(dividend, divisor, prime):
    dividend = trim_mod(dividend, prime)
    divisor = trim_mod(divisor, prime)
    require(divisor != [0], "division by zero polynomial")
    while dividend != [0] and len(dividend) >= len(divisor):
        scale = dividend[-1] * pow(divisor[-1], -1, prime) % prime
        shift = len(dividend) - len(divisor)
        for index, coefficient in enumerate(divisor):
            dividend[index + shift] = (
                dividend[index + shift] - scale * coefficient
            ) % prime
        dividend = trim_mod(dividend, prime)
    return dividend


def gcd_mod(left, right, prime):
    left = trim_mod(left, prime)
    right = trim_mod(right, prime)
    while right != [0]:
        left, right = right, remainder_mod(left, right, prime)
    scale = pow(left[-1], -1, prime)
    return [(coefficient * scale) % prime for coefficient in left]


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


def first_remainder(a, b, prime, modulus=None):
    denominator = 2 * prime + 5
    leading = (2 * prime + 8) * (2 * prime + 7)
    constant_numerator = 2 * (prime + 4) * (prime + 5)
    full = []
    for index in range(prime + 5):
        value = denominator * b[index]
        if index:
            value -= denominator * leading * a[index - 1]
        if index < len(a):
            value += constant_numerator * a[index]
        if modulus is not None:
            value %= modulus
        full.append(value)
    require(full[prime + 4] == full[prime + 3] == 0,
            "first quotient did not cancel its top two coefficients")
    return full[: prime + 3]


def second_parameters(prime):
    denominator = (prime + 4) * (2 * prime + 1) * (24 * prime + 85) ** 2
    linear_numerator = (
        -2
        * (2 * prime + 3)
        * (2 * prime + 5)
        * (2 * prime + 1)
        * (24 * prime + 85)
    )
    constant_numerator = (
        4 * (prime + 5) * (2 * prime + 3) * (28 * prime + 95)
    )
    return denominator, linear_numerator, constant_numerator


def second_remainder(a, first, prime, modulus=None):
    denominator, linear_numerator, constant_numerator = second_parameters(prime)
    full = []
    for index in range(prime + 4):
        value = denominator * a[index]
        if index:
            value -= linear_numerator * first[index - 1]
        if index < len(first):
            value -= constant_numerator * first[index]
        if modulus is not None:
            value %= modulus
        full.append(value)
    require(full[prime + 3] == full[prime + 2] == 0,
            "second quotient did not cancel its top two coefficients")
    return full[: prime + 2]


def endpoint_polynomials(prime):
    endpoint = (
        256 * prime**4
        - 28672 * prime**3
        - 281120 * prime**2
        - 881568 * prime
        - 905545
    )
    penultimate = (
        256 * prime**5
        - 30720 * prime**4
        - 27936 * prime**3
        + 1633888 * prime**2
        + 7109519 * prime
        + 8370560
    )
    return endpoint, penultimate


def expected_high_coefficients(prime):
    endpoint_polynomial, penultimate_polynomial = endpoint_polynomials(prime)
    common = (
        4
        * prime
        * (prime + 1)
        * (prime + 2)
        * (prime + 3)
        * (prime + 4)
        * (2 * prime + 5)
    )
    endpoint = common * endpoint_polynomial * factorial(2 * prime - 2)
    penultimate = (
        -common
        * (prime - 1)
        * penultimate_polynomial
        * factorial(2 * prime - 4)
    )
    return endpoint, penultimate


def valuation(number, prime, cap=None):
    if cap is not None:
        number %= prime**cap
    if number == 0:
        return cap if cap is not None else None
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def lower_hull_from_values(values):
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


def slopes(hull):
    return [
        Fraction(y1 - y0, x1 - x0)
        for (x0, y0), (x1, y1) in zip(hull, hull[1:])
    ]


def midpoint_residue(prime):
    m = (prime - 1) // 2
    return (
        -4284000
        * (-1) ** (m - 3)
        * pow(5, m + 4, prime)
        * pow(m * (m - 1) * (m - 2) * (m - 3), -1, prime)
    ) % prime


def targeted_moment_mod(n, d, j, modulus, factorials):
    return comb(n, j) * sum(
        comb(n - j, ell)
        * pow(d, n - j - ell, modulus)
        * (-1) ** ell
        * factorials[2 * j + ell]
        for ell in range(n - j + 1)
    ) % modulus


def targeted_midpoint_control(prime):
    require(prime >= 11 and is_prime(prime), "midpoint domain drift")
    modulus = prime**3
    n = prime + 3
    d = prime + 5
    k = (prime + 3) // 2
    factorials = [1]
    for index in range(1, 2 * (prime + 4) + 1):
        factorials.append(factorials[-1] * index % modulus)

    def aa(index):
        return targeted_moment_mod(n, d, index, modulus, factorials)

    def bb(index):
        return targeted_moment_mod(n + 1, d, index, modulus, factorials)

    denominator = 2 * prime + 5
    leading = (2 * prime + 8) * (2 * prime + 7)
    constant_numerator = 2 * (prime + 4) * (prime + 5)

    def rr(index):
        value = denominator * bb(index)
        if index:
            value -= denominator * leading * aa(index - 1)
        value += constant_numerator * aa(index)
        return value % modulus

    den2, num1, num0 = second_parameters(prime)
    sk = (den2 * aa(k) - num1 * rr(k - 1) - num0 * rr(k)) % modulus
    require(sk % prime == 0, "midpoint lost first p factor")
    return valuation(sk, prime, 3), (sk // prime) % prime, midpoint_residue(prime)


def moment_coefficients_mod_stream(n, d, modulus):
    powers = [1]
    for _ in range(n):
        powers.append(powers[-1] * d % modulus)
    factorials = [1]
    for index in range(1, 2 * n + 1):
        factorials.append(factorials[-1] * index % modulus)
    row = [1]
    answer = []
    # h=n-j runs upward, so append coefficients in descending j and reverse.
    for h in range(n + 1):
        j = n - h
        total = 0
        for ell, choose in enumerate(row):
            term = choose * powers[h - ell] % modulus
            term = term * factorials[2 * j + ell] % modulus
            total = total - term if ell % 2 else total + term
        answer.append(comb(n, h) * total % modulus)
        if h < n:
            row = [1] + [
                (row[index - 1] + row[index]) % modulus
                for index in range(1, len(row))
            ] + [1]
    answer.reverse()
    return answer


def modular_polygon(prime, cap=5):
    modulus = prime**cap
    a = moment_coefficients_mod_stream(prime + 3, prime + 5, modulus)
    b = moment_coefficients_mod_stream(prime + 4, prime + 5, modulus)
    first = first_remainder(a, b, prime, modulus)
    second = second_remainder(a, first, prime, modulus)
    expected_endpoint, expected_penultimate = expected_high_coefficients(prime)
    require(second[prime + 1] == expected_endpoint % modulus,
            "modular endpoint identity drift")
    require(second[prime] == expected_penultimate % modulus,
            "modular penultimate identity drift")
    return (
        lower_hull_from_values([valuation(x, prime, cap) for x in a]),
        lower_hull_from_values([valuation(x, prime, cap) for x in first]),
        lower_hull_from_values([valuation(x, prime, cap) for x in second]),
    )


def exact_small_atlas():
    records = []
    expected = {
        11: [(0, 0), (1, 0), (7, 1), (12, 2)],
        17: [(0, 0), (1, 0), (18, 2)],
        19: [(0, 0), (1, 0), (11, 1), (20, 2)],
        29: [(0, 0), (1, 0), (16, 1), (30, 2)],
        41: [(0, 0), (1, 0), (22, 1), (42, 2)],
        61: [(0, 0), (1, 0), (32, 1), (61, 2), (62, 3)],
    }
    for prime, wanted in expected.items():
        a = moment_coefficients(prime + 3, prime + 5)
        b = moment_coefficients(prime + 4, prime + 5)
        first = first_remainder(a, b, prime)
        second = second_remainder(a, first, prime)
        if prime in (11, 17):
            modulus = prime**5
            require(
                moment_coefficients_mod_stream(prime + 3, prime + 5, modulus)
                == [coefficient % modulus for coefficient in a],
                "streaming A implementation drift",
            )
            require(
                moment_coefficients_mod_stream(prime + 4, prime + 5, modulus)
                == [coefficient % modulus for coefficient in b],
                "streaming B implementation drift",
            )
        high = expected_high_coefficients(prime)
        require((second[prime + 1], second[prime]) == high,
                "high coefficient identity drift")
        s_hull = lower_hull_from_values([valuation(x, prime) for x in second])
        require(s_hull == wanted, "small exact S polygon drift")
        mid = targeted_midpoint_control(prime)
        require(mid[1] == mid[2], "midpoint residue drift")
        records.append((prime, s_hull, slopes(s_hull), mid))
    return records


def main():
    print("THM-3170 five-step prime resonance Euclidean--Newton separation")
    print("offset_five_normalization=d=p+5,A=A_(p+3),B=A_(p+4)")
    print(
        "first_quotient=Cv+q C=(2p+8)(2p+7) "
        "q=-2(p+4)(p+5)/(2p+5)"
    )
    print(
        "second_quotient_q1=-2(2p+3)(2p+5)/((p+4)(24p+85))"
    )
    print(
        "second_quotient_q0=4(p+5)(2p+3)(28p+95)/"
        "((p+4)(2p+1)(24p+85)^2)"
    )
    print("low_residual_A=10(360v^3+21v+37)")
    print("low_residual_R=-75(544v^2+1216v-307)")
    print("low_residual_S=250(2338491v-482198)")
    print("midpoint_formula=-4284000*(-1)^(m-3)*5^(m+4)/(m(m-1)(m-2)(m-3))")
    print(f"midpoint_constant_factorization={factor_integer(4284000)}")
    records = exact_small_atlas()
    for prime, s_hull, s_slopes, mid in records:
        print(
            f"exact_p={prime} S_hull={s_hull} S_slopes={s_slopes} "
            f"midpoint_valuation={mid[0]} midpoint_residue={mid[1]}"
        )
    require(targeted_midpoint_control(2969)[1:] == (
        midpoint_residue(2969), midpoint_residue(2969)
    ), "p=2969 midpoint drift")
    print(f"endpoint_constant_factorization={factor_integer(905545)}")
    print(f"penultimate_constant_factorization={factor_integer(8370560)}")
    print(f"S_constant_wall_factorization={factor_integer(482198)}")
    print(f"S_linear_wall_factorization={factor_integer(2338491)}")
    for prime in (307, 353, 601, 683, 1297, 2969):
        a_hull, r_hull, s_hull = modular_polygon(prime)
        print(
            f"modular_p={prime} A_hull={a_hull} R_hull={r_hull} "
            f"S_hull={s_hull} S_slopes={slopes(s_hull)}"
        )
    require(is_prime(2969), "2969 primality drift")
    base_left = moment_coefficients(3, 5)
    base_right = moment_coefficients(4, 5)
    delta5 = resultant(base_left, base_right)
    delta5_factors = factor_integer(delta5)
    expected_delta5 = {2: 16, 3: 7, 5: 2, 11: 2, 20747: 1, 249721: 1}
    require(delta5_factors == expected_delta5, "delta5 factorization drift")
    print(f"delta5={delta5} factorization={delta5_factors}")
    for prime, expected_gcd in (
        (11, [9, 1]),
        (20747, [15102, 1]),
        (249721, [26953, 1]),
    ):
        common = gcd_mod(base_left, base_right, prime)
        require(common == expected_gcd, "delta5 exceptional gcd drift")
        print(f"delta5_exception_p={prime} monic_gcd={common}")
    require(is_prime(20749), "20749 primality drift")
    require(is_prime(249727), "249727 primality drift")
    require(20747 + 3 == 20749 + 1 and 20747 + 5 == 20749 + 3,
            "THM3146 transport coordinates drift")
    require(249721 + 3 == 249727 - 3 and 249721 + 5 == 249727 - 1,
            "THM3159 transport coordinates drift")
    print("neighbor_transport_p20747_q20749_prime=True r=q+1 d=q+3")
    print("neighbor_transport_p249721_q249727_prime=True r=q-3 d=q-1")
    fallback_primes = [prime for prime in primes_up_to(197) if prime % 2]
    require(fallback_primes[-1] == 197, "finite fallback endpoint drift")
    require(all(prime + 3 <= 200 for prime in fallback_primes),
            "finite fallback escaped THM3124 range")
    print(
        f"finite_fallback=odd_primes_le_197 count={len(fallback_primes)} "
        f"max_r={fallback_primes[-1] + 3}"
    )
    print("p17_classification=midpoint_wall_and_exact_slope_overlap_2/17")
    print("raised_endpoint_primes=[61,2969]")
    print("all_exact_checks=True")


if __name__ == "__main__":
    main()
