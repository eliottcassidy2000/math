#!/usr/bin/env python3
"""Exact controls for THM-3176's offset-six Euclidean--Newton flag.

The companion uses integer and modular arithmetic only.
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
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [number for number in range(2, limit + 1) if sieve[number]]


def bareiss_determinant(matrix):
    work = [row[:] for row in matrix]
    size = len(work)
    if not size:
        return 1
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]), None
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
                require(numerator % previous == 0, "Bareiss division failed")
                work[row][entry] = numerator // previous
            work[row][column] = 0
        previous = pivot
    return sign * work[-1][-1]


def resultant(left, right):
    left_degree = len(left) - 1
    right_degree = len(right) - 1
    left_descending = list(reversed(left))
    right_descending = list(reversed(right))
    matrix = []
    for shift in range(right_degree):
        matrix.append(
            [0] * shift
            + left_descending
            + [0] * (right_degree - 1 - shift)
        )
    for shift in range(left_degree):
        matrix.append(
            [0] * shift
            + right_descending
            + [0] * (left_degree - 1 - shift)
        )
    return bareiss_determinant(matrix)


def trim_mod(polynomial, prime):
    answer = [coefficient % prime for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def remainder_mod(dividend, divisor, prime):
    dividend = trim_mod(dividend, prime)
    divisor = trim_mod(divisor, prime)
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
    inverse = pow(left[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in left]


def polynomial_add(left, right):
    answer = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] += coefficient
    for index, coefficient in enumerate(right):
        answer[index] += coefficient
    return answer


def polynomial_scale(polynomial, scalar, shift=0):
    return [0] * shift + [scalar * coefficient for coefficient in polynomial]


def polynomial_multiply(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            answer[left_index + right_index] += left_coefficient * right_coefficient
    return answer


def moment_coefficient(n, d, index):
    return comb(n, index) * sum(
        comb(n - index, ell)
        * d ** (n - index - ell)
        * (-1) ** ell
        * factorial(2 * index + ell)
        for ell in range(n - index + 1)
    )


def moment_coefficients_exact(n, d):
    return [moment_coefficient(n, d, index) for index in range(n + 1)]


def moment_coefficients_mod(n_final, d, modulus):
    """Build A_n by the exact three-term recurrence, coefficientwise."""
    if n_final == 0:
        return [1]
    previous = [1]
    current = [(d - 1) % modulus, 2 % modulus]
    d_power = d % modulus
    for n in range(1, n_final):
        next_n = n + 1
        answer = [0] * (n + 2)
        answer[0] = d_power * (d - next_n) % modulus
        shifted_current = 2 * next_n * (2 * n + 1)
        previous_scale = n * next_n
        for index, coefficient in enumerate(current):
            answer[index + 1] = (
                answer[index + 1] + shifted_current * coefficient
            ) % modulus
        for index, coefficient in enumerate(previous):
            answer[index] = (
                answer[index] + previous_scale * coefficient
            ) % modulus
            answer[index + 1] = (
                answer[index + 1]
                - 4 * d * previous_scale * coefficient
            ) % modulus
        previous, current = current, answer
        d_power = d_power * d % modulus
    return current


def first_parameters(prime):
    return (
        2 * prime + 7,
        (2 * prime + 10) * (2 * prime + 9),
        2 * (prime + 5) * (prime + 6),
    )


def first_remainder(a, b, prime, modulus=None):
    denominator, leading, constant = first_parameters(prime)
    answer = []
    for index in range(prime + 6):
        value = denominator * b[index]
        if index:
            value -= denominator * leading * a[index - 1]
        if index < len(a):
            value += constant * a[index]
        if modulus is not None:
            value %= modulus
        answer.append(value)
    require(answer[-1] == answer[-2] == 0, "first top cancellation failed")
    return answer[:-2]


def second_parameters(prime):
    h = 24 * prime + 109
    denominator = (prime + 5) * (2 * prime + 3) * h**2
    linear = -2 * (2 * prime + 5) * (2 * prime + 7) * (2 * prime + 3) * h
    constant = 4 * (prime + 6) * (2 * prime + 5) * (28 * prime + 123)
    return denominator, linear, constant


def second_remainder(a, first, prime, modulus=None):
    denominator, linear, constant = second_parameters(prime)
    answer = []
    for index in range(prime + 5):
        value = denominator * a[index]
        if index:
            value -= linear * first[index - 1]
        if index < len(first):
            value -= constant * first[index]
        if modulus is not None:
            value %= modulus
        answer.append(value)
    require(answer[-1] == answer[-2] == 0, "second top cancellation failed")
    return answer[:-2]


def third_polynomials(prime):
    j_polynomial = (
        256 * prime**4
        - 27648 * prime**3
        - 365600 * prime**2
        - 1528800 * prime
        - 2096649
    )
    k_polynomial = (
        5120 * prime**5
        - 810240 * prime**4
        - 14447872 * prime**3
        - 92004672 * prime**2
        - 256323456 * prime
        - 265142241
    )
    return j_polynomial, k_polynomial


def third_parameters(prime):
    j_polynomial, k_polynomial = third_polynomials(prime)
    denominator = (2 * prime - 1) * (2 * prime + 7) * j_polynomial**2
    linear = (
        -2
        * (2 * prime + 1)
        * (2 * prime + 3)
        * (24 * prime + 109)
        * (2 * prime - 1)
        * j_polynomial
    )
    constant = 4 * (prime + 6) * (2 * prime + 1) * k_polynomial
    return denominator, linear, constant


def third_remainder(first, second, prime, modulus=None):
    denominator, linear, constant = third_parameters(prime)
    answer = []
    for index in range(prime + 4):
        value = denominator * first[index]
        if index:
            value -= linear * second[index - 1]
        if index < len(second):
            value -= constant * second[index]
        if modulus is not None:
            value %= modulus
        answer.append(value)
    require(answer[-1] == answer[-2] == 0, "third top cancellation failed")
    return answer[:-2]


def all_rows_exact(prime):
    d = prime + 6
    a = moment_coefficients_exact(prime + 4, d)
    b = moment_coefficients_exact(prime + 5, d)
    first = first_remainder(a, b, prime)
    second = second_remainder(a, first, prime)
    third = third_remainder(first, second, prime)
    return a, b, first, second, third


def all_rows_mod(prime, cap=5):
    modulus = prime**cap
    d = prime + 6
    a = moment_coefficients_mod(prime + 4, d, modulus)
    b = moment_coefficients_mod(prime + 5, d, modulus)
    first = first_remainder(a, b, prime, modulus)
    second = second_remainder(a, first, prime, modulus)
    third = third_remainder(first, second, prime, modulus)
    return a, b, first, second, third


def valuation(number, prime, cap=5):
    number %= prime**cap
    if number == 0:
        return cap
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def lower_hull(polynomial, prime, cap=5):
    hull = []
    for point in enumerate(valuation(value, prime, cap) for value in polynomial):
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


def endpoint_polynomials(prime):
    endpoint = (
        327680 * prime**8
        + 91422720 * prime**7
        - 1525587968 * prime**6
        - 58319874048 * prime**5
        - 605420542720 * prime**4
        - 3106619397120 * prime**3
        - 8698849881696 * prime**2
        - 12762150724080 * prime
        - 7697077854849
    )
    penultimate = (
        327680 * prime**9
        + 87818240 * prime**8
        - 2581766144 * prime**7
        - 38490294272 * prime**6
        + 127672001792 * prime**5
        + 4526613136640 * prime**4
        + 30758409758112 * prime**3
        + 98688343652784 * prime**2
        + 157343795925183 * prime
        + 100505677023531
    )
    return endpoint, penultimate


def expected_third_high(prime):
    endpoint_polynomial, penultimate_polynomial = endpoint_polynomials(prime)
    common = (
        prime
        * (prime - 3)
        * (prime - 2)
        * (prime - 1)
        * (prime + 1)
        * (prime + 2)
        * (prime + 3)
        * (prime + 4)
        * (prime + 5)
        * (2 * prime - 7)
        * (2 * prime + 3)
        * (2 * prime + 7)
        * (24 * prime + 109) ** 2
        * factorial(2 * prime - 8)
    )
    endpoint = -32 * common * (2 * prime - 5) * endpoint_polynomial
    penultimate = 16 * common * penultimate_polynomial
    return endpoint, penultimate


def evaluate_polynomial(polynomial, value, modulus=None):
    answer = 0
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
        if modulus is not None:
            answer %= modulus
    return answer


def derivative(polynomial):
    return [index * coefficient for index, coefficient in enumerate(polynomial)][1:]


def first_witt_wedge(left, right, residue, lift, prime):
    modulus = prime * prime
    left_value = evaluate_polynomial(left, lift, modulus)
    right_value = evaluate_polynomial(right, lift, modulus)
    require(left_value % prime == right_value % prime == 0, "non-root lift")
    left_quotient = left_value // prime % prime
    right_quotient = right_value // prime % prime
    left_derivative = evaluate_polynomial(derivative(left), residue, prime)
    right_derivative = evaluate_polynomial(derivative(right), residue, prime)
    return (
        left_quotient * right_derivative
        - right_quotient * left_derivative
    ) % prime


def scalar_window_first_jet(prime, offset, residue):
    """Return the resonant A,B first jet modulo p^2 in O(p) time."""
    modulus = prime * prime
    d = prime + offset
    first_index = prime + offset - 2
    second_index = first_index + 1
    previous = 1
    current = (d - 1 + 2 * residue) % modulus
    derivative_previous = 0
    derivative_current = 2
    d_power = d % modulus
    a = b = derivative_a = derivative_b = None
    for n in range(1, second_index):
        next_n = n + 1
        shifted = 2 * next_n * (2 * n + 1)
        previous_scale = n * next_n
        affine = (1 - 4 * d * residue) % modulus
        following = (
            d_power * (d - next_n)
            + shifted * residue * current
            + previous_scale * affine * previous
        ) % modulus
        derivative_following = (
            shifted * (current + residue * derivative_current)
            + previous_scale
            * (-4 * d * previous + affine * derivative_previous)
        ) % modulus
        previous, current = current, following
        derivative_previous, derivative_current = (
            derivative_current,
            derivative_following,
        )
        d_power = d_power * d % modulus
        if next_n == first_index:
            a, derivative_a = current, derivative_current
        if next_n == second_index:
            b, derivative_b = current, derivative_current
    return a, b, derivative_a, derivative_b


def scalar_moment_and_derivative(prime, residue):
    """Offset-six compatibility wrapper."""
    return scalar_window_first_jet(prime, 6, residue)


def selected_third_coefficients(prime, indices, cap=4):
    """Compute selected T coefficients without constructing a degree-p row.

    The requested indices must be less than p.  Each moment coefficient is
    evaluated from the exact coefficient sum in O(p) modular operations.  A
    shared factorial/inverse bank and memoization make the midpoint wall at
    p=232961 practical while remaining independent of the polynomial-row
    recurrence used for the finite controls.
    """
    require(all(0 <= index < prime for index in indices), "selected index range")
    modulus = prime**cap
    d = prime + 6
    factorials = [1] * (2 * prime + 21)
    for number in range(1, len(factorials)):
        factorials[number] = factorials[number - 1] * number % modulus
    inverses = [0] * prime
    inverses[1] = 1
    for number in range(2, prime):
        inverses[number] = (
            -(modulus // number) * inverses[modulus % number]
        ) % modulus
    inverse_d = pow(d, -1, modulus)
    cache = {}

    def moment(n, index):
        if index < 0 or index > n:
            return 0
        key = (n, index)
        if key in cache:
            return cache[key]
        remaining = n - index
        outer_numerator = 1
        for number in range(n - index + 1, n + 1):
            outer_numerator = outer_numerator * number % modulus
        outer = (
            outer_numerator * pow(factorials[index], -1, modulus) % modulus
        )
        inner_binomial = 1
        d_power = pow(d, remaining, modulus)
        falling_factorial = factorials[2 * index]
        total = 0
        sign = 1
        for ell in range(remaining + 1):
            total = (
                total + sign * inner_binomial * d_power * falling_factorial
            ) % modulus
            if ell < remaining:
                inner_binomial = (
                    inner_binomial
                    * (remaining - ell)
                    * inverses[ell + 1]
                ) % modulus
                d_power = d_power * inverse_d % modulus
                falling_factorial = (
                    falling_factorial * (2 * index + ell + 1) % modulus
                )
                sign = -sign
        cache[key] = outer * total % modulus
        return cache[key]

    def a(index):
        return moment(prime + 4, index)

    def b(index):
        return moment(prime + 5, index)

    denominator_1, leading_1, constant_1 = first_parameters(prime)

    def first(index):
        return (
            denominator_1 * b(index)
            - denominator_1 * leading_1 * a(index - 1)
            + constant_1 * a(index)
        ) % modulus

    denominator_2, linear_2, constant_2 = second_parameters(prime)

    def second(index):
        return (
            denominator_2 * a(index)
            - linear_2 * first(index - 1)
            - constant_2 * first(index)
        ) % modulus

    denominator_3, linear_3, constant_3 = third_parameters(prime)

    def third(index):
        return (
            denominator_3 * first(index)
            - linear_3 * second(index - 1)
            - constant_3 * second(index)
        ) % modulus

    return {index: third(index) for index in indices}


def main():
    # A deterministic DVR control for lift-independence, polynomial-frame
    # covariance, and the sharp determinant-wall boundary.
    control_prime = 101
    control_root = 7
    control_left = [
        -control_root + control_prime * 2,
        1 + control_prime * 3,
    ]
    control_right = [
        -2 * control_root + control_prime * 5,
        2,
        control_prime,
    ]
    control_wedges = [
        first_witt_wedge(
            control_left,
            control_right,
            control_root,
            control_root + control_prime * shift,
            control_prime,
        )
        for shift in [-3, -1, 0, 2, 5]
    ]
    require(len(set(control_wedges)) == 1, "first-Witt lift dependence")
    control_wedge = control_wedges[0]
    matrix = [
        [[1, 1], [2]],
        [[3, 0, 1], [4, -1]],
    ]
    transformed_left = polynomial_add(
        polynomial_multiply(matrix[0][0], control_left),
        polynomial_multiply(matrix[0][1], control_right),
    )
    transformed_right = polynomial_add(
        polynomial_multiply(matrix[1][0], control_left),
        polynomial_multiply(matrix[1][1], control_right),
    )
    matrix_values = [
        [evaluate_polynomial(entry, control_root, control_prime) for entry in row]
        for row in matrix
    ]
    matrix_determinant = (
        matrix_values[0][0] * matrix_values[1][1]
        - matrix_values[0][1] * matrix_values[1][0]
    ) % control_prime
    transformed_wedge = first_witt_wedge(
        transformed_left,
        transformed_right,
        control_root,
        control_root,
        control_prime,
    )
    require(
        transformed_wedge == matrix_determinant * control_wedge % control_prime,
        "first-Witt frame covariance",
    )
    determinant_wall_wedge = first_witt_wedge(
        polynomial_scale(control_left, control_prime),
        control_right,
        control_root,
        control_root,
        control_prime,
    )
    require(determinant_wall_wedge == 0, "determinant wall did not kill wedge")

    residual_a = [744, 384, 864, -2880, 40320]
    residual_b = [4056, 2160, 1440, 57600, -604800, 3628800]
    residual_a = polynomial_scale(residual_a, 6)
    residual_b = polynomial_scale(residual_b, 6)

    denominator_1, leading_1, constant_1 = first_parameters(0)
    residual_r = polynomial_add(
        polynomial_add(
            polynomial_scale(residual_b, denominator_1),
            polynomial_scale(residual_a, -denominator_1 * leading_1, 1),
        ),
        polynomial_scale(residual_a, constant_1),
    )
    while residual_r[-1] == 0:
        residual_r.pop()

    denominator_2, linear_2, constant_2 = second_parameters(0)
    residual_s = polynomial_add(
        polynomial_add(
            polynomial_scale(residual_a, denominator_2),
            polynomial_scale(residual_r, -linear_2, 1),
        ),
        polynomial_scale(residual_r, -constant_2),
    )
    while residual_s[-1] == 0:
        residual_s.pop()

    denominator_3, linear_3, constant_3 = third_parameters(0)
    residual_t = polynomial_add(
        polynomial_add(
            polynomial_scale(residual_r, denominator_3),
            polynomial_scale(residual_s, -linear_3, 1),
        ),
        polynomial_scale(residual_s, -constant_3),
    )
    while residual_t[-1] == 0:
        residual_t.pop()

    expected_a = polynomial_scale([31, 16, 36, -120, 1680], 144)
    expected_r = polynomial_scale([3043, -17940, -7500, -13080], 144)
    expected_s = polynomial_scale([-375143, 3212382, -2795532], 15120)
    expected_t = polynomial_scale(
        [-51108355843, 392547973190, -174357330840], 970059888
    )
    require(residual_a == expected_a, "A residual mismatch")
    require(residual_r == expected_r, "R residual mismatch")
    require(residual_s == expected_s, "S residual mismatch")
    require(residual_t == expected_t, "T residual mismatch")

    for prime in [5, 7, 11, 13, 17]:
        *_, third = all_rows_exact(prime)
        expected_endpoint, expected_penultimate = expected_third_high(prime)
        require(third[-1] == expected_endpoint, "endpoint formula mismatch")
        require(third[-2] == expected_penultimate, "penultimate formula mismatch")

    midpoint_constant = 474570855293572800
    midpoint_factors = 2**6 * 3**7 * 5**2 * 7**2 * 109**2 * 232961
    require(midpoint_constant == midpoint_factors, "midpoint factorization mismatch")

    polygon_rows = {}
    for prime in [109, 197, 211]:
        a, _, first, second, third = all_rows_mod(prime)
        m = (prime - 1) // 2
        midpoint = m + 3
        denominator = 1
        for shift in range(5):
            denominator = denominator * (m - shift) % prime
        if prime != 109:
            midpoint_expected = (
                -midpoint_constant
                * pow(-1, m, prime)
                * pow(6, m + 5, prime)
                * pow(denominator, -1, prime)
            ) % prime
            require(
                third[midpoint] // prime % prime == midpoint_expected,
                "midpoint formula mismatch",
            )
        polygon_rows[prime] = tuple(
            lower_hull(row, prime) for row in [a, first, second, third]
        )

    expected_197 = (
        [(0, 0), (4, 0), (201, 2)],
        [(0, 0), (3, 0), (200, 2)],
        [(0, 0), (2, 0), (199, 2)],
        [(0, 0), (2, 0), (101, 1), (198, 2)],
    )
    require(polygon_rows[197] == expected_197, "p=197 polygon mismatch")
    selected_197 = selected_third_coefficients(197, [100, 101], cap=5)
    full_197_third = all_rows_mod(197)[-1]
    require(
        selected_197 == {100: full_197_third[100], 101: full_197_third[101]},
        "selected coefficient engine disagrees with full recurrence",
    )
    require(
        polygon_rows[109][-1]
        == [(0, 2), (2, 2), (57, 3), (110, 4)],
        "p=109 global-shift polygon mismatch",
    )

    frame_wall = 232961
    frame_midpoint = (frame_wall + 5) // 2
    frame_selected = selected_third_coefficients(
        frame_wall, [frame_midpoint - 1, frame_midpoint], cap=4
    )
    require(
        valuation(frame_selected[frame_midpoint - 1], frame_wall, 4) == 1,
        "frame-wall pre-midpoint valuation",
    )
    require(
        valuation(frame_selected[frame_midpoint], frame_wall, 4) == 2,
        "frame-wall midpoint valuation",
    )

    endpoint_wall = 1067703961
    require(is_prime(endpoint_wall), "endpoint wall is not prime")
    endpoint_value, penultimate_value = endpoint_polynomials(endpoint_wall)
    require(endpoint_value % endpoint_wall**2 != 0, "endpoint wall not simple")
    require(endpoint_value % endpoint_wall == 0, "endpoint wall missing")
    require(penultimate_value % endpoint_wall == 567766679, "penultimate wall")
    require(
        endpoint_value // endpoint_wall % endpoint_wall == 114714544,
        "endpoint wall quotient",
    )

    delta_six = resultant(
        [744, 384, 864, -2880, 40320],
        [4056, 2160, 1440, 57600, -604800, 3628800],
    )
    expected_delta = 44965855750876894470144533873830133760000
    require(delta_six == expected_delta, "offset-six resultant mismatch")
    require(
        delta_six == 2**47 * 3**15 * 5**4 * 7 * 139 * 3767 * 12041 * 807241,
        "offset-six factorization mismatch",
    )
    exceptional_gcds = {
        prime: gcd_mod(
            [744, 384, 864, -2880, 40320],
            [4056, 2160, 1440, 57600, -604800, 3628800],
            prime,
        )
        for prime in [139, 3767, 12041, 807241]
    }
    require(exceptional_gcds[139] == [26, 1], "p=139 gcd")
    require(exceptional_gcds[3767] == [2126, 1], "p=3767 gcd")
    require(exceptional_gcds[12041] == [3037, 1], "p=12041 gcd")
    require(exceptional_gcds[807241] == [784752, 1], "p=807241 gcd")

    require(is_prime(3769) and is_prime(12043), "redundant p+2 transport primality")
    require(is_prime(51108355847), "redundant constant-wall transport primality")

    hard_prime = 807241
    hard_root = 22489
    a_value, b_value, a_derivative, b_derivative = scalar_moment_and_derivative(
        hard_prime, hard_root
    )
    require(a_value % hard_prime == b_value % hard_prime == 0, "hard root")
    a_quotient = a_value // hard_prime % hard_prime
    b_quotient = b_value // hard_prime % hard_prime
    a_derivative %= hard_prime
    b_derivative %= hard_prime
    require(
        (a_quotient, b_quotient, a_derivative, b_derivative)
        == (341729, 683313, 364516, 409175),
        "hard first-jet data mismatch",
    )
    hard_lift_a = -a_quotient * pow(a_derivative, -1, hard_prime) % hard_prime
    hard_lift_b = -b_quotient * pow(b_derivative, -1, hard_prime) % hard_prime
    hensel_wedge = (
        a_quotient * b_derivative - b_quotient * a_derivative
    ) % hard_prime
    require((hard_lift_a, hard_lift_b, hensel_wedge) == (414823, 557512, 439007),
            "hard Hensel split mismatch")

    residual_a_base = [744, 384, 864, -2880, 40320]
    residual_b_base = [4056, 2160, 1440, 57600, -604800, 3628800]
    residual_a_derivative = derivative(residual_a_base)
    residual_b_derivative = derivative(residual_b_base)
    require(evaluate_polynomial(residual_a_base, hard_root, hard_prime) == 0,
            "hard residual A root")
    require(evaluate_polynomial(residual_b_base, hard_root, hard_prime) == 0,
            "hard residual B root")
    require(
        evaluate_polynomial(residual_a_derivative, hard_root, hard_prime)
        == 329833,
        "hard residual A derivative",
    )
    require(
        evaluate_polynomial(residual_b_derivative, hard_root, hard_prime)
        == 202736,
        "hard residual B derivative",
    )
    require(a_derivative == 6 * 329833 % hard_prime, "A derivative scale")
    require(b_derivative == 6 * 202736 % hard_prime, "B derivative scale")

    offset_six_hensel_controls = {
        (139, 113): (16, 34, 32, 28, 55, 69, 108),
        (3767, 1641): (1566, 2007, 2749, 2460, 131, 2466, 1042),
        (12041, 9004): (6042, 5530, 1522, 7790, 10951, 1515, 463),
        (807241, 22489): (
            341729,
            683313,
            364516,
            409175,
            439007,
            414823,
            557512,
        ),
    }
    for (prime, root), expected in offset_six_hensel_controls.items():
        value_a, value_b, derivative_a, derivative_b = scalar_window_first_jet(
            prime, 6, root
        )
        require(value_a % prime == value_b % prime == 0, "offset-six root")
        quotient_a = value_a // prime % prime
        quotient_b = value_b // prime % prime
        derivative_a %= prime
        derivative_b %= prime
        wedge = (
            quotient_a * derivative_b - quotient_b * derivative_a
        ) % prime
        lift_a = -quotient_a * pow(derivative_a, -1, prime) % prime
        lift_b = -quotient_b * pow(derivative_b, -1, prime) % prime
        require(
            (
                quotient_a,
                quotient_b,
                derivative_a,
                derivative_b,
                wedge,
                lift_a,
                lift_b,
            )
            == expected,
            "offset-six Hensel control mismatch",
        )

    prior_hensel_controls = {
        (4, 4547, 3243): (2552, 3588, 4280, 1741, 3739, 3671, 3132),
        (5, 20747, 5645): (6993, 17070, 8474, 18464, 7075, 13673, 8568),
        (5, 249721, 222768): (
            129614,
            136893,
            167574,
            84554,
            39549,
            155348,
            15040,
        ),
    }
    for (offset, prime, root), expected in prior_hensel_controls.items():
        value_a, value_b, derivative_a, derivative_b = scalar_window_first_jet(
            prime, offset, root
        )
        require(value_a % prime == value_b % prime == 0, "prior residual root")
        data = (
            value_a // prime % prime,
            value_b // prime % prime,
            derivative_a % prime,
            derivative_b % prime,
        )
        wedge = (data[0] * data[3] - data[1] * data[2]) % prime
        lift_a = -data[0] * pow(data[2], -1, prime) % prime
        lift_b = -data[1] * pow(data[3], -1, prime) % prime
        require(
            data + (wedge, lift_a, lift_b) == expected,
            "prior Hensel control mismatch",
        )

    bounded_primes = [prime for prime in primes_up_to(193) if prime % 2]

    print("OFFSET-SIX EUCLIDEAN--NEWTON / HENSEL-WEDGE EXACT CONTROL")
    print(
        "first_witt_control = "
        f"omega:{control_wedge},det:{matrix_determinant},"
        f"transformed:{transformed_wedge},wall:{determinant_wall_wedge}"
    )
    print("low_degrees A,R,S,T = 4,3,2,2")
    print("third_J_constant = -2096649 = -3^2*232961")
    print("third_K_constant = -265142241 = -3^4*7*13^2*2767")
    print("T_low_constants = 51108355843,392547973190,174357330840")
    print("T_linear_walls = 22447,1748777 (horizontal face unchanged)")
    print("frame_midpoint_wall = 232961; direct normalized hull checked")
    print("T_constant_wall = 51108355843; positive slopes unchanged")
    print(f"midpoint_constant = {midpoint_constant}")
    print("midpoint_factors = 2^6*3^7*5^2*7^2*109^2*232961")
    print(f"generic_hulls_p197 = {polygon_rows[197]}")
    print(f"generic_hulls_p211 = {polygon_rows[211]}")
    print(f"global_shift_hulls_p109 = {polygon_rows[109]}")
    print(
        "frame_wall_selected = "
        f"T(k-1)/p:{frame_selected[frame_midpoint - 1] // frame_wall % frame_wall},"
        "T(k)/p^2:"
        f"{frame_selected[frame_midpoint] // frame_wall**2 % frame_wall}"
    )
    print(
        "frame_wall_hull = "
        "[(0,0),(1,0),((p+3)/2,1),(p+1,2)]"
    )
    print("endpoint_wall_prime = 1067703961")
    print("endpoint_wall_quotient = 114714544")
    print("endpoint_wall_penultimate = 567766679")
    print("endpoint_wall_hull = [(0,0),(2,0),((p+5)/2,1),(p,2),(p+1,3)]")
    print(f"delta6 = {delta_six}")
    print("delta6_factors = 2^47*3^15*5^4*7*139*3767*12041*807241")
    print(f"exceptional_gcds = {exceptional_gcds}")
    print("redundant_transports = 3767->3769,12041->12043")
    print("constant_wall = initial low-coordinate rise only; positive slopes unchanged")
    print(f"hard_prime_root = ({hard_prime},{hard_root})")
    print(
        "hard_first_jet = "
        f"A/p:{a_quotient},B/p:{b_quotient},A':{a_derivative},B':{b_derivative}"
    )
    print(f"hard_lift_digits = A:{hard_lift_a},B:{hard_lift_b}")
    print(f"hard_hensel_wedge = {hensel_wedge}")
    print(f"offset_six_hensel_controls = {offset_six_hensel_controls}")
    print(f"prior_hensel_controls = {prior_hensel_controls}")
    print(f"bounded_odd_primes_through_193 = {len(bounded_primes)}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
