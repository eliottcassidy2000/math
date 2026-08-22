#!/usr/bin/env python3
"""Exact controls for THM-3175's first-Witt/Hensel--Pluecker theorem."""

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
    if not matrix:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for column in range(len(work) - 1):
        pivot_row = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, len(work)):
            for entry in range(column + 1, len(work)):
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
    require(divisor != [0], "polynomial division by zero")
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


def moment_coefficient(n, d, index):
    return comb(n, index) * sum(
        comb(n - index, ell)
        * d ** (n - index - ell)
        * (-1) ** ell
        * factorial(2 * index + ell)
        for ell in range(n - index + 1)
    )


def fixed_endpoint(index, offset):
    return [moment_coefficient(index, offset, j) for j in range(index + 1)]


def evaluate(polynomial, value, modulus=None):
    answer = 0
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
        if modulus is not None:
            answer %= modulus
    return answer


def derivative(polynomial):
    return [index * coefficient for index, coefficient in enumerate(polynomial)][1:]


def polynomial_add(left, right):
    answer = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        answer[index] += coefficient
    for index, coefficient in enumerate(right):
        answer[index] += coefficient
    return answer


def polynomial_scale(polynomial, scalar):
    return [scalar * coefficient for coefficient in polynomial]


def polynomial_multiply(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return answer


def first_witt_wedge(left, right, lift, prime):
    modulus = prime**2
    left_value = evaluate(left, lift, modulus)
    right_value = evaluate(right, lift, modulus)
    require(left_value % prime == right_value % prime == 0, "not a residual root")
    f0 = left_value // prime % prime
    g0 = right_value // prime % prime
    f1 = evaluate(derivative(left), lift, prime)
    g1 = evaluate(derivative(right), lift, prime)
    return (f0 * g1 - g0 * f1) % prime


def scalar_window_first_jet(prime, offset, residue):
    """Evaluate the resonant pair and its v-derivative modulo p^2 in O(p)."""
    modulus = prime**2
    d = prime + offset
    first_index = prime + offset - 2
    second_index = first_index + 1
    previous = 1
    current = (d - 1 + 2 * residue) % modulus
    derivative_previous = 0
    derivative_current = 2
    d_power = d % modulus
    answer = {}
    for n in range(1, second_index):
        next_n = n + 1
        shifted = 2 * next_n * (2 * n + 1)
        previous_scale = n * next_n
        discriminant = 1 - 4 * d * residue
        following = (
            d_power * (d - next_n)
            + shifted * residue * current
            + previous_scale * discriminant * previous
        ) % modulus
        derivative_following = (
            shifted * (current + residue * derivative_current)
            + previous_scale
            * (-4 * d * previous + discriminant * derivative_previous)
        ) % modulus
        previous, current = current, following
        derivative_previous, derivative_current = (
            derivative_current,
            derivative_following,
        )
        d_power = d_power * d % modulus
        if next_n in (first_index, second_index):
            answer[next_n] = (current, derivative_current)
    return answer[first_index] + answer[second_index]


DELTAS = {
    3: 116,
    4: 586_672_128,
    5: 2_246_282_972_173_187_481_600,
    6: 44_965_855_750_876_894_470_144_533_873_830_133_760_000,
}

DELTA_FACTORS = {
    3: {2: 2, 29: 1},
    4: {2: 11, 3: 2, 7: 1, 4547: 1},
    5: {2: 16, 3: 7, 5: 2, 11: 2, 20747: 1, 249721: 1},
    6: {2: 47, 3: 15, 5: 4, 7: 1, 139: 1, 3767: 1, 12041: 1, 807241: 1},
}

# (offset, prime): (monic residual gcd, root,
#                    A/p, B/p, A', B', omega, z_A, z_B)
CASES = {
    (3, 29): ([1, 1], 28, 19, 24, 6, 1, 20, 21, 5),
    (4, 7): ([2, 1], 5, 3, 5, 3, 3, 1, 6, 3),
    (4, 4547): ([1304, 1], 3243, 2552, 3588, 4280, 1741, 3739, 3671, 3132),
    (5, 11): ([9, 1], 2, 10, 6, 4, 4, 5, 3, 4),
    (5, 20747): ([15102, 1], 5645, 6993, 17070, 8474, 18464, 7075, 13673, 8568),
    (5, 249721): ([26953, 1], 222768, 129614, 136893, 167574, 84554, 39549, 155348, 15040),
    (6, 139): ([26, 1], 113, 16, 34, 32, 28, 55, 69, 108),
    (6, 3767): ([2126, 1], 1641, 1566, 2007, 2749, 2460, 131, 2466, 1042),
    (6, 12041): ([3037, 1], 9004, 6042, 5530, 1522, 7790, 10951, 1515, 463),
    (6, 807241): ([784752, 1], 22489, 341729, 683313, 364516, 409175, 439007, 414823, 557512),
}


def frame_controls():
    prime = 101
    root = 7
    left = [-root + 2 * prime, 1 + 3 * prime]
    right = [-2 * root + 5 * prime, 2, prime]
    wedges = [
        first_witt_wedge(left, right, root + shift * prime, prime)
        for shift in (-3, -1, 0, 2, 5)
    ]
    require(len(set(wedges)) == 1, "lift independence failed")
    matrix = [
        [[1, 1], [2]],
        [[3, 0, 1], [4, -1]],
    ]
    transformed_left = polynomial_add(
        polynomial_multiply(matrix[0][0], left),
        polynomial_multiply(matrix[0][1], right),
    )
    transformed_right = polynomial_add(
        polynomial_multiply(matrix[1][0], left),
        polynomial_multiply(matrix[1][1], right),
    )
    values = [
        [evaluate(entry, root, prime) for entry in row]
        for row in matrix
    ]
    determinant = (values[0][0] * values[1][1] - values[0][1] * values[1][0]) % prime
    transformed = first_witt_wedge(
        transformed_left, transformed_right, root, prime
    )
    require(transformed == determinant * wedges[0] % prime, "frame covariance failed")
    wall = first_witt_wedge(polynomial_scale(left, prime), right, root, prime)
    require(wall == 0 and wedges[0] != 0, "determinant-wall hostile failed")
    return wedges[0], determinant, transformed, wall


def main():
    frame = frame_controls()
    records = []
    for offset in range(3, 7):
        left = fixed_endpoint(offset - 2, offset)
        right = fixed_endpoint(offset - 1, offset)
        delta = resultant(left, right)
        require(delta == DELTAS[offset], f"delta_{offset} drift")
        require(factor_integer(delta) == DELTA_FACTORS[offset], f"delta_{offset} factors")

    for (offset, prime), expected in CASES.items():
        gcd_expected, root, *jet_expected = expected
        require(is_prime(prime), "exception is not prime")
        left = fixed_endpoint(offset - 2, offset)
        right = fixed_endpoint(offset - 1, offset)
        common = gcd_mod(left, right, prime)
        require(common == gcd_expected, "residual gcd drift")
        require((-common[0]) % prime == root, "residual root drift")
        a_value, a_derivative, b_value, b_derivative = scalar_window_first_jet(
            prime, offset, root
        )
        require(a_value % prime == b_value % prime == 0, "residual lift drift")
        data = (
            a_value // prime % prime,
            b_value // prime % prime,
            a_derivative % prime,
            b_derivative % prime,
        )
        wedge = (data[0] * data[3] - data[1] * data[2]) % prime
        lift_a = -data[0] * pow(data[2], -1, prime) % prime
        lift_b = -data[1] * pow(data[3], -1, prime) % prime
        require(data + (wedge, lift_a, lift_b) == tuple(jet_expected), "first jet drift")
        require(wedge != 0 and lift_a != lift_b, "Hensel branches failed to split")
        require(
            a_derivative % prime
            == offset * evaluate(derivative(left), root, prime) % prime,
            "left Frobenius derivative scale drift",
        )
        require(
            b_derivative % prime
            == offset * evaluate(derivative(right), root, prime) % prime,
            "right Frobenius derivative scale drift",
        )
        records.append((offset, prime, root, data, wedge, lift_a, lift_b))

    print("FIRST-WITT HENSEL--PLUECKER COVARIANCE EXACT CONTROL")
    print(
        "frame_control="
        f"omega:{frame[0]},det:{frame[1]},transformed:{frame[2]},wall:{frame[3]}"
    )
    for offset in range(3, 7):
        print(
            f"delta_{offset}={DELTAS[offset]} "
            f"factorization={DELTA_FACTORS[offset]}"
        )
    for offset, prime, root, data, wedge, lift_a, lift_b in records:
        print(
            f"s={offset} p={prime} root={root} "
            f"A_over_p={data[0]} B_over_p={data[1]} "
            f"A_derivative={data[2]} B_derivative={data[3]} "
            f"omega={wedge} z_A={lift_a} z_B={lift_b}"
        )
    print(f"exception_count={len(records)} all_omega_nonzero=True")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
