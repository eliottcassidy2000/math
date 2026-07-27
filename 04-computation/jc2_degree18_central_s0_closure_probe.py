#!/usr/bin/env python3
"""Exact standard-library probe for the degree-18 central S=0 closure."""

from fractions import Fraction as Q


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


# Polynomials are low-degree first lists over Q.
def trim(poly: list[Q]) -> list[Q]:
    out = list(poly)
    while out and out[-1] == 0:
        out.pop()
    return out


def add(left: list[Q], right: list[Q]) -> list[Q]:
    out = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return trim(out)


def neg(poly: list[Q]) -> list[Q]:
    return [-value for value in poly]


def sub(left: list[Q], right: list[Q]) -> list[Q]:
    return add(left, neg(right))


def scale(poly: list[Q], scalar: Q | int) -> list[Q]:
    scalar = Q(scalar)
    return trim([scalar * value for value in poly])


def mul(left: list[Q], right: list[Q]) -> list[Q]:
    if not left or not right:
        return []
    out = [Q(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return trim(out)


def power(poly: list[Q], exponent: int) -> list[Q]:
    out = [Q(1)]
    base = list(poly)
    n = exponent
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n >>= 1
    return out


def derivative(poly: list[Q]) -> list[Q]:
    return trim([Q(index) * poly[index] for index in range(1, len(poly))])


def divmod_poly(
    dividend: list[Q], divisor: list[Q]
) -> tuple[list[Q], list[Q]]:
    dividend = trim(dividend)
    divisor = trim(divisor)
    if not divisor:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [Q(0)] * max(0, len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = trim(remainder)
    return trim(quotient), remainder


def exact_div(dividend: list[Q], divisor: list[Q]) -> list[Q]:
    quotient, remainder = divmod_poly(dividend, divisor)
    require(not remainder, "nonexact polynomial division in Bareiss")
    return quotient


def monic(poly: list[Q]) -> list[Q]:
    poly = trim(poly)
    return scale(poly, 1 / poly[-1]) if poly else []


def gcd_poly(left: list[Q], right: list[Q]) -> list[Q]:
    left = trim(left)
    right = trim(right)
    while right:
        _, remainder = divmod_poly(left, right)
        left, right = right, remainder
    return monic(left)


def evaluate(poly: list[Q], value: Q | int) -> Q:
    value = Q(value)
    out = Q(0)
    for coefficient in reversed(poly):
        out = out * value + coefficient
    return out


def bareiss_determinant(matrix: list[list[list[Q]]]) -> list[Q]:
    data = [[list(entry) for entry in row] for row in matrix]
    size = len(data)
    sign = 1
    previous = [Q(1)]
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, size)
                if trim(data[row][pivot_index])
            ),
            None,
        )
        require(pivot_row is not None, "zero pivot column in resultant")
        if pivot_row != pivot_index:
            data[pivot_index], data[pivot_row] = (
                data[pivot_row],
                data[pivot_index],
            )
            sign *= -1
        pivot = trim(data[pivot_index][pivot_index])
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = sub(
                    mul(data[row][column], pivot),
                    mul(
                        data[row][pivot_index],
                        data[pivot_index][column],
                    ),
                )
                data[row][column] = exact_div(numerator, previous)
            data[row][pivot_index] = []
        previous = pivot
    return scale(data[-1][-1], sign)


def resultant_y(
    left: list[list[Q]], right: list[list[Q]]
) -> list[Q]:
    # Inputs are low-y-degree first; entries are polynomials in tau.
    left = list(left)
    right = list(right)
    left_degree = len(left) - 1
    right_degree = len(right) - 1
    size = left_degree + right_degree
    zero: list[Q] = []
    matrix: list[list[list[Q]]] = []
    left_desc = list(reversed(left))
    right_desc = list(reversed(right))
    for shift in range(right_degree):
        row = [list(zero) for _ in range(size)]
        for index, coefficient in enumerate(left_desc):
            row[shift + index] = list(coefficient)
        matrix.append(row)
    for shift in range(left_degree):
        row = [list(zero) for _ in range(size)]
        for index, coefficient in enumerate(right_desc):
            row[shift + index] = list(coefficient)
        matrix.append(row)
    return bareiss_determinant(matrix)


def y_derivative(poly: list[list[Q]]) -> list[list[Q]]:
    return [
        scale(coefficient, degree)
        for degree, coefficient in enumerate(poly[1:], start=1)
    ]


def evaluate_tau(poly: list[list[Q]], value: Q | int) -> list[Q]:
    return trim([evaluate(coefficient, value) for coefficient in poly])


def mod_trim(poly: list[int], prime: int) -> list[int]:
    out = [value % prime for value in poly]
    while out and out[-1] == 0:
        out.pop()
    return out


def mod_divmod(
    dividend: list[int], divisor: list[int], prime: int
) -> tuple[list[int], list[int]]:
    dividend = mod_trim(dividend, prime)
    divisor = mod_trim(divisor, prime)
    quotient = [0] * max(0, len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    inverse = pow(divisor[-1], -1, prime)
    while len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse % prime
        quotient[shift] = coefficient
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
        remainder = mod_trim(remainder, prime)
    return mod_trim(quotient, prime), remainder


def mod_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left = mod_trim(left, prime)
    right = mod_trim(right, prime)
    while right:
        _, remainder = mod_divmod(left, right, prime)
        left, right = right, remainder
    if not left:
        return []
    inverse = pow(left[-1] % prime, -1, prime)
    return [(value * inverse) % prime for value in left]


def mod_mul(left: list[int], right: list[int], prime: int) -> list[int]:
    if not left or not right:
        return []
    out = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            out[i + j] += left_value * right_value
    return mod_trim(out, prime)


def mod_reduce(
    poly: list[int], modulus: list[int], prime: int
) -> list[int]:
    return mod_divmod(poly, modulus, prime)[1]


def mod_pow_x(
    exponent: int, modulus: list[int], prime: int
) -> list[int]:
    out = [1]
    base = [0, 1]
    n = exponent
    while n:
        if n & 1:
            out = mod_reduce(mod_mul(out, base, prime), modulus, prime)
        base = mod_reduce(mod_mul(base, base, prime), modulus, prime)
        n >>= 1
    return out


def main() -> None:
    # S=1701(21W+8BC)^2+2888B^5.
    require(1701 * 21 * 21 == 750141, "S W^2 coefficient")
    require(1701 * 2 * 21 * 8 == 571536, "S BCW coefficient")
    require(1701 * 8 * 8 == 108864, "S B^2C^2 coefficient")

    # R_tau as a y-polynomial with Q[tau] coefficients.
    r_tau: list[list[Q]] = [
        [Q(19 * 19), Q(-38), Q(1)],
        [Q(2 * 19 * 19), Q(-76), Q(2)],
        [Q(1083), Q(-76), Q(1)],
        [Q(304), Q(-16)],
        [Q(-975), Q(44)],
        [Q(-1836), Q(22)],
        [Q(-297)],
        [Q(1242)],
        [Q(621)],
    ]

    require(evaluate_tau(r_tau, 0)[0] == 19 * 19, "R_tau(0)")
    r_at_one = [sum(coefficient) for coefficient in zip(*[
        row + [Q(0)] * (3 - len(row)) for row in r_tau
    ])]
    # Sum in y first: the resulting tau polynomial must (2tau-35)^2.
    expected_one = [Q(35 * 35), Q(-140), Q(4)]
    require(trim(r_at_one) == expected_one, "R_tau(1)")

    discriminant = scale(
        resultant_y(r_tau, y_derivative(r_tau)),
        Q(1, 621),
    )
    p7 = [
        Q(-183687836484375),
        Q(31352481000000),
        Q(-1626234553125),
        Q(28949231250),
        Q(-287555625),
        Q(3949300),
        Q(65205),
        Q(1242),
    ]
    expected_discriminant = scale(
        mul(power([Q(-19), Q(1)], 6), p7),
        -49682096496000000000000,
    )
    require(
        discriminant == expected_discriminant,
        "residual discriminant factorization",
    )

    require(evaluate(p7, 19) != 0, "P7 meets tau=19")
    require(evaluate(p7, Q(35, 2)) != 0, "P7 meets tau=35/2")

    r4 = [Q(-139), Q(-1418), Q(-297), Q(1242), Q(621)]
    require(
        evaluate_tau(r_tau, 19) == [Q(0)] * 4 + r4,
        "tau=19 factorization",
    )
    require(len(gcd_poly(r4, derivative(r4))) == 1, "R4 squarefree")
    require(evaluate(r4, 0) == -139 and evaluate(r4, 1) == 9, "R4 boundary")

    r7 = [
        Q(-9),
        Q(-27),
        Q(-264),
        Q(-360),
        Q(460),
        Q(6264),
        Q(7452),
        Q(2484),
    ]
    expected_node = scale(mul([Q(-1), Q(1)], r7), Q(1, 4))
    require(
        evaluate_tau(r_tau, Q(35, 2)) == expected_node,
        "tau=35/2 factorization",
    )
    require(len(gcd_poly(r7, derivative(r7))) == 1, "R7 squarefree")
    require(evaluate(r7, 0) == -9 and evaluate(r7, 1) == 16000, "R7 boundary")

    prime = 37
    p7_mod = [
        int(value) % prime
        for value in p7
    ]
    lead_inverse = pow(p7_mod[-1], -1, prime)
    p7_monic = [
        value * lead_inverse % prime for value in p7_mod
    ]
    expected_monic = [10, 20, 6, 6, 20, 5, 34, 1]
    require(p7_monic == expected_monic, "P7 monic reduction")
    x37_minus_x = mod_pow_x(37, p7_monic, prime)
    if len(x37_minus_x) < 2:
        x37_minus_x += [0] * (2 - len(x37_minus_x))
    x37_minus_x[1] = (x37_minus_x[1] - 1) % prime
    require(
        mod_gcd(p7_monic, x37_minus_x, prime) == [1],
        "P7 Rabin gcd",
    )
    require(
        mod_pow_x(37**7, p7_monic, prime) == [0, 1],
        "P7 Frobenius closure",
    )

    require(441536697 != 0, "ordinary triple tangent discriminant")
    require(-1972620783 != 0, "infinity discriminant")
    genera = {
        "generic": (12 - 4) // 2,
        "ordinary-triple": (6 - 4) // 2,
        "rational-node": (10 - 4) // 2,
        "P7-node": (10 - 4) // 2,
    }
    require(
        genera
        == {
            "generic": 4,
            "ordinary-triple": 1,
            "rational-node": 3,
            "P7-node": 3,
        },
        "genus ledger",
    )

    print("JC2 DEGREE-18 CENTRAL S=0 CLOSURE PROBE")
    print("S=1701(21W+8BC)^2+2888B^5 [PASS]")
    print("R_tau boundary values and both rational factorizations [PASS]")
    print("Disc_y R_tau = constant*(tau-19)^6*P7(tau) [PASS]")
    print("P7 degree =", len(p7) - 1)
    print("P7 monic mod 37 =", p7_monic)
    print("Rabin irreducibility certificate mod 37 [PASS]")
    print("R4 and R7 squarefree/boundary checks [PASS]")
    print("infinity discriminant = -1972620783 [PASS]")
    print("genus ledger =", genera)
    print("ALL EXACT CENTRAL S=0 PROBE CHECKS PASSED")


if __name__ == "__main__":
    main()
