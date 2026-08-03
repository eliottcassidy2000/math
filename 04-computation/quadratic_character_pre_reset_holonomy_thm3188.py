#!/usr/bin/env python3
"""Exact controls for THM-3188's quadratic-character pre-reset holonomy."""

from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matrix_product(left, right, modulus=None):
    answer = [
        [sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]
    if modulus is not None:
        answer = [[entry % modulus for entry in row] for row in answer]
    return answer


def transfer(n, d, v):
    return [
        [-(n + 1), 2 * v * (n + 1), d],
        [
            -(n + 1) * (2 * n + 3 + 2 * d),
            (n + 1) * (1 + 2 * v * (2 * n + 3)),
            d * (2 * n + 3),
        ],
        [0, 0, d],
    ]


def transfer_product(length, d, v, modulus=None):
    answer = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for n in range(length):
        answer = matrix_product(transfer(n, d, v), answer, modulus)
    return answer


def exterior_square(matrix):
    pairs = [(0, 1), (0, 2), (1, 2)]
    return [
        [
            (matrix[row_a][column_a] * matrix[row_b][column_b]
             - matrix[row_a][column_b] * matrix[row_b][column_a])
            for column_a, column_b in pairs
        ]
        for row_a, row_b in pairs
    ]


def reduce_matrix(matrix, prime):
    return [[entry % prime for entry in row] for row in matrix]


def odd_double_factorial(number):
    require(number >= 1 and number % 2 == 1, "positive odd double factorial")
    answer = 1
    for value in range(1, number + 1, 2):
        answer *= value
    return answer


def closed_continuant_coefficients(length):
    coefficients = {}
    for dimers in range(length // 2 + 1):
        numerator = odd_double_factorial(2 * length - 2 * dimers + 1)
        denominator = odd_double_factorial(2 * dimers + 1)
        require(numerator % denominator == 0, "double-factorial quotient")
        coefficient = (
            comb(length - dimers, dimers) * numerator // denominator
        )
        coefficients[(length - 2 * dimers, dimers)] = coefficient
    return coefficients


def recursive_continuant_coefficients(length):
    # K_m for diagonal entries (3x,5x,...,(2m+1)x), dimer weight Delta.
    previous_previous = {(0, 0): 1}
    if length == 0:
        return previous_previous
    previous = {(1, 0): 3}
    if length == 1:
        return previous
    for index in range(2, length + 1):
        answer = {}
        diagonal = 2 * index + 1
        for (x_degree, delta_degree), coefficient in previous.items():
            key = (x_degree + 1, delta_degree)
            answer[key] = answer.get(key, 0) + diagonal * coefficient
        for (x_degree, delta_degree), coefficient in previous_previous.items():
            key = (x_degree, delta_degree + 1)
            answer[key] = answer.get(key, 0) + coefficient
        previous_previous, previous = previous, answer
    return previous


def continuant(start, end, x, discriminant, modulus):
    # Determinant/matching recurrence with diagonal (2j+1)x and dimer Delta.
    if start > end:
        return 1
    previous_previous = 1
    previous = (2 * start + 1) * x % modulus
    if start == end:
        return previous
    for index in range(start + 1, end + 1):
        answer = (
            (2 * index + 1) * x * previous
            + discriminant * previous_previous
        ) % modulus
        previous_previous, previous = previous, answer
    return previous


def check_closed_continuant_formula():
    for length in range(13):
        require(
            recursive_continuant_coefficients(length)
            == closed_continuant_coefficients(length),
            "closed continuant coefficient formula",
        )
    return 13


def check_holonomy_and_walls():
    unit_checks = 0
    wall_checks = 0
    for prime in [3, 5, 7, 11, 13, 17, 19]:
        for residue in range(1, prime):
            for v in range(prime):
                discriminant = (1 - 4 * residue * v) % prime
                pre_reset = transfer_product(prime - 1, residue, v, prime)
                top = [row[:2] for row in pre_reset[:2]]
                if discriminant == 0:
                    require(top == [[0, 0], [0, 0]],
                            "discriminant-wall homogeneous collapse")
                    require(exterior_square(pre_reset) == [[0, 0, 0]] * 3,
                            "discriminant-wall exterior collapse")
                    wall_checks += 1
                    continue

                character = pow(discriminant, (prime - 1) // 2, prime)
                alpha = (-character) % prime
                require(
                    top == [[alpha, 0], [(-alpha) % prime, alpha]],
                    "quadratic-character homogeneous shear",
                )
                require(pre_reset[2] == [0, 0, 1], "unit bottom holonomy")
                beta = pre_reset[0][2]
                gamma = pre_reset[1][2]
                expected = [
                    [alpha, 0, beta],
                    [(-alpha) % prime, alpha, gamma],
                    [0, 0, 1],
                ]
                require(pre_reset == expected, "full pre-reset extension form")

                compound = reduce_matrix(exterior_square(pre_reset), prime)
                expected_compound = [
                    [1, alpha * (beta + gamma) % prime,
                     -alpha * beta % prime],
                    [0, alpha, 0],
                    [0, -alpha % prime, alpha],
                ]
                require(compound == expected_compound,
                        "exterior pre-reset holonomy")
                require([row[0] for row in compound] == [1, 0, 0],
                        "determinant wedge fixed")

                if prime >= 5:
                    x = 2 * v % prime
                    require(
                        continuant(1, prime - 3, x, discriminant, prime)
                        == pow(discriminant, (prime - 3) // 2, prime),
                        "surviving endpoint continuant",
                    )
                    require(
                        continuant(1, prime - 2, x, discriminant, prime) == 0,
                        "full odd continuant vanishing",
                    )
                    require(
                        continuant(2, prime - 3, x, discriminant, prime) == 0,
                        "reflected shifted continuant vanishing",
                    )
                unit_checks += 1
    return unit_checks, wall_checks


def check_full_block_exterior_flags():
    checks = 0
    for prime in [3, 5, 7, 11, 13]:
        for residue in range(1, prime):
            for v in range(prime):
                discriminant = 1 - 4 * residue * v
                if discriminant % prime == 0:
                    continue
                pre_reset = transfer_product(prime - 1, residue, v)
                reset = transfer(prime - 1, residue, v)
                full_block = matrix_product(reset, pre_reset)
                reset_compound = exterior_square(reset)
                pre_compound = exterior_square(pre_reset)
                block_compound = exterior_square(full_block)
                require(
                    block_compound == matrix_product(
                        reset_compound, pre_compound
                    ),
                    "compound multiplicativity",
                )
                require(
                    all(entry % prime == 0
                        for row in block_compound for entry in row),
                    "full-block exterior first divisibility",
                )
                first_layer = reduce_matrix(
                    [[entry // prime for entry in row]
                     for row in block_compound],
                    prime,
                )
                require([row[0] for row in first_layer] == [0, 0, 0],
                        "full-block right missing wedge")
                require(
                    [
                        (first_layer[0][column] - first_layer[1][column]
                         + first_layer[2][column]) % prime
                        for column in range(3)
                    ] == [0, 0, 0],
                    "full-block left conormal flag",
                )
                first_layer_minors = reduce_matrix(
                    exterior_square(first_layer), prime
                )
                require(any(entry for row in first_layer_minors for entry in row),
                        "full-block exterior rank two")
                require(
                    all(block_compound[row][0] % prime**2 == 0
                        for row in range(3)),
                    "full-block squared missing-column divisibility",
                )
                return_column = [
                    block_compound[row][0] // prime**2 % prime
                    for row in range(3)
                ]
                require(
                    return_column == [-discriminant % prime, 0, 0],
                    "full-block squared-scale return",
                )
                checks += 1
    return checks


def main():
    formula_checks = check_closed_continuant_formula()
    holonomy_checks, wall_checks = check_holonomy_and_walls()
    flag_checks = check_full_block_exterior_flags()
    print("QUADRATIC-CHARACTER PRE-RESET HOLONOMY EXACT CONTROL")
    print(f"closed_continuant_formula_checks={formula_checks}")
    print(f"unit_holonomy_parameter_checks={holonomy_checks}")
    print(f"discriminant_wall_parameter_checks={wall_checks}")
    print(f"full_block_exterior_flag_checks={flag_checks}")
    print("homogeneous_holonomy=-chi(Delta)*[[1,0],[-1,1]]")
    print("determinant_wedge_fixed=1")
    print("wall_homogeneous_rank=0")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
