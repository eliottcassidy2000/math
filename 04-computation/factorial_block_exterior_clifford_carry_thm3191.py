#!/usr/bin/env python3
"""Exact controls for THM-3191's exterior Clifford/carry law."""


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


def matrix_power(matrix, exponent, modulus):
    answer = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    while exponent:
        if exponent & 1:
            answer = matrix_product(answer, matrix, modulus)
        matrix = matrix_product(matrix, matrix, modulus)
        exponent //= 2
    return answer


def matrix_vector(matrix, vector, modulus):
    return [
        sum(matrix[i][j] * vector[j] for j in range(3)) % modulus
        for i in range(3)
    ]


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


def transfer_interval(start, stop, d, v, modulus=None):
    answer = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for n in range(start, stop):
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


def determinant_three(matrix):
    return (
        matrix[0][0]
        * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1]
        * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2]
        * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def valuation(number, prime):
    require(number != 0, "valuation of zero")
    answer = 0
    while number % prime == 0:
        number //= prime
        answer += 1
    return answer


def minimum_entry_valuation(matrix, prime):
    entries = [entry for row in matrix for entry in row if entry]
    require(entries, "nonzero valuation matrix")
    return min(valuation(entry, prime) for entry in entries)


def factorial_valuation(number, prime):
    answer = 0
    while number:
        number //= prime
        answer += number
    return answer


def unit_part(number, prime):
    while number % prime == 0:
        number //= prime
    return number % prime


def reduce_divided(matrix, divisor, prime):
    require(all(entry % divisor == 0 for row in matrix for entry in row),
            "exact layer divisibility")
    return [[entry // divisor % prime for entry in row] for row in matrix]


def one_block_operator(prime, residue, v):
    discriminant = (1 - 4 * residue * v) % prime
    require(discriminant, "unit discriminant block")
    character = pow(discriminant, (prime - 1) // 2, prime)
    alpha = -character % prime
    exterior_tail = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    pre_reset = transfer_interval(0, prime - 1, residue, v, prime)
    exterior_tail = [[entry % prime for entry in row]
                     for row in exterior_square(pre_reset)]
    beta = pre_reset[0][2]
    gamma = pre_reset[1][2]
    expected_tail = [
        [1, alpha * (beta + gamma) % prime, -alpha * beta % prime],
        [0, alpha, 0],
        [0, -alpha % prime, alpha],
    ]
    require(exterior_tail == expected_tail, "THM-3188 tail layer")

    reset_compound = exterior_square(transfer(prime - 1, residue, v))
    reset_layer = reduce_divided(reset_compound, prime, prime)
    operator = matrix_product(reset_layer, exterior_tail, prime)
    core = [
        [0, 2 * residue + 1, -1],
        [0, -(1 + 2 * v), 2 * v],
        [0, -2 * (residue + v + 1), 2 * v + 1],
    ]
    expected = [
        [alpha * residue * entry % prime for entry in row] for row in core
    ]
    require(operator == expected, "extension-free one-block operator")
    return operator, discriminant


def check_one_block_clifford_law():
    parameter_checks = 0
    chart_checks = 0
    for prime in [3, 5, 7, 11, 13]:
        for residue in range(1, prime):
            for v in range(prime):
                discriminant = (1 - 4 * residue * v) % prime
                if not discriminant:
                    continue
                operator, _ = one_block_operator(prime, residue, v)
                cube = matrix_power(operator, 3, prime)
                expected_cube = [
                    [residue**2 * discriminant * entry % prime for entry in row]
                    for row in operator
                ]
                require(cube == expected_cube, "cubic Clifford law")
                require([row[0] for row in operator] == [0, 0, 0],
                        "right missing-wedge kernel")
                require(
                    [
                        (operator[0][column] - operator[1][column]
                         + operator[2][column]) % prime
                        for column in range(3)
                    ] == [0, 0, 0],
                    "left conormal kernel",
                )
                require(any(entry % prime for row in exterior_square(operator)
                            for entry in row), "one-block rank two")
                for coordinate_02 in range(prime):
                    for coordinate_12 in range(prime):
                        if coordinate_02 == coordinate_12 == 0:
                            continue
                        image = matrix_vector(
                            operator, [0, coordinate_02, coordinate_12], prime
                        )
                        require(image[1] or image[2],
                                "complementary projective charts")
                        chart_checks += 1
                parameter_checks += 1
    return parameter_checks, chart_checks


def check_carried_block_layers():
    checks = 0
    for prime in [3, 5, 7]:
        for residue in range(1, prime):
            for v in range(prime):
                discriminant = (1 - 4 * residue * v) % prime
                if not discriminant:
                    continue
                operator, _ = one_block_operator(prime, residue, v)
                for block in range(1, prime + 3):
                    start = (block - 1) * prime
                    stop = block * prime
                    exact_block = transfer_interval(start, stop, residue, v)
                    exact_compound = exterior_square(exact_block)
                    thickness = 1 + valuation(block, prime) if block % prime == 0 else 1
                    divisor = prime**thickness
                    layer = reduce_divided(exact_compound, divisor, prime)
                    expected = [
                        [unit_part(block, prime) * entry % prime for entry in row]
                        for row in operator
                    ]
                    require(layer == expected, "carried block unit layer")
                    checks += 1
    return checks


def check_global_block_products():
    checks = 0
    for prime in [3, 5, 7]:
        for residue in range(1, prime):
            for v in range(prime):
                discriminant = (1 - 4 * residue * v) % prime
                if not discriminant:
                    continue
                operator, _ = one_block_operator(prime, residue, v)
                exact_product = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                unit_factorial = 1
                for blocks in range(1, 8):
                    exact_block = transfer_interval(
                        (blocks - 1) * prime, blocks * prime, residue, v
                    )
                    exact_product = matrix_product(exact_block, exact_product)
                    unit_factorial = (
                        unit_factorial * unit_part(blocks, prime)
                    ) % prime
                    total_thickness = factorial_valuation(blocks * prime, prime)
                    require(
                        total_thickness
                        == blocks + factorial_valuation(blocks, prime),
                        "global carry thickness",
                    )
                    expected_state = [
                        [0, 0, pow(residue, blocks, prime)],
                        [0, 0, pow(residue, blocks, prime)],
                        [0, 0, pow(residue, blocks, prime)],
                    ]
                    require(
                        [[entry % prime for entry in row]
                         for row in exact_product] == expected_state,
                        "global rank-one state reduction",
                    )

                    exact_compound = exterior_square(exact_product)
                    layer = reduce_divided(
                        exact_compound, prime**total_thickness, prime
                    )
                    expected_layer = [
                        [unit_factorial * entry % prime for entry in row]
                        for row in matrix_power(operator, blocks, prime)
                    ]
                    require(layer == expected_layer,
                            "global Clifford/carry layer")
                    require(any(entry for row in exterior_square(layer)
                                for entry in row), "global layer rank two")

                    require(minimum_entry_valuation(exact_product, prime) == 0,
                            "global state first divisor")
                    require(
                        minimum_entry_valuation(exact_compound, prime)
                        == total_thickness,
                        "global state second divisor",
                    )
                    require(
                        valuation(determinant_three(exact_product), prime)
                        == 2 * total_thickness,
                        "global state determinant",
                    )
                    require(
                        minimum_entry_valuation(
                            exterior_square(exact_compound), prime
                        ) == 2 * total_thickness,
                        "global exterior second divisor",
                    )
                    require(
                        valuation(determinant_three(exact_compound), prime)
                        == 4 * total_thickness,
                        "global exterior determinant",
                    )

                    require(
                        all(exact_compound[row][0] % prime**(2 * total_thickness)
                            == 0 for row in range(3)),
                        "global squared missing-column divisibility",
                    )
                    return_column = [
                        exact_compound[row][0]
                        // prime**(2 * total_thickness) % prime
                        for row in range(3)
                    ]
                    expected_return = [
                        unit_factorial**2
                        * pow(-discriminant, blocks, prime) % prime,
                        0,
                        0,
                    ]
                    require(return_column == expected_return,
                            "global squared-scale return")
                    checks += 1
    return checks


def main():
    parameter_checks, chart_checks = check_one_block_clifford_law()
    carried_checks = check_carried_block_layers()
    global_checks = check_global_block_products()
    print("FACTORIAL-BLOCK EXTERIOR CLIFFORD CARRY EXACT CONTROL")
    print(f"one_block_parameter_checks={parameter_checks}")
    print(f"complementary_chart_vector_checks={chart_checks}")
    print(f"carried_local_block_checks={carried_checks}")
    print(f"global_fixed_parameter_product_checks={global_checks}")
    print("one_block_cubic=D^3=s^2*Delta*D")
    print("global_p_primary_state_smith=1,p^H,p^H")
    print("global_p_primary_exterior_smith=p^H,p^H,p^(2H)")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
