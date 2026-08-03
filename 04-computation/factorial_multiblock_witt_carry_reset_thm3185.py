#!/usr/bin/env python3
"""Exact controls for THM-3185's multiblock Witt-carry hierarchy."""

from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(row):
    answer = list(row)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def direct_state(n, d):
    moment = [0] * (n + 1)
    weighted = [0] * (n + 1)
    for j in range(n + 1):
        for ell in range(n - j + 1):
            multinomial = factorial(n) // (
                factorial(j) * factorial(ell) * factorial(n - j - ell)
            )
            scalar = multinomial * d ** (n - j - ell) * (-1) ** ell
            moment[j] += scalar * factorial(2 * j + ell)
            weighted[j] += scalar * factorial(2 * j + ell + 1)
    return moment, weighted, d**n


def reduce_row(row, prime):
    return trim([value % prime for value in row])


def valuation(number, prime):
    require(number != 0, "valuation of zero")
    answer = 0
    while number % prime == 0:
        number //= prime
        answer += 1
    return answer


def matrix_product(left, right, modulus):
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(3)) % modulus
            for j in range(3)
        ]
        for i in range(3)
    ]


def transfer_matrix(n, d, v):
    return [
        [-(n + 1), 2 * v * (n + 1), d],
        [
            -(n + 1) * (2 * n + 3 + 2 * d),
            (n + 1) * (1 + 2 * v * (2 * n + 3)),
            d * (2 * n + 3),
        ],
        [0, 0, d],
    ]


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


def minimum_entry_valuation(matrix, prime):
    entries = [entry for row in matrix for entry in row if entry]
    require(entries, "nonzero determinantal-divisor matrix")
    return min(valuation(entry, prime) for entry in entries)


def transfer_product(length, d, v, prime):
    answer = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for n in range(length):
        answer = matrix_product(transfer_matrix(n, d, v), answer, prime)
    return answer


def evaluate_row(row, v, prime):
    return sum(coefficient * pow(v, degree, prime)
               for degree, coefficient in enumerate(row)) % prime


def check_block_propagators():
    checks = 0
    for prime in [5, 7, 11]:
        for residue in range(1, min(prime, 5)):
            for v in [1, 2]:
                block = transfer_product(prime, residue, v, prime)
                expected_block = [[0, 0, residue % prime]] * 3
                require(block == expected_block, "rank-one full-block propagator")
                checks += 1
                for d in [residue, 2 * prime + residue]:
                    for blocks in [1, 2, 3]:
                        for a in range(min(prime, 5)):
                            product = transfer_product(blocks * prime + a,
                                                       d, v, prime)
                            small = direct_state(a, residue)
                            state = [
                                evaluate_row(small[0], v, prime),
                                evaluate_row(small[1], v, prime),
                                small[2] % prime,
                            ]
                            scale = pow(residue, blocks, prime)
                            expected = [
                                [0, 0, scale * coordinate % prime]
                                for coordinate in state
                            ]
                            require(product == expected,
                                    "iterated rank-one block propagator")
                            checks += 1
    return checks


def check_multiblock_descent():
    checks = 0
    for prime in [5, 7, 11]:
        for residue in range(1, min(prime, 5)):
            for d in [residue, 2 * prime + residue, 4 * prime + residue]:
                for blocks in range(4):
                    for a in range(min(prime, 5)):
                        n = blocks * prime + a
                        large = direct_state(n, d)
                        small = direct_state(a, residue)
                        scale = pow(residue, blocks, prime)
                        for large_row, small_row in zip(large[:2], small[:2]):
                            require(
                                reduce_row(large_row, prime)
                                == reduce_row(
                                    [scale * value for value in small_row], prime
                                ),
                                "multiblock state descent",
                            )
                        require(
                            large[2] % prime == scale * small[2] % prime,
                            "multiblock power descent",
                        )
                        checks += 1
    return checks


def check_resonant_pairs():
    checks = 0
    for prime, residue in [(7, 2), (11, 3), (13, 5), (17, 6)]:
        require(prime > 2 * (residue - 1), "degree-preservation fixture")
        for blocks in [1, 2, 4]:
            d = blocks * prime + residue
            scale = pow(residue, blocks, prime)
            for n, small_n in [(d - 2, residue - 2), (d - 1, residue - 1)]:
                large = direct_state(n, d)[0]
                small = direct_state(small_n, residue)[0]
                require(
                    reduce_row(large, prime)
                    == reduce_row([scale * value for value in small], prime),
                    "resonant pair descent",
                )
                require(len(reduce_row(small, prime)) == small_n + 1,
                        "residual degree preservation")
                checks += 1
    return checks


def check_unit_resets():
    fixtures = [
        (5, 1, 7, 1),
        (5, 5, 7, 1),
        (7, 2, 9, 2),
        (7, 7, 9, 2),
        (11, 11, 14, 2),
        (11, 121, 14, 2),
    ]
    for prime, ell, d, v in fixtures:
        n_block = ell * prime
        discriminant = 1 - 4 * d * v
        require(d % prime and discriminant % prime, "unit reset fixture")
        h = valuation(n_block, prime)

        weighted = transfer_matrix(n_block - 1, d, v)
        compound = exterior_square(weighted)
        expected_compound = [
            [-n_block**2 * discriminant, 2 * n_block * d**2,
             -n_block * d],
            [0, -n_block * d, 2 * n_block * v * d],
            [0, -n_block * d * (2 * n_block + 1 + 2 * d),
             n_block * d * (1 + 2 * v * (2 * n_block + 1))],
        ]
        require(compound == expected_compound,
                "exact exterior-square reconstruction")
        require(
            all(entry == 0 or valuation(entry, prime) >= h
                for row in compound for entry in row),
            "all weighted 2x2 minors contain reset thickness",
        )

        x_minor = compound[0][2]
        x_determinant = -d * n_block**2 * discriminant
        require(valuation(x_minor, prime) == h, "X second divisor")
        require(valuation(x_determinant, prime) == 2 * h,
                "X determinant divisor")

        scalar_unit_minor = n_block - d
        scalar_determinant = -(n_block - 1) * n_block * discriminant * d
        require(valuation(scalar_unit_minor, prime) == 0,
                "scalar unit minor")
        require(valuation(scalar_determinant, prime) == h,
                "scalar determinant divisor")
        require(valuation((n_block - 1) * discriminant, prime) == 0,
                "input gauge index")
        require(valuation(n_block * discriminant, prime) == h,
                "output gauge index")

        require(all(entry % n_block == 0 for row in compound for entry in row),
                "exact N divisibility of exterior layer")
        first_layer = [
            [entry // n_block for entry in row]
            for row in compound
        ]
        expected_first_layer = [
            [0, 2 * d**2, -d],
            [0, -d, 2 * v * d],
            [0, -d * (1 + 2 * d), d * (1 + 2 * v)],
        ]
        require(
            [[entry % prime for entry in row] for row in first_layer]
            == [[entry % prime for entry in row]
                for row in expected_first_layer],
            "derived normalized exterior first layer",
        )
        minors = []
        for row_a, row_b in [(0, 1), (0, 2), (1, 2)]:
            minors.append(
                (first_layer[row_a][1] * first_layer[row_b][2]
                 - first_layer[row_a][2] * first_layer[row_b][1]) % prime
            )
        expected = (-d**2 * discriminant) % prime
        require(minors == [expected] * 3, "exterior first-layer minors")
        require(expected != 0, "normalized exterior rank two")
        require([row[0] % prime for row in first_layer] == [0, 0, 0],
                "normalized exterior right kernel")
        require(
            [
                (first_layer[0][column] - first_layer[1][column]
                 + first_layer[2][column]) % prime
                for column in range(3)
            ] == [0, 0, 0],
            "normalized exterior left kernel",
        )
        require(
            all(compound[row][0] % n_block**2 == 0 for row in range(3)),
            "exact N-squared missing-column divisibility",
        )
        second_layer_column = [
            compound[row][0] // n_block**2 % prime for row in range(3)
        ]
        require(second_layer_column == [-discriminant % prime, 0, 0],
                "derived exterior second-layer return")
    return len(fixtures)


def check_discriminant_walls():
    # Exact orders t=1,2 for p=7,d=9.
    fixtures = [(7, 1, 9, 22, 1), (7, 7, 9, 211, 2)]
    for prime, ell, d, v, t in fixtures:
        n_block = ell * prime
        discriminant = 1 - 4 * d * v
        h = valuation(n_block, prime)
        require(valuation(discriminant, prime) == t,
                "discriminant-wall order")

        weighted = transfer_matrix(n_block - 1, d, v)
        weighted_minors = exterior_square(weighted)
        weighted_determinant = determinant_three(weighted)
        require(minimum_entry_valuation(weighted, prime) == 0,
                "X wall first divisor")
        require(minimum_entry_valuation(weighted_minors, prime) == h,
                "X wall second determinantal divisor")
        require(valuation(weighted_determinant, prime) == 2 * h + t,
                "X wall determinant")

        scalar = [
            [2 * n_block * (2 * n_block - 1) * v,
             (n_block - 1) * n_block * discriminant,
             d - n_block],
            [1, 0, 0],
            [0, 0, d],
        ]
        scalar_minors = exterior_square(scalar)
        scalar_determinant = determinant_three(scalar)
        require(minimum_entry_valuation(scalar, prime) == 0,
                "scalar wall first divisor")
        require(minimum_entry_valuation(scalar_minors, prime) == 0,
                "scalar wall second determinantal divisor")
        require(valuation(scalar_determinant, prime) == h + t,
                "scalar wall determinant")

        exterior_minors = exterior_square(weighted_minors)
        exterior_determinant = determinant_three(weighted_minors)
        require(minimum_entry_valuation(weighted_minors, prime) == h,
                "exterior wall first divisor")
        require(minimum_entry_valuation(exterior_minors, prime) == 2 * h + t,
                "exterior wall second determinantal divisor")
        require(valuation(exterior_determinant, prime) == 4 * h + 2 * t,
                "exterior wall determinant")
        require(valuation((n_block - 1) * discriminant, prime) == t,
                "wall input gauge")
        require(valuation(n_block * discriminant, prime) == h + t,
                "wall output gauge")
    return len(fixtures)


def factorial_valuation(n, prime):
    answer = 0
    while n:
        n //= prime
        answer += n
    return answer


def check_global_ledgers():
    checks = 0
    for prime in [5, 7, 11]:
        for blocks in [1, 2, prime - 1, prime, prime + 2]:
            endpoint = blocks * prime
            factorial_length = factorial_valuation(endpoint, prime)
            require(
                factorial_length == blocks + factorial_valuation(blocks, prime),
                "Legendre block identity",
            )
            reset_sum = sum(
                1 + valuation(ell, prime) if ell % prime == 0 else 1
                for ell in range(1, blocks + 1)
            )
            require(reset_sum == factorial_length, "reset thickness sum")
            weighted_length = 2 * factorial_length
            scalar_length = (
                factorial_valuation(endpoint - 1, prime) + factorial_length
            )
            gauge_endpoint = valuation(endpoint, prime)
            require(weighted_length == scalar_length + gauge_endpoint,
                    "global gauge telescope")
            checks += 1
    return checks


def check_primitive_colour_cycles():
    checks = 0
    for prime, primitive in [(5, 2), (7, 3), (11, 2)]:
        colours = {pow(primitive, block, prime) for block in range(prime - 1)}
        require(colours == set(range(1, prime)), "primitive colour cycle")
        checks += 1
    return checks


def main():
    descent_checks = check_multiblock_descent()
    propagator_checks = check_block_propagators()
    pair_checks = check_resonant_pairs()
    reset_checks = check_unit_resets()
    wall_checks = check_discriminant_walls()
    ledger_checks = check_global_ledgers()
    colour_checks = check_primitive_colour_cycles()
    print("FACTORIAL MULTIBLOCK WITT-CARRY RESET EXACT CONTROL")
    print(f"multiblock_state_descent_checks={descent_checks}")
    print(f"rank_one_block_propagator_checks={propagator_checks}")
    print(f"resonant_pair_descent_checks={pair_checks}")
    print(f"carry_thickness_reset_checks={reset_checks}")
    print(f"discriminant_wall_checks={wall_checks}")
    print(f"global_carry_ledger_checks={ledger_checks}")
    print(f"primitive_colour_cycle_checks={colour_checks}")
    print("X_smith=1,p^h,p^h")
    print("scalar_smith=1,1,p^h")
    print("exterior_smith=p^h,p^h,p^(2h)")
    print("wall_X_smith=1,p^h,p^(h+t)")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
