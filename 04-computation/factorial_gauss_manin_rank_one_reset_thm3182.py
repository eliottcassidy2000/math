#!/usr/bin/env python3
"""Exact controls for THM-3182's factorial Gauss--Manin reset."""

from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(row):
    answer = list(row)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def add(left, right):
    answer = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return trim(answer)


def scale(row, scalar):
    return trim([scalar * value for value in row])


def shift_v(row):
    return [0] + row


def direct_state(n, d):
    """Coefficient rows in v for M_n, X_n, and the scalar d^n."""
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


def transfer(n, d, state):
    moment, weighted, power = state
    next_moment = add(
        add(scale(moment, -(n + 1)),
            scale(shift_v(weighted), 2 * (n + 1))),
        [d * power],
    )
    next_weighted = add(
        add(
            scale(moment, -(n + 1) * (2 * n + 3 + 2 * d)),
            add(
                scale(weighted, n + 1),
                scale(shift_v(weighted), 2 * (n + 1) * (2 * n + 3)),
            ),
        ),
        [d * (2 * n + 3) * power],
    )
    return next_moment, next_weighted, d * power


def scalar_recurrence(n, d, previous, current):
    return add(
        add(
            [d**n * (d - n - 1)],
            scale(shift_v(current), 2 * (n + 1) * (2 * n + 1)),
        ),
        add(
            scale(previous, n * (n + 1)),
            scale(shift_v(previous), -4 * d * n * (n + 1)),
        ),
    )


def determinant_formula(n, d, v):
    matrix = [
        [-(n + 1), 2 * v * (n + 1), d],
        [-(n + 1) * (2 * n + 3 + 2 * d),
         (n + 1) * (1 + 2 * v * (2 * n + 3)), d * (2 * n + 3)],
        [0, 0, d],
    ]
    determinant = (
        matrix[0][0] * matrix[1][1] * matrix[2][2]
        - matrix[0][1] * matrix[1][0] * matrix[2][2]
    )
    expected = -d * (n + 1) ** 2 * (1 - 4 * d * v)
    require(determinant == expected, "transfer determinant mismatch")
    return determinant


def reduce_row(row, prime):
    return trim([value % prime for value in row])


def valuation(number, prime):
    require(number != 0, "valuation of zero")
    answer = 0
    while number % prime == 0:
        number //= prime
        answer += 1
    return answer


def compound_matrix(matrix):
    pairs = [(0, 1), (0, 2), (1, 2)]
    answer = []
    for row_a, row_b in pairs:
        row = []
        for column_a, column_b in pairs:
            row.append(
                matrix[row_a][column_a] * matrix[row_b][column_b]
                - matrix[row_a][column_b] * matrix[row_b][column_a]
            )
        answer.append(row)
    return answer


def main():
    transfer_checks = 0
    scalar_checks = 0
    gauge_checks = 0
    determinant_checks = 0
    for d in [3, 8, 17]:
        for n in range(8):
            state = direct_state(n, d)
            require(transfer(n, d, state) == direct_state(n + 1, d),
                    "transfer/direct mismatch")
            transfer_checks += 1
            for v in [-3, 0, 2, 5]:
                determinant_formula(n, d, v)
                determinant_checks += 1
            if n >= 1:
                previous = direct_state(n - 1, d)[0]
                current = state[0]
                require(
                    scalar_recurrence(n, d, previous, current)
                    == direct_state(n + 1, d)[0],
                    "eliminated scalar recurrence mismatch",
                )
                scalar_checks += 1
                gauge_left = scale(shift_v(state[1]), 2)
                gauge_right = add(
                    add(
                        current,
                        scale(shift_v(current), 2 * (2 * n + 1)),
                    ),
                    add(
                        add(scale(previous, n),
                            scale(shift_v(previous), -4 * d * n)),
                        [-state[2]],
                    ),
                )
                require(gauge_left == gauge_right,
                        "weighted/scalar gauge identity")
                gauge_checks += 1

    descent_checks = 0
    reset_checks = 0
    for prime in [7, 11, 13, 17]:
        for offset in range(2, min(prime, 8)):
            d = prime + offset
            reset = [
                [-prime, 2 * prime, d],
                [-prime * (2 * prime + 1 + 2 * d),
                 prime * (1 + 2 * (2 * prime + 1)), d * (2 * prime + 1)],
                [0, 0, d],
            ]
            require(
                [[entry % prime for entry in row] for row in reset]
                == [[0, 0, offset], [0, 0, offset], [0, 0, offset]],
                "rank-one reset mismatch",
            )
            reset_checks += 1
            for a in range(7):
                large = direct_state(prime + a, d)
                small = direct_state(a, offset)
                for large_entry, small_entry in zip(large[:2], small[:2]):
                    require(
                        reduce_row(large_entry, prime)
                        == reduce_row(scale(small_entry, offset), prime),
                        "state Frobenius descent mismatch",
                    )
                require(
                    large[2] % prime == offset * small[2] % prime,
                    "power Frobenius descent mismatch",
                )
                descent_checks += 1

    smith_checks = 0
    compound_checks = 0
    scalar_smith_checks = 0
    for prime, offset, v in [(7, 2, 2), (11, 3, 2), (13, 5, 4), (17, 6, 3)]:
        d = prime + offset
        discriminant = 1 - 4 * d * v
        require(d % prime != 0 and discriminant % prime != 0,
                "nonunit Smith fixture")
        determinant = -d * prime**2 * discriminant
        minor = -prime * d
        require(valuation(d, prime) == 0, "unit-entry divisor")
        require(valuation(minor, prime) == 1, "second divisor")
        require(valuation(determinant, prime) == 2, "determinant divisor")
        smith_checks += 1

        reset_matrix = [
            [-prime, 2 * prime * v, d],
            [-prime * (2 * prime + 1 + 2 * d),
             prime * (1 + 2 * v * (2 * prime + 1)), d * (2 * prime + 1)],
            [0, 0, d],
        ]
        compound = compound_matrix(reset_matrix)
        expected_compound = [
            [prime**2 * (4 * prime * v + 4 * offset * v - 1),
             2 * prime * d**2, -prime * d],
            [0, -prime * d, 2 * prime * v * d],
            [0, -prime * d * (4 * prime + 2 * offset + 1),
             prime * d * (4 * prime * v + 2 * v + 1)],
        ]
        require(compound == expected_compound, "exterior-square matrix")
        first_layer = [
            [(entry // prime) % prime for entry in row]
            for row in compound
        ]
        expected_first = [
            [0, 2 * offset**2 % prime, -offset % prime],
            [0, -offset % prime, 2 * offset * v % prime],
            [0, -offset * (2 * offset + 1) % prime,
             offset * (2 * v + 1) % prime],
        ]
        require(first_layer == expected_first, "exterior first layer")
        left_kernel = [1, -1, 1]
        require(
            all(
                sum(
                    left_kernel[row] * first_layer[row][column]
                    for row in range(3)
                ) % prime == 0
                for column in range(3)
            ),
            "exterior left kernel",
        )
        require(all(first_layer[row][0] == 0 for row in range(3)),
                "exterior right kernel")
        nonzero_minors = []
        for row_a, row_b in [(0, 1), (0, 2), (1, 2)]:
            value = (
                first_layer[row_a][1] * first_layer[row_b][2]
                - first_layer[row_a][2] * first_layer[row_b][1]
            ) % prime
            nonzero_minors.append(value)
        expected_minor = (-offset**2 * (1 - 4 * offset * v)) % prime
        require(nonzero_minors == [expected_minor] * 3,
                "exterior nonzero minors")
        require(
            (compound[0][0] // prime**2) % prime
            == (-(1 - 4 * offset * v)) % prime,
            "exterior second-layer return",
        )
        compound_checks += 1

        # The scalar companion lattice Y=(M_n,M_(n-1),D_n) has reset
        # Smith type (1,1,p).  The gauge to X has determinant
        # n*Delta/(2v), unit at n=p-1 and p-divisible at n=p.
        scalar_b = (prime - 1) * prime * discriminant
        scalar_determinant = -scalar_b * d
        scalar_unit_minor = -(d - prime)
        require(valuation(d, prime) == 0, "scalar first divisor")
        require(valuation(scalar_unit_minor, prime) == 0,
                "scalar unit 2x2 minor")
        require(valuation(scalar_determinant, prime) == 1,
                "scalar companion determinant")
        require(valuation((prime - 1) * discriminant, prime) == 0,
                "input gauge determinant")
        require(valuation(prime * discriminant, prime) == 1,
                "output gauge determinant")
        scalar_smith_checks += 1

    # Simple discriminant wall: p=7,d=9,v=1 has Delta=-35.
    wall_prime, wall_d, wall_v = 7, 9, 1
    wall_discriminant = 1 - 4 * wall_d * wall_v
    wall_determinant = -wall_d * wall_prime**2 * wall_discriminant
    require(valuation(-wall_prime * wall_d, wall_prime) == 1,
            "wall second divisor")
    require(valuation(wall_determinant, wall_prime) == 3,
            "simple discriminant-wall thickening")

    # If p|s then d is nonunit and the reset reduction is zero, not rank one.
    hostile_prime, hostile_offset = 7, 7
    hostile_d = hostile_prime + hostile_offset
    require(hostile_d % hostile_prime == 0, "p|s hostile setup")
    require(
        all(
            entry % hostile_prime == 0
            for row in [
                [-hostile_prime, 2 * hostile_prime, hostile_d],
                [-hostile_prime * (2 * hostile_prime + 1 + 2 * hostile_d),
                 hostile_prime * (1 + 2 * (2 * hostile_prime + 1)),
                 hostile_d * (2 * hostile_prime + 1)],
                [0, 0, hostile_d],
            ]
            for entry in row
        ),
        "p|s reset is not zero",
    )

    print("FACTORIAL GAUSS--MANIN RANK-ONE RESET EXACT CONTROL")
    print(f"transfer_direct_checks={transfer_checks}")
    print(f"scalar_recurrence_checks={scalar_checks}")
    print(f"weighted_scalar_gauge_checks={gauge_checks}")
    print(f"determinant_checks={determinant_checks}")
    print(f"rank_one_reset_checks={reset_checks}")
    print(f"full_state_descent_checks={descent_checks}")
    print(f"smith_type_1_p_p_checks={smith_checks}")
    print(f"exterior_square_p_p_p2_checks={compound_checks}")
    print(f"scalar_companion_smith_1_1_p_checks={scalar_smith_checks}")
    print("simple_discriminant_wall_smith=1,p,p^2")
    print("p_divides_s_hostile=zero_reset")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
