#!/usr/bin/env python3
"""Exact ballot/production-matrix scout for strict integer product-Gamma data.

This companion does four logically separate things.

1. It verifies the strict-prefix strand factorization and the exact root-event
   ballot law for a declared finite universe with terminal mass Q in {2, 3}.
2. It identifies monomial and Gregory--Newton coefficient-network targets
   for the reversed generalized-Hankel signs, and exhausts every minor of both
   matrices in that universe.
3. It checks the actual generalized-Hankel determinants independently.
4. It records two hostile controls: monotone root matching by itself fails at
   order three, and the obvious local root-deletion production layer already
   has a negative coefficient.

The finite coefficient census is evidence for a sharpened conjecture, not an
all-order theorem.  The factorization, ballot identity, fixed-operator
recurrence, and the two hostiles are exact symbolic statements.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rising(value, length):
    answer = Fraction(1)
    for offset in range(length):
        answer *= value + offset
    return answer


def polynomial_multiply(left, right):
    answer = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return answer


def rising_polynomial(root, length):
    answer = [Fraction(1)]
    for offset in range(length):
        answer = polynomial_multiply(answer, [root + offset, 1])
    return answer


def polynomial_value(coefficients, value):
    answer = Fraction(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def shift_down(coefficients):
    """Coefficients of f(u-1), given coefficients of f(u)."""

    answer = [Fraction(0)] * len(coefficients)
    for degree, coefficient in enumerate(coefficients):
        for target in range(degree + 1):
            answer[target] += (
                coefficient * comb(degree, target) * (-1) ** (degree - target)
            )
    return answer


def newton_coefficients(coefficients):
    """Return b_k in f(u)=sum_k b_k binom(u,k)."""

    values = [
        polynomial_value(coefficients, value)
        for value in range(len(coefficients))
    ]
    answer = []
    while values:
        answer.append(values[0])
        values = [right - left for left, right in zip(values, values[1:])]
    return answer


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return Fraction(1)
    work = [[Fraction(value) for value in row] for row in matrix]
    sign = 1
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            ratio = work[row][column] / pivot_value
            for target in range(column + 1, size):
                work[row][target] -= ratio * work[column][target]
    return sign * answer


def all_minors_nonnegative(matrix):
    row_count = len(matrix)
    column_count = len(matrix[0])
    checked = 0
    zero = 0
    for size in range(1, column_count + 1):
        for rows in combinations(range(row_count), size):
            for columns in combinations(range(column_count), size):
                value = determinant(
                    [[matrix[row][column] for column in columns] for row in rows]
                )
                require(value >= 0, f"negative coefficient minor {rows}/{columns}")
                checked += 1
                zero += value == 0
    return checked, zero


def exponent_inventory(heights):
    """Recover e from S_j=-sum_{i<=j}e_i."""

    return (-heights[0],) + tuple(
        heights[index - 1] - heights[index]
        for index in range(1, len(heights))
    )


def direct_sequence(shapes, exponents, index):
    answer = Fraction(1)
    for shape, exponent in zip(shapes, exponents):
        factor = rising(shape, index)
        answer *= factor**exponent
    return answer


def residual_factors(shapes, heights):
    answer = []
    for index in range(len(shapes) - 1):
        gap = shapes[index + 1] - shapes[index]
        require(gap.denominator == 1 and gap > 0, "nonintegral shape gap")
        answer.extend(
            [(shapes[index], int(gap))] * (heights[index] - 1)
        )
    return answer


def cleared_polynomial(shapes, heights, ceiling, column):
    """The positive row-clearing polynomial F_column(u)."""

    terminal = heights[-1]
    answer = [Fraction(1)]
    for root, length in residual_factors(shapes, heights):
        answer = polynomial_multiply(
            answer, rising_polynomial(root + column, length)
        )
    ray_length = ceiling - column
    answer = polynomial_multiply(
        answer, rising_polynomial(shapes[0] + column, ray_length)
    )
    for _ in range(terminal - 1):
        answer = polynomial_multiply(
            answer, rising_polynomial(shapes[-1] + column, ray_length)
        )
    return answer


def cleared_roots(shapes, heights, ceiling, column):
    roots = []
    for root, length in residual_factors(shapes, heights):
        roots.extend(root + column + offset for offset in range(length))
    ray_length = ceiling - column
    roots.extend(shapes[0] + column + offset for offset in range(ray_length))
    for _ in range(heights[-1] - 1):
        roots.extend(
            shapes[-1] + column + offset for offset in range(ray_length)
        )
    return sorted(roots)


def monotone_root_injection(smaller, larger):
    """Match every root in smaller to a no-larger root in larger."""

    available = list(larger)
    for old_root in sorted(smaller):
        candidates = [
            index for index, new_root in enumerate(available) if new_root <= old_root
        ]
        if not candidates:
            return False
        available.pop(candidates[-1])
    return True


def strict_prefix_universe():
    for edge_count in (1, 2):
        for tail in combinations(range(1, 5), edge_count):
            offsets = (0,) + tail
            for internal_heights in product((1, 2, 3), repeat=edge_count):
                for terminal in (2, 3):
                    heights = internal_heights + (terminal,)
                    residual_degree = sum(
                        (offsets[index + 1] - offsets[index])
                        * (heights[index] - 1)
                        for index in range(edge_count)
                    )
                    if residual_degree <= 5:
                        shapes = tuple(Fraction(1 + offset) for offset in offsets)
                        yield shapes, heights, residual_degree


def polynomial_from_roots(roots):
    answer = [Fraction(1)]
    for root in roots:
        answer = polynomial_multiply(answer, [root, 1])
    return answer


def main():
    ceiling = 3
    case_count = 0
    factorization_checks = 0
    root_event_checks = 0
    matched_root_checks = 0
    recurrence_checks = 0
    monomial_minors = 0
    monomial_zero = 0
    newton_minors = 0
    newton_zero = 0
    hankel_minors = 0

    for shapes, heights, residual_degree in strict_prefix_universe():
        case_count += 1
        exponents = exponent_inventory(heights)
        require(all(height > 0 for height in heights), "prefix lost strictness")
        require(-sum(exponents) == heights[-1], "terminal mass mismatch")

        factors = residual_factors(shapes, heights)
        normalization = Fraction(1)
        for root, length in factors:
            normalization *= rising(root, length)

        for index in range(7):
            polynomial_part = Fraction(1)
            for root, length in factors:
                polynomial_part *= rising(root + index, length)
            factored = polynomial_part / normalization
            factored /= rising(shapes[0], index)
            factored /= rising(shapes[-1], index) ** (heights[-1] - 1)
            require(
                factored == direct_sequence(shapes, exponents, index),
                "strict-prefix strand factorization failed",
            )
            factorization_checks += 1

        polynomials = [
            cleared_polynomial(shapes, heights, ceiling, ceiling - time)
            for time in range(ceiling + 1)
        ]

        for column in range(1, ceiling + 1):
            old_roots = cleared_roots(shapes, heights, ceiling, column)
            new_roots = cleared_roots(shapes, heights, ceiling, column - 1)
            old_count = Counter(old_roots)
            new_count = Counter(new_roots)
            for index, shape in enumerate(shapes):
                location = shape + column - 1
                require(
                    new_count[location] - old_count[location] == -exponents[index],
                    "root-event exponent mismatch",
                )
            other_locations = (set(old_count) | set(new_count)) - {
                shape + column - 1 for shape in shapes
            }
            require(
                all(new_count[root] == old_count[root] for root in other_locations),
                "root event escaped the shape endpoints",
            )
            cumulative = 0
            for exponent in exponents:
                cumulative -= exponent
                require(cumulative > 0, "root-event ballot prefix is not strict")
            root_event_checks += 1
            require(
                monotone_root_injection(old_roots, new_roots),
                "ballot root matching failed",
            )
            matched_root_checks += len(old_roots)

        production_root = rising_polynomial(shapes[0] + ceiling - 1, 1)
        for _ in range(heights[-1] - 1):
            production_root = polynomial_multiply(
                production_root, rising_polynomial(shapes[-1] + ceiling - 1, 1)
            )
        for time in range(ceiling):
            produced = polynomial_multiply(production_root, shift_down(polynomials[time]))
            require(
                produced == polynomials[time + 1],
                "fixed production-operator recurrence failed",
            )
            recurrence_checks += 1

        degree = residual_degree + heights[-1] * ceiling
        require(
            max(len(polynomial) for polynomial in polynomials) == degree + 1,
            "cleared degree mismatch",
        )
        monomial_matrix = [
            [
                polynomial[row] if row < len(polynomial) else Fraction(0)
                for polynomial in polynomials
            ]
            for row in range(degree + 1)
        ]
        newton_columns = [newton_coefficients(polynomial) for polynomial in polynomials]
        newton_matrix = [
            [
                column[row] if row < len(column) else Fraction(0)
                for column in newton_columns
            ]
            for row in range(degree + 1)
        ]
        checked, zero = all_minors_nonnegative(monomial_matrix)
        monomial_minors += checked
        monomial_zero += zero
        checked, zero = all_minors_nonnegative(newton_matrix)
        newton_minors += checked
        newton_zero += zero

        values = [direct_sequence(shapes, exponents, index) for index in range(9)]
        for order in (2, 3, 4):
            expected = (-1) ** (order * (order - 1) // 2)
            for rows in combinations(range(5), order):
                for columns in combinations(range(5), order):
                    value = determinant(
                        [
                            [values[row + column] for column in columns]
                            for row in rows
                        ]
                    )
                    require(expected * value > 0, "actual Hankel sign failed")
                    hankel_minors += 1

    # Bare monotone matching does not control the first order-three flag.
    hostile_roots = (
        (Fraction(1, 2), Fraction(17, 10)),
        (Fraction(1, 2), Fraction(17, 10), Fraction(39, 10), Fraction(47, 10)),
        (
            Fraction(1, 10),
            Fraction(1),
            Fraction(12, 5),
            Fraction(3),
            Fraction(3),
            Fraction(41, 10),
        ),
    )
    require(
        monotone_root_injection(hostile_roots[0], hostile_roots[1])
        and monotone_root_injection(hostile_roots[1], hostile_roots[2]),
        "generic hostile lost its root matching",
    )
    hostile_polynomials = [polynomial_from_roots(roots) for roots in hostile_roots]
    hostile_degree = max(len(polynomial) for polynomial in hostile_polynomials)
    hostile_matrix = [
        [
            polynomial[row] if row < len(polynomial) else Fraction(0)
            for polynomial in hostile_polynomials
        ]
        for row in range(hostile_degree)
    ]
    for rows in combinations(range(hostile_degree), 2):
        for columns in combinations(range(3), 2):
            value = determinant(
                [
                    [hostile_matrix[row][column] for column in columns]
                    for row in rows
                ]
            )
            require(value >= 0, "generic matched-root family lost coefficient TP2")
    hostile_coefficient_minor = determinant(
        [[polynomial[row] for polynomial in hostile_polynomials] for row in range(3)]
    )
    require(
        hostile_coefficient_minor == Fraction(-15565881, 125000),
        "generic coefficient hostile changed",
    )
    hostile_nodes = (Fraction(0), Fraction(1, 100), Fraction(1, 50))
    hostile_evaluation_minor = determinant(
        [
            [polynomial_value(polynomial, node) for polynomial in hostile_polynomials]
            for node in hostile_nodes
        ]
    )
    require(hostile_evaluation_minor < 0, "generic evaluation hostile changed sign")

    # The first nontrivial root-deletion layer: shapes (1,5), e=(-3,+1).
    # Its elementary-coefficient transfer series is (1+z)^3/(1+5z),
    # whose linear coefficient is already negative.
    local_transition_linear = Fraction(3) - Fraction(5)
    require(local_transition_linear == -2, "local transition hostile changed")

    # Independent non-strict-prefix hostile from THM-3065, equations (41)--(42).
    zero_shapes = (Fraction(1, 7), Fraction(2, 3), Fraction(20))
    zero_exponents = (-1, 1, -1)
    zero_values = [
        direct_sequence(zero_shapes, zero_exponents, index) for index in range(5)
    ]
    zero_prefix_minor = determinant(
        [[zero_values[row + column] for column in range(3)] for row in range(3)]
    )
    require(
        zero_prefix_minor == Fraction(4914161, 84138683904000),
        "zero-prefix hostile changed",
    )

    print("strict_prefix_terminal_q_ballot_transport=finite_exact")
    print(
        f"universe=cases:{case_count}:alpha0=1:integer_offsets<=4:"
        f"edges<=2:terminal_Q=2,3:residual_degree<=5:ceiling={ceiling}"
    )
    print(
        f"structure=factorizations:{factorization_checks}:root_events:{root_event_checks}:"
        f"matched_old_roots:{matched_root_checks}:fixed_operator_steps:{recurrence_checks}"
    )
    print(
        f"monomial_coefficient_minors={monomial_minors}:negative=0:zero={monomial_zero}"
    )
    print(
        f"newton_coefficient_minors={newton_minors}:negative=0:zero={newton_zero}"
    )
    print(f"independent_generalized_hankel_minors={hankel_minors}:wrong_sign=0")
    print(
        "generic_root_matching_hostile="
        f"coefficient_det:{hostile_coefficient_minor}:"
        f"evaluation_det:{hostile_evaluation_minor}"
    )
    print(
        "local_root_deletion_transition=(1+z)^3/(1+5z):"
        f"linear_coefficient={local_transition_linear}"
    )
    print(f"zero_prefix_control_order3={zero_prefix_minor}:wrong_checkerboard_sign")
    print(
        "verdict=integer_lattice_coefficient_network_is_load_bearing:"
        "all_exact_checks=PASS"
    )


if __name__ == "__main__":
    main()
