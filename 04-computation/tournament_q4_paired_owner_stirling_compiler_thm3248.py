#!/usr/bin/env python3
"""Exact companion for THM-3248's Q4 paired-owner Stirling compiler."""

from fractions import Fraction
from itertools import permutations
from math import comb, factorial


Q4_ARCS = (
    (0, 1),
    (1, 2),
    (2, 0),
    (2, 3),
    (3, 0),
    (3, 1),
)
ZERO_EXPONENT = (0, 0, 0, 0)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def polynomial_add(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, 0) + coefficient
        if result[exponent] == 0:
            del result[exponent]
    return result


def polynomial_scale(polynomial, scalar):
    return {
        exponent: scalar * coefficient
        for exponent, coefficient in polynomial.items()
        if scalar * coefficient != 0
    }


def polynomial_multiply(left, right):
    result = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index] for index in range(4)
            )
            result[exponent] = (
                result.get(exponent, 0) + left_coefficient * right_coefficient
            )
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def polynomial_power(polynomial, exponent):
    result = {ZERO_EXPONENT: 1}
    base = dict(polynomial)
    remaining = exponent
    while remaining:
        if remaining % 2:
            result = polynomial_multiply(result, base)
        remaining //= 2
        if remaining:
            base = polynomial_multiply(base, base)
    return result


def variable(index):
    exponent = [0, 0, 0, 0]
    exponent[index] = 1
    return {tuple(exponent): 1}


def determinant(matrix):
    size = len(matrix)
    total = {}
    for permutation in permutations(range(size)):
        inversions = sum(
            1
            for left in range(size)
            for right in range(left + 1, size)
            if permutation[left] > permutation[right]
        )
        term = {ZERO_EXPONENT: -1 if inversions % 2 else 1}
        for row, column in enumerate(permutation):
            term = polynomial_multiply(term, matrix[row][column])
        total = polynomial_add(total, term)
    return total


def q4_adjacency():
    adjacency = [[False for _ in range(4)] for _ in range(4)]
    for source, target in Q4_ARCS:
        adjacency[source][target] = True
    for left in range(4):
        for right in range(left + 1, 4):
            require(
                adjacency[left][right] != adjacency[right][left],
                "Q4 arc list does not define a tournament",
            )
    return adjacency


def derive_walk_rational_function():
    adjacency = q4_adjacency()
    variables = [variable(index) for index in range(4)]
    identity_minus_ax = []
    determinant_lemma_matrix = []
    for row in range(4):
        matrix_row = []
        lemma_row = []
        for column in range(4):
            entry = {ZERO_EXPONENT: 1} if row == column else {}
            if adjacency[row][column]:
                entry = polynomial_add(entry, polynomial_scale(variables[column], -1))
            matrix_row.append(entry)
            lemma_row.append(polynomial_add(entry, variables[column]))
        identity_minus_ax.append(matrix_row)
        determinant_lemma_matrix.append(lemma_row)
    denominator = determinant(identity_minus_ax)
    numerator = polynomial_add(
        determinant(determinant_lemma_matrix), polynomial_scale(denominator, -1)
    )
    return denominator, numerator


def paired_factorization_polynomials():
    one = {ZERO_EXPONENT: 1}
    x1, x2, x3, x4 = (variable(index) for index in range(4))
    y = polynomial_add(
        polynomial_add(x1, x4), polynomial_multiply(x1, x4)
    )
    z = polynomial_add(
        polynomial_add(x2, x3), polynomial_multiply(x2, x3)
    )
    c = polynomial_multiply(x2, x3)
    cy = polynomial_multiply(c, y)
    denominator = polynomial_add(one, polynomial_scale(cy, -1))
    numerator = polynomial_add(y, z)
    numerator = polynomial_add(numerator, polynomial_multiply(y, z))
    numerator = polynomial_add(numerator, polynomial_scale(cy, 2))
    return denominator, numerator


def permute_polynomial(polynomial, permutation):
    result = {}
    for exponent, coefficient in polynomial.items():
        new_exponent = [0, 0, 0, 0]
        for old_index, new_index in enumerate(permutation):
            new_exponent[new_index] = exponent[old_index]
        result[tuple(new_exponent)] = coefficient
    return result


def owner_group(polynomial):
    return tuple(
        permutation
        for permutation in permutations(range(4))
        if permute_polynomial(polynomial, permutation) == polynomial
    )


def group_orbits(group):
    unseen = set(range(4))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {permutation[seed] for permutation in group}
        orbits.append(tuple(sorted(index + 1 for index in orbit)))
        unseen -= orbit
    return tuple(orbits)


def stirling_second_table(maximum):
    table = [[0 for _ in range(maximum + 1)] for _ in range(maximum + 1)]
    table[0][0] = 1
    for n in range(1, maximum + 1):
        for k in range(1, n + 1):
            table[n][k] = table[n - 1][k - 1] + k * table[n - 1][k]
    return table


def falling_surjection(table, n, k):
    if k < 0 or k > n:
        return 0
    return factorial(k) * table[n][k]


def hamiltonian_stirling_formula(r, table):
    total = 0
    for k in range(r + 1):
        a_k = falling_surjection(table, 2 * r, k)
        a_next = falling_surjection(table, 2 * r, k + 1)
        b_k = falling_surjection(table, r, k)
        b_next = falling_surjection(table, r, k + 1)
        total += a_next * b_k * b_k
        total += 2 * (a_k + a_next) * b_k * b_next
        total += (a_k + 3 * a_next) * b_next * b_next
    return total


def finite_difference_power(n, k, shift):
    return sum(
        (-1 if (k - j) % 2 else 1) * comb(k, j) * (shift + j) ** n
        for j in range(k + 1)
    )


def abstract_multiply(left, right):
    result = {}
    for (a1, b1, c1), coefficient1 in left.items():
        for (a2, b2, c2), coefficient2 in right.items():
            exponent = (a1 + a2, b1 + b2, c1 + c2)
            result[exponent] = result.get(exponent, 0) + coefficient1 * coefficient2
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def abstract_power(polynomial, exponent):
    result = {(0, 0, 0): 1}
    for _ in range(exponent):
        result = abstract_multiply(result, polynomial)
    return result


def path_cover_stirling_formula(r, depth):
    numerator = {
        (1, 0, 0): 1,
        (0, 1, 0): 1,
        (1, 1, 0): 1,
        (1, 0, 1): 2,
    }
    numerator_power = abstract_power(numerator, depth)
    total = 0
    for (a, b, c), coefficient in numerator_power.items():
        for k in range(2 * r + 1):
            paired_14 = finite_difference_power(2 * r, k + a, 0)
            if paired_14 == 0:
                continue
            paired_23 = 0
            for shift in range(b + 1):
                sign = -1 if (b - shift) % 2 else 1
                difference = finite_difference_power(r, k + c, shift)
                paired_23 += sign * comb(b, shift) * difference * difference
            total += (
                coefficient
                * comb(depth + k - 1, k)
                * paired_14
                * paired_23
            )
    return total


def transitive_blowup_adjacency(r):
    quotient = q4_adjacency()
    size = 4 * r
    adjacency = [[False for _ in range(size)] for _ in range(size)]
    for left in range(size):
        left_block, left_position = divmod(left, r)
        for right in range(size):
            if left == right:
                continue
            right_block, right_position = divmod(right, r)
            if left_block == right_block:
                adjacency[left][right] = left_position < right_position
            else:
                adjacency[left][right] = quotient[left_block][right_block]
    return adjacency


def held_karp_hamiltonian_count(adjacency):
    size = len(adjacency)
    full_mask = (1 << size) - 1
    dynamic = [[0 for _ in range(size)] for _ in range(1 << size)]
    for vertex in range(size):
        dynamic[1 << vertex][vertex] = 1
    for mask in range(1, full_mask + 1):
        for terminal in range(size):
            count = dynamic[mask][terminal]
            if count == 0:
                continue
            available = full_mask ^ mask
            while available:
                bit = available & -available
                successor = bit.bit_length() - 1
                if adjacency[terminal][successor]:
                    dynamic[mask | bit][successor] += count
                available -= bit
    return sum(dynamic[full_mask])


def backward_adjacency_histogram(adjacency):
    size = len(adjacency)
    histogram = [0 for _ in range(size)]
    for ordering in permutations(range(size)):
        backward = sum(
            1
            for index in range(size - 1)
            if not adjacency[ordering[index]][ordering[index + 1]]
        )
        histogram[backward] += 1
    require(sum(histogram) == factorial(size), "permutation universe was incomplete")
    return histogram


def path_cover_count_from_histogram(histogram, depth):
    size = len(histogram)
    total = 0
    for backward, multiplicity in enumerate(histogram):
        optional_cuts = depth - 1 - backward
        if 0 <= optional_cuts <= size - 1 - backward:
            total += multiplicity * comb(size - 1 - backward, optional_cuts)
    return total


def convolved_finite_difference_row(n, shift):
    values = [Fraction((shift + index) ** n, factorial(index)) for index in range(n + 2)]
    alternating = [
        Fraction(-1 if index % 2 else 1, factorial(index))
        for index in range(n + 2)
    ]
    convolution = []
    for k in range(n + 2):
        convolution.append(
            sum(values[j] * alternating[k - j] for j in range(k + 1))
        )
    return tuple(convolution[k] * factorial(k) for k in range(n + 2))


def polynomial_table(polynomial):
    return tuple(sorted((exponent, coefficient) for exponent, coefficient in polynomial.items()))


def one_based_group(group):
    return tuple(tuple(index + 1 for index in permutation) for permutation in group)


def main():
    derived_denominator, derived_numerator = derive_walk_rational_function()
    expected_denominator, expected_numerator = paired_factorization_polynomials()
    require(
        derived_denominator == expected_denominator,
        "matrix determinant did not reproduce 1-CY",
    )
    require(
        derived_numerator == expected_numerator,
        "matrix determinant lemma did not reproduce Y+Z+YZ+2CY",
    )

    denominator_owners = owner_group(derived_denominator)
    numerator_owners = owner_group(derived_numerator)
    require(denominator_owners == numerator_owners, "D and W have different owners")
    require(len(denominator_owners) == 4, "paired owner group does not have order four")
    orbits = group_orbits(denominator_owners)
    require(orbits == ((1, 4), (2, 3)), "owner group has unexpected coordinate orbits")

    table = stirling_second_table(20)
    hamiltonian_values = tuple(
        hamiltonian_stirling_formula(r, table) for r in range(1, 11)
    )
    held_karp_values = tuple(
        held_karp_hamiltonian_count(transitive_blowup_adjacency(r))
        for r in range(1, 5)
    )
    require(
        held_karp_values == hamiltonian_values[:4],
        "Held-Karp control disagreed with the Stirling contraction",
    )

    path_cover_controls = []
    for r in (1, 2):
        adjacency = transitive_blowup_adjacency(r)
        histogram = backward_adjacency_histogram(adjacency)
        for depth in range(1, 5):
            direct = path_cover_count_from_histogram(histogram, depth)
            compiled = path_cover_stirling_formula(r, depth)
            require(direct == compiled, "path-cover compiler disagreed with permutation control")
            path_cover_controls.append((r, depth, direct))

    transform_cells = 0
    for n in range(13):
        for shift in range(5):
            convolved = convolved_finite_difference_row(n, shift)
            for k, value in enumerate(convolved):
                require(value.denominator == 1, "binomial transform produced a noninteger")
                require(
                    value.numerator == finite_difference_power(n, k, shift),
                    "binomial-transform row disagreed with direct finite differences",
                )
                transform_cells += 1

    denominator_only_r1 = sum(
        finite_difference_power(2, k, 0)
        * finite_difference_power(1, k, 0)
        * finite_difference_power(1, k, 0)
        for k in range(3)
    )
    require(denominator_only_r1 == 1, "denominator-only hostile control changed")
    require(hamiltonian_values[0] == 5, "full Q4 Hamiltonian control changed")

    print("THM-3248 exact companion")
    print("determinant_terms=" + repr(polynomial_table(derived_denominator)))
    print("walk_numerator_terms=" + repr(polynomial_table(derived_numerator)))
    print("paired_factorization=D=1-CY;N=Y+Z+YZ+2CY")
    print("owner_group=" + repr(one_based_group(denominator_owners)))
    print("owner_orbits=" + repr(orbits))
    print("hamiltonian_r1_to_r10=" + repr(hamiltonian_values))
    print("held_karp_r1_to_r4=" + repr(held_karp_values))
    print("path_cover_controls=" + repr(tuple(path_cover_controls)))
    print("binomial_transform_cells=" + str(transform_cells))
    print(
        "hostile_denominator_only_r1="
        + str(denominator_only_r1)
        + ";full_walk_r1="
        + str(hamiltonian_values[0])
    )
    print("scope=fixed-depth unit-cost arithmetic; no C-finiteness or P-recursiveness claim")


if __name__ == "__main__":
    main()
