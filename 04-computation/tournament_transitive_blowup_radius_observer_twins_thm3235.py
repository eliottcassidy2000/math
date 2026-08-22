#!/usr/bin/env python3
"""Exact companion for THM-3235.

All mathematical checks use explicit RuntimeError gates, so ``python -O``
executes the same verification path.  The only noninteger arithmetic is exact
``fractions.Fraction`` arithmetic in the general block-determinant control.
"""

from __future__ import annotations

import ast
from array import array
from fractions import Fraction
from itertools import permutations
from math import comb, factorial
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def exact_div(numerator, denominator):
    require(denominator != 0, "Bareiss division by zero")
    if isinstance(numerator, int) and isinstance(denominator, int):
        require(numerator % denominator == 0, "nonexact Bareiss division")
        return numerator // denominator
    return numerator / denominator


def bareiss(matrix):
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "determinant matrix not square")
    if n == 0:
        return 1
    if n == 1:
        return matrix[0][0]
    work = [list(row) for row in matrix]
    previous = 1
    sign = 1
    for pivot_index in range(n - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, n) if work[row][pivot_index] != 0),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                work[row][column] = exact_div(numerator, previous)
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def matrix_for_determinant(adjacency, variables):
    n = len(adjacency)
    require(len(variables) == n, "wrong variable count")
    return [
        [
            (1 if row == column else 0)
            - adjacency[row][column] * variables[column]
            for column in range(n)
        ]
        for row in range(n)
    ]


def determinant_value(adjacency, variables):
    return bareiss(matrix_for_determinant(adjacency, variables))


def solve_exact(matrix, right_hand_side):
    n = len(matrix)
    require(len(right_hand_side) == n, "linear-system dimension mismatch")
    augmented = [
        [Fraction(value) for value in matrix[row]] + [Fraction(right_hand_side[row])]
        for row in range(n)
    ]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if augmented[row][column] != 0),
            None,
        )
        require(pivot is not None, "singular exact linear system")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            multiplier = augmented[row][column]
            if multiplier != 0:
                augmented[row] = [
                    augmented[row][entry] - multiplier * augmented[column][entry]
                    for entry in range(n + 1)
                ]
    return [augmented[row][-1] for row in range(n)]


def walk_observer_value(adjacency, variables):
    system = matrix_for_determinant(adjacency, variables)
    solution = solve_exact(system, [1] * len(adjacency))
    return sum(Fraction(variables[index]) * solution[index] for index in range(len(adjacency)))


def adjacency_from_mask(order: int, mask: int):
    adjacency = [[0] * order for _ in range(order)]
    bit = 0
    for lower in range(order):
        for upper in range(lower + 1, order):
            if (mask >> bit) & 1:
                adjacency[lower][upper] = 1
            else:
                adjacency[upper][lower] = 1
            bit += 1
    return adjacency


def principal_determinant_coefficients(adjacency):
    n = len(adjacency)
    coefficients = [0] * (1 << n)
    coefficients[0] = 1
    for subset in range(1, 1 << n):
        vertices = [vertex for vertex in range(n) if (subset >> vertex) & 1]
        principal = [[adjacency[row][column] for column in vertices] for row in vertices]
        coefficients[subset] = (-1) ** len(vertices) * bareiss(principal)
    return tuple(coefficients)


def permute_subset(subset: int, permutation) -> int:
    image = 0
    for vertex, target in enumerate(permutation):
        if (subset >> vertex) & 1:
            image |= 1 << target
    return image


def determinant_symmetries(coefficients):
    n = (len(coefficients) - 1).bit_length()
    symmetries = []
    for permutation in permutations(range(n)):
        if all(
            coefficients[subset] == coefficients[permute_subset(subset, permutation)]
            for subset in range(1 << n)
        ):
            symmetries.append(permutation)
    return tuple(symmetries)


def group_vertex_orbits(group, order):
    unseen = set(range(order))
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {permutation[seed] for permutation in group}
        orbits.append(tuple(sorted(orbit)))
        unseen -= orbit
    return tuple(orbits)


def cyclic_subset_orbit(seed: int, generator):
    orbit = []
    current = seed
    while current not in orbit:
        orbit.append(current)
        current = permute_subset(current, generator)
    return tuple(orbit)


def subset(vertices_one_based):
    answer = 0
    for vertex in vertices_one_based:
        answer |= 1 << (vertex - 1)
    return answer


def univariate_coefficients(multivariate_coefficients, order):
    answer = [0] * (order + 1)
    for monomial, coefficient in enumerate(multivariate_coefficients):
        answer[monomial.bit_count()] += coefficient
    return tuple(answer)


def polynomial_product(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] += left_coefficient * right_coefficient
    return tuple(answer)


def isomorphic(first, second) -> bool:
    n = len(first)
    for permutation in permutations(range(n)):
        if all(
            first[row][column] == second[permutation[row]][permutation[column]]
            for row in range(n)
            for column in range(n)
        ):
            return True
    return False


def cycle_is_directed(adjacency, cycle_one_based) -> bool:
    cycle = [vertex - 1 for vertex in cycle_one_based]
    return all(
        adjacency[cycle[index]][cycle[(index + 1) % len(cycle)]] == 1
        for index in range(len(cycle))
    )


def backward_histogram(adjacency):
    n = len(adjacency)
    histogram = [0] * n
    for ordering in permutations(range(n)):
        backward = sum(
            adjacency[ordering[index + 1]][ordering[index]]
            for index in range(n - 1)
        )
        histogram[backward] += 1
    return tuple(histogram)


def ordered_path_cover_profile(histogram):
    n = len(histogram)
    profile = []
    for depth in range(1, n + 1):
        count = sum(
            histogram[backward] * comb(n - 1 - backward, depth - 1 - backward)
            for backward in range(depth)
        )
        profile.append(count)
    return tuple(profile)


def falling_factorial_polynomial(depth):
    polynomial = (1,)
    for root in range(depth):
        polynomial = polynomial_product(polynomial, (-root, 1))
    return polynomial


def path_colour_polynomial(unordered_profile):
    answer = [0] * (len(unordered_profile) + 1)
    for depth, count in enumerate(unordered_profile, start=1):
        basis = falling_factorial_polynomial(depth)
        for degree, coefficient in enumerate(basis):
            answer[degree] += count * coefficient
    return tuple(answer)


def transitive_tournament(order):
    return [[1 if row < column else 0 for column in range(order)] for row in range(order)]


def transitive_blowup(base, multiplicity):
    base_order = len(base)
    order = base_order * multiplicity
    adjacency = [[0] * order for _ in range(order)]
    for base_row in range(base_order):
        for local_row in range(multiplicity):
            row = base_row * multiplicity + local_row
            for base_column in range(base_order):
                for local_column in range(multiplicity):
                    column = base_column * multiplicity + local_column
                    if base_row == base_column:
                        adjacency[row][column] = int(local_row < local_column)
                    else:
                        adjacency[row][column] = base[base_row][base_column]
    return adjacency


def block_diagonal_law_control(base, factors, variable_blocks):
    require(len(base) == len(factors) == len(variable_blocks), "block-law arity mismatch")
    total_order = sum(len(factor) for factor in factors)
    substitution = [[0] * total_order for _ in range(total_order)]
    offsets = []
    cursor = 0
    for factor in factors:
        offsets.append(cursor)
        cursor += len(factor)
    for block_row, factor_row in enumerate(factors):
        for local_row in range(len(factor_row)):
            row = offsets[block_row] + local_row
            for block_column, factor_column in enumerate(factors):
                for local_column in range(len(factor_column)):
                    column = offsets[block_column] + local_column
                    if block_row == block_column:
                        substitution[row][column] = factor_row[local_row][local_column]
                    else:
                        substitution[row][column] = base[block_row][block_column]
    flat_variables = [value for block in variable_blocks for value in block]
    left = determinant_value(substitution, flat_variables)

    determinants = []
    observers = []
    for factor, variables in zip(factors, variable_blocks):
        determinants.append(determinant_value(factor, variables))
        observers.append(walk_observer_value(factor, variables))
    right = Fraction(1)
    for determinant in determinants:
        right *= determinant
    right *= determinant_value(base, observers)
    require(Fraction(left) == right, "general block determinant law failed")
    return Fraction(left), tuple(determinants), tuple(observers)


def uniform_determinant_controls(base, multiplicity, seed):
    lifted = transitive_blowup(base, multiplicity)
    variables = [((seed + 3 * index) % 7) + 1 for index in range(len(lifted))]
    block_variables = [
        variables[block * multiplicity : (block + 1) * multiplicity]
        for block in range(len(base))
    ]
    effective = []
    for block in block_variables:
        product = 1
        for value in block:
            product *= 1 + value
        effective.append(product - 1)
    left = determinant_value(lifted, variables)
    right = determinant_value(base, effective)
    require(left == right, "uniform transitive determinant control failed")
    return left


def hamiltonian_path_count(adjacency):
    n = len(adjacency)
    require(n <= 18, "companion Held--Karp memory guard exceeded")
    state_count = 1 << n
    dynamic = [array("Q", [0]) * state_count for _ in range(n)]
    incoming_masks = []
    for terminal in range(n):
        incoming = 0
        for predecessor in range(n):
            if adjacency[predecessor][terminal]:
                incoming |= 1 << predecessor
        incoming_masks.append(incoming)
        dynamic[terminal][1 << terminal] = 1
    for visited in range(1, state_count):
        terminals = visited
        while terminals:
            terminal_bit = terminals & -terminals
            terminal = terminal_bit.bit_length() - 1
            previous = visited ^ terminal_bit
            if previous:
                candidates = previous & incoming_masks[terminal]
                count = 0
                while candidates:
                    predecessor_bit = candidates & -candidates
                    predecessor = predecessor_bit.bit_length() - 1
                    count += dynamic[predecessor][previous]
                    candidates ^= predecessor_bit
                dynamic[terminal][visited] = count
            terminals ^= terminal_bit
    full = state_count - 1
    return sum(dynamic[terminal][full] for terminal in range(n))


def main() -> None:
    source = Path(__file__).read_text(encoding="utf-8")
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source)))
    require(assertion_nodes == 0, "optimized mode could remove an exact check")

    q408 = adjacency_from_mask(6, 408)
    q1332 = adjacency_from_mask(6, 1332)
    expected_q408 = (
        (0, 0, 0, 0, 1, 1),
        (1, 0, 0, 0, 1, 1),
        (1, 1, 0, 0, 0, 0),
        (1, 1, 1, 0, 0, 0),
        (0, 0, 1, 1, 0, 0),
        (0, 0, 1, 1, 1, 0),
    )
    expected_q1332 = (
        (0, 0, 0, 1, 0, 1),
        (1, 0, 1, 0, 0, 1),
        (1, 0, 0, 0, 1, 0),
        (0, 1, 1, 0, 0, 0),
        (1, 1, 0, 1, 0, 0),
        (0, 0, 1, 1, 1, 0),
    )
    require(tuple(map(tuple, q408)) == expected_q408, "mask 408 convention mismatch")
    require(tuple(map(tuple, q1332)) == expected_q1332, "mask 1332 convention mismatch")
    cycle_408 = (1, 5, 3, 2, 6, 4)
    cycle_1332 = (1, 4, 2, 6, 3, 5)
    require(cycle_is_directed(q408, cycle_408), "mask 408 strong-cycle witness failed")
    require(cycle_is_directed(q1332, cycle_1332), "mask 1332 strong-cycle witness failed")

    determinant_408 = principal_determinant_coefficients(q408)
    determinant_1332 = principal_determinant_coefficients(q1332)

    formula_408 = [0] * 64
    formula_408[0] = 1
    pair_choices = (
        (subset((1,)), subset((2,)), subset((1, 2))),
        (subset((3,)), subset((4,)), subset((3, 4))),
        (subset((5,)), subset((6,)), subset((5, 6))),
    )
    for first in pair_choices[0]:
        for second in pair_choices[1]:
            for third in pair_choices[2]:
                formula_408[first | second | third] = -1
    require(determinant_408 == tuple(formula_408), "mask 408 determinant factorization failed")

    generator = (1, 2, 5, 4, 0, 3)
    orbit_seeds_and_weights = (
        (subset((1, 2, 4)), 1),
        (subset((1, 3, 4)), 1),
        (subset((1, 2, 3, 4)), 1),
        (subset((2, 3, 4, 5)), 1),
        (subset((1, 2, 3, 4, 5)), 2),
        (subset((1, 2, 3, 4, 5, 6)), 4),
    )
    formula_1332 = [0] * 64
    formula_1332[0] = 1
    orbit_sizes = []
    for seed, weight in orbit_seeds_and_weights:
        orbit = cyclic_subset_orbit(seed, generator)
        orbit_sizes.append(len(orbit))
        for monomial in orbit:
            formula_1332[monomial] -= weight
    require(tuple(orbit_sizes) == (6, 2, 6, 3, 6, 1), "mask 1332 orbit sizes failed")
    require(determinant_1332 == tuple(formula_1332), "mask 1332 orbit determinant failed")

    symmetries_408 = determinant_symmetries(determinant_408)
    symmetries_1332 = determinant_symmetries(determinant_1332)
    orbits_408 = group_vertex_orbits(symmetries_408, 6)
    orbits_1332 = group_vertex_orbits(symmetries_1332, 6)
    require(len(symmetries_408) == 48, "mask 408 determinant-owner group order failed")
    require(len(symmetries_1332) == 6, "mask 1332 determinant-owner group order failed")
    require(orbits_408 == ((0, 1, 2, 3, 4, 5),), "mask 408 owner group not transitive")
    require(orbits_1332 == ((0, 1, 2, 3, 4, 5),), "mask 1332 owner group not transitive")
    require(not isomorphic(q408, q1332), "observer twins unexpectedly isomorphic")

    diagonal_408 = univariate_coefficients(determinant_408, 6)
    diagonal_1332 = univariate_coefficients(determinant_1332, 6)
    expected_diagonal_408 = (1, 0, 0, -8, -12, -6, -1)
    expected_diagonal_1332 = (1, 0, 0, -8, -9, -12, -4)
    require(diagonal_408 == expected_diagonal_408, "mask 408 diagonal determinant failed")
    require(diagonal_1332 == expected_diagonal_1332, "mask 1332 diagonal determinant failed")
    factor_408 = tuple(-coefficient for coefficient in polynomial_product(
        (-1, 2, 1), (1, 2, 5, 4, 1)
    ))
    squared_factor = polynomial_product((1, 1, 2), (1, 1, 2))
    factor_1332 = tuple(-coefficient for coefficient in polynomial_product((-1, 2, 1), squared_factor))
    require(factor_408 == expected_diagonal_408, "mask 408 exact factorization failed")
    require(factor_1332 == expected_diagonal_1332, "mask 1332 exact factorization failed")

    histogram_408 = backward_histogram(q408)
    histogram_1332 = backward_histogram(q1332)
    expected_histogram = (45, 117, 198, 198, 117, 45)
    require(histogram_408 == expected_histogram, "mask 408 backward histogram failed")
    require(histogram_1332 == expected_histogram, "mask 1332 backward histogram failed")
    ordered_408 = ordered_path_cover_profile(histogram_408)
    ordered_1332 = ordered_path_cover_profile(histogram_1332)
    expected_ordered = (45, 342, 1116, 1944, 1800, 720)
    require(ordered_408 == expected_ordered, "mask 408 ordered profile failed")
    require(ordered_1332 == expected_ordered, "mask 1332 ordered profile failed")
    unordered_408 = tuple(ordered_408[index] // factorial(index + 1) for index in range(6))
    unordered_1332 = tuple(ordered_1332[index] // factorial(index + 1) for index in range(6))
    expected_unordered = (45, 171, 186, 81, 15, 1)
    require(unordered_408 == expected_unordered, "mask 408 unordered profile failed")
    require(unordered_1332 == expected_unordered, "mask 1332 unordered profile failed")
    colour_408 = path_colour_polynomial(unordered_408)
    colour_1332 = path_colour_polynomial(unordered_1332)
    expected_colour = (0, 0, 28, 0, 16, 0, 1)
    require(colour_408 == expected_colour, "mask 408 path-colour polynomial failed")
    require(colour_1332 == expected_colour, "mask 1332 path-colour polynomial failed")

    c3 = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    q4 = (
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (1, 0, 0, 1),
        (1, 1, 0, 0),
    )
    require(
        isomorphic(transitive_blowup(c3, 2), q408),
        "mask 408 is not C3[T2,T2,T2]",
    )
    uniform_values = []
    for base in (c3, q4):
        for multiplicity in (2, 3):
            for seed in range(1, 5):
                uniform_values.append(uniform_determinant_controls(base, multiplicity, seed))
    require(len(uniform_values) == 16, "wrong uniform determinant control count")

    general_value, internal_determinants, internal_observers = block_diagonal_law_control(
        c3,
        (c3, transitive_tournament(2), q4),
        ((1, 2, 3), (2, 4), (1, 1, 2, 1)),
    )
    require(general_value.denominator == 1, "general block determinant did not cancel integrally")

    for base in (c3, q4):
        for multiplicity, repetition in ((2, 3), (3, 2)):
            iterated = transitive_blowup(transitive_blowup(base, multiplicity), repetition)
            direct = transitive_blowup(base, multiplicity * repetition)
            require(iterated == direct, "transitive substitution associativity failed")

    generic_variables = (1, 2, 3, 4, 5, 6)
    observer_408 = walk_observer_value(q408, generic_variables)
    observer_1332 = walk_observer_value(q1332, generic_variables)
    require(observer_408 != observer_1332, "generic walk observers failed to split twins")

    response_408 = tuple(
        hamiltonian_path_count(transitive_blowup(q408, repetition))
        for repetition in (1, 2, 3)
    )
    response_1332 = tuple(
        hamiltonian_path_count(transitive_blowup(q1332, repetition))
        for repetition in (1, 2, 3)
    )
    expected_response_408 = (45, 421425, 73948095585)
    expected_response_1332 = (45, 411513, 71347853961)
    require(response_408 == expected_response_408, "mask 408 lift response failed")
    require(response_1332 == expected_response_1332, "mask 1332 lift response failed")
    response_minor = response_408[0] * response_1332[1] - response_408[1] * response_1332[0]
    require(response_minor == -446040, "rank-two response minor failed")

    print("THM-3235 exact companion")
    print("assertion-independent AST Assert nodes:", assertion_nodes)
    print("mask convention: vertices 1..6; pairs lexicographic; bit 0=(1,2); bit 1 means i->j")
    print("strong Hamilton cycles: 408", cycle_408, "; 1332", cycle_1332)
    print("determinant orbit sizes for 1332:", tuple(orbit_sizes))
    print("determinant-owner group orders: 408", len(symmetries_408), "; 1332", len(symmetries_1332))
    display_orbits_408 = tuple(tuple(vertex + 1 for vertex in orbit) for orbit in orbits_408)
    display_orbits_1332 = tuple(tuple(vertex + 1 for vertex in orbit) for orbit in orbits_1332)
    print("determinant-owner vertex orbits: 408", display_orbits_408, "; 1332", display_orbits_1332)
    print("diagonal determinants: 408", diagonal_408, "; 1332", diagonal_1332)
    print("common positive factor: t^2+2t-1; alpha=sqrt(2)-1; rho=1+sqrt(2)")
    print("backward-adjacency histogram (both):", expected_histogram)
    print("ordered path-cover profile F (both):", expected_ordered)
    print("unordered path-cover profile p (both):", expected_unordered)
    print("path-colour ordinary coefficients [t^0..t^6] (both):", expected_colour)
    print("uniform transitive determinant evaluations checked:", len(uniform_values))
    print("C3[T2,T2,T2] isomorphism for mask 408: PASS")
    print("general block determinant value:", general_value)
    print("general block internal determinants:", internal_determinants)
    print("general block internal observers:", internal_observers)
    print("transitive substitution associativity controls: 4")
    print("generic walk observers at X=(1,2,3,4,5,6): 408", observer_408, "; 1332", observer_1332)
    print("uniform transitive lift Hamilton responses r=1,2,3:")
    print("  408:", response_408)
    print("  1332:", response_1332)
    print("response minor using r=1,2:", response_minor)
    print("all exact checks: PASS")


if __name__ == "__main__":
    main()
