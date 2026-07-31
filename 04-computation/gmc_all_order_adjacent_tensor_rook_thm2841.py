#!/usr/bin/env python3
"""Exact referee for THM-2841 all-order rook positivity."""

from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from random import Random


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def derangement(number):
    return sum(
        (-1) ** index * factorial(number) // factorial(index)
        for index in range(number + 1)
    )


def elementary(values, degree):
    total = 0
    for selection in combinations_with_replacement(range(len(values)), degree):
        if len(set(selection)) != degree:
            continue
        term = 1
        for index in selection:
            term *= values[index]
        total += term
    return total


def elementary_subset(values, degree):
    if degree == 0:
        return 1
    total = 0
    for mask in range(1 << len(values)):
        if mask.bit_count() != degree:
            continue
        term = 1
        for index, value in enumerate(values):
            if (mask >> index) & 1:
                term *= value
        total += term
    return total


def rook_numerator(a_values):
    size = sum(a_values)
    return sum(
        (-1) ** degree
        * elementary_subset(a_values, degree)
        * factorial(size - degree)
        for degree in range(len(a_values) + 1)
    )


def tensor_formula(indices):
    a_values = tuple(index + 1 for index in indices)
    return Fraction(
        rook_numerator(a_values),
        prod(factorial(value) for value in a_values),
    )


def direct_tensor(indices):
    terms = {0: Fraction(1)}
    for index in indices:
        next_terms = {}
        for degree, coefficient in terms.items():
            next_terms[degree + index] = (
                next_terms.get(degree + index, Fraction(0))
                - coefficient / factorial(index)
            )
            next_terms[degree + index + 1] = (
                next_terms.get(degree + index + 1, Fraction(0))
                + coefficient / factorial(index + 1)
            )
        terms = next_terms
    return sum(coefficient * factorial(degree) for degree, coefficient in terms.items())


def forbidden_permutation_count(a_values):
    size = sum(a_values)
    blocks = []
    start = 0
    for block_size in a_values:
        blocks.append(set(range(start, start + block_size)))
        start += block_size
    count = 0
    for assignment in permutations(range(size)):
        if all(assignment[row] not in blocks[row] for row in range(len(a_values))):
            count += 1
    return count


def multinomial(indices):
    total = sum(indices)
    return factorial(total) // prod(factorial(index) for index in indices)


formula_cells = 0
positive_cells = 0
equality_cells = 0
minimum_ratio = None

for order in range(2, 8):
    for indices in combinations_with_replacement(range(7), order):
        value = tensor_formula(indices)
        require(value == direct_tensor(indices), "direct tensor formula")
        lower = derangement(order) * multinomial(indices)
        require(value >= lower, "derangement lower bound")
        require(value > 0, "strict tensor positivity")
        expected_equality = order == 2 or all(index == 0 for index in indices)
        require((value == lower) == expected_equality, "lower-bound equality")
        ratio = value / lower
        if minimum_ratio is None or ratio < minimum_ratio:
            minimum_ratio = ratio
        formula_cells += 1
        positive_cells += 1
        equality_cells += value == lower

permutation_cells = 0
for order in range(2, 6):
    for a_values in product(range(1, 5), repeat=order):
        if sum(a_values) > 8:
            continue
        require(
            rook_numerator(a_values) == forbidden_permutation_count(a_values),
            "forbidden-board permutation count",
        )
        permutation_cells += 1

rng = Random(2841)
random_cells = 0
for _ in range(2000):
    order = rng.randrange(2, 10)
    indices = tuple(rng.randrange(0, 30) for _ in range(order))
    value = tensor_formula(indices)
    lower = derangement(order) * multinomial(indices)
    require(value == direct_tensor(indices), "random direct tensor")
    require(value >= lower > 0, "random positivity")
    random_cells += 1

subfactorials = []
for order in range(2, 11):
    value = tensor_formula((0,) * order)
    require(value == derangement(order), "subfactorial boundary")
    subfactorials.append(value)

require(elementary((2, 3, 5), 2) == elementary_subset((2, 3, 5), 2), "elementary control")

print("THM-2841 ALL-ORDER ADJACENT TENSOR ROOK REFEREE")
print("status=PROOF-COMPLETE CANDIDATE; VERIFIED-EXACT")
print(f"formula_cells={formula_cells}")
print(f"positive_cells={positive_cells}")
print(f"lower_bound_equalities={equality_cells}")
print(f"minimum_ratio_to_derangement_bound={minimum_ratio}")
print(f"permutation_cells={permutation_cells}")
print(f"random_cells={random_cells}")
print("subfactorials_k2_to_k10=" + ",".join(str(value) for value in subfactorials))
print("all_exact_controls=PASS")
