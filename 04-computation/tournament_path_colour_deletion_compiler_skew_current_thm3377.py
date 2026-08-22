#!/usr/bin/env python3
"""Exact companion for THM-3377.

Audits the full subset-deletion lift of THM-3166's path-colour polynomial,
its specialization to THM-3372, coefficientwise negative-colour reciprocity,
and the skew ordered-join current.  Arithmetic is integral.  The two imported
modules are the audited standard-library companions for THM-3166 and THM-3372.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import sys
from functools import lru_cache
from pathlib import Path


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))

import tournament_multiaffine_deletion_transform_skew_current_thm3372 as deck
import tournament_order_join_falling_factorial_transform_thm3166 as colour


EXPECTED_SEMANTIC_SHA256 = "4403a78a56ee1c57dcb0e64865318e4f09187b26fa6a83697f852308d8e612ba"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(polynomial):
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim(tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ))


def poly_scale(polynomial, scalar):
    return trim(tuple(scalar * value for value in polynomial))


def poly_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            result[left_degree + right_degree] += left_value * right_value
    return trim(result)


def poly_value(polynomial, value):
    result = 0
    for coefficient in reversed(polynomial):
        result = result * value + coefficient
    return result


def compiler_trim(compiler):
    values = [trim(polynomial) for polynomial in compiler]
    while len(values) > 1 and values[-1] == (0,):
        values.pop()
    return tuple(values)


def compiler_add(left, right):
    size = max(len(left), len(right))
    return compiler_trim(tuple(
        poly_add(
            left[index] if index < len(left) else (0,),
            right[index] if index < len(right) else (0,),
        )
        for index in range(size)
    ))


def compiler_multiply(left, right):
    result = [(0,)] * (len(left) + len(right) - 1)
    for left_degree, left_polynomial in enumerate(left):
        for right_degree, right_polynomial in enumerate(right):
            degree = left_degree + right_degree
            result[degree] = poly_add(
                result[degree],
                poly_multiply(left_polynomial, right_polynomial),
            )
    return compiler_trim(result)


def compiler_derivative(compiler):
    if len(compiler) == 1:
        return ((0,),)
    return compiler_trim(tuple(
        poly_scale(compiler[degree], degree)
        for degree in range(1, len(compiler))
    ))


def compiler_value_in_colour(compiler, value):
    return trim(tuple(poly_value(polynomial, value) for polynomial in compiler))


def compiler_flatten(compiler):
    return {
        (deletion_degree, colour_degree): coefficient
        for deletion_degree, polynomial in enumerate(compiler)
        for colour_degree, coefficient in enumerate(polynomial)
        if coefficient
    }


def current_add(*terms):
    result = {}
    for term, scalar in terms:
        for monomial, coefficient in term.items():
            result[monomial] = result.get(monomial, 0) + scalar * coefficient
    return {
        monomial: coefficient
        for monomial, coefficient in sorted(result.items())
        if coefficient
    }


def current_outer(left, right, scalar=1):
    left_flat = compiler_flatten(left)
    right_flat = compiler_flatten(right)
    return {
        (left_deletion, left_colour, right_deletion, right_colour):
            scalar * left_value * right_value
        for (left_deletion, left_colour), left_value in left_flat.items()
        for (right_deletion, right_colour), right_value in right_flat.items()
        if scalar * left_value * right_value
    }


def current_times_separated(current, left, right):
    left_flat = compiler_flatten(left)
    right_flat = compiler_flatten(right)
    result = {}
    for monomial, coefficient in current.items():
        a_degree, z_degree, b_degree, w_degree = monomial
        for (extra_a, extra_z), left_value in left_flat.items():
            for (extra_b, extra_w), right_value in right_flat.items():
                target = (
                    a_degree + extra_a,
                    z_degree + extra_z,
                    b_degree + extra_b,
                    w_degree + extra_w,
                )
                result[target] = result.get(target, 0) + (
                    coefficient * left_value * right_value
                )
    return {
        monomial: coefficient
        for monomial, coefficient in sorted(result.items())
        if coefficient
    }


def current_swap(current):
    return {
        (b_degree, w_degree, a_degree, z_degree): coefficient
        for (a_degree, z_degree, b_degree, w_degree), coefficient
        in current.items()
    }


def all_tournaments(order):
    yield from range(1 << math.comb(order, 2))


def converse_code(order, code):
    return ((1 << math.comb(order, 2)) - 1) ^ code


def relabel_code(order, code, permutation):
    result = 0
    for left in range(order):
        for right in range(left + 1, order):
            if colour.beats(order, code, permutation[left], permutation[right]):
                result |= 1 << colour.edge_index(order, left, right)
    return result


def canonical_code(order, code):
    return min(
        relabel_code(order, code, permutation)
        for permutation in itertools.permutations(range(order))
    )


def out_rows(order, code):
    return tuple(
        sum(
            1 << right
            for right in range(order)
            if left != right and colour.beats(order, code, left, right)
        )
        for left in range(order)
    )


@lru_cache(maxsize=None)
def path_colour_polynomial(order, code):
    if order == 0:
        require(code == 0, (order, code))
        return (1,)
    return trim(colour.power_polynomial(colour.path_cover_profile(order, code)))


@lru_cache(maxsize=None)
def deletion_compiler(order, code):
    """Return coefficients in a of F_T(a,t)=sum_X Q_(T-X)(t)a^|X|."""
    require(0 <= code < (1 << math.comb(order, 2)), (order, code))
    result = [(0,)] * (order + 1)
    full = (1 << order) - 1
    deleted = full
    while True:
        vertices = tuple(
            vertex for vertex in range(order)
            if not deleted & (1 << vertex)
        )
        induced = colour.induced_code(order, code, vertices) if vertices else 0
        degree = deleted.bit_count()
        result[degree] = poly_add(
            result[degree], path_colour_polynomial(len(vertices), induced)
        )
        if deleted == 0:
            break
        deleted = (deleted - 1) & full
    return compiler_trim(result)


@lru_cache(maxsize=None)
def marked_compilers(order, code):
    result = []
    for deleted_vertex in range(order):
        vertices = tuple(
            vertex for vertex in range(order) if vertex != deleted_vertex
        )
        induced = colour.induced_code(order, code, vertices) if vertices else 0
        result.append(deletion_compiler(order - 1, induced))
    return tuple(result)


@lru_cache(maxsize=None)
def skew_current(order, code):
    responses = marked_compilers(order, code)
    result = {}
    for left in range(order):
        for right in range(order):
            if left == right:
                continue
            sign = 1 if colour.beats(order, code, left, right) else -1
            result = current_add(
                (result, 1),
                (current_outer(responses[left], responses[right]), sign),
            )
    return result


def shift_one(polynomial):
    """Return coefficients after y=1+a."""
    result = [0] * len(polynomial)
    for degree, coefficient in enumerate(polynomial):
        for target in range(degree + 1):
            result[target] += coefficient * math.comb(degree, target)
    return trim(result)


def feed(digest, *items):
    digest.update(repr(items).encode("ascii"))
    digest.update(b"\n")


def audit_tournaments():
    digest = hashlib.sha256()
    labelled = 0
    derivative_cells = 0
    reciprocity_cells = 0
    relabel_cells = 0
    nonzero = []
    banks = {}

    for order in range(1, 6):
        banks[order] = tuple(all_tournaments(order))
        order_nonzero = 0
        for code in banks[order]:
            labelled += 1
            compiler = deletion_compiler(order, code)
            marked = marked_compilers(order, code)
            marked_sum = ((0,),)
            for response in marked:
                marked_sum = compiler_add(marked_sum, response)
            require(marked_sum == compiler_derivative(compiler),
                    ("derivative", order, code))
            derivative_cells += order

            out = out_rows(order, code)
            hamiltonian = deck.hamiltonian_counts(out)
            diagonal = deck.diagonal_deletion_transform((1 << order) - 1,
                                                         hamiltonian)
            require(compiler_value_in_colour(compiler, 1) == shift_one(diagonal),
                    ("D-specialization", order, code))

            current = skew_current(order, code)
            require(current_add((current, 1), (current_swap(current), 1)) == {},
                    ("skew", order, code))
            opposite = skew_current(order, converse_code(order, code))
            require(current_add((current, 1), (opposite, 1)) == {},
                    ("converse", order, code))
            if current:
                order_nonzero += 1

            permutations = (
                tuple(range(1, order)) + (0,),
                tuple(reversed(range(order))),
            ) if order >= 2 else ((0,),)
            for permutation in permutations:
                relabelled = relabel_code(order, code, permutation)
                require(deletion_compiler(order, relabelled) == compiler,
                        ("compiler relabel", order, code, permutation))
                require(skew_current(order, relabelled) == current,
                        ("current relabel", order, code, permutation))
                relabel_cells += 1

            full = (1 << order) - 1
            deleted = full
            while True:
                vertices = tuple(
                    vertex for vertex in range(order)
                    if not deleted & (1 << vertex)
                )
                induced = colour.induced_code(order, code, vertices) if vertices else 0
                remaining_order = len(vertices)
                polynomial = path_colour_polynomial(remaining_order, induced)
                for m in range(1, 6):
                    if remaining_order == 0:
                        positive = 1
                    else:
                        backward = colour.backward_distribution(remaining_order, induced)
                        positive = colour.negative_colour_value(backward, m)
                    require(
                        (-1) ** remaining_order * poly_value(polynomial, -m)
                        == positive >= 0,
                        ("negative colour", order, code, deleted, m),
                    )
                    reciprocity_cells += 1
                if deleted == 0:
                    break
                deleted = (deleted - 1) & full

            for m in range(1, 6):
                signed_layers = tuple(
                    (-1) ** (order + degree) * poly_value(polynomial, -m)
                    for degree, polynomial in enumerate(compiler)
                )
                require(all(value >= 0 for value in signed_layers),
                        ("compiler positivity", order, code, m, signed_layers))
            feed(digest, order, code, compiler, current)
        nonzero.append(order_nonzero)

    return (
        banks,
        digest,
        labelled,
        derivative_cells,
        reciprocity_cells,
        relabel_cells,
        tuple(nonzero),
    )


def audit_ordered_joins(banks, digest):
    pair_cells = 0
    for first_order in range(1, 6):
        for second_order in range(1, 7 - first_order):
            for first in banks[first_order]:
                for second in banks[second_order]:
                    joined_order = first_order + second_order
                    joined = colour.join_code(
                        first_order, first, second_order, second
                    )
                    first_compiler = deletion_compiler(first_order, first)
                    second_compiler = deletion_compiler(second_order, second)
                    joined_compiler = deletion_compiler(joined_order, joined)
                    require(
                        joined_compiler
                        == compiler_multiply(first_compiler, second_compiler),
                        ("compiler product", first_order, first,
                         second_order, second),
                    )

                    first_current = skew_current(first_order, first)
                    second_current = skew_current(second_order, second)
                    first_derivative = compiler_derivative(first_compiler)
                    second_derivative = compiler_derivative(second_compiler)
                    predicted = current_add(
                        (current_times_separated(
                            first_current, second_compiler, second_compiler
                        ), 1),
                        (current_times_separated(
                            second_current, first_compiler, first_compiler
                        ), 1),
                        (current_outer(
                            compiler_multiply(second_compiler, first_derivative),
                            compiler_multiply(first_compiler, second_derivative),
                        ), 1),
                        (current_outer(
                            compiler_multiply(first_compiler, second_derivative),
                            compiler_multiply(second_compiler, first_derivative),
                        ), -1),
                    )
                    actual = skew_current(joined_order, joined)
                    require(actual == predicted,
                            ("current product", first_order, first,
                             second_order, second))
                    pair_cells += 1
                    feed(digest, first_order, first, second_order, second,
                         joined_compiler, actual)
    return pair_cells


def scalar_current(current):
    return {
        (z_degree, w_degree): coefficient
        for (a_degree, z_degree, b_degree, w_degree), coefficient
        in current.items()
        if a_degree == b_degree == 0
    }


def audit_hostiles(digest):
    cycle = 5  # 0->1->2->0 in the THM-3166 pair-bit convention.
    forward = colour.join_code(1, 0, 3, cycle)
    reverse = colour.join_code(3, cycle, 1, 0)
    forward_compiler = deletion_compiler(4, forward)
    reverse_compiler = deletion_compiler(4, reverse)
    require(forward_compiler == reverse_compiler, "hostile compiler mismatch")
    forward_current = scalar_current(skew_current(4, forward))
    reverse_current = scalar_current(skew_current(4, reverse))
    expected = {(1, 3): 6, (3, 1): -6}
    require(forward_current == expected, (forward_current, expected))
    require(current_add(
        ({(0, z, 0, w): value for (z, w), value in forward_current.items()}, 1),
        ({(0, z, 0, w): value for (z, w), value in reverse_current.items()}, 1),
    ) == {}, "reverse hostile")

    sidecars = []
    for code in (16, 83):
        compiler = deletion_compiler(6, code)
        joined = colour.join_code(6, code, 1, 0)
        sidecars.append((
            path_colour_polynomial(6, code),
            compiler[1],
            canonical_code(6, code),
            colour.scc_components(6, code),
            scalar_current(skew_current(6, code)),
            scalar_current(skew_current(7, joined)),
        ))
    require(sidecars[0][0] == sidecars[1][0]
            == (0, 0, 8, 0, 8, 0, 1), sidecars)
    require(sidecars[0][1] == (0, 8, 0, 24, 0, 6), sidecars)
    require(sidecars[1][1] == (0, 4, 0, 24, 0, 6), sidecars)
    require(sidecars[0][2] != sidecars[1][2], sidecars)
    require(sidecars[0][3] == sidecars[1][3]
            == ((0, 1, 2, 3, 4, 5),), sidecars)
    require(sidecars[0][4] == sidecars[1][4] == {}, sidecars)
    require(sidecars[0][5] != sidecars[1][5], sidecars)
    feed(digest, forward_compiler, forward_current, sidecars)
    return forward_compiler, forward_current, tuple(sidecars)


def main():
    (
        banks,
        digest,
        labelled,
        derivative_cells,
        reciprocity_cells,
        relabel_cells,
        nonzero,
    ) = audit_tournaments()
    pair_cells = audit_ordered_joins(banks, digest)
    hostile_compiler, hostile_current, sidecars = audit_hostiles(digest)
    semantic = digest.hexdigest()

    print("THM-3377 PATH-COLOUR DELETION COMPILER AND SKEW CURRENT")
    print(
        f"labelled_tournaments={labelled};compiler_D_specialization=PASS;"
        f"marked_derivative_cells={derivative_cells}"
    )
    print(
        f"negative_colour_deletion_cells={reciprocity_cells};"
        "coefficientwise_positive_reciprocity=PASS"
    )
    print(
        f"relabel_cells={relabel_cells};converse_and_variable_skew=PASS;"
        "Psi_nonzero_labelled_orders_1_to_5=" + ",".join(map(str, nonzero))
    )
    print(
        f"ordered_join_pair_cells={pair_cells};compiler_product=PASS;"
        "full_Psi_law=PASS"
    )
    print(
        f"sharp_hostile_compiler={hostile_compiler};"
        f"scalar_Psi={tuple(sorted(hostile_current.items()))};reverse=-forward"
    )
    print(
        "sidecar_hostile_masks=16,83;common_Q="
        f"{sidecars[0][0]};q_left={sidecars[0][1]};q_right={sidecars[1][1]};"
        "common_scalar_Psi=0;K1_join_currents_differ=PASS;"
        "both_strong_nonisomorphic=PASS"
    )
    print(f"semantic_sha256={semantic}")
    print("truth_gates=THM3166_HP_cover_DP+THM3372_deletion_DP+direct_full_current")
    print(
        "scope=orders_1_to_5;ordered_factor_total_at_most_6;"
        "negative_colours_1_to_5;sidecar_hostile_order_7"
    )
    print("all_exact_controls=PASS")
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            (semantic, EXPECTED_SEMANTIC_SHA256))


if __name__ == "__main__":
    main()
