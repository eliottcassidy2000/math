#!/usr/bin/env python3
"""Exact companion for THM-3372.

Audits the multiaffine deletion transform against the odd-cycle-packing
formula, every Taylor/deletion identity, the sharp variance certificate, and
the skew ordered-join current.  All arithmetic is integral and the file uses
only the Python standard library.
"""

from __future__ import annotations

import hashlib
import itertools
import math
from collections import defaultdict


EXPECTED_SEMANTIC_SHA256 = "c66025a0fa9e9e36bc8f9e2576f95998ff6c3e6035eaa871f97c8997dbea7688"


def trim(polynomial):
    values = list(polynomial)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def polynomial_add(left, right):
    size = max(len(left), len(right))
    return trim(tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ))


def polynomial_scale(polynomial, scalar):
    return trim(tuple(scalar * coefficient for coefficient in polynomial))


def polynomial_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            result[left_degree + right_degree] += (
                left_coefficient * right_coefficient
            )
    return trim(result)


def polynomial_derivative(polynomial, order=1):
    result = tuple(polynomial)
    for _step in range(order):
        result = tuple(
            degree * result[degree] for degree in range(1, len(result))
        ) or (0,)
    return trim(result)


def shifted_monomial_power(power):
    return tuple(
        math.comb(power, degree) * (-1) ** (power - degree)
        for degree in range(power + 1)
    )


def bivariate_add(*terms):
    result = defaultdict(int)
    for term in terms:
        for monomial, coefficient in term.items():
            result[monomial] += coefficient
    return {
        monomial: coefficient
        for monomial, coefficient in sorted(result.items())
        if coefficient
    }


def bivariate_scale(polynomial, scalar):
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in polynomial.items()
        if scalar * coefficient
    }


def outer_product(left, right, scalar=1):
    return {
        (left_degree, right_degree): scalar * left_coefficient * right_coefficient
        for left_degree, left_coefficient in enumerate(left)
        for right_degree, right_coefficient in enumerate(right)
        if scalar * left_coefficient * right_coefficient
    }


def bivariate_times_separated(bivariate, left, right):
    result = defaultdict(int)
    for (z_degree, w_degree), coefficient in bivariate.items():
        for left_degree, left_coefficient in enumerate(left):
            for right_degree, right_coefficient in enumerate(right):
                result[(z_degree + left_degree, w_degree + right_degree)] += (
                    coefficient * left_coefficient * right_coefficient
                )
    return {
        monomial: coefficient
        for monomial, coefficient in sorted(result.items())
        if coefficient
    }


def all_tournaments(order):
    pairs = tuple(itertools.combinations(range(order), 2))
    for bits in range(1 << len(pairs)):
        out = [0] * order
        for index, (left, right) in enumerate(pairs):
            if bits & (1 << index):
                out[left] |= 1 << right
            else:
                out[right] |= 1 << left
        yield tuple(out)


def converse(out):
    full = (1 << len(out)) - 1
    return tuple(full ^ (1 << vertex) ^ out[vertex]
                 for vertex in range(len(out)))


def canonical_code(out):
    order = len(out)
    return min(
        tuple(
            1 if out[permutation[left]] & (1 << permutation[right]) else 0
            for left, right in itertools.combinations(range(order), 2)
        )
        for permutation in itertools.permutations(range(order))
    )


def ordered_join(left, right):
    left_order = len(left)
    right_order = len(right)
    right_vertices = ((1 << right_order) - 1) << left_order
    return tuple(
        left[vertex] | right_vertices for vertex in range(left_order)
    ) + tuple(
        right[vertex] << left_order for vertex in range(right_order)
    )


def hamiltonian_counts(out):
    order = len(out)
    full = (1 << order) - 1
    paths = [[0] * order for _mask in range(full + 1)]
    for vertex in range(order):
        paths[1 << vertex][vertex] = 1
    for mask in range(1, full + 1):
        for endpoint, count in enumerate(paths[mask]):
            if not count:
                continue
            available = out[endpoint] & (full ^ mask)
            while available:
                bit = available & -available
                available -= bit
                successor = bit.bit_length() - 1
                paths[mask | bit][successor] += count
    counts = [1]
    counts.extend(sum(paths[mask]) for mask in range(1, full + 1))
    return tuple(counts)


def diagonal_deletion_transform(universe, hamiltonian):
    result = (0,)
    deleted = universe
    while True:
        remaining = universe ^ deleted
        term = polynomial_scale(
            shifted_monomial_power(deleted.bit_count()),
            hamiltonian[remaining],
        )
        result = polynomial_add(result, term)
        if deleted == 0:
            break
        deleted = (deleted - 1) & universe
    return result


def directed_odd_cycles(out):
    order = len(out)
    cycles = []
    for length in range(3, order + 1, 2):
        for support in itertools.combinations(range(order), length):
            anchor = support[0]
            for tail in itertools.permutations(support[1:]):
                cycle = (anchor,) + tail
                if all(
                    out[cycle[index]] & (1 << cycle[(index + 1) % length])
                    for index in range(length)
                ):
                    cycles.append(sum(1 << vertex for vertex in support))
    return tuple(cycles)


def packing_coefficients(out):
    full = (1 << len(out)) - 1
    states = {0: 1}
    for cycle_mask in directed_odd_cycles(out):
        additions = defaultdict(int)
        for used, weight in tuple(states.items()):
            if not used & cycle_mask:
                additions[used | cycle_mask] += 2 * weight
        for used, weight in additions.items():
            states[used] = states.get(used, 0) + weight
    return {
        full ^ used: weight for used, weight in states.items()
    }


def multiaffine_from_deletions(out, hamiltonian):
    full = (1 << len(out)) - 1
    coefficients = defaultdict(int)
    deleted = full
    while True:
        value = hamiltonian[full ^ deleted]
        monomial = deleted
        while True:
            coefficients[monomial] += (
                value * (-1) ** (deleted.bit_count() - monomial.bit_count())
            )
            if monomial == 0:
                break
            monomial = (monomial - 1) & deleted
        if deleted == 0:
            break
        deleted = (deleted - 1) & full
    return {
        monomial: coefficient
        for monomial, coefficient in coefficients.items()
        if coefficient
    }


def diagonal_from_multiaffine(coefficients):
    degree = max((monomial.bit_count() for monomial in coefficients), default=0)
    result = [0] * (degree + 1)
    for monomial, coefficient in coefficients.items():
        result[monomial.bit_count()] += coefficient
    return trim(result)


def skew_current(out, deletion_polynomials):
    terms = []
    for left in range(len(out)):
        for right in range(len(out)):
            if left == right:
                continue
            sign = 1 if out[left] & (1 << right) else -1
            terms.append(outer_product(
                deletion_polynomials[left], deletion_polynomials[right], sign
            ))
    return bivariate_add(*terms)


def tournament_invariants(out):
    order = len(out)
    full = (1 << order) - 1
    hamiltonian = hamiltonian_counts(out)
    diagonal = diagonal_deletion_transform(full, hamiltonian)
    deletions = tuple(
        diagonal_deletion_transform(full ^ (1 << vertex), hamiltonian)
        for vertex in range(order)
    )
    return hamiltonian, diagonal, deletions, skew_current(out, deletions)


def has_directed_triangle(out):
    for left, middle, right in itertools.combinations(range(len(out)), 3):
        if (
            out[left] & (1 << middle)
            and out[middle] & (1 << right)
            and out[right] & (1 << left)
        ) or (
            out[left] & (1 << right)
            and out[right] & (1 << middle)
            and out[middle] & (1 << left)
        ):
            return True
    return False


def audit_tournaments():
    labelled = derivative_cells = variance_zero = 0
    nonzero_current_by_order = []
    collision_by_order = []
    order_four_types = defaultdict(int)
    digest = hashlib.sha256()
    banks = {}
    for order in range(1, 6):
        banks[order] = tuple(all_tournaments(order))
        nonzero_current = 0
        signature_classes = defaultdict(set)
        for out in banks[order]:
            labelled += 1
            full = (1 << order) - 1
            hamiltonian, diagonal, _deletions, current = tournament_invariants(out)
            multiaffine = multiaffine_from_deletions(out, hamiltonian)
            packing = packing_coefficients(out)
            assert multiaffine == packing
            assert diagonal_from_multiaffine(multiaffine) == diagonal

            for derivative_order in range(order + 1):
                deletion_sum = (0,)
                for deleted_vertices in itertools.combinations(
                    range(order), derivative_order
                ):
                    deleted = sum(1 << vertex for vertex in deleted_vertices)
                    deletion_sum = polynomial_add(
                        deletion_sum,
                        diagonal_deletion_transform(full ^ deleted, hamiltonian),
                    )
                assert polynomial_derivative(diagonal, derivative_order) == (
                    polynomial_scale(deletion_sum, math.factorial(derivative_order))
                )
                derivative_cells += 1

            base = hamiltonian[full]
            first = sum(hamiltonian[full ^ (1 << vertex)]
                        for vertex in range(order))
            second = sum(
                hamiltonian[full ^ (1 << left) ^ (1 << right)]
                for left, right in itertools.combinations(range(order), 2)
            )
            delta = base * (2 * second + first) - first * first
            assert delta >= 0
            assert (delta == 0) == (not has_directed_triangle(out))
            variance_zero += delta == 0

            opposite = tournament_invariants(converse(out))
            assert opposite[1] == diagonal
            assert opposite[3] == bivariate_scale(current, -1)
            nonzero_current += bool(current)
            if order == 4:
                order_four_types[(diagonal, delta, bool(current))] += 1
            signature_classes[
                (diagonal, tuple(sorted(current.items())))
            ].add(canonical_code(out))
            digest.update(repr((out, multiaffine, diagonal, delta, current)).encode())
        nonzero_current_by_order.append(nonzero_current)
        collision_by_order.append(any(
            len(classes) > 1 for classes in signature_classes.values()
        ))
    assert nonzero_current_by_order == [0, 0, 0, 16, 320]
    assert collision_by_order == [False, False, False, False, True]
    assert variance_zero == sum(math.factorial(order) for order in range(1, 6))
    expected_order_four = {
        ((0, 0, 0, 0, 1), 0, False): 24,
        ((0, 2, 0, 0, 1), 18, True): 16,
        ((0, 4, 0, 0, 1), 36, False): 24,
    }
    assert dict(order_four_types) == expected_order_four
    return (
        banks,
        labelled,
        derivative_cells,
        variance_zero,
        tuple(nonzero_current_by_order),
        tuple(collision_by_order),
        digest,
    )


def audit_ordered_joins(banks, digest):
    pair_cells = 0
    for left_order in range(1, 6):
        for right_order in range(1, 7 - left_order):
            if right_order not in banks:
                banks[right_order] = tuple(all_tournaments(right_order))
            for left in banks[left_order]:
                left_data = tournament_invariants(left)
                for right in banks[right_order]:
                    right_data = tournament_invariants(right)
                    joined = ordered_join(left, right)
                    joined_data = tournament_invariants(joined)
                    left_diagonal = left_data[1]
                    right_diagonal = right_data[1]
                    assert joined_data[1] == polynomial_multiply(
                        left_diagonal, right_diagonal
                    )
                    internal_left = bivariate_times_separated(
                        left_data[3], right_diagonal, right_diagonal
                    )
                    internal_right = bivariate_times_separated(
                        right_data[3], left_diagonal, left_diagonal
                    )
                    left_marked = polynomial_multiply(
                        right_diagonal, polynomial_derivative(left_diagonal)
                    )
                    right_marked = polynomial_multiply(
                        left_diagonal, polynomial_derivative(right_diagonal)
                    )
                    predicted = bivariate_add(
                        internal_left,
                        internal_right,
                        outer_product(left_marked, right_marked),
                        outer_product(right_marked, left_marked, -1),
                    )
                    assert joined_data[3] == predicted
                    pair_cells += 1
                    digest.update(repr((left, right, joined_data[1], predicted)).encode())
    assert pair_cells == 2553
    return pair_cells


def audit_sharp_hostile():
    singleton = (0,)
    cycle_three = (1 << 1, 1 << 2, 1 << 0)
    forward = tournament_invariants(ordered_join(singleton, cycle_three))
    reverse = tournament_invariants(ordered_join(cycle_three, singleton))
    expected_diagonal = (0, 2, 0, 0, 1)
    expected_current = {(0, 3): 6, (3, 0): -6}
    assert forward[1] == reverse[1] == expected_diagonal
    assert forward[3] == expected_current
    assert reverse[3] == bivariate_scale(expected_current, -1)
    return expected_diagonal, expected_current


def audit_loss_hostile():
    order_five = tuple(all_tournaments(5))
    left = order_five[8]
    right = order_five[10]
    left_data = tournament_invariants(left)
    right_data = tournament_invariants(right)
    common_diagonal = (2, 0, 6, 0, 0, 1)
    assert left_data[1] == right_data[1] == common_diagonal
    assert left_data[3] == right_data[3] == {}
    left_deletion_types = tuple(sorted(left_data[2]))
    right_deletion_types = tuple(sorted(right_data[2]))
    assert left_deletion_types != right_deletion_types
    assert canonical_code(left) != canonical_code(right)
    assert canonical_code(left) == canonical_code(converse(left))
    assert canonical_code(right) == canonical_code(converse(right))
    return common_diagonal, left_deletion_types, right_deletion_types


def main():
    print("THM-3372 MULTIAFFINE DELETION TRANSFORM AND SKEW CURRENT")
    (
        banks,
        labelled,
        derivative_cells,
        variance_zero,
        nonzero_current,
        collision_by_order,
        digest,
    ) = audit_tournaments()
    print(
        f"labelled_tournaments={labelled};multiaffine_OCF=PASS;"
        f"taylor_deletion_cells={derivative_cells}"
    )
    print(
        f"variance_cells={labelled};zero_exactly_transitive={variance_zero};"
        "sharp_C3_delta=18"
    )
    print(
        "order4_D_types=y4:24,y4+2y:16,y4+4y:24;"
        "middle_type_xi_sign_split=8+8"
    )
    print(
        "xi_nonzero_labelled_orders_1_to_5="
        + ",".join(map(str, nonzero_current))
    )
    print(
        "D_xi_cross_iso_collision_orders_1_to_5="
        + ",".join("1" if value else "0" for value in collision_by_order)
    )
    pair_cells = audit_ordered_joins(banks, digest)
    print(f"ordered_join_pair_cells={pair_cells};D_product=PASS;xi_law=PASS")
    diagonal, current = audit_sharp_hostile()
    print(
        f"sharp_hostile_D={diagonal};"
        f"forward_xi={tuple(sorted(current.items()))};reverse=-forward"
    )
    loss_diagonal, left_types, right_types = audit_loss_hostile()
    print(
        f"loss_hostile_masks=8,10;common_D={loss_diagonal};common_xi=0;"
        f"marked_deletion_types_differ={left_types}!={right_types}"
    )
    semantic = digest.hexdigest()
    print(f"semantic_sha256={semantic}")
    print("truth_gates=integer_HP_DP+literal_multiaffine_OCF+all_derivatives")
    print("scope=orders_1_to_5;ordered_factor_total_at_most_6;no_injectivity_claim")
    print("all_exact_controls=PASS")
    assert semantic == EXPECTED_SEMANTIC_SHA256


if __name__ == "__main__":
    main()
