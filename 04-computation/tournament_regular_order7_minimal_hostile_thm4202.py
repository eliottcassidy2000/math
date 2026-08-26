#!/usr/bin/env python3
"""Standalone exact audit of THM-4202's first regular-score hostile.

This maintained certificate imports no repository module.  It exhausts all
labelled tournaments of orders below seven to show that every regular one is
vertex-transitive, constructs a regular non-vertex-transitive order-seven
tournament by a directed-triangle reversal, and rebuilds its Hamilton count,
exposure capacities, rooted parity states, ordinal child, and symmetry-defect
ledger by literal enumeration.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction
from functools import cache


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def tournament_from_bits(order: int, bits: int) -> tuple[int, ...]:
    out = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            if bits & (1 << cursor):
                out[left] |= 1 << right
            else:
                out[right] |= 1 << left
            cursor += 1
    return tuple(out)


def pair_label(out: tuple[int, ...]) -> str:
    return "".join(
        "1" if out[left] & (1 << right) else "0"
        for left in range(len(out))
        for right in range(left + 1, len(out))
    )


def relabelled_pair_label(out: tuple[int, ...], permutation: tuple[int, ...]) -> str:
    return "".join(
        "1" if out[permutation[left]] & (1 << permutation[right]) else "0"
        for left in range(len(out))
        for right in range(left + 1, len(out))
    )


def canonical_label(out: tuple[int, ...]) -> str:
    return min(
        relabelled_pair_label(out, permutation)
        for permutation in itertools.permutations(range(len(out)))
    )


def automorphisms(out: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    target = pair_label(out)
    return tuple(
        permutation
        for permutation in itertools.permutations(range(len(out)))
        if relabelled_pair_label(out, permutation) == target
    )


def orbit_partition(
    order: int, automorphism_group: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(order))
    answer = []
    while unseen:
        seed = min(unseen)
        orbit = tuple(sorted({permutation[seed] for permutation in automorphism_group}))
        answer.append(orbit)
        unseen.difference_update(orbit)
    return tuple(answer)


def circulant(order: int, connection: set[int]) -> tuple[int, ...]:
    return tuple(
        sum(
            1 << right
            for right in range(order)
            if right != left and (right - left) % order in connection
        )
        for left in range(order)
    )


def is_directed_triangle(
    out: tuple[int, ...], triangle: tuple[int, int, int]
) -> bool:
    return all(
        sum(
            bool(out[vertex] & (1 << other))
            for other in triangle
            if other != vertex
        )
        == 1
        for vertex in triangle
    )


def reverse_triangle(
    out: tuple[int, ...], triangle: tuple[int, int, int]
) -> tuple[int, ...]:
    answer = list(out)
    for left, right in itertools.combinations(triangle, 2):
        if answer[left] & (1 << right):
            answer[left] ^= 1 << right
            answer[right] |= 1 << left
        else:
            answer[right] ^= 1 << left
            answer[left] |= 1 << right
    return tuple(answer)


def path_valid(out: tuple[int, ...], path: tuple[int, ...]) -> bool:
    return all(
        out[path[index]] & (1 << path[index + 1])
        for index in range(len(path) - 1)
    )


def gate(out: tuple[int, ...], capacities: tuple[tuple[int, ...], ...]) -> int:
    order = len(out)
    degrees = tuple(sum(row) for row in capacities)
    currents = tuple(
        sum(
            capacities[vertex][other]
            if out[vertex] & (1 << other)
            else -capacities[vertex][other]
            for other in range(order)
            if other != vertex
        )
        for vertex in range(order)
    )
    mass = sum(
        capacities[left][right]
        for left in range(order)
        for right in range(left + 1, order)
    )
    squares = sum(
        capacities[left][right] ** 2
        for left in range(order)
        for right in range(left + 1, order)
    )
    disjoint = (
        mass * mass + squares - sum(degree * degree for degree in degrees)
    ) // 2
    current = sum(degrees[index] * currents[index] for index in range(order))
    return disjoint + 2 * current


@cache
def exposed_data(
    out: tuple[int, ...],
) -> tuple[int, tuple[tuple[int, ...], ...], int, int]:
    """Build H, the exposure tensor, W, and G_+ from marked words."""
    order = len(out)
    capacities = [[0] * order for _ in range(order)]
    hamilton = 0
    for permutation in itertools.permutations(range(order)):
        bad = tuple(
            index
            for index in range(order - 1)
            if not out[permutation[index]] & (1 << permutation[index + 1])
        )
        if not bad:
            hamilton += 1
            marked_gaps = range(order - 1)
        elif len(bad) == 1:
            marked_gaps = bad
        else:
            continue
        for gap in marked_gaps:
            left, right = permutation[gap], permutation[gap + 1]
            capacities[left][right] += 1
            capacities[right][left] += 1
    tensor = tuple(tuple(row) for row in capacities)
    mass = sum(
        tensor[left][right]
        for left in range(order)
        for right in range(left + 1, order)
    )
    return hamilton, tensor, mass, gate(out, tensor)


@cache
def rooted_states(
    out: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Build U and V by simple paths and complement Hamilton paths."""
    order = len(out)
    full = (1 << order) - 1
    complement_hamilton = [0] * (1 << order)
    complement_hamilton[0] = 1
    for mask in range(1, 1 << order):
        vertices = tuple(
            vertex for vertex in range(order) if mask & (1 << vertex)
        )
        complement_hamilton[mask] = sum(
            path_valid(out, permutation)
            for permutation in itertools.permutations(vertices)
        )
    starts = [[0, 0] for _ in range(order)]
    ends = [[0, 0] for _ in range(order)]
    for length in range(1, order + 1):
        parity = length & 1
        for path in itertools.permutations(range(order), length):
            if not path_valid(out, path):
                continue
            mask = sum(1 << vertex for vertex in path)
            weight = complement_hamilton[full ^ mask]
            starts[path[0]][parity] += weight
            ends[path[-1]][parity] += weight
    return tuple(map(tuple, starts)), tuple(map(tuple, ends))


def ordinal_out(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    left_order = len(left)
    right_mask = ((1 << len(right)) - 1) << left_order
    return tuple(row | right_mask for row in left) + tuple(
        row << left_order for row in right
    )


def direct_remainder(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    left_h, _, _, left_gate = exposed_data(left)
    right_h, _, _, right_gate = exposed_data(right)
    _, _, _, child_gate = exposed_data(ordinal_out(left, right))
    return child_gate - right_h**2 * left_gate - left_h**2 * right_gate


def transferred_remainder(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    """Separate rank-two assembly, checked against the literal order-ten child."""
    left_h, left_cap, _, left_gate = exposed_data(left)
    right_h, right_cap, _, right_gate = exposed_data(right)
    left_starts, _ = rooted_states(left)
    _, right_ends = rooted_states(right)
    left_order, right_order = len(left), len(right)
    child_order = left_order + right_order
    capacities = [[0] * child_order for _ in range(child_order)]
    for left_vertex in range(left_order):
        for right_vertex in range(left_vertex + 1, left_order):
            value = right_h * left_cap[left_vertex][right_vertex]
            capacities[left_vertex][right_vertex] = value
            capacities[right_vertex][left_vertex] = value
    for left_vertex in range(right_order):
        for right_vertex in range(left_vertex + 1, right_order):
            value = left_h * right_cap[left_vertex][right_vertex]
            capacities[left_order + left_vertex][left_order + right_vertex] = value
            capacities[left_order + right_vertex][left_order + left_vertex] = value
    for left_vertex in range(left_order):
        for right_vertex in range(right_order):
            value = 2 * (
                left_starts[left_vertex][0] * right_ends[right_vertex][0]
                + left_starts[left_vertex][1] * right_ends[right_vertex][1]
            )
            capacities[left_vertex][left_order + right_vertex] = value
            capacities[left_order + right_vertex][left_vertex] = value
    tensor = tuple(map(tuple, capacities))
    child_gate = gate(ordinal_out(left, right), tensor)
    return child_gate - right_h**2 * left_gate - left_h**2 * right_gate


def uniform_formula(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
    m, n = len(left), len(right)
    h, _, w, _ = exposed_data(left)
    s, _, v, _ = exposed_data(right)
    k = w * v + w * s + h * v + 2 * h * s
    return (
        h * s * w * v
        + k * s * w * Fraction(m + 2, m)
        + k * h * v * Fraction(n - 6, n)
        + Fraction(k * k, 2)
        * (1 + Fraction(1, m * n) + Fraction(3, m) - Fraction(5, n))
    )


def main() -> None:
    digest = hashlib.sha256()

    minimality_rows = []
    for order in range(1, 7):
        regular = []
        for bits in range(1 << (order * (order - 1) // 2)):
            out = tournament_from_bits(order, bits)
            if len({row.bit_count() for row in out}) == 1:
                regular.append(out)
        classes = {canonical_label(out) for out in regular}
        every_vt = True
        for out in regular:
            group = automorphisms(out)
            if len({permutation[0] for permutation in group}) != order:
                every_vt = False
                break
        need(every_vt, f"regular non-VT tournament below seven at order {order}")
        row = (order, len(regular), len(classes), every_vt)
        minimality_rows.append(row)
        digest.update(("minimality|" + repr(row) + "\n").encode("ascii"))

    seed = circulant(7, {1, 2, 3})
    triangle = (0, 2, 4)
    need(is_directed_triangle(seed, triangle), "selected triangle is not directed")
    hostile = reverse_triangle(seed, triangle)
    score = tuple(row.bit_count() for row in hostile)
    need(score == (3,) * 7, "triangle reversal did not preserve regularity")
    group = automorphisms(hostile)
    orbits = orbit_partition(7, group)
    need(len(orbits) > 1, "order-seven hostile stayed vertex-transitive")

    hamilton, capacities, mass, _ = exposed_data(hostile)
    starts, ends = rooted_states(hostile)
    degrees = tuple(sum(row) for row in capacities)
    incoming = tuple(
        sum(
            capacities[vertex][other]
            for other in range(7)
            if hostile[other] & (1 << vertex)
        )
        for vertex in range(7)
    )

    cycle = circulant(3, {1})
    cycle_starts, _ = rooted_states(cycle)
    need(len(set(cycle_starts)) == 1, "C3 rooted states not uniform")
    cross = tuple(
        2
        * (
            cycle_starts[0][0] * ends[vertex][0]
            + cycle_starts[0][1] * ends[vertex][1]
        )
        for vertex in range(7)
    )
    mean_cross = Fraction(sum(cross), 7)
    mean_degree = Fraction(sum(degrees), 7)
    mean_incoming = Fraction(sum(incoming), 7)
    variance = sum((Fraction(value) - mean_cross) ** 2 for value in cross)
    covariance = sum(
        (Fraction(cross[vertex]) - mean_cross)
        * (
            Fraction(degrees[vertex])
            - mean_degree
            + 4 * (Fraction(incoming[vertex]) - mean_incoming)
        )
        for vertex in range(7)
    )
    cycle_h, _, _, _ = exposed_data(cycle)
    penalty = 3 * cycle_h * covariance + Fraction(3 * (5 * 3 - 1), 2) * variance

    transferred = transferred_remainder(cycle, hostile)
    direct = direct_remainder(cycle, hostile)
    false_uniform = uniform_formula(cycle, hostile)
    need(transferred == direct, "literal child and transferred remainder disagree")
    need(false_uniform.denominator == 1, "uniform formula is nonintegral")
    discrepancy = false_uniform - direct
    need(discrepancy == penalty == 3816, "variance-covariance ledger mismatch")

    hostile_row = (
        pair_label(hostile),
        score,
        len(group),
        orbits,
        hamilton,
        mass,
        starts,
        ends,
        degrees,
        incoming,
        cross,
        direct,
        false_uniform.numerator,
        variance,
        covariance,
        discrepancy,
    )
    digest.update(("hostile|" + repr(hostile_row) + "\n").encode("ascii"))

    print("REGULAR_ORDER7_MINIMAL_HOSTILE_ACCEPT")
    for row in minimality_rows:
        print(
            "below7",
            "order",
            row[0],
            "regular_labelled",
            row[1],
            "regular_classes",
            row[2],
            "all_vertex_transitive",
            row[3],
        )
    print("seed_connection", (1, 2, 3), "reversed_triangle", triangle)
    print("pair_label", hostile_row[0])
    print("score", hostile_row[1])
    print("automorphism_size", hostile_row[2], "orbits", hostile_row[3])
    print("hamilton", hostile_row[4], "mass", hostile_row[5])
    print("rooted_starts", hostile_row[6])
    print("rooted_ends", hostile_row[7])
    print("capacity_degrees", hostile_row[8])
    print("incoming_masses", hostile_row[9])
    print("cross_row", hostile_row[10])
    print(
        "remainder_literal",
        hostile_row[11],
        "remainder_transfer",
        transferred,
        "false_uniform",
        hostile_row[12],
    )
    print(
        "defect",
        "variance",
        hostile_row[13],
        "covariance",
        hostile_row[14],
        "penalty",
        hostile_row[15],
    )
    print("semantic_sha256", digest.hexdigest())


if __name__ == "__main__":
    main()
