#!/usr/bin/env python3
"""Independent literal-permutation referee for THM-4202.

No tournament computation module is imported.  Factor Hamilton counts,
exposed-gap capacities, and rooted path states are rebuilt by literal
permutations.  Small ordinal children are rebuilt directly rather than via the
rank-two transfer.  The order-nine score-regular hostile uses an independently
implemented transfer only after its factor sidecars have been rebuilt.
"""

from __future__ import annotations

import hashlib
import itertools
from functools import cache
from fractions import Fraction


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def transitive(order: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << right for right in range(left + 1, order))
        for left in range(order)
    )


def circulant(order: int, connection: set[int]) -> tuple[int, ...]:
    return tuple(
        sum(
            1 << right
            for right in range(order)
            if right != left and (right - left) % order in connection
        )
        for left in range(order)
    )


def circulant_bank(order: int) -> tuple[tuple[int, ...], ...]:
    if order == 1:
        return (transitive(1),)
    pairs = tuple((step, order - step) for step in range(1, (order + 1) // 2))
    return tuple(
        circulant(order, {pairs[i][bit] for i, bit in enumerate(choice)})
        for choice in itertools.product((0, 1), repeat=len(pairs))
    )


def label(out: tuple[int, ...]) -> str:
    return "".join(
        "1" if out[left] & (1 << right) else "0"
        for left in range(len(out))
        for right in range(left + 1, len(out))
    )


def path_valid(out: tuple[int, ...], path: tuple[int, ...]) -> bool:
    return all(out[path[i]] & (1 << path[i + 1]) for i in range(len(path) - 1))


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
    disjoint = (mass * mass + squares - sum(value * value for value in degrees)) // 2
    current = sum(degrees[i] * currents[i] for i in range(order))
    return disjoint + 2 * current


@cache
def exposed_data(
    out: tuple[int, ...],
) -> tuple[int, tuple[tuple[int, ...], ...], int, int]:
    """Literal marked-word construction of H and c."""
    order = len(out)
    capacities = [[0] * order for _ in range(order)]
    hamilton = 0
    for permutation in itertools.permutations(range(order)):
        good = tuple(
            bool(out[permutation[i]] & (1 << permutation[i + 1]))
            for i in range(order - 1)
        )
        bad = tuple(i for i, value in enumerate(good) if not value)
        if not bad:
            hamilton += 1
            marked = range(order - 1)
        elif len(bad) == 1:
            marked = bad
        else:
            continue
        for gap in marked:
            left, right = permutation[gap], permutation[gap + 1]
            capacities[left][right] += 1
            capacities[right][left] += 1
    capacity_tuple = tuple(tuple(row) for row in capacities)
    mass = sum(
        capacities[left][right]
        for left in range(order)
        for right in range(left + 1, order)
    )
    return hamilton, capacity_tuple, mass, gate(out, capacity_tuple)


@cache
def rooted_states(
    out: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Literal simple-path/complement-H construction of U and V."""
    order = len(out)
    full = (1 << order) - 1
    subset_h = [0] * (1 << order)
    subset_h[0] = 1
    for mask in range(1, 1 << order):
        vertices = tuple(vertex for vertex in range(order) if mask & (1 << vertex))
        subset_h[mask] = sum(
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
            weight = subset_h[full ^ mask]
            starts[path[0]][parity] += weight
            ends[path[-1]][parity] += weight
    return tuple(tuple(row) for row in starts), tuple(tuple(row) for row in ends)


def ordinal_out(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    nleft = len(left)
    nright = len(right)
    right_mask = ((1 << nright) - 1) << nleft
    return tuple(row | right_mask for row in left) + tuple(
        row << nleft for row in right
    )


def direct_remainder(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    left_h, _, _, left_gate = exposed_data(left)
    right_h, _, _, right_gate = exposed_data(right)
    _, _, _, child_gate = exposed_data(ordinal_out(left, right))
    return child_gate - right_h**2 * left_gate - left_h**2 * right_gate


def transferred_remainder(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    """Independent rank-two assembly, used only for the order-12 hostile child."""
    left_h, left_cap, _, left_gate = exposed_data(left)
    right_h, right_cap, _, right_gate = exposed_data(right)
    left_starts, _ = rooted_states(left)
    _, right_ends = rooted_states(right)
    nleft, nright = len(left), len(right)
    order = nleft + nright
    capacities = [[0] * order for _ in range(order)]
    for i in range(nleft):
        for j in range(i + 1, nleft):
            value = right_h * left_cap[i][j]
            capacities[i][j] = capacities[j][i] = value
    for i in range(nright):
        for j in range(i + 1, nright):
            value = left_h * right_cap[i][j]
            capacities[nleft + i][nleft + j] = capacities[nleft + j][nleft + i] = value
    for i in range(nleft):
        for j in range(nright):
            value = 2 * (
                left_starts[i][0] * right_ends[j][0]
                + left_starts[i][1] * right_ends[j][1]
            )
            capacities[i][nleft + j] = capacities[nleft + j][i] = value
    capacity_tuple = tuple(tuple(row) for row in capacities)
    child_gate = gate(ordinal_out(left, right), capacity_tuple)
    return child_gate - right_h**2 * left_gate - left_h**2 * right_gate


def vt_formula(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
    m, n = len(left), len(right)
    h, _, w, _ = exposed_data(left)
    s, _, v, _ = exposed_data(right)
    cross_mass = w * v + w * s + h * v + 2 * h * s
    return (
        h * s * w * v
        + cross_mass * s * w * Fraction(m + 2, m)
        + cross_mass * h * v * Fraction(n - 6, n)
        + Fraction(cross_mass * cross_mass, 2)
        * (1 + Fraction(1, m * n) + Fraction(3, m) - Fraction(5, n))
    )


def incoming_masses(
    out: tuple[int, ...], capacities: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    order = len(out)
    return tuple(
        sum(capacities[i][j] for j in range(order) if out[j] & (1 << i))
        for i in range(order)
    )


def one_sided_vt_formula(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[Fraction, tuple[int, ...]]:
    m = len(left)
    h, _, w, _ = exposed_data(left)
    s, right_cap, v, _ = exposed_data(right)
    left_starts, _ = rooted_states(left)
    _, right_ends = rooted_states(right)
    cross = tuple(
        2
        * (
            left_starts[0][0] * right_ends[j][0]
            + left_starts[0][1] * right_ends[j][1]
        )
        for j in range(len(right))
    )
    total = sum(cross)
    squares = sum(value * value for value in cross)
    degrees = tuple(sum(row) for row in right_cap)
    incoming = incoming_masses(right, right_cap)
    linear = sum(
        cross[j] * (v - degrees[j] - 4 * incoming[j])
        for j in range(len(right))
    )
    value = (
        h * s * w * v
        + (m + 2) * s * w * total
        + m * h * linear
        + Fraction(m, 2)
        * ((m + 3) * total * total - (5 * m - 1) * squares)
    )
    return value, cross


def odd_symmetry_defect(
    left: tuple[int, ...], right: tuple[int, ...], cross: tuple[int, ...]
) -> tuple[Fraction, Fraction, Fraction]:
    m, n = len(left), len(right)
    h, _, _, _ = exposed_data(left)
    _, right_cap, _, _ = exposed_data(right)
    degrees = tuple(sum(row) for row in right_cap)
    incoming = incoming_masses(right, right_cap)
    mean_cross = Fraction(sum(cross), n)
    mean_degree = Fraction(sum(degrees), n)
    mean_incoming = Fraction(sum(incoming), n)
    variance = sum((Fraction(value) - mean_cross) ** 2 for value in cross)
    covariance = sum(
        (Fraction(cross[j]) - mean_cross)
        * (
            Fraction(degrees[j])
            - mean_degree
            + 4 * (Fraction(incoming[j]) - mean_incoming)
        )
        for j in range(n)
    )
    penalty = m * h * covariance + Fraction(m * (5 * m - 1), 2) * variance
    return variance, covariance, penalty


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


def phase_tournament(positions: tuple[Fraction, ...]) -> tuple[int, ...]:
    order = len(positions)
    out = [0] * order
    for left, right in itertools.combinations(range(order), 2):
        delta = (positions[right] - positions[left]) % 1
        winner = left if delta == 0 or delta == Fraction(1, 2) or delta < Fraction(1, 2) else right
        loser = right if winner == left else left
        out[winner] |= 1 << loser
    return tuple(out)


def canonical(out: tuple[int, ...], fixed_zero: bool) -> str:
    order = len(out)
    permutations = (
        ((0,) + tail for tail in itertools.permutations(range(1, order)))
        if fixed_zero
        else itertools.permutations(range(order))
    )
    return min(
        "".join(
            "1" if out[permutation[i]] & (1 << permutation[j]) else "0"
            for i in range(order)
            for j in range(order)
        )
        for permutation in permutations
    )


def main() -> None:
    digest = hashlib.sha256()
    factor_banks = {order: circulant_bank(order) for order in (1, 3, 5)}

    # Directly recover the uniform coordinates used by the proof.
    uniform_checks = 0
    for order, bank in factor_banks.items():
        for out in bank:
            hamilton, capacities, mass, _ = exposed_data(out)
            starts, ends = rooted_states(out)
            expected = (
                Fraction(mass, 2 * order),
                Fraction(mass + 2 * hamilton, 2 * order),
            )
            need(all(tuple(map(Fraction, row)) == expected for row in starts), "start uniformity")
            need(all(tuple(map(Fraction, row)) == expected for row in ends), "end uniformity")
            degrees = tuple(sum(row) for row in capacities)
            outgoing = tuple(
                sum(capacities[i][j] for j in range(order) if out[i] & (1 << j))
                for i in range(order)
            )
            incoming = tuple(
                sum(capacities[i][j] for j in range(order) if out[j] & (1 << i))
                for i in range(order)
            )
            need(len(set(degrees)) == len(set(outgoing)) == len(set(incoming)) == 1, "capacity uniformity")
            need(degrees[0] == Fraction(2 * mass, order), "degree value")
            need(outgoing[0] == incoming[0] == Fraction(mass, order), "current value")
            uniform_checks += 1

    lefts = tuple(out for bank in factor_banks.values() for out in bank)
    rights = factor_banks[3] + factor_banks[5]
    direct_checks = 0
    minimum = None
    for left in lefts:
        for right in rights:
            if len(left) + len(right) > 8:
                continue
            direct = direct_remainder(left, right)
            formula = vt_formula(left, right)
            need(formula.denominator == 1 and formula.numerator == direct, "direct formula mismatch")
            need(direct > 0, "direct nonpositive row")
            row = (direct, len(left), len(right), label(left), label(right))
            if minimum is None or row < minimum:
                minimum = row
            digest.update(("direct|" + "|".join(map(str, row)) + "\n").encode("ascii"))
            direct_checks += 1

    # One-sided identity on all eight labelled order-three right factors.
    cycle = circulant(3, {1})
    one_sided_checks = 0
    for bits in range(8):
        right_rows = [0] * 3
        cursor = 0
        for left_vertex in range(3):
            for right_vertex in range(left_vertex + 1, 3):
                if bits & (1 << cursor):
                    right_rows[left_vertex] |= 1 << right_vertex
                else:
                    right_rows[right_vertex] |= 1 << left_vertex
                cursor += 1
        right = tuple(right_rows)
        direct = direct_remainder(cycle, right)
        formula, cross = one_sided_vt_formula(cycle, right)
        need(formula.denominator == 1 and formula.numerator == direct, "one-sided direct formula")
        variance, covariance, penalty = odd_symmetry_defect(cycle, right, cross)
        need(vt_formula(cycle, right) - penalty == direct, "one-sided defect")
        one_sided_checks += 1

    # Rebuild the primary order-nine score-regular hostile literally.
    seed = circulant(9, {1, 2, 3, 4})
    triangle = (0, 2, 5)
    need(
        all(
            sum(bool(seed[v] & (1 << u)) for u in triangle if u != v) == 1
            for v in triangle
        ),
        "hostile seed triangle",
    )
    hostile = reverse_triangle(seed, triangle)
    need({row.bit_count() for row in hostile} == {4}, "hostile regular score")
    hostile_h, hostile_cap, _, _ = exposed_data(hostile)
    _, hostile_ends = rooted_states(hostile)
    hostile_degrees = tuple(sum(row) for row in hostile_cap)
    hostile_incoming = tuple(
        sum(hostile_cap[i][j] for j in range(9) if hostile[j] & (1 << i))
        for i in range(9)
    )
    need(len(set(hostile_ends)) > 1, "hostile rooted states stayed uniform")
    need(len(set(hostile_degrees)) > 1, "hostile degrees stayed uniform")
    need(len(set(hostile_incoming)) > 1, "hostile incoming masses stayed uniform")
    hostile_direct = transferred_remainder(cycle, hostile)
    false_uniform = vt_formula(cycle, hostile)
    need(false_uniform.denominator == 1 and hostile_direct != false_uniform, "hostile separation")
    hostile_cross = tuple(2 * (row[0] + 2 * row[1]) for row in hostile_ends)
    hostile_variance, hostile_covariance, hostile_penalty = odd_symmetry_defect(
        cycle, hostile, hostile_cross
    )
    need(hostile_penalty == false_uniform - hostile_direct, "hostile defect penalty")
    hostile_row = (
        triangle,
        label(hostile),
        hostile_h,
        hostile_ends,
        hostile_degrees,
        hostile_incoming,
        hostile_direct,
        false_uniform.numerator,
        false_uniform.numerator - hostile_direct,
    )
    digest.update(("hostile|" + repr(hostile_row) + "\n").encode("ascii"))

    # Independent N=4 LRC projection hostile.
    lrc_rows = []
    for time in (Fraction(41, 96), Fraction(7, 48)):
        positions = tuple((speed * time) % 1 for speed in (0, 1, 3, 4))
        distances = tuple(min(point, 1 - point) for point in positions[1:])
        safe = all(distance >= Fraction(1, 4) for distance in distances)
        tournament = phase_tournament(positions)
        lrc_rows.append(
            (
                time,
                positions,
                distances,
                safe,
                canonical(tournament, False),
                canonical(tournament, True),
            )
        )
    need(lrc_rows[0][3] and not lrc_rows[1][3], "independent LRC safe flip")
    need(lrc_rows[0][4:] == lrc_rows[1][4:], "independent LRC class mismatch")
    digest.update(("lrc|" + repr(lrc_rows) + "\n").encode("ascii"))

    print("VT_ORDINAL_REMAINDER_INDEPENDENT_ACCEPT")
    print("uniform_factor_checks", uniform_checks)
    print("direct_literal_child_checks", direct_checks)
    print("direct_minimum", minimum)
    print("one_sided_odd_defect_checks", one_sided_checks)
    print("regular_hostile_label", hostile_row[1])
    print("regular_hostile_hamilton", hostile_row[2])
    print("regular_hostile_end_states", hostile_row[3])
    print("regular_hostile_capacity_degrees", hostile_row[4])
    print("regular_hostile_incoming_masses", hostile_row[5])
    print(
        "regular_hostile_remainder",
        hostile_row[6],
        "false_uniform",
        hostile_row[7],
        "difference",
        hostile_row[8],
    )
    print(
        "regular_hostile_symmetry_defect",
        "cross",
        hostile_cross,
        "variance",
        hostile_variance,
        "covariance",
        hostile_covariance,
        "penalty",
        hostile_penalty,
    )
    for row in lrc_rows:
        print("lrc_projection_row", row)
    print("lrc_same_unmarked_and_pointed", lrc_rows[0][4:] == lrc_rows[1][4:])
    print("semantic_sha256", digest.hexdigest())


if __name__ == "__main__":
    main()
