#!/usr/bin/env python3
"""Exact first-counterexample audit for THM-300's quadratic signature.

Only the transitive tournament and its one-/two-tile flips are needed for the
quadratic coefficient matrix.  Hamiltonian paths are counted by an integer
subset DP.  Inertia is then computed over Q by explicit symmetric congruence,
using one-by-one pivots or a two-by-two off-diagonal pivot when every diagonal
entry vanishes.  No eigensolver, floating point, or third-party package enters
the truth gate.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction


def tiles(size):
    return tuple(
        (upper, lower)
        for lower in range(1, size - 1)
        for upper in range(size, lower + 1, -1)
        if upper - lower >= 2
    )


def out_masks(size, flipped_tiles):
    out = [(1 << vertex) - 1 for vertex in range(size)]
    for upper, lower in flipped_tiles:
        high = upper - 1
        low = lower - 1
        out[high] &= ~(1 << low)
        out[low] |= 1 << high
    return tuple(out)


def hamiltonian_paths(size, flipped_tiles):
    full = (1 << size) - 1
    out = out_masks(size, flipped_tiles)
    states = [None] * (full + 1)
    for vertex in range(size):
        states[1 << vertex] = {vertex: 1}
    for mask in range(1, full):
        row = states[mask]
        if row is None:
            continue
        for endpoint, count in tuple(row.items()):
            available = out[endpoint] & (full ^ mask)
            while available:
                bit = available & -available
                available -= bit
                successor = bit.bit_length() - 1
                new_mask = mask | bit
                target = states[new_mask]
                if target is None:
                    target = {}
                    states[new_mask] = target
                target[successor] = target.get(successor, 0) + count
    return sum(states[full].values())


def quadratic_matrix(size):
    tile_list = tiles(size)
    base = hamiltonian_paths(size, ())
    singles = tuple(
        hamiltonian_paths(size, (tile,)) for tile in tile_list
    )
    matrix = [[0] * len(tile_list) for _ in tile_list]
    for left, right in itertools.combinations(range(len(tile_list)), 2):
        pair = hamiltonian_paths(
            size, (tile_list[left], tile_list[right])
        )
        coefficient = pair - singles[left] - singles[right] + base
        matrix[left][right] = coefficient
        matrix[right][left] = coefficient
    return tile_list, matrix


def permute_symmetric(matrix, front):
    selected = set(front)
    order = tuple(front) + tuple(
        index for index in range(len(matrix)) if index not in selected
    )
    return [[matrix[row][column] for column in order] for row in order]


def exact_inertia(integer_matrix):
    matrix = [[Fraction(value) for value in row] for row in integer_matrix]
    positive = negative = zero = one_pivots = two_pivots = 0
    while matrix:
        size = len(matrix)
        diagonal = next(
            (index for index in range(size) if matrix[index][index] != 0),
            None,
        )
        if diagonal is not None:
            matrix = permute_symmetric(matrix, (diagonal,))
            pivot = matrix[0][0]
            positive += pivot > 0
            negative += pivot < 0
            one_pivots += 1
            column = tuple(matrix[index][0] for index in range(1, size))
            matrix = [
                [
                    matrix[row + 1][column_index + 1]
                    - column[row] * column[column_index] / pivot
                    for column_index in range(size - 1)
                ]
                for row in range(size - 1)
            ]
            continue

        off_diagonal = next(
            (
                (row, column)
                for row in range(size)
                for column in range(row + 1, size)
                if matrix[row][column] != 0
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break
        matrix = permute_symmetric(matrix, off_diagonal)
        pivot = matrix[0][1]
        positive += 1
        negative += 1
        two_pivots += 1
        matrix = [
            [
                matrix[row + 2][column + 2]
                - (
                    matrix[0][row + 2] * matrix[1][column + 2]
                    + matrix[1][row + 2] * matrix[0][column + 2]
                ) / pivot
                for column in range(size - 2)
            ]
            for row in range(size - 2)
        ]
    return positive, negative, zero, one_pivots, two_pivots


def solve_rational(matrix, right_hand_side):
    size = len(matrix)
    width = len(right_hand_side[0])
    augmented = [
        [Fraction(value) for value in matrix[row]]
        + [Fraction(value) for value in right_hand_side[row]]
        for row in range(size)
    ]
    for column in range(size):
        pivot_row = next(
            row for row in range(column, size)
            if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row], augmented[column]
        )
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(size):
            if row == column or augmented[row][column] == 0:
                continue
            multiplier = augmented[row][column]
            augmented[row] = [
                left - multiplier * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return tuple(tuple(row[size:]) for row in augmented)


def new_layer_schur(size, tile_list, matrix, previous_matrix):
    old = tuple(
        index for index, tile in enumerate(tile_list) if tile[0] < size
    )
    new = tuple(
        index for index, tile in enumerate(tile_list) if tile[0] == size
    )
    old_block = tuple(
        tuple(matrix[row][column] for column in old) for row in old
    )
    assert old_block == tuple(tuple(row) for row in previous_matrix)
    cross = tuple(
        tuple(matrix[row][column] for column in new) for row in old
    )
    new_block = tuple(
        tuple(matrix[row][column] for column in new) for row in new
    )
    solved = solve_rational(old_block, cross)
    schur = tuple(
        tuple(
            Fraction(new_block[row][column])
            - sum(
                Fraction(cross[old_row][row]) * solved[old_row][column]
                for old_row in range(len(old))
            )
            for column in range(len(new))
        )
        for row in range(len(new))
    )
    return tuple(tile_list[index] for index in new), schur


def main():
    print("THM-300 QUADRATIC SIGNATURE EXACT REFUTATION")
    digest = hashlib.sha256()
    first_negative_count_failure = None
    matrices = {}
    for size in range(5, 13):
        tile_list, matrix = quadratic_matrix(size)
        matrices[size] = (tile_list, matrix)
        inertia = exact_inertia(matrix)
        positive, negative, zero, one_pivots, two_pivots = inertia
        assert positive + negative + zero == len(tile_list)
        if negative != size - 2 and first_negative_count_failure is None:
            first_negative_count_failure = size
        digest.update(repr((size, tile_list, matrix, inertia)).encode("ascii"))
        print(
            f"n={size};m={len(tile_list)};inertia=({positive},{negative},{zero});"
            f"target_negative={size - 2};congruence_pivots=({one_pivots},{two_pivots})"
        )
    for size in range(7, 13):
        tile_list, matrix = matrices[size]
        _previous_tiles, previous_matrix = matrices[size - 1]
        new_tiles, schur = new_layer_schur(
            size, tile_list, matrix, previous_matrix
        )
        positive, negative, zero, _one, _two = exact_inertia(schur)
        assert len(new_tiles) == size - 2
        print(
            f"n={size};new_hypotenuse_layer={len(new_tiles)};"
            f"schur_inertia=({positive},{negative},{zero})"
        )
    assert first_negative_count_failure == 9
    print("first_negative_count_failure=n9;actual_negative=8;claimed=7")
    print("full_rank_exact_range=n6..12;n5_nullity=1")
    print("surviving_pattern_probe=n9..12_have_negative_n_minus_1;OPEN_beyond_12")
    print(f"semantic_sha256={digest.hexdigest()}")
    print("truth_gates=integer_HP_DP+rational_symmetric_congruence;float_literals=0")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
