#!/usr/bin/env python3
"""Exact bounded controls for THM-3729's switching-class mean envelope."""

from __future__ import annotations

from hashlib import sha256
from math import factorial


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    work = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (row for row in range(pivot_index + 1, size)
                 if work[row][pivot_index]),
                None,
            )
            require(swap is not None, ("singular pivot", pivot_index))
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(numerator % previous == 0,
                        ("Bareiss division", pivot_index, numerator, previous))
                work[row][column] = numerator // previous
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def representative(order, tail_mask):
    matrix = [[0] * order for _ in range(order)]
    for right in range(1, order):
        matrix[0][right] = 1
        matrix[right][0] = -1
    bit = 0
    for left in range(1, order):
        for right in range(left + 1, order):
            sign = 1 if (tail_mask >> bit) & 1 else -1
            matrix[left][right] = sign
            matrix[right][left] = -sign
            bit += 1
    return matrix


def switched(matrix, switch_mask):
    order = len(matrix)
    signs = [1] + [
        -1 if (switch_mask >> (vertex - 1)) & 1 else 1
        for vertex in range(1, order)
    ]
    return [
        [signs[row] * matrix[row][column] * signs[column]
         for column in range(order)]
        for row in range(order)
    ]


def identity_plus(matrix):
    return [
        [matrix[row][column] + (row == column)
         for column in range(len(matrix))]
        for row in range(len(matrix))
    ]


def delete_vertex(matrix, vertex):
    return [
        [entry for column, entry in enumerate(row) if column != vertex]
        for row_index, row in enumerate(matrix) if row_index != vertex
    ]


def rooted_extension(matrix):
    order = len(matrix)
    return [row + [1] for row in matrix] + [[-1] * order + [0]]


def hamiltonian_paths(matrix):
    order = len(matrix)
    total_masks = 1 << order
    dynamic = [[0] * order for _ in range(total_masks)]
    for vertex in range(order):
        dynamic[1 << vertex][vertex] = 1
    for mask in range(1, total_masks):
        for last in range(order):
            if not (mask >> last) & 1:
                continue
            previous_mask = mask ^ (1 << last)
            if not previous_mask:
                continue
            dynamic[mask][last] = sum(
                dynamic[previous_mask][previous]
                for previous in range(order)
                if (previous_mask >> previous) & 1
                and matrix[previous][last] == 1
            )
    return sum(dynamic[-1])


def stable_digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def main():
    rows = []
    semantic = []
    total_classes = 0
    total_members = 0
    for order in range(2, 7):
        class_size = 1 << (order - 1)
        class_count = 1 << ((order - 1) * (order - 2) // 2)
        target_h_sum = factorial(order)
        even_slacks = []
        odd_slacks = []
        maximum_disc = 0
        maximum_deletion_sum = 0

        for tail_mask in range(class_count):
            base = representative(order, tail_mask)
            base_det = determinant(identity_plus(base))
            require(base_det % class_size == 0,
                    ("disc normalization", order, tail_mask, base_det))
            disc = base_det // class_size
            deletion_sum = 0
            for vertex in range(order):
                deleted_det = determinant(identity_plus(
                    delete_vertex(base, vertex)
                ))
                deletion_denominator = 1 << (order - 2)
                require(deleted_det % deletion_denominator == 0,
                        ("deletion normalization", order, tail_mask, vertex))
                deletion_sum += deleted_det // deletion_denominator

            h_sum = 0
            odd_sum = 0
            seen = set()
            for switch_mask in range(class_size):
                member = switched(base, switch_mask)
                key = tuple(tuple(row) for row in member)
                require(key not in seen,
                        ("nonfree switching", order, tail_mask, switch_mask))
                seen.add(key)
                member_det = determinant(identity_plus(member))
                require(member_det == base_det,
                        ("disc not switching invariant", order, tail_mask))
                h_sum += hamiltonian_paths(member)
                augmented_det = determinant(identity_plus(
                    rooted_extension(member)
                ))
                numerator = augmented_det - member_det
                require(numerator % class_size == 0,
                        ("odd normalization", order, tail_mask, switch_mask))
                odd_sum += numerator // class_size

            require(len(seen) == class_size,
                    ("class size", order, tail_mask, len(seen)))
            require(h_sum == target_h_sum,
                    ("Hamiltonian class mean", order, tail_mask, h_sum))
            require(odd_sum == (class_size // 2) * deletion_sum,
                    ("rooted deletion mean", order, tail_mask,
                     odd_sum, deletion_sum))

            even_slack = target_h_sum - class_size * disc
            odd_slack = target_h_sum - odd_sum
            require(even_slack >= 0,
                    ("Hadamard even envelope", order, tail_mask, even_slack))
            require(odd_slack >= 0,
                    ("Hadamard odd envelope", order, tail_mask, odd_slack))
            even_slacks.append(even_slack)
            odd_slacks.append(odd_slack)
            maximum_disc = max(maximum_disc, disc)
            maximum_deletion_sum = max(maximum_deletion_sum, deletion_sum)

        row = (
            order,
            class_count,
            class_count * class_size,
            min(even_slacks),
            sum(slack == 0 for slack in even_slacks),
            min(odd_slacks),
            sum(slack == 0 for slack in odd_slacks),
            maximum_disc,
            maximum_deletion_sum,
        )
        rows.append(row)
        semantic.append((row, stable_digest(tuple(even_slacks)),
                         stable_digest(tuple(odd_slacks))))
        total_classes += class_count
        total_members += class_count * class_size

    print("tournament switching-class mean envelope")
    print("orders=2..6")
    print("row_fields=(order,classes,members,min_even_slack,even_equalities,"
          "min_odd_slack,odd_equalities,max_disc,max_deletion_sum)")
    print("rows=" + repr(tuple(rows)))
    print("total_classes=" + repr(total_classes))
    print("total_members=" + repr(total_members))
    print("semantic_sha256=" + stable_digest(tuple(semantic)))
    print("scope=finite_controls_for_all_scale_switching_and_hadamard_proof")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
