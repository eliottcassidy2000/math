#!/usr/bin/env python3
"""Independent row-rank audit for THM-4275.

This path never enumerates recurrence words or observer fibres. It builds the
coordinate functionals and row-reduces them over F4/F13.
"""

from collections import Counter
from itertools import combinations
import sys


def require(condition, message="independent audit failure"):
    if not condition:
        raise RuntimeError(message)


def f4_add(left, right):
    return left ^ right


def f4_mul(left, right):
    a0, a1 = left & 1, (left >> 1) & 1
    b0, b1 = right & 1, (right >> 1) & 1
    return (a0 * b0 ^ a1 * b1) | (
        (a0 * b1 ^ a1 * b0 ^ a1 * b1) << 1
    )


def f4_pow(value, exponent):
    answer = 1
    for _ in range(exponent):
        answer = f4_mul(answer, value)
    return answer


def f4_inverse(value):
    require(value != 0, "zero inverse")
    return f4_pow(value, 2)


def rank_f4(matrix):
    rows = [list(row) for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = f4_inverse(rows[rank][column])
        rows[rank] = [f4_mul(inverse, entry) for entry in rows[rank]]
        for index in range(len(rows)):
            if index == rank or rows[index][column] == 0:
                continue
            scalar = rows[index][column]
            rows[index] = [
                f4_add(entry, f4_mul(scalar, pivot_entry))
                for entry, pivot_entry in zip(rows[index], rows[rank])
            ]
        rank += 1
    return rank


OMEGA = 2
ALPHA = f4_mul(OMEGA, OMEGA)
require(f4_mul(OMEGA, ALPHA) == 1, "bad F4")


def ambient_row(index):
    scale = f4_pow(ALPHA, index)
    return tuple(
        f4_mul(scale, entry)
        for entry in (
            1,
            index & 1,
            (index * (index - 1) // 2) & 1,
        )
    )


def visible_row(index):
    scale = f4_pow(OMEGA, index // 2)
    return (scale, 0) if index % 2 == 0 else (0, scale)


ambient_rows = tuple(ambient_row(index) for index in range(12))
visible_rows = tuple(visible_row(index) for index in range(12))

ambient_counts = Counter()
visible_counts = Counter()
ambient_rank_sum = 0
visible_rank_sum = 0
for mask in range(1 << 12):
    sample = [index for index in range(12) if mask & (1 << index)]
    ambient_rank = rank_f4([ambient_rows[index] for index in sample])
    visible_rank = rank_f4([visible_rows[index] for index in sample])
    require(
        ambient_rank == min(len({index % 4 for index in sample}), 3),
        "ambient residue formula",
    )
    require(
        visible_rank == len({index % 2 for index in sample}),
        "visible parity formula",
    )
    ambient_counts[ambient_rank] += 1
    visible_counts[visible_rank] += 1
    ambient_rank_sum += ambient_rank
    visible_rank_sum += visible_rank

ambient_pair_counts = Counter()
ambient_triple_counts = Counter()
visible_pair_counts = Counter()
visible_triple_counts = Counter()
opposite_pairs = 0
same_pairs = 0
gap_counts = Counter()
for sample in combinations(range(12), 2):
    ambient_pair_counts[rank_f4([ambient_rows[index] for index in sample])] += 1
    rank = rank_f4([visible_rows[index] for index in sample])
    visible_pair_counts[rank] += 1
    if (sample[0] - sample[1]) & 1:
        require(rank == 2, "opposite pair rank")
        gap = sample[1] - sample[0]
        gap_counts[min(gap, 12 - gap)] += 1
        opposite_pairs += 1
    else:
        require(rank == 1, "same pair rank")
        same_pairs += 1

for sample in combinations(range(12), 3):
    ambient_triple_counts[rank_f4([ambient_rows[index] for index in sample])] += 1
    visible_triple_counts[rank_f4([visible_rows[index] for index in sample])] += 1


def rank_mod_prime(matrix, prime):
    rows = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    columns = len(rows[0]) if rows else 0
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(inverse * entry) % prime for entry in rows[rank]]
        for index in range(len(rows)):
            if index == rank or rows[index][column] == 0:
                continue
            scalar = rows[index][column]
            rows[index] = [
                (entry - scalar * pivot_entry) % prime
                for entry, pivot_entry in zip(rows[index], rows[rank])
            ]
        rank += 1
    return rank


prime = 13
zeta = 2
contact_rank_counts = Counter()
contact_matrix_cells = 0
for parameter in range(prime):
    matrix = [
        [
            pow(parameter * pow(zeta, row, prime) % prime, column, prime)
            for column in range(12)
        ]
        for row in range(12)
    ]
    rank = rank_mod_prime(matrix, prime)
    require(rank == (1 if parameter == 0 else 12), "contact rank")
    contact_rank_counts[rank] += 1
    contact_matrix_cells += 144

require(ambient_counts == Counter({3: 3773, 2: 294, 1: 28, 0: 1}))
require(visible_counts == Counter({2: 3969, 1: 126, 0: 1}))
require(ambient_pair_counts == Counter({2: 54, 1: 12}))
require(ambient_triple_counts == Counter({3: 108, 2: 108, 1: 4}))
require(visible_pair_counts == Counter({2: 36, 1: 30}))
require(visible_triple_counts == Counter({2: 180, 1: 40}))
require(opposite_pairs == 36 and same_pairs == 30)
require(gap_counts == Counter({1: 12, 3: 12, 5: 12}))


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    print("jc23_confluent_sample_matroids_rank_audit_thm4275")
    print("method=direct_F4_coordinate_row_reduction_no_word_enumeration")
    print(f"ambient_subset_rank_counts={dict(sorted(ambient_counts.items()))}")
    print(f"visible_subset_rank_counts={dict(sorted(visible_counts.items()))}")
    print(f"ambient_subset_rank_sum={ambient_rank_sum}")
    print(f"visible_subset_rank_sum={visible_rank_sum}")
    print(f"ambient_pair_rank_counts={dict(sorted(ambient_pair_counts.items()))}")
    print(f"ambient_triple_rank_counts={dict(sorted(ambient_triple_counts.items()))}")
    print(f"visible_pair_rank_counts={dict(sorted(visible_pair_counts.items()))}")
    print(f"visible_triple_rank_counts={dict(sorted(visible_triple_counts.items()))}")
    print(f"opposite_parity_pairs={opposite_pairs}")
    print(f"opposite_parity_gap_counts={dict(sorted(gap_counts.items()))}")
    print(f"same_parity_pairs={same_pairs}")
    print(f"contact_rank_counts_F13={dict(sorted(contact_rank_counts.items()))}")
    print(f"contact_matrix_cells={contact_matrix_cells}")
    print("status=PASS")


if __name__ == "__main__":
    main()
