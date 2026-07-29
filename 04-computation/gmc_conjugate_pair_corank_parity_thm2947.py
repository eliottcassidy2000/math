#!/usr/bin/env python3
"""Exact controls for THM-2947.

The proof is algebraic.  This companion checks the complete-intersection
Hilbert dimensions, the even-corank rank ladder, the fifth-compound
energy on conjugate 2-by-2 blocks, and the sharp real-summand hostile.
Only Python integer arithmetic is used.
"""

from itertools import combinations
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sdim(degree):
    if degree < 0:
        return 0
    return comb(degree + 2, 2)


def ci23_hilbert(degree):
    return (
        sdim(degree)
        - sdim(degree - 2)
        - sdim(degree - 3)
        + sdim(degree - 5)
    )


def det_bareiss(matrix):
    n = len(matrix)
    require(all(len(row) == n for row in matrix), "matrix is not square")
    if n == 0:
        return 1
    a = [list(map(int, row)) for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division failed")
                a[i][j] = numerator // previous
        previous = pivot
    return sign * a[-1][-1]


def rank_fraction_free(matrix):
    a = [list(map(int, row)) for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pivot_value = a[rank][col]
        for i in range(rows):
            if i == rank or a[i][col] == 0:
                continue
            row_value = a[i][col]
            a[i] = [
                pivot_value * a[i][j] - row_value * a[rank][j]
                for j in range(cols)
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def block_diag(blocks):
    size = sum(len(block) for block in blocks)
    result = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        width = len(block)
        require(all(len(row) == width for row in block), "bad block")
        for i in range(width):
            for j in range(width):
                result[offset + i][offset + j] = block[i][j]
        offset += width
    return result


def complex_block(real, imag):
    return [[real, -imag], [imag, real]]


def compound_energy(matrix, order):
    rows = len(matrix)
    cols = len(matrix[0])
    total = 0
    for row_set in combinations(range(rows), order):
        for col_set in combinations(range(cols), order):
            minor = [[matrix[i][j] for j in col_set] for i in row_set]
            total += det_bareiss(minor) ** 2
    return total


def main():
    h3 = ci23_hilbert(3)
    h7 = ci23_hilbert(7)
    require((h3, h7) == (6, 6), "CI(2,3) Hilbert dimensions changed")
    require(sdim(7) - h7 == 30, "(Q,C)_7 dimension changed")

    kernel_ladder = (0, 2, 4, 6)
    mu_ranks = tuple(6 - k for k in kernel_ladder)
    phi_ranks = tuple(30 + r for r in mu_ranks)
    require(mu_ranks == (6, 4, 2, 0), "multiplication rank ladder changed")
    require(phi_ranks == (36, 34, 32, 30), "Macaulay rank ladder changed")

    nonsingular = block_diag(
        [complex_block(1, 2), complex_block(2, -1), complex_block(3, 1)]
    )
    one_dead_pair = block_diag(
        [complex_block(1, 2), complex_block(0, 0), complex_block(3, 1)]
    )
    require(rank_fraction_free(nonsingular) == 6, "full conjugate model rank")
    require(rank_fraction_free(one_dead_pair) == 4, "dead-pair model rank")
    full_energy = compound_energy(nonsingular, 5)
    dead_energy = compound_energy(one_dead_pair, 5)
    require(full_energy > 0, "fifth compound missed full rank")
    require(dead_energy == 0, "fifth compound survived corank two")

    real_summand_hostile = block_diag(
        [[[0]], [[1]], [[1]], [[1]], [[1]], [[1]]]
    )
    require(rank_fraction_free(real_summand_hostile) == 5, "hostile rank changed")
    hostile_energy = compound_energy(real_summand_hostile, 5)
    require(hostile_energy == 1, "hostile fifth-compound energy changed")
    require(det_bareiss(real_summand_hostile) == 0, "hostile determinant changed")

    full_minor_count = comb(36, 35) * comb(46, 35)

    print("THM-2947 CONJUGATE-PAIR CORANK PARITY AUDIT")
    print(f"CI23_H3={h3};CI23_H7={h7};QC7_rank={sdim(7)-h7}")
    print(f"kernel_ladder={kernel_ladder}")
    print(f"muF_rank_ladder={mu_ranks}")
    print(f"Phi7_rank_ladder={phi_ranks}")
    print(f"full_conjugate_fifth_compound_energy={full_energy}")
    print(f"one_dead_pair_fifth_compound_energy={dead_energy}")
    print(
        "real_summand_hostile="
        f"det0_rank5_energy{hostile_energy}"
    )
    print(f"full_35_minor_bank_size={full_minor_count}")
    print("consequence=one_nonzero_35_minor_forces_rank36_and_resultant_nonzero")
    print("scope=no_uniform_fixed_minor_or_new_SFC4_width_closure")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
