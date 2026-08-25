#!/usr/bin/env python3
"""Independent exact audit for THM-4080.

This path does not use DVR pivoting.  Small cases are reconstructed from all
determinantal divisors.  Larger cases audit the proof's weighted-minor lower
bound and independently build a full-rank confluent-matroid witness modulo p.
"""

from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd
import sys


sys.stdout.reconfigure(newline="\n")

EXPECTED_SEMANTIC_SHA256 = "fb55e3c1dd445b524333b4bef3e5cc1825ab52c53033840d1bf9f1424bf385a9"
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def primes_through(limit):
    answer = []
    for candidate in range(2, limit + 1):
        if all(candidate % prime for prime in answer if prime * prime <= candidate):
            answer.append(candidate)
    return tuple(answer)


def vp_integer(value, prime):
    require(value != 0, "valuation needs nonzero integer")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def twojet_matrix(nodes):
    size = 2 * len(nodes)
    return [
        [comb(degree, order) * node ** (degree - order)
         if degree >= order else 0
         for degree in range(size)]
        for node in nodes
        for order in (0, 1)
    ]


def bareiss_det(matrix):
    size = len(matrix)
    require(size >= 1 and all(len(row) == size for row in matrix), "square matrix")
    work = [list(row) for row in matrix]
    if size == 1:
        return work[0][0]
    previous = 1
    sign = 1
    for pivot_index in range(size - 1):
        pivot_row = next((row for row in range(pivot_index, size)
                          if work[row][pivot_index]), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                require(numerator % previous == 0, "Bareiss exact division")
                work[row][column] = numerator // previous
        for row in range(pivot_index + 1, size):
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def submatrix(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def determinantal_divisors(matrix):
    size = len(matrix)
    answer = [1]
    for rank in range(1, size + 1):
        divisor = 0
        for rows in combinations(range(size), rank):
            for columns in combinations(range(size), rank):
                divisor = gcd(divisor, abs(bareiss_det(submatrix(matrix, rows, columns))))
                if divisor == 1:
                    break
            if divisor == 1:
                break
        require(divisor > 0, "full-rank determinantal divisor")
        answer.append(divisor)
    return tuple(answer)


def smith_diagonal(divisors):
    diagonal = []
    for left, right in zip(divisors, divisors[1:]):
        require(right % left == 0, "determinantal divisor chain")
        diagonal.append(right // left)
    return tuple(diagonal)


def single_scale_profile(size, exponent):
    return tuple(sorted(
        (0, 0)
        + tuple(exponent * value for value in range(1, size))
        + tuple(exponent * value for value in range(size + 1, 2 * size))
    ))


def rank_mod_prime(matrix, prime):
    if not matrix:
        return 0
    work = [[entry % prime for entry in row] for row in matrix]
    row_count = len(work)
    column_count = len(work[0])
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(row_count):
            if row == rank or work[row][column] == 0:
                continue
            multiplier = work[row][column]
            work[row] = [
                (left - multiplier * right) % prime
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


def residual_rows(units, column_count):
    rows = []
    for unit in units:
        rows.append(tuple(pow(unit, degree) for degree in range(column_count)))
        rows.append(tuple(
            degree * pow(unit, degree - 1) if degree else 0
            for degree in range(column_count)
        ))
    return tuple(rows)


def predicted_minor_valuation(rank, size, exponent):
    if rank <= size:
        return exponent * (rank - 1) * (rank - 2) // 2
    return exponent * (rank * (rank - 1) // 2 - size)


def main():
    require(sys.argv[1:] == [], "no command-line arguments")

    exhaustive_cases = []
    for prime, exponent, units, translate in (
        (3, 1, (0, 1), 0),
        (3, 3, (1, 2), -4),
        (5, 1, (0, 1, 2), 7),
        (5, 2, (1, 3, 4), -9),
        (5, 1, (0, 1, 2, 3), 0),
        (7, 2, (0, 2, 3, 5), 3),
    ):
        scale = prime ** exponent
        nodes = tuple(translate + scale * unit for unit in units)
        divisors = determinantal_divisors(twojet_matrix(nodes))
        diagonal = smith_diagonal(divisors)
        actual = tuple(0 if entry == 1 else vp_integer(entry, prime)
                       for entry in diagonal)
        expected = single_scale_profile(len(units), exponent)
        require(actual == expected,
                ("all-minor profile", prime, exponent, units, actual, expected))
        exhaustive_cases.append((prime, exponent, units, divisors, actual))

    matroid_cases = []
    witness_minors = 0
    for prime in primes_through(19):
        if prime == 2:
            continue
        for size in range(1, prime):
            for variant, units in enumerate((
                tuple(range(size)),
                tuple((2 * index + 1) % prime for index in range(size)),
            )):
                require(len(set(units)) == size, "distinct residual units")
                rows = residual_rows(units, 2 * size)
                for rank in range(1, 2 * size + 1):
                    columns = tuple(range(rank))
                    if rank <= size:
                        selected = (0,) + tuple(2 * index + 1 for index in range(rank - 1))
                    else:
                        selected_list = [2 * index + 1 for index in range(size)]
                        current_rank = rank_mod_prime(
                            [[rows[row][column] for column in columns]
                             for row in selected_list], prime
                        )
                        require(current_rank == size, "derivative-row independence")
                        for index in range(size):
                            candidate = 2 * index
                            trial = selected_list + [candidate]
                            trial_rank = rank_mod_prime(
                                [[rows[row][column] for column in columns]
                                 for row in trial], prime
                            )
                            if trial_rank > current_rank:
                                selected_list.append(candidate)
                                current_rank = trial_rank
                            if len(selected_list) == rank:
                                break
                        require(len(selected_list) == rank and current_rank == rank,
                                "value-row augmentation")
                        selected = tuple(selected_list)
                    witness = submatrix(rows, selected, columns)
                    require(bareiss_det(witness) % prime != 0,
                            "unit residual witness minor")
                    witness_minors += 1

                    for exponent in (1, 2, 4):
                        lower_candidates = []
                        for derivative_count in range(min(rank, size) + 1):
                            if derivative_count == rank:
                                column_sum = rank * (rank + 1) // 2
                            else:
                                column_sum = rank * (rank - 1) // 2
                            lower_candidates.append(
                                exponent * (column_sum - derivative_count)
                            )
                        lower = min(lower_candidates)
                        require(lower == predicted_minor_valuation(
                            rank, size, exponent
                        ), "weighted-minor lower envelope")
                matroid_cases.append((prime, size, variant))

    consecutive_formula_cases = []
    for prime in primes_through(19):
        if prime == 2:
            continue
        for node_count in range(1, prime * (prime - 1) + 1):
            quotient, remainder = divmod(node_count, prime)
            sizes = tuple(
                quotient + (1 if residue < remainder else 0)
                for residue in range(prime)
                if quotient + (1 if residue < remainder else 0) > 0
            )
            require(max(sizes) < prime, "quadratic cluster bound")
            profile = tuple(sorted(
                exponent
                for size in sizes
                for exponent in single_scale_profile(size, 1)
            ))
            determinant_valuation = 4 * sum(
                vp_integer(distance, prime)
                for distance in range(1, node_count)
                for _ in range(node_count - distance)
            )
            require(len(profile) == 2 * node_count, "global profile length")
            require(sum(profile) == determinant_valuation,
                    "global determinant valuation")
            consecutive_formula_cases.append((prime, node_count, sizes, profile))

    boundary_hostiles = []
    for prime in (2, 3):
        nodes = tuple(prime * index for index in range(prime))
        divisors = determinantal_divisors(twojet_matrix(nodes))
        diagonal = smith_diagonal(divisors)
        actual = tuple(0 if entry == 1 else vp_integer(entry, prime)
                       for entry in diagonal)
        false_extension = single_scale_profile(prime, 1)
        require(actual != false_extension, "sharp size boundary")
        boundary_hostiles.append((prime, actual, false_extension))

    for prime in (3, 5, 7, 11):
        size = prime
        rows = residual_rows(tuple(range(size)), 2 * size - 1)
        derivative_rows = [rows[2 * index + 1] for index in range(size)]
        require(rank_mod_prime(derivative_rows, prime) == size - 1,
                "boundary derivative rank drop")

    semantic = {
        "exhaustive": [
            [prime, exponent, list(units), list(divisors), list(profile)]
            for prime, exponent, units, divisors, profile in exhaustive_cases
        ],
        "matroid_cases": [list(row) for row in matroid_cases],
        "witness_minors": witness_minors,
        "consecutive": [
            [prime, count, list(sizes), list(profile)]
            for prime, count, sizes, profile in consecutive_formula_cases
        ],
        "hostiles": [
            [prime, list(actual), list(false_extension)]
            for prime, actual, false_extension in boundary_hostiles
        ],
    }
    semantic_digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("CONFLUENT_TWOJET_SINGLE_SCALE_SMITH_THM4080_INDEPENDENT_AUDIT")
    print("status=ALL_MINORS_PLUS_MOD_P_MATROID;PASS")
    print(f"all_determinantal_divisor_cases={len(exhaustive_cases)}")
    print(f"matroid_cases={len(matroid_cases)};unit_witness_minors={witness_minors}")
    print(f"quadratic_range_formula_cases={len(consecutive_formula_cases)}")
    for prime, actual, false_extension in boundary_hostiles:
        print(f"hostile_s_eq_p(p={prime})=actual:{actual};false_extension:{false_extension}")
    print("boundary_mechanism=derivative_row_rank_drops_from_p_to_p_minus_1")
    print(f"semantic_sha256={semantic_digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
