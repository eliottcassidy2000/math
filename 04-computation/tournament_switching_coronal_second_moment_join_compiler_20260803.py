#!/usr/bin/env python3
"""Exact switching-cube second moments and order-join coronal compiler.

For a tournament centered adjacency C=2A-J, put

    P(z)=det(I-zC), B(z)=adj(I-zC), N_d(z)=d^T B(z) d.

The Rademacher second moment of N_d is checked coefficientwise against an
exact formula consisting of two P-determined trace kernels and the Gram
kernel of the diagonal cofactors B_ii(z).  The same script finds and freezes
the first P-collision at which that principal-minor sidecar differs.

It also checks the exact order-join multiplication law for the pairs (P,N)
attached to selected switched tournaments.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
from collections import Counter
from pathlib import Path


EXPECTED_LABELLED_MOMENT_COUNTS = (
    (1, 1, 1),
    (2, 2, 4),
    (3, 8, 32),
    (4, 64, 512),
    (5, 1024, 16_384),
)
EXPECTED_P_PROFILE_COUNTS = (
    (1, 1, 1),
    (2, 2, 1),
    (3, 8, 1),
    (4, 64, 2),
    (5, 1024, 2),
    (6, 32_768, 6),
)
EXPECTED_JOIN_CASES = 1369
EXPECTED_TRIANGLE = (
    (1, 3, 6, 4),
    ((1, 2, 2), (1, 2, 2), (1, 2, 2)),
    (((3, 6, 4), 3), ((3, 6, 12), 1)),
)
EXPECTED_COLLISION = (
    4164,
    5122,
    (1, 7, 42, 140, 360, 576, 576, 256),
    (
        (1, 6, 30, 80, 152, 160, 64),
        (1, 6, 30, 80, 152, 160, 64),
        (1, 6, 30, 80, 144, 144, 64),
        (1, 6, 30, 80, 144, 144, 64),
        (1, 6, 30, 80, 152, 160, 64),
        (1, 6, 30, 80, 168, 192, 128),
        (1, 6, 30, 80, 168, 192, 128),
    ),
    (
        (1, 6, 30, 80, 144, 144, 64),
        (1, 6, 30, 80, 152, 160, 64),
        (1, 6, 30, 80, 144, 144, 64),
        (1, 6, 30, 80, 160, 176, 96),
        (1, 6, 30, 80, 152, 160, 64),
        (1, 6, 30, 80, 168, 192, 128),
        (1, 6, 30, 80, 160, 176, 96),
    ),
    (6, 5, 3, 2, 2, 2, 1),
    (5, 5, 3, 3, 2, 2, 1),
    (0, 0, 0, 0, 1, 2, 4),
)
EXPECTED_JOIN_HOSTILE = (
    (1, 4, 12, 16, 16),
    (4, 12, 24, 16),
    (3, 1, 1, 1),
    (2, 2, 2, 0),
)
EXPECTED_MOMENT_SHA256 = "32541a343cc90cc1a0573f3e9851b7b7700d9e8e05f911a4471b1506262e5726"
EXPECTED_COLLISION_SHA256 = "a8b584345412058b846292e376e52e922ddfce1930f4990d49aaea79963d1dad"
EXPECTED_JOIN_SHA256 = "fab670744525cfeecf0bd2f0b43a09a86b60f1d90b8e8c065dd8ae747026ff63"
EXPECTED_SEMANTIC_SHA256 = "40427da1f4e5a2d18a971f25c3a81ccaac48314b81d16b2e1cc4f409142ba345"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(order):
    return tuple(
        tuple(int(row == column) for column in range(order))
        for row in range(order)
    )


def zero_matrix(order):
    return tuple(tuple(0 for _column in range(order)) for _row in range(order))


def matrix_product(left, right):
    rows = len(left)
    middle_size = len(right)
    columns = len(right[0]) if middle_size else 0
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(middle_size))
            for column in range(columns)
        )
        for row in range(rows)
    )


def matrix_add_scalar_identity(matrix, scalar):
    return tuple(
        tuple(
            value + (scalar if row == column else 0)
            for column, value in enumerate(matrix[row])
        )
        for row in range(len(matrix))
    )


def matrix_scale(matrix, scalar):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def matrix_transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix)))


def trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def trace_product(left, right):
    order = len(left)
    return sum(
        left[row][column] * right[column][row]
        for row in range(order)
        for column in range(order)
    )


def frobenius_product(left, right):
    order = len(left)
    return sum(
        left[row][column] * right[row][column]
        for row in range(order)
        for column in range(order)
    )


def polynomial_add(left, right):
    size = max(len(left), len(right))
    return tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def polynomial_scale(polynomial, scalar):
    return tuple(scalar * value for value in polynomial)


def polynomial_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return tuple(result)


def polynomial_shift(polynomial, amount):
    return (0,) * amount + tuple(polynomial)


def polynomial_pad(polynomial, size):
    require(len(polynomial) <= size, (polynomial, size))
    return tuple(polynomial) + (0,) * (size - len(polynomial))


def outer(left, right):
    return tuple(tuple(a * b for b in right) for a in left)


def bivariate_add(*matrices):
    rows = max(len(matrix) for matrix in matrices)
    columns = max(len(matrix[0]) if matrix else 0 for matrix in matrices)
    return tuple(
        tuple(
            sum(
                matrix[row][column]
                if row < len(matrix) and column < len(matrix[row]) else 0
                for matrix in matrices
            )
            for column in range(columns)
        )
        for row in range(rows)
    )


def bivariate_scale(matrix, scalar):
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def bivariate_shift(matrix, row_shift, column_shift):
    rows = len(matrix) + row_shift
    columns = (len(matrix[0]) if matrix else 0) + column_shift
    result = [[0] * columns for _row in range(rows)]
    for row in range(len(matrix)):
        for column in range(len(matrix[row])):
            result[row + row_shift][column + column_shift] = matrix[row][column]
    return tuple(tuple(row) for row in result)


def bivariate_factor(matrix, terms):
    max_row_shift = max(term[0] for term in terms)
    max_column_shift = max(term[1] for term in terms)
    result = [
        [0] * (len(matrix[0]) + max_column_shift)
        for _row in range(len(matrix) + max_row_shift)
    ]
    for row_shift, column_shift, scalar in terms:
        for row in range(len(matrix)):
            for column in range(len(matrix[row])):
                result[row + row_shift][column + column_shift] += (
                    scalar * matrix[row][column]
                )
    return tuple(tuple(row) for row in result)


def char_polynomial_i_minus(matrix):
    """Return det(I-zM) from Newton traces, with exact division checks."""
    order = len(matrix)
    if order == 0:
        return (1,)
    power = identity(order)
    traces = []
    for _degree in range(1, order + 1):
        power = matrix_product(power, matrix)
        traces.append(trace(power))
    coefficients = [1]
    for degree in range(1, order + 1):
        numerator = -sum(
            traces[index - 1] * coefficients[degree - index]
            for index in range(1, degree + 1)
        )
        require(numerator % degree == 0, (matrix, degree, numerator))
        coefficients.append(numerator // degree)
    return tuple(coefficients)


def determinant_polynomial_permutation(matrix):
    """Independent direct determinant of I-zM by permutation expansion."""
    order = len(matrix)
    result = (0,) * (order + 1)
    for permutation in itertools.permutations(range(order)):
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(order)
            for right in range(left + 1, order)
        )
        sign = -1 if inversions % 2 else 1
        term = (1,)
        for row, column in enumerate(permutation):
            factor = (
                (1, -matrix[row][column])
                if row == column else (0, -matrix[row][column])
            )
            term = polynomial_multiply(term, factor)
        result = polynomial_add(result, polynomial_scale(term, sign))
    return result


def adjugate_coefficients(matrix, denominator):
    """Coefficient matrices of adj(I-zM), using (I-zM)B=P I."""
    order = len(matrix)
    require(len(denominator) == order + 1, denominator)
    coefficients = [identity(order)]
    for degree in range(1, order):
        coefficients.append(
            matrix_add_scalar_identity(
                matrix_product(matrix, coefficients[-1]),
                denominator[degree],
            )
        )
    terminal = matrix_product(matrix, coefficients[-1])
    expected_terminal = matrix_scale(identity(order), -denominator[-1])
    require(terminal == expected_terminal, (matrix, denominator, terminal))
    return tuple(coefficients)


def adjacency(order, mask):
    matrix = [[0] * order for _row in range(order)]
    edge = 0
    for left in range(order):
        for right in range(left + 1, order):
            orientation = (mask >> edge) & 1
            matrix[left][right] = orientation
            matrix[right][left] = 1 - orientation
            edge += 1
    return tuple(tuple(row) for row in matrix)


def centered(matrix):
    order = len(matrix)
    return tuple(
        tuple(2 * matrix[row][column] - 1 for column in range(order))
        for row in range(order)
    )


def validate_tournament(matrix):
    order = len(matrix)
    return all(
        matrix[row][row] == 0
        and all(
            matrix[row][column] in (0, 1)
            and matrix[row][column] + matrix[column][row] == 1
            for column in range(order)
            if row != column
        )
        for row in range(order)
    )


def delete_vertex(matrix, vertex):
    order = len(matrix)
    return tuple(
        tuple(matrix[row][column] for column in range(order) if column != vertex)
        for row in range(order)
        if row != vertex
    )


def switch(matrix, signs):
    order = len(matrix)
    return tuple(
        tuple(
            0 if row == column else (
                matrix[row][column]
                if signs[row] == signs[column]
                else 1 - matrix[row][column]
            )
            for column in range(order)
        )
        for row in range(order)
    )


def relabel(matrix, permutation):
    return tuple(
        tuple(matrix[permutation[row]][permutation[column]]
              for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def order_join(left, right):
    left_order = len(left)
    right_order = len(right)
    order = left_order + right_order
    return tuple(
        tuple(
            left[row][column]
            if row < left_order and column < left_order
            else right[row - left_order][column - left_order]
            if row >= left_order and column >= left_order
            else 1
            if row < left_order and column >= left_order
            else 0
            for column in range(order)
        )
        for row in range(order)
    )


def quadratic_form(matrix, vector):
    order = len(matrix)
    return sum(
        vector[row] * matrix[row][column] * vector[column]
        for row in range(order)
        for column in range(order)
    )


def observer_numerator(adjugate, signs):
    return tuple(quadratic_form(coefficient, signs) for coefficient in adjugate)


def principal_deck(adjugate):
    order = len(adjugate[0])
    return tuple(
        tuple(adjugate[degree][vertex][vertex] for degree in range(order))
        for vertex in range(order)
    )


def deck_gram(deck):
    degree_count = len(deck[0])
    return tuple(
        tuple(
            sum(polynomial[left] * polynomial[right] for polynomial in deck)
            for right in range(degree_count)
        )
        for left in range(degree_count)
    )


def spectral_kernels(adjugate):
    degree_count = len(adjugate)
    same = tuple(
        tuple(trace_product(adjugate[left], adjugate[right])
              for right in range(degree_count))
        for left in range(degree_count)
    )
    transpose = tuple(
        tuple(frobenius_product(adjugate[left], adjugate[right])
              for right in range(degree_count))
        for left in range(degree_count)
    )
    return same, transpose


def second_moment_packet(matrix, direct_determinant):
    order = len(matrix)
    core = centered(matrix)
    require(
        tuple(
            tuple(core[row][column] + core[column][row]
                  for column in range(order))
            for row in range(order)
        ) == matrix_scale(identity(order), -2),
        (order, "centered transpose relation"),
    )
    denominator = char_polynomial_i_minus(core)
    if direct_determinant:
        require(
            denominator == determinant_polynomial_permutation(core),
            (matrix, denominator, "direct determinant"),
        )
    adjugate = adjugate_coefficients(core, denominator)
    mean = tuple(trace(coefficient) for coefficient in adjugate)
    expected_mean = tuple(
        (order - degree) * denominator[degree]
        for degree in range(order)
    )
    require(mean == expected_mean, (matrix, mean, expected_mean))

    deck = principal_deck(adjugate)
    direct_deck = tuple(
        char_polynomial_i_minus(delete_vertex(core, vertex))
        for vertex in range(order)
    )
    require(deck == direct_deck, (matrix, deck, direct_deck))
    gram = deck_gram(deck)
    same_kernel, transpose_kernel = spectral_kernels(adjugate)
    covariance = bivariate_add(
        same_kernel,
        transpose_kernel,
        bivariate_scale(gram, -2),
    )
    raw_second = bivariate_add(outer(mean, mean), covariance)

    same_left = bivariate_factor(same_kernel, ((1, 0, 1), (0, 1, -1)))
    same_right = bivariate_add(
        bivariate_shift(outer(mean, denominator), 1, 0),
        bivariate_scale(bivariate_shift(outer(denominator, mean), 0, 1), -1),
    )
    require(same_left == same_right, (matrix, "same resolvent kernel"))

    transpose_left = bivariate_factor(
        transpose_kernel,
        ((1, 0, 1), (0, 1, 1), (1, 1, 2)),
    )
    transpose_right = bivariate_add(
        bivariate_shift(outer(mean, denominator), 1, 0),
        bivariate_shift(outer(denominator, mean), 0, 1),
    )
    require(
        transpose_left == transpose_right,
        (matrix, "transpose resolvent kernel"),
    )

    cut_count = 1 << (order - 1)
    numerator_sum = [0] * order
    second_sum = [[0] * order for _row in range(order)]
    determinant_sum = [0] * (order + 1)
    determinant_second_sum = [[0] * (order + 1) for _row in range(order + 1)]
    cut_records = []
    for tail in itertools.product((1, -1), repeat=order - 1):
        signs = (1, *tail)
        numerator = observer_numerator(adjugate, signs)
        determinant = tuple(
            denominator[degree]
            - (numerator[degree - 1] if degree else 0)
            for degree in range(order + 1)
        )
        direct_switched_determinant = char_polynomial_i_minus(
            matrix_scale(switch(matrix, signs), 2)
        )
        require(
            determinant == direct_switched_determinant,
            (matrix, signs, determinant, direct_switched_determinant),
        )
        for left in range(order):
            numerator_sum[left] += numerator[left]
            for right in range(order):
                second_sum[left][right] += numerator[left] * numerator[right]
        for left in range(order + 1):
            determinant_sum[left] += determinant[left]
            for right in range(order + 1):
                determinant_second_sum[left][right] += (
                    determinant[left] * determinant[right]
                )
        cut_records.append((signs, numerator, determinant))

    require(
        tuple(numerator_sum) == tuple(cut_count * value for value in mean),
        (matrix, "numerator mean"),
    )
    require(
        tuple(tuple(row) for row in second_sum)
        == bivariate_scale(raw_second, cut_count),
        (matrix, "numerator second moment"),
    )

    determinant_mean = tuple(
        denominator[degree] - (mean[degree - 1] if degree else 0)
        for degree in range(order + 1)
    )
    determinant_covariance = bivariate_shift(covariance, 1, 1)
    determinant_raw_second = bivariate_add(
        outer(determinant_mean, determinant_mean),
        determinant_covariance,
    )
    require(
        tuple(determinant_sum)
        == tuple(cut_count * value for value in determinant_mean),
        (matrix, "determinant mean"),
    )
    require(
        tuple(tuple(row) for row in determinant_second_sum)
        == bivariate_scale(determinant_raw_second, cut_count),
        (matrix, "determinant second moment"),
    )
    return {
        "core": core,
        "P": denominator,
        "B": adjugate,
        "M": mean,
        "deck": deck,
        "gram": gram,
        "Ksame": same_kernel,
        "Ktranspose": transpose_kernel,
        "covariance": covariance,
        "raw_second": raw_second,
        "determinant_mean": determinant_mean,
        "determinant_raw_second": determinant_raw_second,
        "cuts": tuple(cut_records),
    }


def relabel_control(matrix, packet):
    order = len(matrix)
    permutation = tuple(reversed(range(order)))
    relabelled = relabel(matrix, permutation)
    relabelled_packet = second_moment_packet(relabelled, False)
    require(relabelled_packet["P"] == packet["P"], (matrix, "relabel P"))
    require(relabelled_packet["gram"] == packet["gram"], (matrix, "relabel Gram"))
    require(
        relabelled_packet["deck"]
        == tuple(packet["deck"][permutation[index]] for index in range(order)),
        (matrix, "relabel deck"),
    )
    for signs, numerator, _determinant in packet["cuts"]:
        relabelled_signs = tuple(signs[permutation[index]] for index in range(order))
        require(
            observer_numerator(relabelled_packet["B"], relabelled_signs)
            == numerator,
            (matrix, signs, "relabel observer"),
        )


def triangle_control():
    transitive = (
        (0, 1, 1),
        (0, 0, 1),
        (0, 0, 0),
    )
    cycle = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    packets = tuple(second_moment_packet(matrix, True) for matrix in (transitive, cycle))
    summaries = []
    for packet in packets:
        census = tuple(sorted(Counter(record[1] for record in packet["cuts"]).items()))
        summaries.append((packet["P"], packet["deck"], census))
    require(summaries[0] == summaries[1] == EXPECTED_TRIANGLE, summaries)
    return summaries[0]


def score_sequence(matrix):
    return tuple(sorted((sum(row) for row in matrix), reverse=True))


def minimal_sidecar_collision():
    digest = hashlib.sha256()
    profile_counts = []
    for order in range(1, 7):
        tournament_count = 1 << (order * (order - 1) // 2)
        by_denominator = {}
        for mask in range(tournament_count):
            matrix = adjacency(order, mask)
            core = centered(matrix)
            denominator = char_polynomial_i_minus(core)
            deck = tuple(
                char_polynomial_i_minus(delete_vertex(core, vertex))
                for vertex in range(order)
            )
            gram = deck_gram(deck)
            if denominator in by_denominator:
                require(
                    by_denominator[denominator] == gram,
                    (order, mask, denominator, "premature P/Gram collision"),
                )
            else:
                by_denominator[denominator] = gram
            digest.update(f"{order}|{mask}|{denominator}|{gram}\n".encode())
        profile_counts.append((order, tournament_count, len(by_denominator)))
    profile_counts = tuple(profile_counts)
    require(profile_counts == EXPECTED_P_PROFILE_COUNTS, profile_counts)

    first_mask, second_mask = EXPECTED_COLLISION[:2]
    packets = tuple(
        second_moment_packet(adjacency(7, mask), True)
        for mask in (first_mask, second_mask)
    )
    first, second = packets
    collision_packet = (
        first_mask,
        second_mask,
        first["P"],
        first["deck"],
        second["deck"],
        score_sequence(adjacency(7, first_mask)),
        score_sequence(adjacency(7, second_mask)),
        EXPECTED_COLLISION[-1],
    )
    require(collision_packet == EXPECTED_COLLISION, collision_packet)
    require(first["P"] == second["P"], "collision P")
    require(first["gram"] != second["gram"], "collision Gram must differ")

    factor = EXPECTED_COLLISION[-1]
    expected_gram_difference = bivariate_scale(outer(factor, factor), 128)
    actual_gram_difference = bivariate_add(
        first["gram"], bivariate_scale(second["gram"], -1)
    )
    require(
        actual_gram_difference == expected_gram_difference,
        (actual_gram_difference, expected_gram_difference),
    )
    expected_second_difference = bivariate_scale(outer(factor, factor), -256)
    actual_second_difference = bivariate_add(
        first["raw_second"], bivariate_scale(second["raw_second"], -1)
    )
    require(
        actual_second_difference == expected_second_difference,
        (actual_second_difference, expected_second_difference),
    )
    digest.update(repr(collision_packet).encode())
    digest.update(repr(expected_gram_difference).encode())
    digest.update(repr(expected_second_difference).encode())
    return profile_counts, collision_packet, digest.hexdigest()


def pair_from_tournament(matrix):
    core = centered(matrix)
    denominator = char_polynomial_i_minus(core)
    adjugate = adjugate_coefficients(core, denominator)
    ones = (1,) * len(matrix)
    numerator = observer_numerator(adjugate, ones)
    return denominator, numerator


def order_join_audit():
    marked = {}
    for order in range(1, 5):
        cases = []
        for mask in range(1 << (order * (order - 1) // 2)):
            matrix = adjacency(order, mask)
            core = centered(matrix)
            denominator = char_polynomial_i_minus(core)
            adjugate = adjugate_coefficients(core, denominator)
            for tail in itertools.product((1, -1), repeat=order - 1):
                signs = (1, *tail)
                selected = switch(matrix, signs)
                numerator = observer_numerator(adjugate, signs)
                require(
                    pair_from_tournament(selected) == (denominator, numerator),
                    (order, mask, signs, "selected target pair"),
                )
                cases.append((mask, signs, selected, denominator, numerator))
        marked[order] = tuple(cases)

    digest = hashlib.sha256()
    join_cases = 0
    for left_order in range(1, 5):
        for right_order in range(1, 5):
            if left_order + right_order > 5:
                continue
            for left_case in marked[left_order]:
                for right_case in marked[right_order]:
                    left_matrix, left_p, left_n = left_case[2:]
                    right_matrix, right_p, right_n = right_case[2:]
                    joined = order_join(left_matrix, right_matrix)
                    direct_p, direct_n = pair_from_tournament(joined)
                    predicted_p = polynomial_add(
                        polynomial_multiply(left_p, right_p),
                        polynomial_shift(polynomial_multiply(left_n, right_n), 2),
                    )
                    predicted_n = polynomial_add(
                        polynomial_multiply(left_n, right_p),
                        polynomial_multiply(left_p, right_n),
                    )
                    require(
                        (direct_p, direct_n) == (predicted_p, predicted_n),
                        (left_case[:2], right_case[:2], direct_p, direct_n),
                    )
                    direct_q = tuple(
                        direct_p[degree]
                        - (direct_n[degree - 1] if degree else 0)
                        for degree in range(len(direct_p))
                    )
                    left_q = tuple(
                        left_p[degree]
                        - (left_n[degree - 1] if degree else 0)
                        for degree in range(len(left_p))
                    )
                    right_q = tuple(
                        right_p[degree]
                        - (right_n[degree - 1] if degree else 0)
                        for degree in range(len(right_p))
                    )
                    require(
                        direct_q == polynomial_multiply(left_q, right_q),
                        (left_case[:2], right_case[:2], "Q multiplicativity"),
                    )
                    digest.update(
                        f"{left_case[:2]}|{right_case[:2]}|{direct_p}|{direct_n}\n".encode()
                    )
                    join_cases += 1
    require(join_cases == EXPECTED_JOIN_CASES, join_cases)

    singleton = ((0,),)
    cycle = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    source_join = order_join(singleton, cycle)
    sink_join = order_join(cycle, singleton)
    source_pair = pair_from_tournament(source_join)
    sink_pair = pair_from_tournament(sink_join)
    hostile = (
        source_pair[0],
        source_pair[1],
        score_sequence(source_join),
        score_sequence(sink_join),
    )
    require(source_pair == sink_pair, (source_pair, sink_pair))
    require(hostile == EXPECTED_JOIN_HOSTILE, hostile)
    digest.update(repr(hostile).encode())
    return join_cases, hostile, digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    moment_digest = hashlib.sha256()
    labelled_counts = []
    total_tournaments = 0
    total_cuts = 0
    for order in range(1, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        order_cuts = 0
        for mask in range(tournament_count):
            matrix = adjacency(order, mask)
            require(validate_tournament(matrix), (order, mask))
            packet = second_moment_packet(matrix, True)
            relabel_control(matrix, packet)
            record = (
                order,
                mask,
                packet["P"],
                packet["M"],
                packet["deck"],
                packet["gram"],
                packet["raw_second"],
                packet["determinant_raw_second"],
            )
            moment_digest.update(repr(record).encode())
            order_cuts += len(packet["cuts"])
            total_tournaments += 1
        labelled_counts.append((order, tournament_count, order_cuts))
        total_cuts += order_cuts
    labelled_counts = tuple(labelled_counts)
    require(labelled_counts == EXPECTED_LABELLED_MOMENT_COUNTS, labelled_counts)
    triangle = triangle_control()
    moment_digest.update(repr(triangle).encode())
    moment_sha = moment_digest.hexdigest()

    profile_counts, collision, collision_sha = minimal_sidecar_collision()
    join_cases, join_hostile, join_sha = order_join_audit()

    if EXPECTED_MOMENT_SHA256 is not None:
        require(moment_sha == EXPECTED_MOMENT_SHA256, moment_sha)
    if EXPECTED_COLLISION_SHA256 is not None:
        require(collision_sha == EXPECTED_COLLISION_SHA256, collision_sha)
    if EXPECTED_JOIN_SHA256 is not None:
        require(join_sha == EXPECTED_JOIN_SHA256, join_sha)

    semantic_packet = (
        labelled_counts,
        total_tournaments,
        total_cuts,
        triangle,
        profile_counts,
        collision,
        join_cases,
        join_hostile,
        moment_sha,
        collision_sha,
        join_sha,
    )
    semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)

    lines = [
        "Tournament switching coronal second moment and order-join compiler",
        (
            "objects=P(z)=det(I-zC);B(z)=adj(I-zC);N_d(z)=dT_B(z)_d;"
            "M(z)=tr_B(z)=nP-zPprime;D(z,w)=sum_i_Bii(z)_Bii(w)"
        ),
        (
            "second_moment=E_d[N_d(z)N_d(w)]=M(z)M(w)+Ksame(z,w)+"
            "Ktranspose(z,w)-2D(z,w),where_Ksame=tr(B(z)B(w))_and_"
            "Ktranspose=tr(B(z)B(w)^T)"
        ),
        (
            "spectral_closure=(z-w)Ksame=zM(z)P(w)-wP(z)M(w);"
            "(z+w+2zw)Ktranspose=zM(z)P(w)+wP(z)M(w);"
            "therefore_only_D_beyond_P_survives_the_second_average"
        ),
        (
            "variance=Ksame(z,z)+Ktranspose(z,z)-2D(z,z);"
            "Ksame(z,z)=P(M+zMprime)-zMPprime;"
            "Ktranspose(z,z)=PM/(1+z)"
        ),
        (
            "determinant_second_moment=for_Q_d=P-zN_d,"
            "E_Q_is_(1-nz)P+z^2Pprime_and_Cov(Q_d(z),Q_d(w))="
            "zw_Cov(N_d(z),N_d(w))"
        ),
        (
            "sidecar_type=D_is_the_Gram_kernel_of_the_vertex_deleted_"
            "centered_characteristic_polynomials_Bii;it_is_switching_and_"
            "relabel_invariant_but_forgets_the_individual_deck_assignment_"
            "individual_degree_two_Walsh_coefficients_and_higher_distribution_moments"
        ),
        (
            f"labelled_audit=orders_1_to_5;tournaments:{total_tournaments};"
            f"cuts:{total_cuts};counts:{labelled_counts};moment_sha256:{moment_sha}"
        ),
        f"transitive_C3_control={triangle}",
        (
            f"minimal_P_without_D_collision=no_collision_all_labelled_orders_1_to_6:"
            f"{profile_counts};first_order_7_masks:{collision[0]},{collision[1]};"
            f"common_P:{collision[2]};score_sequences:{collision[5]},{collision[6]};"
            "D_first_minus_second=128*z^4*w^4*(1+2z+4z^2)*(1+2w+4w^2);"
            "second_moment_first_minus_second_is_minus_twice_this;"
            f"collision_sha256:{collision_sha}"
        ),
        (
            "order_join_law=for_selected_targets_S_i=T_i^d_i,"
            "P_(S1_join_S2)=P1P2+z^2N1N2_and_"
            "N_(S1_join_S2)=N1P2+P1N2;equivalently_"
            "P_plus_or_minus_zN_multiplies"
        ),
        (
            f"join_audit=all_labelled_marked_pairs_total_order_at_most_5;"
            f"cases:{join_cases};join_sha256:{join_sha}"
        ),
        (
            f"join_loss_hostile=singleton_join_C3_versus_C3_join_singleton_have_"
            f"common_P:{join_hostile[0]};common_N:{join_hostile[1]};"
            f"different_scores:{join_hostile[2]}_versus_{join_hostile[3]};"
            "the_commutative_pair_law_forgets_SCC_order_and_source_versus_sink"
        ),
        (
            "nonconsequences=neither_second_moment_nor_join_pair_restores_cut_"
            "labels_higher_distribution_moments_Hamiltonian_or_path_cover_profiles_SCC_"
            "order_substitution_run_content_or_chronology"
        ),
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic_sha}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
