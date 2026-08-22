#!/usr/bin/env python3
"""Exact order-join compiler for tournament vertex-deletion Gram responses.

For a tournament T put

    s_T(z) = (P_T(z), z N_T(z))^T,
    r_{T,v}(z) = s_{T-v}(z),
    Gamma_T(z,w) = sum_v r_{T,v}(z) r_{T,v}(w)^T.

The split product (a,b)*(c,d)=(ac+bd,ad+bc) is the order-join law for s.
Deleting a vertex before joining proves the corresponding conjugation law
for Gamma.  Its (0,0) entry is the deletion-polynomial Gram D_T.

The audit is assertion-independent and uses exact integer arithmetic only.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
from pathlib import Path


EXPECTED_JOIN_CASES = 2553
EXPECTED_REDUCED_PROFILE_COUNTS = (
    (1, 1, 1),
    (2, 2, 1),
    (3, 8, 2),
    (4, 64, 3),
    (5, 1024, 9),
)
EXPECTED_JOIN_SHA256 = "d78fb46977680cb7b8dba9d0d1c0f77614c866c3768a084321f011e23cc3e8b8"
EXPECTED_HOSTILE_SHA256 = "24654c5a8dc4d41ec641cbd074c8afad1c51a03a6cfdf35ad7e22b9045b8a9e3"
EXPECTED_SEMANTIC_SHA256 = "e4e973e83afebe69d28f7da009a4cb5e75c87c53d2d7668b0b24b71666f42fe8"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def identity(order):
    return tuple(
        tuple(int(row == column) for column in range(order))
        for row in range(order)
    )


def matrix_product(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0]) if middle else 0
    return tuple(
        tuple(
            sum(left[row][index] * right[index][column]
                for index in range(middle))
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


def trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def polynomial_add(*polynomials):
    size = max((len(polynomial) for polynomial in polynomials), default=0)
    return tuple(
        sum(polynomial[index] if index < len(polynomial) else 0
            for polynomial in polynomials)
        for index in range(size)
    )


def polynomial_scale(polynomial, scalar):
    return tuple(scalar * value for value in polynomial)


def polynomial_multiply(left, right):
    if not left or not right:
        return ()
    result = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            result[left_degree + right_degree] += left_value * right_value
    return tuple(result)


def polynomial_pad(polynomial, size):
    require(len(polynomial) <= size, (polynomial, size))
    return tuple(polynomial) + (0,) * (size - len(polynomial))


def polynomial_z(polynomial):
    return (0,) + tuple(polynomial)


def outer(left, right):
    return tuple(tuple(a * b for b in right) for a in left)


def kernel_add(*kernels):
    rows = max((len(kernel) for kernel in kernels), default=0)
    columns = max(
        (len(kernel[0]) if kernel else 0 for kernel in kernels),
        default=0,
    )
    return tuple(
        tuple(
            sum(
                kernel[row][column]
                if row < len(kernel) and column < len(kernel[row]) else 0
                for kernel in kernels
            )
            for column in range(columns)
        )
        for row in range(rows)
    )


def kernel_scale(kernel, scalar):
    return tuple(tuple(scalar * value for value in row) for row in kernel)


def kernel_separable_multiply(kernel, left, right):
    if not kernel:
        return ()
    result = [
        [0] * (len(kernel[0]) + len(right) - 1)
        for _row in range(len(kernel) + len(left) - 1)
    ]
    for row in range(len(kernel)):
        for column in range(len(kernel[row])):
            value = kernel[row][column]
            for left_degree, left_value in enumerate(left):
                for right_degree, right_value in enumerate(right):
                    result[row + left_degree][column + right_degree] += (
                        value * left_value * right_value
                    )
    return tuple(tuple(row) for row in result)


def char_polynomial_i_minus(matrix):
    """Return det(I-zM) from Newton traces with exact divisions."""
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
    """Independent permutation expansion of det(I-zM)."""
    order = len(matrix)
    if order == 0:
        return (1,)
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
    """Coefficient matrices of adj(I-zM), including a CH terminal check."""
    order = len(matrix)
    if order == 0:
        return ()
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
    require(
        terminal == matrix_scale(identity(order), -denominator[-1]),
        (matrix, denominator, terminal),
    )
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


def delete_vertex(matrix, vertex):
    order = len(matrix)
    return tuple(
        tuple(matrix[row][column] for column in range(order) if column != vertex)
        for row in range(order)
        if row != vertex
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


def score_sequence(matrix):
    return tuple(sorted((sum(row) for row in matrix), reverse=True))


def relabel(matrix, permutation):
    return tuple(
        tuple(matrix[permutation[row]][permutation[column]]
              for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def canonical_matrix(matrix):
    return min(
        relabel(matrix, permutation)
        for permutation in itertools.permutations(range(len(matrix)))
    )


def pair_from_tournament(matrix, independent=False):
    order = len(matrix)
    if order == 0:
        return (1,), (0,)
    core = centered(matrix)
    denominator = char_polynomial_i_minus(core)
    if independent:
        require(
            denominator == determinant_polynomial_permutation(core),
            (matrix, denominator, "independent determinant"),
        )
    adjugate = adjugate_coefficients(core, denominator)
    numerator = tuple(
        sum(coefficient[row][column]
            for row in range(order) for column in range(order))
        for coefficient in adjugate
    )
    return denominator, polynomial_z(numerator)


def split_product(left, right):
    return (
        polynomial_add(
            polynomial_multiply(left[0], right[0]),
            polynomial_multiply(left[1], right[1]),
        ),
        polynomial_add(
            polynomial_multiply(left[0], right[1]),
            polynomial_multiply(left[1], right[0]),
        ),
    )


def response_packet(matrix, independent=False):
    order = len(matrix)
    response = pair_from_tournament(matrix, independent)
    size = order + 1
    response = tuple(polynomial_pad(polynomial, size) for polynomial in response)
    deletion_responses = []
    for vertex in range(order):
        deleted = pair_from_tournament(delete_vertex(matrix, vertex), independent)
        deletion_responses.append(
            tuple(polynomial_pad(polynomial, size) for polynomial in deleted)
        )
    gram = tuple(
        tuple(
            kernel_add(
                *(outer(deleted[left], deleted[right])
                  for deleted in deletion_responses)
            )
            for right in range(2)
        )
        for left in range(2)
    )
    return response, tuple(deletion_responses), gram


def predicted_gram(left_gram, right_gram, left_response, right_response):
    """The two L_s Gamma L_s^T summands, coefficientwise in z,w."""
    result = []
    for output_left in range(2):
        row = []
        for output_right in range(2):
            terms = []
            for inner_left in range(2):
                for inner_right in range(2):
                    terms.append(
                        kernel_separable_multiply(
                            left_gram[inner_left][inner_right],
                            right_response[output_left ^ inner_left],
                            right_response[output_right ^ inner_right],
                        )
                    )
                    terms.append(
                        kernel_separable_multiply(
                            right_gram[inner_left][inner_right],
                            left_response[output_left ^ inner_left],
                            left_response[output_right ^ inner_right],
                        )
                    )
            row.append(kernel_add(*terms))
        result.append(tuple(row))
    return tuple(result)


def direct_join_audit():
    tournaments = {
        order: tuple(
            (mask, adjacency(order, mask))
            for mask in range(1 << (order * (order - 1) // 2))
        )
        for order in range(1, 6)
    }
    packets = {
        (order, mask): response_packet(matrix)
        for order, cases in tournaments.items()
        for mask, matrix in cases
    }
    digest = hashlib.sha256()
    cases = 0
    reverse_matches = 0
    for left_order in range(1, 6):
        for right_order in range(1, 6):
            if left_order + right_order > 6:
                continue
            for left_mask, left in tournaments[left_order]:
                left_response, left_deletions, left_gram = packets[
                    (left_order, left_mask)
                ]
                for right_mask, right in tournaments[right_order]:
                    right_response, right_deletions, right_gram = packets[
                        (right_order, right_mask)
                    ]
                    joined = order_join(left, right)
                    direct_response, direct_deletions, direct_gram = response_packet(
                        joined,
                        independent=True,
                    )
                    expected_response = split_product(left_response, right_response)
                    expected_deletions = tuple(
                        split_product(deleted, right_response)
                        for deleted in left_deletions
                    ) + tuple(
                        split_product(left_response, deleted)
                        for deleted in right_deletions
                    )
                    expected_gram = predicted_gram(
                        left_gram,
                        right_gram,
                        left_response,
                        right_response,
                    )
                    require(
                        direct_response == expected_response,
                        (left_order, left_mask, right_order, right_mask, "response"),
                    )
                    require(
                        direct_deletions == expected_deletions,
                        (left_order, left_mask, right_order, right_mask, "deletions"),
                    )
                    require(
                        direct_gram == expected_gram,
                        (left_order, left_mask, right_order, right_mask, "gram"),
                    )
                    reverse_packet = response_packet(order_join(right, left))
                    require(
                        (direct_response, direct_gram)
                        == (reverse_packet[0], reverse_packet[2]),
                        (left_order, left_mask, right_order, right_mask, "order loss"),
                    )
                    reverse_matches += 1
                    digest.update(
                        repr((
                            left_order,
                            left_mask,
                            right_order,
                            right_mask,
                            direct_response,
                            direct_gram,
                        )).encode()
                    )
                    cases += 1
    require(cases == EXPECTED_JOIN_CASES, cases)
    require(reverse_matches == EXPECTED_JOIN_CASES, reverse_matches)
    return cases, reverse_matches, digest.hexdigest()


def reduced_interface_hostile():
    profile_counts = []
    for order in range(1, 6):
        profiles = {}
        tournament_count = 1 << (order * (order - 1) // 2)
        for mask in range(tournament_count):
            response, _deletions, gram = response_packet(adjacency(order, mask))
            reduced = (response, gram[0][0], gram[0][1], gram[1][0])
            if reduced in profiles:
                require(
                    profiles[reduced] == gram[1][1],
                    (order, mask, "premature reduced-interface collision"),
                )
            else:
                profiles[reduced] = gram[1][1]
        profile_counts.append((order, tournament_count, len(profiles)))
    profile_counts = tuple(profile_counts)
    if EXPECTED_REDUCED_PROFILE_COUNTS is not None:
        require(
            profile_counts == EXPECTED_REDUCED_PROFILE_COUNTS,
            profile_counts,
        )

    first_mask, second_mask = 73, 83
    first_matrix = adjacency(6, first_mask)
    second_matrix = adjacency(6, second_mask)
    first = response_packet(first_matrix, independent=True)
    second = response_packet(second_matrix, independent=True)
    first_reduced = (first[0], first[2][0][0], first[2][0][1], first[2][1][0])
    second_reduced = (
        second[0], second[2][0][0], second[2][0][1], second[2][1][0]
    )
    require(first_reduced == second_reduced, "hostile reduced interfaces")
    require(first[2][1][1] != second[2][1][1], "hostile missing block")
    factor = (0, 0, 0, 1, 2, 4, 0)
    expected_missing_difference = kernel_scale(outer(factor, factor), 128)
    actual_missing_difference = kernel_add(
        first[2][1][1], kernel_scale(second[2][1][1], -1)
    )
    require(
        actual_missing_difference == expected_missing_difference,
        actual_missing_difference,
    )

    singleton = ((0,),)
    joined = tuple(
        response_packet(order_join(adjacency(6, mask), singleton), independent=True)
        for mask in (first_mask, second_mask)
    )
    joined_d_difference = kernel_add(
        joined[0][2][0][0], kernel_scale(joined[1][2][0][0], -1)
    )
    joined_factor = (0, 0, 0, 0, 1, 2, 4, 0)
    expected_joined_difference = kernel_scale(outer(joined_factor, joined_factor), 128)
    require(joined_d_difference == expected_joined_difference, joined_d_difference)

    cycle = (
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 0),
    )
    source = order_join(singleton, cycle)
    sink = order_join(cycle, singleton)
    source_packet = response_packet(source, independent=True)
    sink_packet = response_packet(sink, independent=True)
    require(
        (source_packet[0], source_packet[2])
        == (sink_packet[0], sink_packet[2]),
        "order-loss interface",
    )
    require(
        tuple(sorted(source_packet[1])) == tuple(sorted(sink_packet[1])),
        "order-loss marked deletion multiset",
    )
    require(score_sequence(source) != score_sequence(sink), "order-loss scores")
    lower_order_reverse_cases = 0
    for left_order in range(1, 3):
        for right_order in range(1, 3):
            if left_order + right_order > 3:
                continue
            for left_mask in range(1 << (left_order * (left_order - 1) // 2)):
                for right_mask in range(1 << (right_order * (right_order - 1) // 2)):
                    left = adjacency(left_order, left_mask)
                    right = adjacency(right_order, right_mask)
                    require(
                        canonical_matrix(order_join(left, right))
                        == canonical_matrix(order_join(right, left)),
                        (left_order, left_mask, right_order, right_mask),
                    )
                    lower_order_reverse_cases += 1

    hostile = (
        profile_counts,
        first_mask,
        second_mask,
        first[0],
        score_sequence(first_matrix),
        score_sequence(second_matrix),
        factor,
        joined_factor,
        source_packet[0],
        source_packet[2],
        score_sequence(source),
        score_sequence(sink),
        lower_order_reverse_cases,
    )
    hostile_sha = hashlib.sha256(repr(hostile).encode()).hexdigest()
    if EXPECTED_HOSTILE_SHA256 is not None:
        require(hostile_sha == EXPECTED_HOSTILE_SHA256, hostile_sha)
    return hostile, hostile_sha


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

    join_cases, reverse_matches, join_sha = direct_join_audit()
    hostile, hostile_sha = reduced_interface_hostile()
    profile_counts = hostile[0]
    semantic_packet = (
        join_cases,
        reverse_matches,
        join_sha,
        hostile,
        hostile_sha,
    )
    semantic_sha = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_JOIN_SHA256 is not None:
        require(join_sha == EXPECTED_JOIN_SHA256, join_sha)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_sha == EXPECTED_SEMANTIC_SHA256, semantic_sha)

    lines = [
        "Tournament deletion-response Gram order-join compiler",
        (
            "objects=s_T(z)=(P_T(z),zN_T(z));r_Tv(z)=s_(T-v)(z)="
            "(p_v(z),znu_v(z));Gamma_T(z,w)=sum_v_r_Tv(z)_r_Tv(w)^T;"
            "D_T=Gamma_T[0,0]"
        ),
        (
            "split_product=(a,b)*(c,d)=(ac+bd,ad+bc);"
            "s_(X_join_Y)=s_X*s_Y"
        ),
        (
            "deletion_response=for_v_in_X,r_(X_join_Y),v=r_Xv*s_Y;"
            "for_v_in_Y,r_(X_join_Y),v=s_X*r_Yv"
        ),
        (
            "gram_law=Gamma_(X_join_Y)=L_sY(z)_Gamma_X(z,w)_L_sY(w)^T+"
            "L_sX(z)_Gamma_Y(z,w)_L_sX(w)^T;"
            "L_(a,b)=((a,b),(b,a))"
        ),
        (
            "minimal_interface=two_response_coordinates_and_three_Gram_kernels:"
            "D=sum_pp;E=sum_p_times_znu;F=sum_znu_times_znu;"
            "Gamma=((D,E),(E_transposed,F))"
        ),
        (
            f"labelled_join_audit=all_ordered_factor_pairs_total_order_at_most_6;"
            f"cases:{join_cases};direct_deletion_vectors_and_full_Gram_match;"
            f"reverse_order_response_and_Gram_matches:{reverse_matches};"
            f"join_sha256:{join_sha}"
        ),
        (
            f"sharp_missing_coordinate_audit=no_reduced_interface_collision_orders_1_to_5:"
            f"{profile_counts};first_order_6_masks:73,83;common_s_D_E_Etranspose;"
            "F_73_minus_F_83=128*z^3*w^3*(1+2z+4z^2)*(1+2w+4w^2);"
            "after_joining_K1,D_difference=128*z^4*w^4*(1+2z+4z^2)*"
            "(1+2w+4w^2)"
        ),
        (
            f"missing_coordinate_scores={hostile[4]}_versus_{hostile[5]};"
            "therefore_(P,N,D)_and_even_the_mixed_deletion_Gram_do_not_compose_D"
        ),
        (
            f"sharp_order_loss=K1_join_C3_versus_C3_join_K1_have_equal_s_and_Gamma;"
            f"scores:{hostile[10]}_versus_{hostile[11]};"
            f"all_reverse_joins_total_order_at_most_3_isomorphic_cases:{hostile[12]}"
        ),
        (
            "lost_information=the_interface_forgets_vertex_assignment_of_deletion_"
            "responses_and_the_ordered_SCC_condensation;no_isomorphism_Hamiltonian_"
            "path_cover_substitution_or_chronology_consequence"
        ),
        f"hostile_sha256={hostile_sha}",
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
