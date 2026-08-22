#!/usr/bin/env python3
"""Exact audit for the skew deletion-response orientation current (THM-3369).

The scalar/Gram ordered-join interface of THM-3324 is commutative.  This
companion contracts the same marked deletion responses against the intrinsic
skew adjacency K=A-A^T:

    Omega_T(z,w) = sum_(u,v) K_uv r_Tu(z) r_Tv(w)^T.

It verifies the Euler first-response identity, the noncommutative two-factor
and repeated-factor ordered-join laws, relabel/reversal covariance, and the
sharp K1/C3 order hostile.  Arithmetic is exact and assertion-independent.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
from pathlib import Path

from tournament_deletion_response_gram_order_join_compiler_20260803 import (
    adjacency,
    delete_vertex,
    kernel_add,
    kernel_scale,
    kernel_separable_multiply,
    order_join,
    outer,
    pair_from_tournament,
    polynomial_add,
    polynomial_multiply,
    polynomial_pad,
    require,
    response_packet,
    score_sequence,
    split_product,
)


EXPECTED_JOIN_CASES = 2553
EXPECTED_TRIPLE_CASES = 339
EXPECTED_JOIN_SHA256 = "02db9efeb5ad40e9629c61ce5e4c8cc4e5487d5894f9b1989cdd3884798e3a79"
EXPECTED_COVARIANCE_SHA256 = "2742c1fd54b4f34b6c894ce4d22c01888cab4254cf0a8e41496830a9af699fb1"
EXPECTED_HOSTILE_SHA256 = "1a69962011e3f18906f6b1aea516c840460302d4d1c1037a1b15b329d4e3ce0d"
EXPECTED_SEMANTIC_SHA256 = "cd70e45bc435acb281084bdee703dcf36e7b53b6877003b8260a5ff40671f12e"


def zero_kernel(size):
    return tuple(tuple(0 for _column in range(size)) for _row in range(size))


def zero_kernel_matrix(size):
    zero = zero_kernel(size)
    return ((zero, zero), (zero, zero))


def kernel_matrix_add(*matrices):
    return tuple(
        tuple(kernel_add(*(matrix[row][column] for matrix in matrices))
              for column in range(2))
        for row in range(2)
    )


def kernel_matrix_scale(matrix, scalar):
    return tuple(
        tuple(kernel_scale(matrix[row][column], scalar) for column in range(2))
        for row in range(2)
    )


def vector_outer(left, right):
    return tuple(
        tuple(outer(left[row], right[column]) for column in range(2))
        for row in range(2)
    )


def skew_adjacency(matrix):
    order = len(matrix)
    return tuple(
        tuple(matrix[row][column] - matrix[column][row]
              for column in range(order))
        for row in range(order)
    )


def omega_from_deletions(matrix, deletions):
    size = len(matrix) + 1
    skew = skew_adjacency(matrix)
    result = [[zero_kernel(size) for _column in range(2)] for _row in range(2)]
    for row in range(2):
        for column in range(2):
            terms = []
            for left in range(len(matrix)):
                for right in range(len(matrix)):
                    if skew[left][right]:
                        terms.append(
                            kernel_scale(
                                outer(deletions[left][row], deletions[right][column]),
                                skew[left][right],
                            )
                        )
            result[row][column] = kernel_add(*terms) if terms else zero_kernel(size)
    return tuple(tuple(row) for row in result)


def direct_packet(matrix, independent=False):
    response, deletions, _gram = response_packet(matrix, independent=independent)
    omega = omega_from_deletions(matrix, deletions)
    return response, deletions, omega


def deletion_sum(deletions):
    if not deletions:
        return ((0,), (0,))
    return tuple(
        tuple(sum(deletion[coordinate][degree] for deletion in deletions)
              for degree in range(len(deletions[0][coordinate])))
        for coordinate in range(2)
    )


def euler_first_response(order, response):
    """n*s-z*s', coefficientwise."""
    return tuple(
        tuple((order - degree) * value for degree, value in enumerate(polynomial))
        for polynomial in response
    )


def h_apply_vector(response, vector):
    return tuple(
        polynomial_add(
            *(polynomial_multiply(response[output ^ inner], vector[inner])
              for inner in range(2))
        )
        for output in range(2)
    )


def conjugate_omega(omega, response):
    result = []
    for output_left in range(2):
        row = []
        for output_right in range(2):
            row.append(
                kernel_add(
                    *(
                        kernel_separable_multiply(
                            omega[inner_left][inner_right],
                            response[output_left ^ inner_left],
                            response[output_right ^ inner_right],
                        )
                        for inner_left in range(2)
                        for inner_right in range(2)
                    )
                )
            )
        result.append(tuple(row))
    return tuple(result)


def predicted_two_factor(left_packet, right_packet):
    left_response, left_deletions, left_omega = left_packet
    right_response, right_deletions, right_omega = right_packet
    left_q = deletion_sum(left_deletions)
    right_q = deletion_sum(right_deletions)
    left_in_join = h_apply_vector(right_response, left_q)
    right_in_join = h_apply_vector(left_response, right_q)
    internal = kernel_matrix_add(
        conjugate_omega(left_omega, right_response),
        conjugate_omega(right_omega, left_response),
    )
    cross = kernel_matrix_add(
        vector_outer(left_in_join, right_in_join),
        kernel_matrix_scale(vector_outer(right_in_join, left_in_join), -1),
    )
    return kernel_matrix_add(internal, cross)


def split_identity():
    return ((1,), (0,))


def split_product_many(responses):
    result = split_identity()
    for response in responses:
        result = split_product(result, response)
    return result


def predicted_many_factor(packets):
    responses = [packet[0] for packet in packets]
    qs = [deletion_sum(packet[1]) for packet in packets]
    transformed_qs = []
    terms = []
    for index, packet in enumerate(packets):
        others = split_product_many(
            response for other, response in enumerate(responses) if other != index
        )
        transformed_qs.append(h_apply_vector(others, qs[index]))
        terms.append(conjugate_omega(packet[2], others))
    for left in range(len(packets)):
        for right in range(left + 1, len(packets)):
            terms.append(vector_outer(transformed_qs[left], transformed_qs[right]))
            terms.append(
                kernel_matrix_scale(
                    vector_outer(transformed_qs[right], transformed_qs[left]), -1
                )
            )
    return kernel_matrix_add(*terms)


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix)))


def converse(matrix):
    return transpose(matrix)


def relabel(matrix, permutation):
    return tuple(
        tuple(matrix[permutation[row]][permutation[column]]
              for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def omega_alternating(omega):
    for row in range(2):
        for column in range(2):
            left = omega[row][column]
            right = omega[column][row]
            if left != tuple(tuple(-right[j][i] for j in range(len(right)))
                             for i in range(len(right[0]) if right else 0)):
                return False
    return True


def tournaments_through(max_order):
    return {
        order: tuple(
            (mask, adjacency(order, mask))
            for mask in range(1 << (order * (order - 1) // 2))
        )
        for order in range(1, max_order + 1)
    }


def covariance_audit(tournaments, packets):
    digest = hashlib.sha256()
    cases = 0
    relabel_cases = 0
    for order, entries in tournaments.items():
        for mask, matrix in entries:
            response, deletions, omega = packets[(order, mask)]
            require(deletion_sum(deletions) == euler_first_response(order, response),
                    (order, mask, "Euler first response"))
            require(omega_alternating(omega), (order, mask, "alternating"))
            opposite_response, _opposite_deletions, opposite_omega = direct_packet(
                converse(matrix)
            )
            require(opposite_response == response, (order, mask, "converse response"))
            require(opposite_omega == kernel_matrix_scale(omega, -1),
                    (order, mask, "converse omega"))
            digest.update(repr((order, mask, response, omega)).encode())
            cases += 1
            if order <= 4:
                for permutation in itertools.permutations(range(order)):
                    relabelled = direct_packet(relabel(matrix, permutation))
                    require(relabelled[0] == response, (order, mask, permutation, "s"))
                    require(relabelled[2] == omega, (order, mask, permutation, "omega"))
                    relabel_cases += 1
    return cases, relabel_cases, digest.hexdigest()


def two_factor_audit(tournaments, packets):
    digest = hashlib.sha256()
    cases = 0
    reverse_equal_small = 0
    for left_order in range(1, 6):
        for right_order in range(1, 6):
            if left_order + right_order > 6:
                continue
            for left_mask, left in tournaments[left_order]:
                left_packet = packets[(left_order, left_mask)]
                for right_mask, right in tournaments[right_order]:
                    right_packet = packets[(right_order, right_mask)]
                    joined = order_join(left, right)
                    direct = direct_packet(joined)
                    require(
                        direct[0] == split_product(left_packet[0], right_packet[0]),
                        (left_order, left_mask, right_order, right_mask, "response"),
                    )
                    expected = predicted_two_factor(left_packet, right_packet)
                    require(direct[2] == expected,
                            (left_order, left_mask, right_order, right_mask, "omega"))
                    reverse = direct_packet(order_join(right, left))
                    if left_order + right_order <= 3:
                        require(reverse[2] == direct[2],
                                (left_order, left_mask, right_order, right_mask,
                                 "small reverse"))
                        reverse_equal_small += 1
                    digest.update(
                        repr((left_order, left_mask, right_order, right_mask,
                              direct[0], direct[2], reverse[2])).encode()
                    )
                    cases += 1
    return cases, reverse_equal_small, digest.hexdigest()


def triple_audit(tournaments, packets):
    digest = hashlib.sha256()
    cases = 0
    for orders in itertools.product(range(1, 5), repeat=3):
        if sum(orders) > 6:
            continue
        banks = [tournaments[order] for order in orders]
        for entries in itertools.product(*banks):
            matrices = [entry[1] for entry in entries]
            factor_packets = [packets[(orders[index], entries[index][0])]
                              for index in range(3)]
            joined = order_join(order_join(matrices[0], matrices[1]), matrices[2])
            direct = direct_packet(joined)
            require(
                direct[0] == split_product_many(packet[0] for packet in factor_packets),
                (orders, entries, "triple response"),
            )
            expected = predicted_many_factor(factor_packets)
            require(direct[2] == expected, (orders, entries, "triple omega"))
            digest.update(repr((orders, tuple(entry[0] for entry in entries),
                                direct[2])).encode())
            cases += 1
    return cases, digest.hexdigest()


def hostile_audit():
    singleton = adjacency(1, 0)
    cycle = adjacency(3, 5)
    forward = order_join(singleton, cycle)
    backward = order_join(cycle, singleton)
    forward_packet = direct_packet(forward, independent=True)
    backward_packet = direct_packet(backward, independent=True)
    require(forward_packet[0] == backward_packet[0], "hostile common response")
    require(forward_packet[2] == kernel_matrix_scale(backward_packet[2], -1),
            "hostile opposite omega")
    require(forward_packet[2] != backward_packet[2], "hostile nonzero omega")
    factor_z = (1, 3, 6, 4)
    factor_w = (0, 0, 0, -24)
    expected_01 = outer(polynomial_pad(factor_z, 5), polynomial_pad(factor_w, 5))
    require(forward_packet[2][0][1] == expected_01,
            (forward_packet[2][0][1], expected_01, "hostile factor"))
    packet = (
        forward_packet[0],
        score_sequence(forward),
        score_sequence(backward),
        forward_packet[2],
    )
    digest = hashlib.sha256(repr(packet).encode()).hexdigest()
    return packet, digest


def main():
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    tournaments = tournaments_through(5)
    packets = {
        (order, mask): direct_packet(matrix)
        for order, entries in tournaments.items()
        for mask, matrix in entries
    }
    covariance_cases, relabel_cases, covariance_sha = covariance_audit(
        tournaments, packets
    )
    join_cases, reverse_equal_small, join_sha = two_factor_audit(
        tournaments, packets
    )
    triple_cases, triple_sha = triple_audit(tournaments, packets)
    hostile, hostile_sha = hostile_audit()

    require(join_cases == EXPECTED_JOIN_CASES, join_cases)
    require(triple_cases == EXPECTED_TRIPLE_CASES, triple_cases)
    semantic = (
        covariance_cases,
        relabel_cases,
        covariance_sha,
        join_cases,
        reverse_equal_small,
        join_sha,
        triple_cases,
        triple_sha,
        hostile,
        hostile_sha,
    )
    semantic_sha = hashlib.sha256(repr(semantic).encode()).hexdigest()
    for actual, expected in (
        (join_sha, EXPECTED_JOIN_SHA256),
        (covariance_sha, EXPECTED_COVARIANCE_SHA256),
        (hostile_sha, EXPECTED_HOSTILE_SHA256),
        (semantic_sha, EXPECTED_SEMANTIC_SHA256),
    ):
        if expected is not None:
            require(actual == expected, (actual, expected))

    lines = [
        "THM-3369 SKEW DELETION RESPONSE ORDERED-JOIN ORIENTATION CURRENT",
        (
            "object=Omega_T(z,w)=sum_uv(A_uv-A_vu)r_Tu(z)r_Tv(w)^T;"
            "r_Tv=s_(T-v);s_T=(P_T,zN_T)^T"
        ),
        "first_response=sum_v_r_Tv=n*s_T-z*s_T_prime",
        (
            "two_factor_law=two_internal_H-conjugates+"
            "(H_Y*q_X)(z)(H_X*q_Y)(w)^T-(H_X*q_Y)(z)(H_Y*q_X)(w)^T"
        ),
        "many_factor_law=sum_internal_conjugates+ordered_pairwise_transformed_q_wedges",
        (
            f"covariance_audit=labelled_tournaments_orders_1_to_5:{covariance_cases};"
            f"relabels_through_order_4:{relabel_cases};converse_negates_Omega;"
            f"Omega(w,z)^T=-Omega(z,w);sha256:{covariance_sha}"
        ),
        (
            f"join_audit=all_ordered_factor_pairs_total_order_at_most_6:{join_cases};"
            f"reverse_equal_total_order_at_most_3:{reverse_equal_small};sha256:{join_sha}"
        ),
        (
            f"three_factor_audit=all_ordered_labelled_triples_total_order_at_most_6:"
            f"{triple_cases};sha256:{triple_sha}"
        ),
        (
            "sharp_hostile=K1_join_C3_versus_C3_join_K1;common_s;"
            f"scores:{hostile[1]}_versus_{hostile[2]};Omega_reverse=-Omega_forward;"
            "Omega_forward_01=-24*w^3*(z+1)*(4*z^2+2*z+1)"
        ),
        f"hostile_sha256={hostile_sha}",
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic_sha}",
        "scope=exact_operation_response_and_one_sharp_order_bit;no_injectivity_or_chronology_claim",
        "all_exact_controls=PASS",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    main()
