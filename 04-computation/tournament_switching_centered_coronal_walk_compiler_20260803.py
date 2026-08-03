#!/usr/bin/env python3
"""Exact centered-coronal compiler for tournament switching walk sequences.

For a tournament adjacency matrix A, put C=2A-J.  Reversing every arc across
the cut encoded by d in {+-1}^n sends C to D C D.  If

    P(z)=det(I-zC),
    U_d(z)=d^T(I-zC)^(-1)d=N_d(z)/P(z),

then the switched tournament T^d has total directed-walk series

    G_(T^d)(2z)=N_d(z)/(P(z)-zN_d(z)).

Thus one switching-invariant denominator and one observer numerator compile
the complete walk sequence.  The script derives every object over the
integers and exhausts all labelled tournaments through order five.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
from collections import Counter
from pathlib import Path


EXPECTED_TOURNAMENTS = 1099
EXPECTED_SWITCH_CASES = 16_933
EXPECTED_ORDER_COUNTS = (
    (1, 1, 1),
    (2, 2, 4),
    (3, 8, 32),
    (4, 64, 512),
    (5, 1024, 16_384),
)
EXPECTED_HOSTILE = (
    (1, 3, 6, 4),
    (3, 6, 4),
    (3, 6, 12),
    (1, 0, 0, 0),
    (1, 0, 0, -8),
    (3, 3, 1, 0, 0, 0, 0, 0),
    (3, 3, 3, 3, 3, 3, 3, 3),
)
EXPECTED_CASE_SHA256 = "f5edbe1308c1f43503d0f790acbdabb9b09dbd741d6eed7c3df2ce1eb0057ad4"
EXPECTED_AVERAGE_SHA256 = "8289c5c82f0194084af60048cf055219ada6739fb6d76887603f5477c7fcde6f"
EXPECTED_SEMANTIC_SHA256 = "282a0c09458e9b539e349135f3800cd8bfc850971b4e8782fc14de75ea36c1be"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def add_polynomials(left, right):
    size = max(len(left), len(right))
    return tuple(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def multiply_polynomials(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return tuple(result)


def derivative(polynomial):
    if len(polynomial) <= 1:
        return (0,)
    return tuple(index * polynomial[index] for index in range(1, len(polynomial)))


def determinant_polynomial(matrix, scale=1):
    """Return det(I-z*scale*matrix) by the permutation definition."""
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
            if row == column:
                factor = (1, -scale * matrix[row][column])
            else:
                factor = (0, -scale * matrix[row][column])
            term = multiply_polynomials(term, factor)
        result = add_polynomials(
            result,
            tuple(sign * coefficient for coefficient in term),
        )
    require(result[0] == 1 and len(result) == order + 1, result)
    return result


def adjacency(order, mask):
    matrix = [[0] * order for _ in range(order)]
    edge = 0
    for left in range(order):
        for right in range(left + 1, order):
            orientation = (mask >> edge) & 1
            matrix[left][right] = orientation
            matrix[right][left] = 1 - orientation
            edge += 1
    return tuple(tuple(row) for row in matrix)


def centered(matrix):
    return tuple(
        tuple(2 * matrix[row][column] - 1 for column in range(len(matrix)))
        for row in range(len(matrix))
    )


def switch(matrix, signs):
    order = len(matrix)
    return tuple(
        tuple(
            0 if row == column
            else (
                matrix[row][column]
                if signs[row] == signs[column]
                else 1 - matrix[row][column]
            )
            for column in range(order)
        )
        for row in range(order)
    )


def conjugate(matrix, signs):
    return tuple(
        tuple(
            signs[row] * matrix[row][column] * signs[column]
            for column in range(len(matrix))
        )
        for row in range(len(matrix))
    )


def matrix_vector(matrix, vector):
    return tuple(
        sum(entry * value for entry, value in zip(row, vector))
        for row in matrix
    )


def matrix_product(left, right):
    order = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(order))
            for column in range(order)
        )
        for row in range(order)
    )


def identity(order):
    return tuple(
        tuple(int(row == column) for column in range(order))
        for row in range(order)
    )


def bilinear_moments(matrix, left, right, count):
    current = tuple(right)
    moments = []
    for _index in range(count):
        moments.append(sum(a * b for a, b in zip(left, current)))
        current = matrix_vector(matrix, current)
    return tuple(moments)


def traces_of_powers(matrix, count):
    power = identity(len(matrix))
    traces = []
    for _index in range(count):
        traces.append(sum(power[index][index] for index in range(len(matrix))))
        power = matrix_product(power, matrix)
    return tuple(traces)


def observer_numerator(denominator, moments, order):
    numerator = tuple(
        sum(
            denominator[index] * moments[degree - index]
            for index in range(degree + 1)
        )
        for degree in range(order)
    )
    for degree in range(order, len(moments)):
        coefficient = sum(
            denominator[index] * moments[degree - index]
            for index in range(len(denominator))
        )
        require(coefficient == 0, (degree, coefficient))
    return numerator


def subtract_z_numerator(denominator, numerator):
    return tuple(
        denominator[degree]
        - (numerator[degree - 1] if 1 <= degree <= len(numerator) else 0)
        for degree in range(len(denominator))
    )


def compiled_scaled_walks(observer_moments):
    """Coefficients of U/(1-zU) from the coefficients of U."""
    result = []
    for degree, moment in enumerate(observer_moments):
        value = moment
        if degree:
            value += sum(
                observer_moments[left] * result[degree - 1 - left]
                for left in range(degree)
            )
        result.append(value)
    return tuple(result)


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


def case_packet(order, mask, matrix, signs, depth, denominator):
    core = centered(matrix)
    switched = switch(matrix, signs)
    switched_core = centered(switched)
    require(validate_tournament(switched), (order, mask, signs, "not tournament"))
    require(switched_core == conjugate(core, signs), (order, mask, signs, "core"))
    negative_signs = tuple(-value for value in signs)
    require(switch(matrix, negative_signs) == switched, (order, mask, signs, "gauge"))

    moments = bilinear_moments(core, signs, signs, depth)
    negative_moments = bilinear_moments(
        core, negative_signs, negative_signs, depth
    )
    require(moments == negative_moments, (order, mask, signs, "observer gauge"))
    numerator = observer_numerator(denominator, moments, order)
    switched_denominator = subtract_z_numerator(denominator, numerator)
    direct_denominator = determinant_polynomial(switched, scale=2)
    require(
        switched_denominator == direct_denominator,
        (order, mask, signs, switched_denominator, direct_denominator),
    )

    ones = (1,) * order
    direct_walks = bilinear_moments(switched, ones, ones, depth)
    direct_scaled = tuple(
        (2 ** degree) * value for degree, value in enumerate(direct_walks)
    )
    compiled = compiled_scaled_walks(moments)
    require(direct_scaled == compiled, (order, mask, signs, "walk compiler"))
    for degree in range(order, depth):
        recurrence = sum(
            switched_denominator[index] * direct_scaled[degree - index]
            for index in range(order + 1)
        )
        require(recurrence == 0, (order, mask, signs, degree, recurrence))
    return (
        signs,
        numerator,
        switched_denominator,
        direct_scaled[: order + 2],
    ), moments


def hostile_control():
    matrix = (
        (0, 1, 1),
        (0, 0, 1),
        (0, 0, 0),
    )
    core = centered(matrix)
    denominator = determinant_polynomial(core)
    source_signs = (-1, 1, 1)
    middle_signs = (1, -1, 1)
    records = []
    for signs in (source_signs, middle_signs):
        switched = switch(matrix, signs)
        moments = bilinear_moments(core, signs, signs, 10)
        numerator = observer_numerator(denominator, moments, 3)
        adjacency_denominator = subtract_z_numerator(denominator, numerator)
        walks = bilinear_moments(switched, (1, 1, 1), (1, 1, 1), 8)
        records.append((numerator, adjacency_denominator, walks))
    source, middle = records
    packet = (
        denominator,
        source[0],
        middle[0],
        source[1],
        middle[1],
        source[2],
        middle[2],
    )
    require(packet == EXPECTED_HOSTILE, packet)
    require(
        sum(value < 0 for value in source_signs)
        == sum(value < 0 for value in middle_signs) == 1,
        "hostile cut sizes",
    )
    return packet


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

    case_digest = hashlib.sha256()
    average_records = []
    order_records = []
    total_tournaments = 0
    total_cases = 0
    total_moment_cells = 0
    distinct_target_profiles = Counter()

    for order in range(1, 6):
        tournament_count = 1 << (order * (order - 1) // 2)
        cut_count = 1 << (order - 1)
        order_cases = 0
        for mask in range(tournament_count):
            matrix = adjacency(order, mask)
            require(validate_tournament(matrix), (order, mask))
            core = centered(matrix)
            denominator = determinant_polynomial(core)
            depth = 2 * order + 5
            traces = traces_of_powers(core, depth)
            numerator_sum = [0] * order
            denominator_sum = [0] * (order + 1)
            moment_sums = [0] * depth
            target_profiles = set()

            for tail in itertools.product((1, -1), repeat=order - 1):
                signs = (1, *tail)
                record, moments = case_packet(
                    order, mask, matrix, signs, depth, denominator
                )
                for index, value in enumerate(record[1]):
                    numerator_sum[index] += value
                for index, value in enumerate(record[2]):
                    denominator_sum[index] += value
                for index, value in enumerate(moments):
                    moment_sums[index] += value
                target_profiles.add(record[3])
                case_digest.update(f"{order}|{mask}|{record}\n".encode())
                total_cases += 1
                order_cases += 1
                total_moment_cells += depth

            derivative_denominator = derivative(denominator)
            expected_average_numerator = tuple(
                order * denominator[degree]
                - (
                    derivative_denominator[degree - 1]
                    if 1 <= degree <= len(derivative_denominator)
                    else 0
                )
                for degree in range(order)
            )
            expected_average_denominator = subtract_z_numerator(
                denominator, expected_average_numerator
            )
            require(
                tuple(numerator_sum)
                == tuple(cut_count * value for value in expected_average_numerator),
                (order, mask, "numerator average"),
            )
            require(
                tuple(denominator_sum)
                == tuple(cut_count * value for value in expected_average_denominator),
                (order, mask, "denominator average"),
            )
            require(
                tuple(moment_sums)
                == tuple(cut_count * value for value in traces),
                (order, mask, "moment average"),
            )
            average_record = (
                order,
                mask,
                denominator,
                expected_average_numerator,
                expected_average_denominator,
                traces[: order + 2],
                len(target_profiles),
            )
            average_records.append(average_record)
            distinct_target_profiles[(order, len(target_profiles))] += 1
            total_tournaments += 1

        order_records.append((order, tournament_count, order_cases))

    order_records = tuple(order_records)
    require(order_records == EXPECTED_ORDER_COUNTS, order_records)
    require(total_tournaments == EXPECTED_TOURNAMENTS, total_tournaments)
    require(total_cases == EXPECTED_SWITCH_CASES, total_cases)
    case_sha = case_digest.hexdigest()
    average_records = tuple(average_records)
    average_sha = hashlib.sha256(repr(average_records).encode()).hexdigest()
    if EXPECTED_CASE_SHA256 is not None:
        require(case_sha == EXPECTED_CASE_SHA256, case_sha)
    if EXPECTED_AVERAGE_SHA256 is not None:
        require(average_sha == EXPECTED_AVERAGE_SHA256, average_sha)

    hostile = hostile_control()
    profile_census = tuple(sorted(distinct_target_profiles.items()))
    semantic_packet = (
        order_records,
        total_tournaments,
        total_cases,
        total_moment_cells,
        case_sha,
        average_sha,
        profile_census,
        hostile,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "Tournament switching centered-coronal walk compiler",
        (
            "operation=intrinsic_cut_switching_reverses_every_arc_between_"
            "the_two_sign_parts;centered_adjacency_C=2A-J_transforms_by_"
            "C_to_DCD;switch_sign_d_is_retained_up_to_global_negation"
        ),
        (
            "identity=P(z)=det(I-zC);U_d(z)=dT(I-zC)^(-1)d=N_d(z)/P(z);"
            "G_(T^d)(2z)=sum_k_2^k_walk_k(T^d)z^k="
            "U_d/(1-zU_d)=N_d/(P-zN_d)"
        ),
        (
            "characteristic=det(I-2zA_(T^d))=P-zN_d;the_complete_scaled_"
            "walk_sequence_has_a_constant_coefficient_recurrence_from_this_"
            "degree_at_most_n_denominator"
        ),
        (
            "coefficient_compiler=s_0=u_0=n;s_k=u_k+sum_(i=0)^(k-1)_"
            "u_i_s_(k-1-i),where_u_k=dT_C^k_d_and_s_k=2^k_walk_k(T^d)"
        ),
        (
            "switching_class_average=average_d_U_d=tr((I-zC)^(-1));"
            "average_d_N_d=nP-zPprime;average_d_det(I-2zA_(T^d))="
            "(1-nz)P+z^2Pprime"
        ),
        (
            f"exhaustion=orders_1_to_5;tournaments:{total_tournaments};"
            f"switch_cases:{total_cases};moment_cells:{total_moment_cells};"
            f"order_counts:{order_records};case_sha256:{case_sha};"
            f"average_sha256:{average_sha}"
        ),
        f"target_profile_census={profile_census}",
        (
            "hostile=transitive_T3_two_singleton_cuts_have_the_same_cut_size_"
            f"and_centered_denominator_P:{hostile[0]};source_N:{hostile[1]};"
            f"middle_N:{hostile[2]};source_adjacency_denominator:{hostile[3]};"
            f"middle_adjacency_denominator:{hostile[4]};source_walks:"
            f"{hostile[5]};middle_walks:{hostile[6]}"
        ),
        (
            "first_failed_implication=centered_spectrum_or_switching_class_"
            "plus_cut_cardinality_does_not_determine_the_walk_sequence;the_"
            "observer_numerator_N_d_is_load_bearing"
        ),
        (
            "source_target_map=source_is_the_centered_switching_class_with_"
            "a_signed_cut_observer_d;map_is_diagonal_conjugation_then_one_"
            "rank_one_J_update;target_is_the_total_directed_walk_sequence_"
            "and_adjacency_recurrence_of_the_selected_switched_tournament"
        ),
        (
            "preserved_and_lost=P_is_switching_invariant_and_N_d_restores_"
            "the_selected_total_walk_series;forgetting_d_or_N_d_loses_the_"
            "fixed_all_ones_observer;the_interface_does_not_restore_"
            "Hamiltonian_paths_path_cover_profiles_SCC_order_or_substitution_"
            "run_content"
        ),
        (
            "cheapest_decisive_test=the_three_vertex_transitive_tournament_"
            "switched_at_a_source_versus_at_the_middle_vertex"
        ),
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
