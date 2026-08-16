#!/usr/bin/env python3
"""Exact typing audit for the r=5 tent-location / hypothetical C4 cospan.

This script deliberately separates four coefficient sorts:

* the six-coordinate exception-location section;
* the arc-reversal-odd coefficient register;
* a divergence-free chain current on a formally closed C4; and
* an edge cochain representing graph H^1.

The same tuples can occur in these spaces, but no response-amplitude or
physical-current transport between them is assumed.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import lrc_r5_pointed_bidirected_p4_cocycle_sidecar_20260816 as P4SIDE


H = (12, 12, 9, 3, 0, 0)
ARC_REVERSAL = (1, 0, 3, 2, 5, 4)
J_MID = (0, 0, 1, -1, 0, 0)
P13 = 13
RESPONSE_PRIME = 755373809845391722745761

EXPECTED_SEMANTIC_SHA256 = (
    "d662e9ffb12e35f60b34fcf86d5cc60788aded7de28fdcc132cfe1fc693290c3"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank_mod(rows, prime: int) -> int:
    matrix = [[value % prime for value in row] for row in rows]
    if not matrix:
        return 0
    columns = len(matrix[0])
    require(all(len(row) == columns for row in matrix), "ragged matrix")
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, prime)
        matrix[pivot_row] = [value * inverse % prime for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % prime
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def mat_vec(matrix, vector):
    return tuple(sum(left * right for left, right in zip(row, vector))
                 for row in matrix)


def in_column_space(matrix, vector, prime: int) -> bool:
    base_rank = rank_mod(matrix, prime)
    augmented = tuple(tuple(row) + (value,)
                      for row, value in zip(matrix, vector))
    return rank_mod(augmented, prime) == base_rank


def main() -> None:
    reversed_h = tuple(H[index] for index in ARC_REVERSAL)
    sums = tuple(left + right for left, right in zip(H, reversed_h))
    differences = tuple(left - right for left, right in zip(H, reversed_h))
    require(all(value % 2 == 0 for value in sums + differences),
            "integral half projection")
    h_even = tuple(value // 2 for value in sums)
    h_odd = tuple(value // 2 for value in differences)
    require(h_even == (12, 12, 6, 6, 0, 0), h_even)
    require(h_odd == tuple(3 * value for value in J_MID), h_odd)
    require(tuple(left + right for left, right in zip(h_even, h_odd)) == H,
            "decomposition")
    require(tuple(h_even[index] for index in ARC_REVERSAL) == h_even,
            "even parity")
    require(tuple(h_odd[index] for index in ARC_REVERSAL)
            == tuple(-value for value in h_odd), "odd parity")

    # Vertices are A,B,C,D.  Path edges are AB, BC, CD; the formal closure is
    # DA, all oriented around A->B->C->D->A.
    boundary_p4 = (
        (-1, 0, 0),
        (1, -1, 0),
        (0, 1, -1),
        (0, 0, 1),
    )
    boundary_c4 = tuple(row + (closure,)
                        for row, closure in zip(
                            boundary_p4, (1, 0, 0, -1)))
    require(rank_mod(boundary_p4, P13) == 3, "P4 boundary rank")
    require(rank_mod(boundary_c4, P13) == 3, "C4 boundary rank F13")
    require(rank_mod(boundary_c4, 2) == 3, "C4 boundary rank F2")

    # In the odd-basis coordinate convention, h_odd has middle coefficient 3.
    # Under literal directed-arc chain realization, e_BC-e_CB maps to 2[BC].
    odd_basis_coefficients = (0, 3, 0)
    normalized_path_chain = odd_basis_coefficients
    literal_arc_path_chain = tuple(2 * value for value in odd_basis_coefficients)
    require(mat_vec(boundary_p4, normalized_path_chain) == (0, -3, 3, 0),
            "normalized path boundary")
    require(mat_vec(boundary_p4, literal_arc_path_chain) == (0, -6, 6, 0),
            "literal arc boundary")

    normalized_cycle = (3, 3, 3, 3)
    literal_cycle = (6, 6, 6, 6)
    require(mat_vec(boundary_c4, normalized_cycle) == (0, 0, 0, 0),
            "normalized cycle")
    require(mat_vec(boundary_c4, literal_cycle) == (0, 0, 0, 0),
            "literal cycle")

    # Divergence zero on C4 forces all four coefficients equal, so evaluation
    # on BC is an isomorphism from the cycle line to the coefficient line.
    for coefficient, cycle in ((3, normalized_cycle), (6, literal_cycle)):
        require(cycle[1] == coefficient, "middle evaluation")
        require(all(value == coefficient for value in cycle), "unique extension")

    restriction_normalized = normalized_cycle[:3]
    restriction_literal = literal_cycle[:3]
    require(restriction_normalized != normalized_path_chain,
            "whole normalized section is not extended")
    require(restriction_literal != literal_arc_path_chain,
            "whole literal section is not extended")

    # Cochains: delta is the transpose of the boundary matrix.
    coboundary_c4 = tuple(
        tuple(boundary_c4[vertex][edge] for vertex in range(4))
        for edge in range(4)
    )
    f2_alternating_potential = (0, 1, 0, 1)
    f2_all_one = (1, 1, 1, 1)
    require(tuple(value % 2 for value in
                  mat_vec(coboundary_c4, f2_alternating_potential))
            == f2_all_one, "F2 all-one cochain is exact")
    require(in_column_space(coboundary_c4, f2_all_one, 2),
            "F2 cochain containment")
    require(not in_column_space(coboundary_c4, normalized_cycle, P13),
            "F13 constant-three cochain is nonexact")

    seam_integer = sum(normalized_cycle)
    seam_f13 = seam_integer % P13
    seam_f2 = seam_integer % 2
    literal_seam_f13 = sum(literal_cycle) % P13
    require((seam_integer, seam_f13, seam_f2, literal_seam_f13)
            == (12, 12, 0, 11), "seam ledger")

    # As chains, both constant vectors remain nonzero H1 generators modulo
    # their coefficient fields; there are no graph 2-cells whose boundary
    # could make them exact.  The F2 exactness above is cohomological only.
    require(tuple(value % 2 for value in normalized_cycle) == f2_all_one,
            "integral cycle base change")
    require(mat_vec(boundary_c4, f2_all_one) == (0, 0, 0, 0),
            "F2 chain cycle")

    # There is no unital ring map F13 -> F2: 13*1=0 in F13, while its image
    # would be 13*1=1 in F2.  The lawful comparison is the coefficient cospan
    # Z -> F13 and Z -> F2 applied to the integral cochain of seam 12.
    characteristic_hostile = (
        "no_unital_ring_map_F13_to_F2",
        "no_nonzero_additive_map_F13_to_split_response_field",
        "integer_12_equals_minus_1_only_after_F13_reduction_not_in_response_field",
        "lawful_common_integral_lift_Z_to_F13_and_Z_to_F2",
        "arc_reversal_plus_minus_split_unavailable_in_characteristic_2",
    )
    require(RESPONSE_PRIME != P13 and RESPONSE_PRIME > 13,
            "response characteristic hostile")

    formal_cospan = (
        "Loc_exception_threshold",
        "middle_odd_coefficient_3_in_Z_address",
        "MISSING_tau_address_to_response_current",
        "ev_BC:Z1(C4;Z_current)_to_Z_current_is_iso",
        "formal_normalized_completion_3*(1,1,1,1)",
    )
    missing_gates = (
        "lawful_same-copy_D_to_A_closure_edge",
        "address-location_to-response-amplitude_coefficient_transport",
        "marked_chain_to_cochain_constitutive_identification",
        "semantic_C4_to_THM3496_C7_chart_map",
    )
    preserved = (
        "arc_reversal_middle_support",
        "middle_coefficient_after_explicit_normalization",
        "C4_divergence_zero_after_formal_completion",
        "integral_seam_and_its_separate_base_changes",
    )
    lost = (
        "exception_set_semantics",
        "even_location_baseline",
        "outer_edge_coefficients",
        "source_cells_and_endpoint_response",
        "clock_temporal_copy_and_physical_current",
    )

    # Pay the cheapest actual-response hostile.  The accepted diagonal profile
    # has an odd middle-arc response on exactly the source-middle address
    # support.  If the location coefficient 3 transported without retaining
    # r1, this response would be constant on that support.  Rebuild the pinned
    # parent and test the ten values at r0=3 and r0=9 directly.
    _cells, tensor, _source_digest, _tensor_digest = P4SIDE.build_actual_tensor()
    profile, _kernels, profile_digest = P4SIDE.extract_kernels(tensor)
    require(profile_digest == P4SIDE.EXPECTED_PROFILE_SHA256,
            "accepted profile drift")
    nonroot_digits = tuple(r1 for r1 in range(13) if r1 not in (0, 6, 12))
    middle_response_rows = tuple(
        tuple((profile[r0][r1][2] - profile[r0][r1][3]) % RESPONSE_PRIME
              for r1 in nonroot_digits)
        for r0 in (3, 9)
    )
    require(all(value for row in middle_response_rows for value in row),
            "middle response support")
    middle_distinct_counts = tuple(len(set(row)) for row in middle_response_rows)
    require(middle_distinct_counts == (8, 8), middle_distinct_counts)
    require(P4SIDE.M.rank_rows(middle_response_rows) == 2,
            "middle response row rank")
    middle_reflection = all(
        middle_response_rows[1][nonroot_digits.index(r1)]
        == -middle_response_rows[0][nonroot_digits.index(12 - r1)]
           % RESPONSE_PRIME
        for r1 in nonroot_digits
    )
    require(middle_reflection, "middle response reflection")
    response_transport_hostile = (
        profile_digest,
        nonroot_digits,
        middle_distinct_counts,
        P4SIDE.M.rank_rows(middle_response_rows),
        middle_reflection,
        digest_json(middle_response_rows),
        tuple(row[0] for row in middle_response_rows),
        "no_r1_blind_scalar_transport",
        "r1_dependent_section_still_requires_a_typed_map",
    )

    record = (
        H,
        ARC_REVERSAL,
        RESPONSE_PRIME,
        h_even,
        h_odd,
        odd_basis_coefficients,
        normalized_path_chain,
        literal_arc_path_chain,
        normalized_cycle,
        literal_cycle,
        restriction_normalized,
        restriction_literal,
        (seam_integer, seam_f13, seam_f2, literal_seam_f13),
        f2_alternating_potential,
        characteristic_hostile,
        formal_cospan,
        missing_gates,
        preserved,
        lost,
        response_transport_hostile,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 tent-location / hypothetical C4 cospan hostile audit ==")
    print(f"tent_(h,reversal,even,odd)={(H, reversed_h, h_even, h_odd)}")
    print(f"odd_identity=h_odd=3*j_mid; j_mid={J_MID}")
    print(f"odd_basis_P4_coefficients={odd_basis_coefficients}")
    print("normalization_record=(normalized_half_difference,literal_directed_arc)")
    print(f"P4_chains_and_boundaries={((normalized_path_chain, mat_vec(boundary_p4, normalized_path_chain)), (literal_arc_path_chain, mat_vec(boundary_p4, literal_arc_path_chain)))}")
    print(f"formal_C4_cycles={((normalized_cycle, literal_cycle), (rank_mod(boundary_c4, P13), rank_mod(boundary_c4, 2)))}")
    print(f"restriction_hostile_(cycle_restriction,observed_path)={((restriction_normalized, normalized_path_chain), (restriction_literal, literal_arc_path_chain))}")
    print(f"cochain_seams_(Z,F13,F2,literal_F13)={(seam_integer, seam_f13, seam_f2, literal_seam_f13)}")
    print(f"F2_exact_cochain_potential={f2_alternating_potential};delta_f={f2_all_one}")
    print("F2_chain_status=all-one is the nonzero H1 cycle; exactness above is H^1/cochain only")
    print(f"characteristic_hostile={characteristic_hostile}")
    print(f"formal_shadow_cospan={formal_cospan}")
    print(f"missing_typed_gates={missing_gates}")
    print(f"preserved={preserved}")
    print(f"lost={lost}")
    print("response_transport_record=(profile_sha,nonroot_r1,distinct_counts,row_rank,reflection,row_digest,first_values,scalar_verdict,remaining_gate)")
    print(f"response_transport_hostile={response_transport_hostile}")
    print(f"semantic_sha256={semantic}")
    print("verdict=useful marked candidate interface, blocked first at the closure edge and then at coefficient transport; no response amplitude, current, or D5 class")
    print("PASS")


if __name__ == "__main__":
    main()
