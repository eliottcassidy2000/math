#!/usr/bin/env python3
"""Primary exact companion for the THM-2806 tau-twelve addendum.

The universe is deliberately narrow: the three positive THM-2807 cylinders
on the fixed rail-eight sheet

    (physical clock, sigma, tau) = (1, 0, 12).

Before integration, write A for the source one-sided carrier, B for the
pulled-back target one-sided carrier, M=A intersect B, L=A\\B, and R=B\\A.
The script proves that M is already empty before any address cylinder is
selected, while each of the three selected cylinders is one whole positive
R-piece.  It also reconstructs the separate 81-label common bank so that its
clock support is not conflated with the tau-twelve right cofiber.

This closes only item 1 in THM-2807's proposed allocation test for this
specified simplex.  It constructs no endpoint origin, determinant current,
root-deck boundary, row exclusion, or LRC(14) proof.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

P = 13
Q = 7
N0 = 3_454_614
N_PLUS = N0 + P
N_A = N0 + 689_364
CENTRAL_INDICES = (0, 1, 53_028)
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
EXPECTED_COMMON_CLOCKS = (0, 81, 56, 56, 0, 0, 0)

PINNED = {
    "lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    "lrc14_positive_graded_address_two_simplex_thm2807.py":
        "11cdbe3c6cc7f9d5b6b24863ced71eb91cc84adc67fe38a3f8a3e637362453fb",
    "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected in PINNED.items():
    require(
        lf_hash(COMP / dependency) == expected,
        f"pinned dependency changed: {dependency}",
    )


import lrc14_semantic_arm_right_wing_central_digit_thm2782 as arm


def shift_weighted(pieces, amount):
    result = []
    for left, right, weight in pieces:
        length = right - left
        start = (left + amount) % arm.T
        stop = start + length
        if stop <= arm.T:
            result.append((start, stop, weight))
        else:
            result.extend(
                ((start, arm.T, weight), (0, stop - arm.T, weight))
            )
    return tuple(sorted(result))


def restricted_allocation(context, objects, key, central_index):
    """Apply exactly the THM-2807 target-cylinder filters to one face."""
    module, delayed, present, _source, _clocks, source_weight, _rail = context
    source_safe = arm.marked.complement(present[1, 7])
    target_safe = arm.marked.complement(
        arm.marked.shift_union(present[1, 7], arm.marked.SHIFT)
    )
    carrier = arm.marked.intersect(objects[key], source_safe)
    carrier = arm.marked.intersect(carrier, target_safe)
    carrier = arm.marked.shift_union(carrier, -arm.marked.SHIFT)
    weighted = arm.marked.restrict_weighted(source_weight, carrier)
    weighted = arm.marked.private.old.intersect_weighted_comb(
        weighted, module.C3, 182, 1, 13
    )
    weighted = arm.marked.restrict_weighted(
        weighted,
        arm.target_cylinder(arm.EXPECTED_BASES[1], central_index),
    )
    mass = arm.weighted_mass(weighted)
    coefficient = arm.marked.private.delayed_carry_pair(
        weighted, delayed[1], {}
    )[6][1]
    return tuple(weighted), mass, coefficient


def main():
    require(
        arm.EXPECTED_BASES[1] == N0
        and tuple(N0 + P * index for index in CENTRAL_INDICES)
        == (N0, N_PLUS, N_A),
        "THM-2807 selected addresses changed",
    )
    context = arm.build_context()
    module, _delayed, _present, source, clocks, _weight, rail_common = context

    objects12 = arm.carrier_objects(
        module, source, clocks, rail_common, 1, 0, 12
    )
    raw_piece_counts = tuple(
        len(objects12[key]) for key in ("A", "B", "M", "L", "R")
    )
    require(
        raw_piece_counts == (240, 241, 0, 240, 241)
        and not objects12["M"]
        and arm.marked.intersect(objects12["A"], objects12["B"])
        == objects12["M"]
        and objects12["L"] == objects12["A"]
        and objects12["R"] == objects12["B"],
        "tau-twelve physical allocation decomposition changed",
    )

    selected_rows = []
    right_pieces = []
    for address, central_index in zip(
        (N0, N_PLUS, N_A), CENTRAL_INDICES
    ):
        allocations = {
            key: restricted_allocation(
                context, objects12, key, central_index
            )
            for key in ("A", "B", "M", "L", "R")
        }
        require(
            all(
                allocations[key] == ((), Fraction(0), 0)
                for key in ("A", "M", "L")
            )
            and allocations["B"] == allocations["R"]
            and allocations["R"][1:]
            == (EXPECTED_MASS, EXPECTED_COEFFICIENT)
            and len(allocations["R"][0]) == 1
            and allocations["R"][0][0][2] == EXPECTED_WEIGHT,
            f"selected tau-twelve allocation face changed at {address}",
        )
        right_pieces.append(allocations["R"][0])
        selected_rows.append((
            address,
            central_index,
            tuple(
                (key, len(value[0]), value[1], value[2])
                for key, value in allocations.items()
            ),
        ))

    for index in (1, 2):
        address_gap = (N_PLUS, N_A)[index - 1] - N0
        coordinate_shift = address_gap * arm.EXPECTED_SHIFT
        require(
            shift_weighted(right_pieces[0], coordinate_shift)
            == right_pieces[index],
            "whole right-cofiber cylinder stopped translating",
        )

    # Positive and hostile controls: the tau-three diagonal survives on the
    # first and third vertices, while tau twelve fills all three.
    objects3 = arm.carrier_objects(
        module, source, clocks, rail_common, 1, 0, 3
    )
    require(
        objects3["M"] and not objects12["M"],
        "common-bank positive / tau-twelve hostile control changed",
    )
    table3 = tuple(
        restricted_allocation(context, objects3, "R", index)[1:]
        for index in CENTRAL_INDICES
    )
    table12 = tuple(
        restricted_allocation(context, objects12, "R", index)[1:]
        for index in CENTRAL_INDICES
    )
    positive = (EXPECTED_MASS, EXPECTED_COEFFICIENT)
    require(
        table3 == (positive, (Fraction(0), 0), positive)
        and table12 == (positive, positive, positive),
        "THM-2807 selected positive table changed",
    )

    # Reconstruct exactly the separate bank from THM-2806 Section 2.
    common_clock_counts = []
    right_clock_counts = []
    paired_support_checks = 0
    disjoint_checks = 0
    common_piece_histograms = []
    for clock in range(Q):
        common_count = 0
        right_count = 0
        piece_histogram = Counter()
        for sigma in COMMON_S:
            for tau in COMMON_T:
                objects = arm.carrier_objects(
                    module, source, clocks, rail_common,
                    clock, sigma, tau,
                )
                common_present = bool(objects["M"])
                right_present = bool(objects["R"])
                common_count += common_present
                right_count += right_present
                paired_support_checks += common_present == right_present
                piece_histogram[len(objects["M"])] += 1
                require(
                    not arm.marked.intersect(objects["M"], objects["R"]),
                    "common carrier met the right cofiber",
                )
                disjoint_checks += 1
        common_clock_counts.append(common_count)
        right_clock_counts.append(right_count)
        common_piece_histograms.append(tuple(sorted(piece_histogram.items())))

    common_clock_counts = tuple(common_clock_counts)
    right_clock_counts = tuple(right_clock_counts)
    require(
        common_clock_counts == EXPECTED_COMMON_CLOCKS
        and right_clock_counts == EXPECTED_COMMON_CLOCKS
        and paired_support_checks == 81 * Q
        and disjoint_checks == 81 * Q
        and 12 not in COMMON_T,
        "separate common-bank/right-cofiber census changed",
    )

    omega_digits = tuple(
        (address % P, address // P % P)
        for address in (N0, N_PLUS, N_A)
    )
    require(
        omega_digits == ((7, 6), (7, 7), (7, 7))
        and tuple(address % (P**2)
                  for address in (N0, N_PLUS, N_A))
        == (85, 98, 98),
        "two-value endpoint-decoder address projection changed",
    )

    source_text = Path(__file__).read_text()
    tree = ast.parse(source_text)
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "optimized-mode validity gate contains assert",
    )

    print("THM-2806 TAU-TWELVE SIMPLEX ALLOCATION-SUPPORT ADDENDUM")
    print("status=FINITE-EXACT primary; independent hostile audit required")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(
        "universe=rail8;(clock,sigma,tau)=(1,0,12);"
        f"addresses=({N0},{N_PLUS},{N_A});"
        "both_relative_safeties;target_root1;delayed_carry6"
    )
    print(
        "typing=A=source_one_sided;B=pullback(target_one_sided);"
        "M=A_intersect_B;L=A_minus_B;R=B_minus_A"
    )
    print(
        f"tau12_raw_piece_counts_(A,B,M,L,R)={raw_piece_counts};"
        "M_empty_before_address_restriction=yes;"
        "therefore_address_independent_on_this_fixed_sheet=yes"
    )
    for row in selected_rows:
        print(
            f"tau12_vertex=address:{row[0]};central_index:{row[1]};"
            f"allocation_rows:{row[2]}"
        )
    print(
        f"selected_right_cofiber_tables=tau3:{table3};tau12:{table12};"
        "whole_cylinder_translation=3/3"
    )
    print(
        f"separate_81_label_common_bank_clock_counts="
        f"{common_clock_counts};"
        f"right_cofiber_clock_counts={right_clock_counts};"
        f"cellwise_nonempty_support_agreement={paired_support_checks}/567"
    )
    print(
        f"common_piece_histograms={tuple(common_piece_histograms)};"
        f"common_bank_vs_right_cofiber_disjoint={disjoint_checks}/567;"
        "tau12_not_in_common_target_bank=yes"
    )
    print(
        f"Omega_digits_n0_nplus_na={omega_digits};"
        "projected_mod169=(85,98,98)"
    )
    print(
        "CLOSES=THM2807_item1_only_for_the_specified_tau12_simplex:"
        "the_requested_common_allocation_atom_is_empty"
    )
    print(
        "FIRST_FAILED_IMPLICATION="
        "positive_same_ancestry_address_simplex_does_not_imply_one_common_"
        "physical_allocation_atom"
    )
    print(
        "DECODER_BOUNDARY="
        "THM2803_two_value_ratio_requires_one_shared_physical_current_"
        "and_scalar;this_tau12_simplex_is_target_only_before_endpoint_"
        "origin_or_Fourier_data_exist"
    )
    print(
        "NEXT_TEST="
        "construct_a_typed_boundary_map_from_the_separate_common_bank_to_"
        "the_tau12_right_cofiber_or_find_a_positive_address_simplex_inside_"
        "the_common_bank_before_endpoint_allocation"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
