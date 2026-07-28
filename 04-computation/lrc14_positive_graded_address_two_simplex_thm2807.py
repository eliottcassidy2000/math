#!/usr/bin/env python3
"""Exact companion for THM-2807's positive graded-address two-simplex.

The script rebuilds the inherited THM-2782 carrier once, restricts the
tau=3 and tau=12 profiles to the three claimed cylinders, and checks the
lower-central and affine-address arithmetic.  It deliberately does not
claim factorwise covariance or endpoint allocation.
"""

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_semantic_arm_right_wing_central_digit_thm2782 as arm
import lrc14_full_arm_orbit_lower_central_chord_thm2791 as chord


P = 13
ADDRESS_MODULUS = P**6
N0 = 3_454_614
N_PLUS = N0 + 13
N_A = N0 + 689_364
J0 = 0
J_PLUS = 1
J_A = 53_028
EXPECTED_WEIGHT = 27_581_135_604
EXPECTED_MASS = Fraction(60_781_651_775_958_960, 371_293)
EXPECTED_COEFFICIENT = 790_161_473_087_466_480


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cylinder(base, central_index):
    center = arm.arm_target_center(base, central_index)
    half = arm.relative.Q_RADIUS * arm.T
    coordinate = center * arm.T
    require(0 < coordinate - half < coordinate + half < arm.T,
            "selected cylinder crossed the circle seam")
    return ((coordinate - half, coordinate + half),)


def restricted_piece(weighted, base, central_index):
    selected = cylinder(base, central_index)
    pieces = arm.marked.restrict_weighted(weighted, selected)
    require(
        len(pieces) == 1 and pieces[0][:2] == selected[0],
        "selected positive cylinder stopped being one whole weighted piece",
    )
    return pieces[0]


def translate_piece(piece, shift):
    left = (piece[0] + shift * arm.T) % arm.T
    right = (piece[1] + shift * arm.T) % arm.T
    require(left < right, "translated cylinder wrapped or reversed")
    return left, right, piece[2]


def typed_metadata(base, central_index):
    center = arm.arm_target_center(base, central_index)
    sigma_labels, tau_labels = arm.relative.full_target_labels(center)
    return (
        arm.relative.semantic.semantic_record(center),
        arm.predecessor_carry(center),
        sigma_labels,
        tau_labels,
    )


def main():
    print("THM-2807 POSITIVE GRADED-ADDRESS TWO-SIMPLEX")
    print("status=VERIFIED-EXACT candidate; no endpoint allocation or LRC14")

    require(
        arm.EXPECTED_BASES[1] == N0
        and N_PLUS == N0 + P
        and N_A == N0 + P * J_A,
        "selected absolute addresses changed",
    )
    context = arm.build_context()
    base = arm.EXPECTED_BASES[1]
    delayed3, weighted3 = chord.target_weighted_by_tau(context, 3)
    delayed12, weighted12 = chord.target_weighted_by_tau(context, 12)

    cells3 = tuple(
        chord.exact_cell(delayed3, weighted3, base, index)
        for index in (J0, J_PLUS, J_A)
    )
    cells12 = tuple(
        chord.exact_cell(delayed12, weighted12, base, index)
        for index in (J0, J_PLUS, J_A)
    )
    positive_cell = (EXPECTED_MASS, EXPECTED_COEFFICIENT)
    require(
        cells3 == (positive_cell, (0, 0), positive_cell),
        "tau=3 diagonal/deleted-vertex profile changed",
    )
    require(
        cells12 == (positive_cell, positive_cell, positive_cell),
        "tau=12 three-vertex filler profile changed",
    )

    pieces12 = tuple(
        restricted_piece(weighted12, base, index)
        for index in (J0, J_PLUS, J_A)
    )
    require(
        tuple(piece[2] for piece in pieces12)
        == (EXPECTED_WEIGHT,) * 3,
        "tau=12 cylinder weights stopped agreeing",
    )
    pieces3 = (
        restricted_piece(weighted3, base, J0),
        restricted_piece(weighted3, base, J_A),
    )
    require(
        pieces3 == (pieces12[0], pieces12[2])
        and arm.marked.restrict_weighted(
            weighted3, cylinder(base, J_PLUS)
        ) == (),
        "tau=3 diagonal is no longer the tau=12 face with one vertex deleted",
    )

    metadata = tuple(
        typed_metadata(base, index) for index in (J0, J_PLUS, J_A)
    )
    require(
        all(
            semantic == (3, (1, 2))
            and carry == 6
            and 0 in sigma_labels
            and 12 in tau_labels
            for semantic, carry, sigma_labels, tau_labels in metadata
        )
        and 3 in metadata[0][3]
        and 3 in metadata[2][3],
        "semantic/carry/target typing changed",
    )

    pure_gap = N_PLUS - N0
    vertical_gap = N_A - N_PLUS
    diagonal_gap = N_A - N0
    require(
        pure_gap == P
        and vertical_gap == P**2 * 4_079
        and diagonal_gap == P * 53_028
        and 53_028 == 1 + P * 4_079
        and 4_079 % P == 10,
        "lower-central edge arithmetic changed",
    )
    require(
        N0 % (P**2) == 85
        and N_PLUS % (P**2) == N_A % (P**2) == 98,
        "endpoint-quotient horn changed",
    )

    shift_pure = Fraction(7 * pure_gap, ADDRESS_MODULUS)
    shift_vertical = Fraction(7 * vertical_gap, ADDRESS_MODULUS)
    shift_diagonal = Fraction(7 * diagonal_gap, ADDRESS_MODULUS)
    require(
        (shift_pure, shift_vertical, shift_diagonal)
        == (
            Fraction(7, 371_293),
            Fraction(28_553, 28_561),
            Fraction(371_196, 371_293),
        )
        and (shift_pure + shift_vertical) % 1 == shift_diagonal,
        "physical translation fractions changed",
    )
    require(
        translate_piece(pieces12[0], shift_pure) == pieces12[1]
        and translate_piece(pieces12[1], shift_vertical) == pieces12[2]
        and translate_piece(pieces12[0], shift_diagonal) == pieces12[2],
        "whole-cylinder triangle stopped commuting",
    )

    depth5 = P**5
    logarithms = tuple(
        exponent
        for exponent in range(P**4)
        if pow(14, exponent, depth5) == 53_028
    )
    require(
        logarithms == (23_098,),
        "14-adic logarithm stopped being unique modulo 13^4",
    )
    exponent = logarithms[0]
    multiplier = pow(14, exponent, ADDRESS_MODULUS)
    translation = (N0 - multiplier * N0) % ADDRESS_MODULUS
    require(
        multiplier == 2_652_079
        and translation == 352_469
        and (multiplier * N0 + translation) % ADDRESS_MODULUS == N0
        and (multiplier * N_PLUS + translation) % ADDRESS_MODULUS == N_A
        and multiplier * pure_gap % ADDRESS_MODULUS == diagonal_gap,
        "bare affine address conjugacy changed",
    )

    print(
        f"addresses=({N0},{N_PLUS},{N_A}); "
        f"quotient_mod169=({N0 % 169},{N_PLUS % 169},{N_A % 169})"
    )
    print(
        "tau3_cells="
        f"{cells3}; tau12_cells={cells12}"
    )
    print(
        "whole_cylinder_weight="
        f"{EXPECTED_WEIGHT}; tau3_equals_tau12_on_diagonal=yes"
    )
    print(
        "edge_gaps="
        f"({pure_gap},{vertical_gap},{diagonal_gap})="
        "(Z1,Z2^4079,Z1*Z2^4079); transgression=4079_mod13=10"
    )
    print(
        "physical_shifts="
        f"({shift_pure},{shift_vertical},{shift_diagonal}); "
        "triangle_commutes=yes"
    )
    print(
        "affine_address_conjugacy="
        f"k={exponent}; 14^k_mod13^6={multiplier}; "
        f"g(n)=14^k*n+{translation}; g(n0)=n0; g(nplus)=na"
    )
    print(
        "scope=positive collapsed-carrier/address nerve simplex only; "
        "no factor covariance, endpoint origin, K4 allocation, Cech map, "
        "row exclusion, or LRC14"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
