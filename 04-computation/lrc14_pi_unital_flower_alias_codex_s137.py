#!/usr/bin/env python3
"""
Numerical guardrail for the pi/unital/flower prompt.

This is not a proof script.  It fixes the normalization facts that should be
kept attached before importing the numerology into the LRC14 carrier stack.
"""

import math
from fractions import Fraction


def signed_error(x: float, target: float = math.pi) -> float:
    return x - target


def best_closure(theta: float, qmax: int):
    """Return best q<=qmax for q*theta near an integer number of turns."""
    best = None
    for q in range(1, qmax + 1):
        turns = q * theta / (2 * math.pi)
        k = round(turns)
        residual = q * theta - k * 2 * math.pi
        item = (abs(residual), q, k, residual)
        if best is None or item < best:
            best = item
    abs_residual, q, k, residual = best
    return q, k, residual, abs_residual


def print_pi_approximants():
    print("PI APPROXIMANTS")
    approximants = [
        ("22/7", 22 / 7, "rational/Farey approximant"),
        ("cuberoot(31)", 31 ** (1 / 3), "cubic algebraic approximant"),
    ]
    for name, value, role in approximants:
        err = signed_error(value)
        print(f"{name}: {value:.15f}")
        print(f"  role: {role}")
        print(f"  signed_error_vs_pi: {err:+.15e}")
        print(f"  abs_error_vs_pi: {abs(err):.15e}")
        print(f"  relative_error_ppm: {abs(err) / math.pi * 1_000_000:.6f}")
    ratio = abs(signed_error(22 / 7)) / abs(signed_error(31 ** (1 / 3)))
    print(f"cuberoot(31)_is_better_than_22/7_by_abs_error_factor: {ratio:.6f}")
    print(f"pi_cubed_minus_31: {math.pi ** 3 - 31:+.15e}")
    print()


def print_flower_normalization():
    print("FLOWER ROTATION NORMALIZATION")

    theta_literal = 1 / math.pi
    print("literal_radian_step: theta = 1/pi radians")
    print(f"  theta_rad: {theta_literal:.15f}")
    print(f"  turn_fraction: {theta_literal / (2 * math.pi):.15f}")
    print(f"  local_steps_per_turn: {2 * math.pi / theta_literal:.15f}")
    q, k, residual, abs_residual = best_closure(theta_literal, 100)
    print(f"  best_closure_q<=100: q={q}, turns={k}, residual_rad={residual:+.15e}")
    q22_residual = 22 * theta_literal - round(22 * theta_literal / (2 * math.pi)) * 2 * math.pi
    print(f"  q=22_residual_rad: {q22_residual:+.15e}")
    print("  readout: literal 1/pi radians does not make a 22-period orbit.")

    theta_turn = 2.0
    print("fractional_turn_step: step = 1/pi of a full turn = 2 radians")
    print(f"  theta_rad: {theta_turn:.15f}")
    print(f"  exact_turn_fraction: {1 / math.pi:.15f}")
    print(f"  22/7_pi_gives_turn_fraction: {Fraction(7, 22)}")
    residual22 = 22 * theta_turn - 7 * 2 * math.pi
    print(f"  q=22_vs_7_turns_residual_rad: {residual22:+.15e}")
    print(f"  q=22_vs_7_turns_residual_turns: {residual22 / (2 * math.pi):+.15e}")
    q, k, residual, abs_residual = best_closure(theta_turn, 100)
    print(f"  best_closure_q<=100: q={q}, turns={k}, residual_rad={residual:+.15e}")
    print("  readout: the 22-family flower belongs to the full-turn normalization,")
    print("           because 1/pi turns is approximated by 7/22 turns.")
    print()


def print_unital_guardrail():
    print("UNITAL TERMINOLOGY GUARDRAIL")
    q = 3
    v = q ** 3 + 1
    block_size = q + 1
    lam = 1
    blocks = v * (v - 1) // (block_size * (block_size - 1))
    replication = (v - 1) // (block_size - 1)
    print("block_design_unital_q=3:")
    print(f"  points: {v}")
    print(f"  block_size: {block_size}")
    print(f"  lambda: {lam}")
    print(f"  blocks: {blocks}")
    print(f"  point_replication: {replication}")
    print("  LRC use: secondary pair-frame candidate for C(8,2)=28 slots.")
    print("algebraic_unital:")
    print("  meaning: an identity/unit is present, or a map preserves identity.")
    print("  LRC use: metaphor for label-preserving quotient maps, not a design.")
    print("residue_unit_shell:")
    print("  meaning: gcd-class visibility such as g1 versus g3/g9 in C=27 shells.")
    print("  LRC use: primary finite shell label in HYP-2937, not the word unital.")
    print()


def print_carrier_order():
    print("PROOF-CARRIER ORDER")
    print("angle normalization")
    print("> exact M/Farey branch")
    print("> C=27 unit/nonunit shell transfer")
    print("> 22/7 and cuberoot(31) pi-approximant labels")
    print("> q=3 unital pair-frame after category-1 labelling")
    print("> algebraic unital identity-preservation metaphor")
    print("> raw flower/family visual count")
    print()
    print("ASSUMPTION CHALLENGED")
    print("A visible 22-family flower should not be treated as an exact")
    print("literal-radian invariant.  It is an aliasing witness for the")
    print("1/pi-of-a-turn normalization and for the rational approximant 7/22.")


def main():
    print_pi_approximants()
    print_flower_normalization()
    print_unital_guardrail()
    print_carrier_order()


if __name__ == "__main__":
    main()
