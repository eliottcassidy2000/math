#!/usr/bin/env python3
"""Exact parity atlas and minimal physical-row hostile for LRC word currents.

This standard-library companion checks two symbolic mechanisms.

First, at a common boundary of a source comb bank and its K-clock pulled
copy, safe/danger bits obey an XOR alignment law.  The law belongs to the
combined collision packet; separately assigned Leibniz terms can cancel.

Second, an exact u=2 interval model realizes complementary full-set and
word-current spectra in both physical rows.  It preserves a literal
T^2 pullback and the arithmetic row data, but is deliberately not claimed
to be a full canonical THM-2305 Boolean word or global LRC cover.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bit_value(bit: int, orientation: int, side: int) -> int:
    """Safe=0, danger=1; side and orientation are each +/-1."""

    return int(((-1) ** bit) * orientation * side > 0)


def main() -> None:
    # The XOR germ and the repaired combined-Leibniz provenance.
    for clock_parity in (0, 1):
        require(
            pow(13, clock_parity, 14)
            == (1 if clock_parity == 0 else 13),
            "13-clock orientation parity changed",
        )
        for source_bit in (0, 1):
            for target_bit in (0, 1):
                aligned = (
                    target_bit - clock_parity - source_bit
                ) % 2 == 0
                for epsilon in (-1, 1):
                    source_germ = [
                        bit_value(source_bit, epsilon, side)
                        for side in (-1, 1)
                    ]
                    target_germ = [
                        bit_value(
                            target_bit,
                            ((-1) ** clock_parity) * epsilon,
                            side,
                        )
                        for side in (-1, 1)
                    ]
                    product_germ = [
                        source_germ[index] * target_germ[index]
                        for index in range(2)
                    ]
                    require(
                        product_germ
                        == (source_germ if aligned else [0, 0]),
                        "safe/danger XOR germ failed",
                    )

                    # Specialize to a safe source S.  Only the sum of the
                    # right-sided Leibniz pieces has the XOR law.
                    if source_bit == 0:
                        sigma = (-1) ** (clock_parity + target_bit)
                        source_plus = Fraction(1 + epsilon, 2)
                        source_minus = Fraction(1 - epsilon, 2)
                        target_plus = Fraction(1 + sigma * epsilon, 2)
                        delta_source = Fraction(epsilon)
                        delta_target = Fraction(sigma * epsilon)
                        source_piece = target_plus * delta_source
                        target_piece = source_minus * delta_target
                        require(
                            source_piece + target_piece
                            == Fraction(epsilon * (1 + sigma), 2),
                            "combined collision provenance changed",
                        )

    # Exact arithmetic of the minimal u=2 hostile.
    lam = 1
    clock = lam + 1
    b = 3
    shell = 2 * 13**5
    multiplier = 170
    a_prime = 13
    a_value = a_prime + multiplier * shell
    step = 13**b
    require(clock == 2, "clock changed")
    require(shell == 742586, "shell changed")
    require(a_value == 126239633, "second row base changed")
    require(gcd(multiplier, 91) == 1, "multiplier lost its 7/13 unit")

    frequencies = {
        0: [a_prime, a_value],
        1: [a_prime + step, a_value + step],
    }
    require(
        frequencies[1] == [2210, 126241830],
        "shifted physical frequencies changed",
    )
    require(
        frequencies[1][1] == 57123 * frequencies[1][0],
        "shifted rows lost the common interval period",
    )

    i_length = Fraction(3, 4420)
    j_length = Fraction(1, 2210)
    require(j_length < i_length < Fraction(1, 169), "bad interval nesting")

    # For E=I union (I+1/2), the half-shift multiplier kills exactly odd
    # frequencies.  A centered interval coefficient vanishes exactly when
    # its nonzero frequency times its length is an integer.
    for address in (0, 1):
        for frequency in frequencies[address]:
            e_zero = frequency % 2 == 1 or (
                frequency * i_length
            ).denominator == 1
            j_zero = (frequency * j_length).denominator == 1
            if address == 0:
                require(e_zero and not j_zero, "k=0 zero pattern changed")
            else:
                require(not e_zero and j_zero, "k=1 zero pattern changed")
                require(
                    (frequency * i_length).denominator == 2,
                    "the surviving I phase is no longer half-integral",
                )

    # In K[X]/(X^2-1), the evaluation patterns are
    # F proportional to 1-X and W proportional to 1+X.
    # Their product is 1-X^2 and hence zero in the quotient.
    full_current = [1, -1]
    word_current = [1, 1]
    product = [
        full_current[0] * word_current[0],
        full_current[0] * word_current[1]
        + full_current[1] * word_current[0],
        full_current[1] * word_current[1],
    ]
    require(product == [1, 0, -1], "u=2 current product changed")
    require(
        [product[0] + product[2], product[1]] == [0, 0],
        "u=2 product is not zero modulo X^2-1",
    )

    # Literal prescribed-clock pullback geometry.  If J is centered at tau,
    # the other T^2-preimages are shifted by r/169.  They miss I and I+1/2.
    same_component_clearance = Fraction(1, 169)
    half_component_clearance = Fraction(1, 338)
    combined_half_width = (i_length + j_length) / 2
    require(
        combined_half_width < half_component_clearance
        < same_component_clearance,
        "pullback components can intersect the wrong interval",
    )
    q_length = 169 * j_length
    require(q_length == Fraction(13, 170), "target word length changed")
    require(q_length < Fraction(1, 7), "target word left one danger arc")

    common_boundaries = 2 * shell
    target_only_boundaries = 2 * shell * (13**clock - 1)
    require(
        common_boundaries + target_only_boundaries
        == 2 * shell * 13**clock,
        "source/target boundary count no longer partitions",
    )

    print("clock_parity_law=delta_equals_K_plus_alpha_mod_2")
    print("parity_scope=combined_source_pulled_target_collision_packet")
    print("separate_Leibniz_provenance=can_cancel_exactly")
    print("canonical_even_fork_source_support=guard_plus_five_units_only")
    print("canonical_odd_fork_source_support=all_three_blockers")
    print("whole_row_incidence=some_FjWj_nonzero_mod_Xu_minus_1")
    print("whole_row_gcd_test=lcm(gcd(Fj,N),gcd(Wj,N))_not_equal_N")
    print("hostile_u=2")
    print("hostile_frequencies_k0=13,126239633")
    print("hostile_frequencies_k1=2210,126241830")
    print("hostile_full_current=1-X")
    print("hostile_word_current=1+X")
    print("hostile_product=0_mod_X2_minus_1")
    print("literal_pullback=E_intersect_Tminus2_Q_equals_J")
    print("target_word_length=13/170_less_than_1/7")
    print(f"common_boundaries={common_boundaries}")
    print(f"target_only_boundaries={target_only_boundaries}")
    print("scope=local_word_control_not_full_THM2305_cover")
    print("status=LRC_WORD_PARITY_AND_U2_HOSTILE_EXACT_AUDIT")


if __name__ == "__main__":
    main()
