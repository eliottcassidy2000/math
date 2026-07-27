#!/usr/bin/env python3
"""Exact hostile for fixed-colour overlap versus a prescribed LRC row.

This strengthens the u=2 hostile from the parity-current packet.  The
shifted intervals satisfy the abstract hypotheses of THM-2323, have an
XOR-aligned common endpoint, and in fact have common Fourier support in
every four consecutive gauges of every primitive colour modulo 91 and
the actual normalized row modulus 169.
Nevertheless the four canonical diagonal samples retain complementary
spectra, and the total current product vanishes in K[C_2].
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def f_zero(frequency: int) -> bool:
    """P_1 1_J is the centered interval of length 1/170."""

    return frequency % 170 == 0


def g_zero(frequency: int) -> bool:
    """P_1 1_E is one interval plus its half translate."""

    return frequency % 2 != 0 or frequency % 340 == 0


def overlap(frequency: int) -> bool:
    return not f_zero(frequency) and not g_zero(frequency)


def main() -> None:
    # Original circle intervals.
    j_left = Fraction(-1, 4420)
    j_right = Fraction(1, 4420)
    i_left = j_left
    i_right = Fraction(1, 2210)
    require(j_left == i_left, "the XOR-aligned left endpoint moved")
    require(j_right < i_right, "the right endpoint is no longer target-only")
    require(
        j_right - j_left == Fraction(1, 2210),
        "J length changed",
    )
    require(
        i_right - i_left == Fraction(3, 4420),
        "I length changed",
    )

    # T^2 is injective on J and its image lies in the centered D_1 arc.
    q_left = 169 * j_left
    q_right = 169 * j_right
    require(
        (q_left, q_right) == (Fraction(-13, 340), Fraction(13, 340)),
        "T^2 J endpoints changed",
    )
    require(q_right < Fraction(1, 14), "Q left D_1")

    # Every other T^2 preimage is displaced by a nonzero r/169.  Since I
    # is not concentric with J, retain its center shift explicitly.  The
    # sum of the two half-widths plus that shift is 3/4420.  It is smaller
    # than both the same-component and half-shift spacings, proving
    # E intersect T^-2 Q = J.
    total_clearance_cost = Fraction(3, 4420)
    require(
        Fraction(1, 169) - total_clearance_cost > 0
        and Fraction(1, 338) - total_clearance_cost > 0,
        "pullback clearance failed",
    )

    # After one Perron step:
    # f=(1/13)1_[-1/340,1/340],
    # g=(1/13)(1_[-1/340,1/170] plus its half translate).
    f_left, f_right = 13 * j_left, 13 * j_right
    g_left, g_right = 13 * i_left, 13 * i_right
    require(
        (f_left, f_right)
        == (Fraction(-1, 340), Fraction(1, 340)),
        "P1 J changed",
    )
    require(
        (g_left, g_right)
        == (Fraction(-1, 340), Fraction(1, 170)),
        "P1 I changed",
    )
    require(
        max(abs(g_left), abs(g_right)) < Fraction(1, 28),
        "support left D_2",
    )
    require(gcd(2, 91) == 1, "D_2 lost its primitive comb coordinate")

    # Stronger-than-THM-2323 thickness: at both the 91-colour comparison
    # modulus and the actual normalized row modulus 169, every primitive
    # colour and every four-gauge block has a common coefficient.  The
    # property is periodic modulo 340.
    primitive_counts = {}
    for modulus in (91, 169):
        primitive_colours = [
            kappa for kappa in range(1, modulus)
            if gcd(kappa, modulus) == 1
        ]
        primitive_counts[modulus] = len(primitive_colours)
        for kappa in primitive_colours:
            for start in range(340):
                block = [
                    kappa + modulus * h
                    for h in range(start, start + 4)
                ]
                require(
                    any(overlap(frequency) for frequency in block),
                    "a primitive four-gauge block lost common support",
                )
    require(
        primitive_counts == {91: 72, 169: 156},
        "primitive colour counts changed",
    )

    # The actual diagonal row samples remain complementary.
    lam = 1
    b = 3
    shell = 2 * 13**5
    multiplier = 170
    a_prime = 13
    a_value = a_prime + multiplier * shell
    step = 13**b
    physical = {
        0: [a_prime, a_value],
        1: [a_prime + step, a_value + step],
    }
    divided = {
        address: [frequency // 13 for frequency in row]
        for address, row in physical.items()
    }
    require(
        divided == {
            0: [1, 9710741],
            1: [170, 9710910],
        },
        "selected divided gauges changed",
    )
    require(
        divided[0][0] % 169 == divided[0][1] % 169 == 1,
        "the two physical row bases left the primitive 169-colour",
    )
    require(
        [value % 91 for row in divided.values() for value in row]
        == [1, 40, 79, 27],
        "selected primitive colours changed",
    )
    require(
        all(
            gcd(value, 91) == 1
            for row in divided.values()
            for value in row
        ),
        "a selected gauge lost primitive colour",
    )
    for frequency in divided[0]:
        require(
            g_zero(frequency) and not f_zero(frequency),
            "row zero lost g-only cancellation",
        )
    for frequency in divided[1]:
        require(
            f_zero(frequency) and not g_zero(frequency),
            "row one lost f-only cancellation",
        )

    # At the two f-zero samples, the common left and target-only right
    # endpoints differ by 1/170.  Their frequency times that difference
    # is an odd integer, so both endpoint phases equal -1 and opposite
    # jumps cancel inside the total word current.
    endpoint_gap = f_right - f_left
    require(endpoint_gap == Fraction(1, 170), "endpoint gap changed")
    for frequency in divided[1]:
        quotient = frequency * endpoint_gap
        require(
            quotient.denominator == 1
            and quotient.numerator % 2 == 1,
            "aligned and target-only endpoint phases no longer coincide",
        )

    # In K[X]/(X^2-1), F is in the - eigenspace and W in the +
    # eigenspace.  Every central idempotent preserves their zero product.
    full_current = [1, -1]
    word_current = [1, 1]
    product = [1, 0, -1]
    require(
        [
            full_current[0] * word_current[0],
            full_current[0] * word_current[1]
            + full_current[1] * word_current[0],
            full_current[1] * word_current[1],
        ]
        == product,
        "row-current product changed",
    )
    require(
        [product[0] + product[2], product[1]] == [0, 0],
        "row-current product is not zero modulo X^2-1",
    )

    print("interval_J=[-1/4420,1/4420]")
    print("interval_I=[-1/4420,1/2210]")
    print("literal_pullback=E_intersect_Tminus2_Q_equals_J")
    print("xor_common_endpoint=-1/4420")
    print("target_only_endpoint=1/4420")
    print("P1_f_support=[-1/340,1/340]")
    print("P1_g_support=[-1/340,1/170]_plus_half_shift")
    print("THM2323_hypotheses=zero_le_f_le_g_support_in_D2")
    print("primitive_modulus=91")
    print("primitive_colours=72")
    print("physical_row_modulus=169")
    print("physical_row_primitive_colours=156")
    print("common_support_gap_at_most=4_gauges")
    print("selected_colours=1,40,79,27")
    print("selected_row0=g_zero_f_nonzero")
    print("selected_row1=f_zero_g_nonzero")
    print("row_algebra=F_in_minus_W_in_plus")
    print("central_idempotents=diagnose_but_cannot_repair_zero_product")
    print("missing_sidecar=gauge_block_or_private_relative_conductor")
    print("scope=local_word_control_not_full_THM2305_cover")
    print("status=LRC_FIXED_COLOUR_GAUGE_THICKNESS_HOSTILE_EXACT")


if __name__ == "__main__":
    main()
