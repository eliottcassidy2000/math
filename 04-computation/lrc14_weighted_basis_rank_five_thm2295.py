#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2295.

Checks the weighted coordinate-basis safe floors, the Machin-certified
squared-Fejer error at N=99, the adjacent N=98 failure of this exact
certificate, and the scalar/fixed-section height arithmetic.  All validity
checks use ``require`` so normal and optimized Python execute the same gate.
"""

from fractions import Fraction
import sys


GUARD_DANGER = Fraction(2, 7)
ORDINARY_DANGER = Fraction(1, 7)
GUARD_SAFE = 1 - GUARD_DANGER
ORDINARY_SAFE = 1 - ORDINARY_DANGER
TOTAL_COORDINATES = 9
TELESCOPE_FACTOR = 2 * TOTAL_COORDINATES

SHARP_PI_UPPER_NUMERATOR = 104_348
SHARP_PI_UPPER_DENOMINATOR = 33_215
JACKSON_N = 99
SCALAR_HEIGHT = 2 * JACKSON_N - 2
FIXED_SECTION_HEIGHT = 2 * SCALAR_HEIGHT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def alternating_arctangent_partial(
    denominator: int,
    last_index: int,
) -> Fraction:
    require(denominator > 1, "arctangent denominator must exceed one")
    require(last_index >= 0, "arctangent index must be nonnegative")
    terms = [
        Fraction(1, (2 * index + 1) * denominator ** (2 * index + 1))
        for index in range(last_index + 1)
    ]
    require(
        all(left > right for left, right in zip(terms, terms[1:])),
        "arctangent terms must strictly decrease",
    )
    return sum(
        (
            term if index % 2 == 0 else -term
            for index, term in enumerate(terms)
        ),
        Fraction(0),
    )


def jackson_coefficient(bandwidth: int, frequency: int) -> int:
    require(bandwidth >= 2, "Jackson bandwidth must be at least two")
    require(
        0 <= frequency <= 2 * bandwidth - 2,
        "Jackson frequency outside support",
    )
    if frequency <= bandwidth:
        return (
            4 * bandwidth**3
            - 6 * bandwidth * frequency**2
            + 2 * bandwidth
            + 3 * frequency**3
            - 3 * frequency
        ) // 6
    reflected = 2 * bandwidth - frequency
    return (reflected**3 - reflected) // 6


def jackson_coefficient_by_convolution(
    bandwidth: int,
    frequency: int,
) -> int:
    return sum(
        (bandwidth - abs(index))
        * (bandwidth - abs(frequency - index))
        for index in range(-(bandwidth - 1), bandwidth)
        if abs(frequency - index) < bandwidth
    )


def jackson_odd_sum(bandwidth: int) -> tuple[int, Fraction]:
    c_zero = bandwidth * (2 * bandwidth * bandwidth + 1) // 3
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(bandwidth, frequency), frequency**2)
            for frequency in range(1, 2 * bandwidth - 2, 2)
        ),
        Fraction(0),
    )
    return c_zero, odd_sum


def jackson_eta_cap(bandwidth: int) -> Fraction:
    c_zero, odd_sum = jackson_odd_sum(bandwidth)
    return Fraction(1, 2) - Fraction(
        4 * SHARP_PI_UPPER_DENOMINATOR**2,
        SHARP_PI_UPPER_NUMERATOR**2 * c_zero,
    ) * odd_sum


def mixed_basis_floor(rank: int) -> Fraction:
    require(0 <= rank <= 4, "weighted basis floor rank outside 0..4")
    return (
        Fraction(5 - rank, 7)
        * ORDINARY_SAFE ** (8 - rank)
    )


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    expected_floors = (
        Fraction(8_398_080, 40_353_607),
        Fraction(1_119_744, 5_764_801),
        Fraction(139_968, 823_543),
        Fraction(15_552, 117_649),
        Fraction(1_296, 16_807),
    )
    floors = tuple(mixed_basis_floor(rank) for rank in range(5))
    require(floors == expected_floors, "weighted floor table changed")
    require(
        all(left > right for left, right in zip(floors, floors[1:])),
        "weighted floors must decrease through rank four",
    )

    for rank in range(5):
        dimension = TOTAL_COORDINATES - rank
        guard_in_basis = (
            GUARD_SAFE * ORDINARY_SAFE ** (dimension - 1)
            - rank
            * ORDINARY_DANGER
            * ORDINARY_SAFE ** (dimension - 1)
        )
        require(
            guard_in_basis == floors[rank],
            f"guard-in-basis ledger failed at rank {rank}",
        )
        if rank >= 1:
            guard_extra = (
                ORDINARY_SAFE**dimension
                - (
                    GUARD_DANGER
                    + (rank - 1) * ORDINARY_DANGER
                )
                * ORDINARY_SAFE ** (dimension - 1)
            )
            require(
                guard_extra == floors[rank],
                f"guard-extra ledger failed at rank {rank}",
            )

    tangent_one_fifth = Fraction(1, 5)
    tangent_double = (
        2 * tangent_one_fifth / (1 - tangent_one_fifth**2)
    )
    tangent_quadruple = (
        2 * tangent_double / (1 - tangent_double**2)
    )
    tangent_one_239 = Fraction(1, 239)
    tangent_difference = (
        (tangent_quadruple - tangent_one_239)
        / (1 + tangent_quadruple * tangent_one_239)
    )
    require(tangent_double == Fraction(5, 12), "Machin double changed")
    require(
        tangent_quadruple == Fraction(120, 119),
        "Machin quadruple changed",
    )
    require(tangent_difference == 1, "Machin tangent identity failed")

    atan_one_fifth_upper = alternating_arctangent_partial(5, 6)
    atan_one_239_lower = alternating_arctangent_partial(239, 1)
    machin_pi_upper = 4 * (
        4 * atan_one_fifth_upper - atan_one_239_lower
    )
    certified_pi_upper = Fraction(
        SHARP_PI_UPPER_NUMERATOR,
        SHARP_PI_UPPER_DENOMINATOR,
    )
    require(
        0 < machin_pi_upper < certified_pi_upper,
        "Machin rational upper bound for pi failed",
    )

    for bandwidth in (98, 99):
        for frequency in range(2 * bandwidth - 1):
            require(
                jackson_coefficient(bandwidth, frequency)
                == jackson_coefficient_by_convolution(
                    bandwidth,
                    frequency,
                ),
                "Jackson closed/convolution coefficient mismatch "
                f"at N={bandwidth}, k={frequency}",
            )

    worst_floor = floors[-1]
    accepted = []
    selected = {}
    for bandwidth in range(2, JACKSON_N + 1):
        eta_cap = jackson_eta_cap(bandwidth)
        require(eta_cap > 0, f"eta cap nonpositive at N={bandwidth}")
        gap = worst_floor - TELESCOPE_FACTOR * eta_cap
        if gap > 0:
            accepted.append(bandwidth)
        if bandwidth in (98, 99):
            selected[bandwidth] = (eta_cap, gap)

    require(
        accepted == [JACKSON_N],
        "N=99 was not the first passing exact certificate",
    )
    eta_98, gap_98 = selected[98]
    eta_99, gap_99 = selected[99]
    require(
        eta_98 > Fraction(43, 10_000) > Fraction(1, 234),
        "N=98 eta failure wrapper changed",
    )
    require(
        eta_99 < Fraction(213, 50_000) < Fraction(1, 234),
        "N=99 eta success wrapper changed",
    )
    require(
        Fraction(-301, 1_000_000)
        < gap_98
        < Fraction(-300, 1_000_000)
        < 0,
        "N=98 failure wrapper changed",
    )
    require(
        0
        < Fraction(481, 1_000_000)
        < gap_99
        < Fraction(482, 1_000_000),
        "N=99 success wrapper changed",
    )
    require(
        gap_99 > Fraction(1, 2_500),
        "N=99 exact margin no longer exceeds 1/2500",
    )
    require(
        eta_99 < Fraction(1, 234),
        "N=99 eta cap no longer below 1/234",
    )

    simple_gap = worst_floor - Fraction(1, 13)
    require(
        simple_gap == Fraction(41, 218_491) > 0,
        "simple rank-four margin changed",
    )
    sharp_simple_gap = (
        worst_floor - TELESCOPE_FACTOR * Fraction(213, 50_000)
    )
    require(
        sharp_simple_gap
        == Fraction(180_981, 420_175_000)
        > 0,
        "sharp simple N=99 margin changed",
    )
    adjacent_failure = (
        worst_floor - TELESCOPE_FACTOR * Fraction(43, 10_000)
    )
    require(
        adjacent_failure
        == Fraction(-24_309, 84_035_000)
        < 0,
        "simple N=98 failure margin changed",
    )
    require(SCALAR_HEIGHT == 196, "scalar height changed")
    require(FIXED_SECTION_HEIGHT == 392, "fixed-section height changed")

    print("THM-2295 weighted basis-safe and scalar rank-five referee")
    print(
        "mixed_basis_floors="
        + ",".join(f"r{rank}:{floor}" for rank, floor in enumerate(floors))
    )
    print(f"worst_rank4_floor={worst_floor}")
    print("Machin_pi_upper=104348/33215")
    print("exact_Jackson_scan_N=2..99:first_accepted=99")
    print("N98_FAIL:-301/1000000<gap<-300/1000000<0")
    print("N99_PASS:1/2500<481/1000000<gap<482/1000000")
    print("eta98>43/10000>1/234;eta99<213/50000<1/234")
    print(
        f"N99_simple_margin={sharp_simple_gap},"
        f"coarse_margin={simple_gap}"
    )
    print(
        f"scalar_rank>=5_height={SCALAR_HEIGHT},"
        f"fixed_section_height={FIXED_SECTION_HEIGHT}"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
