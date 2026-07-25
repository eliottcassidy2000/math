#!/usr/bin/env python3
"""Exact referee for THM-2283.

The theorem proves a positive mixed safe-box floor on every saturated
rank-two scalar relation torus.  This companion checks all finite algebra:

* the cubic danger-count minorant and every sparse-support case floor;
* the numerical inequalities in the delicate six-coordinate core;
* the squared-Fejer/Jackson coefficient formula by direct convolution;
* the elementary Machin upper bound for pi used by the sharp certificate;
* exact Jackson-ledger scans through the first accepted bandwidths; and
* the scalar and fixed-section relation heights.

The geometric projection and sparse-line classification are proved in the
theorem document.  All checks use explicit ``require`` calls so that normal
and optimized Python executions have the same validity gate.
"""

from fractions import Fraction
from math import comb
import sys


GUARD_DANGER = Fraction(2, 7)
ORDINARY_DANGER = Fraction(1, 7)
GUARD_SAFE = 1 - GUARD_DANGER
ORDINARY_SAFE = 1 - ORDINARY_DANGER
ORDINARY_COORDINATES = 8
TOTAL_COORDINATES = 9

MIXED_TORUS_FLOOR = Fraction(72, 16_807)
TELESCOPE_FACTOR = 2 * TOTAL_COORDINATES
ETA_TARGET = MIXED_TORUS_FLOOR / TELESCOPE_FACTOR

COARSE_PI_UPPER_NUMERATOR = 355
COARSE_PI_UPPER_DENOMINATOR = 113
SHARP_PI_UPPER_NUMERATOR = 104_348
SHARP_PI_UPPER_DENOMINATOR = 33_215

SHARP_JACKSON_N = 1_771
COARSE_JACKSON_N = 1_772
SCALAR_RELATION_HEIGHT = 2 * SHARP_JACKSON_N - 2
FIXED_SECTION_HEIGHT = 2 * SCALAR_RELATION_HEIGHT


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def danger_minorant(count: int) -> Fraction:
    """Cubic pointwise minorant of the zero-danger atom."""
    return -Fraction((count - 1) * (count - 3) * (count - 4), 12)


def danger_minorant_binomial(count: int) -> Fraction:
    """The same cubic in the factorial-moment basis."""
    return (
        1
        - count
        + Fraction(5, 6) * comb(count, 2)
        - Fraction(1, 2) * comb(count, 3)
    )


def alternating_arctangent_partial(denominator: int, last_index: int) -> Fraction:
    """Sum atan(1/q) through the indicated alternating-series index."""
    require(denominator > 1, "arctangent denominator must exceed one")
    require(last_index >= 0, "arctangent truncation index must be nonnegative")
    terms = [
        Fraction(1, (2 * index + 1) * denominator ** (2 * index + 1))
        for index in range(last_index + 1)
    ]
    require(
        all(left > right for left, right in zip(terms, terms[1:])),
        "arctangent term magnitudes must strictly decrease",
    )
    return sum(
        (
            term if index % 2 == 0 else -term
            for index, term in enumerate(terms)
        ),
        Fraction(0),
    )


def jackson_coefficient(bandwidth: int, frequency: int) -> int:
    """Closed integer coefficient C_k of the squared Fejer kernel."""
    require(bandwidth >= 2, "Jackson bandwidth must be at least two")
    require(
        0 <= frequency <= 2 * bandwidth - 2,
        "Jackson coefficient frequency outside support",
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


def jackson_coefficient_by_convolution(bandwidth: int, frequency: int) -> int:
    """Independent convolution of the two triangular Fejer coefficient lists."""
    return sum(
        (bandwidth - abs(index))
        * (bandwidth - abs(frequency - index))
        for index in range(-(bandwidth - 1), bandwidth)
        if abs(frequency - index) < bandwidth
    )


def jackson_odd_sum(bandwidth: int) -> tuple[int, Fraction]:
    """Return C_0 and the exact odd-mode sum in the distance moment."""
    c_zero = bandwidth * (2 * bandwidth * bandwidth + 1) // 3
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(bandwidth, frequency), frequency**2)
            for frequency in range(1, 2 * bandwidth - 2, 2)
        ),
        Fraction(0),
    )
    return c_zero, odd_sum


def jackson_eta_pi_cap_from_sum(
    c_zero: int,
    odd_sum: Fraction,
    pi_upper_numerator: int,
    pi_upper_denominator: int,
) -> Fraction:
    """Strict eta upper bound obtained from a strict rational upper bound on pi."""
    return Fraction(1, 2) - Fraction(
        4 * pi_upper_denominator**2,
        pi_upper_numerator**2 * c_zero,
    ) * odd_sum


def main() -> None:
    # Freeze transcript newlines across Windows and POSIX hosts.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    # Cubic minorant on the complete nine-coordinate danger-count range.
    for count in range(TOTAL_COORDINATES + 1):
        require(
            danger_minorant(count) == danger_minorant_binomial(count),
            f"danger-minorant binomial identity failed at count={count}",
        )
        require(
            danger_minorant(count) <= (1 if count == 0 else 0),
            f"danger minorant failed at count={count}",
        )

    p_zero = GUARD_DANGER
    p = ORDINARY_DANGER
    first_moment = p_zero + ORDINARY_COORDINATES * p
    second_moment = (
        ORDINARY_COORDINATES * p_zero * p
        + comb(ORDINARY_COORDINATES, 2) * p**2
    )
    third_moment = (
        comb(ORDINARY_COORDINATES, 2) * p_zero * p**2
        + comb(ORDINARY_COORDINATES, 3) * p**3
    )
    baseline = (
        1
        - first_moment
        + Fraction(5, 6) * second_moment
        - Fraction(1, 2) * third_moment
    )
    require(first_moment == Fraction(10, 7), "mixed first moment changed")
    require(second_moment == Fraction(44, 49), "mixed second moment changed")
    require(third_moment == Fraction(16, 49), "mixed third moment changed")
    require(baseline == Fraction(23, 147), "mixed baseline changed")

    # Case B: one rational sparse line with support three.
    case_b_ordinary = baseline - Fraction(1, 2) * (p**2 - p**3)
    case_b_guard = baseline - Fraction(1, 2) * (p**2 - p_zero * p**2)
    require(
        case_b_ordinary == Fraction(152, 1029),
        "ordinary support-three floor changed",
    )
    require(
        case_b_guard == Fraction(307, 2058),
        "guard support-three floor changed",
    )

    # Case C: one rational sparse line with support two.
    ordinary_pair_coefficient = Fraction(5, 6) - Fraction(1, 2) * (
        p_zero + 6 * p
    )
    guard_pair_coefficient = Fraction(5, 6) - Fraction(1, 2) * (7 * p)
    require(
        ordinary_pair_coefficient == Fraction(11, 42),
        "ordinary-pair correction coefficient changed",
    )
    require(
        guard_pair_coefficient == Fraction(1, 3),
        "guard-pair correction coefficient changed",
    )
    case_c_ordinary = baseline - ordinary_pair_coefficient * p**2
    case_c_guard = baseline - guard_pair_coefficient * p_zero * p
    require(
        case_c_ordinary == Fraction(311, 2058),
        "ordinary support-two floor changed",
    )
    require(
        case_c_guard == Fraction(1, 7),
        "guard support-two floor changed",
    )

    for case_name, case_floor in (
        ("A", baseline),
        ("B ordinary", case_b_ordinary),
        ("B guard", case_b_guard),
        ("C ordinary", case_c_ordinary),
        ("C guard", case_c_guard),
    ):
        require(
            case_floor >= MIXED_TORUS_FLOOR,
            f"sparse case {case_name} fell below the mixed torus target",
        )

    # Case D: two independent sparse lines.  The union bound handles every
    # support-union type except m=6 with the guard in the core.
    union_floor_table: list[tuple[int, Fraction, Fraction]] = []
    for support_union_size in range(3, 7):
        guard_in_core = (
            1 - p_zero - (support_union_size - 1) * p
        ) * ORDINARY_SAFE ** (TOTAL_COORDINATES - support_union_size)
        guard_outside_core = (
            (1 - support_union_size * p)
            * GUARD_SAFE
            * ORDINARY_SAFE ** (ORDINARY_COORDINATES - support_union_size)
        )
        union_floor_table.append(
            (support_union_size, guard_in_core, guard_outside_core)
        )
        if support_union_size < 6:
            require(
                guard_in_core >= MIXED_TORUS_FLOOR,
                f"guard-in union floor failed at m={support_union_size}",
            )
        require(
            guard_outside_core >= MIXED_TORUS_FLOOR,
            f"guard-out union floor failed at m={support_union_size}",
        )

    expected_union_table = (
        (3, Fraction(139_968, 823_543), Fraction(155_520, 823_543)),
        (4, Fraction(15_552, 117_649), Fraction(19_440, 117_649)),
        (5, Fraction(1_296, 16_807), Fraction(2_160, 16_807)),
        (6, Fraction(0), Fraction(180, 2_401)),
    )
    require(
        tuple(union_floor_table) == expected_union_table,
        "case-D union-bound table changed",
    )

    # Delicate m=6 guard core: its danger count has mean one.  Therefore
    # Pr(N=0)=E[(N-1)_+], and the following pointwise inequality prices a
    # single Haar pair.  Rank two cannot contain all five independent star
    # pair relations on six nonzero coordinates, so at least one of the
    # fifteen pair projections is Haar.
    delicate_core_mean = p_zero + 5 * p
    require(delicate_core_mean == 1, "delicate core mean changed")
    for count in range(7):
        positive_part = max(count - 1, 0)
        require(
            Fraction(positive_part) >= Fraction(comb(count, 2), 3),
            f"delicate pointwise hinge failed at count={count}",
        )
    star_pair_rank = 6 - 1
    require(
        star_pair_rank > 2,
        "rank-two plane could contain all six-coordinate star pair relations",
    )
    pair_count = comb(6, 2)
    minimum_haar_pair_mass = min(p**2, p_zero * p)
    delicate_core_floor = minimum_haar_pair_mass / 3
    delicate_full_floor = delicate_core_floor * ORDINARY_SAFE**3
    require(pair_count == 15, "delicate pair count changed")
    require(
        minimum_haar_pair_mass == Fraction(1, 49),
        "minimum Haar-pair danger mass changed",
    )
    require(
        delicate_core_floor == Fraction(1, 147),
        "delicate core safe floor changed",
    )
    require(
        delicate_full_floor == MIXED_TORUS_FLOOR,
        "delicate full mixed safe floor changed",
    )

    # Machin's identity is checked at the rational tangent level:
    # tan(2 atan(1/5))=5/12, tan(4 atan(1/5))=120/119, and subtracting
    # atan(1/239) gives tangent one.  The relevant angle lies in (0,pi/2).
    tangent_one_fifth = Fraction(1, 5)
    tangent_double = (
        2 * tangent_one_fifth / (1 - tangent_one_fifth**2)
    )
    tangent_quadruple = 2 * tangent_double / (1 - tangent_double**2)
    tangent_one_239 = Fraction(1, 239)
    tangent_difference = (
        (tangent_quadruple - tangent_one_239)
        / (1 + tangent_quadruple * tangent_one_239)
    )
    require(tangent_double == Fraction(5, 12), "Machin double angle changed")
    require(
        tangent_quadruple == Fraction(120, 119),
        "Machin quadruple angle changed",
    )
    require(tangent_difference == 1, "Machin tangent identity changed")

    # Alternating-series directions are parity-sensitive.  Index 6 is even,
    # hence an upper bound for atan(1/5); index 1 is odd, hence a lower bound
    # for atan(1/239).  This gives a strict upper bound for pi.
    atan_one_fifth_upper = alternating_arctangent_partial(5, 6)
    atan_one_fifth_lower = alternating_arctangent_partial(5, 7)
    atan_one_239_lower = alternating_arctangent_partial(239, 1)
    atan_one_239_upper = alternating_arctangent_partial(239, 2)
    require(
        atan_one_fifth_lower < atan_one_fifth_upper,
        "atan(1/5) alternating bracket direction failed",
    )
    require(
        atan_one_239_lower < atan_one_239_upper,
        "atan(1/239) alternating bracket direction failed",
    )
    machin_angle_lower = 4 * atan_one_fifth_lower - atan_one_239_upper
    machin_angle_upper = 4 * atan_one_fifth_upper - atan_one_239_lower
    require(
        0 < machin_angle_lower < machin_angle_upper < 1,
        "Machin tangent branch was not confined to the first quadrant",
    )
    machin_truncation_upper = 4 * (
        machin_angle_upper
    )
    sharp_pi_upper = Fraction(
        SHARP_PI_UPPER_NUMERATOR,
        SHARP_PI_UPPER_DENOMINATOR,
    )
    require(
        atan_one_fifth_upper
        == Fraction(2_170_821_043_343, 10_997_314_453_125),
        "atan(1/5) upper truncation changed",
    )
    require(
        atan_one_239_lower == Fraction(171_362, 40_955_757),
        "atan(1/239) lower truncation changed",
    )
    require(
        machin_truncation_upper
        == Fraction(
            471_661_273_023_004_128_472,
            150_134_446_131_591_796_875,
        ),
        "Machin truncation upper bound changed",
    )
    machin_ratio_margin = sharp_pi_upper - machin_truncation_upper
    require(
        machin_ratio_margin
        == Fraction(
            464_759_401_292,
            1_565_687_795_372_314_453_125,
        ),
        "Machin rational-upper margin changed",
    )
    require(machin_ratio_margin > 0, "sharp rational pi bound failed")

    # Independent Jackson coefficient reconstruction at both load-bearing
    # bandwidths.
    for bandwidth in (SHARP_JACKSON_N, COARSE_JACKSON_N):
        for frequency in range(2 * bandwidth - 1):
            closed = jackson_coefficient(bandwidth, frequency)
            convolved = jackson_coefficient_by_convolution(
                bandwidth,
                frequency,
            )
            require(
                closed == convolved,
                "Jackson convolution mismatch "
                f"at N={bandwidth}, frequency={frequency}",
            )
            require(
                closed > 0,
                f"Jackson coefficient vanished at N={bandwidth}, k={frequency}",
            )

    # Exact scans use the same odd sums for the coarse and sharp pi bounds.
    sharp_accepted: list[int] = []
    coarse_accepted: list[int] = []
    selected_caps: dict[tuple[str, int], Fraction] = {}
    for bandwidth in range(2, COARSE_JACKSON_N + 1):
        c_zero, odd_sum = jackson_odd_sum(bandwidth)
        require(
            c_zero == jackson_coefficient(bandwidth, 0),
            f"Jackson C_0 mismatch at N={bandwidth}",
        )
        sharp_cap = jackson_eta_pi_cap_from_sum(
            c_zero,
            odd_sum,
            SHARP_PI_UPPER_NUMERATOR,
            SHARP_PI_UPPER_DENOMINATOR,
        )
        coarse_cap = jackson_eta_pi_cap_from_sum(
            c_zero,
            odd_sum,
            COARSE_PI_UPPER_NUMERATOR,
            COARSE_PI_UPPER_DENOMINATOR,
        )
        require(sharp_cap > 0, f"sharp eta cap nonpositive at N={bandwidth}")
        require(coarse_cap > 0, f"coarse eta cap nonpositive at N={bandwidth}")
        if bandwidth <= SHARP_JACKSON_N and sharp_cap < ETA_TARGET:
            sharp_accepted.append(bandwidth)
        if coarse_cap < ETA_TARGET:
            coarse_accepted.append(bandwidth)
        if bandwidth in (1_770, 1_771, 1_772):
            selected_caps[("sharp", bandwidth)] = sharp_cap
            selected_caps[("coarse", bandwidth)] = coarse_cap

    require(
        sharp_accepted == [SHARP_JACKSON_N],
        "sharp exact scan did not first accept N=1771",
    )
    require(
        coarse_accepted == [COARSE_JACKSON_N],
        "355/113 exact scan did not first accept N=1772",
    )

    sharp_eta_1770 = selected_caps[("sharp", 1_770)]
    sharp_eta_1771 = selected_caps[("sharp", 1_771)]
    coarse_eta_1771 = selected_caps[("coarse", 1_771)]
    coarse_eta_1772 = selected_caps[("coarse", 1_772)]

    require(
        sharp_eta_1770 > Fraction(11_903, 50_000_000) > ETA_TARGET,
        "sharp N=1770 boundary no longer fails",
    )
    require(
        Fraction(11_903, 50_000_000) - ETA_TARGET
        == Fraction(53_721, 840_350_000_000),
        "sharp N=1770 simple failure margin changed",
    )
    require(
        sharp_eta_1771 < Fraction(23_794, 100_000_000),
        "sharp N=1771 simple eta cap changed",
    )
    sharp_simple_gap = (
        MIXED_TORUS_FLOOR
        - TELESCOPE_FACTOR * Fraction(23_794, 100_000_000)
    )
    require(
        sharp_simple_gap == Fraction(424_089, 420_175_000_000),
        "sharp N=1771 simple gap changed",
    )
    require(
        sharp_simple_gap > Fraction(1, 1_000_000),
        "sharp N=1771 gap no longer exceeds one millionth",
    )
    require(
        sharp_simple_gap - Fraction(1, 1_000_000)
        == Fraction(1_957, 210_087_500_000),
        "sharp millionth comparison changed",
    )

    require(
        coarse_eta_1771 > Fraction(119, 500_000) > ETA_TARGET,
        "coarse N=1771 boundary no longer fails",
    )
    require(
        Fraction(119, 500_000) - ETA_TARGET
        == Fraction(33, 8_403_500_000),
        "coarse N=1771 simple failure margin changed",
    )
    require(
        coarse_eta_1772 < Fraction(2_379, 10_000_000),
        "coarse N=1772 simple eta cap changed",
    )
    coarse_simple_gap = (
        MIXED_TORUS_FLOOR
        - TELESCOPE_FACTOR * Fraction(2_379, 10_000_000)
    )
    require(
        coarse_simple_gap == Fraction(145_323, 84_035_000_000),
        "coarse N=1772 simple gap changed",
    )
    require(
        coarse_simple_gap > Fraction(1, 600_000),
        "coarse N=1772 gap no longer exceeds 1/600000",
    )
    require(
        coarse_simple_gap - Fraction(1, 600_000)
        == Fraction(7_897, 126_052_500_000),
        "coarse 1/600000 comparison changed",
    )

    require(ETA_TARGET == Fraction(4, 16_807), "eta target changed")
    require(SCALAR_RELATION_HEIGHT == 3_540, "scalar height changed")
    require(FIXED_SECTION_HEIGHT == 7_080, "fixed-section height changed")

    print("THM-2283 exact mixed safe-torus and scalar rank-three referee")
    print("danger_minorant=-(N-1)(N-3)(N-4)/12")
    print(f"mixed_three_wise_baseline={baseline}")
    print(f"case_B_support3_ordinary={case_b_ordinary}")
    print(f"case_B_support3_guard={case_b_guard}")
    print(f"case_C_support2_ordinary={case_c_ordinary}")
    print(f"case_C_support2_guard={case_c_guard}")
    print(
        "case_D_union_floors="
        + ";".join(
            f"m{size}:guard_in={guard_in},guard_out={guard_out}"
            for size, guard_in, guard_out in union_floor_table
        )
    )
    print(
        f"delicate_core_mean={delicate_core_mean},"
        f"pair_count={pair_count},star_pair_rank={star_pair_rank}"
    )
    print(
        f"delicate_haar_pair_floor={minimum_haar_pair_mass},"
        f"core_safe_floor={delicate_core_floor},"
        f"full_safe_floor={delicate_full_floor}"
    )
    print("Machin_directions=atan(1/5):k6_upper;atan(1/239):k1_lower")
    print(f"Machin_truncation_pi_upper={machin_truncation_upper}")
    print(
        f"Machin_pi_upper<{SHARP_PI_UPPER_NUMERATOR}/"
        f"{SHARP_PI_UPPER_DENOMINATOR},margin={machin_ratio_margin}"
    )
    print("sharp_Jackson_scan_N=2..1771:first_accepted=1771")
    print(
        "sharp_N1770_FAIL:"
        "eta_cap>11903/50000000>4/16807"
    )
    print(
        "sharp_N1771_PASS:"
        "eta_cap<23794/100000000,"
        f"gap>{sharp_simple_gap}>1/1000000"
    )
    print("coarse_355_over_113_scan_N=2..1772:first_accepted=1772")
    print(
        "coarse_N1771_FAIL:"
        "eta_cap>119/500000>4/16807"
    )
    print(
        "coarse_N1772_PASS:"
        "eta_cap<2379/10000000,"
        f"gap>{coarse_simple_gap}>1/600000"
    )
    print(
        f"scalar_relation_height={SCALAR_RELATION_HEIGHT},"
        f"fixed_section_height={FIXED_SECTION_HEIGHT}"
    )
    print("coarse_fallback_scalar_height=3542,fixed_section_height=7084")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
