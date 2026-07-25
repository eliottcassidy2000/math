#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2193.

The theorem document proves the Haar disintegration, slicing, cube inclusion,
and Jackson-kernel inequalities. This companion keeps every numerical factor
live under both ordinary Python and ``python -O``.
"""

from fractions import Fraction


P = Fraction(1, 7)
Q = Fraction(6, 7)
DIMENSION = 7
EXTRA_CHARACTERS = 6
LOCAL_COEFFICIENT_BOUND = 6
LOCAL_L1_BOUND = DIMENSION * LOCAL_COEFFICIENT_BOUND
CUBE_RADIUS = Fraction(1, 686)
CUBE_SIDE = 2 * CUBE_RADIUS
UNIFORM_FLOOR = Fraction(1, 7**21)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    require(CUBE_SIDE == Fraction(1, 343), "cube side changed")
    require(343**7 == 7**21, "cube-volume power identity changed")
    require(
        CUBE_SIDE**DIMENSION == UNIFORM_FLOOR,
        "cube volume is not 7^-21",
    )
    require(
        LOCAL_L1_BOUND * CUBE_RADIUS == Fraction(3, 49),
        "local linear-form variation changed",
    )

    # Case 1: the essential residue-count minimum.
    for h in range(7, 100_001):
        require(
            Fraction(h // 7, h) >= Fraction(1, 13),
            f"residue-count floor failed at h={h}",
        )
    strong_facet_floor = P * Q**6 / 13

    # Case 2: both possible locations of a maximal coefficient.
    for coefficient in range(7, 100_001):
        danger_slice = Fraction(coefficient // 7, coefficient)
        safe_slice = Fraction((6 * coefficient) // 7, coefficient)
        require(
            danger_slice >= Fraction(1, 13),
            f"danger-coordinate slice failed at B={coefficient}",
        )
        require(
            safe_slice >= Fraction(6, 13),
            f"safe-coordinate slice failed at B={coefficient}",
        )
    require(
        Q**6 / 91 == strong_facet_floor,
        "danger-slice domain factor changed",
    )
    require(
        P * Q**5 * Fraction(6, 91) == strong_facet_floor,
        "safe-slice domain factor changed",
    )

    # Case 3: parity centre, displacement, and open-cube margins.
    require(
        Fraction(1, 14) + Fraction(3, 7) == Fraction(1, 2),
        "minimum odd support radius changed",
    )
    for h in range(1, 7):
        target_residual = Fraction(h, 98)
        domain_margin = Fraction(h, 686)
        target_margin = Fraction(h, 14) - target_residual
        require(
            target_margin == Fraction(3 * h, 49),
            f"target margin changed at h={h}",
        )
        require(
            domain_margin >= CUBE_RADIUS,
            f"domain cube no longer fits at h={h}",
        )
        require(
            target_margin >= LOCAL_L1_BOUND * CUBE_RADIUS,
            f"target cube no longer fits at h={h}",
        )
        theta_upper = 1 - Fraction(h, 49)
        require(0 < theta_upper < 1, f"odd shift failed at h={h}")

    local_facet_floor = CUBE_SIDE**DIMENSION / 6
    require(
        strong_facet_floor > local_facet_floor,
        "strong cases no longer dominate local case",
    )

    # The factor six is load-bearing at the critical union cancellation.
    require(Q**7 == EXTRA_CHARACTERS * P * Q**6, "critical cancellation")
    require(
        EXTRA_CHARACTERS * local_facet_floor == UNIFORM_FLOOR,
        "six facet deficits do not sum to 7^-21",
    )
    require(7**21 == 558_545_864_083_284_007, "7^21 changed")

    # Jackson first-moment bound. For the pointwise majorant in the theorem,
    # the two exact pieces are 3/(8N) and 3/(8N)-3/(8N^3).
    for n in range(2, 101):
        c_zero_by_sum = sum(
            (n - abs(k)) ** 2 for k in range(-(n - 1), n)
        )
        c_zero_closed = n * (2 * n * n + 1) // 3
        require(
            c_zero_by_sum == c_zero_closed,
            f"Jackson normalization failed at N={n}",
        )
        low_moment = Fraction(3, 8 * n)
        high_moment = Fraction(3, 8 * n) - Fraction(3, 8 * n**3)
        moment_cap = low_moment + high_moment
        eta_cap = 2 * moment_cap
        require(
            moment_cap
            == Fraction(3, 4 * n) - Fraction(3, 8 * n**3),
            f"Jackson moment split failed at N={n}",
        )
        require(
            eta_cap < Fraction(3, 2 * n),
            f"Jackson eta cap failed at N={n}",
        )

    n_star = 39 * 7**21 + 1
    h_star = 2 * n_star - 2
    require(
        n_star == 21_783_288_699_248_076_274,
        "explicit Jackson N changed",
    )
    require(
        h_star == 78 * 7**21,
        "degree identity changed",
    )
    require(
        h_star == 43_566_577_398_496_152_546,
        "explicit rank-seven height changed",
    )
    require(
        Fraction(39, n_star) < UNIFORM_FLOOR,
        "Jackson error does not clear uniform floor",
    )

    # THM-2190's rank-five floor handles all lower ranks at the same height.
    rank_five_floor = Fraction(279_936, 5_764_801)
    require(
        rank_five_floor > UNIFORM_FLOOR,
        "rank-five floor no longer dominates uniform floor",
    )

    deletion_height = 900 * h_star
    require(
        deletion_height == 39_209_919_658_646_537_291_400,
        "deletion height changed",
    )

    print("THM-2193 exact uniform rank-six safe-torus referee")
    print(f"strong_facet_floor={strong_facet_floor}")
    print(f"local_cube_radius={CUBE_RADIUS}, side={CUBE_SIDE}")
    print(f"local_facet_floor={local_facet_floor}")
    print(f"uniform_rank_six_floor={UNIFORM_FLOOR}=7^-21")
    print("critical_deficit_count=6")
    print("residue_h=7..100000 and slice_B=7..100000 exact sweeps passed")
    print("local_h=1..6 parity/margin audit passed")
    print("Jackson_eta_N<3/(2N)")
    print(f"N_star={n_star}")
    print(f"H_star={h_star}")
    print(f"deletion_rank=6, deletion_height={deletion_height}")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
