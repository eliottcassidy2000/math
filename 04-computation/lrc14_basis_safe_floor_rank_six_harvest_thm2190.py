#!/usr/bin/env python3
"""Exact referee for THM-2190.

The finite-cover and Haar-disintegration arguments are proved in the theorem
document. This companion checks all load-bearing rational constants, rebuilds
the squared-Fejer/Jackson coefficients by two independent formulas, verifies
the inherited N=251 error ledger, and checks the rank-six parity-radius and
deletion ledgers. All checks remain live under ``python -O``.
"""

from fractions import Fraction


P = Fraction(1, 7)
Q = Fraction(6, 7)
JACKSON_N = 251
PI_UPPER_NUMERATOR = 355
PI_UPPER_DENOMINATOR = 113
ETA_CAP = Fraction(21, 12_500)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def safe_floor(rank: int) -> Fraction:
    """Basis-safe floor q^(d-1)(q-rp), where d=13-r."""
    require(0 <= rank <= 6, "rank outside checked range")
    dimension = 13 - rank
    return Q ** (dimension - 1) * (Q - rank * P)


def jackson_coefficient(n: int, k: int) -> int:
    """Closed coefficient of the square of the triangular Fejer list."""
    require(0 <= k <= 2 * n - 2, "Jackson coefficient index")
    if k <= n:
        return (
            4 * n**3
            - 6 * n * k**2
            + 2 * n
            + 3 * k**3
            - 3 * k
        ) // 6
    return ((2 * n - k) ** 3 - (2 * n - k)) // 6


def jackson_coefficient_by_convolution(n: int, k: int) -> int:
    """Independent convolution of the two triangular coefficient lists."""
    return sum(
        (n - abs(a)) * (n - abs(k - a))
        for a in range(-(n - 1), n)
        if abs(k - a) < n
    )


def jackson_eta_pi_cap(n: int) -> Fraction:
    """Strict eta upper bound obtained from pi < 355/113."""
    c_zero = n * (2 * n * n + 1) // 3
    odd_sum = sum(
        (
            Fraction(jackson_coefficient(n, k), k * k)
            for k in range(1, 2 * n - 2, 2)
        ),
        Fraction(0),
    )
    return Fraction(1, 2) - Fraction(
        4 * PI_UPPER_DENOMINATOR**2,
        PI_UPPER_NUMERATOR**2 * c_zero,
    ) * odd_sum


def main() -> None:
    expected_floors = {
        2: Fraction(241_864_704, 1_977_326_743),
        3: Fraction(30_233_088, 282_475_249),
        4: Fraction(3_359_232, 40_353_607),
        5: Fraction(279_936, 5_764_801),
    }
    for rank, expected in expected_floors.items():
        require(safe_floor(rank) == expected, f"rank-{rank} floor changed")

    require(
        all(
            safe_floor(rank) > safe_floor(rank + 1)
            for rank in range(2, 5)
        ),
        "rank-two-through-five floors are not strictly decreasing",
    )
    require(safe_floor(6) == 0, "rank-six union ledger should be critical")
    require(
        Q**7 - 6 * P * Q**6 == 0,
        "rank-six conditional union cancellation changed",
    )

    # The odd-centre branch of the strict facet lemma has at least this
    # much numerator variation before paying the additional h/14 target
    # radius. The even-centre branch is centered exactly.
    minimum_odd_variation = Fraction(1, 14) + Fraction(3, 7)
    require(
        minimum_odd_variation == Fraction(1, 2),
        "strict facet odd-parity variation changed",
    )
    target_radius_per_h = Fraction(1, 14)
    require(
        target_radius_per_h > 0,
        "strict facet target radius is not positive",
    )

    # Rebuild all degree-500 coefficients independently.
    for k in range(2 * JACKSON_N - 1):
        closed = jackson_coefficient(JACKSON_N, k)
        convolved = jackson_coefficient_by_convolution(JACKSON_N, k)
        require(closed == convolved, f"Jackson mismatch at k={k}")
        require(closed > 0, f"Jackson coefficient not positive at k={k}")

    c_zero = JACKSON_N * (2 * JACKSON_N * JACKSON_N + 1) // 3
    require(
        c_zero == jackson_coefficient(JACKSON_N, 0),
        "Jackson zero coefficient mismatch",
    )
    require(c_zero == 10_542_251, "Jackson C_0 changed")

    eta_cap = jackson_eta_pi_cap(JACKSON_N)
    require(eta_cap > 0, "Jackson eta cap should be positive")
    require(eta_cap < ETA_CAP, "N=251 Jackson eta cap changed")

    smallest_floor = safe_floor(5)
    inherited_error = 26 * ETA_CAP
    comparison_gap = smallest_floor - inherited_error
    require(
        comparison_gap == Fraction(175_809_327, 36_030_006_250),
        "rank-five/Jackson comparison gap changed",
    )
    require(comparison_gap > 0, "rank-five/Jackson gap is not positive")

    optimized_thresholds = {
        2: (90, Fraction(47, 10_000)),
        3: (103, Fraction(41, 10_000)),
        4: (132, Fraction(2, 625)),
        5: (226, Fraction(373, 200_000)),
    }
    for rank, (target_n, simple_cap) in optimized_thresholds.items():
        first_passing = next(
            (
                n
                for n in range(2, target_n + 1)
                if 26 * jackson_eta_pi_cap(n) < safe_floor(rank)
            ),
            None,
        )
        require(
            first_passing == target_n,
            f"first passing Jackson N changed at rank={rank}",
        )
        target_eta = jackson_eta_pi_cap(target_n)
        previous_eta = jackson_eta_pi_cap(target_n - 1)
        require(
            target_eta < simple_cap,
            f"simple eta cap failed at rank={rank}",
        )
        require(
            26 * simple_cap < safe_floor(rank),
            f"simple rank-floor comparison failed at rank={rank}",
        )
        require(
            26 * previous_eta >= safe_floor(rank),
            f"adjacent Jackson failure changed at rank={rank}",
        )

    nested_heights = (105, 105, 178, 204, 262, 450)
    require(
        tuple(
            2 * optimized_thresholds[rank][0] - 2
            for rank in range(2, 6)
        )
        == nested_heights[2:],
        "optimized nested heights changed",
    )
    pivot_bounds = tuple(
        2 * nested_heights[index] * nested_heights[-1]
        for index in range(5)
    )
    require(
        pivot_bounds == (94_500, 94_500, 160_200, 183_600, 235_800),
        "optimized deletion pivot bounds changed",
    )
    optimized_deletion_height = max(
        (*pivot_bounds, nested_heights[-2])
    )
    require(
        optimized_deletion_height == 235_800,
        "optimized deletion height changed",
    )
    simple_deletion_height = 2 * 500**2
    require(simple_deletion_height == 500_000, "simple deletion height changed")

    print("THM-2190 exact basis-safe floor and rank-six-harvest referee")
    print(
        "safe_floors_r2_to_r5="
        + ",".join(str(expected_floors[rank]) for rank in range(2, 6))
    )
    print(f"rank_six_union_ledger={safe_floor(6)}")
    print("rank_six_facet_odd_variation=1/2+h/14>1/2")
    print(
        f"Jackson_N={JACKSON_N}, degree={2 * JACKSON_N - 2}, "
        f"C_0={c_zero}"
    )
    print(f"eta_{JACKSON_N}_pi_cap<{ETA_CAP}")
    print(f"rank_five_Jackson_gap={comparison_gap}")
    print("optimized_first_N_r2_to_r5=90,103,132,226")
    print("optimized_heights_r3_to_r6=178,204,262,450")
    print("adjacent_N_r2_to_r5=89,102,131,225 all fail exact cap")
    print("nested_basis_heights=105,105,178,204,262,450")
    print(
        "deletion_rank=5, optimized_height="
        f"{optimized_deletion_height}, simple_height={simple_deletion_height}"
    )
    print("all exact checks passed")


if __name__ == "__main__":
    main()
