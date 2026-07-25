#!/usr/bin/env python3
"""Exact referee for THM-2185.

The theorem's Haar-projection arguments are proved in the accompanying
document.  This companion checks every finite polynomial and rational
constant used there, independently recomputes the Jackson coefficients by
convolution, verifies the N=251 height-500 error ledger, and checks that the
adjacent N=250 rational-pi ledger does not clear the same target.

All validity checks remain active under ``python -O``.
"""

from fractions import Fraction
from math import comb, factorial


DANGER_MASS = Fraction(1, 7)
SAFE_MASS = 1 - DANGER_MASS
AMBIENT_FLOOR = Fraction(15, 343)
JACKSON_N = 251
RELATION_HEIGHT = 2 * JACKSON_N - 2
ETA_SIMPLE_CAP = Fraction(21, 12_500)
PI_UPPER_NUMERATOR = 355
PI_UPPER_DENOMINATOR = 113


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def danger_minorant(count: int) -> Fraction:
    """The cubic pointwise minorant of the zero-danger atom."""
    return -Fraction((count - 1) * (count - 3) * (count - 4), 12)


def danger_minorant_binomial(count: int) -> Fraction:
    """The same cubic in the binomial/factorial-moment basis."""
    return (
        1
        - count
        + Fraction(5, 6) * comb(count, 2)
        - Fraction(1, 2) * comb(count, 3)
    )


def jackson_coefficient(n: int, k: int) -> int:
    """Closed integer coefficient C_k of the squared Fejer kernel."""
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
    """Independent convolution of the triangular Fejer coefficients."""
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
    # The cubic is a pointwise minorant on the complete danger-count range.
    for count in range(14):
        require(
            danger_minorant(count) == danger_minorant_binomial(count),
            f"binomial-basis identity failed at count={count}",
        )
        require(
            danger_minorant(count) <= (1 if count == 0 else 0),
            f"pointwise minorant failed at count={count}",
        )

    p = DANGER_MASS
    baseline = sum(
        (
            coefficient * comb(13, degree) * p**degree
            for degree, coefficient in enumerate(
                (Fraction(1), -Fraction(1), Fraction(5, 6), -Fraction(1, 2))
            )
        ),
        Fraction(0),
    )
    require(baseline == Fraction(18, 343), "three-wise baseline changed")

    support_three_floor = baseline - Fraction(1, 2) * (p**2 - p**3)
    require(
        support_three_floor == AMBIENT_FLOOR,
        "unique support-three floor changed",
    )

    support_two_coefficient = Fraction(5, 6) - Fraction(11, 2) * p
    require(
        support_two_coefficient == Fraction(1, 21),
        "unique support-two coefficient changed",
    )
    support_two_floor = baseline - support_two_coefficient * p**2
    require(
        support_two_floor == Fraction(53, 1029),
        "unique support-two floor changed",
    )

    two_line_floors = []
    for support_union_size in range(3, 7):
        floor = (
            (1 - support_union_size * p)
            * SAFE_MASS ** (13 - support_union_size)
        )
        two_line_floors.append(floor)
        require(
            floor >= AMBIENT_FLOOR,
            f"two-line factor floor failed at m={support_union_size}",
        )
    require(
        all(
            left > right
            for left, right in zip(two_line_floors, two_line_floors[1:])
        ),
        "two-line factor floors are not strictly decreasing",
    )
    require(
        two_line_floors[-1] == Fraction(279936, 5764801),
        "m=6 two-line floor changed",
    )
    require(
        two_line_floors[-1] - AMBIENT_FLOOR
        == Fraction(27831, 5764801),
        "m=6 comparison gap changed",
    )

    # Independent Jackson coefficient path at the load-bearing degree.
    for k in range(2 * JACKSON_N - 1):
        closed = jackson_coefficient(JACKSON_N, k)
        convolved = jackson_coefficient_by_convolution(JACKSON_N, k)
        require(closed == convolved, f"Jackson convolution mismatch at k={k}")
        require(closed > 0, f"Jackson coefficient not positive at k={k}")

    c_zero = JACKSON_N * (2 * JACKSON_N * JACKSON_N + 1) // 3
    require(
        c_zero == jackson_coefficient(JACKSON_N, 0),
        "Jackson zero coefficient formula mismatch",
    )
    require(c_zero == 10_542_251, "Jackson C_0 changed")

    eta_251_cap = jackson_eta_pi_cap(JACKSON_N)
    require(eta_251_cap > 0, "Jackson eta cap should be positive")
    require(
        eta_251_cap < ETA_SIMPLE_CAP,
        "N=251 Jackson eta cap no longer clears 21/12500",
    )
    require(
        26 * ETA_SIMPLE_CAP < AMBIENT_FLOOR,
        "height-500 telescope no longer clears ambient floor",
    )
    final_gap = AMBIENT_FLOOR - 26 * ETA_SIMPLE_CAP
    require(
        final_gap == Fraction(111, 2_143_750),
        "final Jackson/Haar gap changed",
    )

    eta_250_cap = jackson_eta_pi_cap(JACKSON_N - 1)
    require(
        eta_250_cap >= AMBIENT_FLOOR / 26,
        "adjacent N=250 rational-pi ledger unexpectedly clears target",
    )

    # Elementary alternative Fejer/log controls retained as a cross-check.
    e_partial = sum(
        (Fraction(1, factorial(k)) for k in range(5)),
        Fraction(0),
    )
    require(e_partial > Fraction(8, 3), "elementary lower bound for e failed")
    require(8**9 > 5649 * 3**9, "elementary log(5649)<9 check failed")
    require(
        Fraction(247, 5649) < AMBIENT_FLOOR,
        "fallback Fejer comparison changed",
    )

    print("THM-2185 exact rank-two safe-cube and rank-three-harvest referee")
    print("danger_minorant=-(N-1)(N-3)(N-4)/12")
    print(f"three_wise_floor={baseline}")
    print(f"unique_support_three_floor={support_three_floor}")
    print(f"unique_support_two_floor={support_two_floor}")
    print(
        "two_sparse_line_floors_m3_to_m6="
        + ",".join(str(value) for value in two_line_floors)
    )
    print(f"universal_ambient_floor={AMBIENT_FLOOR}")
    print(
        f"Jackson_N={JACKSON_N}, relation_height={RELATION_HEIGHT}, "
        f"C_0={c_zero}"
    )
    print(f"eta_{JACKSON_N}_pi_cap<{ETA_SIMPLE_CAP}")
    print(
        f"N={JACKSON_N - 1} rational-pi ledger does not clear "
        f"{AMBIENT_FLOOR / 26}"
    )
    print(f"final_height_500_gap={final_gap}")
    print("deletion_rank_two_height=105000")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
