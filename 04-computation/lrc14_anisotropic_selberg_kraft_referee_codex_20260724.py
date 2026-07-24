#!/usr/bin/env python3
"""Exact referee for THM-2144 (anisotropic Selberg--Kraft boxes).

The analytic input is the signed Selberg tensor minorant of THM-2085.
This companion checks every rational specialization and the two discrete
certificate-optimality boundaries used by THM-2144.
"""

from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def kraft_term(height: int) -> Fraction:
    """The h=1/14 anisotropic Kraft summand."""
    require(height >= 1, "height must be positive")
    return Fraction(14, 6 * height + 13)


def equal_cap_bound(k: int, h: Fraction, height: int) -> Fraction:
    """Constant coefficient for k equal caps."""
    require(k >= 1 and 0 < h < Fraction(1, 2) and height >= 1, "invalid input")
    alpha = 1 - 2 * h
    epsilon = Fraction(1, height + 1)
    u = alpha + epsilon
    return u ** (k - 1) * (alpha - (2 * k - 1) * epsilon)


def anisotropic_factor(heights: tuple[int, ...]) -> Fraction:
    return 1 - sum((kraft_term(height) for height in heights), Fraction(0))


def balance(total: int, count: int) -> tuple[int, ...]:
    q, r = divmod(total, count)
    return (q,) * (count - r) + (q + 1,) * r


def main() -> None:
    # Equal-cap LRC(14) certificate.
    h_lrc14 = Fraction(1, 14)
    b29 = equal_cap_bound(13, h_lrc14, 29)
    b28 = equal_cap_bound(13, h_lrc14, 28)
    u29 = Fraction(6, 7) + Fraction(1, 30)
    final29 = Fraction(6, 7) - 25 * Fraction(1, 30)
    final28 = Fraction(6, 7) - 25 * Fraction(1, 29)
    require(u29 == Fraction(187, 210), "H=29 upper constant")
    require(final29 == Fraction(1, 42), "H=29 final factor")
    require(final28 == Fraction(-1, 203), "H=28 final factor")
    require(b29 == u29**12 * final29 > 0, "H=29 positivity")
    require(b28 < 0, "H=28 must fail")

    # The relation-free endpoint forced by H=29 is 1/12, obtained as a
    # supremum because the tensor bound itself vanishes at the endpoint.
    require(equal_cap_bound(13, Fraction(1, 12), 29) == 0, "endpoint vanishing")
    require(
        equal_cap_bound(13, Fraction(1, 12) - Fraction(1, 10**6), 29) > 0,
        "strictly sub-endpoint positivity",
    )

    # Total coordinate budget: discrete balancing minimizes the convex Kraft
    # sum.  The exact adjacent boundary is 366 (fail) versus 367 (pass).
    profile366 = balance(366, 13)
    profile367 = balance(367, 13)
    require(profile366 == (28,) * 11 + (29,) * 2, "budget-366 balance")
    require(profile367 == (28,) * 10 + (29,) * 3, "budget-367 balance")
    sum366 = sum((kraft_term(height) for height in profile366), Fraction(0))
    sum367 = sum((kraft_term(height) for height in profile367), Fraction(0))
    require(sum366 == Fraction(33866, 33847) > 1, "budget 366 must fail")
    require(sum367 == Fraction(33782, 33847) < 1, "budget 367 must pass")

    # Exact exchange check behind the balancing argument on the whole relevant
    # range: moving one unit from a high cap to a cap at least two lower
    # strictly decreases the Kraft sum.
    for low in range(1, 368):
        for high in range(low + 2, 369):
            before = kraft_term(low) + kraft_term(high)
            after = kraft_term(low + 1) + kraft_term(high - 1)
            require(after < before, f"balancing failed at low={low}, high={high}")

    # One coordinate of cap 1 and twelve complementary caps.  H=105 is the
    # first uniform complement accepted by this profile.
    rank105 = kraft_term(1) + 12 * kraft_term(105)
    rank104 = kraft_term(1) + 12 * kraft_term(104)
    require(rank105 == Fraction(12194, 12217) < 1, "rank-105 profile")
    require(rank104 == Fraction(1730, 1729) > 1, "rank-104 boundary")

    print("THM-2144 exact anisotropic Selberg--Kraft referee")
    print(f"H=29: u={u29}, final_factor={final29}, B={b29} > 0")
    print(f"H=28: final_factor={final28}, B={b28} < 0")
    print("H=29 endpoint formula: B(h=1/12)=0 and its final factor is positive exactly for h<1/12")
    print(f"budget 366 profile={profile366}: Kraft sum={sum366} > 1")
    print(f"budget 367 profile={profile367}: Kraft sum={sum367} < 1")
    print("balancing exchange: strict decrease verified for 1<=low and low+2<=high<=368")
    print(f"one cap 1 + twelve caps 104: Kraft sum={rank104} > 1")
    print(f"one cap 1 + twelve caps 105: Kraft sum={rank105} < 1")
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
