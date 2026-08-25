#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4089.

This is independent of the external interval certificate.  It checks the
first-chamber stationary point, the p=11 boundary, the complete piecewise
minimization of the idealized p=13,s=3 arithmetic cost, and a rational
exponential-series bound proving that the full margin is negative.

Reproduction:
  python -B 04-computation/hybrid_padic_margin_boundary_thm4089.py
  python -B -O 04-computation/hybrid_padic_margin_boundary_thm4089.py
"""

from fractions import Fraction
from math import factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def xi_star(p: int, s: int) -> Fraction:
    return Fraction(
        12 * (s * (s + 1) - 1),
        12 * s * s + (s - 1) * (p + 1),
    )


def q_chamber_one(x: Fraction) -> Fraction:
    return Fraction(9, 2) * x * x - 11 * x + Fraction(95, 4)


def q_chamber_two(x: Fraction) -> Fraction:
    return x + Fraction(63, 4)


def q_chamber_three(x: Fraction) -> Fraction:
    return 3 * x * x - 8 * x + Fraction(45, 2)


def q_chamber_four(x: Fraction) -> Fraction:
    return 4 * x + Fraction(21, 2)


def q_chamber_five(x: Fraction) -> Fraction:
    return 7 * x + Fraction(3, 2) * (4 - x) ** 2


def main() -> None:
    print("THM-4089 HYBRID P-ADIC MARGIN BOUNDARY AUDIT")
    print()
    print("first-chamber stationary points")
    for p in (2, 3, 5, 7, 11, 13):
        x = xi_star(p, 3)
        sign = "interior" if x > 1 else ("boundary" if x == 1 else "outside")
        require((x > 1) == (p < 11), "p=11 stationary-point boundary failed")
        print(f"p={p:2d}, s=3: xi*={x}, position={sign}")

    # Exact identity behind the general p=11 boundary, tested on a broad box.
    identity_gates = 0
    for p in range(2, 40):
        for s in range(2, 30):
            x = xi_star(p, s)
            denominator = 12 * s * s + (s - 1) * (p + 1)
            require(
                x - 1 == Fraction((s - 1) * (11 - p), denominator),
                "xi*-1 identity failed",
            )
            identity_gates += 1

    # For p=13,s=3 with the most favorable possible one-power Hasse cost,
    # Q(xi)=S_0(xi)+3I^4_xi(xi) has five chambers on 1<xi<4.
    candidates = (
        ("1<xi<4/3", Fraction(11, 9), q_chamber_one(Fraction(11, 9))),
        ("4/3<=xi<3/2", Fraction(4, 3), q_chamber_two(Fraction(4, 3))),
        ("3/2<=xi<2", Fraction(3, 2), q_chamber_three(Fraction(3, 2))),
        ("2<=xi<3", Fraction(2), q_chamber_four(Fraction(2))),
        ("3<=xi<4", Fraction(3), q_chamber_five(Fraction(3))),
    )
    best = min(candidates, key=lambda row: row[2])
    require(best == candidates[0], "ideal p=13 cost has the wrong chamber minimum")
    require(best[1] == Fraction(11, 9), "ideal p=13 optimizer is not 11/9")
    require(best[2] == Fraction(613, 36), "ideal p=13 cost minimum is not 613/36")
    tau_min = Fraction(1, 8) * best[2]
    require(tau_min == Fraction(613, 288), "ideal p=13 tau minimum is not 613/288")

    # For the actual p=13 integrand, c_13=7/6 and xi>1 lies beyond its
    # support, so J_13(xi)=1/(2c_13)=3/7.  The five cost derivatives are
    # 9xi-9, 3, 6xi-6, 6, 3xi-3: all positive on their chambers.
    actual_j = Fraction(3, 7)
    actual_q_inf = (9 - 2 * actual_j) + 3 * Fraction(41, 12)
    actual_tau_inf = Fraction(1, 8) * actual_q_inf
    require(actual_q_inf == Fraction(515, 28), "actual p=13 cost infimum is wrong")
    require(actual_tau_inf == Fraction(515, 224), "actual p=13 tau infimum is wrong")

    print()
    print("ideal p=13,s=3 one-power Hasse chamber minima")
    for chamber, point, value in candidates:
        print(f"{chamber:16s}: point={point:>4}, Q={value}, tau={value / 8}")

    # The positive exponential series proves exp(8/3)>13, hence log(13)<8/3.
    exp_partial = sum(
        (Fraction(8, 3) ** n / factorial(n) for n in range(6)),
        Fraction(0),
    )
    require(exp_partial == Fraction(49621, 3645), "unexpected exponential partial sum")
    require(exp_partial > 13, "partial sum does not prove exp(8/3)>13")
    margin_upper = 8 - 4 * tau_min
    require(margin_upper == Fraction(-37, 72), "unexpected global margin upper bound")
    require(margin_upper < 0, "ideal p=13 margin was not excluded")
    actual_margin_upper = 8 - 4 * actual_tau_inf
    require(actual_margin_upper == Fraction(-67, 56), "unexpected actual margin upper bound")

    # At the exact one-layer analytic stationary point x=13Y, the derivative
    # equation is sqrt(1-x^2)=13/16 and x=sqrt(87)/16.
    require(Fraction(87, 256) + Fraction(169, 256) == 1, "analytic stationary identity failed")

    print()
    print(f"general xi*-1 identity gates = {identity_gates}")
    print(f"exp(8/3) partial sum n<=5   = {exp_partial} > 13")
    print(f"ideal arithmetic tau floor  = {tau_min}")
    print(f"all-Y margin upper bound     = {margin_upper} < 0")
    print(f"actual J_13(xi), xi>1        = {actual_j}")
    print(f"actual tau infimum, xi->1+   = {actual_tau_inf}")
    print(f"actual all-Y margin bound    = {actual_margin_upper} < 0")
    print("analytic stationary point    = 13Y=sqrt(87)/16, sqrt(1-(13Y)^2)=13/16")
    print("EXACT NEXT-CASE OBSTRUCTION: PASS")


if __name__ == "__main__":
    main()
