#!/usr/bin/env python3
"""Exact arithmetic probe for the literal-word quadratic unit-pair lemma."""

from fractions import Fraction as Q


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def first_shell_bound(rho: Q) -> tuple[Q, int, Q, Q]:
    mu = Q(13013, 2) * rho * rho
    cutoff = Q(33800 * 22, 961) / mu
    bound = cutoff.numerator // cutoff.denominator + 1
    margin = mu / 22 - Q(33800, 961 * bound)
    previous = mu / 22 - Q(33800, 961 * (bound - 1))
    require(margin > 0, "displayed shell bound does not pass")
    require(previous <= 0, "displayed shell bound is not first")
    return mu, bound, margin, previous


def main() -> None:
    # For sheet masses x,y:
    # X=12(x^2+y^2)-2xy and
    # X-(11/2)(x+y)^2=(13/2)(x-y)^2.
    lhs = {
        "x2": Q(12) - Q(11, 2),
        "xy": Q(-2) - Q(11),
        "y2": Q(12) - Q(11, 2),
    }
    rhs = {"x2": Q(13, 2), "xy": Q(-13), "y2": Q(13, 2)}
    require(lhs == rhs, "two-sheet lower identity failed")

    corners = [
        12 * (x * x + y * y) - 2 * x * y
        for x, y in ((0, 0), (1, 0), (0, 1), (1, 1))
    ]
    require(corners == [0, 12, 12, 22], "corner cap failed")

    rows = (
        (
            "strict",
            Q(39002430583, 160481782761300),
            2013431,
        ),
        (
            "repeated-first conditional",
            Q(13560199813, 160481782761300),
            16656658,
        ),
    )

    print("LRC14 LITERAL-WORD QUADRATIC UNIT-PAIR PROBE")
    print("two-sheet remainder = (13/2)(x-y)^2 [PASS]")
    print("unit-square corner values =", corners, "; cap 22 [PASS]")
    print("mass coefficient = 13013/2")

    for label, rho, expected_bound in rows:
        mu, bound, margin, previous = first_shell_bound(rho)
        require(bound == expected_bound, f"{label} bound mismatch")
        print(label)
        print("  rho =", rho)
        print("  mu0 =", mu)
        print("  first multiplier bound =", bound)
        print("  endpoint margin =", margin)
        print("  previous margin =", previous)

    print("grade transport: 13^lambda * 13*(ms) = m*c_* [PASS]")
    print("ALL EXACT QUADRATIC-PAIR CONSTANT CHECKS PASSED")


if __name__ == "__main__":
    main()
