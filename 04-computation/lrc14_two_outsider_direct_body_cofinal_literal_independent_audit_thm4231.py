#!/usr/bin/env python3
"""Standalone literal hostile/boundary audit for the Q=1307 body witness."""

from fractions import Fraction
from math import gcd

B_STAR = (170, 176, 190, 193, 240, 252, 264, 286, 290)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    denominator = 1
    for speed in speeds:
        denominator = lcm(denominator, 14 * speed)
    walls = {0, denominator}
    for speed in speeds:
        unit = denominator // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * unit)
            walls.add((14 * tooth + 13) * unit)
    ordered = sorted(walls)
    period_twice = 2 * denominator
    mass = 0
    components = 0
    previous_safe = False
    for left, right in zip(ordered, ordered[1:]):
        midpoint_twice = left + right
        safe = all(
            7 * (speed * midpoint_twice % period_twice) >= denominator
            and 7 * (speed * midpoint_twice % period_twice) <= 13 * denominator
            for speed in speeds
        )
        if safe:
            mass += right - left
            components += not previous_safe
        previous_safe = safe
    return denominator, mass, components


def reduced(numerator: int, denominator: int) -> str:
    value = Fraction(numerator, denominator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    denominator, mass, components = geometry(B_STAR)
    surplus = 45 * mass - 4 * denominator
    numerator = 108 * components * denominator
    divisor = 7 * surplus
    activation = (numerator + divisor - 1) // divisor
    assert (denominator, mass, components, surplus, activation) == (
        18_241_159_416_480,
        4_579_301_272_924,
        618,
        133_103_919_615_660,
        1307,
    )
    left_slack = numerator - divisor * (activation - 1)
    right_slack = divisor * activation - numerator
    assert left_slack == 651_910_967_177_400
    assert right_slack == 279_816_470_132_220

    print("LRC14_COMPOSITE_DESCENT_LITERAL_REFEREE")
    print(
        f"BASE B {B_STAR} D {denominator} MASS {mass} COMPONENTS {components} "
        f"SURPLUS45 {surplus} KAPPA {activation} LEFT_SLACK {left_slack} "
        f"RIGHT_SLACK {right_slack}"
    )
    for pair in ((1306, 1307), (1307, 1308)):
        speeds = tuple(sorted(B_STAR + pair))
        pair_denominator, pair_mass, pair_components = geometry(speeds)
        delta = 63 * pair_mass - 4 * pair_denominator
        assert delta > 0
        print(
            f"PAIR {pair} D {pair_denominator} MASS {pair_mass} "
            f"COMPONENTS {pair_components} DELTA63 {delta} "
            f"MU {reduced(pair_mass, pair_denominator)} "
            f"DELTA {reduced(delta, pair_denominator)}"
        )
    print("VERDICT BASE_CEILING_SHARP LITERAL_ADJACENT_CONTROLS_STRICTLY_SAFE")


if __name__ == "__main__":
    main()
