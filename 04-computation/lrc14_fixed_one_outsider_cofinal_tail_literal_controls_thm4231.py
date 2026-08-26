#!/usr/bin/env python3
"""Exact endpoint-toggle controls for THM-4231's fixed q=1 cofinal tail."""

from fractions import Fraction
from math import gcd


BASE = (1, 170, 176, 190, 193, 240, 252, 264, 286, 290)
D = 18_241_159_416_480
MASS_TICKS = 3_885_436_686_322
COMPONENTS = 528


def lcm(left: int, right: int) -> int:
    return left // gcd(left, right) * right


def geometry(speeds: tuple[int, ...]) -> tuple[int, int, int]:
    common = 1
    for speed in speeds:
        common = lcm(common, 14 * speed)
    delta: dict[int, int] = {}
    for speed in speeds:
        unit = common // (14 * speed)
        for tooth in range(speed):
            leave = (14 * tooth + 1) * unit
            enter = (14 * tooth + 13) * unit
            delta[leave] = delta.get(leave, 0) - 1
            delta[enter] = delta.get(enter, 0) + 1
    failed = len(speeds)
    previous = 0
    previous_safe = False
    mass = 0
    components = 0
    for tick in sorted(delta):
        safe = failed == 0
        if safe:
            mass += tick - previous
            if not previous_safe:
                components += 1
        previous_safe = safe
        failed += delta[tick]
        previous = tick
    safe = failed == 0
    if safe:
        mass += common - previous
        if not previous_safe:
            components += 1
    assert failed == len(speeds)
    return common, mass, components


def show(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    mass = Fraction(MASS_TICKS, D)
    surplus = 54 * MASS_TICKS - 4 * D
    rows = []
    for r in (542, 543):
        gap = Fraction(surplus, 63 * D) - Fraction(6 * COMPONENTS, 49 * r)
        integer_margin = 7 * surplus * r - 54 * COMPONENTS * D
        common, ticks, components = geometry(BASE + (r,))
        literal_mass = Fraction(ticks, common)
        literal_delta = 63 * literal_mass - 4
        assert literal_delta > 0
        rows.append(
            (r, gap, integer_margin, components, literal_mass, literal_delta)
        )
    assert rows[0][1] < 0 < rows[1][1]

    print("LRC14_FIXED_Q1_LITERAL_CONTROLS")
    print(
        f"BASE MASS {show(mass)} COMPONENTS {COMPONENTS} "
        f"SURPLUS54 {surplus} TAIL_K 543"
    )
    for r, gap, margin, components, literal_mass, literal_delta in rows:
        print(
            f"R {r} BOUND_GAP {show(gap)} INTEGER_MARGIN {margin} "
            f"LITERAL_COMPONENTS {components} LITERAL_MASS {show(literal_mass)} "
            f"LITERAL_DELTA63 {show(literal_delta)} POSITIVE YES"
        )
    print("VERDICT CERTIFICATE_BOUNDARY_NOT_LITERAL_COUNTEREXAMPLE")


if __name__ == "__main__":
    main()
