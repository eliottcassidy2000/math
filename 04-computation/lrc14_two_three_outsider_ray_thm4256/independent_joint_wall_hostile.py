#!/usr/bin/env python3
"""Exact full-joint-wall control for the robust-boundary body at g=1015.

This deliberately does not use the endpoint primitive or repair transform.
It constructs the wall arrangement of the literal 24-label constrained set
and sums its safe cells on one common (arbitrary-precision) denominator.
"""

from fractions import Fraction
from math import gcd

POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
D = 18_241_159_416_480
G = 1015
BODY = (88, 95, 170, 193, 240, 252, 264, 286, 290)
REPAIR = (10, 40, 60, 80, 85, 143, 145, 176)
EXPECTED_COCYCLE_MARGIN = 301_703_072_020_224_840


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def main() -> None:
    if set(BODY) & set(REPAIR):
        raise RuntimeError("hostile body and exact witness repair intersect")
    speeds = tuple(x for x in POOL if x not in REPAIR) + (2 * G, 3 * G)
    denominator = 1
    for speed in speeds:
        denominator = lcm(denominator, 14 * speed)

    walls = {0, denominator}
    for speed in speeds:
        quantum = denominator // (14 * speed)
        for tooth in range(speed):
            walls.add((14 * tooth + 1) * quantum)
            walls.add((14 * tooth + 13) * quantum)
    ordered = sorted(walls)

    safe_ticks = 0
    safe_cells = 0
    for left, right in zip(ordered, ordered[1:]):
        midpoint_twice = left + right
        safe = True
        for speed in speeds:
            residue = speed * midpoint_twice % (2 * denominator)
            if not (7 * residue >= denominator and
                    7 * residue <= 13 * denominator):
                safe = False
                break
        if safe:
            safe_ticks += right - left
            safe_cells += 1

    literal_margin = 63 * safe_ticks - 4 * denominator
    if literal_margin < 0:
        raise RuntimeError("literal hostile witness falls below 4/63")
    scaled = Fraction(literal_margin, denominator) * (84 * D * G)
    if scaled.denominator != 1:
        raise RuntimeError("cocycle normalization is unexpectedly nonintegral")
    if scaled.numerator != EXPECTED_COCYCLE_MARGIN:
        raise RuntimeError("literal/cocycle margin mismatch")

    mass = Fraction(safe_ticks, denominator)
    excess = mass - Fraction(4, 63)
    print("LRC_23_RAY_FULL_JOINT_WALL_HOSTILE_V1")
    print("G", G, "LABELS", ",".join(map(str, speeds)))
    print("WALLS", len(ordered), "OPEN_CELLS", len(ordered) - 1,
          "SAFE_CELLS", safe_cells)
    print("MASS", f"{mass.numerator}/{mass.denominator}")
    print("EXCESS_OVER_4_63", f"{excess.numerator}/{excess.denominator}")
    print("COCYCLE_MARGIN_NUM", scaled.numerator)
    print("VERDICT ROBUST_1015_BODY_LITERAL_EXACT_WITNESS_PASS")


if __name__ == "__main__":
    main()
