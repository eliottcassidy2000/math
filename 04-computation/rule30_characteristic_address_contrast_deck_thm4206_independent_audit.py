#!/usr/bin/env python3
"""Independent scalar audit for THM-4206's Rule-30 finite controls.

This path imports nothing from the bit-parallel primary implementation.  It
enumerates ordinary input words, evolves scalar rows, and checks depths one
through eight.
"""

from fractions import Fraction


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def step(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def hostile_census() -> None:
    triple = [0] * 8
    pairs = {"AB": [0] * 4, "AC": [0] * 4, "BC": [0] * 4}
    triple_walsh = 0
    deck = set()
    conditional_support = {0: set(), 1: set()}

    for assignment in range(8):
        x0 = assignment & 1
        x1 = (assignment >> 1) & 1
        x2 = (assignment >> 2) & 1
        a = x1
        b = x0
        c = step(x0, x1, x2)
        code = (a << 2) | (b << 1) | c
        triple[code] += 1
        pairs["AB"][(a << 1) | b] += 1
        pairs["AC"][(a << 1) | c] += 1
        pairs["BC"][(b << 1) | c] += 1
        triple_walsh += (-1) ** (a ^ b ^ c)
        deck.add((x2, a, b ^ c))
        conditional_support[x2].add(code)

    require(triple == [1, 1, 1, 1, 0, 2, 2, 0], "direct triple")
    require(pairs["AB"] == pairs["AC"] == [2, 2, 2, 2], "direct cross pairs")
    require(pairs["BC"] == [1, 3, 3, 1], "direct collision pair")
    require(triple_walsh == 4, "direct third Walsh")
    require(deck == {(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 1)},
            "direct deck")
    require(tuple(len(conditional_support[x2]) for x2 in (0, 1)) == (4, 4),
            "direct conditional supports")
    print(
        "independent_hostile=triple:(1,1,1,1,0,2,2,0);"
        "pairs:(uniform,uniform,(1,3,3,1));third_walsh:1/2;"
        "deck:D=A_or_x2;conditional_supports:(4,4)"
    )


def direct_count(depth: int) -> tuple[int, int]:
    variables = 2 * depth
    ones = 0
    for assignment in range(1 << variables):
        row = [0] + [(assignment >> j) & 1 for j in range(variables)]
        for _ in range(depth):
            row = [step(row[j], row[j + 1], row[j + 2])
                   for j in range(len(row) - 2)]
        require(len(row) == 1, "scalar cone closes")
        ones += row[0]
    return ones, 1 << variables


def slope_one_census() -> None:
    expected = (
        Fraction(-1, 2), Fraction(1, 4), Fraction(-1, 4), Fraction(5, 32),
        Fraction(-5, 64), Fraction(77, 1024), Fraction(-141, 2048),
        Fraction(39, 512),
    )
    observed = []
    counts = []
    for depth in range(1, 9):
        ones, total = direct_count(depth)
        observed.append(Fraction(total - 2 * ones, total))
        counts.append((ones, total))
    require(tuple(observed) == expected, "independent slope-one bank")
    print("independent_slope_one_ones=" + ";".join(f"{a}/{b}" for a, b in counts))
    print("independent_slope_one_correlations=" +
          ",".join(str(value) for value in observed))
    print("finite_scope=depths:1..8;no_extrapolation:true")


def main() -> None:
    hostile_census()
    slope_one_census()
    print("THM4206_INDEPENDENT_PASS")


if __name__ == "__main__":
    main()
