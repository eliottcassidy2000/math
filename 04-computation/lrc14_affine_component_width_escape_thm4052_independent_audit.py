#!/usr/bin/env python3
"""No-import circle-wall audit for THM-4052.

Unlike the primary certificate, this path never uses the THM-4041 component
formula.  It constructs every strict danger wall on the pack-phase circle,
tests the intervening cells literally, and recovers component lengths.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations


GATES = 0


def check(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)


def norm_distance(q):
    residue = q - q.numerator // q.denominator
    return min(residue, 1 - residue)


def spoiled(d, speeds, y):
    for label in range(d):
        x = (y + label) / d
        if not any(norm_distance(speed * x) < Fraction(1, 14)
                   for speed in speeds):
            return False
    return True


def circle_walls(d, speeds):
    walls = {Fraction(0)}
    for speed in speeds:
        for label in range(d):
            for integer in range(-1, speed + 2):
                for sign in (-1, 1):
                    y = Fraction(d, speed) * (integer + Fraction(sign, 14)) - label
                    if 0 < y < 1:
                        walls.add(y)
    return tuple(sorted(walls))


def component_lengths(d, speeds):
    walls = circle_walls(d, speeds)
    endpoints = walls + (Fraction(1),)
    lengths = tuple(endpoints[i + 1] - endpoints[i]
                    for i in range(len(endpoints) - 1))
    flags = tuple(spoiled(d, speeds, (endpoints[i] + endpoints[i + 1]) / 2)
                  for i in range(len(lengths)))

    runs = []
    current = Fraction(0)
    for flag, length in zip(flags, lengths):
        if flag:
            current += length
        elif current:
            runs.append(current)
            current = Fraction(0)
    if current:
        runs.append(current)

    if runs and len(runs) > 1 and flags[0] and flags[-1]:
        runs[0] += runs[-1]
        runs.pop()
    if flags and all(flags):
        runs = [Fraction(1)]
    return tuple(sorted(runs)), len(lengths)


def ftext(q):
    return f"{q.numerator}/{q.denominator}"


def direct_clearance(speeds, x):
    return min(norm_distance(speed * x) for speed in speeds)


def main():
    # Complete literal reconstruction of the d2 component multiset.
    semantic = []
    d2_cells = 0
    for alpha in range(1, 80, 2):
        for beta in range(alpha + 2, 80, 2):
            lengths, cells = component_lengths(2, (alpha, beta))
            d2_cells += cells
            semantic.append(
                f"{alpha},{beta}:" + ",".join(ftext(length) for length in lengths)
            )
    check(len(semantic) == 780, "d2 audit universe changed")
    digest = sha256("\n".join(semantic).encode()).hexdigest()

    # The capacity-equality partitions force every component into one
    # assignment stratum.  The literal wall engine verifies the resulting
    # d/(7E) width bound over the exact inherited hostile boxes.
    d3_profiles = 0
    d3_cells = 0
    d3_ratio = Fraction(0)
    units3 = tuple(value for value in range(1, 24) if value % 3)
    for exceptions in combinations(units3, 3):
        lengths, cells = component_lengths(3, exceptions)
        bound = Fraction(3, 7 * max(exceptions))
        check(all(length <= bound for length in lengths),
              f"d3 width bound failed at {exceptions}")
        if lengths:
            d3_ratio = max(d3_ratio, max(lengths) / bound)
        d3_profiles += 1
        d3_cells += cells
    check(d3_profiles == 560, "d3 audit universe changed")

    d4_profiles = 0
    d4_cells = 0
    d4_ratio = Fraction(0)
    odds = tuple(range(1, 20, 2))
    for r in range(1, 12, 2):
        for first, second in combinations(odds, 2):
            exceptions = (2 * r, first, second)
            lengths, cells = component_lengths(4, exceptions)
            bound = Fraction(4, 7 * max(exceptions))
            check(all(length <= bound for length in lengths),
                  f"d4 width bound failed at {exceptions}")
            if lengths:
                d4_ratio = max(d4_ratio, max(lengths) / bound)
            d4_profiles += 1
            d4_cells += cells
    check(d4_profiles == 270, "d4 audit universe changed")

    # Rebuild the three endpoint controls directly from the full row.
    positive = (
        (2, tuple(range(1, 12)), (1, 133), Fraction(1, 12),
         Fraction(13, 154), Fraction(167, 308), Fraction(1, 14)),
        (3, tuple(range(1, 11)), (1, 110, 23), Fraction(1, 11),
         Fraction(13, 140), Fraction(51, 140), Fraction(1, 14)),
        (4, tuple(range(1, 11)), (2, 185, 11), Fraction(1, 11),
         Fraction(137, 1540), Fraction(4757, 6160), Fraction(137, 1540)),
    )
    for d, pack, exceptions, central, escape, x, expected in positive:
        check(spoiled(d, exceptions, central), "positive central phase changed")
        check(direct_clearance(tuple(d * h for h in pack) + exceptions, x) == expected,
              "positive full-row clearance changed")
        check((d * x - escape).denominator == 1,
              "positive escape left the lift bank")

    # Fixed two-phase endpoint banks are not all-height certificates in d3/d4.
    d3_hostile = (1, 55, 56)
    check(spoiled(3, d3_hostile, Fraction(1, 14)), "d3 bank hostile phase one")
    check(spoiled(3, d3_hostile, Fraction(13, 140)), "d3 bank hostile phase two")
    check(direct_clearance(tuple(3 * h for h in range(1, 11)) + d3_hostile,
                           Fraction(1, 13)) == Fraction(1, 13),
          "d3 hostile is not the declared positive row")

    d4_hostile = (26, 31, 57)
    check(spoiled(4, d4_hostile, Fraction(15, 98)), "d4 bank hostile phase one")
    check(spoiled(4, d4_hostile, Fraction(1, 14)), "d4 bank hostile phase two")
    check(direct_clearance(tuple(4 * h for h in range(1, 11)) + d4_hostile,
                           Fraction(1, 11)) == Fraction(1, 11),
          "d4 hostile is not the declared positive row")

    print("LRC14 AFFINE COMPONENT-WIDTH INDEPENDENT CIRCLE AUDIT")
    print("method=literal_strict_circle_walls+midcell_masks;no_formula_import")
    print(f"d2_pairs=780;cells={d2_cells};semantic_digest={digest}")
    print(f"d3_profiles={d3_profiles};cells={d3_cells};max_width_ratio={ftext(d3_ratio)}")
    print(f"d4_profiles={d4_profiles};cells={d4_cells};max_width_ratio={ftext(d4_ratio)}")
    print("endpoint_controls=d2,d3,d4:PASS")
    print("fixed_bank_hostiles=d3:(1,55,56),period420;d4:(26,31,57),period392")
    print(f"gates={GATES}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
