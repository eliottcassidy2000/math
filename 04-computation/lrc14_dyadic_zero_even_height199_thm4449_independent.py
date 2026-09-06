#!/usr/bin/env python3
"""Independent interval-union height-199 census for THM-4449.

Unlike the main independent audit's wall/midpoint method, this constructs the
two parity-owner tooth lists, intersects opposite colours pairwise, and unions
the resulting open intervals on a common exact integer lattice.
"""

import sys
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def allowed(n):
    return n > 0 and n % 2 == 1 and n % 3 != 0


def owner_intervals(speed, parity, common):
    """Clipped open quotient intervals owned by one tail and one lift."""
    scale = common // speed
    ceiling = 7 * common
    result = []
    for nearest in range(parity, speed + 1, 2):
        left = max(0, (7 * nearest - 1) * scale)
        right = min(ceiling, (7 * nearest + 1) * scale)
        if left < right:
            result.append((left, right))
    return result


def intersect_sorted(first, second):
    i = j = 0
    out = []
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return out


def quotient_stats(triple):
    common = lcm(*triple)
    denominator = 7 * common
    owners = {(s, bit): owner_intervals(s, bit, common) for s in triple for bit in (0, 1)}
    pieces = []
    for a, b in combinations(triple, 2):
        pieces.extend(intersect_sorted(owners[a, 0], owners[b, 1]))
        pieces.extend(intersect_sorted(owners[a, 1], owners[b, 0]))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left < merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    mass = F(sum(right - left for left, right in merged), denominator)
    longest = F(max((right - left for left, right in merged), default=0), denominator)
    return mass, longest, len(merged)


def main():
    height = int(sys.argv[1]) if len(sys.argv) > 1 else 199
    values = [n for n in range(1, height + 1) if allowed(n)]
    maximum = F(-1)
    winners = []
    histogram = {}
    checked = 0
    for triple in combinations(values, 3):
        mass, width, count = quotient_stats(triple)
        checked += 1
        histogram[mass] = histogram.get(mass, 0) + 1
        record = (triple, width, count)
        if mass > maximum:
            maximum = mass
            winners = [record]
        elif mass == maximum:
            winners.append(record)
    require(checked == 47905 if height == 199 else checked >= 0, (checked, "triple count"))
    if height == 199:
        expected = [(t, 7 * t, 11 * t) for t in (1, 5, 7, 11, 13, 17)]
        require(maximum == F(72, 539), (maximum, "height-199 maximum"))
        require([triple for triple, _, _ in winners] == expected, (winners, "height-199 leaders"))
        # These are quotient-y components.  The physical pullback under
        # doubling has twice as many components, each with half the length.
        require(winners[0][1:] == (F(2, 77), 6), (winners[0], "primitive quotient geometry"))
        require(
            all(width == F(2, 77 * t) and count == 6 * t for (triple, width, count), t in zip(winners, (1, 5, 7, 11, 13, 17))),
            (winners, "dilated quotient geometry"),
        )
    print("THM4449_INDEPENDENT_HEIGHT199_CENSUS")
    print(f"height={height} labels={len(values)} triples={checked} distinct_masses={len(histogram)}")
    print(f"maximum={maximum}")
    for triple, width, count in winners:
        print(
            f"winner={triple} primitive_gcd={gcd(gcd(*triple[:2]), triple[2])} "
            f"quotient_longest={width} quotient_components={count} "
            f"physical_longest={width/2} physical_components={2*count}"
        )
    print("PASS")


if __name__ == "__main__":
    main()
