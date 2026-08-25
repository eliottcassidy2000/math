#!/usr/bin/env python3
"""Optimization-safe finite controls for the THM-4028 average-order audit.

The asymptotic itself is proved in REPORT.md by an anisotropic Riemann-sum
sandwich.  This companion checks the constant at high precision, compares two
independent exact cumulative counters, and replays the 2-4-6-9 supercritical
hole using a triangular-value set rather than a floating or integer square
test.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial

import mpmath as mp


DEGREES = (2, 4, 6, 8)
LOWER_BOUNDS = (2, 3, 5, 7)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def indexed_values(degree: int, lower: int, limit: int) -> list[tuple[int, int]]:
    values: list[tuple[int, int]] = []
    top = lower
    while comb(top, degree) <= limit:
        values.append((comb(top, degree), top))
        top += 1
    return values


def triangular_count(limit: int) -> int:
    """Count w>=2 with C(w,2)<=limit, by monotone binary search."""
    if limit < 1:
        return 0
    low, high = 2, 2
    while comb(high, 2) <= limit:
        high *= 2
    while low < high:
        middle = (low + high) // 2
        if comb(middle, 2) <= limit:
            low = middle + 1
        else:
            high = middle
    return low - 2


def cumulative_nested(limit: int) -> int:
    """Sum the exact triangular count over the 4/6/8 coordinates."""
    higher = [indexed_values(r, lower, limit) for r, lower in zip(DEGREES[1:], LOWER_BOUNDS[1:])]
    total = 0
    for value8, _ in higher[2]:
        for value6, _ in higher[1]:
            if value8 + value6 >= limit:
                break
            for value4, _ in higher[0]:
                remainder = limit - value8 - value6 - value4
                if remainder < 1:
                    break
                total += triangular_count(remainder)
    return total


def cumulative_meet_in_middle(limit: int) -> int:
    """Independent pair-sum/bisect count of all four-coordinate tuples."""
    atoms = [
        indexed_values(r, lower, limit)
        for r, lower in zip(DEGREES, LOWER_BOUNDS)
    ]
    left = sorted(a + b for a, _ in atoms[0] for b, _ in atoms[1] if a + b <= limit)
    right = sorted(c + d for c, _ in atoms[2] for d, _ in atoms[3] if c + d <= limit)
    return sum(bisect_right(right, limit - value) for value in left)


def mixed_2469_certificate(target: int) -> tuple[int, list[tuple[int, int, int, int]]]:
    degrees = (4, 6, 9)
    lowers = (3, 5, 8)
    higher = [indexed_values(r, lower, target) for r, lower in zip(degrees, lowers)]
    triangular: dict[int, int] = {0: 0}  # stronger hostile: even allow T=0
    w = 2
    while comb(w, 2) <= target:
        triangular[comb(w, 2)] = w
        w += 1
    triples = 0
    representations: list[tuple[int, int, int, int]] = []
    for value9, z in higher[2]:
        for value6, y in higher[1]:
            if value9 + value6 > target:
                break
            for value4, x in higher[0]:
                total = value9 + value6 + value4
                if total > target:
                    break
                triples += 1
                remainder = target - total
                if remainder in triangular:
                    representations.append((triangular[remainder], x, y, z))
    return triples, representations


def first_positive_2469(target: int) -> tuple[int, int, int, int] | None:
    _, reps = mixed_2469_certificate(target)
    return next((rep for rep in reps if rep[0] >= 2), None)


def main() -> None:
    reciprocal_sum = sum((Fraction(1, r) for r in DEGREES), Fraction(0))
    require(reciprocal_sum == Fraction(25, 24), "reciprocal-degree sum")
    mp.mp.dps = 70
    s_mp = mp.mpf(reciprocal_sum.numerator) / reciprocal_sum.denominator
    volume = mp.fprod(
        mp.gamma(1 + mp.mpf(1) / r) * mp.mpf(factorial(r)) ** (mp.mpf(1) / r)
        for r in DEGREES
    ) / mp.gamma(1 + s_mp)
    shell = s_mp * volume
    expected = mp.mpf("24.31102486226095468290593868402939357063510418733886655384259744")
    require(abs(volume - expected) < mp.mpf("1e-62"), "high-precision volume constant")
    print(f"RECIPROCAL_SUM={reciprocal_sum}")
    print(f"CUMULATIVE_CONSTANT={mp.nstr(volume, 66)}")
    print(f"FORMAL_SHELL_CONSTANT={mp.nstr(shell, 66)}")

    semantic: list[str] = []
    for limit in (10, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000):
        nested = cumulative_nested(limit)
        meet = cumulative_meet_in_middle(limit)
        require(nested == meet, f"counter mismatch at X={limit}")
        ratio = mp.mpf(nested) / mp.power(limit, s_mp)
        row = f"X={limit} A={nested} ratio={mp.nstr(ratio, 18)}"
        semantic.append(row)
        print(row)

    hole = 1_061_619
    triples, representations = mixed_2469_certificate(hole)
    require(triples == 28_037, "2-4-6-9 triple universe")
    require(representations == [], "2-4-6-9 hostile is represented")
    left = first_positive_2469(hole - 1)
    right = first_positive_2469(hole + 1)
    require(left == (1447, 18, 17, 9), f"left control {left}")
    require(right == (1456, 17, 5, 8), f"right control {right}")
    require(sum(comb(v, r) for v, r in zip(left, (2, 4, 6, 9))) == hole - 1, "left equality")
    require(sum(comb(v, r) for v, r in zip(right, (2, 4, 6, 9))) == hole + 1, "right equality")
    print(f"HOSTILE_2469_N={hole} RECIPROCAL_SUM={Fraction(37, 36)} TRIPLES={triples} REPS=0")
    print(f"HOSTILE_CONTROLS={hole-1}:{left};{hole+1}:{right}")
    print(f"SEMANTIC_SHA256={sha256(';'.join(semantic).encode('ascii')).hexdigest()}")
    print("ALL_AVERAGE_ORDER_GATES_PASSED")


if __name__ == "__main__":
    main()
