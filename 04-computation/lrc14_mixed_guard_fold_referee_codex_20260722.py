#!/usr/bin/env python3
"""Exact referee for THM-2080's mixed 1/14--1/7 fold law.

The direct calculation samples one midpoint in every cell of a common exact
boundary grid; it is not a floating mesh.  The terminal replay is the rank-six
arithmetic slice of THM-2078 and tests the necessary first-order union bound.

Tournament Analysis is inapplicable: simultaneous set cover is not pairwise
dominance.  The exact carrier is the signed fold invoice plus endpoint phase.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd


MAX_FORMULA_Q = 40
MAX_FORMULA_H = 39
MAX_CORE_SPEED = 24
FULL_DIVISOR_MASK = (1 << 13) - 1
DIVISOR_MASK = {
    q: sum(1 << (d - 2) for d in range(2, 15) if q % d == 0)
    for q in range(1, MAX_CORE_SPEED + 1)
}


def fold14(r: int) -> int:
    r %= 14
    return r * (14 - r)


def folded_overlap(q: int, h: int) -> Fraction:
    d = gcd(q, h)
    a, b = q // d, h // d
    numerator = 8 * a * b + fold14(b + 2 * a) - fold14(b - 2 * a)
    return Fraction(numerator, 196 * a * b)


def atom_overlap(q: int, h: int) -> Fraction:
    # Every endpoint of either comb lies on the grid of denominator 14*q*h.
    # The predicate is constant on the open cells, so midpoint counting is exact.
    cells = 14 * q * h
    denominator = 2 * cells
    count = 0
    for j in range(cells):
        numerator = 2 * j + 1
        rq = (q * numerator) % denominator
        rh = (h * numerator) % denominator
        rq = min(rq, denominator - rq)
        rh = min(rh, denominator - rh)
        if 14 * rq < denominator and 7 * rh < denominator:
            count += 1
    return Fraction(count, cells)


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def main() -> None:
    formula_checks = 0
    maximum_correction = Fraction(-1)
    maximum_case = None
    for q in range(1, MAX_FORMULA_Q + 1):
        for h in range(1, MAX_FORMULA_H + 1, 2):
            direct = atom_overlap(q, h)
            folded = folded_overlap(q, h)
            assert direct == folded, (q, h, direct, folded)
            d = gcd(q, h)
            correction = folded - Fraction(2, 49)
            assert correction <= Fraction(d * d, 4 * q * h)
            if correction > maximum_correction:
                maximum_correction = correction
                maximum_case = (q, h, correction)
            formula_checks += 1

    core_count = 0
    guard_pairs = 0
    union_bound_survivors = 0
    resonance_failures = 0
    first_survivor = None
    for Q in combinations(range(1, MAX_CORE_SPEED + 1), 6):
        if not divisor_complete(Q) or not hereditarily_primitive(Q):
            continue
        core_count += 1
        # For s=6, THM-2077 gives 7h<14B, equivalently h<2B.
        for h in range(1, 2 * max(Q), 2):
            guard_pairs += 1
            overlap_sum = sum((folded_overlap(q, h) for q in Q), Fraction())
            if overlap_sum < Fraction(2, 7):
                continue
            union_bound_survivors += 1
            if first_survivor is None:
                first_survivor = (Q, h, overlap_sum)
            products = []
            for q in Q:
                d = gcd(q, h)
                products.append((q // d) * (h // d))
            if min(products) > 36:
                resonance_failures += 1

    assert core_count == 6
    assert guard_pairs == 144
    assert union_bound_survivors == 28
    assert resonance_failures == 0

    print("LRC14 MIXED GUARD FOLD REFEREE -- exact atoms and rank-six gate")
    print(f"formula cases: {formula_checks}")
    print(f"maximum positive correction: {maximum_case}")
    print(f"rank-six hereditary divisor-complete cores through 24: {core_count}")
    print(f"allowed odd guard pairs: {guard_pairs}")
    print(f"first-order exact-overlap survivors: {union_bound_survivors}")
    print(f"first survivor: {first_survivor}")
    print(f"survivors violating reduced-product <=36: {resonance_failures}")
    print("PASS")


if __name__ == "__main__":
    main()
