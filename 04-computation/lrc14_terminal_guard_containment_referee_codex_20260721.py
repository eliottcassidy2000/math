#!/usr/bin/env python3
"""Bounded exact referee for the THM-2073/2076 terminal residual.

The sampled bitset is only a lossless prefilter: every true continuous
containment passes every sampled phase.  Every prefilter survivor is then
checked on the complete rational boundary arrangement, including endpoints.

Tournament Analysis is inapplicable: the exact predicate is containment of a
closed multi-frequency safe set in one open guard set, not pairwise dominance.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd


MAX_CORE_SPEED = 24
MIN_CORE_SIZE = 6
MAX_CORE_SIZE = 10
MAX_GUARD = 175
GRID = 8192
FULL_DIVISOR_MASK = (1 << 13) - 1
DIVISOR_MASK = {
    q: sum(1 << (d - 2) for d in range(2, 15) if q % d == 0)
    for q in range(1, MAX_CORE_SPEED + 1)
}


def abs_residue(k: int, modulus: int) -> int:
    r = k % modulus
    return min(r, modulus - r)


def divisor_complete(Q: tuple[int, ...]) -> bool:
    mask = 0
    for q in Q:
        mask |= DIVISOR_MASK[q]
    return mask == FULL_DIVISOR_MASK


def hereditarily_primitive(Q: tuple[int, ...]) -> bool:
    return all(gcd(*(Q[:i] + Q[i + 1 :])) == 1 for i in range(len(Q)))


def safe_at(Q: tuple[int, ...], t: Fraction) -> bool:
    den = t.denominator
    return all(14 * abs_residue(q * t.numerator, den) >= den for q in Q)


def eligible_at(h: int, t: Fraction) -> bool:
    den = t.denominator
    return 7 * abs_residue(h * t.numerator, den) < den


def exact_containment(Q: tuple[int, ...], h: int) -> tuple[bool, Fraction | None]:
    boundaries = {Fraction(0), Fraction(1)}
    for q in Q:
        for k in range(q + 1):
            for sign in (-1, 1):
                t = Fraction(14 * k + sign, 14 * q)
                if 0 <= t <= 1:
                    boundaries.add(t)
    for k in range(h + 1):
        for sign in (-1, 1):
            t = Fraction(7 * k + sign, 7 * h)
            if 0 <= t <= 1:
                boundaries.add(t)

    points = sorted(boundaries)
    tests = set(points)
    tests.update((a + b) / 2 for a, b in zip(points, points[1:]))
    for t in sorted(tests):
        if safe_at(Q, t) and not eligible_at(h, t):
            return False, t
    return True, None


def sampled_masks() -> tuple[list[int], dict[int, int]]:
    safe = [0] * (MAX_CORE_SPEED + 1)
    for q in range(1, MAX_CORE_SPEED + 1):
        mask = 0
        for j in range(GRID):
            if 14 * abs_residue(q * j, GRID) >= GRID:
                mask |= 1 << j
        safe[q] = mask

    eligible: dict[int, int] = {}
    for h in range(1, MAX_GUARD + 1, 2):
        mask = 0
        for j in range(GRID):
            if 7 * abs_residue(h * j, GRID) < GRID:
                mask |= 1 << j
        eligible[h] = mask
    return safe, eligible


def main() -> None:
    print("LRC14 TERMINAL GUARD-CONTAINMENT SEARCH -- sampled prefilter + exact atoms")
    safe_masks, eligible_masks = sampled_masks()
    full = (1 << GRID) - 1
    totals: dict[int, int] = {}
    arithmetic_counts: dict[int, int] = {}
    prefilter_counts: dict[int, int] = {}
    exact_counts: dict[int, int] = {}
    prefilter_pairs = 0
    exact_pairs = 0
    first_witness = None

    for size in range(MIN_CORE_SIZE, MAX_CORE_SIZE + 1):
        total = 0
        arithmetic = 0
        size_prefilter = 0
        size_exact = 0
        for Q in combinations(range(1, MAX_CORE_SPEED + 1), size):
            total += 1
            if not divisor_complete(Q) or not hereditarily_primitive(Q):
                continue
            arithmetic += 1
            safe = full
            for q in Q:
                safe &= safe_masks[q]

            # THM-2073's maximizer interval gives the strict bound
            # (13-s)h < 2(s+1)max(Q) for a terminal guard.
            guard_numerator = 2 * (size + 1) * max(Q)
            guard_denominator = 13 - size
            for h, eligible in eligible_masks.items():
                if guard_denominator * h >= guard_numerator:
                    continue
                if safe & ~eligible:
                    continue
                size_prefilter += 1
                prefilter_pairs += 1
                contains, failure = exact_containment(Q, h)
                if contains:
                    size_exact += 1
                    exact_pairs += 1
                    if first_witness is None:
                        first_witness = (Q, h)
                elif failure is None:
                    raise AssertionError("failed exact pair lacks a witness")
        totals[size] = total
        arithmetic_counts[size] = arithmetic
        prefilter_counts[size] = size_prefilter
        exact_counts[size] = size_exact

    for size in range(MIN_CORE_SIZE, MAX_CORE_SIZE + 1):
        print(
            f"size {size} through {MAX_CORE_SPEED}: {totals[size]} total, "
            f"{arithmetic_counts[size]} hereditary divisor-complete, "
            f"{prefilter_counts[size]} prefilter, {exact_counts[size]} exact"
        )
    print(f"odd guards through {MAX_GUARD}: {(MAX_GUARD + 1) // 2}")
    print(f"sampled prefilter survivors: {prefilter_pairs}")
    print(f"exact containment pairs: {exact_pairs}")
    print(f"first exact pair: {first_witness}")
    print("PASS")


if __name__ == "__main__":
    main()
