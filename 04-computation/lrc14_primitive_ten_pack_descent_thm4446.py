#!/usr/bin/env python3
"""Exact audit for THM-4446's primitive ten-pack descent.

If S=3C union T is primitive, |C|=10, every t in T is a ternary unit,
and gcd(C)>1, choose a prime p|gcd(C), a C/p-safe phase y, and the common
3p lifts x_j=(y+j)/(3p).  This audit checks the open-arc branch caps, all
small-prime margins, and 6,292 literal controls.  Arithmetic is exact.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
from lrc14_bounded_ten_body_geometry_thm4442 import (  # noqa: E402
    safe_components_by_union,
)


DELTA = Q(1, 14)
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(x: Q) -> Q:
    return x - (x.numerator // x.denominator)


def norm(x: Q) -> Q:
    z = frac(x)
    return min(z, 1 - z)


def danger(x: Q) -> bool:
    # Full-row safety is weak, so the complementary danger teeth are open.
    return norm(x) < DELTA


def critical_phase_max(order: int) -> tuple[int, Q]:
    """Maximum points of a translated order-grid in the open 1/7 arc."""
    endpoints = sorted(
        {
            frac(sign * DELTA - Q(j, order))
            for sign in (-1, 1)
            for j in range(order)
        }
    )
    probes = set(endpoints)
    for i, left in enumerate(endpoints):
        right = endpoints[(i + 1) % len(endpoints)]
        if i + 1 == len(endpoints):
            right += 1
        probes.add(frac((left + right) / 2))
    scored = [
        (sum(danger(alpha + Q(j, order)) for j in range(order)), alpha)
        for alpha in probes
    ]
    return max(scored)


def expected_grid_cap(order: int) -> int:
    # This remains exact when 7|order because the danger arc is open.
    return (order + 6) // 7


def tail_catalogue(p: int, divisible_count: int) -> tuple[int, int, int]:
    require(0 <= divisible_count <= 2, "primitivity leaves at most two")
    if p == 3:
        require(divisible_count == 0, "ternary units cannot be 3-divisible")
    multiples = [n for n in range(1, 20 * p + 20) if n % 3 and n % p == 0]
    nonmultiples = [n for n in range(1, 20 * p + 20) if n % 3 and n % p]
    chosen = multiples[:divisible_count] + nonmultiples[: 3 - divisible_count]
    require(len(chosen) == 3 and len(set(chosen)) == 3, (p, divisible_count, chosen))
    return tuple(sorted(chosen))


def literal_row_audit(
    base: tuple[int, ...], p: int, tails: tuple[int, int, int]
) -> tuple[int, int]:
    """Return (number of safe lifts, certified union-bound margin)."""
    body = tuple(p * b for b in base)
    row = tuple(3 * c for c in body) + tails
    require(len(set(row)) == 13, (base, p, tails, "collision"))
    require(gcd(*row) == 1, (base, p, tails, "nonprimitive control"))
    require(all(t % 3 for t in tails), tails)
    require(gcd(*body) > 1, body)

    components = safe_components_by_union(base)
    require(bool(components), (base, "empty good set"))
    left, right = max(components, key=lambda interval: interval[1] - interval[0])
    y = (left + right) / 2
    require(all(norm(b * y) >= DELTA for b in base), (base, y))

    lifts = [frac((y + j) / (3 * p)) for j in range(3 * p)]
    winners = [
        x for x in lifts if all(norm(speed * x) >= DELTA for speed in row)
    ]

    predicted_bad = 0
    for t in tails:
        divisor = gcd(t, 3 * p)
        order = 3 * p // divisor
        predicted_bad += divisor * expected_grid_cap(order)
        actual_bad = sum(danger(t * x) for x in lifts)
        require(
            actual_bad <= divisor * expected_grid_cap(order),
            (p, tails, t, actual_bad),
        )
    margin = 3 * p - predicted_bad
    require(margin > 0, (p, tails, predicted_bad, "nonpositive union margin"))
    require(len(winners) >= margin, (base, p, tails, len(winners), margin))
    return len(winners), margin


def main() -> None:
    print("PRIMITIVE TEN-PACK REDUCTION -- EXACT HOSTILE AUDIT")
    print("threshold=1/14; danger teeth open; safe endpoints retained")

    grid_records = []
    for order in range(1, 92):
        observed, phase = critical_phase_max(order)
        expected = expected_grid_cap(order)
        require(observed == expected, (order, observed, expected, phase))
        grid_records.append((order, observed, phase))
    print("translated order-grid caps: exact ceil(m/7) for m=1..91")
    print(
        "boundary controls: m=7 ->",
        critical_phase_max(7)[0],
        "; m=21 ->",
        critical_phase_max(21)[0],
    )

    print("prime | full-order cap | p-divisible cap | worst bad | margin")
    analytic_margins = []
    for p in PRIMES:
        full = expected_grid_cap(3 * p)
        repeated = p * expected_grid_cap(3)
        worst = 3 * full if p == 3 else 2 * repeated + full
        margin = 3 * p - worst
        require(full < p, (p, full))
        require(margin > 0, (p, worst, margin))
        analytic_margins.append(margin)
        print(f"{p:5d} | {full:14d} | {repeated:15d} | {worst:9d} | {margin:6d}")

    rows = 0
    min_winners = None
    min_margin = None
    digest = hashlib.sha256()
    for base in combinations(range(1, 14), 10):
        for p in PRIMES:
            possible_counts = (0,) if p == 3 else (0, 1, 2)
            for divisible_count in possible_counts:
                tails = tail_catalogue(p, divisible_count)
                winners, margin = literal_row_audit(base, p, tails)
                rows += 1
                min_winners = winners if min_winners is None else min(min_winners, winners)
                min_margin = margin if min_margin is None else min(min_margin, margin)
                digest.update(f"{base}|{p}|{tails}|{winners}|{margin}\n".encode())

    expected_rows = 286 * (1 + 7 * 3)
    require(rows == expected_rows, (rows, expected_rows))
    require(min_winners is not None and min_winners >= 1, min_winners)
    require(min_margin == min(analytic_margins), (min_margin, analytic_margins))
    print(f"literal rows audited: {rows}")
    print(f"minimum actual safe lifts: {min_winners}")
    print(f"minimum certified union-bound margin: {min_margin}")
    print("literal-record sha256:", digest.hexdigest())
    print("PASS")


if __name__ == "__main__":
    main()
