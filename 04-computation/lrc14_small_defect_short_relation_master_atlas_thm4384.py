#!/usr/bin/env python3
"""Explore every parity-compatible one-determinant short relation.

The analytic range is m >= h+1, m+h <= 13, gcd(3,mh)=1.  For odd speeds,
parity further requires m == h+1 (mod 2).
"""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "comb141",
    ROOT / "04-computation/lrc14_signed_one_four_one_comb_exact_measure_thm4382.py",
)
COMB = module_from_spec(SPEC)
SPEC.loader.exec_module(COMB)

R = Fraction(3, 14)
PATTERNS = [
    (h, m)
    for h in range(1, 14)
    for m in range(h + 1, 14 - h)
    if h % 3 and m % 3 and (m - h - 1) % 2 == 0
]


def error(t):
    return sum(
        (t - k for k in range(1, int(t) + 1) if k % 3),
        Fraction(0),
    ) - t * t / 3


def measure(m, p, q):
    a = R * abs(q - p) / m
    b = R * (p + q) / m
    return Fraction(6, 49 * m) + Fraction(2 * m, p * q) * (error(b) - error(a))


def representations(h, m, product_bound):
    for p in range(1, product_bound, 2):
        if p % 3 == 0:
            continue
        max_q = (product_bound - 1) // p
        for q in range(1, max_q + 1, 2):
            if q % 3 == 0:
                continue
            for s, t in product((-1, 1), repeat=2):
                numerator = s * h * p + t * q
                if numerator <= 0 or numerator % m:
                    continue
                b = numerator // m
                if b % 2 == 0 or b % 3 == 0 or len({p, b, q}) != 3:
                    continue
                if gcd(gcd(p, b), q) != 1:
                    continue
                yield p, b, q, s, t


def ceiling(m, product_value):
    return Fraction(6, 49 * m) + Fraction(2 * m, 3 * product_value)


def cutoff_for(m, target):
    gap = target - Fraction(6, 49 * m)
    assert gap > 0
    ratio = Fraction(2 * m, 3) / gap
    return ratio.numerator // ratio.denominator + 1


def main():
    discovery_bound = 20000
    print(f"patterns={PATTERNS}")
    sectors = {}
    for h, m in PATTERNS:
        discovered = {}
        for p, b, q, s, t in representations(h, m, discovery_bound):
            triple = tuple(sorted((p, b, q)))
            value = measure(m, p, q)
            discovered.setdefault(triple, []).append((value, p, b, q, s, t))
        ranked = sorted(
            ((rows[0][0], triple, rows) for triple, rows in discovered.items()),
            reverse=True,
        )
        target = ranked[0][0]
        cutoff = cutoff_for(m, target)
        assert cutoff < discovery_bound
        exact = {}
        for p, b, q, s, t in representations(h, m, cutoff):
            triple = tuple(sorted((p, b, q)))
            value = measure(m, p, q)
            direct = COMB.failure_measure_direct(triple)
            assert value == direct, ((h, m), p, b, q, s, t, value, direct)
            exact.setdefault(triple, []).append((value, p, b, q, s, t))
        exact_ranked = sorted(
            ((rows[0][0], triple, rows) for triple, rows in exact.items()),
            reverse=True,
        )
        assert exact_ranked[0][0] == target
        assert ceiling(m, cutoff) < target
        maximizers = [row for row in exact_ranked if row[0] == target]
        sectors[(h, m)] = set(exact)
        print(
            f"pattern=({h},{m},1) max={target} maximizers={maximizers} "
            f"cutoff={cutoff} ceiling={ceiling(m,cutoff)} "
            f"presentations={sum(len(x[2]) for x in exact_ranked)} "
            f"triples={len(exact_ranked)}"
        )
        print(f"top5={exact_ranked[:5]}")

    for left_index, left in enumerate(PATTERNS):
        for right in PATTERNS[left_index + 1 :]:
            overlap = sectors[left] & sectors[right]
            if overlap:
                print(f"small-product overlap {left}<->{right}: {sorted(overlap)[:20]}")
    print("PASS")


if __name__ == "__main__":
    main()
