#!/usr/bin/env python3
"""Exact probe for p|d with a small quotient d/p."""

import importlib.util
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
SOURCE = ROOT / "04-computation/factorial_multiplace_newton_degree_barcode_thm3152.py"
spec = importlib.util.spec_from_file_location("thm3152", SOURCE)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def main():
    cases = ((22, 11), (214, 107), (999, 37), (1002, 167), (1384, 173))
    print("LARGE DIVISOR FLAG BLOCK EXACT PROBE")
    for d, prime in cases:
        p, q = c.pair(d)
        r = c.first_full_pseudoremainder(p, q, d)
        _, pair_blocks = c.degree_barcode_many((p, q), prime)
        degrees, flag_blocks = c.degree_barcode_many((p, q, r), prime)
        expected = ((c.Fraction(1, prime), d - prime, prime),)
        print(
            "d=%d p=%d m=%d pair=%s flag=%s expected_flag=%s match=%s degree1=%s"
            % (
                d,
                prime,
                d // prime,
                tuple((str(s), cap, den) for s, cap, den in pair_blocks),
                tuple((str(s), cap, den) for s, cap, den in flag_blocks),
                tuple((str(s), cap, den) for s, cap, den in expected),
                flag_blocks == expected,
                1 in degrees,
            )
        )
    hostile = []
    tested = 0
    cache = {}
    for prime in tuple(p for p in range(5, 44) if c.is_prime(p)):
        for multiplier in range(2, min(8, (prime - 1) // 2) + 1):
            d = multiplier * prime
            if d not in cache:
                ppoly, qpoly = c.pair(d)
                rem = c.first_full_pseudoremainder(ppoly, qpoly, d)
                cache[d] = (ppoly, qpoly, rem)
            rows = cache[d]
            _, actual = c.degree_barcode_many(rows, prime)
            expected = ((c.Fraction(1, prime), (multiplier - 1) * prime, prime),)
            tested += 1
            if actual != expected:
                hostile.append(
                    (
                        d,
                        prime,
                        multiplier,
                        tuple((str(s), cap, den) for s, cap, den in expected),
                        tuple((str(s), cap, den) for s, cap, den in actual),
                    )
                )
    print("symbolic-lane probe p>2m, p<=43, 2<=m<=8 tested=%d failures=%d" % (tested, len(hostile)))
    print("symbolic-lane failures=%s" % (hostile,))
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
