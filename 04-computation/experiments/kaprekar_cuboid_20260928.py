#!/usr/bin/env python3
"""Kaprekar's routine on d-digit numbers (descending digits minus ascending digits, leading zeros kept):
cycles for d = 2..6, the divisibility-by-9 invariant, and the digit-multiset factoring (an exact invariant
in the sense of the nine-iteration typology).  Also the smallest Euler bricks (perfect cuboid context).
Session: opus, collatz-poset-dag-20260927 (S19), 2026-09-28.
Run: python 04-computation/experiments/kaprekar_cuboid_20260928.py
"""
from __future__ import annotations

import math
from collections import Counter


def kaprekar(n: int, d: int) -> int:
    s = f"{n:0{d}d}"
    return int("".join(sorted(s, reverse=True))) - int("".join(sorted(s)))


def cycles(d: int):
    seen_cycle = {}
    counts = Counter()
    for n in range(10 ** d):
        if len(set(f"{n:0{d}d}")) == 1:
            continue  # repdigits go to 0
        path = []
        x = n
        visited = {}
        while x not in visited:
            visited[x] = len(path)
            path.append(x)
            x = kaprekar(x, d)
        cyc = path[visited[x]:]
        key = tuple(sorted(cyc))
        seen_cycle[key] = cyc
        counts[key] += 1
    return seen_cycle, counts


if __name__ == "__main__":
    print("== Kaprekar's routine: cycles by digit count ==")
    for d in range(2, 7):
        cyc, counts = cycles(d)
        desc = []
        for key, c in sorted(cyc.items(), key=lambda kv: -counts[kv[0]]):
            desc.append(f"{c} (len {len(c)}, {counts[key]} starts)")
        print(f"   d={d}: " + "; ".join(desc))
    # invariant: every value after one step is divisible by 9 (digit sum preserved by sorting)
    for d in range(2, 7):
        assert all(kaprekar(n, d) % 9 == 0 for n in range(10 ** d))
    print("   every image is divisible by 9 (n - reverse-sorted differs from n by a multiple of 9); the map factors through the digit multiset: a finite state space per d")
    print(f"   6174 = {6174} = 2*3^2*7^3, 495 = 3^2*5*11; the 2-digit cycle 9 -> 81 -> 63 -> 27 -> 45 -> 9 is 9 * (1, 9, 7, 3, 5)")
    # Euler bricks: integer edges with integer face diagonals; perfect cuboid also needs an integer space diagonal
    bricks = []
    for a in range(1, 300):
        for b in range(a + 1, 300):
            ab = a * a + b * b
            r = math.isqrt(ab)
            if r * r != ab:
                continue
            for c in range(b + 1, 300):
                ac, bc = a * a + c * c, b * b + c * c
                if math.isqrt(ac) ** 2 == ac and math.isqrt(bc) ** 2 == bc:
                    sp = a * a + b * b + c * c
                    bricks.append((a, b, c, math.isqrt(sp) ** 2 == sp))
    print(f"   Euler bricks with edges < 300: {bricks[:6]} (last flag: integer space diagonal, i.e. a perfect cuboid; none)")
    print("DONE")
