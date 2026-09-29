#!/usr/bin/env python3
"""3-adic zero lines of the multiplier families: the census behind the ridges of the negative family (S23, 2026-09-29).

A ridge of the negative-power family (section 2d of collatz_five_mirrors_20260929.md) is born at a coincidence
    u * 2^Q = -+1 (mod 3^d)   with u small,
which puts the multiplier family (-+u) 2^j at the negative exponents -Q + j for every level n <= d ("depth" d = the
3-adic valuation of u 2^Q -+ 1).  This is the shape of Chocian's twisted-Bernoulli zero lines (arXiv:2608.08724):
a p-adic coincidence per (character, prime) pair, counted against a Poisson null model, with a depth.  Here the pairs
are (u, Q), u odd, 3 does not divide u, 1 <= Q <= W; exactly one sign has depth >= 1 (u 2^Q is +-1 mod 3), and under
the uniform model P(depth >= d) = 3^-(d-1), P(depth = d) = (2/3) 3^-(d-1).
Reports: the depth histogram against the geometric law, the number of lines of depth >= d against the Poisson
expectation, the deepest lines (with the seed (55, 423, -) of depth 15 of the S22 note), the depth-by-Q profile for
the window Q <= 60 that enters Ntilde_n, and the predicted ridge heights 0.2 M(d) 3^(d/2) of the deepest lines.
Run: python 04-computation/experiments/collatz_three_mirrors_zero_lines_20260929.py [U] [W]
"""
from __future__ import annotations

import math
import sys
from collections import Counter


def v3(x: int) -> int:
    if x == 0:
        return 10 ** 9
    d = 0
    while x % 3 == 0:
        x //= 3
        d += 1
    return d


if __name__ == "__main__":
    U = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    us = [u for u in range(1, U + 1, 2) if u % 3]
    pow2 = [pow(2, Q) for Q in range(W + 1)]
    hist = Counter()
    lines = []          # (depth, u, Q, sign)
    depth_by_Q = {}
    for u in us:
        for Q in range(1, W + 1):
            x = u * pow2[Q]
            if (x - 1) % 3 == 0:
                d, sgn = v3(x - 1), "-"
            else:
                d, sgn = v3(x + 1), "+"
            hist[d] += 1
            if d >= 6:
                lines.append((d, u, Q, sgn))
            if Q <= 60 and d >= 5:
                depth_by_Q.setdefault(Q, []).append((d, u, sgn))
    N = len(us) * W
    print(f"census: {len(us)} multipliers u <= {U} (odd, 3 does not divide u), Q <= {W}: N = {N} pairs")
    print("depth d: observed count, expected N (2/3) 3^-(d-1), ratio; cumulative >= d: observed, expected N 3^-(d-1)")
    dmax = max(hist)
    for d in range(1, dmax + 1):
        obs = hist[d]
        exp = N * (2 / 3) * 3 ** (-(d - 1))
        cum = sum(hist[e] for e in range(d, dmax + 1))
        cexp = N * 3 ** (-(d - 1))
        print(f"   d={d:2d}: {obs:8d}  {exp:12.2f}  {obs / exp if exp > 0 else float('nan'):6.3f}   | >= d: {cum:8d}  {cexp:12.3f}  {cum / cexp if cexp > 0 else float('nan'):6.3f}")
    print("deepest lines (depth, u, Q, sign of u 2^Q -+ 1 = 0 mod 3^depth):")
    for d, u, Q, sgn in sorted(lines, reverse=True)[:25]:
        print(f"   depth {d:2d}: u = {u:5d}, Q = {Q:5d}, sign {sgn}   (2^-{Q} = {'-' if sgn == '+' else ''}{u} mod 3^{d})")
    # the S22 seed and its neighbourhood
    for (u, Q) in ((55, 423), (13, 154), (1, 154), (1, 480), (1, 486)):
        x = u * pow2[Q]
        print(f"   check (u, Q) = ({u}, {Q}): v3(u 2^Q - 1) = {v3(x - 1)}, v3(u 2^Q + 1) = {v3(x + 1)}")
    print("lines of depth >= 5 in the window Q <= 60 (these place a multiplier family at the negative exponents -Q + j entering Ntilde_n for n <= depth):")
    for Q in sorted(depth_by_Q):
        items = sorted(depth_by_Q[Q], reverse=True)[:4]
        print(f"   Q = {Q:2d}: " + ", ".join(f"depth {d} (u = {u}{s})" for d, u, s in items))
    # Poisson-model expectation of the deepest line in a box, and the ridge-height prediction
    print("expected number of lines of depth >= d in this box under the uniform model, and the depth at which it drops below 1:")
    for d in range(8, 20):
        e = N * 3 ** (-(d - 1))
        print(f"   d >= {d:2d}: {e:10.4f}")
    print("DONE")
