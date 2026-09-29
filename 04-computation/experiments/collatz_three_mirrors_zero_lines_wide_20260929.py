#!/usr/bin/env python3
"""Wide 3-adic zero-line census (numpy, residues mod 3^20), S23.

Same objects as collatz_three_mirrors_zero_lines_20260929.py: pairs (u, Q), u odd and prime to 3, 1 <= Q <= W, depth
d(u, Q) = v_3(u 2^Q -+ 1) for the sign that makes it >= 1, capped at 20 (residues mod 3^20 fit int64).  Reports the
depth histogram and the cumulative counts against the uniform model (P(depth >= d) = 3^-(d-1)), Poisson tail
probabilities for the observed deep counts, the deepest lines, and the depth-by-u structure (whether the deep lines
cluster on particular multipliers).
Run: python 04-computation/experiments/collatz_three_mirrors_zero_lines_wide_20260929.py [U] [W]
"""
from __future__ import annotations

import math
import sys

import numpy as np

CAP = 20
MOD = 3 ** CAP


def poisson_tail(lam: float, k: int) -> float:
    """P(Poisson(lam) >= k)."""
    s = 0.0
    term = math.exp(-lam)
    for i in range(k):
        s += term
        term *= lam / (i + 1)
    return max(0.0, 1.0 - s)


if __name__ == "__main__":
    U = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    us = np.array([u for u in range(1, U + 1, 2) if u % 3], dtype=np.int64)
    pw = np.empty(W + 1, dtype=np.int64)
    x = 1
    for Q in range(W + 1):
        pw[Q] = x
        x = (2 * x) % MOD
    N = len(us) * W
    hist = np.zeros(CAP + 2, dtype=np.int64)
    deep = []   # (depth, u, Q, sign)
    per_u_max = np.zeros(len(us), dtype=np.int64)
    for i, u in enumerate(us):
        r = (u * pw[1:]) % MOD                 # u 2^Q mod 3^20, Q = 1..W
        minus = (r - 1) % MOD
        plus = (r + 1) % MOD
        # depth of the sign that works: exactly one of minus, plus is divisible by 3
        sel_minus = (minus % 3 == 0)
        val = np.where(sel_minus, minus, plus)
        d = np.zeros(W, dtype=np.int64)
        v = val.copy()
        alive = np.ones(W, dtype=bool)
        for k in range(1, CAP + 1):
            alive &= (v % 3 == 0)
            d[alive] += 1
            v = np.where(alive, v // 3, v)
        # v == 0 means divisible by 3^20 (depth >= 20): cap
        d[val == 0] = CAP + 1
        hist += np.bincount(d, minlength=CAP + 2)
        per_u_max[i] = d.max()
        for Q in np.nonzero(d >= 11)[0]:
            deep.append((int(d[Q]), int(u), int(Q) + 1, "-" if sel_minus[Q] else "+"))
    print(f"census: {len(us)} multipliers u <= {U}, Q <= {W}: N = {N} pairs; depth capped at {CAP}")
    print("depth d: observed, expected N (2/3) 3^-(d-1), ratio | cumulative >= d: observed, expected, Poisson P(X >= observed)")
    for d in range(1, CAP + 2):
        obs = int(hist[d])
        exp = N * (2 / 3) * 3 ** (-(d - 1))
        cum = int(hist[d:].sum())
        cexp = N * 3 ** (-(d - 1))
        tail = poisson_tail(cexp, cum) if cum > 0 else 1.0
        print(f"   d={d:2d}: {obs:9d} {exp:13.3f} {obs / exp if exp > 0 else float('nan'):7.3f} | >= d: {cum:9d} {cexp:12.4f}  P = {tail:.3g}")
    print("deepest lines:")
    for d, u, Q, sgn in sorted(deep, reverse=True)[:30]:
        print(f"   depth {d:2d}: u = {u:5d}, Q = {Q:6d}, sign {sgn}")
    print("multipliers with a line of depth >= 12 (u : max depth):")
    idx = np.nonzero(per_u_max >= 12)[0]
    print("   " + ", ".join(f"{int(us[i])}:{int(per_u_max[i])}" for i in idx))
    print("DONE")
