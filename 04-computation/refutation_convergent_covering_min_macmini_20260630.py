#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DEFINITIVE refutation: the convergent n/Phi_6 is NOT the covering-min for n>=7.
mac-mini-2026-06-30-S47. Exact M (complete breakpoints {k/2s_i, k/(s_i+-s_j)}, MISTAKE-86 guard)
+ dense-grid cross-check (so the M's are not underestimates). HYP-3725 / CASE-convergent-not-covering-min.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
print = functools.partial(print, flush=True)


def Mexact(S):
    """max_t min_v ||v t||, exact. Complete candidate set t=k/d, d in {2s_i} U {s_i+-s_j}."""
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):   # i==j gives 2*S[i]; i!=j gives s_i+-s_j
                if d > 0:
                    for k in range(1, d):
                        cand.add(F(k, d))
    best = F(0); at = None
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in S)
        if g > best:
            best = g; at = t
    return best, at


def grid_max(S, N=2_000_000):
    best = 0.0; bt = 0.0
    for k in range(1, N):
        t = k/N
        g = min(abs(((v*t+0.5) % 1)-0.5) for v in S)
        if g > best:
            best = g; bt = t
    return best, bt


def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))


def prim(S):
    g = 0
    for v in S:
        g = math.gcd(g, v)
    return g == 1


def main():
    print("REFUTATION: the convergent n/Phi_6 is NOT the covering-min for n>=7.\n")
    cases = [
        (7, [1, 2, 5, 6, 7, 8]),
        (8, [1, 4, 5, 6, 7, 11, 16]),
        (9, [1, 3, 4, 5, 7, 11, 18, 32]),
    ]
    for n, S in cases:
        m, t = Mexact(S)
        gm, gt = grid_max(S)
        con = F(n, n*n-n+1)
        ok = covers(S, n) and prim(S) and len(set(S)) == n-1
        guard = abs(gm-float(m)) < 2e-4
        print(f"n={n}: {S}")
        print(f"   M = {m} = {float(m):.6f}  at t={t}   (grid-max {gm:.6f} -> no underestimate: {guard})")
        print(f"   construction n/Phi_6 = {con} = {float(con):.6f}   ==> beats construction: {m < con}")
        print(f"   valid covering (n-1)-set (covers 2..n, primitive, size n-1): {ok};  M>1/n: {m > F(1, n)}")
    print("\nCONCLUSION: there is no mediant->convergent transition at n=7; the convergent is NOT the")
    print("covering-min. LRC floor untouched (all M>1/n). See HYP-3725, MISTAKE-087, CASE-convergent-not-covering-min.")


if __name__ == "__main__":
    main()
