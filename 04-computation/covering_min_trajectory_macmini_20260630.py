#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is the convergent n/Phi_6 really the covering-min for n>=7? mac-mini-2026-06-30-S47.

REFUTED at n=7,8 (exact): covering-min = MEDIANT 2/(2n-1), achieved by sets living mod (2n-1):
  n=7: {1,2,5,6,7,8}      M=2/13   (t=4/13)
  n=8: {1,4,5,6,7,11,16}  M=2/15   (t=8/15)
both < the convergent n/Phi_6 (7/43, 8/57). This script pins the trajectory at n=9..16:
  (A) exhaustive-ish search (low speeds) for the covering-min;
  (B) targeted mediant construction mod (2n-1): a covering (n-1)-set whose speeds avoid {0,+-1}*k^{-1}.
If the mediant 2/(2n-1) is achievable at n=14, then opus's 14/183 covering-min is WRONG.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
from itertools import combinations
print = functools.partial(print, flush=True)


def Mfloat(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add(k/d)
    return max(min(abs(((v*t+0.5) % 1)-0.5) for v in S) for t in cand)


def Mexact(S):
    S = sorted(set(S)); cand = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    for k in range(1, d):
                        cand.add(F(k, d))
    best = F(0); at = None
    for t in cand:
        g = min(min((v*t) % 1, 1-((v*t) % 1)) for v in S)
        if g > best:
            best = g; at = t
    return best, at


def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))


def prim(S):
    g = 0
    for v in S:
        g = math.gcd(g, v)
    return g == 1


def targeted_mediant(n):
    """Find a covering (n-1)-set, all speeds avoiding {0,+-1}*k^{-1} mod (2n-1); return best M found."""
    m = 2*n-1
    bestM = None; bestS = None
    for k in range(1, m):
        if math.gcd(k, m) != 1:
            continue
        kinv = pow(k, -1, m)
        forbidden = {0, kinv, (m-kinv) % m}
        allowed = lambda v: (v % m) not in forbidden
        # greedy cover q=2..n with smallest allowed multiple of q
        S = set()
        for q in range(2, n+1):
            if any(v % q == 0 for v in S):
                continue
            j = 1
            while not allowed(q*j):
                j += 1
            S.add(q*j)
        # pad to exactly n-1 with smallest allowed speeds not yet in S
        v = 1
        while len(S) < n-1:
            if allowed(v) and v not in S:
                S.add(v)
            v += 1
        if len(S) != n-1 or not covers(S, n) or not prim(S):
            continue
        M = Mfloat(S)
        if bestM is None or M < bestM:
            bestM = M; bestS = sorted(S)
    return bestM, bestS


def main():
    print("COVERING-MIN TRAJECTORY: is the convergent really the covering-min for n>=7?")
    print("(REFUTED exact at n=7: {1,2,5,6,7,8} M=2/13<7/43; n=8: {1,4,5,6,7,11,16} M=2/15<8/57)\n")
    print(f"  {'n':>3} {'1/n':>8} {'mediant 2/(2n-1)':>17} {'conv n/Phi6':>12} | {'search(low) min':>16} | {'targeted-mediant':>16}")
    for n in range(7, 16):
        m = 2*n-1; P = n*n-n+1
        med, con = 2/m, n/P
        # (A) exhaustive-ish over low speeds
        B = 2*n + 4
        best = 1.0; arg = None
        if n <= 11:
            for S in combinations(range(1, B+1), n-1):
                if not covers(S, n):
                    continue
                if not prim(S):
                    continue
                M = Mfloat(S)
                if M < best - 1e-12:
                    best = M; arg = S
            sa = f"{best:.5f}"
        else:
            sa = "(skip)"
        # (B) targeted mediant
        tm, tms = targeted_mediant(n)
        flag = ""
        if tm is not None and tm < con - 1e-9:
            flag = " <-- BEATS convergent!"
        print(f"  {n:>3} {1/n:>8.5f} {med:>17.5f} {con:>12.5f} | {sa:>16} | {tm if tm else 0:>16.5f}{flag}")
        if n <= 11 and arg:
            print(f"        search argmin n={n}: {arg}")
        if tms:
            print(f"        targeted set  n={n}: {tms}")


if __name__ == "__main__":
    main()
