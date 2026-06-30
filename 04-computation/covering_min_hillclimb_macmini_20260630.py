#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hill-climbing search to beat the convergent n/Phi_6 at n=10..14 (esp. n=14 vs opus's 14/183).
mac-mini-2026-06-30-S47. Local search: swap speeds to reduce M while staying a primitive covering (n-1)-set.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
import numpy as np
print = functools.partial(print, flush=True)


class LCG:
    def __init__(self, s): self.s = s & 0x7fffffff
    def r(self): self.s = (1103515245*self.s + 12345) & 0x7fffffff; return self.s


def Mfloat(S):
    S = sorted(set(S))
    diffs = set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]-S[j], S[i]+S[j]):
                if d > 0:
                    diffs.add(d)
    ts = np.array(sorted({k/d for d in diffs for k in range(1, d)}))
    sv = np.array(S)
    A = np.abs(((np.outer(ts, sv)+0.5) % 1.0)-0.5)
    return float(A.min(axis=1).max())


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


def climb(n, restarts, steps, seed, vmax):
    rng = LCG(seed)
    bestM, bestS = 2.0, None
    for _ in range(restarts):
        # random covering start
        S = set(range(2, n+1))  # consecutive covers 2..n, size n-1
        S = set(sorted(S))
        cur = list(S)
        if len(cur) != n-1:
            continue
        curM = Mfloat(cur)
        for _ in range(steps):
            i = rng.r() % len(cur)
            old = cur[i]
            new = rng.r() % vmax + 1
            if new in cur:
                continue
            cand = cur[:i] + [new] + cur[i+1:]
            if not covers(cand, n) or not prim(cand) or len(set(cand)) != n-1:
                continue
            m = Mfloat(cand)
            if m < curM - 1e-12 or (m < curM + 0.003 and rng.r() % 5 == 0):  # accept downhill + occasional sideways
                cur, curM = cand, m
                if curM < bestM:
                    bestM, bestS = curM, sorted(cur)
    return bestM, bestS


def main():
    print("HILL-CLIMB: beat the convergent n/Phi_6? (esp n=14 vs opus 14/183)\n")
    print(f"  {'n':>3} {'1/n':>8} {'med 2/(2n-1)':>13} {'conv n/Phi6':>12} | {'best M':>9}  verdict")
    for n in range(9, 15):
        med, con = 2/(2*n-1), n/(n*n-n+1)
        bM, bS = climb(n, restarts=60, steps=400, seed=999+n, vmax=4*n)
        beats = bM < con - 1e-9
        ex = ""
        if bS is not None and beats:
            em, et = Mexact(bS)
            ex = f"  EXACT M={em}={float(em):.5f} at t={et}"
        v = ("BEATS conv!" + ex) if beats else "conv not beaten"
        print(f"  {n:>3} {1/n:>8.5f} {med:>13.5f} {con:>12.5f} | {bM:>9.5f}  {v}")
        if bS is not None:
            print(f"        set: {bS}")


if __name__ == "__main__":
    main()
