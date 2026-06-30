#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Direct attack: can ANY covering (n-1)-set beat the convergent n/Phi_6 at n=9..14?
mac-mini-2026-06-30-S47. If yes at n=14, opus's 14/183 covering-min is WRONG.

Two engines: (A) targeted mediant mod (2n-1) over all k and small allowed-speed exhaustion;
(B) randomized greedy covering sets. Tracks the global min M per n. numpy-vectorized M.
"""
from __future__ import annotations
import functools, math
from fractions import Fraction as F
from itertools import combinations, product
import numpy as np
print = functools.partial(print, flush=True)

# deterministic PRNG (no Math.random / Date in this env's spirit; use a fixed LCG)
class LCG:
    def __init__(self, s): self.s = s
    def rnd(self): self.s = (1103515245*self.s + 12345) & 0x7fffffff; return self.s
    def choice(self, lst): return lst[self.rnd() % len(lst)]


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
    # min over speeds of ||v t|| for each t
    A = np.outer(ts, sv)
    A = np.abs(((A+0.5) % 1.0)-0.5)
    return float(A.min(axis=1).max())


def covers(S, n):
    return all(any(v % q == 0 for v in S) for q in range(2, n+1))


def prim(S):
    g = 0
    for v in S:
        g = math.gcd(g, v)
    return g == 1


def targeted(n, maxmult=12):
    """For each k mod (2n-1), exhaust small covering sets among allowed speeds; return best (M,S)."""
    m = 2*n-1
    bestM, bestS = 2.0, None
    for k in range(1, m):
        if math.gcd(k, m) != 1:
            continue
        kinv = pow(k, -1, m)
        forb = {0, kinv, (m-kinv) % m}
        # allowed speeds up to maxmult*n
        allowed = [v for v in range(1, maxmult*n+1) if (v % m) not in forb]
        # for each q, the allowed multiples (a few smallest)
        permult = {q: [v for v in allowed if v % q == 0][:6] for q in range(2, n+1)}
        if any(not permult[q] for q in range(2, n+1)):
            continue
        # greedy-ish: cover hardest q first (fewest options), build set, then verify M
        order = sorted(range(2, n+1), key=lambda q: len(permult[q]))
        S = set()
        ok = True
        for q in order:
            if any(v % q == 0 for v in S):
                continue
            added = False
            for v in permult[q]:
                if len(S) < n-1:
                    S.add(v); added = True; break
            if not added:
                ok = False; break
        if not ok:
            continue
        # pad with smallest allowed
        i = 0
        while len(S) < n-1 and i < len(allowed):
            S.add(allowed[i]); i += 1
        if len(S) != n-1 or not covers(S, n) or not prim(S):
            continue
        M = Mfloat(S)
        if M < bestM:
            bestM, bestS = M, sorted(S)
    return bestM, bestS


def randomized(n, iters, seed):
    rng = LCG(seed)
    m = 2*n-1
    bestM, bestS = 2.0, None
    mults = {q: [q*j for j in range(1, 8)] for q in range(2, n+1)}
    for _ in range(iters):
        S = set()
        for q in range(2, n+1):
            if any(v % q == 0 for v in S):
                continue
            S.add(rng.choice(mults[q]))
        # trim/pad to n-1
        Sl = sorted(S)
        while len(Sl) < n-1:
            v = rng.rnd() % (3*n) + 1
            if v not in Sl:
                Sl.append(v)
        Sl = sorted(set(Sl))
        if len(Sl) != n-1 or not covers(Sl, n) or not prim(Sl):
            continue
        M = Mfloat(Sl)
        if M < bestM:
            bestM, bestS = M, Sl
    return bestM, bestS


def main():
    print("DIRECT ATTACK: can any covering (n-1)-set beat the convergent n/Phi_6?\n")
    print(f"  {'n':>3} {'1/n':>8} {'med 2/(2n-1)':>13} {'conv n/Phi6':>12} | {'targeted min':>13} {'rand min':>10} {'GLOBAL min':>11}  verdict")
    for n in range(7, 15):
        med, con = 2/(2*n-1), n/(n*n-n+1)
        tM, tS = targeted(n)
        rM, rS = randomized(n, 40000, 12345 + n)
        gM = min(tM, rM)
        gS = tS if tM <= rM else rS
        verdict = "MEDIANT beats conv" if gM < con-1e-9 and abs(gM-med) < 1e-6 else \
                  ("BEATS conv" if gM < con-1e-9 else "conv not beaten (low search)")
        print(f"  {n:>3} {1/n:>8.5f} {med:>13.5f} {con:>12.5f} | {tM:>13.5f} {rM:>10.5f} {gM:>11.5f}  {verdict}")
        print(f"        best set n={n}: {gS}")


if __name__ == "__main__":
    main()
