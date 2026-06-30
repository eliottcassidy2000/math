#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Smart covering-min search: map the exact margin deviation + hunt odd-n winners.
mac-mini-2026-06-30-S50. Structured seeds (QR/Paley mod m, spread APs, scaled blocks) + local search.
Tracks the best primitive covering (n-1)-set per n; reports M, margin = M - 1/n, and the witness modulus.
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


def fix_to_covering(S, n, rng, vmax):
    """greedily repair S to a primitive covering (n-1)-set."""
    S = set(S)
    for q in range(2, n+1):
        if not any(v % q == 0 for v in S):
            S.add(q * (1 + rng.r() % 3))
    S = sorted(S)
    while len(S) > n-1:
        # drop a redundant element (covering preserved)
        for i in range(len(S)):
            T = S[:i]+S[i+1:]
            if covers(T, n):
                S = T; break
        else:
            break
    while len(S) < n-1:
        v = 1 + rng.r() % vmax
        if v not in S:
            S.append(v); S = sorted(set(S))
    return S if (len(S) == n-1 and covers(S, n) and prim(S)) else None


def search(n, rng, restarts, steps, vmax):
    best = (2.0, None)
    # structured seeds
    seeds = []
    for m in range(2*n-1, 4*n, 1):
        QR = sorted(set((x*x) % m for x in range(1, m)) - {0})
        seeds.append([q if q > 0 else m for q in QR][:n-1])
    seeds.append(list(range(1, n)))               # consecutive
    seeds.append([2*v for v in range(1, n)])       # even block (non-prim, will be repaired)
    for _ in range(restarts):
        if seeds:
            s0 = seeds.pop()
        else:
            s0 = [1+rng.r() % vmax for _ in range(n-1)]
        cur = fix_to_covering(s0, n, rng, vmax)
        if cur is None:
            continue
        curM = Mfloat(cur)
        for _ in range(steps):
            i = rng.r() % len(cur)
            new = 1 + rng.r() % vmax
            cand = sorted(set(cur[:i]+[new]+cur[i+1:]))
            if len(cand) != n-1 or not covers(cand, n) or not prim(cand):
                continue
            m = Mfloat(cand)
            if m < curM - 1e-12 or (m < curM + 0.002 and rng.r() % 4 == 0):
                cur, curM = cand, m
                if curM < best[0]:
                    best = (curM, cur[:])
    return best


def main():
    print("SMART covering-min search: margin deviation + odd-n winners\n")
    print(f"  {'n':>3} {'1/n':>9} {'best M (float)':>14} {'best M exact':>14} {'margin':>10} {'mod':>5}  set")
    rng = LCG(20260630)
    for n in range(7, 15):
        bM, bS = search(n, rng, restarts=120, steps=300, vmax=5*n)
        if bS is None:
            print(f"  {n:>3}  (no set found)")
            continue
        em, et = Mexact(bS)
        mod = et.denominator if et else 0
        print(f"  {n:>3} {1/n:>9.5f} {bM:>14.5f} {str(em):>14} {float(em)-1/n:>10.5f} {mod:>5}  {bS}")


if __name__ == "__main__":
    main()
