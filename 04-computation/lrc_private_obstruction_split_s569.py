#!/usr/bin/env python3
"""Leverage the paired-or-anchored dichotomy (THM-396/397) to split LRC(n):
worry-set (measure-0) has an UNBLOCKED small pair (route 1); residual (block-all)
is positive-measure. Verified across even n. opus-2026-06-02-S569."""
from fractions import Fraction
from math import gcd
from itertools import combinations
import random

def dist(x):
    x = x % 1
    return min(x, 1 - x)

def run(n):
    THR = Fraction(1, n); m = n - 1
    def psp(a, b):
        D = a + b
        return [k for k in range(1, D)
                if dist(Fraction(a*k, D)) >= THR and dist(Fraction(b*k, D)) >= THR]
    def small_pairs(V):
        return [(a, b) for a, b in combinations(sorted(set(V)), 2)
                if (a + b) // gcd(a, b) <= n]
    def pair_unblocked(V, a, b):   # a pinch where ALL runners safe = a witness
        for k in psp(a, b):
            if all(dist(Fraction(c*k, a+b)) >= THR for c in V if c not in (a, b)):
                return True
        return False
    def has_unblocked_small(V):
        return any(pair_unblocked(V, a, b) for a, b in small_pairs(V))
    def positive_measure(V):
        eps = set()
        for v in V:
            for k in range(v+1):
                for s in (-1, 1):
                    eps.add(Fraction(k*n+s, n*v) % 1)
        pts = sorted(eps); L = len(pts)
        for i in range(L):
            a, b = pts[i], pts[(i+1) % L]; ln = (b-a) if b > a else (b-a+1)
            if all(dist(v*((a+ln/2) % 1)) > THR for v in V):
                return True
        return False
    def prim(V):
        g = 0
        for v in V:
            g = gcd(g, v)
        return tuple(sorted(v//g for v in V))
    B = {4: 10, 6: 11, 8: 13, 10: 14, 12: 15, 14: 16}[n]
    cands = [c for c in combinations(range(1, B+1), m) if prim(c) == tuple(c)]
    worry = worry1 = blockall = blockall_pos = 0
    for V in cands:
        pos = positive_measure(V); ub = has_unblocked_small(V)
        if not pos:
            worry += 1; worry1 += ub
        if not ub:
            blockall += 1; blockall_pos += pos
    print(f"n={n:2d}: WORRY(meas-0)={worry}, with unblocked small pair={worry1}/{worry}; "
          f"block-all={blockall}, positive-measure={blockall_pos}/{blockall}")

if __name__ == "__main__":
    for n in [4, 6, 8, 10, 12, 14]:
        run(n)
