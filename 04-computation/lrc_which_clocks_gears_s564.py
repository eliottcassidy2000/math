#!/usr/bin/env python3
"""Which clocks matter: worry/ignore separation + the prime-gear n-clock.
opus-2026-06-02-S564. IGNORE = positive safe-measure (orbit spends an interval in
B, hits trivially; ~all sets, all incommensurate, all low-resonance). WORRY =
measure-0 (resonance-maximal tight family). For WORRY the only clock is t=j/n,
which factors via CRT into prime gears (mod 2 x mod 7 for n=14): a runner is at the
observer iff it aligns on EVERY gear (is a multiple of n)."""
from fractions import Fraction
from math import gcd
import random
n = 14; thr = Fraction(1, n)

def positive_measure(V):
    eps = set()
    for v in V:
        for k in range(v + 1):
            for s in (-1, 1): eps.add(Fraction(k*n+s, n*v) % 1)
    pts = sorted(eps); L = len(pts)
    for i in range(L):
        a, b = pts[i], pts[(i+1) % L]; ln = (b-a) if b > a else (b-a+1)
        mid = (a + ln/2) % 1
        if all(min((v*mid) % 1, 1-(v*mid) % 1) > thr for v in V): return True
    return False

def n_gear(V):
    for j in range(1, n):
        if gcd(j, n) == 1 and all((v*j) % n != 0 for v in V): return j
    return None

def prim(V):
    g = 0
    for v in V: g = gcd(g, v)
    return tuple(sorted(v//g for v in V))

if __name__ == "__main__":
    rng = random.Random(1)
    samp = [prim(tuple(rng.sample(range(1, 60), 13))) for _ in range(300)]
    samp += [tuple(range(1, 14)), (1,2,3,4,5,6,7,8,9,10,11,13,24)]
    ig = wo = ng = 0
    for V in samp:
        if positive_measure(V): ig += 1
        else:
            wo += 1; ng += 1 if n_gear(V) else 0
    print(f"{len(samp)} configs: IGNORE(pos measure)={ig}, WORRY(measure 0)={wo}; "
          f"worry caught by n-gear={ng}/{wo}")
    V = tuple(range(1, 14))
    for j in [1, 3, 5, 9, 11, 13]:
        bad = [v for v in V if (v*j) % 14 == 0]
        print(f"  AP gear j={j}: doubly-aligned={bad} -> {'CLEAR' if not bad else 'blocked'}")
