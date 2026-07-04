#!/usr/bin/env python3
"""
Extend THM-618: is {1..12} the TIGHTEST base? For any 12-runner base B covering {2..12} + single killer 182,
is M(B u {182}) >= 14/183? If yes => ALL single-killer families >= 14/183 (whole stratum). (mac-mini-S44)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np, random
def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 8000)/8000.0
def approxM(sp):
    v = np.array(sp, float); ph = np.outer(v, _G) % 1.0
    return np.minimum(ph, 1.0-ph).min(axis=0).max()
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
def covers_range(sp, lo, hi): return all(any(v % q == 0 for v in sp) for q in range(lo, hi+1))
if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    cmin = F(14,183); rng = random.Random(44)
    out(f"single-killer families (12-runner base covers {{2..12}} + killer 182): M >= 14/183={float(cmin):.6f}?")
    tested = below = 0; minM = 1.0; minfam = None
    for _ in range(20000):
        B = sorted(set(rng.sample(range(1, rng.choice([13,20,36])), 12)))
        if len(B) != 12 or not covers_range(B, 2, 12): continue
        S = B + [182]
        if len(set(S)) != 13 or gcd_all(S) != 1: continue
        tested += 1
        aM = approxM(S)
        if aM < minM: minM = aM; minfam = B
        if aM < float(cmin) - 5e-4:  # float-flag possible below; confirm exact
            if M_exact(S) < cmin: below += 1; out(f"   BELOW 14/183: {B}+[182] M={M_exact(S)}")
    out(f"tested {tested}; families flagged below 14/183 (exact-confirmed): {below}; min approxM = {minM:.6f}")
    out(f"deep well {{1..12,182}}: M = {float(M_exact(list(range(1,13))+[182])):.6f}")
    out("=> if 0 below and min ~14/183: {1..12} is the TIGHTEST base => WHOLE single-killer stratum >= 14/183.")
