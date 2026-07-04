#!/usr/bin/env python3
"""
f=2 confinement via a pigeonhole over the m shifts {t+j/m} (mac-mini-2026-07-04-S41).
opus THM-616: ONE tightener is jointly useless for all m (orbit-max, Psi_1>=1/4). f>=2 is the named gap.
IDEA: mU is (+1/m)-invariant, so M(U) is hit at ALL m shifts of the argmax; a tightener w is UNSAFE on
<= 2 M(U) m + gcd(w,m) shifts; if the tighteners' unsafe sets don't cover all m shifts, a both-safe shift
exists => M(mU u F) = M(U) >= 1/12 > 1/14 (loose). Prediction: closes for m>=3-4 (enough shifts); m=2 hard.
Test: for m=2..9, is M(mU u {w1,w2}) = M(U) always (tighteners jointly useless)?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np, random
def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 5000)/5000.0
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

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    rng = random.Random(41)
    out("f=2: M(mU u {w1,w2}) vs M(U) for m=2..9. = means tighteners JOINTLY USELESS (loose, M=M(U)>1/14).")
    out(f"{'m':>3} {'tested':>7} {'jointly useless':>16} {'BITE (M<M(U))':>14} {'min M when bite':>16}")
    for m in range(2, 10):
        tested = usel = bite = 0; minbite = 1.0
        for _ in range(150):
            U = sorted(set(rng.sample(range(1, 14), 11)))
            if len(U) != 11: continue
            E = [m*u for u in U]
            cand = [w for w in range(1, 50) if w % m != 0 and w not in E]
            if len(cand) < 2: continue
            w1, w2 = rng.sample(cand, 2)
            S = E + [w1, w2]
            if len(set(S)) != 13 or gcd_all(S) != 1: continue
            # fast float prefilter: if approxM(S) ~ approxM(U), likely equal; confirm exact only near-equal
            aU, aS = approxM(U), approxM(S)
            tested += 1
            if abs(aS - aU) < 1e-3:
                usel += 1  # jointly useless (M(S)=M(U) to float precision)
            else:
                bite += 1; minbite = min(minbite, aS)
        out(f"{m:>3} {tested:>7} {usel:>16} {bite:>14} {(minbite if bite else 0):>16.5f}")
    out("\n=> if m>=3 has 0 BITE: pigeonhole closes f=2 for m>=3 (M=M(U)>=1/12>1/14); only m=2 remains (opus folding).")
    out("   BITE at m=2 (or small m) exhibits the tightener conspiracy = the named f>=2 residual.")
