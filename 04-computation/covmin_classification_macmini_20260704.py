#!/usr/bin/env python3
"""
Map the covering-min residual (mac-mini-2026-07-04-S45): classify covering families by which theorem covers
them, verify M>=14/183 holistically, identify the EXACT residual toward closing the proof.
 THM-617 (shift-pigeonhole): m>=3 large-scale, f<=6. THM-618: single-killer {1..12}-base ladder.
 Report: any M<14/183? which families are RESIDUAL (not covered by 617/618)?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import numpy as np, random
def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 9000)/9000.0
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
def is_cov(sp, n): return all(any(v % q == 0 for v in sp) for q in range(2, n+1))
def scale_m(sp):
    """largest m such that the m-divisible runners = m*U (i.e. all m-divisible v have v/m forming U); return the
       max m>=2 with |m-divisible runners| = 13 - f, f = # non-m-divisible; here just report min f over m>=2."""
    best=(99,1)
    for m in range(2,15):
        E=[v for v in sp if v%m==0]; f=len(sp)-len(E)
        if 2<=len(E) and f<best[0]: best=(f,m)  # f tighteners at scale m
    return best  # (f, m)
if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    cmin=F(14,183); rng=random.Random(45)
    n=14; tested=below=residual=cov617=cov618=0; minM=1.0; worst=None
    for _ in range(30000):
        hi=rng.choice([15,18,25,40,190])
        sp=sorted(set(rng.sample(range(1,hi),13)))
        if len(sp)!=13 or gcd_all(sp)!=1 or not is_cov(sp,n): continue
        tested+=1
        aM=approxM(sp)
        if aM<minM: minM=aM; worst=sp
        # classify
        f,m=scale_m(sp)
        covered_617 = (m>=3 and f<=6 and f*2*(14/183)*m + sum(gcd(w,m) for w in sp if w%m!=0) < m)
        base=[v for v in sp if v not in (182,) and v<=12]  # rough single-killer {1..12} check
        covered_618 = (182 in sp and set(range(1,13)).issubset(set(sp)))
        if covered_617: cov617+=1
        elif covered_618: cov618+=1
        else: residual+=1
        if aM < float(cmin)-5e-4:
            if M_exact(sp)<cmin: below+=1; out(f"   !! M<14/183: {sp} M={M_exact(sp)}")
    out(f"tested {tested} primitive covering families; M<14/183 (exact): {below}; min approxM={minM:.6f} (14/183={float(cmin):.6f})")
    out(f"  covered by THM-617 (large-scale m>=3): {cov617}")
    out(f"  covered by THM-618 (single-killer {{1..12}}): {cov618}")
    out(f"  RESIDUAL (neither): {residual}  ({100*residual//max(tested,1)}%)")
    out(f"  min-M family: {worst}")
    out("=> the residual = the strata still needing closure (small-m non-single-killer, split-killer, other bases).")
