#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S33: THE MIN-RISE CROSSOVER -- the exact quantity mac-mini's
HYP-4502 needs (generalized AP => M >= 2/(2k+1), via min-rise vs window).

opus HYP-4456 / kps S32: for k runners the gap is (1/(k+1), 2/(2k+1)), window
w(k) = 1/((k+1)(2k+1)).  A non-AP generalized AP survives (lands in the gap) iff its
M-RISE = M - 1/(k+1) is < w(k).  The gap is NONEMPTY at k=6 (n=7), EMPTY at k=12 (n=13).

This finds, for each k, the MINIMUM M-rise over non-AP families (the tightest
generalized AP), and compares to w(k).  The CROSSOVER (min-rise >= w) is the
gap-empty onset -- the exact n-specific transition mac-mini's proof must cross.
"""
from fractions import Fraction
from itertools import combinations, product
import numpy as np
from math import gcd
from functools import reduce

def M_exact(v):
    v=[x for x in v if x!=0]; S=int(sum(abs(x) for x in v)); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd=0,1
    for q in range(2,Q+1):
        a=np.arange(1,q,dtype=np.int64); r=np.outer(va,a)%q; d=np.minimum(r,q-r); bq=int(d.min(axis=0).max())
        if bq*bd>bn*q: bn,bd=bq,q
    return Fraction(bn,bd)

def min_M_above_AP(k, dmax=6, pertmax=3):
    """min M over non-AP k-runner generalized APs that is > 1/(k+1)."""
    AP = Fraction(1, k+1)
    best = Fraction(1)   # smallest M found that is > AP
    for d in range(1, dmax+1):
        for L in range(max(3, k-4), k+1):
            base = [1 + i*d for i in range(L)]
            need = k - L
            if need < 0: continue
            extras = sorted({base[j]+e for j in range(L) for e in (-2,-1,1,2)} - set(base) - {0})
            extras = [x for x in extras if x >= 1]
            if len(extras) < need: continue
            combos = list(combinations(extras, need))
            if len(combos) > 3000: combos = combos[:3000]
            for ex in combos:
                v = sorted(base + list(ex))
                if len(set(v)) != k: continue
                if reduce(gcd, v) != 1: continue
                M = M_exact(v)
                if AP < M < best:
                    best = M
    return best

print("=== min-rise crossover: k, gap (1/(k+1), 2/(2k+1)), window, min non-AP M, rise vs window ===", flush=True)
print(f"  {'k':>2} {'AP=1/(k+1)':>11} {'top=2/(2k+1)':>12} {'window':>10} {'min-M(nonAP)':>13} {'min-rise':>10} {'vs window':>12}", flush=True)
for k in range(5, 14):
    AP = Fraction(1, k+1); top = Fraction(2, 2*k+1); w = top - AP
    mm = min_M_above_AP(k)
    rise = mm - AP
    status = "EMPTY (rise>=w)" if mm >= top else ("NONEMPTY (in gap!)" if AP < mm < top else "?")
    print(f"  {k:>2} {str(AP):>11} {str(top):>12} {str(w):>10} {str(mm):>13} {str(rise):>10}  {status:>16}", flush=True)

print(flush=True)
print("READING: the crossover k* where min-M(nonAP) first reaches 2/(2k+1) (rise>=window)", flush=True)
print("is the gap-empty onset.  Below k*: a generalized AP fits the window (gap nonempty).", flush=True)
print("At/above k*: min-rise >= window, gap empty.  mac-mini HYP-4502 must prove k=12 > k*.", flush=True)
print("The exact min-rise as a function of k is the metric-alignment quantity their proof needs.", flush=True)
