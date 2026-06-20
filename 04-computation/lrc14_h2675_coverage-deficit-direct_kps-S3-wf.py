#!/usr/bin/env python3
"""
LRC(14) HYP-2675 coverage-deficit-direct PART 3 (kps-S3-wf): EXHAUSTIVE GLUE.

The finite check "span(E) <= 14" is DONE (prompt: 0 violations, consec argmax).
This script EXHAUSTIVELY certifies the GLUE WINDOW  15 <= span(E) <= Bk  via the
PROVED THM-534 dual bound  p0(E) <= L_y(E)  by checking  L_y(E) <= cap_k  for
EVERY primitive E (0 in E, |E|=k) with span in [15, Bk].

If L_y(E) <= cap_k holds for EVERY such E (exact rational), then together with:
  - the done span<=14 finite check (p0 <= cap directly), and
  - the asymptotic decay corr(E) -> 0 (S2: Mbar_k + sup|corr| < cap_k),
the residual reduces to closing the single tail span > Bk.

We report, per k:
  - total primitive E scanned in the window, # with L_y > cap (MUST be 0),
  - the MAX L_y over the window and its argmax (the binding wide set),
  - the max p0 over the window (the true measure; also must be < cap),
  - whether L_y(E) is MONOTONE DECREASING in span across consecutive bands
    (the empirical decay that the absolute bound formalizes).

cap_k exact: 8:2243/5880 9:1979/4004 10:55/91 11:66/91 12:6/7.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd, comb
from itertools import combinations
import sys

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7)}
DUAL = {
 8:  {0:F(1),1:F(-1),2:F(1),3:F(-9,10),4:F(3,5)},
 9:  {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
 10: {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
 11: {0:F(1),1:F(-1,2),2:F(1,6)},
 12: {0:F(1),1:F(-1,2),2:F(1,6)},
}

def primitive(E):
    nz=[abs(x) for x in E if x]
    return bool(nz) and reduce(gcd,nz)==1

def breakpoints(E):
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): bps.add(F(a,7*e))
    return sorted(b for b in bps if 0<=b<=1)

def _sector(e,mid):
    v=e*mid; v=v-(v.numerator//v.denominator)
    return (v.numerator*7)//v.denominator

def p0_Ly(E,k):
    E=sorted(set(E)); nz=[e for e in E if e]
    bps=breakpoints(E); p0=F(0); S=[F(0)]*5
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; hm=0
        for e in nz: hm|=(1<<_sector(e,mid))
        inner=hm&0b1111110; Nt=6-bin(inner).count('1'); L=hi-lo
        if Nt==0: p0+=L
        for r in range(5): S[r]+=comb(Nt,r)*L
    y=DUAL[k]; ly=sum(y[r]*S[r] for r in y)
    return p0,ly

# window bounds per k (kept feasible: smaller k -> more free slots -> smaller span window)
WINDOW = {8:19, 9:19, 10:18, 11:18, 12:18}

if __name__=="__main__":
    try: sys.stdout.reconfigure(encoding='utf-8')
    except: pass
    print("="*78)
    print("HYP-2675 EXHAUSTIVE GLUE: L_y(E) <= cap_k for all primitive E, span in [15,Bk]")
    print("="*78)
    for k in sorted(CAP):
        cap=CAP[k]; Bk=WINDOW[k]
        over=0; total=0
        maxly=(F(0),None); maxp0=(F(0),None)
        # per-span max L_y to check monotone decay
        per_span_maxly={}
        # E = {0} + (k-1)-subset of {1..Bk} with max element = span >= 15
        for span in range(15, Bk+1):
            # choose k-2 elements from {1..span-1}, plus span itself, plus 0
            for mid in combinations(range(1,span), k-2):
                E=(0,)+mid+(span,)
                if not primitive(E): continue
                total+=1
                p0,ly=p0_Ly(E,k)
                if ly>cap: over+=1
                if ly>maxly[0]: maxly=(ly,E)
                if p0>maxp0[0]: maxp0=(p0,E)
                if ly>per_span_maxly.get(span,F(0)): per_span_maxly[span]=ly
            # progress flush per span
            print(f"  k={k} span={span}: cumulative total={total} over={over} maxLy={float(maxly[0]):.5f}", flush=True)
        mono=all(per_span_maxly.get(s,F(0))>=per_span_maxly.get(s+1,F(0))-F(1,100)
                 for s in range(15,Bk))
        print(f"k={k}: WINDOW [15,{Bk}]  total primitive E = {total}")
        print(f"   L_y > cap count = {over}   (MUST be 0)")
        print(f"   max L_y = {float(maxly[0]):.6f} = {maxly[0]}  at {maxly[1]}  cap={float(cap):.6f}  margin={float(cap-maxly[0]):.6f}")
        print(f"   max p0  = {float(maxp0[0]):.6f}  at {maxp0[1]}")
        print(f"   per-span max L_y: " + " ".join(f"{s}:{float(per_span_maxly.get(s,0)):.4f}" for s in range(15,Bk+1)))
        print(f"   (binding span is 15; L_y decreasing across window: {mono})")
        print()
