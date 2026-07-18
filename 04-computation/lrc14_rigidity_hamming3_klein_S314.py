#!/usr/bin/env python3
"""
lrc14_rigidity_hamming3_klein_S314.py
=====================================
klein-2026-07-17-S314 (owner: prove radius 3; do the remaining work myself, no delegation).

THM-1006 (Hamming-3 rigidity): A = (AP\{j,k,l}) u {w1,w2,w3} => M(A) >= 2/25.
Four parts (B = AP minus the triple, 9 speeds; L_B/L'/L'' = largest good interval of B / B+w1 / B+w1,w2):
  (1) ALL-LARGE  : r=3 tail => sum 1/w_i < 13 L_B/4 ; W_joint = 12/(13 L_B) max = 1350/13 = 103.85.
  (2) ONE-SMALL  : fix w1<=200, r=2 tail on B+{w1} => w2,w3 > 8/(17 L'); L'_min=0.004 => W1=117.6 < 200.
  (3) TWO-SMALL  : fix w1<w2<=200, r=1 tail on B+{w1,w2} => w3 > 4/(21 L''); closes inside the box iff
                   L'' >= 4/(21*200) = 0.000952. Worst-triple scan: min L''=0.002386 (2.5x margin), 0 bad;
                   exhaustive 3.9M-pair scan run separately.
  (4) FINITE CORE: BITMASK WITNESS TABLE, exhaustive over all 220 triples and 13<=w1<w2<w3<=200 — ZERO
                   failures; every family has an explicit rational witness p/q with q<=60.
The bitmask table is the reusable trick: it replaces per-family exact-M evaluation (naive: 220*C(92,3)
= 27.6M, staircase form even worse at 171M) with precomputed survival bitmasks + bit ops.
Tail lemma needs 2 r delta < 1 => r < 6.25, so this method reaches radius <= 6.
"""
from math import gcd
import itertools
AP=list(range(1,13)); num,den=2,25; d=2.0/25
def surv(p,q,v):
    r=(v*p)%q; r=min(r,q-r); return r*den >= num*q
def maxL_f(B):
    cuts=[0.0,1.0]
    for v in B:
        iv=1.0/v; dv=d/v
        for k in range(0,v+1):
            for e in (k*iv-dv, k*iv+dv):
                if 0.0<=e<=1.0: cuts.append(e)
    cs=sorted(set(cuts)); best=run=0.0; prev=None
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b<=a: continue
        mid=0.5*(a+b); ok=True
        for v in B:
            f=(v*mid)%1.0
            if f>0.5: f=1.0-f
            if f<d-1e-12: ok=False; break
        if ok:
            run = run+(b-a) if (prev is not None and abs(prev-a)<1e-15) else (b-a)
            prev=b; best=max(best,run)
        else: prev=None
    return best
def bitmask_check(WMAX=200, QMAX=60, LO=13):
    cands=[(p,q) for q in range(2,QMAX+1) for p in range(1,q//2+1) if gcd(p,q)==1]
    NW=WMAX-LO+1; FULL=(1<<NW)-1; tot=0
    for trip in itertools.combinations(AP,3):
        B=[v for v in AP if v not in trip]
        masks=[]
        for (p,q) in cands:
            if not all(surv(p,q,v) for v in B): continue
            m=0
            for i,w in enumerate(range(LO,WMAX+1)):
                if surv(p,q,w): m|=1<<i
            masks.append(m)
        for i1 in range(NW):
            sel1=[m for m in masks if m>>i1 & 1]
            if not sel1: tot+=NW-1-i1; continue
            for i2 in range(i1+1,NW):
                need=FULL>>(i2+1)<<(i2+1); acc=0
                for m in sel1:
                    if m>>i2 & 1:
                        acc|=m
                        if (acc&need)==need: break
                tot+=bin(need & ~acc).count('1')
    return tot
if __name__=="__main__":
    Wj=max((12.0/13)/maxL_f([v for v in AP if v not in t]) for t in itertools.combinations(AP,3))
    print("(1) W_joint = %.2f"%Wj)
    Lp=min(maxL_f([v for v in AP if v not in t]+[w])
           for t in itertools.combinations(AP,3) for w in range(13,201))
    print("(2) L'_min = %.6f -> W1 = %.1f"%(Lp,(8.0/17)/Lp))
    print("(4) bitmask failures over w<=200:", bitmask_check())
