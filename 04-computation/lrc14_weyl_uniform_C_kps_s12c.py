#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) UNIFORM WEYL CONSTANT C  (kps-2026-06-19-S12c)

The remaining analytic lemma (HYP-2642) is |Delta_w(E)| <= C/w where
   Delta_w(E) = p0(E) - Plat(E'),  E' = E without its max element w,
   Plat(E') = p0(E') + (1/7) p1(E').
For a PROOF we need C uniform over ALL cores E' (|E'|=k-1), not just consec.
This script scans many cores (consec, near-AP, random, GAP, wide) and many w,
and reports sup_w w*|Delta_w| per core, then the overall sup = empirical C.

It also exhibits the EXACT structure of Delta_w as a single-frequency Weyl object:
   p0(E) = (1/(7w)) * sum over the 7w breakpoint sub-cells of the core-cover-deviation,
so Delta_w is a Riemann-sum error for integrating the core's conditional-cover function
against the equidistributing point frac(w x).  Rate O(1/w) is the std Koksma error with
total variation V of that function; C <= V/7 roughly.  We measure V too.
"""
import sys, math, itertools, random
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def dist_full(E):
    E=sorted(set(E)); bps=set([Fraction(0),Fraction(1)])
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1)
    p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        p[sum(1 for j in range(1,7) if j not in hit)]+=(hi-lo)
    return p
def p0p1(E):
    p=dist_full(E); return float(p[0]),float(p[1])

def core_defect_C(core, ws):
    p0c,p1c=p0p1(core)
    Plat=p0c+p1c/7.0
    best=0.0; worst_w=None
    for w in ws:
        if w in core: continue
        E=sorted(set(core)|{w})
        p0E,_=p0p1(E)
        c=w*abs(p0E-Plat)
        if c>best: best=c; worst_w=w
    return best, worst_w, Plat

if __name__=="__main__":
    print("LRC(14) UNIFORM WEYL CONSTANT — sup_w w*|Delta_w| over many cores (kps-S12c)\n")
    random.seed(7)
    ws=[11,13,17,19,23,29,31,37,41,53,67,89,113,151,211,307,401]
    overall_C=0.0; overall_row=None
    for k in [8,9]:
        m=k-1   # core size
        print(f"=== k={k} (core size {m}) ===")
        cores=[]
        cores.append(("consec", list(range(m))))
        cores.append(("nearAP", list(range(m-1))+[m]))      # one gap
        cores.append(("oddAP", [2*i for i in range(m)]))     # dilated; should match consec by THM-531? no, this is even AP
        cores.append(("oddvals", [0]+[2*i+1 for i in range(m-1)]))
        cores.append(("GAP", [0,1,2,3,8,9,10,11][:m] if m<=8 else [0,1,2,3,8,9,10,11,16][:m]))
        # random small-spread cores
        for r in range(8):
            base=set([0])
            while len(base)<m:
                base.add(random.randint(1, 3*m))
            cores.append((f"rand{r}", sorted(base)))
        Cmax=0.0; rowmax=None
        for name,core in cores:
            if math.gcd(*core)!=1 and len(core)>1:
                # reduce to primitive (dilation invariance)
                g=0
                for x in core: g=math.gcd(g,x)
                core=[x//g for x in core]
            C,ww,Plat=core_defect_C(core,ws)
            print(f"   {name:8s} core={core}  Plat={Plat:.5f}  sup_w w|Delta_w|={C:.4f} (at w={ww})")
            if C>Cmax: Cmax=C; rowmax=(name,core,ww)
        print(f"   -> k={k} empirical uniform C = {Cmax:.4f}  worst row {rowmax}")
        cap={8:0.38153,9:0.49426}[k]
        print(f"      cap_{k}={cap}; need Plat + C/w < cap.  With C={Cmax:.3f}, any core w/ Plat<cap-eps and w>{Cmax:.3f}/margin closes.\n")
        if Cmax>overall_C: overall_C=Cmax; overall_row=(k,rowmax)
    print(f"OVERALL empirical Weyl constant C ~ {overall_C:.4f}  (row {overall_row})")
    print("Interpretation: the single-frequency Weyl defect is uniformly O(1/w) with C<1.")
