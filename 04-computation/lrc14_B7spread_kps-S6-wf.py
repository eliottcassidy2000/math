#!/usr/bin/env python3
"""
lrc14_B7spread_kps-S6-wf.py
Test whether the concurrent session's 7-arc NET-MINORANT B_7(E) (a RIGOROUS lower
bound on mu_{1/7}(E)) admits a provable SPREAD behaviour:  does B_7(E) stay >= thr_8
as spread grows?  If B_7 is large for large spread (equidistribution) and we already
have it >= thr_8 for bounded spread (exhaustive), then B_7 >= thr_8 UNIFORMLY, which
PROVES mu_{1/7}(E) >= thr_8 for ALL E (k=8) -> closes the binding case.

B_7(E) = meas{x: SOME of 7 fixed arcs [c_i, c_i+1/7), c_i=(2i+1)/14, is empty of
              {frac(e x): e in E}}.   (empty length-1/7 arc => maxgap>=1/7 => in mu a.e.)

We compute the MIN of B_7 over primitive k=8 E by spread, pushing spread high, to see
if the min stays >= thr_8 and whether it is monotone / bounded below away from thr_8.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(123)
from math import gcd
from functools import reduce
W=F(1,7)
THR8=F(3637,5880)
CENTERS=[F(2*i+1,14) for i in range(7)]  # 7 arcs at (2i+1)/14

def union_arc_measure(E, centers=CENTERS, w=W):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for c in centers:
            j=0
            while True:
                lo=(c+j)/e; hi=(c+w+j)/e
                if lo>=1: break
                bps.add(max(F(0),lo)); bps.add(min(F(1),hi)); j+=1
                if j>e+2: break
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid=(a+b)/2; fracs=[(e*mid)%1 for e in E]; good=False
        for c in centers:
            empty=True
            for fr in fracs:
                if c<=fr<c+w: empty=False; break
            if empty: good=True; break
        if good: tot+=b-a
    return tot

if __name__=="__main__":
    print(f"thr_8={THR8}={float(THR8):.6f}; 7 fixed arcs at (2i+1)/14")
    print("B_7(E) is a RIGOROUS lower bound on mu_{1/7}(E).  Min over primitive k=8 by spread:")
    gmin=F(10); gminE=None
    for Wmax in range(7,16):
        best=F(10); bestE=None; cnt=0
        for body in itertools.combinations(range(1,Wmax+1),7):
            if body[-1]!=Wmax: continue
            E=(0,)+body
            if reduce(gcd,E)!=1: continue
            cnt+=1
            m=union_arc_measure(list(E))
            if m<best: best=m; bestE=E
        if best<gmin: gmin=best; gminE=bestE
        print(f"  spread={Wmax:2d}: min B_7={float(best):.6f} ({best})  E={bestE}  "
              f"[{cnt} sets] {'>=thr_8' if best>=THR8 else '*** < thr_8 ***'}", flush=True)
    # random large spread
    print("\n  random primitive large-spread min B_7:")
    for sp in [20,30,50,80,120,200]:
        best=F(10); bestE=None
        for _ in range(2000):
            E=sorted(set([0]+random.sample(range(1,sp+1),7)))
            if len(E)<8: continue
            if reduce(gcd,E)!=1: continue
            m=union_arc_measure(E)
            if m<best: best=m; bestE=E
        if best<gmin: gmin=best; gminE=bestE
        print(f"  spread<={sp:3d}: min B_7={float(best):.6f} ({best})  E={bestE}  "
              f"{'>=thr_8' if best>=THR8 else '*** < thr_8 ***'}", flush=True)
    print(f"\n  GLOBAL min B_7 found = {float(gmin):.6f} at {gminE}; thr_8={float(THR8):.6f}; "
          f"{'B_7>=thr_8 throughout' if gmin>=THR8 else 'B_7 dipped below thr_8'}")
