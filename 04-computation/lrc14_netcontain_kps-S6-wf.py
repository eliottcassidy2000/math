#!/usr/bin/env python3
"""
lrc14_netcontain_kps-S6-wf.py
Rigorous containment lemma candidates for the direct net bound.

CANDIDATE A:  For x in N(E), and for EVERY consecutive difference d=e_{i+1}-e_i
   of sorted E, we have  ||d x|| <= 1/7.
   Reason (to be proved): the two E-points e_i, e_{i+1} ... NO, not adjacent on circle.
   TEST whether this containment actually holds.

CANDIDATE B (the provable one):  For x in N(E), the SMALLEST positive point
   among {frac(e x): e in E\{0}} ... since 0 in E and its clockwise neighbor is
   within 1/7, there is e with frac(e x) in (0,1/7].  And reflection gives one in
   [6/7,1).  Test: N(E) subset { x : exists e in E, ||e x|| <= 1/7 } -- trivially
   true but measure too big.

CANDIDATE C (strong, provable):  Let d = gcd-witness: since gcd(E)=1, the
   differences generate Z, so there is a SET of differences D (a Z-basis-ish) with
   N(E) subset intersection_{d in D} {||d x|| <= (k-1)/7 * something}.

We TEST candidate A first (cleanest if true), then measure the intersection bound.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(11)
from math import gcd, floor
from functools import reduce
ONE7=F(1,7)

def net_intervals(E):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    good=[]
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(iv): return sum((b-a for a,b in iv),F(0))

def norm_band_meas(d, thr):
    """meas{ x in [0,1): ||d x|| <= thr } = 2*thr (exactly, for any positive int d)."""
    return min(F(1), 2*thr)

def max_norm_over_net(N, d):
    """max ||d x|| over x in N (exact via cell midpoints)."""
    mx=F(0)
    for a,b in N:
        pts=set([a,b])
        m0=floor(a*d*2); m1=floor(b*d*2)+1
        for m in range(m0,m1+1):
            xc=F(m,2*d)
            if a<xc<b: pts.add(xc)
        pts=sorted(pts)
        for u,v in zip(pts,pts[1:]):
            mid=(u+v)/2; dx=(d*mid)%1; nx=min(dx,1-dx)
            if nx>mx: mx=nx
    return mx

if __name__=="__main__":
    print("[CANDIDATE A] For x in N(E), does ||d x||<=1/7 for every CONSECUTIVE diff d?")
    viol=0; tot=0
    for trial in range(200):
        k=random.choice([8,9])
        E=sorted(set([0]+random.sample(range(1,16),k-1)))
        if len(E)<k: continue
        N=net_intervals(E)
        if not N: continue
        Es=sorted(E)
        for i in range(len(Es)-1):
            d=Es[i+1]-Es[i]
            mx=max_norm_over_net(N,d); tot+=1
            if mx>ONE7:
                viol+=1
                if viol<=5: print(f"   A-VIOLATION E={E} consec-d={d}: max||dx||={float(mx):.4f}>1/7")
    print(f"   candidate A: {tot} checks, {viol} violations => {'HOLDS' if viol==0 else 'FALSE'}")

    print("\n[CANDIDATE A'] does ||d x||<=1/7 hold for the SMALLEST diff d=min consecutive gap?")
    viol=0; tot=0
    for trial in range(300):
        k=random.choice([8,9,10])
        E=sorted(set([0]+random.sample(range(1,20),k-1)))
        if len(E)<k: continue
        N=net_intervals(E)
        if not N: continue
        Es=sorted(E); dmin=min(Es[i+1]-Es[i] for i in range(len(Es)-1))
        mx=max_norm_over_net(N,dmin); tot+=1
        if mx>ONE7:
            viol+=1
            if viol<=5: print(f"   A'-VIOLATION E={E} dmin={dmin}: max||dx||={float(mx):.4f}")
    print(f"   candidate A': {tot} checks, {viol} violations => {'HOLDS' if viol==0 else 'FALSE'}")

    print("\n[CANDIDATE D] ALL pairwise differences d: is max ||dx|| over N <= (some multiple)/7?")
    # Report, per random E, the max over consecutive diffs of (max||dx|| over N) in units of 1/7
    print("   sample: for each E, max_d (max_{x in N} ||dx||) / (1/7):")
    for trial in range(8):
        k=8
        E=sorted(set([0]+random.sample(range(1,16),k-1)))
        if len(E)<k: continue
        N=net_intervals(E)
        if not N: continue
        Es=sorted(E)
        ratios=[]
        for i in range(len(Es)-1):
            d=Es[i+1]-Es[i]; mx=max_norm_over_net(N,d); ratios.append(float(mx)/(1/7))
        print(f"   E={E} netmeas={float(meas(N)):.4f}: consec-diff ratios max||dx||/(1/7) = {[round(r,2) for r in ratios]}")
