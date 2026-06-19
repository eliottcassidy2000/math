#!/usr/bin/env python3
"""
lrc14_netspread_kps-S6-wf.py
Probe the SPREAD BOUND for the net measure.  For primitive (gcd=1) k=8 sets,
scan increasing max(E) and confirm max net stays at 44/735 (= consec).  We also
test a candidate provable spread bound:
    net(E) <= sum over the SMALLEST consecutive difference structure ...
Specifically, candidate: net(E) <= 1/(largest gap of E in [0,max])? no.

Real goal here: find the SPREAD CUTOFF S(k) beyond which net(E) < consec, so the
proof = (scaling -> WLOG primitive) + (exhaustive primitive max(E)<=S(k)) +
(spread bound net<consec for max(E)>S(k)).

We measure: max net over primitive k=8 sets, bucketed by max(E), up to as large
as feasible; and the structural decay.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd
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
def netmeas(E): return meas(net_intervals(E))

if __name__=="__main__":
    consec=F(44,735)
    print(f"consec_8 net = {consec} = {float(consec):.6f}; cap_8=0.381463")
    print("\nMax net over PRIMITIVE k=8 sets, by spread W=max(E):")
    overall=-1; overallE=None
    for W in range(7,17):
        best=-1; bestE=None
        for body in itertools.combinations(range(1,W),6):
            E=(0,)+body+(W,)
            if reduce(gcd,E)!=1: continue
            nm=netmeas(list(E))
            if nm>best: best=nm; bestE=E
        if best>overall: overall=best; overallE=bestE
        flag='= consec' if best==consec else ('> CONSEC!!' if best>consec else '< consec')
        print(f"   W={W:2d}: max net = {float(best):.6f} ({best}) {flag}  argmax={bestE}")
    print(f"\n   overall max (W<=16) = {float(overall):.6f} ({overall}) at {overallE}")
    print(f"   {'CONFIRMS consec is the max up to W=16' if overall<=consec else 'CONSEC EXCEEDED'}")
