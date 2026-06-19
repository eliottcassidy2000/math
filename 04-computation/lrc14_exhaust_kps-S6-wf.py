#!/usr/bin/env python3
"""
lrc14_exhaust_kps-S6-wf.py
Strongest available CERTIFIED statement for the binding case:
  EXACT exhaustion of net(E)=meas{maxgap<=1/7} over ALL primitive k=8 E with
  max(E)<=Wmax, confirming max net = 44/735 (=consec) < cap_8, with margin.
By the SCALING THEOREM (proved), every E reduces to a primitive one with the same
net, so this certifies ALL E with a primitive representative of spread<=Wmax.

Outputs the exact global max and its argmax per Wmax, and the running global max.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
ONE7=F(1,7)
CONSEC=F(44,735); CAP8=F(2243,5880)

def net_meas(E):
    E=sorted(set(E)); n=len(E); good=[]
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
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
    good.sort(); out=[]; tot=F(0)
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    for a,b in out: tot+=b-a
    return tot

if __name__=="__main__":
    print(f"EXACT exhaustion of net(E), primitive k=8.  consec=44/735={float(CONSEC):.6f}, cap_8={float(CAP8):.6f}")
    gmax=F(0); gmaxE=None
    import time
    for Wmax in range(7,20):
        best=F(0); bestE=None; cnt=0; t0=time.time()
        # all 8-subsets of {0..Wmax} with 0 and max=Wmax, primitive
        for body in itertools.combinations(range(1,Wmax+1),7):
            if body[-1]!=Wmax: continue
            E=(0,)+body
            if reduce(gcd,E)!=1: continue
            cnt+=1
            nm=net_meas(list(E))
            if nm>best: best=nm; bestE=E
        if best>gmax: gmax=best; gmaxE=bestE
        dt=time.time()-t0
        over = '*** > consec ***' if best>CONSEC else ''
        print(f"  Wmax={Wmax:2d}: max net (spread=Wmax) = {float(best):.6f} ({best})  "
              f"argmax={bestE}  [{cnt} sets, {dt:.1f}s] {over}", flush=True)
    print(f"\n  GLOBAL max net over ALL primitive k=8, spread<=19 = {float(gmax):.6f} ({gmax}) at {gmaxE}")
    print(f"  consec=44/735; cap_8=2243/5880")
    print(f"  RESULT: max net {'== consec (44/735)' if gmax==CONSEC else '!= consec'}, "
          f"{'<= cap_8' if gmax<=CAP8 else '> cap_8'}, margin to cap_8 = {float(CAP8-gmax):.6f}")
