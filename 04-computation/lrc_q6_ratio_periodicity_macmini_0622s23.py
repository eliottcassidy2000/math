#!/usr/bin/env python3
"""
lrc_q6_ratio_periodicity_macmini_0622s23.py  (mac-mini-2026-06-22-S23)
THE q6-RATIO BOUND via THM-563 periodicity (the flagged gK8-concentration gap).
q6(E)=meas{x: frac(e x) in [0,1/7) for all e in E} (all runners in sector 0). For E=B u {f}:
  q6(B u {f}) = int_{Q6(B)} 1_{[0,1/7)}(f x) dx,  Q6(B)={x: all of B in sector 0}.
Deviation delta_f = q6(B u {f}) - (1/7)q6(B) = int_{Q6(B)}[1_{[0,1/7)}(fx)-1/7]dx.
By the THM-563 identity: f*delta_f = Sum_{endpoints t of Q6(B)} +-S_0(frac(f t)), S_0=centered sawtooth
of (1_{[0,1/7)}-1/7). Endpoints of Q6(B) depend ONLY on B => f*delta_f periodic in f, period 7*lcm(B).
So q6-ratio: q6(B u {f})/q6(B) <= 1/7 + periodmax_q6(B)/(f*q6(B)). Verify the identity + the ratio < 1.
"""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def Q6_arcs(B):
    """arcs of Q6(B) = {x in [0,1): frac(b x) in [0,1/7) for all b in B}."""
    B=sorted(set(b for b in B if b>0)); bps=set([F(0),F(1)])
    for b in B:
        for m in range(0,b+1): bps.add(F(m,b)); bps.add(F(7*m+1,7*b))  # boundaries of [0,1/7) under b
    # cleaner: frac(b x)<1/7  <=> x in union over m of [m/b, m/b + 1/(7b))
    bps=set([F(0),F(1)])
    for b in B:
        for m in range(0,b):
            bps.add(F(m,b)); bps.add(F(m,b)+F(1,7*b))
    bps=sorted(z for z in bps if 0<=z<=1)
    arcs=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all((b*xm)%1 < F(1,7) for b in B): arcs.append((x0,x1))
    # merge adjacent
    merged=[]
    for a,b in arcs:
        if merged and merged[-1][1]==a: merged[-1]=(merged[-1][0],b)
        else: merged.append((a,b))
    return merged
def q6_direct(E):
    E=sorted(set(e for e in E if e>0)); 
    if not E: return F(1,7)  # only 0
    bps=set([F(0),F(1)])
    for e in E:
        for m in range(0,e): bps.add(F(m,e)); bps.add(F(m,e)+F(1,7*e))
    bps=sorted(z for z in bps if 0<=z<=1); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        if all((e*((x0+x1)/2))%1 < F(1,7) for e in E): tot+=x1-x0
    return tot
# centered sawtooth S_0 of (1_{[0,1/7)} - 1/7)
def S0(t):
    t=t-int(t)
    if t<0: t+=1
    if t<=F(1,7): return F(6,7)*t
    else: return F(6,7)*F(1,7)-F(1,7)*(t-F(1,7))
def meanS0():
    pts=[F(0),F(1,7),F(1)]; v=[S0(x) for x in pts]; I=F(0)
    for i in range(2): I+=(v[i]+v[i+1])/2*(pts[i+1]-pts[i])
    return I
M0=meanS0()
def Sc0(t): return S0(t)-M0
print("VERIFY q6 identity + ratio for B=consec_{k-1}:")
for k in (9,10):
    B=tuple(range(1,k))  # base elements 1..k-1 (drop 0; 0 always in sector 0)
    arcs=Q6_arcs(B); q6B=q6_direct(B)
    L=1
    for e in B: L=lcm(L,e)
    P=7*L
    # verify f*delta_f == sawtooth sum, a few f
    print(f"\n k={k} B(1..{k-1}): q6(B)={q6B}={float(q6B):.5f}  #arcs(Q6)={len(arcs)} period P={P}")
    for f in (15,20,50):
        q6Bf=q6_direct(B+(f,))
        delta=q6Bf-q6B*F(1,7)
        saw=sum(Sc0(f*b)-Sc0(f*a) for (a,b) in arcs)
        ratio=q6Bf/q6B
        print(f"   f={f}: q6(Bu f)={float(q6Bf):.5f} ratio={float(ratio):.4f} | f*delta={float(delta*f):.5f} sawtooth={float(saw):.5f} {'OK' if delta*f==saw else 'X'}")
    # period-max of f*delta over one period
    if P<=20000:
        mx=F(-10);mn=F(10)
        for f in range(15,15+P):
            v=sum(Sc0(f*b)-Sc0(f*a) for (a,b) in arcs)
            if v>mx: mx=v
            if v<mn: mn=v
        # ratio bound: q6(Bu f)/q6B = 1/7 + delta/q6B <= 1/7 + mx/(f q6B)
        worst_ratio_15 = F(1,7)+mx/(15*q6B)
        print(f"   period-max(f*delta)={float(mx):.4f} (min {float(mn):.4f}); worst q6-ratio at f>=15 <= 1/7+{float(mx):.3f}/(15*{float(q6B):.4f}) = {float(worst_ratio_15):.4f}  {'<1 OK (q6 strictly decreases)' if worst_ratio_15<1 else 'CHECK'}")
