#!/usr/bin/env python3
"""
lrc_dedekind_identity_macmini_0621s6.py  (mac-mini-2026-06-21-S6)
VERIFY the Dedekind identity for the signed one-far deviation:
  Delta_w * w = Sum_j Sum_{endpoints t of A_j} (+-) S_j(frac(w t)),
where A_j={x: B misses EXACTLY sector j}, t are its arc endpoints (rationals k/(7e)),
S_j = centered sawtooth antiderivative of (1_{sector_j} - 1/7) (period 1, |S_j|<=3/49).
Then DECOMPOSE by base element e and exhibit the cancellation that keeps Delta_w*w ~ 1 (bounded).
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def p1(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if 7-len(set(sector_of(e*((x0+x1)/2)) for e in E))==1: tot+=x1-x0
    return tot
def Aj_arcs(E):
    """arcs of A_j = {B misses exactly sector j}, j=1..6 (sector 0 always hit). Return list of (j,a,b)."""
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); arcs=[]; cur={}
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E)
        miss=[j for j in range(1,7) if j not in secs]
        mj=miss[0] if (len(miss)==1) else None
        # track contiguous A_j arcs
        for j in range(1,7):
            active=(mj==j)
            if active and j not in cur: cur[j]=x0
            if (not active) and j in cur: arcs.append((j,cur.pop(j),x0))
    for j,a in cur.items(): arcs.append((j,a,F(1)))
    return arcs
def S_j(t,j):
    """centered sawtooth: antiderivative of 1_{[j/7,(j+1)/7)}-1/7 over [0,1), value at frac(t)."""
    t=t-int(t); 
    if t<0: t+=1
    t=F(t) if not isinstance(t,F) else t
    a=F(j,7); b=F(j+1,7)
    if t<=a: return -t/7
    elif t<=b: return -a/7 + F(6,7)*(t-a)
    else: return -a/7 + F(6,7)*(b-a) - F(1,7)*(t-b)
# centering: subtract mean of S_j over [0,1]
def mean_Sj(j):
    a=F(j,7);b=F(j+1,7)
    # integral of piecewise-linear S_j (raw) over [0,1]
    pts=[F(0),a,b,F(1)]; vals=[S_j(x,j) for x in pts]
    I=F(0)
    for i in range(3): I+=(vals[i]+vals[i+1])/2*(pts[i+1]-pts[i])
    return I
MEAN=[None]+[mean_Sj(j) for j in range(1,7)]
def Sc(t,j): return S_j(t,j)-MEAN[j]

print(f"{'k':>3}{'w':>4}{'Delta_w*w direct':>18}{'sawtooth-sum':>16}{'match':>7}")
for k in (8,9):
    B=tuple(range(k-1)); Phi=measS7(B)+p1(B)/7; arcs=Aj_arcs(B)
    for w in (17,18,25,50):
        direct=(measS7(B+(w,))-Phi)*w
        saw=F(0)
        for (j,a,b) in arcs:
            saw += Sc(w*b,j) - Sc(w*a,j)
        print(f"{k:>3}{w:>4}{str(direct):>18}{str(saw):>16}{'OK' if direct==saw else 'X':>7}")
print("\n=> if match: Delta_w*w = signed sum of sawtooths S_j(w*k/(7e)) = generalized Dedekind sum (exact).")
