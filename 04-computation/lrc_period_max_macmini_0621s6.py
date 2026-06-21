#!/usr/bin/env python3
"""
lrc_period_max_macmini_0621s6.py  (mac-mini-2026-06-21-S6)
THE PERIODICITY CRACK: Delta_w*w = Sum_j Sum_{endpoints t of A_j(B)} +-S_j(frac(w t)) (Dedekind id,
verified exact). Endpoints depend only on B => Delta_w*w periodic in w, period P=7*lcm(B).
Compute the EXACT max over one period via the FAST sawtooth sum (no measS7), for consec_{k-1},
and check max < 15*margin_k  =>  Delta_w < margin for all w>=15  =>  SINGLE-FAR CLOSED.
"""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def sector_of(p): return int((p%1)*7)
def Aj_arcs(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); arcs=[]; cur={}
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E)
        miss=[j for j in range(1,7) if j not in secs]; mj=miss[0] if len(miss)==1 else None
        for j in range(1,7):
            active=(mj==j)
            if active and j not in cur: cur[j]=x0
            if (not active) and j in cur: arcs.append((j,cur.pop(j),x0))
    for j,a in cur.items(): arcs.append((j,a,F(1)))
    return arcs
# centered sawtooth S_j(frac(t))
def Sj_raw(t,j):
    t=t-int(t)
    if t<0: t+=1
    a=F(j,7);b=F(j+1,7)
    if t<=a: return -t/7
    elif t<=b: return -a/7+F(6,7)*(t-a)
    else: return -a/7+F(6,7)*(b-a)-F(1,7)*(t-b)
def meanSj(j):
    a=F(j,7);b=F(j+1,7);pts=[F(0),a,b,F(1)];v=[Sj_raw(x,j) for x in pts];I=F(0)
    for i in range(3): I+=(v[i]+v[i+1])/2*(pts[i+1]-pts[i])
    return I
MEAN={j:meanSj(j) for j in range(1,7)}
def Sc(t,j): return Sj_raw(t,j)-MEAN[j]
def Dw_w(arcs,w):  # = Delta_w * w via sawtooth (fast, exact)
    tot=F(0)
    for (j,a,bb) in arcs: tot+=Sc(w*bb,j)-Sc(w*a,j)
    return tot
# margin_k = cap_k - Plat(consec_{k-1}); compute Plat exactly via measS7+p1 on the base
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
def p1f(E):
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
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
for k in (8,9,10):
    B=tuple(range(k-1)); arcs=Aj_arcs(B)
    L=1
    for e in B[1:]: L=lcm(L,e)
    P=7*L
    plat=measS7(B)+p1f(B)/7; margin=caps[k]-plat; thr=15*margin
    # verify periodicity: Dw_w(15)==Dw_w(15+P)
    per_ok = (Dw_w(arcs,15)==Dw_w(arcs,15+P))
    # max over one period [15,15+P)
    mx=F(-10); mw=0
    for w in range(15,15+P):
        v=Dw_w(arcs,w)
        if v>mx: mx=v; mw=w
    closes = mx<thr
    print(f"k={k}: period P={P}  periodic? {per_ok}  period-max(Delta_w*w)={float(mx):.5f}={mx} at w={mw}")
    print(f"      margin={float(margin):.5f}  threshold 15*margin={float(thr):.5f}  =>  max<thr? {closes}  {'SINGLE-FAR CLOSED for all w>=15' if closes else 'NOT closed'}")
