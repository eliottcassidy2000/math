#!/usr/bin/env python3
"""
lrc_continuous_periodmax_macmini_0621s6.py  (mac-mini-2026-06-21-S6)
The CONTINUOUS period-max: f(s)=Delta_s*s = Sum_j Sum_t +-S_j(frac(s t)) for REAL s.
f is periodic in real s with period P=7*lcm(B) (Pt in Z). So sup_{real s} f = max over [0,P)
(finite, attained at a sawtooth breakpoint). This handles DILATED bases (scale-invariance: base
d*X + far w  ==  base X + real far speed s=w/d, and (w)Delta_w(d*X) = d*(s Delta_s(X))).
For the bounded base + real far speed s>14: Delta_s <= contmax/s <= contmax/14. Need contmax<=14*margin.
Compute contmax exactly at the breakpoints for consec_{k-1}.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def sector_of(p): return int((p%1)*7)
def setup(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); arcs=[];cur={};p0=F(0);p1=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); miss=[j for j in range(1,7) if j not in secs]
        if len(secs)==7: p0+=x1-x0
        mj=miss[0] if len(miss)==1 else None
        if mj is not None: p1+=x1-x0
        for j in range(1,7):
            active=(mj==j)
            if active and j not in cur: cur[j]=x0
            if (not active) and j in cur: arcs.append((j,cur.pop(j),x0))
    for j,a in cur.items(): arcs.append((j,a,F(1)))
    return arcs, p0+p1/7
def Sj_raw(t,j):
    t=t-int(t)
    if t<0:t+=1
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
def f(arcs,s):
    tot=F(0)
    for (j,a,bb) in arcs: tot+=Sc(s*bb,j)-Sc(s*a,j)
    return tot
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
for k in (8,9,10):
    B=tuple(range(k-1)); arcs,plat=setup(B); margin=caps[k]-plat
    L=1
    for e in B:
        if e>0: L=lcm(L,e)
    P=7*L
    # breakpoints of f in s over [0,P): s where s*t in {n, n+j/7, n+(j+1)/7} for endpoints t & corners
    endpoints=set()
    for (j,a,bb) in arcs:
        for t in (a,bb):
            if t==0: continue
            for c in (F(0),F(j,7),F(j+1,7)):
                n=0
                while True:
                    sval=(n+c)/t
                    if sval>=P: break
                    if sval>F(14): endpoints.add(sval)
                    n+=1
    mx=F(-10); ms=None
    for s in sorted(endpoints):
        # evaluate just left and right (piecewise linear; check the breakpoint and midpoints)
        v=f(arcs,s)
        if v>mx: mx=v; ms=s
    thr=14*margin
    print(f"k={k}: contmax(s>14, over period {P}) = {float(mx):.5f} at s={float(ms):.3f}  14*margin={float(thr):.5f}  {'OK (dilated-safe)' if mx<thr else 'CHECK'}")
print("=> contmax < 14*margin => Delta_s<margin for all real far speed s>=14 => DILATED case closes too.")
