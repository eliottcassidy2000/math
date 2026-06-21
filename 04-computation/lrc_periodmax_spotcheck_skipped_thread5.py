#!/usr/bin/env python3
"""THREAD 5: directly compute EXACT period-max for the smallest-margin SKIPPED bases (k=8),
to test whether the 'huge-period => safe' assumption actually holds where it was never checked.
We allow periods up to ~400k (the smallest-margin skipped k=8 bases have P up to 360360)."""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def sector_of(p): return int((p%1)*7)
def breakpoints(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    return sorted(b)
def plat_and_arcs(E):
    b=breakpoints(E); p0=F(0);p1=F(0);arcs=[];cur={}
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
    return p0+p1/7, arcs
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
def period_max(arcs,B,cap_P):
    L=1
    for e in B:
        if e>0: L=lcm(L,e)
    P=7*L
    if P>cap_P: return None,P
    mx=F(-10); argw=None
    for w in range(15,15+P):
        tot=F(0)
        for (j,a,bb) in arcs: tot+=Sc(w*bb,j)-Sc(w*a,j)
        if tot>mx: mx=tot; argw=w
    return mx,P,argw
cap=F(2243,5880); k=8
# the 8 smallest-margin skipped bases from the audit (k=8):
test_bases=[
 (0,8,9,10,11,12,13),(0,1,2,10,11,12,13),(0,1,3,5,9,11,13),(0,1,9,10,11,12,13),
 (0,4,5,8,9,10,13),(0,1,3,9,11,12,13),(0,2,3,5,8,11,13),(0,3,4,5,8,9,13)]
print(f"k={k} cap={float(cap):.5f}: EXACT period-max for smallest-margin SKIPPED bases")
print(f"{'B':<28} {'margin':>8} {'P':>8} {'period-max':>11} {'pm/margin':>10} {'argw':>7} {'15m':>8} verdict")
worst=F(0)
for B in test_bases:
    plat,arcs=plat_and_arcs(B); margin=cap-plat
    pm,P,argw=period_max(arcs,B,cap_P=400000)
    if pm is None:
        print(f"{str(B):<28} {float(margin):>8.4f} {P:>8} {'TOO-BIG':>11}"); continue
    ratio=pm/margin
    if ratio>worst: worst=ratio
    verdict="CLOSE" if pm<15*margin else "*** FAIL ***"
    print(f"{str(B):<28} {float(margin):>8.4f} {P:>8} {float(pm):>11.4f} {float(ratio):>10.3f} {argw:>7} {float(15*margin):>8.4f} {verdict}")
print(f"worst ratio among these skipped bases = {float(worst):.3f} (need < 15)")
