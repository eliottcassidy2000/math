#!/usr/bin/env python3
"""period-max(B) <= 15*margin(B) for the DANGEROUS (high-Plat, near-consec) bounded bases B subset [0,14]."""
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
    t=t-int(t); 
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
def period_max(arcs,B):
    L=1
    for e in B:
        if e>0: L=lcm(L,e)
    P=7*L
    if P>30000: return None,P   # skip huge-period (low-Plat, easy) bases
    mx=F(-10)
    for w in range(15,15+P):
        tot=F(0)
        for (j,a,bb) in arcs: tot+=Sc(w*bb,j)-Sc(w*a,j)
        if tot>mx: mx=tot
    return mx,P
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}
for k in (8,9):
    cap=caps[k]; rows=[]
    for combo in itertools.combinations(range(1,15),k-2):
        B=(0,)+combo; plat,arcs=plat_and_arcs(B); rows.append((plat,B,arcs))
    rows.sort(reverse=True)  # by Plat (most dangerous first)
    print(f"\nk={k}: top-20 highest-Plat bounded bases, period-max(B) vs 15*margin(B):")
    worst_ratio=F(0); fails=0; checked=0
    for plat,B,arcs in rows[:20]:
        margin=cap-plat; pm,P=period_max(arcs,B)
        if pm is None: 
            print(f"  B={B} P={P} (skipped, huge period)"); continue
        checked+=1; thr=15*margin; ratio=pm/margin if margin>0 else F(999)
        if pm>=thr: fails+=1
        if ratio>worst_ratio: worst_ratio=ratio
        print(f"  B={str(B):<26} Plat={float(plat):.4f} margin={float(margin):.4f} P={P} period-max={float(pm):.4f} 15m={float(thr):.4f} {'CLOSE' if pm<thr else 'FAIL'}")
    print(f"  => checked {checked}, fails={fails}, worst period-max/margin = {float(worst_ratio):.3f} (need <15)")
