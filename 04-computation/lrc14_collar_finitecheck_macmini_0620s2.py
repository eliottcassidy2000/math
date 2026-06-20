#!/usr/bin/env python3
"""
lrc14_collar_finitecheck_macmini_0620s2.py  (mac-mini-2026-06-20-S2)

THM-547 case 2 verification: boundary-collar configs E = E' u {w}, E' subset {0..14} size k-1,
14 < w <= w*(k).  Check p0(E) <= cap_k EXACTLY.  Focus on highest-risk E' (top Plat / top V) x
all w in (14, w*], plus a deterministic stride-sample of the rest.  Report worst margin.
w*: k=8 ->54, k=9 ->90, k=10 ->103.  caps: cap_8=2243/5880, cap_9=1979/4004, cap_10=4/7 (floor).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def plat_V(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); p0=F(0);p1=F(0);arc=[0]*7;inB=[False]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(sector_of(e*xm) for e in E); w=x1-x0
        if len(secs)==7: p0+=w
        miss=[j for j in range(7) if j not in secs]; mj=miss[0] if len(miss)==1 else None
        if mj is not None: p1+=w
        for j in range(7):
            if j==mj:
                if not inB[j]: arc[j]+=1; inB[j]=True
            else: inB[j]=False
    return p0+SEV*p1, sum(arc[1:7])
from math import gcd
def g(t):
    r=0
    for x in t: r=gcd(r,x)
    return r

caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7)}; wstar={8:54,9:90,10:103}
for k in (8,9,10):
    r=k-1; base=list(range(1,15))
    # rank E' by Plat (risk) once; take top-200 plus a stride sample
    scored=[]
    for combo in itertools.combinations(base, r-1):
        Ep=(0,)+combo; pl,V=plat_V(Ep); scored.append((pl,Ep))
    scored.sort(reverse=True)
    topEps=[ep for _,ep in scored[:200]]
    strideEps=[ep for _,ep in scored[200::max(1,len(scored)//300)]]
    testEps=topEps+strideEps
    worst=None; nchk=0; viol=0
    for Ep in testEps:
        for w in range(15, wstar[k]+1):
            E=Ep+(w,)
            if g(E)!=1: continue
            p0=measS7(E); m=caps[k]-p0; nchk+=1
            if m<0: viol+=1
            if worst is None or m<worst[0]: worst=(m,E)
    print(f"k={k}: checked {nchk} collar configs ({len(testEps)} risky E' x w in (14,{wstar[k]}]); "
          f"violations={viol}; worst margin={float(worst[0]):.5f} at E={worst[1]}")
print("\nDONE  (violations=0 in all rows => boundary-collar finite check passes on the risk-ranked sample).")
