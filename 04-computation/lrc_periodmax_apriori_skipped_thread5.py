#!/usr/bin/env python3
"""THREAD 5: a-priori period-max bound for ALL bounded bases, NO period scan.
THM-563 identity: w*Delta_w = sum over one-miss arcs A_j=[a,b] of [Sc_j(w*b)-Sc_j(w*a)].
Each centered sawtooth Sc_j has range R_j = max(Sc_j)-min(Sc_j). The signed difference
[Sc_j(w*b)-Sc_j(w*a)] is bounded in absolute value by R_j. So
    period-max(B) <= sum_arcs R_j   (a-priori, holds for ALL w, no scan).
If sum_arcs R_j < 15*margin(B) for EVERY bounded base, single-far is CLOSED with no scan
and no skipping. This is the honest universal certificate the broad-scan was approximating.

We compute R_j exactly for each centered sawtooth, count one-miss arcs per sector, and
test sum_arcs R_j < 15*margin over ALL C(14,k-2) bases for k=8..12. Report worst ratio.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def breakpoints(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    return sorted(b)
# centered sawtooth Sc_j = Sj_raw - mean; need its exact range R_j over [0,1).
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
# Sj_raw is piecewise linear with breakpoints {0, a=j/7, b=(j+1)/7, 1}. Extremes at breakpoints.
def range_Sc(j):
    a=F(j,7);b=F(j+1,7);pts=[F(0),a,b,F(1)]
    vals=[Sj_raw(x,j)-MEAN[j] for x in pts]
    return max(vals)-min(vals)
R={j:range_Sc(j) for j in range(1,7)}
print("centered-sawtooth ranges R_j:", {j:(R[j], float(R[j])) for j in range(1,7)})
def arcs_per_sector(E):
    """count one-miss arcs by sector j (each arc contributes one b and one a endpoint diff,
    bounded by R_j)."""
    b=breakpoints(E); cur={}; cnt={j:0 for j in range(1,7)}
    prevmiss=None
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); miss=[j for j in range(1,7) if j not in secs]
        mj=miss[0] if len(miss)==1 else None
        for j in range(1,7):
            active=(mj==j)
            if active and j not in cur: cur[j]=True
            if (not active) and j in cur:
                cnt[j]+=1; del cur[j]
    for j in list(cur): cnt[j]+=1
    return cnt
def plat(E):
    b=breakpoints(E); p0=F(0);p1=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); miss=[j for j in range(1,7) if j not in secs]
        if len(secs)==7: p0+=x1-x0
        if len(miss)==1: p1+=x1-x0
    return p0+p1/7
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(66,91),12:F(6,7)}
for k in (8,9,10,11,12):
    cap=caps[k]; worst=F(0); wB=None; fails=0; nb=0; worst_apriori_pm=None
    for combo in itertools.combinations(range(1,15),k-2):
        B=(0,)+combo; nb+=1
        margin=cap-plat(B)
        cnt=arcs_per_sector(B)
        apriori_pm=sum(cnt[j]*R[j] for j in range(1,7))  # crude upper bound on period-max
        if margin<=0:
            fails+=1; continue
        ratio=apriori_pm/margin
        if apriori_pm>=15*margin: fails+=1
        if ratio>worst:
            worst=ratio; wB=B; worst_apriori_pm=apriori_pm
    print(f"k={k}: {nb} bases, a-priori (sum R_j over arcs) fails(>=15*margin)={fails}, "
          f"WORST ratio={float(worst):.3f} at B={wB} (apriori_pm={float(worst_apriori_pm):.3f})")
