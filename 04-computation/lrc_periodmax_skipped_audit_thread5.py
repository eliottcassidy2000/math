#!/usr/bin/env python3
"""THREAD 5 CRITIC: Audit the 267 SKIPPED (huge-period P>30000) bounded bases.
The general-script CLAIM is they are 'low-Plat, easy'. Test that: for k=8, over ALL
C(14,6)=3003 bases (0 in B, span<=14), compute Plat(B) and margin, and check whether
the period>30000 (skipped) bases are actually all low-margin/easy, OR whether any skipped
base has small margin (=dangerous) where period-max could threaten cap.

We do NOT compute period-max for the huge ones (too slow); instead we verify the
PREMISE that high-Plat <=> small period. If a skipped base has small margin, the
'easy' claim is FALSE and the broad check has a hole.
"""
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
def plat(E):
    b=breakpoints(E); p0=F(0);p1=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        secs=set(sector_of(e*((x0+x1)/2)) for e in E); miss=[j for j in range(1,7) if j not in secs]
        if len(secs)==7: p0+=x1-x0
        if len(miss)==1: p1+=x1-x0
    return p0+p1/7
def period_of(B):
    L=1
    for e in B:
        if e>0: L=lcm(L,e)
    return 7*L
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(66,91),12:F(6,7)}
for k in (8,9,10,11,12):
    cap=caps[k]
    skipped=[]; small=[]
    nbases=0
    for combo in itertools.combinations(range(1,15),k-2):
        B=(0,)+combo; nbases+=1
        P=period_of(B)
        pl=plat(B); margin=cap-pl
        rec=(margin,B,P,pl)
        if P>30000: skipped.append(rec)
        else: small.append(rec)
    skipped.sort()  # smallest margin first
    small.sort()
    # Worst-case period-max bound for a base: |period-max| <= sum over arcs of 2*max|Sc|.
    # max|Sc| over one sawtooth Sj <= ~6/49 + |mean|; but the THM-563 claim is |S_j|<=3/49.
    # A crude a-priori upper bound on period-max = (#one-miss arcs)*2*(3/49). We report margin only.
    print(f"\n=== k={k}: {nbases} bases, {len(skipped)} skipped(P>30000), {len(small)} small-period ===")
    print(f"  cap={float(cap):.5f}")
    print(f"  SMALLEST-margin SKIPPED bases (the dangerous-if-any among huge-period):")
    for margin,B,P,pl in skipped[:8]:
        print(f"    B={str(B):<28} Plat={float(pl):.5f} margin={float(margin):.5f} P={P}")
    print(f"  SMALLEST-margin SMALL-period bases (these WERE checked if in top-20):")
    for margin,B,P,pl in small[:5]:
        print(f"    B={str(B):<28} Plat={float(pl):.5f} margin={float(margin):.5f} P={P}")
    # Is the overall smallest-margin base skipped?
    allrec=sorted(skipped+small)
    worst_margin,wB,wP,wpl = allrec[0]
    print(f"  GLOBAL smallest-margin base: B={wB} margin={float(worst_margin):.6f} P={wP} {'(SKIPPED!)' if wP>30000 else '(checked-region)'}")
