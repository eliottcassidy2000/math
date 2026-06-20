#!/usr/bin/env python3
"""
lrc14_boundary_collar_cutoff_macmini_0620s2.py  (mac-mini-2026-06-20-S2)

BOUNDARY-COLLAR CLOSURE (codex HYP-2675's collar sub-case) via THM-546 sharpened.
Boundary collar: E = E' u {w}, E' subset [0,14] (second-largest <= 14), w=max(E) > 14.
Since E' is BOUNDED, both Plat(E')=p0(E')+(1/7)p1(E') and V(E')=arc-complexity are bounded by
FINITE maxima over E' subset {0..14}, |E'|=k-1, 0 in E'.  Then:
   p0(E) = Plat(E') + Delta_w  <=  Qb(k-1) + (6/49)V_max/w   [THM-546 sharpened].
For w > w* := (6/49) V_max / (cap_k - Qb(k-1)),  p0(E) < cap_k  => COLLAR CLOSED (gapped).
For 14 < w <= w*: a FINITE check (E' bounded, w bounded).  This script computes Qb(k-1),
V_max, and w* for k=8,9,10, establishing the collar reduces to a finite computation.
caps: cap_8=2243/5880, cap_9=1979/4004, cap_10 from THM-535 (>= (k-6)/7).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)
def sector_of(p): return int((p%1)*7)
def stats(E):
    """return (p0, p1, V=sum_j #arcs(B_j))."""
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); p0=F(0); p1=F(0); arccount=[0]*7; inBj=[False]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(sector_of(e*xm) for e in E); w=x1-x0
        nmiss=7-len(secs)
        if nmiss==0: p0+=w
        missed=[j for j in range(7) if j not in secs]; mj=missed[0] if nmiss==1 else None
        if mj is not None: p1+=w
        for j in range(7):
            if j==mj:
                if not inBj[j]: arccount[j]+=1; inBj[j]=True
            else: inBj[j]=False
    return p0,p1,sum(arccount[1:7])

caps={8:F(2243,5880), 9:F(1979,4004), 10:F(4,7)}  # cap_10 >= 4/7 (THM-535 lower bd; use 4/7 as floor)
print("Boundary-collar cutoff w* (E' subset {0..14}, |E'|=k-1, 0 in E'):")
print(f"{'k':>3}{'#E-prime':>10}{'Qb(k-1)':>12}{'cap_k':>10}{'margin':>10}{'V_max':>7}{'w*':>9}")
for k in (8,9,10):
    r=k-1  # |E'|
    base=list(range(1,15))  # {1..14}; E' = {0} u (r-1 chosen from base)
    bestPlat=F(0); Vmax=0; cnt=0
    for combo in itertools.combinations(base, r-1):
        Ep=(0,)+combo; p0,p1,V=stats(Ep); plat=p0+SEV*p1
        if plat>bestPlat: bestPlat=plat
        if V>Vmax: Vmax=V
        cnt+=1
    margin=caps[k]-bestPlat
    wstar = (F(6,49)*Vmax/margin) if margin>0 else None
    ws = float(wstar) if wstar else float('inf')
    print(f"{k:>3}{cnt:>10}{float(bestPlat):>12.5f}{float(caps[k]):>10.5f}{float(margin):>10.5f}{Vmax:>7}{ws:>9.0f}")
print()
print("If margin>0 and w* is finite: boundary collar = {finite check 14<w<=w*} + {THM-546 gapped w>w*}.")
print("Qb(k-1) = max Plat over bounded E' (the finite-half extremal, kps-S19 style); margin=cap_k-Qb.")
print("NOTE: w* here uses the closed-form (6/49)V_max bound; the full signed sum (5-76x tighter)")
print("shrinks w* by that factor, making the finite intermediate check far smaller.")
print("\nDONE.")
