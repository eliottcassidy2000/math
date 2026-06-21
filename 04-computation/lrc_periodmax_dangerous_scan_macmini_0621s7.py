#!/usr/bin/env python3
"""
lrc_periodmax_dangerous_scan_macmini_0621s7.py  (mac-mini-2026-06-21-S7)
THREAD 1 (THM-563 general bounded-base finite check), STAGE 1: SCAN.

For ALL primitive bounded bases B = {0} u combo, combo subset of [1,14], |B|=k-1, k=9..13:
  - Plat(B) = decorrelated plateau = measS7(B) + p1(B)/7  (canonical breakpoint method)
  - margin_k = cap_k - Plat(B)   (skip base if margin<=0: that base is over-cap and not a valid plateau)
  - V = #arc-endpoints of the A_j arcs (THM-546)
  - absolute bound period-max <= (6/49)*V   (THM-546 / HYP-2784)
  - TRIVIAL PASS if (6/49)*V <= 15*margin  (no exact period-max needed)
  - else DANGEROUS: report B, Plat, margin, V, absbound, 15*margin, period P=7*lcm(B)

Goal: count dangerous bases per k and their periods, to size the exact computation in stage 2.
caps (canonical, from orchestrator): cap_8..13 = 2243/5880,1979/4004,55/91,66/91,6/7,1.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def lcm(a,b): return a*b//gcd(a,b)
def sector_of(p): return int((p%1)*7)

def setup(E):
    """Return (Plat, arcs, V) where arcs = list of (j,a,b) maximal sector-j-missing intervals,
    V = number of distinct arc endpoints, Plat = measS7 + p1/7 (canonical breakpoint plateau)."""
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
    # V = distinct endpoints (exclude t=0 since S_j(0)-S_j(0) contributes nothing periodic-wise? keep all)
    eps=set()
    for (j,a,bb) in arcs:
        eps.add(a); eps.add(bb)
    V=len(eps)
    return p0+p1/7, arcs, V

caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
ABS=F(6,49)  # per-endpoint absolute bound |S_j|<=3/49, two contributions; THM-546 uses (6/49)*V

for k in range(9,14):
    cap=caps[k]
    n_total=0; n_skip_overcap=0; n_trivial=0; dangerous=[]
    for combo in itertools.combinations(range(1,15),k-2):
        B=(0,)+combo
        plat,arcs,V=setup(B)
        n_total+=1
        margin=cap-plat
        if margin<=0:
            n_skip_overcap+=1
            continue
        absbound=ABS*V
        thr=15*margin
        if absbound<=thr:
            n_trivial+=1
        else:
            L=1
            for e in B:
                if e>0: L=lcm(L,e)
            P=7*L
            dangerous.append((B,plat,margin,V,absbound,thr,P))
    print(f"\n=== k={k}  cap={cap}={float(cap):.5f}  total bases C(14,{k-2})={n_total} ===")
    print(f"  over-cap (Plat>=cap, skipped): {n_skip_overcap}")
    print(f"  TRIVIAL pass via (6/49)V <= 15*margin: {n_trivial}")
    print(f"  DANGEROUS (need exact period-max): {len(dangerous)}")
    if dangerous:
        dangerous.sort(key=lambda r:r[5]-r[4])  # most dangerous: smallest (15m - absbound)
        Pmax=max(r[6] for r in dangerous)
        print(f"    max period P among dangerous = {Pmax}")
        for (B,plat,margin,V,ab,thr,P) in dangerous[:25]:
            print(f"    B={str(B):<34} Plat={float(plat):.5f} margin={float(margin):.5f} V={V} (6/49)V={float(ab):.4f} 15m={float(thr):.4f} P={P}")
        if len(dangerous)>25: print(f"    ... ({len(dangerous)-25} more)")
