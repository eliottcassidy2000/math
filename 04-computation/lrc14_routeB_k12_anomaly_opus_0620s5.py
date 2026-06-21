#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_k12_anomaly_opus_0620s5.py   (opus-2026-06-20-S5)

S4 found: at k=12, cut (I) "consec maximizes p0" FAILS at E=[0..10,12].
But the survival certificate only needs the WEIGHTED sum U4=1-G_1+G_5+4G_6 to be
consec-max, not each cut.  Two questions:
  (Q1) Does E=[0..10,12] (and any other) actually BEAT consec on U4 itself at k=12?
       (If U4(consec) is still bank-max, the cut-by-cut proof fails but the
        TARGET HYP-2693 route survives -- it's genuinely joint at k=12.)
  (Q2) Characterize: when cut (I) fails, do (II)+(III) over-compensate so U4 still
       favors consec?  Print the full survival vectors for the offenders.
  (Q3) Re-examine k=8,9,10 (the ONLY tight rows per THM-535) with a MAXIMAL exact bank
       (every primitive shape up to the largest span the exact engine can afford) to
       see if the survival certificate is clean exactly where it MATTERS.
       THM-535: cap_k>=(k-6)/7 and consec<(k-6)/7 for k>=9, so the binding rows are
       k=8,9,10.  At k>=11 the floor route (codex) handles it; consec need NOT be the
       global U4-max there.  So the survival certificate is the RIGHT tool precisely
       on k=8,9,10 where it is CLEAN.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = (e*mid) % 1
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p
def survival(p): return [sum(p[t] for t in range(a,7)) for a in range(7)]
def U4(p): return p[0]+p[5]+5*p[6]
def consec(k): return list(range(k))
def primitive(E):
    g=reduce(gcd,[e for e in E if e>0]); return tuple(e//g for e in E)

print("="*78); print("(Q1/Q2) k=12: does anyone beat consec on U4 even though cut(I) fails?")
print("="*78)
k=12; C=consec(k); pc=dist_p(C); Gc=survival(pc); Uc=U4(pc)
print(f"consec_12: U4={Uc}={float(Uc):.6f}  p0={float(pc[0]):.6f}  G=[{','.join(f'{float(Gc[b]):.4f}' for b in range(1,7))}]")
span=14
seen=set(); beats=[]; cutI_fail=[]
for rest in itertools.combinations(range(1,span+1),k-1):
    E=[0]+list(rest); pr=primitive(E)
    if pr in seen: continue
    seen.add(pr)
    p=dist_p(E); G=survival(p); U=U4(p)
    if U>Uc+F(1,10**15): beats.append((E,U,p,G))
    if G[1]<Gc[1]-F(1,10**15): cutI_fail.append((E,U,p,G))
print(f"bank={len(seen)} primitive shapes (span<= {span})")
print(f"shapes BEATING consec on U4: {len(beats)}")
for E,U,p,G in beats[:10]:
    print(f"   {E}: U4={float(U):.6f} (consec {float(Uc):.6f})  p0={float(p[0]):.5f}")
print(f"\nshapes where cut(I) fails (p0>consec) but U4 check:")
for E,U,p,G in cutI_fail[:6]:
    rel = "BEATS U4" if U>Uc+F(1,10**15) else "still <= consec U4 (joint compensation)"
    print(f"   {E}: p0={float(p[0]):.5f}(>{float(pc[0]):.5f}) U4={float(U):.6f} -> {rel}")
    print(f"       G_b={[float(G[b]) for b in range(1,7)]}  vs consec {[float(Gc[b]) for b in range(1,7)]}")

print("\n"+"="*78)
print("(Q3) MAXIMAL exact bank on the TIGHT rows k=8,9,10 (THM-535 binding rows)")
print("="*78)
for k in (8,9,10):
    C=consec(k); pc=dist_p(C); Gc=survival(pc); Uc=U4(pc)
    span = 16 if k==8 else (15 if k==9 else 14)
    seen=set(); n=0; v1=v5=v6=0; beats=0
    worst1=F(0); t1=None
    for rest in itertools.combinations(range(1,span+1),k-1):
        E=[0]+list(rest); pr=primitive(E)
        if pr in seen: continue
        seen.add(pr); n+=1
        p=dist_p(E); G=survival(p); U=U4(p)
        if U>Uc+F(1,10**15): beats+=1
        if G[1]<Gc[1]-F(1,10**15):
            v1+=1
            if G[1]-Gc[1]<worst1: worst1=G[1]-Gc[1]; t1=E
        if G[5]>Gc[5]+F(1,10**15): v5+=1
        if G[6]>Gc[6]+F(1,10**15): v6+=1
    print(f"k={k}: bank={n} (span<= {span})  U4 beats={beats}  "
          f"cut(I) fail={v1}  cut(II) fail={v5}  cut(III) fail={v6}"
          + (f"  worst(I) by {float(worst1):.2e} at {t1}" if t1 else ""))
    print(f"      => survival certificate {'CLEAN (all 3 cuts + U4)' if (v1==v5==v6==beats==0) else 'NOT clean'}"
          f" on the binding row k={k}")
