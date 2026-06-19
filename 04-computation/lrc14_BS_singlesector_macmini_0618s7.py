#!/usr/bin/env python3
"""
lrc14_BS_singlesector_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A baseline

The SIMPLEST rigorous upper bound on meas(S7(E)), as a baseline before the Vaaler refinement:
For ANY fixed sector j in {1..6}:   1_{S7}(x) <= 1 - g_j(x),   g_j= "sector j empty".
=>  meas(S7(E)) <= 1 - P_j,   P_j := meas{x: sector j empty} = int prod_i (1-psi_j(e_i x)) dx.
Take the best (max) P_j:  meas(S7) <= 1 - max_j P_j  =: U1(E).
Also the 2-sector inclusion-exclusion upper bound (Bonferroni even order):
  1_{S7} = prod(1-g_j); but a rigorous *upper* IE truncation of meas(S7) directly:
  meas(S7) <= 1 - sum_j P_j + sum_{j<l} P_{jl}   (even Bonferroni on the union of empties is
  a LOWER bound on the union meas, i.e. UPPER bound on its complement meas(S7) -- WRONG parity).
  The correct Bonferroni: meas(union empties) >= sum P_j - sum P_{jl}  (lower bd on union)
  => meas(S7)=1-meas(union) <= 1 - sum_j P_j + sum_{j<l} P_{jl}.  We compute this too.

All P_j, P_{jl} are EXACT rationals via the breakpoint method (the empty-sector set is a finite
union of rational intervals).  This gives a fully rigorous, exact-rational certificate.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def meas_empty(E, J):
    """meas{x in [0,1): for every e in E, frac(e x) avoids every sector in J}.
       J = list of sector indices. 'sector j empty' for all j in J simultaneously.
       = meas{x: no e in E has frac(e x) in union_{j in J}[j/7,(j+1)/7)}."""
    E=sorted(set(E)); Jset=set(J)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        empty=True
        for e in E:
            if e==0: continue
            y=(e*xm)%1; sec=int(y*7)
            if sec in Jset: empty=False; break
        # also e=0: frac(0)=0 in sector 0; if 0 in J then sector0 occupied -> but J subset {1..6}
        if empty: total+=x1-x0
    return total

def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add(int(y*7))
        if len(secs)==7: total+=x1-x0
    return total

cap_float={8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}
cap8=F(2243,5880)

shapes_by_k={
 8:[("consec{0..7}",list(range(8))),
    ("perf{0,2,3,4,5,6,7,9}",[0,2,3,4,5,6,7,9]),
    ("dissoc 2^i",[0,1,3,7,15,31,63,127]),
    ("Sidon",[0,1,3,7,12,20,30,44]),
    ("spread 0..3,40..43",[0,1,2,3,40,41,42,43])],
 9:[("consec{0..8}",list(range(9)))],
 10:[("consec{0..9}",list(range(10)))],
 11:[("consec{0..10}",list(range(11)))],
}

print("="*96)
print("ANGLE A baseline: rigorous single-sector + 2-Bonferroni upper bounds on meas(S7)")
print("="*96)
print(f"{'shape':<26}{'k':>3}{'measS7':>10}{'U1=1-maxPj':>12}{'U2(Bonf2)':>11}{'cap_k':>9}{'U1<=cap?':>10}{'U2<=cap?':>10}")
print("-"*96)
for k in sorted(shapes_by_k):
    capk = float(cap8) if k==8 else cap_float[k]
    for name,E in shapes_by_k[k]:
        s7=measS7_geom(E)
        Pj=[meas_empty(E,[j]) for j in range(1,7)]
        U1=F(1)-max(Pj)
        # 2-Bonferroni: meas(S7) <= 1 - sum_j P_j + sum_{j<l} P_{jl}
        sP=sum(Pj,F(0))
        sPP=F(0)
        for j,l in itertools.combinations(range(1,7),2):
            sPP+=meas_empty(E,[j,l])
        U2=F(1)-sP+sPP
        f1="OK" if float(U1)<=capk else "VIOLATE"
        f2="OK" if float(U2)<=capk else "VIOLATE"
        print(f"{name:<26}{k:>3}{float(s7):>10.4f}{float(U1):>12.4f}{float(U2):>11.4f}{capk:>9.4f}{f1:>10}{f2:>10}")
print("-"*96)
print("U1 = best single-sector bound (rigorous, exact). U2 = even-Bonferroni (rigorous, exact).")
print("If U2 <= cap_k for consec at k=8..11, the 2-term certificate already closes the dangerous rows.")
