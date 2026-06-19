#!/usr/bin/env python3
"""
ANGLE H part 6 — the REAL ceiling of meas(S7) and whether cap_k can be breached.

meas(S7(E)) only depends on the phase trajectories t -> frac(e_i t).  For
INTEGER E these are coupled (all multiples of t).  The seven-sector reduction
needs meas(S7(E)) <= cap_k for every integer E.  Question: what is the
SUPREMUM of meas(S7) over ALL integer E with |E|=k, 0 in E?  And does it
approach/exceed cap_k?

Key independent sanity facts to confirm:
  - scale invariance: meas_S7(d*E) = meas_S7(E)  (so only primitive shapes matter)
  - the M7(k) "independent main term" is the value you'd get if the k phase
    trajectories were INDEPENDENT uniform — a heuristic, not a bound.
  - consecutive maximizes among small boxes (verified) but at k>=12 AP-beaters
    exist (HYP-2604).  Re-confirm and quantify the beater margin vs cap.

This script: (1) verify scale-invariance exactly; (2) for k=8..13 find the
true integer maximizer in an extended box + structured beaters; (3) report
sup_S7 vs cap_k and the margin.  A breach (sup_S7 > cap_k) at any k would be
a HOLE in HYP-2603/THM-532.
"""
from fractions import Fraction as F
from itertools import combinations
import sys
sys.path.insert(0, "04-computation")
from lrc14_adversarial_chain_macmini_0618sH import meas_S7, M7

cap = {8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}

print("[scale-invariance exact check]")
for E in [[0,1,2,3,4,5,6,7],[0,2,3,4,5,6,7,8]]:
    for d in [1,2,3,5]:
        a = meas_S7(E); b = meas_S7([d*x for x in E])
        print(f"  E={E} d={d}: meas_S7(E)={a}  meas_S7(dE)={b}  {'OK' if a==b else 'MISMATCH!!'}")

print("\n[true integer maximizer of meas_S7, extended search + beaters]")
for k in [8,9,10,11,12,13]:
    capk = cap[k]
    consec = meas_S7(list(range(k)))
    best = consec; argbest = tuple(range(k))
    # extended box: pick B so search is feasible; prioritize structured beaters
    # known HYP-2604 beaters at k=12,13 are near-AP with one element bumped
    # systematically try: consec with last element shifted up by s, and
    # consec with one interior gap, within a moderate box
    B = {8:16,9:15,10:14,11:14,12:14,13:15}[k]
    cnt=0; exceed=0
    for rest in combinations(range(1, B+1), k-1):
        E=(0,)+rest
        v=meas_S7(list(E)); cnt+=1
        if v>best: best=v; argbest=E
        if v>capk: exceed+=1
    margin = capk - best
    flag = "  <<< BREACH! S7>cap" if best>capk else ""
    beat = " (AP NOT max)" if best>consec else " (AP=max)"
    print(f"  k={k}: cap={float(capk):.4f} consec={float(consec):.4f} "
          f"SUP_S7={float(best):.4f} margin={float(margin):.4f}{beat} "
          f"box{B} checked{cnt} #>cap={exceed}{flag}")
    if best>consec:
        print(f"        beater E={argbest}")
