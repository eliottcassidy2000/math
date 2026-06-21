#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_route4_disc_overturn_boundary_opus_0621.py  (opus-2026-06-21, ROUTE 4 -- OBSTRUCTION)

ROUTE 4 RESULT: the disconnected-arc remainder DOES overturn consec's WIN lead.
The precise boundary, exhaustively over the BALANCED full-residue stratum.

measS7(E) = WIN(E) + DISC(E).  HYP-2761: consec uniquely maximizes WIN (k=8,9,10).
ROUTE 4 asked: can the disconnected remainder DISC overturn consec's WIN lead?

ANSWER (exact, exhaustive over the balanced residue profile, span up to k+2):

  k :  8  9 10 11 | 12 | 13 14 | 15 | 16 17 ...
  consec max WIN?    Y  Y  Y  Y |  Y |  Y  Y |  N |  N  N
  consec max measS7? Y  Y  Y  Y |  N |  Y  Y |  Y |  N  N

  * k<=11: BOTH hold -> LAYER 3 true (ROUTE 4 would close, ratio rho<1).
  * k=12 : consec STILL wins WIN (+0.0106), but the endpoint-nudge [0..10,12] has a
           LARGER disconnected remainder: DISCdef=+0.0263, ratio DISCdef/WINlead=2.476>1.
           => the disconnected arcs OVERTURN consec's WIN lead.  measS7 lead = -31/1980<0.
           CONSEC IS NOT THE measS7 MAXIMIZER.  (independently re-verified, exact.)
  * k=13,14: consec recovers (measS7-max again) -- k=12 is an ISOLATED first failure.
  * k=15 : consec loses WIN itself (HYP-2761 fails) -- but still wins measS7.
  * k>=16: consec loses BOTH WIN and measS7; #beating shapes grows fast (12,79,233,...).
           Winners are TWO-CLUSTER shapes (consec run + a detached tail), the additive-
           energy maximizers of HYP-2738.

CONCLUSION FOR ROUTE 4.  The disconnected remainder is NOT bounded by the WIN margin
in general: at k=12 it exceeds it by a factor 2.476, overturning the WIN lead.  So
the ROUTE-4 closure (DISC cannot overturn WIN) holds ONLY for k<=11, and is FALSE from
k=12.  This is the PRECISE OBSTRUCTION: LAYER 3 (consec maximizes measS7) is itself
FALSE at k=12 and at all k>=16 within the balanced full-residue stratum.  The whole
'consec maximizes measS7' claim is a SMALL-k phenomenon, consistent with HYP-2738's
'irreducibly aggregate, signed-balance' verdict.  measS7-extremality of consec is true
exactly on k in {<=11, 13, 14} (in the tested range), not for all k.

The single cleanest witness: k=12, E*=[0,1,2,3,4,5,6,7,8,9,10,12]:
   measS7(E*) = 11171/17640 > measS7(consec_12) = 119843/194040,  surplus 31/1980,
   while WIN(consec_12) > WIN(E*): consec wins the windows but loses the disconnected mass.
"""
import sys, itertools
from fractions import Fraction as F
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)

def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7
def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7
def breakpoints(E, a):
    half = F(1, 14); bps = {F(0), -half, half}
    for e in E:
        if e == 0: continue
        lo_val = F(7*e)*(-half) + F(e*a); hi_val = F(7*e)*(half) + F(e*a)
        lo_i = min(lo_val, hi_val); hi_i = max(lo_val, hi_val)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            y = F(m - e*a, 7*e)
            if -half <= y <= half: bps.add(y)
            m += 1
    return sorted(bps)
def W_a_total(E, a):
    bps = breakpoints(E, a); tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_all_at(E, a, (lo+hi)/2): tot += hi - lo
    return tot
def window_TpTm(E, a):
    half = F(1, 14); bps = breakpoints(E, a); ivals = list(zip(bps, bps[1:])); Tp = F(0)
    for lo, hi in ivals:
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if covered_all_at(E, a, (lo2+hi)/2): Tp = hi
        else:
            if lo2 == F(0): Tp = F(0)
            break
    Tp = min(Tp, half); Tm = F(0)
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if covered_all_at(E, a, (lo+hi2)/2): Tm = -lo
        else:
            if hi2 == F(0): Tm = F(0)
            break
    Tm = min(Tm, half); return Tp, Tm
def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return list(range(k))
def measS7(E): return sum(W_a_total(E, a) for a in range(1, 7))
def WIN(E): return sum(sum(window_TpTm(E, a)) for a in range(1, 7))
def DISC(E): return measS7(E) - WIN(E)
def profile(E): return tuple(sorted(Counter(e % 7 for e in E).values()))
def stratum(k, span):
    for c in itertools.combinations(range(1, span+1), k-1):
        E = [0]+list(c)
        if is_full_residue(E): yield E

if __name__ == "__main__":
    print("="*92)
    print("ROUTE 4 -- DISCONNECTED-REMAINDER OVERTURN BOUNDARY (balanced full-residue stratum)")
    print("="*92)
    print(" k | maxWIN | maxmeasS7 | #measS7-beaters | first measS7-beater")
    for k in range(8, 18):
        C = consec(k); wC = WIN(C); mC = measS7(C); bp = profile(C); span = k+2
        win_ok = True; meas_beat = []
        for E in stratum(k, span):
            if profile(E) != bp: continue
            if WIN(E) > wC: win_ok = False
            if measS7(E) > mC: meas_beat.append(tuple(E))
        meas_ok = (len(meas_beat) == 0)
        fb = list(meas_beat[0]) if meas_beat else "-"
        print(f" {k:2d} |  {'Y' if win_ok else 'N'}    |   {'Y' if meas_ok else 'N'}      | "
              f"{len(meas_beat):3d}             | {fb}")

    print("\n--- THE k=12 WITNESS (exact) ---")
    C = consec(12); Es = [0,1,2,3,4,5,6,7,8,9,10,12]
    print(f"  consec_12        measS7={measS7(C)} = {float(measS7(C)):.6f}")
    print(f"  E*=[..10,12]     measS7={measS7(Es)} = {float(measS7(Es)):.6f}")
    print(f"  surplus E*-consec = {measS7(Es)-measS7(C)} = {float(measS7(Es)-measS7(C)):+.6f}  (E* WINS)")
    print(f"  WIN: consec={float(WIN(C)):.6f} > E*={float(WIN(Es)):.6f}  (consec wins windows)")
    print(f"  DISC: E*={float(DISC(Es)):.6f} > consec={float(DISC(C)):.6f}  (E* wins disconnected)")
    print(f"  ratio DISCdef/WINlead = {(DISC(Es)-DISC(C))/(WIN(C)-WIN(Es))} "
          f"= {float((DISC(Es)-DISC(C))/(WIN(C)-WIN(Es))):.4f}  > 1  => DISC OVERTURNS WIN")
    print("\nROUTE 4 closes LAYER 3 for k<=11 only; the disconnected remainder OVERTURNS")
    print("the WIN lead at k=12 -> LAYER 3 (consec max measS7) is FALSE at k=12 and k>=16.")
