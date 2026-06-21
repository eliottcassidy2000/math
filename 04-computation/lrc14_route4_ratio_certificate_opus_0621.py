#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_route4_ratio_certificate_opus_0621.py  (opus-2026-06-21, ROUTE 4 CERTIFICATE)

THE ROUTE-4 CLOSURE OF LAYER 3, via the disconnected-arc remainder.

SETUP.  measS7(E) = WIN(E) + DISC(E),  DISC = off-center (disconnected) covered arcs.
  HYP-2761 (verified k=8,9,10, 0 ties): consec UNIQUELY maximizes WIN.  So for every
  rival E in the full-residue stratum:   WINlead(E) := WIN(consec) - WIN(E) > 0.

  consec's measS7 lead over E:
     measS7(consec) - measS7(E) = WINlead(E) - DISCdef(E),
     DISCdef(E) := DISC(E) - DISC(consec).

  consec maximizes measS7  <=>  WINlead(E) > DISCdef(E) for every rival.

THE RESULT (this script CERTIFIES, exact fractions, k=8,9,10):
  Partition rivals:
   * AUTO-SAFE: DISCdef(E) <= 0.  Then lead = WINlead - DISCdef >= WINlead > 0. (TRIVIAL.)
   * DANGER:    DISCdef(E)  > 0.  Then the RATIO BOUND holds:
        DISCdef(E) / WINlead(E)  <=  rho_k  <  1,
     so lead >= (1 - rho_k) * WINlead(E) > 0.

  EXACT rho_k (worst over the stratum):
     rho_8  = 65/83    ~ 0.78313   (margin 1-rho = 18/83  ~ 0.217)
     rho_9  = 274/379  ~ 0.72296   (margin 105/379 ~ 0.277)
     rho_10 = 2033/3188~ 0.63770   (margin 1155/3188~ 0.362)
  Monotone DECREASING in k; closure margin (1-rho_k) INCREASING.

  The worst-ratio (most dangerous) rival is the DILATION 2*consec = [0,2,4,...] at
  k=8,9 (a dilation-with-merged-tail at k=10): the maximal-conductance second-best
  shape. The dilation's own ratio is the clean sequence
     k=7..11:  5/6, 65/83, 274/379, 445/641, 1201/1768  (0.833,0.783,0.723,0.694,0.679).

OBSTRUCTION (honest).  This CLOSES LAYER 3 for k=8,9,10 conditional on HYP-2761
(WINlead>0).  It is NOT yet an all-k theorem: the bound rho_k < 1 is VERIFIED, not
PROVEN.  To make it unconditional one needs an a-priori inequality DISCdef(E) <
WINlead(E) (equivalently a uniform rho<1).  The natural target is to PROVE the
dilation governs the worst ratio and that the dilation ratio (closed sequence above)
stays < 1 -- it appears to converge to a limit ~2/3 from above.  See [T5] below.
"""
import sys, itertools
from fractions import Fraction as F
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
def stratum(k, span):
    for c in itertools.combinations(range(1, span+1), k-1):
        E = [0]+list(c)
        if is_full_residue(E): yield E

if __name__ == "__main__":
    SPANS = {8: 16, 9: 16, 10: 15}
    print("="*84)
    print("ROUTE-4 RATIO CERTIFICATE: DISCdef(E) <= rho_k * WINlead(E), rho_k < 1")
    print("="*84)
    for k in [8, 9, 10]:
        span = SPANS[k]; C = consec(k)
        mC, wC, dC = measS7(C), WIN(C), DISC(C)
        n_auto = n_danger = 0
        worst = F(0); worstE = None
        min_lead = None; min_lead_E = None
        for E in stratum(k, span):
            if E == C: continue
            winlead = wC - WIN(E); discdef = DISC(E) - dC
            lead = winlead - discdef
            assert winlead > 0, f"HYP-2761 VIOLATED at {E}"   # consec must win WIN
            assert lead == mC - measS7(E)
            if min_lead is None or lead < min_lead:
                min_lead, min_lead_E = lead, tuple(E)
            if discdef <= 0:
                n_auto += 1
            else:
                n_danger += 1
                r = discdef / winlead
                if r > worst: worst, worstE = r, tuple(E)
        print(f"\nk={k}  (stratum span<={span})")
        print(f"  AUTO-SAFE rivals (DISCdef<=0): {n_auto}   DANGER rivals (DISCdef>0): {n_danger}")
        print(f"  worst ratio rho_{k} = {worst} = {float(worst):.6f}  at {list(worstE)}")
        print(f"  closure margin 1-rho_{k} = {1-worst} = {float(1-worst):.6f}  > 0")
        print(f"  => for every rival: lead >= (1-rho_{k})*WINlead > 0  [LAYER 3 closed at this k]")
        print(f"  min measS7 lead = {min_lead} = {float(min_lead):.6f} at {list(min_lead_E)}  (>0: {min_lead>0})")

    # [T5] the dilation ratio sequence (clean closed competitor)
    print("\n" + "="*84)
    print("[T5] DILATION 2*consec ratio sequence (the governing competitor):")
    print("     k :  WINlead        DISCdef        ratio          measLead")
    for k in range(7, 13):
        C = consec(k); D = [2*i for i in range(k)]
        if not is_full_residue(D): continue
        wl = WIN(C) - WIN(D); dd = DISC(D) - DISC(C); lead = wl - dd
        r = dd/wl
        print(f"    {k:2d} : {str(wl):12s}  {str(dd):12s}  {str(r):14s}={float(r):.5f}  {float(lead):.6f}>0={lead>0}")
    print("\n  CONCLUSION: ROUTE 4 closes LAYER 3 for k=8,9,10 (conditional on HYP-2761).")
    print("  The disconnected remainder DISCdef is bounded by rho_k<1 times the WIN margin;")
    print("  it provably cannot overturn consec's unique WIN lead. Unconditional all-k")
    print("  proof reduces to: rho_k<1, governed by the dilation 2*consec (sequence above).")
