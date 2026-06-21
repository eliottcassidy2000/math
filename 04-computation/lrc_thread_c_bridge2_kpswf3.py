#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_thread_c_bridge2_kpswf3.py   (kind-pasteur 2026-06-21, THREAD C lead-gen)

Bridge1 result: p0_sectorcover(C) <= 1-L(C) always (in tested cases), i.e. L(C) <= 1-p0(C).
This is an UPPER bound on L, the WRONG direction for OPEN-Q-108 (which wants a LOWER bound).
=> the naive "p0<=cap => L>=1-cap" transfer is REFUTED.

WHY: sector-cover and danger-cover are different events on tau.
  - danger(v,tau): ||v tau||<=1/14   (tau near a grid multiple a/v) -- the OPEN-Q-108 / lonely event
  - sectorhit(v,tau): floor(7 frac(v tau)) = some s -- which 7-sector v lands in
G_C = {tau: NO v in danger}. The complement is "SOME v in danger".
The sector-cover event is "the v's TOGETHER paint all 7 sectors". These need not nest.

This script:
  (A) Establishes the EXACT set relation. Claim: danger(v,tau) <=> sector(v tau) in {0,6}? NO.
      Danger band ||vt||<=1/14 = [0,1/14]u[13/14,1) (mod 1) of (v tau). Sector0=[0,1/7),
      sector6=[6/7,1). 1/14<1/7 and 13/14>6/7, so DANGER is a STRICT SUBSET of sector0 u sector6.
      meas(danger)=1/7; meas(sector0 u sector6)=2/7. Verify exactly.
  (B) The RIGHT bridge object: define the SECTOR route at the SAME 1/14 resolution.
      Use 14 half-sectors (mod 14, width 1/14) instead of 7. Then danger = exactly half-sector
      {0} u {13}? danger=[0,1/14]u[13/14,1) = halfsector 0 (=[0,1/14)) plus halfsector 13.
      So at resolution 14, "lonely" = "frac(v tau) avoids the danger half-sector".
      L(C) = meas{tau: union over v of (v lands in danger-halfsector) is EMPTY}
           = meas{tau: NO v lands in halfsector-0-region}.  This is an AVOIDANCE, not a cover.
  (C) DUALITY TEST: is L(C) = meas(avoid) related to a COVER on the COMPLEMENT speeds/structure?
      The sector-cover p0 is a COVER (hit all 7); the lonely L is an AVOID (miss the danger band).
      These are genuinely different (cover of 7 vs avoid 1). Quantify the asymptotic decorrelated
      values of BOTH to see if EITHER is the binding constraint.
"""
import itertools
from fractions import Fraction as Fr

P = 7

def main():
    print("="*82)
    print("THREAD C BRIDGE 2: the EXACT set relation danger vs sector, and which route binds")
    print("="*82)

    # (A) danger band vs sectors, exact measures
    print("\n[A] danger band ||theta||<=1/14 vs 7-sector partition:")
    print("    meas(danger) = 2*(1/14) = 1/7   (one band of total width 1/7)")
    print("    danger = [0,1/14] u [13/14,1).  sector0=[0,1/7), sector6=[6/7,1).")
    print("    1/14=%.4f < 1/7=%.4f ; 13/14=%.4f > 6/7=%.4f" % (1/14,1/7,13/14,6/7))
    print("    => DANGER is a STRICT SUBSET of (sector0 u sector6). NOT equal.")
    print("    => 'avoid danger' (lonely) is WEAKER than 'avoid sectors 0,6'.")

    # The decorrelated (q->inf) value of each route, per single speed:
    # P(v in danger) = 1/7 (band measure). P(v hits a given sector) = 1/7.
    print("\n[B] DECORRELATED single-speed probabilities (independent uniform model):")
    print("    P(one speed in danger) = 1/7 = %.4f" % (1/7))
    print("    P(k indep speeds ALL avoid danger) = (6/7)^k  [lonely decorrelated floor]")
    for k in [7,8,12,13]:
        print("       k=%2d: (6/7)^k = %.5f   <- decorrelated L floor (no relations)" % (k,(6/7)**k))
    print("    These are POSITIVE and bounded below for fixed k. The decorrelated lonely")
    print("    measure NEVER vanishes; L=0 requires arithmetic RELATIONS (the tight locus).")

    # (C) KEY: the two routes attack DIFFERENT 'k'. Sector-cover route works on the FULL 14-set
    # E (it asks if E's 13 speeds cover all 7 sectors at the relevant tau, blocking loneliness).
    # The lonely route works on a 12-CORE C and asks meas(G_C)>=c.
    # The logical connection (THM-525): L(S) = meas(G_C) - meas(G_C cap D_w).
    print("\n[C] THE LOGICAL LINK (THM-525, exact identity):")
    print("    For S=C u {w}:  L(S) = meas(G_C) - meas(G_C cap D_w).")
    print("    meas(G_C) is the OPEN-Q-108 floor (route 2).  D_w is w's danger comb (meas 1/7).")
    print("    The sector-cover p0(S) bounds when the FULL union of danger combs covers [0,1):")
    print("    p0 relates to meas(union D_v over ALL of S). Test the chain numerically below.")

    # numeric: for S=C u {w}, show L(S), meas(G_C), and the sector-cover p0 of S
    from fractions import Fraction as Fr
    def danger_intervals(v):
        v=int(v); ivs=[]
        for k in range(0,v+1):
            lo=max(Fr(14*k-1,14*v),Fr(0)); hi=min(Fr(14*k+1,14*v),Fr(1))
            if lo<hi: ivs.append((lo,hi))
        return ivs
    def umeas(intervals):
        if not intervals: return Fr(0)
        iv=sorted(intervals); tot=Fr(0); cl,ch=iv[0]
        for lo,hi in iv[1:]:
            if lo>ch: tot+=ch-cl; cl,ch=lo,hi
            else: ch=max(ch,hi)
        tot+=ch-cl; return tot
    def L_meas(C):
        ai=[];
        for v in C: ai+=danger_intervals(v)
        return Fr(1)-umeas(ai)

    print("\n    S = C u {w}:  L(S) vs meas(G_C) vs L(C):")
    cases = [
        (list(range(1,12))+[13], 36),  # core {1..11,13}, stranger 36 -> the 1/1260 extremizer
        (list(range(1,12))+[13], 24),  # -> sporadic tight L=0
        (list(range(1,12))+[13], 84),  # lcm parking
        ([1,2,3,4,5,6,7,8,9,10,11,12], 14),  # AP core, w=14
    ]
    print(f"      {'core C':<16}{'w':>5}{'meas(G_C)':>12}{'L(S=Cuw)':>12}{'drop':>12}")
    for C, w in cases:
        gC = L_meas(C)
        LS = L_meas(C+[w])
        label = "{1..11,13}" if C[-1]==13 else "{1..12}"
        print(f"      {label:<16}{w:>5}{float(gC):>12.5f}{float(LS):>12.5f}{float(gC-LS):>12.5f}")
    print("\n    => meas(G_C) is the dominant term; the stranger w only shaves off")
    print("       meas(G_C cap D_w) <= meas(D_w)=1/7. So meas(G_C)>=c & c>1/7 => L(S)>0 directly.")
    print("       The whole game IS lower-bounding meas(G_C). p0/sector-cover is the WRONG functional")
    print("       for THIS (it's an upper bound on the danger-cover). BUT the L7 closure's MACHINERY")
    print("       (resonance atlas + torus discrepancy) may transfer. See bridge3.")

if __name__ == "__main__":
    main()
