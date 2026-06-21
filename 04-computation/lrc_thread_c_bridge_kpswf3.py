#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_thread_c_bridge_kpswf3.py   (kind-pasteur 2026-06-21, THREAD C lead-gen)

GOAL: pin down the EXACT relationship between the two LRC(14) routes, then test
3-4 concrete hypotheses connecting them.

  ROUTE 1 (sector-cover / L7): p0(E) = meas{ tau : {sector(v tau): v in E} = Z/7 }
     where sector(y)=floor(7 frac(y)). Cap inequality p0<=cap_k => survival.
  ROUTE 2 (lonely measure / OPEN-Q-108): meas(G_C) = L(C) = meas{tau: ||v tau||>1/14 all v in C}
     = int_0^1 prod_{v in C} 1_safe(v tau) dtau,  1_safe(theta)=1_{||theta||>1/14}.

KEY BRIDGE FACT to verify: the danger band ||theta||<=1/14 = exactly the points whose
7-sector label is "0" (the band [-1/14,1/14] = [0,1/14] u [13/14,1) = sector 0 split? NO:
sector(y)=floor(7y); sector 0 = [0,1/7). Danger = [0,1/14]u[13/14,1). These DIFFER.
So the "sector" partition (mod 7) and the "danger band" (mod 1, width 1/7) are RELATED but
NOT identical. We test the precise link below.

Tests:
  T1: Is meas(G_C) = 1 - meas{tau: danger combs of C cover [0,1)}?  (trivially: G_C = complement
      of union of danger combs; verify it's literally 1 - meas(union D_v).)
  T2: BRIDGE: does sector-cover p0 (route 1) lower-bound 1-meas(G_C), i.e. is there an inequality
      meas(G_C) >= 1 - p0_sectorcover(C) or similar? Compute both exactly on small C.
  T3: TIGHT-LOCUS <-> RESONANCE-ATLAS: is the set of C with L(C)=0 (tight) related to the finite
      resonance atlas (small-denom ratios)? Compute tight configs and their ratio structure.
"""
import itertools
from fractions import Fraction as Fr

P = 7
HALF = Fr(1,14)  # gap threshold

def sector(yfrac):
    return int(P*yfrac)

# ---------- exact lonely measure L(C) = meas(G_C) via danger-comb union ----------
def danger_intervals(v):
    """v's danger comb: { tau : ||v tau|| <= 1/14 } = union_k [(14k-1)/(14v),(14k+1)/(14v)] cap [0,1)."""
    v = int(v)
    ivs = []
    # ||v tau||<=1/14  <=> exists integer k with |v tau - k| <= 1/14 <=> tau in [ (k-1/14)/v, (k+1/14)/v ]
    for k in range(0, v+1):
        lo = Fr(14*k-1, 14*v); hi = Fr(14*k+1, 14*v)
        lo = max(lo, Fr(0)); hi = min(hi, Fr(1))
        if lo < hi:
            ivs.append((lo, hi))
    return ivs

def union_measure(intervals):
    if not intervals: return Fr(0)
    iv = sorted(intervals)
    tot = Fr(0); cl, ch = iv[0]
    for lo, hi in iv[1:]:
        if lo > ch:
            tot += ch-cl; cl, ch = lo, hi
        else:
            ch = max(ch, hi)
    tot += ch-cl
    return tot

def L_measure(C):
    """meas(G_C) = 1 - meas(union of danger combs). EXACT."""
    allint = []
    for v in C:
        allint += danger_intervals(v)
    return Fr(1) - union_measure(allint)

# ---------- exact sector-cover p0(C): meas{tau: sectors of C cover Z/7} ----------
def p0_sectorcover(C):
    """meas{ tau in [0,1): {floor(7 frac(v tau)): v in C} = {0,...,6} }. EXACT via breakpoints."""
    C = [int(v) for v in C]
    bp = {Fr(0), Fr(1)}
    for v in C:
        for t in range(0, P*v+1):
            bp.add(Fr(t, P*v))
    xs = sorted(b for b in bp if Fr(0)<=b<=Fr(1))
    tot = Fr(0)
    for a, b in zip(xs, xs[1:]):
        if b<=a: continue
        mid = (a+b)/2
        S = set(sector((v*mid) % 1) for v in C)
        if len(S) == P:
            tot += b-a
    return tot

def main():
    print("="*82)
    print("THREAD C BRIDGE: lonely-measure (route 2) vs sector-cover (route 1), EXACT")
    print("="*82)

    # ---- T1 sanity: G_C is literally complement of danger-comb union ----
    print("\n[T1] L(C)=meas(G_C)=1-meas(union danger combs). Spot checks:")
    for C in [[1,2,3], [1,2,3,4,5], [1,3,4,5,7], list(range(1,14))]:
        L = L_measure(C)
        print(f"   C={C if len(C)<6 else 'AP{1..13}'}: L(C)={L}  ={float(L):.5f}")

    # ---- T2: BRIDGE between sector-cover p0 and lonely measure L ----
    print("\n[T2] BRIDGE: compare p0_sectorcover(C) and 1-L(C)=meas(union danger):")
    print(f"   {'C':<26}{'p0_sectcov':>12}{'1-L(C)':>12}{'L(C)':>10}{'p0 vs 1-L':>12}")
    test_sets = [
        [1,2,3], [1,2,3,4], [1,2,3,4,5], [2,3,4,5,6,7],
        [1,3,4,5,7], [1,2,3,5,8], list(range(1,8)), list(range(1,14)),
        [1,2,3,4,5,6,7,8,9,10,11,13,36],   # the 1/1260 extremizer (drop 12, add 36)
        [1,2,3,4,5,6,7,8,9,10,11,13,24],   # the sporadic tight L=0
    ]
    for C in test_sets:
        L = L_measure(C)
        p0 = p0_sectorcover(C)
        oneminusL = Fr(1)-L
        rel = "p0<=1-L" if p0<=oneminusL else "p0>1-L !!"
        label = ("AP1-13" if C==list(range(1,14)) else
                 "1-11,13,36" if C[-1]==36 else
                 "1-11,13,24" if C[-1]==24 else str(C))
        print(f"   {label:<26}{float(p0):>12.5f}{float(oneminusL):>12.5f}{float(L):>10.5f}{rel:>12}")

    print("\n   Interpretation: if p0(C)<=1-L(C) ALWAYS, then L(C)>=1-p0(C):")
    print("   a sector-cover UPPER bound on p0 would give a lonely-measure LOWER bound (the transfer!).")

if __name__ == "__main__":
    main()
