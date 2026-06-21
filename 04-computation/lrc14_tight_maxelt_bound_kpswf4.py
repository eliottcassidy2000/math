#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_tight_maxelt_bound_kpswf4.py  (kind-pasteur, 2026-06-21, THREAD D / AP-saturation finiteness)

GOAL: turn the empirical "tight locus = {AP{1..13}, {1..11,13,24}}" into an a-priori
STRUCTURAL finiteness statement by bounding the LARGEST element of any tight 13-set.

DEFS (exact, rational arc sweep; matches lrc14_lonely_helpers):
  danger(v) = { tau in [0,1) : ||v tau|| < 1/14 } = v teeth, each width 1/(7v), total 1/7.
  L(S) = 1 - measure( union_v danger(v) ).  S is TIGHT iff L(S)=0 (dangers cover [0,1)).

CORE STRUCTURAL LEMMA (the finiteness driver), tested here:
  Let V = max(S), and S' = S \ {V} (12 speeds).  The danger set danger(V) is a union of
  V teeth of width 1/(7V) centred at k/V.  Its complement-within-[0,1) GAPS are V arcs,
  each of width 6/(7V), centred at (k+1/2)/V.  For L(S)=0 these gaps must be covered by
  the OTHER 12 dangers.  In particular EVERY gap-CENTER c_k=(2k+1)/(2V) must lie in some
  danger(v), v in S' -- i.e. ||v c_k|| < 1/14 for some v<V.
  => For each k in {0..V-1}: min_{v in S'} || v (2k+1)/(2V) || < 1/14.
  This is EXACTLY the LRC(14) gap condition M(S')< 1/14 sampled at the V gap-centers of V.

  The 12 dangers have TOTAL measure 12/7.  They must cover at least the V gaps of total
  measure 6/7.  Necessary: 12/7 >= 6/7 (always true) -- but the gaps of V are PERIODIC with
  period 1/V, while danger(v) is periodic with period 1/v.  A single danger(v) can cover the
  gap-centers of V over a full period only if v and V resonate.  We test the strongest
  necessary covering inequalities and look for an a-priori V-bound.

OUTPUTS:
  (1) Exact tight check over a wide window + the per-gap-center covering census.
  (2) The "uncoverable gap-center" obstruction: smallest tau-gap of V that NO small v can reach.
  (3) Empirical V-bound evidence: largest V achievable as a tight max, vs a derived bound.
"""
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

# ---------- exact interval / danger / lonely ----------
def danger_arcs(v, h=F(1, 14)):
    A = []
    w = h / v
    for k in range(v + 1):
        c = F(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(F(0), hi), (1 + lo, F(1))]
        elif hi > 1:
            A += [(lo, F(1)), (F(0), hi - 1)]
        else:
            A.append((lo, hi))
    return A

def union_measure(A):
    A = sorted((a, b) for a, b in A if b > a)
    tot = F(0); cl = ch = None
    for a, b in A:
        if ch is None: cl, ch = a, b
        elif a <= ch: ch = max(ch, b)
        else: tot += ch - cl; cl, ch = a, b
    if ch is not None: tot += ch - cl
    return tot

def L_exact(S, h=F(1, 14)):
    A = []
    for v in S: A += danger_arcs(v, h)
    return 1 - union_measure(A)

def primitive(S): return reduce(gcd, S) == 1

def frac_norm(x):
    f = x - (x.numerator // x.denominator)
    return f if f <= F(1, 2) else 1 - f

# ---------- gap-center covering census ----------
def gap_centers(V):
    """Centers of the V complement-gaps of danger(V): (2k+1)/(2V), k=0..V-1."""
    return [F(2 * k + 1, 2 * V) for k in range(V)]

def covered_by(c, Sp, h=F(1, 14)):
    """Is gap-center c in some danger(v), v in Sp?  i.e. ||v c|| < h."""
    return any(frac_norm(v * c) < h for v in Sp)

def uncovered_gap_centers(S, h=F(1, 14)):
    """For the largest V=max(S): which gap-centers of V are NOT covered by the others?"""
    V = max(S); Sp = [v for v in S if v != V]
    bad = [c for c in gap_centers(V) if not covered_by(c, Sp, h)]
    return V, bad

# ---------- main ----------
if __name__ == "__main__":
    print("=" * 78)
    print("LRC(14) TIGHT-LOCUS MAX-ELEMENT BOUND  (AP-saturation finiteness, kpswf4)")
    print("=" * 78)
    AP = tuple(range(1, 14))
    spor = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)

    # sanity
    print(f"\nL(AP) = {L_exact(AP)} (tight={L_exact(AP)==0})")
    print(f"L(spor{spor[-3:]}...) = {L_exact(spor)} (tight={L_exact(spor)==0})")

    # ---- (1) gap-center obstruction for the two known tight configs ----
    print("\n--- (1) gap-center covering of the LARGEST speed in the known tight configs ---")
    for S in (AP, spor):
        V, bad = uncovered_gap_centers(S)
        print(f"  S(max={V}): {len(bad)} uncovered gap-centers of V (must be 0 for tight): {bad[:5]}")

    # ---- (2) the structural max-element bound ----
    # Each gap of V has width 6/(7V). To cover gap-center c=(2k+1)/(2V), need v with
    # ||v c|| < 1/14, i.e. v c within 1/14 of an integer. The smallest v that can hit a
    # GENERIC center needs v/(2V) near-integer-ish. We measure: over all tight configs in a
    # window, what is the max V, and how does it relate to lcm / the 12 helper speeds.
    print("\n--- (2) wide exhaustive census: max tight element vs window ---")
    for top in (18, 20, 24, 28, 34):
        tights = []
        for S in combinations(range(1, top + 1), 13):
            if S[0] != 1:  # primitivity needs 1 or gcd 1; quick: require min=1 OR check gcd
                if not primitive(S):
                    continue
            if L_exact(S) == 0:
                tights.append(S)
        maxV = max((max(S) for S in tights), default=None)
        print(f"  window [1,{top}]: {len(tights)} tight, max element among them = {maxV}")
        for S in tights:
            print(f"      {S}  lcm={reduce(lambda a,b:a*b//gcd(a,b),S)}")

    # ---- (3) far-outlier sweep: can ONE large element ever be tight? ----
    # Replace one AP element by a large w; or keep {1..11,13} and add huge w; test up to big W.
    print("\n--- (3) far-outlier sweep: largest tight max element over big single replacement ---")
    best_far = None
    for j in range(13):
        kept = [AP[i] for i in range(13) if i != j]; ks = set(kept)
        for w in range(14, 4001):
            if w in ks: continue
            S = tuple(sorted(kept + [w]))
            if not primitive(S): continue
            if L_exact(S) == 0:
                if best_far is None or max(S) > max(best_far):
                    best_far = S
                if max(S) >= 25:
                    print(f"      TIGHT far config max={max(S)}: {S}")
    print(f"  largest single-replacement tight max element (w<=4000): "
          f"{max(best_far) if best_far else None}  config={best_far}")

    # ---- (4) the covering deficit: for a candidate large-V config, how big is L? ----
    # This quantifies WHY large V fails: the gaps of V become too fine for 12 fixed helpers.
    print("\n--- (4) covering deficit grows with V (why finiteness holds): ---")
    print("  Take helpers = {1..11,13} (the densest 12-helper set), add w; report L(w).")
    helpers = [1,2,3,4,5,6,7,8,9,10,11,13]
    rows = []
    for w in [24, 25, 26, 30, 36, 48, 72, 100, 200, 500, 1000]:
        if w in helpers: continue
        S = tuple(sorted(helpers + [w]))
        Lv = L_exact(S)
        rows.append((w, Lv))
        print(f"      w={w:5d}: L={Lv}={float(Lv):.6f}  {'TIGHT' if Lv==0 else ''}")
    print("\n  Interpretation: L(w)=0 only at the resonant w (24); generic large w has L>0 and")
    print("  L is bounded BELOW by the asymptotic covering deficit of 12 fixed helpers (V->inf).")
