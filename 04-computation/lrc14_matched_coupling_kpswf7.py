#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- the MATCHED (permutation) coupling + t-profile of weak-maj failures.

Findings so far:
 - SORTED weak-majorization of the W-cell vector implies the wall (t=6 entry) and
   holds for ALL rivals at k=8,11 but FAILS at exactly ONE rival each at k=9,10,
   ALWAYS at t=1 (a single rival cell beats consec's single best cell).
   => sorted-majorization is STRICTLY STRONGER than the wall; it's not the right
      certificate (the wall is irreducibly about the SUM).

This script tests the genuinely JOINT coupling candidates:

 (M1) t-PROFILE: for EVERY rival, at which t does the partial-sum gap
      D_t = sum_{i<=t}(WC_i - WE_i) first go negative? Is it ALWAYS t=1? If the
      only failures are t=1 and D_t recovers (D_t>=0) for all t>=2 AND at t=6
      (the wall), then we have: "the deficit is confined to the single top cell;
      cells 2..6 of consec already dominate the partial sums" -- a 2-part bound:
      (top cell handled by HYP-2761 WIN extremality?) + (rest by majorization).

 (M2) MATCHED coupling: instead of sorting, MATCH cell a of consec to cell sigma(a)
      of E by a FIXED bijection sigma (e.g. identity, or the reflection a->7-a, or
      the dilation a->c*a mod 7). Does there exist a single sigma (per rival, or
      universal) with W_a(consec) >= W_{sigma(a)}(E) for all a? Universal sigma =
      a clean coupling. Test identity, all 6 reflections/dilations, and best-match.

 (M3) The HALL/RADO transportation test: is there a doubly-stochastic-style
      coupling? Equivalent to sorted-majorization for the SUM, so we instead test
      the WIN-vector (binding-speed) matched coupling which mac-mini found sharper.
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def W_a(E, a):
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7 * abs(e)
        j0 = math.floor(lo * d); j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l + h) / 2; hit = set()
        for e in E:
            v = e * xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += h - l
    return tot

def Wraw(E): return [W_a(E,a) for a in range(1,7)]  # indexed a=1..6
def measS7(E): return sum(Wraw(E))
def full_residue(E): return set(e % 7 for e in E) == set(range(7))
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return tuple(range(k))

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*88)
    print("MATCHED COUPLING + t-PROFILE of W-cell partial-sum deficits")
    print("="*88)
    for k, box in [(8,16),(9,15),(10,14),(11,13)]:
        print(f"\n{'#'*84}\n# k={k}, box=[0,{box}]\n{'#'*84}")
        C=consec(k); WCraw=Wraw(C); mC=measS7(C)
        WCs=sorted(WCraw,reverse=True)
        S=stratum(k,box)
        # ---- M1: t-profile of sorted weak-maj deficits ----
        first_neg_t = {1:0,2:0,3:0,4:0,5:0,6:0}
        wall_ok=True
        t1_fail_recovers=0; t1_fail_total=0
        for E in S:
            if E==C: continue
            WEs=sorted(Wraw(E),reverse=True); mE=measS7(E)
            if mE>mC: wall_ok=False
            pu=pv=F(0); fneg=None
            negs=[]
            for t in range(6):
                pu+=WCs[t]; pv+=WEs[t]
                if pu<pv:
                    negs.append(t+1)
                    if fneg is None: fneg=t+1
            if fneg is not None:
                first_neg_t[fneg]+=1
                if fneg==1:
                    t1_fail_total+=1
                    # recovers if D_t>=0 for all t>=2 (and t=6 = wall)
                    pu=pv=F(0); rec=True
                    for t in range(6):
                        pu+=WCs[t]; pv+=WEs[t]
                        if t>=1 and pu<pv: rec=False
                    if rec: t1_fail_recovers+=1
        print(f"  wall holds (no measS7 rival > consec)? {wall_ok}")
        print(f"  weak-maj first-negative-t histogram: {first_neg_t}")
        print(f"  among t=1 failures ({t1_fail_total}): partial-sums recover for t>=2? {t1_fail_recovers}/{t1_fail_total}")
        # Interpretation: if all failures are t=1 AND all recover at t>=2,
        # then: consec's TOP-5 cells already dominate (sum_{top5} consec >= sum_{top5} E),
        # and consec's full sum dominates (wall). Only the single largest cell can be beaten.

        # ---- M2: matched coupling under dilation a->c*a mod 7 ----
        # For each rival, is there c in Z_7* with W_a(C) >= W_{(c*a)%7}(E) for all a=1..6?
        cstar=[1,2,3,4,5,6]
        univ_c_works={c:0 for c in cstar}
        any_c=0
        for E in S:
            if E==C: continue
            WEraw=Wraw(E)  # a=1..6 at index a-1
            found=False
            for c in cstar:
                ok=True
                for a in range(1,7):
                    b=(c*a)%7
                    if b==0: ok=False; break
                    if WCraw[a-1] < WEraw[b-1]: ok=False; break
                if ok:
                    univ_c_works[c]+=1; found=True
            if found: any_c+=1
        print(f"  matched dilation coupling: EXISTS some c in Z_7* with W_a(C)>=W_(ca)(E) for all a: {any_c}/{len(S)-1}")
        print(f"     per-c success counts: {univ_c_works}")
