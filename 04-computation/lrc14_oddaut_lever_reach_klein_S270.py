#!/usr/bin/env python3
"""
lrc14_oddaut_lever_reach_klein_S270.py
======================================
klein-2026-07-12-S270  (owner: explore the odd-|Aut| tie-resolution lever as a proof route)

THE LEVER (S269): at t=1/14 the runners+observer are the 14-grid, which carries the order-2
antipodal symmetry x↦x+1/2; a tournament has odd |Aut| so it can't carry an order-2
automorphism ⟹ the 7 tied antipodal arcs resolve ⟹ "you can't beat 1/14 at base 14".
mac-mini cont.56 twin: the AP is pinned at 1/14 because it CONTAINS its own distance-1
landers {1,13} mod 14.

QUESTION: is this a LOWER-BOUND engine (proves M≥1/14 for all S = LRC) or only an
EQUALITY/UNIQUENESS tool (why the AP achieves 1/14 exactly)?

This script tests the lever's REACH:
 (T1) distance-1-lander: AP@14 landers {1,13} present→1/14; deep well@183 landers {13,170}
      absent→14/183 (reconfirm mac-mini cont.56).
 (T2) DOES "base 14 caps at 1/14" hold for ALL families?  NO: a family avoiding ±1 mod 14
      beats 1/14 at base 14. So the lever pins ONLY the full-grid (AP) family — an EQUALITY
      statement, not a universal lower bound. And covering families give 0 at base 14 (blocked).
 (T3) antipodal parity (klein-S55): on an odd base q the good-multiplier set is closed under
      a↦q−a (the involution ι), so its size is EVEN — a constraint (shadow), not a detector.
 (T4) the LOWER bound is elsewhere: M≥1/14 needs SOME base ≥1/14; for covering families that
      base is large (up to 183), the lander-DODGE crux — which the base-14 lever does NOT give.
"""
import math
from fractions import Fraction

def dist0(r,q): return min(r,q-r)
def best_at_base(v,q):
    """best_q = max_a min_i dist(v_i a mod q,0)/q, and the argmax a."""
    best=-1; arga=None
    for a in range(1,q):
        m=min(dist0((vi*a)%q,q) for vi in v)
        if m>best: best=m; arga=a
    return Fraction(best,q), arga
def exact_M(v):
    best=Fraction(0); argq=None; arga=None
    qs=set()
    n=len(v)
    for i in range(n):
        for j in range(i,n): qs.add(v[i]+v[j])
    for q in sorted(qs):
        b,a=best_at_base(v,q)
        if b>best: best=b; argq=q; arga=a
    return best,argq,arga
def covering(v): return all(any(x%d==0 for x in v) for d in range(2,15))

AP=list(range(1,14))
deep=list(range(1,13))+[182]
# a family avoiding residues {0,1,13} mod 14 (so it BEATS 1/14 at base 14) -- non-covering
avoid=[2,3,4,5,6,8,9,10,11,12,16,17,18]   # all in {2..12} mod 14
covfam=[2,3,4,6,7,8,9,12,13,14,21,22,26]  # has a mult of 14 -> covering-ish

print("="*74)
print("(T1) distance-1-lander mechanism (mac-mini cont.56)")
print("="*74)
for nm,v,q,a in [("AP {1..13}",AP,14,1),("deep well {1..12,182}",deep,183,14)]:
    res=sorted((vi*a)%q for vi in v)
    landers=[(pow(a,-1,q))%q, (q-pow(a,-1,q))%q]   # v with v*a ≡ ±1
    present=[L for L in landers if L in v or any((vi%q)==L for vi in v)]
    nearest=min(dist0(r,q) for r in res)
    print(f"  {nm} @ base {q}, a={a}: residues min-dist={nearest} => M contribution {Fraction(nearest,q)}={float(Fraction(nearest,q)):.4f}")
    print(f"     distance-1 landers (v≡±a^-1 mod {q}) = {sorted(set(landers))}; present in family? {sorted(set(present))}")

print()
print("="*74)
print("(T2) does 'base 14 caps at 1/14' hold for ALL families? (the lever's reach)")
print("="*74)
for nm,v in [("AP (full grid)",AP),("avoids ±1 mod 14 (non-cov)",avoid),("covering-ish",covfam)]:
    b14,a14=best_at_base(v,14)
    cov=covering(v)
    tag=""
    if b14>Fraction(1,14): tag="  <-- BEATS 1/14 at base 14! (lever's cap is NOT universal)"
    if b14==0: tag="  <-- 0 (a runner ≡0 mod 14: base 14 BLOCKED)"
    if b14==Fraction(1,14): tag="  <-- exactly 1/14 (pinned: full grid, landers ±1 present)"
    print(f"  {nm:28s} cov={str(cov):5s} best@base14 = {str(b14):>6s}={float(b14):.4f} (a={a14}){tag}")
print("  => the lever pins EXACTLY 1/14 only for the FULL-GRID (AP) family. Others beat it (non-cov,")
print("     dispatched anyway) or are blocked (covering). So it is an EQUALITY statement, not a")
print("     universal lower bound. 'Base 14 caps at 1/14' is FALSE off the full grid.")

print()
print("="*74)
print("(T3) antipodal parity: good-multiplier set is ι-closed (a↦q−a) => size EVEN")
print("="*74)
for nm,v,q in [("deep well",deep,183),("AP",AP,13)]:
    bstar,_=best_at_base(v,q)
    # good multipliers achieving the max
    good=[a for a in range(1,q) if Fraction(min(dist0((vi*a)%q,q) for vi in v),q)==bstar]
    iota_closed=all(((q-a)%q) in good for a in good)
    print(f"  {nm} @ base {q}: max margin {bstar}; #good multipliers={len(good)} (even? {len(good)%2==0}); ι-closed? {iota_closed}")
print("  => witnesses come in ι-pairs {a,q−a} (klein-S55: 'loneliness is antipodally paired').")
print("     This is a CONSTRAINT (parity), not a detector: it says the shape of holes, not that")
print("     none exist (drop-11 @ D=41 was lonely at radius 2, a perfect ι-pair — klein-S55).")

print()
print("="*74)
print("(T4) the LOWER bound is NOT the base-14 lever: M≥1/14 needs a GOOD base, and")
print("     for covering families that base is LARGE (the lander-dodge crux, open)")
print("="*74)
for nm,v in [("AP",AP),("deep well",deep)]:
    M,argq,arga=exact_M(v)
    print(f"  {nm}: M={M}={float(M):.4f} achieved at base {argq} (a={arga}), cov={covering(v)}")
# for the deep well: is base 183 really forced? best_q for q<183
belowbest=Fraction(0); bq=None
for q in range(14,183):
    b,_=best_at_base(deep,q)
    if b>belowbest: belowbest=b; bq=q
print(f"  deep well: best over ALL bases q<183 = {belowbest}={float(belowbest):.4f} (at q={bq}) < 14/183={float(Fraction(14,183)):.4f}")
print("  => the covering floor 14/183 is realized ONLY at base 183; NO base<183 reaches it.")
print("     THAT lander-exclusion-up-to-183 is the crux, and the base-14 odd-|Aut| lever gives NONE of it.")
print()
print("VERDICT: the odd-|Aut| tie-resolution lever is an EQUALITY/UNIQUENESS tool (why the full-grid")
print("AP is pinned at exactly 1/14), NOT a lower-bound engine. The inequality M≥1/14 lives in the")
print("lander-DODGE / lander-EXCLUSION count over coarser bases (open) or the global antipodal-degree")
print("(Borsuk-Ulam odd index, open); the parity is a shadow of that degree, not a detector.")
print("\ndone.")
