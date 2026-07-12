"""
opus-2026-07-11-S239: the three-gap "disjunction" for the spread bulk is FAVORABLE, not a wall -- and it
UNIFIES the clean-ruler route with the density route (both reduce to "the AP is the unique good coverer").

Owner: work the three-gap disjunction. Reading mac-mini cont.44 (three-gap regularity = why the AP is the
best coverer; coverage advantage 24x k=7 .. 5x k=10) reframes my spread-bulk residual.

THE REFRAMING. Clearing at a modulus q (a lonely multiplier p, bandCount=0) = the danger arc is EMPTY at p =
the family FAILS to cover at p. mac-mini's three-gap coverage advantage: the AP (Steinhaus {k*alpha} orbit)
is the UNIQUE best coverer (covers 7-24x more than iid); every non-AP (spread) family covers like iid = BADLY.
So a bad coverer has MANY lonely/clearing multipliers. Divisor-complete FORCES spread (S237: 99% longest-AP<=7,
structural: needs mult of 8,9,11,13,14, incompatible with the tight AP). Hence:

  divisor-complete  =>  spread  =>  BAD coverer  =>  MANY clearing multipliers  =>  clears easily.

VERIFIED (window = non-14 moduli in [15,31]):
  - AP {1..13} (good coverer, bucket A): window clearing-count = 0. It clears at NO non-14 modulus -- only at
    multiples of 14 (t=1/14). The AP is the HARD coverer.
  - SPREAD divisor-complete (bad coverers): window clearing-count min=8, mean=32, max=56; each clears at
    min 2, mean ~7 of the 16 window moduli. 0/800 have window-clear-count <= 3. They clear ROBUSTLY.

RECONCILIATION (corrects S237/S238's "spread is the wall"): the spread bulk is the FAVORABLE case (bad
coverers, min 8 clearing multipliers). The genuinely hard object is the AP (good coverer, 0 window-clearing),
which is NOT divisor-complete (bucket A) and is dispatched by the elementary t=1/14 witness. kps's "window-hard
= spread" means the spread families clear at VARYING q (wide window, no small shortcut, S238) -- NOT that they
are hard to clear.

THE UNIFICATION. Both LRC(14) routes now reduce to the SAME crux:
  - DENSITY route (mac-mini/klein): the moment-ladder base needs "consec/AP is the extremal" -- the AP maximizes
    coverage (good coverer). Hard part = proving the extremal HAS {k*alpha} structure (the inverse theorem).
  - CLEAN-RULER route (kps/opus): the residual needs "divisor-complete (spread) clears" -- spread = bad coverer
    => clears. Hard part = proving spread => not-the-AP-good-coverer (the SAME inverse theorem).
So the clean-ruler and density routes have the SAME wall -- the AP inverse theorem ("the good coverer is the
{k*alpha} orbit") -- approached from opposite sides. The residual is the favorable, quantitatively-robust
(min 8) direction: spread => bad coverer => clears.
"""
from math import gcd
from functools import reduce
import random
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def clearing_count(v,q): return sum(1 for p in range(1,q) if bandCount(v,q,p)==0)
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def longest_AP(v):
    s=set(v); best=1
    for a in v:
        for d in range(1,max(v)//2+1):
            L=1;x=a+d
            while x in s:L+=1;x+=d
            if L>best:best=L
    return best
WIN=[q for q in range(15,32) if q%14!=0]

def main():
    ap=list(range(1,14))
    print(f"AP {{1..13}} (good coverer, bucket A): window clearing-count = {sum(clearing_count(ap,q) for q in WIN)} "
          f"(clears at NO non-14 modulus => hard coverer, dispatched by t=1/14)")
    random.seed(1); pool=[]; tries=0
    while len(pool)<800 and tries<300000:
        tries+=1
        v=sorted(random.sample(range(1,120),13))
        if primitive(v) and divisor_complete(v) and longest_AP(v)<=7: pool.append(v)
    wcs=[sum(clearing_count(v,q) for q in WIN) for v in pool]
    nq=[sum(1 for q in WIN if clearing_count(v,q)>0) for v in pool]
    print(f"SPREAD divisor-complete (n={len(pool)}, bad coverers): window clearing-count min={min(wcs)}, "
          f"mean={sum(wcs)//len(wcs)}, max={max(wcs)}; clears at min {min(nq)}, mean {sum(nq)//len(nq)} of {len(WIN)} window q")
    print(f"  #with window-clear-count<=3: {sum(1 for w in wcs if w<=3)}/{len(pool)} => spread clears ROBUSTLY")
    print("\n=> spread bulk is the FAVORABLE case; the AP (good coverer) is the wall, dispatched by t=1/14.")
    print("   Both LRC(14) routes reduce to the SAME AP inverse theorem (good coverer = {k*alpha} orbit).")

if __name__=='__main__':
    main()
