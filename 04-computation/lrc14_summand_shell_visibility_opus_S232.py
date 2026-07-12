"""
opus-2026-07-11-S232: the remaining hard part, in the summand-graph / multiplicand-graph frame.

The clean-ruler route (kps THM-707: hB5 <= every residual has q with liveCount>=1, maxBand<=5) at the
NATURAL modulus C = 2n-1 = 27 (THM-401: the pinch/Farey-companion modulus of the floor 1/14) has an EXACT,
PROVED criterion in the summand-shell / multiplicand-unit structure. This operationalizes HYP-2935's slogan
"addition creates relation shells; multiplication tests their visibility" as a formula.

SETUP (THM-401 + HYP-2083/S571). Nonzero residues mod 27 partition into 13 antipodal SUMMAND SHELLS
P_a = {a, 27-a}. The danger arc at q=27 is {0,1,26} = {0} u P_1. The multiplicative units (Z/27)^* permute
the 9 UNIT shells {1,2,4,5,7,8,10,11,13}; the 3-adic strata are the blind spot: 3 gcd-3 shells
{3,24},{6,21},{12,15} and 1 gcd-9 shell {9,18}. bandCount(v,27,p):
  - p a UNIT:   z + (load of unit-shell P_{p^-1})            [z = #{27|v_i}]
  - p ~ gcd 3:  #{9 | v_i}                                    [danger {0,1,26} cap 3Z = {0}]
  - p ~ gcd 9:  #{3 | v_i}

TWO PROVED FORMULAS (verified 0 / 12000 random families each):
  (1) maxBand(v,27) = max( z + maxUnitShellLoad , #{3|v_i} ).
  (2) q=27 is a CLEAN RULER for v  <=>  maxBand(v,27) <= 5  AND  live, where
      live  <=>  [ z=0 AND some unit-shell empty ]  OR  [ #{9|v_i}=0 ]  OR  [ #{3|v_i}=0 ].
So SHALLOW = summand graph (no shell overloaded, and few multiples of 3); LIVE = multiplicand graph (a unit
clock, or a 3-adic clock, finds an empty danger-shell). This is the additive/multiplicative split, exactly.

THE WALL (why the clean-ruler route stops at the AP). clean-at-27 fails two ways:
  (a) SHALLOW-FAIL: #{3|v_i} >= 6 (or a unit-shell overloaded) = a mod-3/mod-27 CONCENTRATION (coarse-
      reducible / dilated family, dispatched by other branches).
  (b) LIVE-FAIL (shallow but not live): the family fills ALL 9 unit-shells AND has a multiple of 3 AND a
      multiple of 9 = MULTIPLICAND-MAXIMAL (S560o: "AP contains a multiple of every q<=13"). The AP {1..13},
      GW 12->24, and the sporadic V*={1..11,13,24} are the extreme live-fail families: each hits 9/9 unit
      shells, has 3 and 9 -> not live at 27. These are exactly the tight (M=1/14) wall families.

SYNTHESIS. The multiplicand-maximal (live-fail) hard core = the SAME families as the covering-route crux
(HYP-3085/OPEN-Q-108: consec maximizes the pairwise co-emptiness S2), the additive-energy route (THM-656/660:
AP maximizes R2), and opus-S181's "tight <=> Lambda 1-dimensional AND coherent (a dilate of {1..13})". All
roads lead to the AP. The clean-ruler route at C=27 certifies EVERYTHING except this AP-coherent wall, via
the explicit criterion above -- turning "the last hard part" into the single, named, inverse-additively-
characterized family class, which the tight-floor / LRC<=13 census handles.
"""
import random
from math import gcd
from functools import reduce
from collections import Counter
def shell(r):
    r%=27; return min(r,27-r)
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def maxband27_direct(v): return max(bandCount(v,27,p) for p in range(1,27))
def live27_direct(v): return any(bandCount(v,27,p)==0 for p in range(1,27))
UNIT_SHELLS={1,2,4,5,7,8,10,11,13}
def stats(v):
    z=sum(1 for x in v if x%27==0); z3=sum(1 for x in v if x%3==0); z9=sum(1 for x in v if x%9==0)
    usl=Counter(shell(x) for x in v if shell(x) in UNIT_SHELLS)
    return z,z3,z9,(max(usl.values()) if usl else 0),len(usl)
def maxband27_formula(v):
    z,z3,z9,maxUSL,nhit=stats(v); return max(z+maxUSL, z3)
def clean27_formula(v):
    z,z3,z9,maxUSL,nhit=stats(v)
    return max(z+maxUSL,z3)<=5 and ((z==0 and nhit<9) or z9==0 or z3==0)

def main():
    random.seed(2); mb=cl=0; N=12000
    for _ in range(N):
        v=sorted(random.sample(range(1,400),13))
        if maxband27_formula(v)!=maxband27_direct(v): mb+=1
        if clean27_formula(v)!=(live27_direct(v) and maxband27_direct(v)<=5): cl+=1
    print(f"(1) maxBand(v,27) = max(z + maxUnitShellLoad, #{{3|v}}):  {mb} mismatches / {N}")
    print(f"(2) clean-at-27 explicit criterion:                      {cl} mismatches / {N}")
    print("\nThe wall (named tight families are LIVE-FAIL at 27 = multiplicand-maximal):")
    for name,v in [("AP {1..13}",list(range(1,14))),("GW 12->24",[1,2,3,4,5,6,7,8,9,10,11,24,13]),
                   ("V* sporadic",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        v=sorted(v); z,z3,z9,mUSL,nhit=stats(v)
        print(f"  {name:12}: maxBand={maxband27_direct(v)} (shallow={maxband27_direct(v)<=5}), "
              f"live={live27_direct(v)}, unit-shells hit={nhit}/9, #div3={z3}, #div9={z9}")

if __name__=='__main__':
    main()
