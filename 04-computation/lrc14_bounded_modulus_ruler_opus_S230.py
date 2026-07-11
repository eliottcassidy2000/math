"""
opus-2026-07-11-S230: the pair-sum shallow lemma is FALSE unconditionally; the correct clean-ruler
instrument for the hard core is a BOUNDED SMALL MODULUS (diameter-free), not a pair-sum (∝ diameter).

CONTEXT. kps THM-707 reduces hB5 (the single LRC(14) obligation) to: every RESIDUAL covering family has a
CLEAN RULER q (liveCount(q)>=1 and maxBand(v,q)<=5). kps's instrument was the PAIR-SUM q=v_a+v_b; kps
cont.32 flagged "some pair-sum has maxBand<=5" as the remaining rigorous gap ("LRC-hard").

THE DISSOCIATION/ANTI-CONCENTRATION ANGLE (owner-requested) redirects the instrument:

  (1) PAIR-SUM IS FALSE UNCONDITIONALLY. Because q=v_a+v_b ~ 2*Vmax, the number of expected 6-clusters in
      the danger arc grows ∝ Vmax; the min-over-pairsums maxBand climbs with diameter and exceeds 5.
      Verified: on prime-rich dissociated residuals, pair-sum fails (min-over-pairsums maxBand>=6) for
      0,0,1,3,36 of 120 families as Vmax ~ 60,120,250,500,900. NOT just coarse-reducible families -- genuine
      dissociated residuals at large Vmax. THM-701's peel bounds the RATIO (no far element = w<90191*Σe'),
      NOT Vmax, so its bounded-spread core CAN have large Vmax -> pair-sums are insufficient there.

  (2) BOUNDED SMALL MODULI ARE THE RIGHT INSTRUMENT. Cleanness at q depends only on {v_i mod q}, so a
      bounded q is DIAMETER-FREE. Verified: every dissociated prime-rich residual has a clean ruler with
      q<=37, and the bound does NOT grow with diameter (max smallest-clean-q = 37 over Vmax up to 5000, 0
      failures). A fixed 8-element set {8,17,19,23,27,29,35,37} is a clean ruler for 100% of a 2000-family
      diverse pool.

  (3) THE PROVED SMALL-PRIME CRITERION (the anti-concentration statement, made exact). For q in {17,19,23}
      the danger set is EXACTLY {0,1,q-1}={0,±1} (true for all q in [15,28]). With z=#{i:v_i≡0 mod q} and
      c_r=#{i:v_i≡±r mod q} for r=1..(q-1)/2:
          maxBand(v,q) = z + max_r c_r,   liveCount(v,q)=2*#{r:c_r=0} if z=0 else 0.
      Hence q is a CLEAN RULER  <=>  z=0  AND  max_r c_r <= 5  AND  some c_r = 0.
      (Proof: q prime => x->px bijects Z/q fixing 0; v_i p ≡ ±1 <=> v_i ≡ ±p^{-1}; as p ranges, p^{-1}
      ranges over all nonzero residues, each fold-class twice.) Verified 0 mismatches / 90000.
      => "max_r c_r <= 5" = NO antipodal ±-pair mod q holds >=6 of the 13 speeds -- an anti-concentration
      (Sidon-type) condition a dissociated family satisfies for most small primes. The remaining OPEN
      statement (cleanly posed, diameter-free): every dissociated prime-rich residual has SOME small prime
      q in the window with folded multiplicities <=5 and an empty ±-class.
"""
import random
from math import gcd
from functools import reduce

def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def maxband(v,q): return max(bandCount(v,q,p) for p in range(1,q))
def livecount(v,q): return sum(1 for p in range(1,q) if bandCount(v,q,p)==0)
def is_clean(v,q): return livecount(v,q)>=1 and maxband(v,q)<=5
def primitive(v): return reduce(gcd,v)==1
def prime_rich(v): return all(any(x%p==0 for x in v) for p in [2,3,5,7,11,13])
def longest_AP(v):
    s=set(v); best=1
    for a in v:
        for d in range(1,max(v)//2+1):
            L=1; x=a+d
            while x in s: L+=1; x+=d
            if L>best: best=L
    return best
def is_residual(v): return primitive(v) and prime_rich(v) and longest_AP(v)<=7
def minmax_pairsum(v):
    ps={v[i]+v[j] for i in range(13) for j in range(i+1,13) if v[i]+v[j]>=14}
    return min(maxband(v,q) for q in ps)
def smallest_clean(v,Q):
    for q in range(8,Q+1):
        if is_clean(v,q): return q
    return None
def crit_1719_23(v,q):
    """proved criterion for q in {17,19,23}: z=0, max fold-mult <=5, some fold empty."""
    if any(vi%q==0 for vi in v): return False
    c={}
    for vi in v:
        r=vi%q; r=min(r,q-r); c[r]=c.get(r,0)+1
    cr=[c.get(r,0) for r in range(1,(q-1)//2+1)]
    return max(cr)<=5 and min(cr)==0

def main():
    random.seed(11)
    print("(1) PAIR-SUM shallow degrades with diameter (prime-rich dissociated residuals):")
    for hi in [60,250,900]:
        tested=psfail=tries=0
        while tested<120 and tries<200000:
            tries+=1
            v=sorted(random.sample(range(1,hi+1),13))
            if v[-1]<hi*0.6 or not is_residual(v): continue
            tested+=1
            if minmax_pairsum(v)>=6: psfail+=1
        print(f"    Vmax~{hi:4d}: pair-sum fails (min-over-pairsums maxBand>=6) for {psfail}/{tested}")

    print("\n(2) BOUNDED small-modulus ruler is diameter-free (smallest clean q in {8..80}):")
    random.seed(2026); om=0
    for hi in [100,800,5000]:
        tested=mx=none=tries=0
        while tested<200 and tries<600000:
            tries+=1
            v=sorted(random.sample(range(1,hi+1),13))
            if v[-1]<hi*0.5 or not is_residual(v): continue
            tested+=1
            sc=smallest_clean(v,80)
            if sc is None: none+=1
            else: mx=max(mx,sc); om=max(om,sc)
        print(f"    Vmax~{hi:5d}: max(smallest clean q)={mx}, #none<=80={none} (of {tested})")
    print(f"    => bounded window: overall max smallest-clean-q = {om} (does not grow with Vmax)")

    print("\n(3) PROVED small-prime criterion vs is_clean over q in {17,19,23}:")
    random.seed(5); mism=checks=0
    for _ in range(30000):
        v=sorted(random.sample(range(1,200),13))
        for q in [17,19,23]:
            checks+=1
            if is_clean(v,q)!=crit_1719_23(v,q): mism+=1
    print(f"    {mism} mismatches / {checks} checks  (criterion: z=0 & max fold-mult<=5 & some fold empty)")

if __name__=='__main__':
    main()
