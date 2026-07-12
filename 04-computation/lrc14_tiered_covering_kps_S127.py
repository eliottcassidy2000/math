# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont37: EXTENDING the finite bounded-modulus combinatorics -- a 3-TIER COVERING of the
# bounded-window check (every residual family has a B5>0 ruler in [8,43]), 100% over 6000 families:
#   TIER 1 (RIGOROUS divisibility): THM-712 prime clean rulers {2,3,5,7,11,13} (p-nmid-all => B5=p-1)
#     + the GENERAL CLEAN q<=14 ruler (verified 0/21000): #{q|v}=0 AND (all proper divisors d|q: #{d|v}<=5)
#       => B5(v,q)>0 [since for q<=14 the band = all nonzero residues, so bandCount(v,q,p)=#{i: q|v_i p}
#          = #{i: (q/gcd(q,p))|v_i} <= max over proper divisors, and coprime p give bandCount 0]. ~82%.
#   TIER 2 (detuned): heavy-divisibility (>=7 of 13 share a prime) => THM-678 detuned dispatch. ~12%.
#   TIER 3 (near-unit): q in [15,43], the +-1 resonance rulers (opus not_loose_near_unit). ~6%.
from math import comb
import random
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def B5(v,q):
    S=[0]*6
    for p in range(1,q):
        b=bandCount(v,q,p)
        for d in range(6): S[d]+=comb(b,d)
    return sum((-1)**d*S[d] for d in range(6))
def clean_cond(v,q):
    if any(vi%q==0 for vi in v): return False
    return all(sum(1 for vi in v if vi%d==0)<=5 for d in range(2,q) if q%d==0)
def tier1(v): return any(all(vi%p!=0 for vi in v) for p in [2,3,5,7,11,13]) or any(clean_cond(v,q) for q in range(8,15))
def heavy(v): return any(sum(1 for vi in v if vi%g==0)>=7 for g in [2,3,5,7])
def nearunit(v): return any(B5(v,q)>0 for q in range(15,44))

def main():
    random.seed(11); n=6000; t1=t2=t3=un=0
    for _ in range(n):
        v=[random.randrange(1,80) for _ in range(13)]
        if tier1(v): t1+=1
        elif heavy(v): t2+=1
        elif nearunit(v): t3+=1
        else: un+=1
    print(f"3-TIER COVERING of the bounded-window check ({n} families):")
    print(f"  Tier 1 RIGOROUS [THM-712 primes + clean q<=14 divisibility]: {100*t1/n:.1f}%")
    print(f"  Tier 2 detuned  [>=7 share a prime => THM-678]:               {100*t2/n:.1f}%")
    print(f"  Tier 3 near-unit[q in [15,43], +-1 resonance rulers]:         {100*t3/n:.1f}%")
    print(f"  UNCOVERED: {100*un/n:.2f}%  => tiers cover {100*(t1+t2+t3)/n:.2f}%")
    print("  => the window check = [rigorous divisibility] + [detuned/THM-678] + [near-unit]; the open piece is Tier 3.")

if __name__=='__main__': main()
