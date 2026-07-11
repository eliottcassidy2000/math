# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont31: coverage of the GENERAL prime clean ruler (Lean b5_pos_of_prime_ndvd).
# For any prime P<=13 with no speed divisible by P: q=P is a perfectly clean ruler (bandCount=0 everywhere,
# liveCount=P-1) => B5(v,P)=P-1>0 => hB5 discharged. Generalizes opus-S227 THM-709 (q=13) to {2,3,5,7,11,13}.
import random
from functools import reduce
from math import gcd, comb
primes=[2,3,5,7,11,13]
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def B5(v,q):
    S=[0]*6
    for pp in range(1,q):
        b=bandCount(v,q,pp)
        for d in range(6): S[d]+=comb(b,d)
    return sum((-1)**d*S[d] for d in range(6))
def dispatched_by(v):  # the smallest prime <=13 dividing NO speed (the clean ruler), or None
    for P in primes:
        if all(x%P!=0 for x in v): return P
    return None
def generic_resid(v):
    if len(set(v))!=13 or max(v)<23 or reduce(gcd,v)!=1: return False
    return all(sum(1 for x in v if x%P==0)<=6 for P in primes)

def main():
    print("VERIFY the lemma: q=P (prime<=13, P-nmid-all) => B5=P-1:")
    v=[1,2,4,8,9,16,17,23,25,27,29,31,37]  # avoids 5,7,11,13 (not 2,3)
    for P in primes:
        av=all(x%P!=0 for x in v)
        if av: print(f"    P={P}: B5(v,{P})={B5(v,P)} (= P-1={P-1}) [maxBand 0, liveCount {P-1}]")
    random.seed(1); disc=tot=0; byP={P:0 for P in primes}
    for _ in range(6000):
        w=sorted(random.sample(range(1,50),13))
        if not generic_resid(w): continue
        tot+=1; P=dispatched_by(w)
        if P: disc+=1; byP[P]+=1
    print(f"\nCOVERAGE over {tot} generic residuals:")
    print(f"  prime dispatch discharges {disc}/{tot} = {100*disc/tot:.0f}%  (opus q=13 alone ~ {100*byP[13]/tot:.0f}% + smaller primes)")
    print(f"  per-prime first-hit: "+', '.join(f'{P}:{100*byP[P]/tot:.0f}%' for P in primes))
    print(f"  remaining CORE (prime-rich, every prime<=13 divides some speed): {100*(tot-disc)/tot:.0f}% -- AP-like, handed to composite/pair-sum/moment-ladder")

if __name__=='__main__': main()
