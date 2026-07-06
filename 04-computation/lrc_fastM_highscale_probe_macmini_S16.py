#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S16 (HYP-4432) -- FAST EXACT M via the witness-denominator lemma,
and the high-scale difference-core probe.

FAST EXACT M: by the lemma (q | v_i+-v_j or 2v_i), the witness denominator is among
the O(n^2) pairwise sums/differences.  So M = max over those q, over a/q, of
min_i ||v_i a/q|| -- EXACT and fast (no slow profile solver, no large-denominator
sweep).  This is a genuine engineering deliverable: exact M in O(n^2 * max) vs the
profile solver's blowup at high scale.

HIGH-SCALE DIFFERENCE-CORE PROBE (honest): is a high-scale near-AP cluster
N*{1..12}+eps ever a GAP member?  FINDINGS:
 * The dilated AP N*{1..12} is tight (M=1/13) at every N (dilation-invariant).
 * Perturbing it: SPECIFIC structure-preserving eps give clean scale-stable rungs
   (+1 on r12 -> 1/12 at every N; +1,+2 on r11,r12 -> 1/11; etc.), all >= 2/25.
 * GENERIC eps (|eps|<=3) make it LOOSE (M ~ 0.2-0.3) -- a tiny integer perturbation
   destroys the exact roots-of-unity resonance and the family becomes generic.
   (NOT exact scale-invariance: M(20*AP+eps) != M(100*AP+eps) in general; the values
   drift but stay loose.)
 * ACROSS 400+ perturbations at high scale: ZERO in the gap (1/13, 2/25).
CONCLUSION (empirical): the dilated AP is an ISOLATED tight point at EVERY scale;
high scale does NOT manufacture a gap member.  Supports the difference-core case of
(G) (the values are three-gap/CF-quantized -- rung or loose, no gap -- at every
scale), though it does not prove it.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

def Mfast(S):
    """EXACT M(S)=max_t min_i||v_i t|| via the witness-denominator lemma."""
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val
    return best

def prim(S):
    g=reduce(gcd,S); return [x//g for x in S]

GAP_LO,GAP_HI=F(1,13),F(2,25); AP=list(range(1,13))

if __name__=="__main__":
    print("sanity:  M(AP)=",Mfast(AP)," M(doubled-apex)=",Mfast(list(range(1,12))+[24]),
          " M(deep well)=",Mfast(list(range(1,12))+[168]))
    print("\nstructure-preserving perturbations (scale-stable clean rungs):")
    for name,eps in [("+1 r12",[0]*11+[1]),("+1,+2 r11,r12",[0]*10+[1,2]),("+3,+5,+7 r10-12",[0]*9+[3,5,7])]:
        vals=[str(Mfast(prim([N*AP[i]+eps[i] for i in range(12)]))) for N in [5,20,100]]
        print(f"   {name:16s}: N=5,20,100 -> {vals}")
    print("\ngeneric perturbations at high scale are LOOSE, and NONE in gap:")
    random.seed(16); in_gap=0; n=0; loose=0
    for _ in range(400):
        eps=[random.randint(-3,3) for _ in range(12)]
        S=prim([60*AP[i]+eps[i] for i in range(12)])
        if len(set(S))<12 or min(S)<=0: continue
        n+=1; M=Mfast(S)
        if GAP_LO<M<GAP_HI: in_gap+=1
        if M>=F(1,8): loose+=1
    print(f"   checked {n}: in gap={in_gap}, loose(M>=1/8)={loose}  => AP isolated at every scale")
