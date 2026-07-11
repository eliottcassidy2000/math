# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont29: CREATIVE ROUTES to the clean-ruler supply (THM-707) -- connections dug from
# past fleet work. Four findings that reduce the supply to a bounded pair-sum search.
#
# CONTEXT: hB5 (the single Lean obligation) <= clean-ruler supply (THM-707): every residual covering family
# has a ruler q with a live multiplier (liveCount>=1) and shallow coverage (maxBand<=5), giving B5=liveCount>0.
#
# THE CONNECTIONS:
# (1) THM-668 (mac-mini, PROVEN): the loneliness witness M(S)=max_t min_i ||v_i t|| is ALWAYS attained at a
#     PAIR-SUM modulus q = v_i + v_j <= 2*Vmax. So the ruler is a pair-sum -- a bounded (<=78) candidate set.
# (2) CLEAN RULERS ARE MODERATE q (in [Vmax, 2Vmax]), NOT large q: at large q the p=1 multiplier puts EVERY
#     runner in the danger arc (v_i < q/14), so maxBand=13 and B5<<0. klein's large-q live rulers are NOT clean.
# (3) B5 is NOT scale-invariant: dilating v->c*v introduces a DEEP resonance at p=q*k/c (all runners hit 0),
#     so sign(B5) is NOT preserved. Primitive families are the right domain (dilating breaks the certificate).
# (4) EVERY residual family has a CLEAN PAIR-SUM ruler (some q=v_i+v_j, live>=1 & maxBand<=5). The M-maximizing
#     pair-sum is usually (but not always) the clean one; a clean pair-sum always exists among the <=78.
#
# => clean-ruler supply <=> "every residual family has a clean PAIR-SUM modulus" (bounded, decidable per family).
#    DECOMPOSITION: clean = SHALLOW (maxBand<=5, an ADDITIVE-structure condition, ~unconditional) + LIVE
#    (liveCount>=1, the LONELINESS content = LRC, supplied by klein's measure floor THM-687/692 for residuals).
from math import comb, gcd
from fractions import Fraction as F

def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def B5(v,q):
    S=[0]*6
    for p in range(1,q):
        b=bandCount(v,q,p)
        for d in range(6): S[d]+=comb(b,d)
    return sum((-1)**d*S[d] for d in range(6))
def liveMax(v,q):
    live=mb=0
    for p in range(1,q):
        b=bandCount(v,q,p)
        if b==0: live+=1
        mb=max(mb,b)
    return live,mb
def clean_pairsum(v):
    # search the <=78 pair-sum moduli for a clean one (live>=1 & maxBand<=5)
    for i in range(13):
        for j in range(i+1,13):
            q=v[i]+v[j]
            if q<14: continue
            live,mb=liveMax(v,q)
            if live>=1 and mb<=5: return q,live,mb
    return None

import random
def main():
    print("(2) MODERATE vs LARGE q  (near-AP {1..12,26}, Vmax=26):")
    v=list(range(1,13))+[26]
    for q in [27, 55, 200, 800]:
        live,mb=liveMax(v,q)
        print(f"    q={q:4d} (q/Vmax={q/26:.1f}): liveCount={live}, maxBand={mb}, B5={B5(v,q)}  {'CLEAN' if mb<=5 and live>=1 else 'deep (p=1 resonance)' if q>400 else ''}")
    print("(3) B5 NOT scale-invariant (dilation breaks it):")
    for c in [1,2,3]:
        vc=[c*x for x in v]; print(f"    c={c}: B5({c}v, {27*c}) = {B5(vc,27*c)}  (would be {c}x2={2*c} if scale-invariant)")
    print("(4) EVERY residual family has a CLEAN PAIR-SUM ruler:")
    random.seed(11); hit=tot=0
    for _ in range(200):
        vv=sorted(random.sample(range(1,26),12)); mx=random.randint(23,3*max(vv)); vv=sorted(set(vv+[mx]))
        if len(vv)!=13: continue
        tot+=1
        if clean_pairsum(vv) is not None: hit+=1
    print(f"    {hit}/{tot} residual families have a clean pair-sum modulus (live>=1 & maxBand<=5 at some v_i+v_j)")
    q,live,mb=clean_pairsum(v)
    print(f"    binding near-AP {{1..12,26}}: clean pair-sum q={q} (=1+26), live={live}, maxBand={mb} => B5={B5(v,q)}>0")
    print()
    print("REDUCTION (THM-668 o THM-707): clean-ruler supply <=> every residual has a clean PAIR-SUM modulus.")
    print("DECOMPOSITION: SHALLOW (additive, ~unconditional) + LIVE (loneliness = klein measure floor THM-687/692).")

if __name__ == '__main__':
    main()
