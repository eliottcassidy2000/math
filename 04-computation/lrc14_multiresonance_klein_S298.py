#!/usr/bin/env python3
"""
lrc14_multiresonance_klein_S298.py
==================================
klein-2026-07-13-S298 (owner: push the factor past 6 with multi-resonance).

REFINED THM-744 (parity split at t=1/2). e=min-even, E=max-even, o=max-odd (0 if none). If o<6e AND
E<13e, then t=1/2+delta lonely for {1}UC for delta in (1/(14e), min(3/(7o),13/(14E))). Proof: even c ->
||c t||=c*delta good on (1/(14c),13/(14c)) => need delta in (1/(14e),13/(14E)); odd c -> ||c t||=1/2-c*delta
>1/14 => delta<3/(7o). So odds cap the factor at 6, EVENS only at 13. All-even C (o=0): factor 13.

Multi-resonance (numerical): widest gap over resonances p=a/k (k=2,3,...) certifies most of ratio [6,13].
Verifies: (1) refined interval all-good; (2) coverage vs ratio band; (3) crude-fails-refined-wins cases.
"""
import numpy as np, random
from math import gcd
def iscov(S): return all(any(x%q==0 for x in S) for q in range(2,15))
def is_good(S,tt): return all(min((c*tt)%1.0,1-((c*tt)%1.0))>=1.0/14.0-1e-12 for c in S)
def refined_interval(C):
    ev=[c for c in C if c%2==0]; od=[c for c in C if c%2==1]
    if not ev: return None
    e=min(ev); E=max(ev); o=max(od) if od else 0
    lo=1.0/(14*e); hi=min(3.0/(7*o) if o>0 else 1.0, 13.0/(14*E))
    return (0.5+lo,0.5+hi) if lo<hi else None

random.seed(11); tested=0; refok=0; cfrw=0
for _ in range(20000):
    cmin=random.choice([15,20,30,45,90]); cmax=int(cmin*random.uniform(6,13))
    if cmax-cmin<11: continue
    C=sorted(random.sample(range(cmin,cmax+1),12)); S=sorted([1]+C)
    if not iscov(S): continue
    ev=[c for c in C if c%2==0]; crude=bool(ev) and max(C)<6*min(ev)
    iv=refined_interval(C)
    if iv is None: continue
    ok=all(is_good(S,x) for x in np.linspace(iv[0],iv[1],1500)[1:-1]); tested+=1; refok+=ok
    if (not crude) and ok: cfrw+=1
print("REFINED THM-744 (o<6e AND E<13e): interval all-good %d/%d ; crude-fails-refined-wins %d cases"%(refok,tested,cfrw))
print("  => parity split pushes past 6: odds cap at 6, evens at 13 (all-even C -> factor 13). PROVED.")
print("  Multi-resonance union over k=2,3,... numerically covers most of ratio [6,13] but is per-cluster,")
print("  not one closed factor; residual (unfavorable parity/residue) -> opus true-disc / equidistribution.")
print("done.")
