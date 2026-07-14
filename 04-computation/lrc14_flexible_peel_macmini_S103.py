#!/usr/bin/env python3
"""mac-mini-S103: does THM-751 with FLEXIBLE PEELING (peel ANY element, not just the largest) close the
unsafe/loose stratum? For each covering family, try peeling EACH element v; classify v at the core's
tight point as aligned / non-aligned-safe (clean => THM-751 reduces to M(core)) or unsafe. If EVERY
covering family has SOME clean peel (=> recurse), THM-751 + sieve closes the covering case elementarily."""
from fractions import Fraction as F
from math import gcd
import numpy as np, random
def M_t0(S,dens=1_200_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def peel_clean(S, v):
    """is peeling v CLEAN (aligned or non-aligned-safe) at the core's tight point?"""
    core=[u for u in S if u!=v]
    if len(core)<2: return True
    mu0,t0=M_t0(core, dens=600_000)
    r=(v*t0)%1.0; vn=min(r,1-r)
    return (vn < 1.5/max(v,1)) or (vn >= mu0-2e-4)   # aligned OR safe
def has_clean_peel(S): return any(peel_clean(S,v) for v in S)
def covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
rng=random.Random(2024)
# focus on the UNSAFE/loose stratum: high-diameter spread covering families
tot=0; all_have_clean=0; no_clean=[]
for _ in range(30000):
    S=set()
    for q in range(2,15):
        if any(u%q==0 for u in S): continue
        S.add(q*rng.randint(1,20))
    while len(S)<13: S.add(rng.randint(1,350))
    S=sorted(set(list(S)[:13]))
    if len(S)!=13 or 1 not in S or not covering(S): continue
    # only look at spread/unsafe: diameter > 20
    if max(S)/min([u for u in S if u>1]) < 8: continue
    tot+=1
    if has_clean_peel(S): all_have_clean+=1
    else: no_clean.append(sorted(S))
    if tot>=1200: break
print(f"spread covering families tested: {tot}")
print(f"with SOME clean (aligned/safe) peel: {all_have_clean} ({100*all_have_clean/max(tot,1):.1f}%)")
print(f"with NO clean peel (all peels unsafe -- genuine density-floor residual): {len(no_clean)}")
for S in no_clean[:6]:
    Mm,_=M_t0(S); print(f"   M={Mm:.4f}: {S}")
print()
if not no_clean:
    print("=> EVERY spread covering family has a clean peel => THM-751 + recursion + sieve closes the")
    print("   loose stratum too (peel a clean element, reduce to M(core), recurse to non-covering core).")
    print("   The covering case would be ELEMENTARY: THM-751 (any clean peel) + THM-366 sieve, no density floor.")
else:
    print("=> some families have NO clean peel = the genuine density-floor residual (opus THM-745/746).")
    print("   THM-751 flexible peeling reduces but does not eliminate the loose stratum; the no-clean-peel")
    print("   families are the irreducible core of opus's floor.")
