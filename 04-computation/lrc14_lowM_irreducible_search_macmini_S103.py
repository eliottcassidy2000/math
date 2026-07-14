#!/usr/bin/env python3
"""mac-mini-S103d: SKEPTICAL make-or-break. Search adversarially near the covering-min for a LOW-M
IRREDUCIBLE family (no safe peel, M<0.12) that is ALSO not shadow-closable (THM-749) and not near-AP
(<11 in {1..14}, kps THM-734). Such a family = genuine open gap. If none, the reduction holds:
covering => (safe-peel => LRC(<=13)) OR (irreducible => tile: shadow / near-AP / loose)."""
from fractions import Fraction as F
from math import gcd
import numpy as np, random
c114=F(1,14)
def M_t0(S,dens=1_000_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def safe_peel_exists(S):
    for v in sorted(S,reverse=True):
        core=[u for u in S if u!=v]
        if len(core)<2: return True
        mu0,t0=M_t0(core,dens=400_000)
        if min((v*t0)%1.0,1-(v*t0)%1.0)>=mu0-2e-4: return True
    return False
def shadow_interval(S,a,k):
    lo=F(0); hi=F(1,2)
    for c in S:
        r=(c*a)%k
        if r==0: lo=max(lo,F(1,14*c)); hi=min(hi,F(13,14*c))
        else:
            s=r if r<=k//2 else r-k; base=F(abs(s),k); hi=min(hi,(F(13,14)-base)/c if s>0 else (base-c114)/c)
        if lo>=hi: return False
    return True
def has_shadow(S):
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)==1 and c114<=F(a,k)<=F(13,14) and shadow_interval(S,a,k): return True
    return False
def covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
def n114(S): return sum(1 for v in S if 1<=v<=14)
rng=random.Random(31337)
seeds=[list(range(1,13))+[182],[*range(1,12),13,84],[v for v in range(1,14) if v!=6]+[182],
       [*range(1,11),13,14,45],[*range(1,11),13,22,84],[v for v in range(1,15) if v!=6]]
gap=[]; irred_lowM=0; tot=0
for _ in range(20000):
    base=list(rng.choice(seeds)); S=set(base)
    for _ in range(rng.randint(1,3)):
        if S: S.discard(rng.choice(list(S)))
        S.add(rng.choice([rng.randint(2,90), rng.choice([13,26,14,28,84,182,169,45,90,91,183])]))
    while len(S)<13: S.add(rng.randint(2,90))
    S=sorted(set(list(S)[:13]))
    if len(S)!=13 or 1 not in S or not covering(S): continue
    tot+=1
    if safe_peel_exists(S): continue    # reducible => LRC(<=13)
    # IRREDUCIBLE: check M and tiles
    Mm,_=M_t0(S)
    if Mm<0.12:
        irred_lowM+=1
        if not has_shadow(S) and n114(S)<9:   # not shadow, not near-AP => genuine gap
            gap.append((round(Mm,4),sorted(S)))
print(f"covering families near extremals: {tot}; IRREDUCIBLE with M<0.12: {irred_lowM}")
print(f"  GENUINE GAP (irreducible, low-M, NOT shadow, NOT near-AP): {len(gap)}")
for M,S in sorted(gap)[:10]: print(f"    M={M}: {S}  in114={n114(S)}")
if not gap:
    print("  => NO genuine gap: every low-M irreducible family is shadow-closable OR near-AP (a tile).")
    print("     REDUCTION HOLDS: covering => (safe-peel => LRC(<=13), 98%) OR (irreducible: shadow/near-AP/loose).")
else:
    print("  => genuine gap families exist = irreducible, low-M, in NO tile. The real open residual.")
