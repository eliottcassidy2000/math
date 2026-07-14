#!/usr/bin/env python3
"""mac-mini-S103c: THE KEY REDUCTION. Safe-peel recursion preserves M; a terminal core of size <=12
(<=12 moving speeds) is a SETTLED LRC(<=13) instance (M>=1/13>1/14). So families reducing to size<=12
are CLOSED. The residual = IRREDUCIBLE families (size 13, NO safe peel). Count + characterize them:
are they the known extremals (near-AP/single-killer, already closed) or a new open set?"""
from fractions import Fraction as F
import numpy as np, random
def M_t0(S,dens=800_000):
    x=np.arange(1,dens)/dens; mn=None
    for v in S:
        r=(v*x)%1.0; d=np.minimum(r,1-r); mn=d if mn is None else np.minimum(mn,d)
    j=int(mn.argmax()); return mn[j],x[j]
def safe_peel_exists(S):
    for v in sorted(S,reverse=True):
        core=[u for u in S if u!=v]
        if len(core)<2: return v
        mu0,t0=M_t0(core,dens=400_000)
        vn=min((v*t0)%1.0, 1-(v*t0)%1.0)
        if vn>=mu0-2e-4: return v
    return None
def recurse_size(S):
    S=list(S)
    while len(S)>2:
        v=safe_peel_exists(S)
        if v is None: break
        S=[u for u in S if u!=v]
    return len(S)
def covering(S): return all(any(u%q==0 for u in S) for q in range(2,15))
def n_in_114(S): return sum(1 for v in S if 1<=v<=14)
rng=random.Random(7777)
tot=0; closed_le12=0; irred=[]
for _ in range(30000):
    S=set()
    for q in range(2,15):
        if any(u%q==0 for u in S): continue
        S.add(q*rng.randint(1,20))
    while len(S)<13: S.add(rng.randint(1,350))
    S=sorted(set(list(S)[:13]))
    if len(S)!=13 or 1 not in S or not covering(S): continue
    tot+=1
    sz=recurse_size(S)
    if sz<=12: closed_le12+=1     # LRC(<=13) settles the terminal core => M>=1/13
    else: irred.append(sorted(S)) # size 13, no safe peel = irreducible residual
    if tot>=600: break
print(f"covering families tested (all diameters): {tot}")
print(f"  reduce to size<=12 => CLOSED by safe-peel recursion + LRC(<=13): {closed_le12} ({100*closed_le12/tot:.1f}%)")
print(f"  IRREDUCIBLE (size 13, NO safe peel = residual): {len(irred)} ({100*len(irred)/tot:.1f}%)")
if irred:
    print("  irreducible families -- are they near-AP (>=11 in {1..14}, kps THM-734) or single-killer?")
    nAP=sum(1 for S in irred if n_in_114(S)>=11)
    print(f"    >=11 in {{1..14}} (kps near-AP tile): {nAP}/{len(irred)}")
    Ms=[M_t0(S)[0] for S in irred[:60]]
    print(f"    irreducible M: min={min(Ms):.4f} median={sorted(Ms)[len(Ms)//2]:.4f} (1/14={1/14:.4f})")
    for S in irred[:6]:
        Mm,_=M_t0(S); print(f"      M={Mm:.4f} in114={n_in_114(S)}: {S}")
print()
print("=> if the IRREDUCIBLE residual is small AND = the known tiles (near-AP kps / single-killer /")
print("   covering-min), then THM-751 safe-peel + LRC(<=13) + those tiles CLOSE the covering case.")
