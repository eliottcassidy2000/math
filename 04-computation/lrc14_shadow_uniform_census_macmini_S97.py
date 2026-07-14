#!/usr/bin/env python3
"""mac-mini-S97b: does the k<=13 shadow-witness close the WHOLE covering case? Broad census of covering
13-sets {1}uC (single-killer dilates, multi-killer, near-AP, random). For each, does SOME a/k (k<=13,
a/k in middle) give an exact-rational shadow witness (verified lonely)? If 100%, the shadow route is a
UNIFORM elementary alternative to disc_v; the open part = the uniform PROOF (residue-pattern)."""
from fractions import Fraction as F
from math import gcd
import random
c114=F(1,14)
def lonely_at(S,t):
    for c in S:
        r=(c*t)%1
        if min(r,1-r)<c114: return False
    return True
def shadow_interval(S,a,k):
    lo=F(0); hi=F(1,2)
    for c in S:
        r=(c*a)%k
        if r==0: lo=max(lo,F(1,14*c)); hi=min(hi,F(13,14*c))
        else:
            s=r if r<=k//2 else r-k; base=F(abs(s),k)
            hi=min(hi, (F(13,14)-base)/c if s>0 else (base-c114)/c)
        if lo>=hi: return None
    return (lo,hi)
def has_shadow_witness(C):
    for k in range(2,14):
        for a in range(1,k):
            if gcd(a,k)!=1 or not (c114<=F(a,k)<=F(13,14)): continue
            iv=shadow_interval(C,a,k)
            if iv and lonely_at([1]+C,(F(a,k)+ (iv[0]+iv[1])/2)):
                return (k,a,iv)
    return None
def is_cov13(S):  # S a 13-set containing 1, covering 2..14, primitive
    return len(set(S))==13 and 1 in S and all(any(v%q==0 for v in S) for q in range(2,15))
def gen_random_covering(rng, vmax=200):
    # fill each modulus 2..14 with a random multiple, pad with random speeds, force size 13 incl 1
    S={1}
    for q in range(2,15):
        if any(v%q==0 for v in S): continue
        mult=q*rng.randint(1, vmax//q); S.add(mult)
    while len(S)<13: S.add(rng.randint(2,vmax))
    S=set(list(S)[:13])
    while len(S)<13: S.add(rng.randint(2,vmax))
    return sorted(S)

fams=[]
# single-killer dilates {1..12,182m} and dilated cores
for m in range(1,9): fams.append((f"1..12,182*{m}", list(range(2,13))+[182*m]))
for cc in [2,3,5,7,11]:
    core=[cc*i for i in range(1,13)]; 
    if 1 in core: continue
    fams.append((f"{cc}*(1..12)+182*{cc}", [cc]+core[1:]+[182*cc]) if cc*1!=1 else (None,None))
# multi-killer (kps families)
for mk in [[2,3,4,5,6,7,8,9,10,11,13,84],[2,3,4,5,6,7,8,9,10,13,22,84],[2,3,4,5,6,7,8,9,13,14,45,88]]:
    fams.append((f"MK {mk[:3]}..", mk))
# near-AP drop-x + far
for x in range(2,8):
    base=[v for v in range(2,15) if v!=x]  # {2..14}\{x}, 12 elts, {1}u this = near-AP
    fams.append((f"drop-{x}", base))
# random covering
rng=random.Random(20260714)
for i in range(120):
    S=gen_random_covering(rng)
    if is_cov13(S): fams.append((f"rand{i}", [v for v in S if v!=1]))

tot=0; win=0; kdist={}; fails=[]
for name,C in fams:
    if C is None: continue
    if not all(any(v%q==0 for v in C) for q in range(2,15)): continue  # cluster must cover (1 covers nothing)
    tot+=1; w=has_shadow_witness(C)
    if w: win+=1; kdist[w[0]]=kdist.get(w[0],0)+1
    else: fails.append((name,C))
print(f"covering families tested: {tot}")
print(f"with a k<=13 shadow witness (verified lonely): {win}/{tot} = {100*win/tot:.1f}%")
print(f"min-k distribution: {dict(sorted(kdist.items()))}")
if fails:
    print(f"FAILURES (no k<=13 shadow witness) -- the hard cases: {len(fails)}")
    for name,C in fails[:8]: print(f"   {name}: {C}")
else:
    print("ZERO failures: EVERY covering family tested has a k<=13 shadow witness.")
    print("=> the shadow-witness route (exact residue condition) closes the covering case empirically,")
    print("   INCLUDING isolated-far (deep well) -- a UNIFORM elementary alternative to disc_v.")
    print("   OPEN: the uniform PROOF (some k<=13 works for every covering residue pattern) = klein-S299's")
    print("   resonance-grid equidistribution, now with my EXACT decidable residue-mod-k condition.")
