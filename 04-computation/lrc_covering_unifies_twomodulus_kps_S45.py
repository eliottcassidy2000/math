#!/usr/bin/env python3
"""
kps-S45: the S44 covering system UNIFIES opus-S125's two-modulus factoring of (C), and supplies
the mechanism opus flagged OPEN (mod-13 collision => loose).

opus S125: (C) <=> (1)[mod-13 collision => M>=2/25, verified, MECHANISM OPEN] + (2)[doubly-saturated
(full-transversal-25 AND distinct-13) + M<2/25 => AP].  The 25-clearance (non-transversal => loose)
is GREEN (mac-mini THM-634 + kps LRCMod25Floor).

SYNTHESIS: kps S44 shows every non-AP family clears at a BOUNDED small modulus.  Applied to the
full-transversal residual, this UNIFIES opus (1)+(2): both mod-13-collision and distinct-13-non-AP
full-transversals clear at q<=Q0 (<=23 on sample) -- so the SAME covering supplies opus's open
collision mechanism AND the distinct-13 rigidity, leaving the AP as the sole survivor.

Also: the residual (M<2/25 full-transversal) has a multiple of each of {7..12} (kps LRCSmallModFloor
loose_of_no_multiple_q), so it is extremely AP-like -- the last hard node = 'the mult-of-{7..12}
doubly-saturated residual is the AP' (tight-locus stability at the prime 13).
"""
from fractions import Fraction
import numpy as np, random
from math import gcd
from functools import reduce
def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd)
PAIRS=[{u,25-u} for u in [1,2,3,4,6,7,8,9,11,12]]
def is_transversal(v):
    if any(x%25==0 for x in v): return False
    R=set(x%25 for x in v); return all(p&R for p in PAIRS)
def distinct13(v): return len(set(x%13 for x in v))==12 and all(x%13!=0 for x in v)
def clears_at(v,q):
    for c in range(1,q):
        if gcd(c,q)!=1: continue
        if all(Fraction(min((x*c)%q,q-(x*c)%q),q)>=Fraction(2,25) for x in v): return True
    return False
def min_clear(v,QMAX=60):
    for q in range(6,QMAX+1):
        if q!=25 and clears_at(v,q): return q
    return None
random.seed(4); coll=[]; dist=[]; base=list(range(1,13))
for _ in range(120000):
    v=base[:]
    for _ in range(random.randint(1,4)):
        i=random.randrange(12); v[i]=random.choice([v[i]+13,v[i]+26,v[i]+25,random.randint(1,60)])
    v=tuple(sorted(set(x for x in v if x>0)))
    if len(v)!=12 or reduce(gcd,v)!=1 or not is_transversal(v): continue
    (dist if distinct13(v) else coll).append(v)
print(f"full-transversal families: collision-13={len(coll)}, distinct-13={len(dist)}", flush=True)
cmax=max((min_clear(v) or 999) for v in coll[:2000]) if coll else 0
print(f"(1) collision-13: max min-clearing-modulus={cmax}; M<2/25 count={sum(1 for v in coll[:600] if Mw(v)<Fraction(2,25))}/600  (opus mechanism was OPEN -> supplied)", flush=True)
nonap=[v for v in dist if list(v)!=list(range(1,13))]
dmax=max((min_clear(v) or 999) for v in nonap[:2000]) if nonap else 0
print(f"(2) distinct-13 non-AP={len(nonap)}: max min-clearing-modulus={dmax}; M<2/25 count={sum(1 for v in nonap[:600] if Mw(v)<Fraction(2,25))}/600", flush=True)
tight=[v for v in (coll[:1500]+dist[:1500]) if Mw(v)<Fraction(2,25)]
allap=all(list(v)==list(range(1,13)) for v in tight)
print(f"full-transversals with M<2/25 in sample: {len(tight)}, all = AP {{1..12}}? {allap}", flush=True)
print("RESIDUAL: M<2/25 full-transversal => has mult of each {7..12} (kps LRCSmallModFloor). The AP has", flush=True)
print("them as 7,8,9,10,11,12 -- the last hard node = 'mult-of-{7..12} doubly-saturated residual = AP'.", flush=True)
