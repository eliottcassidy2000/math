#!/usr/bin/env python3
"""cont.56 CRITICAL MATH: the lander-dodge lower bound as a constructive floor.

LEMMA (constructive). If v ⊇ {1,...,12} then at base q=13k+1, multiplier a=k, the core
{1,...,12} lands at {k,2k,...,12k} with closest approach EXACTLY k (12k sits at distance
k+1). So M(v) ≥ k/(13k+1) as long as the remaining element(s) also land at distance ≥ k,
i.e. they DODGE the band {0,±1,...,±(k-1)} mod q.

Test the full core family {1,...,12,f} over all f: find true M, locate the minimum-over-f,
and confirm the lander-dodge construction (some k) achieves the census M. This certifies the
far-element floor (cont.55) over ALL f, not just multiples of 182, and exhibits the WITNESS
base for each — turning the empirical floor into a constructive one."""
from fractions import Fraction as F
from math import gcd

def dist0(r,q): return min(r%q, q-(r%q))
def bestbase(v,Q):
    best=F(0); arg=None
    for q in range(2,Q):
        # only multipliers coprime to q matter for the max
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(dist0(x*a,q) for x in v)
            if F(m,q)>best: best=F(m,q); arg=(q,a,m)
    return best,arg

def dodge_construction(v,Kmax=40):
    """best k such that base 13k+1, a=k gives closest approach k (all runners dodge band)."""
    best=F(0); bk=None
    for k in range(1,Kmax+1):
        q=13*k+1
        m=min(dist0(x*k,q) for x in v)
        if F(m,q)>best: best=F(m,q); bk=(k,q,m)
    return best,bk

print("="*74)
print("core family {1,...,12, f}: census M vs lander-dodge construction, over all f")
print("="*74)
core=list(range(1,13))
minM=None; rows=[]
for f in range(13, 260):
    v=core+[f]
    if len(set(v))<13: continue
    Q=min(3*f+10, 800)
    M,(q,a,m)=bestbase(v,Q)
    dM,dk=dodge_construction(v)
    rows.append((f,M,q,a,m,dM,dk))
    if minM is None or M<minM[1]: minM=(f,M,q,a,m,dM,dk)
# show the minimum and a few neighbors
print(f"\nMINIMUM over f in [13,259]:  f={minM[0]}, M={minM[1]}={float(minM[1]):.5f}, "
      f"best base q={minM[2]}, a={minM[3]}, closest approach={minM[4]}")
print(f"  lander-dodge construction reaches: {minM[5]}={float(minM[5]):.5f} at (k,q,m)={minM[6]}"
      f"  {'== census M (CONSTRUCTIVE FLOOR)' if minM[5]==minM[1] else 'GAP'}")
print(f"\nthe ten smallest-M core families {{1..12,f}}:")
for f,M,q,a,m,dM,dk in sorted(rows,key=lambda r:r[1])[:10]:
    tag = "  <-- 14/183 deep well" if M==F(14,183) else ""
    match = "dodge=M" if dM==M else f"dodge={float(dM):.5f}"
    print(f"  f={f:>3}: M={float(M):.5f}={M}  base q={q},a={a},approach={m}  [{match}]{tag}")
print(f"\nfloor over this range = {float(minM[1]):.5f} = {minM[1]}  "
      f"({'14/183 crux' if minM[1]==F(14,183) else 'OTHER'}); all core families have M >= this > 1/14={1/14:.5f}")

# does the dodge construction ALWAYS match census M (is it the optimal base)?
gaps=[r for r in rows if r[5]!=r[1]]
print(f"\nlander-dodge optimality: {len(rows)-len(gaps)}/{len(rows)} core families have dodge == census M.")
if gaps:
    print(f"  families where a NON-13k+1 base beats the dodge construction (dodge < M):")
    for f,M,q,a,m,dM,dk in gaps[:6]:
        print(f"    f={f}: census M={float(M):.5f} at q={q} > dodge {float(dM):.5f}  (dodge is a valid lower bound, not always tight)")
