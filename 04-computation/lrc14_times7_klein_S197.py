#!/usr/bin/env python3
"""
klein-2026-07-09-S197: extend LEM-012 to L=k-6 via the "x7 collapse".

For longest-AP L = k-6 (m=6 stray), cluster the L-AP at dilation j (Dirichlet) to
a super-point at p. If the config is BAD (maxgap<=V/7), then the 7 clumps (super-point
+ 6 stray) have gaps in [V/7-delta, V/7] -- FORCED to a ~V/7-grid {p + cV/7}. Then at
dilation 7j the whole grid collapses: 7*(p+cV/7) = 7p + cV = 7p mod V -- ALL clumps map
to ~7p, one super-point, huge gap. So 7j is GOOD. => j* <= 7*Q, Q the clustering bound.

Test: build hard E = AP_{k-6} u {6 stray}; at the Dirichlet-cluster j, if bad, check 7j good.
Also explore m>=7 (deeply dissociated) to see if a small multiple works.
"""
import numpy as np
from math import ceil, gcd
rng=np.random.default_rng(1970)
INV7=1/7
def maxgap(pts,V):
    p=np.sort(np.array(pts)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()/V
def is_good(E,j,V): return maxgap([(e*j)%V for e in E],V)>INV7+1e-12
def jstar(E,V,Jmax=300):
    for j in range(1,min(Jmax,V)):
        if is_good(E,j,V): return j
    return None
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best
def dirichlet_j(d,V,Q):
    best=None;bv=1.0
    for j in range(1,Q+1):
        v=abs(((j*d/V)+0.5)%1.0-0.5)
        if v<bv: bv=v;best=j
    return best

print("(1) x7 collapse for m=6 (L=k-6): at Dirichlet-cluster j, if BAD, is 7j GOOD?")
print(f"{'k':>3} {'V':>5} {'#hard,L=k-6':>12} {'#bad@clj':>9} {'#7j good when bad':>18} {'maxj*':>6}")
for k in (8,9,11,13):
    L=k-6; m=6
    for V in (91,200,400):
        built=0; badcl=0; sevengood=0; mxj=0; tries=0
        while built<120 and tries<200000:
            tries+=1
            dmin=max(1,ceil(6*V/(7*(L-1)))) if L>=2 else 1
            dmax=max(dmin,(V-1)//max(L-1,1))
            if dmin>dmax: continue
            d=int(rng.integers(dmin,dmax+1))
            ap=[(i*d)%V for i in range(L)]
            if len(set(ap))<L: continue
            pool=[x for x in range(1,V) if x not in set(ap)]
            if len(pool)<m: continue
            extra=[int(x) for x in rng.choice(pool,m,replace=False)]
            E=sorted(set(ap)|set(extra))
            if len(E)!=k or max(E)<6*V/7: continue
            if is_good(E,1,V): continue
            if longest_AP(E)!=L: continue     # exactly L
            built+=1
            Q=ceil(49*(L-1)/6) if L>=2 else 7
            jc=dirichlet_j(d,V,min(Q,V-1))
            js=jstar(E,V);
            if js: mxj=max(mxj,js)
            if not is_good(E,jc,V):
                badcl+=1
                if (7*jc)%V!=0 and is_good(E,(7*jc)%V,V): sevengood+=1
        if built>0:
            print(f"{k:>3} {V:>5} {built:>12} {badcl:>9} {sevengood:>18} {mxj:>6}")

print("\n(2) deeply dissociated m>=7 (L<=k-7): j* and does x7 (or small mult) help?")
for k in (11,13):
    for L in (k-7, k-8):
        if L<2: continue
        m=k-L; V=200
        built=0; mxj=0; tries=0
        while built<100 and tries<300000:
            tries+=1
            dmin=max(1,ceil(6*V/(7*(L-1)))); dmax=max(dmin,(V-1)//(L-1))
            if dmin>dmax: continue
            d=int(rng.integers(dmin,dmax+1))
            ap=[(i*d)%V for i in range(L)]
            if len(set(ap))<L: continue
            pool=[x for x in range(1,V) if x not in set(ap)]
            if len(pool)<m: continue
            extra=[int(x) for x in rng.choice(pool,m,replace=False)]
            E=sorted(set(ap)|set(extra))
            if len(E)!=k or max(E)<6*V/7 or is_good(E,1,V): continue
            if longest_AP(E)>L: continue
            built+=1; js=jstar(E,V)
            if js: mxj=max(mxj,js)
        print(f"  k={k} L={L} (m={m}): {built} hard sets, max j* = {mxj}")
