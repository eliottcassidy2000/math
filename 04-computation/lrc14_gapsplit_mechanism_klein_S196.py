#!/usr/bin/env python3
"""
klein-2026-07-09-S196: verify the gap-splitting MECHANISM on constructed long-AP
clusters E = AP_L u {m extra points}, m=k-L, L in {k-5,...,k}.

Theorem claim: cluster AP_L (Dirichlet j<=Q=ceil(7(L-1)/(L-k+6)), span<S=(L-k+6)/7)
=> complement gap > 1-S=(k-L+1)/7 => m extra points split into <=m+1 => max>1/7.
So the Dirichlet-cluster j (and hence j*) is good. Verify directly.
"""
import numpy as np
from math import ceil, gcd
rng=np.random.default_rng(19623)
INV7=1/7
def maxgap(pts,V):
    p=np.sort(np.array(pts)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()/V
def is_good(E,j,V): return maxgap([(e*j)%V for e in E],V)>INV7+1e-12
def jstar(E,V,Jmax):
    for j in range(1,Jmax+1):
        if is_good(E,j,V): return j
    return None
def dirichlet_j(d,V,Q):
    best=None;bv=1.0
    for j in range(1,Q+1):
        v=abs(((j*d/V)+0.5)%1.0-0.5)
        if v<bv: bv=v;best=j
    return best

print("Mechanism check: E = AP_L (step d) u {m=k-L extra}, hard (spread>=6V/7).")
print("Claim: Dirichlet-cluster j (<=Q) is a good period (>1/7 gap). Also report j*.")
print(f"{'k':>3} {'L':>3} {'m':>2} {'V':>5} {'Q':>4} {'#built':>7} {'#Dcluster good':>15} {'#j*<=Q':>8} {'maxj*':>6}")
for k in (11,12,13):
    for L in range(k, k-6, -1):   # L = k .. k-5
        m=k-L
        S=(L-k+6)/7.0
        if S<=0: continue
        Q=ceil(7*(L-1)/(L-k+6))
        for V in (91, 200, 400):
            built=0; dgood=0; jok=0; mxj=0; tries=0
            while built<80 and tries<40000:
                tries+=1
                # AP_L with step d chosen so spread>=6V/7: put AP spanning most of [0,V)
                # d ~ 6V/(7(L-1)) up to V/(L-1); random d in a range
                dmin=max(1, ceil(6*V/(7*(L-1))))
                dmax=max(dmin, (V-1)//(L-1))
                if dmin>dmax: continue
                d=rng.integers(dmin,dmax+1)
                ap=[(i*d) for i in range(L)]
                if max(ap)>=V: continue
                # m extra points random, distinct, not in AP; keep total spread>=6V/7
                pool=[x for x in range(1,V) if x not in set(ap)]
                if len(pool)<m: continue
                extra=list(rng.choice(pool,m,replace=False)) if m>0 else []
                E=sorted(set(ap)|set(int(x) for x in extra))
                if len(E)!=k: continue
                if max(E)<6*V/7: continue
                if is_good(E,1,V): continue   # want hard
                built+=1
                jc=dirichlet_j(d,V,Q)
                if is_good(E,jc,V): dgood+=1
                js=jstar(E,V,min(Q,V-1))
                if js is not None: jok+=1; mxj=max(mxj,js)
            if built>0:
                print(f"{k:>3} {L:>3} {m:>2} {V:>5} {Q:>4} {built:>7} {dgood:>15} {jok:>8} {mxj:>6}")
print("\n=> if #Dcluster good == #built (mechanism) and #j*<=Q == #built for ALL L>=k-5,")
print("   the gap-split argument is CONFIRMED. (m=0 row = mac-mini's exact-AP case.)")
