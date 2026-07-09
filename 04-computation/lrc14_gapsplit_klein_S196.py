#!/usr/bin/env python3
"""
klein-2026-07-09-S196: THE GAP-SPLITTING ARGUMENT for general j*=O(k) (elementary).

If E contains an AP of length L (common difference d), then at a dilation j that
clusters that AP into a circular span < S := (L-k+6)/7 (Dirichlet: exists j <=
ceil((L-1)/S) = ceil(7(L-1)/(L-k+6)), UNCONDITIONAL in d,Vmax), the AP occupies a
<S arc so its complement is one gap of length > 1-S = (k-L+1)/7. The remaining
m=k-L points split that gap into <= m+1 pieces => max piece > (1-S)/(m+1) = 1/7.
=> a good period at j <= ceil(7(L-1)/(L-k+6)) = O(k), whenever L >= k-5 (so S>0).
Combined with kps-S91 (longest-AP <= k-3 => j*<=3), ALL L are covered.

VERIFY: for hard E, find longest-AP L; if L>=k-5, check a good period exists at
j <= ceil(7(L-1)/(L-k+6)); report the Dirichlet-clustering j and the actual j*.
"""
import numpy as np
from math import ceil, gcd
rng=np.random.default_rng(1962)
INV7=1/7

def maxgap(pts,V):
    p=np.sort(np.array(pts)%V); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()/V
def is_good(E,j,V): return maxgap([(e*j)%V for e in E],V)>INV7+1e-12
def jstar(E,V,Jmax=None):
    Jmax=Jmax or V
    for j in range(1,Jmax+1):
        if is_good(E,j,V): return j
    return None

def longest_AP(E):
    """longest arithmetic progression (integer common difference) that is a subset of E."""
    Es=sorted(set(E)); n=len(Es); Eset=set(Es); best=(1,None,None)
    if n<2: return (n,None,None)
    for i in range(n):
        for jj in range(i+1,n):
            a,d=Es[i],Es[jj]-Es[i]
            L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            if L>best[0]: best=(L,a,d)
    return best  # (length, start a, diff d)

def dirichlet_cluster_j(d,V,Q):
    """Dirichlet: smallest j in 1..Q with ||j d / V|| < 1/Q  (exists by pigeonhole)."""
    best=None; bestval=1.0
    for j in range(1,Q+1):
        val=abs(((j*d/V)+0.5)%1.0-0.5)
        if val<bestval: bestval=val; best=j
    return best,bestval

def sample_hard(k,V,n,tries=60000):
    out=[];t=0
    while len(out)<n and t<tries:
        t+=1
        rest=rng.choice(np.arange(1,V),k-1,replace=False)
        E=tuple(sorted([0]+[int(x) for x in rest]))
        if max(E)<6*V/7: continue
        if is_good(E,1,V): continue
        out.append(E)
    return out

print("Verify the gap-splitting bound: longest-AP L>=k-5 => j* <= ceil(7(L-1)/(L-k+6)).")
print(f"{'k':>3} {'V':>5} {'#hard':>6} {'#(L>=k-5)':>10} {'bound holds':>12} {'max(j*/bound)':>13} {'#(L<k-5)':>9}")
for k in (8,10,11,12,13):
    for V in (2*k+3, 91, 200):
        hs=sample_hard(k,V,300)
        if not hs: continue
        nge=0; nlt=0; ok=0; ratios=[]
        for E in hs:
            L,a,d=longest_AP(E)
            m=k-L
            if L>=k-5:
                nge+=1
                S=(L-k+6)/7.0
                Q=ceil(7*(L-1)/(L-k+6))
                # does a good period exist at j<=Q? (the theorem claim)
                js=jstar(E,V,Jmax=min(Q,V-1))
                if js is not None:
                    ok+=1; ratios.append(js/Q)
            else:
                nlt+=1
        mr=max(ratios) if ratios else 0
        print(f"{k:>3} {V:>5} {len(hs):>6} {nge:>10} {str(ok==nge):>12} {mr:>13.3f} {nlt:>9}")

print("\nSpot-check: the Dirichlet-clustering j actually leaves a >1/7 gap (the mechanism):")
for k in (11,13):
    V=91
    hs=sample_hard(k,V,400)
    checked=0; mech_ok=0
    for E in hs:
        L,a,d=longest_AP(E)
        if L<k-5: continue
        Q=ceil(7*(L-1)/(L-k+6))
        jc,val=dirichlet_cluster_j(d,V,Q)
        checked+=1
        if is_good(E,jc,V): mech_ok+=1
    print(f"  k={k}: over {checked} clusters with L>=k-5, the Dirichlet-cluster j leaves a >1/7 gap: {mech_ok}/{checked}")
print("\n=> if 'bound holds'=True everywhere and mechanism works, general j*=O(k) is PROVED")
print("   (gap-split for L>=k-5) + (kps-S91 dissociated for L<=k-3); ranges OVERLAP => all L covered.")
