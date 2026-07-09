#!/usr/bin/env python3
"""
klein-2026-07-09-S202: is the MULTIVARIATE grid-residual bound tight enough to close the HARD boundary?
The TV bound R_grid_abs <= TV(Phi')/(12Vmax^2) is ~30-50x loose (absolute, ignores that b0(n),c(sigma)
~0.13 << 1/pi). Compute R_grid_abs DIRECTLY from LEM-011 exact |What(n)|, summed over GRID resonances
(n.e = mVmax, m!=0), enumerated by support r and height H, + geometric tail. If this multivariate sum
< E[W] at the hard boundary (spread ~ Vmax), the transport is TIGHT (closes it a-priori). If it too
overshoots, the barely-covers wall is real and the transport is partial (large-ruler only).

LEM-011: |What(n)| = (6/7)^{k-1-r} prod_{n_i!=0}|b0(n_i)| * min(6/7,|c(sigma)|),
 b0(m)=|sin(pi m/7)|/(pi|m|), c(s)=|sin(pi s/7)|/(pi|s|); What(n)=0 if 7|n_i any i, or 7|sigma & sigma!=0.
"""
import numpy as np
from math import sin, pi, gcd
from functools import reduce
from itertools import combinations, product
INV7=1/7
def b0(m):
    m=abs(int(m));  return 0.0 if m%7==0 else abs(sin(pi*m/7))/(pi*m)
def cabs(s):
    s=int(s)
    if s==0: return 6/7
    return 0.0 if s%7==0 else abs(sin(pi*s/7))/(pi*abs(s))
def what_abs(nz, k):  # nz = list of nonzero coeffs; returns |What(n)|
    r=len(nz); sig=sum(nz)
    p=(6/7)**(k-1-r)
    for v in nz: p*=b0(v)
    if p==0: return 0.0
    return p*min(6/7, cabs(sig))
def EW_cont(E,N=1<<17):
    E=np.array(sorted(set(E)),float); y=np.arange(N)/N
    ph=np.mod(np.outer(y,E),1.0); ph.sort(axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return np.maximum(g-INV7,0).sum(axis=1).mean()
def grid_resonant_abs(E, Vmax, Rmax=4, Hmax=10):
    """Sum 2|What(n)| over n with support<=Rmax, height<=Hmax, Vmax | n.e, n.e != 0."""
    E=[int(x) for x in sorted(set(E))]; k=len(E); idx=list(range(1,k))  # e_0=0 pinned, vars n_1..n_{k-1}
    ev=[E[i] for i in idx]
    tot=0.0; count=0
    for r in range(1,Rmax+1):
        for pos in combinations(range(len(idx)), r):
            es=[ev[p] for p in pos]
            # nonzero coeff vectors with sum|.|<=Hmax, each in [-Hmax,Hmax]\{0}
            for vals in product([v for v in range(-Hmax,Hmax+1) if v!=0], repeat=r):
                if sum(abs(v) for v in vals)>Hmax: continue
                ne=sum(v*e for v,e in zip(vals,es))
                if ne==0 or ne%Vmax!=0: continue
                w=what_abs(list(vals), k)
                if w>0: tot+=2*w; count+=1
    return tot,count

clusters={
 'tightAP{0..12}@V=15(Q+1,hasGP)': (list(range(13)),15),
 'tightAP{0..12}@V=13(noGP)':      (list(range(13)),13),
 'near-AP d=3 @boundary':          ([3*i for i in range(13)], 43),
 '7-struct(M128)@V=91':            ([0,7,14,21,26,29,37,44,51,58,67,75,82],91),
 'dissoc@boundary':                ([0,1,3,7,12,20,30,44,65,80,96,112,129],131),
}
print("MULTIVARIATE grid-residual sum (exact LEM-011 |What|) vs E[W]. tight if < E[W].\n")
print(f"{'cluster':>34} {'Vmax':>5} {'E[W]':>7} {'R_grid_mv':>10} {'#res':>6} {'ratio/EW':>8} {'<E[W]?':>6}")
for nm,(E,V) in clusters.items():
    E=[e-min(E) for e in sorted(set(E))]
    ew=EW_cont(E); rg,cnt=grid_resonant_abs(E,V,Rmax=4,Hmax=9)
    print(f"{nm:>34} {V:>5} {ew:>7.4f} {rg:>10.5f} {cnt:>6} {rg/ew:>8.3f} {str(rg<ew):>6}")
print("\nIf R_grid_mv < E[W] at the HARD boundary (esp. 7-struct/dissoc), the multivariate transport is")
print("TIGHT enough to close a-priori. tightAP@V=13 (no GP) should have R_grid_mv >= E[W]-(6/7)/V region.")
