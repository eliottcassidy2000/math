#!/usr/bin/env python3
"""
klein-2026-07-09-S199: investigate the INTERLOCK GAP at longest-AP L=5,6 (k=13).

Route (c) FAILS there (c=#arcs/spread=1.2-1.7 > rho*). These are 2-BLOCK-like sets
(a 6-AP + a far cluster): high arc-count but L<7 so LEM-012 doesn't apply. Questions:
 (1) Do good periods EXIST (j* small)? => proof-route gap, not a real counterexample.
 (2) Does route (c) work with the ACTUAL rho*=mu (not the D3 lower bound)?
 (3) MECHANISM: are these handled by clustering the 6-AP (a LEM-012-like small-j)?
     A 2-block has a long AP in EACH block; clustering the LONGER block's AP should work.
"""
import numpy as np
from math import gcd, ceil
rng=np.random.default_rng(19999)
INV7=1/7
def maxgap1(E,j,V):
    p=np.sort(np.mod(np.array(E)*j,V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()/V
def jstar(E,V,Jmax=200):
    for j in range(1,min(Jmax,V)):
        if maxgap1(E,j,V)>INV7+1e-12: return j
    return None
def mu(E,V,Nx=200000):
    E=np.array(E); x=(np.arange(Nx)+0.5)/Nx*V
    ph=np.sort(np.mod(np.outer(x,E),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return (g.max(axis=1)>INV7+1e-12).mean()
def arcs(E,V,Nx=200000):
    E=np.array(E); x=(np.arange(Nx)+0.5)/Nx*V
    ph=np.sort(np.mod(np.outer(x,E),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    gi=(g.max(axis=1)>INV7+1e-12).astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((ed==1).sum())
    if gi.all():nc=1
    return nc
def longest_AP_with_d(E):
    Es=sorted(set(E)); Eset=set(Es); best=(1,None)
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            if L>best[0]: best=(L,d)
    return best

k=13
print("Worst high-c L=5,6 sets (2-block-like): do good periods exist? route(c) w/ actual mu?")
print(f"{'L':>3} {'V':>5} {'#arcs':>6} {'c':>6} {'mu(rho*)':>8} {'c<mu?':>6} {'j*':>4} {'clusterAP j':>11}")
found={5:[],6:[]}
for _ in range(30000):
    L=rng.choice([5,6]); d=int(rng.integers(1,10))
    ap=[i*d for i in range(L)]; base=max(ap)
    extras=[base+8+int(rng.integers(0,10)) for _ in range(k-L)]  # far cluster (2-block)
    E=sorted(set(ap)|set(extras))
    if len(E)!=k: continue
    Lact,dd=longest_AP_with_d(E)
    if Lact!=L: continue
    E=[e-min(E) for e in E]; sp=max(E)
    if sp<7: continue
    V=int(sp/0.85)+2
    c=arcs(E,V)/sp
    if c<0.9: continue    # want the high-c (route-c-failing) ones
    m=mu(E,V); js=jstar(E,V)
    # cluster the longest AP (LEM-012-style Dirichlet) -- smallest good j via that
    Q=ceil(49*(L-1)/3)
    jcl=None
    for j in range(1,Q+1):
        if maxgap1(E,j,V)>INV7+1e-12: jcl=j; break
    found[L].append((V,arcs(E,V),c,m,js,jcl))
    if all(len(found[l])>=4 for l in (5,6)): break
for L in (5,6):
    for (V,nc,c,m,js,jcl) in found[L][:4]:
        print(f"{L:>3} {V:>5} {nc:>6} {c:>6.2f} {m:>8.3f} {str(c<m):>6} {str(js):>4} {str(jcl):>11}")
print("\n=> if j* is small (<=~14) for ALL these, good periods EXIST -- the gap is a PROOF-ROUTE")
print("   gap, resolved by clustering the block's AP (LEM-012 extended to 2-blocks: cluster the")
print("   LONGER block, a >=ceil(k/2)-AP, which IS >= k-6 for k<=13 within that block).")
