#!/usr/bin/env python3
"""
klein-2026-07-09-S199: is the LEM-012 / route-(c) INTERLOCK airtight?

LEM-012 covers longest-AP L >= k-6 (Dirichlet gap-split + x7). Route (c) covers
L <= k-7 IF c(L):=max #arcs/spread < rho_min(L) >= D3_inf^{(L)}. The 2-block warns
that #arcs can be LARGE (c=2.45) for coherent structure -- but a k=13 2-block has
L>=7 (LEM-012's domain). CHECK: for L <= 6 (route c's domain, k=13), is the WORST
c(L) still < D3_inf^{(L)}?  Adversarially construct high-arc-count sets with
longest-AP EXACTLY L (L-AP + clustered/structured extras that keep L fixed).
"""
import numpy as np
from math import gcd
rng=np.random.default_rng(1999)
INV7=1/7
D3L={2:0.86,3:0.85,4:0.84,5:0.81,6:0.76,7:0.68}  # opus-S158 D3_inf^{(L)}
def arcs(E,V,Nx=200000):
    E=np.array(sorted(set(E)))-min(E); x=(np.arange(Nx)+0.5)/Nx*V
    ph=np.sort(np.mod(np.outer(x,E),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    gi=(g.max(axis=1)>INV7+1e-12).astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((ed==1).sum());
    if gi.all():nc=1
    return nc
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best

k=13
print(f"k={k}: worst c(L)=max #arcs/spread over longest-AP=L sets, vs D3_inf^(L) (route c needs c<D3)")
print(f"{'L':>3} {'#tested':>8} {'max c(L)':>9} {'D3_inf^L':>9} {'c<D3?':>7} {'margin':>7} {'worst shape (rel)':>22}")
for L in (2,3,4,5,6):
    maxc=0; worst=None; ntest=0
    for _ in range(6000):
        # construct: an L-AP (step d) + (k-L) extras that DON'T extend to an (L+1)-AP,
        # biased toward CLUSTERED/structured extras (adversarial for #arcs)
        d=int(rng.integers(1,12))
        ap=[i*d for i in range(L)]
        mode=rng.integers(0,3)
        base=max(ap)
        extras=[]
        # clustered extras (2-3 tight blocks) to maximize coherence
        for _ in range(k-L):
            if mode==0: extras.append(int(rng.integers(1, base+120)))
            elif mode==1: extras.append(base+10+int(rng.integers(0,8)))  # one far cluster
            else: extras.append(int(rng.integers(0, base+40)))
        E=sorted(set(ap)|set(extras))
        if len(E)!=k: continue
        if longest_AP(E)!=L: continue   # exactly L
        sp=max(E)-min(E)
        if sp<7: continue
        V=int(sp/0.85)+2
        c=arcs(E,V)/sp
        ntest+=1
        if c>maxc: maxc=c; worst=E
    d3=D3L[L]
    ws=str([e-min(worst) for e in worst][:7])+"..." if worst else "-"
    print(f"{L:>3} {ntest:>8} {maxc:>9.3f} {d3:>9.3f} {str(maxc<d3):>7} {d3-maxc:>7.3f} {ws:>22}")

print("\n=> if max c(L) < D3_inf^(L) for ALL L in 2..6, the interlock is AIRTIGHT:")
print("   [L>=7: LEM-012 elementary] + [L<=6: route (c) c(L)<D3<=rho*] cover all L, a-priori.")
print("   If some c(L) >= D3, that L is a GAP -- needs LEM-012 extended or a tighter arc bound.")
