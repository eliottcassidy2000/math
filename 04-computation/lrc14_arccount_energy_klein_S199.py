#!/usr/bin/env python3
"""
klein-2026-07-09-S199: the a-priori ARC-COUNT bound via ADDITIVE ENERGY (Freiman theme).

Route (c) closes the Sidon branch IF #arcs(Good_E) < rho*·Vmax, i.e. c(L):=#arcs/spread
< rho*. The a-priori Davenport-Schinzel bound #arcs=O(k^3 spread) is too weak (c=O(k^3)).
FREIMAN/HADWIGER theme (kps-S93): the AP/Sidon dichotomy = additive density=>structure.
So #arcs should be governed by the ADDITIVE ENERGY E(A) (or longest-AP), NOT k^3:
low energy (Sidon) => few coherent resonances => few arcs => small c.

Test: for sets ranging Sidon -> AP (k=13), compute #arcs/spread, additive energy E(A),
longest-AP L, exact-relation count Z, and find the a-priori law #arcs <= f(structure).
"""
import numpy as np
from math import gcd
from collections import Counter
from itertools import product
rng=np.random.default_rng(199)
INV7=1/7; M=6/7

def arcs(E,V,Nx=200000):
    x=(np.arange(Nx)+0.5)/Nx*V; e=np.array(E)
    ph=np.sort(np.mod(np.outer(x,e),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    gi=(g.max(axis=1)>INV7+1e-12).astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((ed==1).sum());
    if gi.all():nc=1
    return nc
def add_energy(E):
    # E(A)=#{(a,b,c,d): a+b=c+d}, a,b,c,d in E
    s=Counter(); Ea=list(E)
    for a in Ea:
        for b in Ea: s[a+b]+=1
    return sum(v*v for v in s.values())
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best
def Zcount(E,H=3):
    # #{n: support<=3, |n_i|<=H, n!=0, sum n_i e_i = 0}  (exact additive relations)
    z=0; E=list(E); k=len(E)
    from itertools import combinations
    for nz in (2,3):
        for combo in combinations(range(k),nz):
            for vals in product(range(-H,H+1),repeat=nz):
                if any(v==0 for v in vals): continue
                if sum(E[i]*v for i,v in zip(combo,vals))==0: z+=1
    return z
def good_j1(E,V):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()>V/7+1e-9

k=13
print(f"k={k}: #arcs/spread vs additive energy, longest-AP, Z (Freiman/energy law)")
print(f"{'kind':>10} {'spread':>7} {'#arcs':>6} {'c=#arcs/sp':>11} {'E(A)':>6} {'longAP':>7} {'Z':>5} {'D3_inf^L(ref)':>13}")
D3L={2:0.86,3:0.85,4:0.84,5:0.81,6:0.76,7:0.68,8:0.60,9:0.52,10:0.465}  # opus-S158
def report(E,V,kind):
    E=sorted(set(E)); sp=max(E)-min(E)
    nc=arcs([e-min(E) for e in E],V); ea=add_energy(E); L=longest_AP(E); z=Zcount([e-min(E) for e in E])
    d3=D3L.get(min(L,10),0.465)
    print(f"{kind:>10} {sp:>7} {nc:>6} {nc/sp:>11.3f} {ea:>6} {L:>7} {z:>5} {d3:>13.3f}")
    return nc/sp, ea, L
# exact AP
V=200; d=13
report([i*d for i in range(k)], V, "exactAP")
# perfect-difference / Sidon set (k=13): use a known Sidon set
sidon=[0,1,3,7,12,20,30,44,65,80,96,122,147]  # roughly Sidon (distinct diffs), fit in ~150
# ensure distinct diffs; if not, perturb -- just use a random-greedy Sidon
def greedy_sidon(k,hi):
    S=[0]; diffs=set()
    x=1
    while len(S)<k and x<hi:
        ok=all(abs(x-s) not in diffs for s in S)
        if ok:
            for s in S: diffs.add(abs(x-s))
            S.append(x)
        x+=1
    return S
sd=greedy_sidon(k,300); Vs=int(max(sd)/0.85)
report(sd, Vs, "Sidon")
# intermediate: perturbed AP / 2-block / random
for tag,E,V in [("2block",[0,1,2,3,4,5,100,101,102,103,104,105,106],130),
                ("randlo",sorted(rng.choice(range(1,180),k-1,replace=False).tolist()+[0]),210),
                ("perturbAP",[i*13+int(rng.integers(-3,4)) for i in range(k)],200)]:
    report(E,V,tag)

print("\n=> if c=#arcs/spread TRACKS additive energy/longest-AP (low for Sidon, high for AP)")
print("   and c < D3_inf^L (which is HIGH for low L), then route (c) closes a-priori via")
print("   the additive-energy law -- the Freiman 'density=>structure' shadow, made quantitative.")
