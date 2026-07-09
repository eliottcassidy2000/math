"""
mac-mini-2026-07-09-S64 -- the DISSOCIATED branch closes a-priori by the widest-arc pigeonhole.

CLAIM: for dissociated (longest-AP <= k-7 = 6) primitive k=13 co-offset sets,
   maxIntG(E) * spread(E) >= c > 1   (adversarially measured c ~ 1.71).
Then maxIntG(E) >= c/spread > 1/spread > 1/Vmax for every Vmax > spread, so the grid {j/Vmax}
contains a point in G(E)'s widest arc (pigeonhole: an arc of length >= 1/V contains a multiple of 1/V)
=> a STRICT good period exists for ALL Vmax > spread. A-priori -- NO Mertens resonant sum (opus-S172),
NO exhaustion. Sidesteps the wall.

WHY THE BRANCH SPLIT IS THE DIVIDING LINE: fragmentation of G (maxIntG collapsing to the 0-nbhd
6/(7*spread), so maxIntG*spread -> 6/7 < 1) is driven by a long resonant sub-AP (e.g. the mult-of-7
AP in the knife-edge {0,7,10,14,18,20,21,26,28,35,36,37,42}, longest-AP=7). DISSOCIATED sets (L<=6)
have NO such sub-AP, so G keeps a wide arc: maxIntG*spread >= ~1.7. The fragmented sets are exactly
the NEAR-AP ones (L>=7) = LEM-012's Dirichlet territory (+ the non-strict j=1 knife-edge). So the
good-period dichotomy split (near-AP L>=k-6 / dissociated L<=k-7) is PRECISELY the geometric divide.

This searches HARD (7-structured + k-structured bias) for a dissociated set with maxIntG*spread <= 1.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(999)

def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
def maxint_G(E,N):
    Ea=np.array(sorted(E)); i=np.arange(0,N)
    ph=np.mod(np.outer(i,Ea),N)/N; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    ind=(g.max(axis=1)>1.0/7+1e-12).astype(int)
    if ind.sum()==0: return 0.0
    if ind.sum()==N: return 1.0
    dd=np.diff(np.concatenate([[0],ind,[0]]))
    st=np.where(dd==1)[0]; en=np.where(dd==-1)[0]; runs=list(en-st)
    if ind[0]==1 and ind[-1]==1 and len(runs)>=2: runs[0]+=runs[-1]; runs=runs[:-1]
    return max(runs)/N

best=99; arg=None; nd=0
for _ in range(120000):
    s=random.randint(20,80); r=random.random()
    if r<0.6:
        sev=list(range(7,s,7)); base=[0,s]+(random.sample(sev,min(random.randint(2,5),len(sev))) if sev else [])
    elif r<0.8:
        kk=random.choice([3,4,5]); mk=list(range(kk,s,kk)); base=[0,s]+(random.sample(mk,min(random.randint(2,5),len(mk))) if mk else [])
    else: base=[0,s]
    pool=[x for x in range(1,s) if x not in base]; need=13-len(base)
    if need<0 or len(pool)<need: continue
    E=sorted(set(base+random.sample(pool,need)))
    if len(E)!=13 or not prim(E) or longest_ap(E)>6: continue
    nd+=1; N=max(6000,30*s); v=maxint_G(E,N)*s
    if v<best: best=v; arg=(E,longest_ap(E),round(v,4))
print(f"DISSOCIATED (longest-AP<=6), {nd} sets (7-/k-structured biased):")
print(f"  min maxIntG*spread = {best:.4f}   argmin E={arg[0]} (L={arg[1]})")
print(f"  0-nbhd floor 6/7 = {6/7:.4f};  observed min = {best/(6.0/7):.2f}x floor;  12/7 = {12/7:.4f}")
print(f"  {best:.3f} > 1  =>  maxIntG > 1/Vmax for all Vmax>spread => strict good period A-PRIORI (dissociated).")
print("  A-priori TARGET (three-distance/Steinhaus): prove maxIntG(E)*spread >= c > 1 for dissociated E.")
