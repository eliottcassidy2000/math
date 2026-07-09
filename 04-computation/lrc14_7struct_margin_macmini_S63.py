"""
mac-mini-2026-07-09-S63 -- the 7-structured existence margin (mu - c).

MISTAKE-128: 7-structured co-offsets (many elements =0 mod 7) spike #arcs (c high) but mu is also
high; moments can't certify (B_6=0.77<c=0.88). The good period still EXISTS iff c<mu. This finds
the MIN margin (mu - c) over PRIMITIVE 7-structured dissociated (L<=k-6) clusters, adversarially,
vs spread -- does it stay positive and grow? If so, dissociated closes via direct existence
[exhaustive small spread] + [margin-floor large spread].
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(63)

def cmu(E, GRID):
    x=(np.arange(GRID)+0.5)/GRID; Ea=np.array(sorted(E),float)
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    mg=g.max(axis=1); good=(mg>1/7)
    arcs=int(np.sum((good.astype(int)-np.roll(good.astype(int),1))==1))
    return arcs/(max(E)-min(E)), float(good.mean())
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
def make_7struct(k, s, n7):
    # n7 elements =0 mod 7 (a step-7 family) + (k-n7) non-mult-of-7, spread ~ s, primitive
    sevens=sorted(random.sample(range(0, s+1, 7), min(n7, s//7+1)))
    if len(sevens) < n7: return None
    sevens=sevens[:n7]
    rest_pool=[x for x in range(1, s) if x%7!=0]
    if len(rest_pool) < k-n7: return None
    rest=random.sample(rest_pool, k-n7)
    E=sorted(set(sevens+[0,s]+rest))
    if len(E)!=k: return None
    return E

print("min existence margin (mu - c) over PRIMITIVE 7-structured dissociated (L<=k-6), by spread:\n")
for k in (13,):
    diss=k-6
    print(f"k={k} (L<={diss}):  {'spread':>7} {'min(mu-c)':>10} {'at c':>7} {'at mu':>7} {'#tested':>8} {'n7 worst':>9}")
    for s in [50, 84, 150, 300, 700]:
        GRID=max(2_000_000, 100*s); best=(9,0,0,0); n=0
        for n7 in range(2, min(9, s//7+1)):
            for _ in range(400):
                E=make_7struct(k, s, n7)
                if E is None or not prim(E) or longest_ap(E)>diss: continue
                if max(E)-min(E) < 6*s//7: continue    # hard (spread-dense)
                n+=1; c,mu=cmu(E,GRID); m=mu-c
                if m<best[0]: best=(m,c,mu,n7)
        print(f"          {s:>7} {best[0]:>10.4f} {best[1]:>7.3f} {best[2]:>7.3f} {n:>8} {best[3]:>9}")
    print("   => if min(mu-c) stays positive and grows with spread, 7-structured closes via existence.")
