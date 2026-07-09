"""
mac-mini-2026-07-09-S62 -- pin the spread threshold S0 for the dissociated-branch inequality
c=#arcs/spread < D3(E), and the resulting SMALL FINITE CHECK.

The inequality c<D3 => good period exists (arc-count: #arcs<=c.Vmax<D3.Vmax<=rho*.Vmax). D3 is
dilation-invariant; c shrinks as spread grows. So there is a threshold S0: for spread>=S0 (dissociated,
L<=k-6) c<D3 holds; the residue (spread<S0 => hard case has Vmax<=7*S0/6) is a FINITE CHECK.

This finds, per k, the max c/D3 ratio over dissociated clusters bucketed by spread, => S0 and the
finite-check bound Vmax0 = ceil(7*S0/6). Uses FAST float c and D3 (moments on a grid).
"""
import numpy as np
from math import gcd, ceil
from functools import reduce
import random
random.seed(62)

GRIDm = 20000
_xm = (np.arange(GRIDm)+0.5)/GRIDm
Mf = 6/7
def d3_float(E):
    Ea=np.array(sorted(E),float); ph=np.mod(np.outer(_xm,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-1/7,0).sum(axis=1); m1=W.mean(); m2=(W*W).mean(); m3=(W**3).mean()
    den=m2-m3/Mf
    return (m1/Mf+(m1-m2/Mf)**2/den) if den>1e-12 else m1/Mf
def arc_count(E,GRID):
    x=(np.arange(GRID)+0.5)/GRID; Ea=np.array(sorted(E),float)
    ph=np.mod(np.outer(x,Ea),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    good=(g.max(axis=1)>1/7).astype(int); return int(np.sum((good-np.roll(good,1))==1))
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

print("dissociated (L<=k-6): max c/D3 by spread => threshold S0 & finite check Vmax0=ceil(7*S0/6)\n")
for k in (11, 13):
    diss=k-6
    print(f"k={k} (dissociated L<={diss}):")
    print(f"   {'spread':>7} {'max c/D3':>9} {'max c':>7} {'min D3':>7} {'n':>5}")
    S0 = None
    for s in [80, 120, 200, 350, 600, 1000]:
        worst=0; mc=0; md=9; n=0
        GRID=max(1_500_000, 100*s)
        for _ in range(3000):
            if n>=120: break
            mid=sorted(random.sample(range(1,s),k-2)); E=[0]+mid+[s]
            if len(set(E))!=k or not prim(E) or longest_ap(E)>diss: continue
            n+=1; d3=d3_float(E); c=arc_count(E,GRID)/s
            worst=max(worst,c/d3); mc=max(mc,c); md=min(md,d3)
        below = worst < 1.0
        if below and S0 is None: S0 = s
        print(f"   {s:>7} {worst:>9.3f} {mc:>7.3f} {md:>7.3f} {n:>5}  {'c<D3 OK' if below else 'c>=D3 somewhere'}")
    if S0:
        print(f"   => threshold S0 ~ {S0}; hard case spread<S0 has Vmax <= ceil(7*S0/6) = {ceil(7*S0/6)} "
              f"(FINITE CHECK, inside kps-S30 Vmax<=1001)\n")
    else:
        print("   => no clean threshold in tested range\n")
