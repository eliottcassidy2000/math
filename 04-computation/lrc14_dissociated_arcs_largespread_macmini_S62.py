"""
mac-mini-2026-07-09-S62 -- confirm the LARGE-SPREAD arc-count for dissociated closes:
c = #arcs/spread stays well below rho* (~0.96) as spread -> inf. opus-S168: #arcs ~ spread^0.92
(sublinear => c -> 0). Either sublinear (c->0) or linear with c<=0.51 closes c<rho*. This fits
the exponent alpha in #arcs ~ C*spread^alpha for DISSOCIATED (L<=k-6) clusters at large spread.
"""
import numpy as np
from math import gcd, log
from functools import reduce
import random
random.seed(622)

def arc_count(E, GRID):
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

print("dissociated (L<=k-6): mean #arcs and c=#arcs/spread vs spread (large-spread arc-count)\n")
for k in (13,):
    diss=k-6
    print(f"k={k} (L<={diss}):  {'spread':>7} {'mean #arcs':>11} {'mean c':>8} {'max c':>7}")
    pts=[]
    for s in [200, 500, 1200, 3000, 8000]:
        arcs=[]; cs=[]; n=0; GRID=min(40_000_000, 60*s)
        tries=0
        while n<30 and tries<600:
            tries+=1
            mid=sorted(random.sample(range(1,s),k-2)); E=[0]+mid+[s]
            if len(set(E))!=k or not prim(E) or longest_ap(E)>diss: continue
            na=arc_count(E,GRID); arcs.append(na); cs.append(na/s); n+=1
        ma=np.mean(arcs); mc=np.mean(cs); xc=max(cs)
        pts.append((s,ma)); print(f"          {s:>7} {ma:>11.1f} {mc:>8.4f} {xc:>7.4f}")
    # fit exponent alpha: #arcs ~ C spread^alpha
    xs=np.array([log(s) for s,_ in pts]); ys=np.array([log(a) for _,a in pts])
    alpha=np.polyfit(xs,ys,1)[0]
    print(f"   => #arcs ~ spread^{alpha:.3f}  ({'SUBLINEAR (c->0)' if alpha<1 else 'linear'}); "
          f"c stays << rho*~0.96 => large-spread c<rho* CLOSES\n")
