# Extend the EXHAUSTIVE existence base to spread s<=24 (kps-S95), shrinking LEM-013's band.
import numpy as np
from math import gcd, floor
from functools import reduce
from itertools import combinations
import sys
def prim(E):
    return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            if L>best: best=L
    return best
def margin(E,s):
    Ea=np.array(E); worst=9.9
    for Vmax in range(s+1,floor(7*s/6)+1):
        j=np.arange(1,Vmax)
        ph=np.mod(np.outer(j,Ea),Vmax); ph.sort(axis=1)
        g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+Vmax-ph[:,-1])[:,None]],axis=1)
        b=7*int(g.max())/Vmax
        if b<worst: worst=b
    return worst
gmin=(9.9,None); n=0; fails=0
for s in range(23,25):
    for mid in combinations(range(1,s),11):
        E=[0]+list(mid)+[s]
        if not prim(E): continue
        if longest_ap(E)>7: continue
        n+=1; m=margin(E,s)
        if m<=1.0: fails+=1; print("FAIL",E,s,m); sys.stdout.flush()
        if m<gmin[0]: gmin=(m,(tuple(E),s))
    print(f"  done s={s}: cumulative {n} clusters, min margin {gmin[0]:.4f}", flush=True)
print(f"EXHAUSTIVE s=23..24 (L<=7): {n} clusters, min margin {gmin[0]:.4f} at {gmin[1]}, failures={fails}")
