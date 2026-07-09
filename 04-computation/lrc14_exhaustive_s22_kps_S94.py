# Push the EXHAUSTIVE existence base to spread s<=22 (kps-S94), pure-dissociated L<=6 AND L<=7.
import numpy as np
from math import gcd, floor
from functools import reduce
from itertools import combinations
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
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
def cluster_margin(E,s):
    Ea=np.array(sorted(E)); worst=9.9
    for Vmax in range(s+1,floor(7*s/6)+1):
        j=np.arange(1,Vmax)
        ph=np.mod(np.outer(j,Ea),Vmax); ph.sort(axis=1)
        g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+Vmax-ph[:,-1])[:,None]],axis=1)
        b=7*g.max(axis=1).max()/Vmax
        if b<worst: worst=b
    return worst
gmin6=(9.9,None); gmin7=(9.9,None); n6=0; n7=0; fails=0
for s in range(13,23):
    for mid in combinations(range(1,s),11):
        E=[0]+list(mid)+[s]
        if not prim(E): continue
        L=longest_ap(E)
        if L>7: continue
        m=cluster_margin(E,s)
        if m<=1.0: fails+=1; print("FAIL",E,s,m)
        n7+=1
        if m<gmin7[0]: gmin7=(m,(tuple(E),s,L))
        if L<=6:
            n6+=1
            if m<gmin6[0]: gmin6=(m,(tuple(E),s,L))
print(f"EXHAUSTIVE s<=22: L<=7: {n7} clusters, min margin {gmin7[0]:.4f} at {gmin7[1]}")
print(f"                  L<=6: {n6} clusters, min margin {gmin6[0]:.4f} at {gmin6[1]}")
print(f"existence failures (margin<=1): {fails}")
