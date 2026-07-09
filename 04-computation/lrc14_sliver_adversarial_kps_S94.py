# Adversarial MINIMUM good-period-existence margin over dissociated k=13 clusters (kps-S94).
# margin(E) := min over critical Vmax in [s+1, floor(7s/6)] of  max_j 7*maxgap(E,j,Vmax)/Vmax.
# Good period exists iff margin>1. Hill-climb to MINIMIZE margin (adversary). If min stays >>1,
# the small-spread dissociated regime is robustly safe (sliver = certificate artifact, no real gap).
# Plus: EXHAUSTIVE over small spreads s<=20 as a rigorous base.
import numpy as np
from math import gcd, ceil, floor
from functools import reduce
from itertools import combinations
import random
random.seed(94099)
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
    # min over critical Vmax of (max_j 7*maxgap/Vmax). Existence iff >1.
    Ea=np.array(sorted(E)); worst=9.9
    lo=s+1; hi=floor(7*s/6)
    for Vmax in range(lo,hi+1):
        j=np.arange(1,Vmax)
        ph=np.mod(np.outer(j,Ea),Vmax); ph.sort(axis=1)
        g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+Vmax-ph[:,-1])[:,None]],axis=1)
        best=7*g.max(axis=1).max()/Vmax          # best period at this Vmax
        if best<worst: worst=best                # adversary picks worst Vmax
    return worst
# ---- (A) EXHAUSTIVE small spread ----
print("(A) EXHAUSTIVE dissociated k=13, spread s<=20 (rigorous base):")
gmin_exh=(9.9,None); nexh=0
for s in range(13,21):
    for mid in combinations(range(1,s),11):
        E=[0]+list(mid)+[s]
        if not prim(E) or longest_ap(E)>7: continue
        nexh+=1; m=cluster_margin(E,s)
        if m<gmin_exh[0]: gmin_exh=(m,(tuple(E),s))
        if m<=1.0: print("   *** EXISTENCE FAILURE",E,s,m)
print(f"    exhausted {nexh} dissociated clusters s<=20; MIN existence margin = {gmin_exh[0]:.4f} at {gmin_exh[1]}")
print(f"    (margin>1 everywhere => good period exists for ALL of them)\n")
# ---- (B) ADVERSARIAL hill-climb min-margin across larger spreads ----
print("(B) ADVERSARIAL min-margin (hill-climb, minimize existence margin) per spread:")
def rand_diss(s):
    for _ in range(200):
        mid=sorted(random.sample(range(1,s),11)); E=[0]+mid+[s]
        if len(set(E))==13 and prim(E) and longest_ap(E)<=7: return E
    return None
for s in (50,70,90,120,160,200):
    best=(9.9,None)
    for _restart in range(60):
        E=rand_diss(s)
        if E is None: continue
        cur=cluster_margin(E,s)
        improved=True
        while improved:
            improved=False
            for idx in range(1,12):
                for delta in (-2,-1,1,2):
                    nv=E[idx]+delta
                    if nv<=0 or nv>=s or nv in E: continue
                    E2=sorted(E[:idx]+[nv]+E[idx+1:])
                    if not prim(E2) or longest_ap(E2)>7: continue
                    m2=cluster_margin(E2,s)
                    if m2<cur: cur=m2; E=E2; improved=True; break
                if improved: break
        if cur<best[0]: best=(cur,tuple(E))
    tag="  <-- APPROACHES 1!" if best[0]<1.1 else ""
    print(f"   s={s:>3}: adversarial MIN existence margin = {best[0]:.4f}{tag}   worst E={best[1]}")
print("\n=> if all mins stay comfortably >1, good-period existence is ROBUST for dissociated clusters:")
print("   the mac-mini c>=D3 sliver is a CERTIFICATE artifact, NOT a gap in the covering leg.")
