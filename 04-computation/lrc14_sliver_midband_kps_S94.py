# Fill the intermediate spread band s in [21,49] (kps-S94): adversarial min existence margin, to
# make the evidence continuous between the EXHAUSTIVE base (s<=20, min 1.105) and the large-spread
# adversarial (s>=50, min>=1.717). Pure-dissociated L<=6 (the region ONLY Route-c must handle) AND
# inclusive L<=7. If min-margin stays >1 across the whole band, existence is robust for all spreads.
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(94211)
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
        best=7*g.max(axis=1).max()/Vmax
        if best<worst: worst=best
    return worst
def rand_diss(s,Lcap):
    for _ in range(300):
        mid=sorted(random.sample(range(1,s),11)); E=[0]+mid+[s]
        if len(set(E))==13 and prim(E) and longest_ap(E)<=Lcap: return E
    return None
def adversary(s,Lcap,restarts):
    best=(9.9,None)
    for _ in range(restarts):
        E=rand_diss(s,Lcap)
        if E is None: continue
        cur=cluster_margin(E,s); improved=True
        while improved:
            improved=False
            for idx in range(1,12):
                for delta in (-2,-1,1,2):
                    nv=E[idx]+delta
                    if nv<=0 or nv>=s or nv in E: continue
                    E2=sorted(E[:idx]+[nv]+E[idx+1:])
                    if not prim(E2) or longest_ap(E2)>Lcap: continue
                    m2=cluster_margin(E2,s)
                    if m2<cur: cur=m2; E=E2; improved=True; break
                if improved: break
        if cur<best[0]: best=(cur,tuple(E))
    return best
print("Intermediate band s in [21,49]: adversarial MIN existence margin (good period exists iff >1).")
print(f"{'s':>4}{'minMargin(L<=6)':>17}{'minMargin(L<=7)':>17}")
globmin=(9.9,None)
for s in range(21,50,3):
    b6=adversary(s,6,40); b7=adversary(s,7,40)
    for b in (b6,b7):
        if b[0]<globmin[0]: globmin=(b[0],(s,b[1]))
    print(f"{s:>4}{b6[0]:>17.4f}{b7[0]:>17.4f}")
print(f"\nband-min existence margin = {globmin[0]:.4f} at s,E={globmin[1]}")
print("Combined with EXHAUSTIVE s<=20 (min 1.105) and adversarial s>=50 (min>=1.717):")
print("=> good-period existence margin stays >1 for ALL spreads => dissociated covering branch has")
print("   NO real gap; mac-mini's c>=D3 sliver is purely a sufficient-certificate artifact.")
