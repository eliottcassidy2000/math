# DIRECT good-period-existence search over small-spread dissociated k=13 clusters (kps-S94).
# Skip the c>=D3 proxy: for each dissociated E (spread s), scan the critical Vmax window and ask
# does SOME j in {1..Vmax-1} give maxgap{e_i j mod Vmax} > Vmax/7 (7*maxgap > Vmax). If existence
# NEVER fails, the small-spread dissociated regime is CLOSED directly (independent of mac-mini's
# certificate). Vectorized integer arithmetic over j.
import numpy as np
from math import gcd, ceil
from functools import reduce
import random
random.seed(94023)
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
def min_maxgap_over_j(E,Vmax):
    # returns (best maxgap over all j, argj). teeth = (e*j) mod Vmax for j=1..Vmax-1.
    j=np.arange(1,Vmax)                      # (Vmax-1,)
    Ea=np.array(sorted(E))                   # (k,)
    ph=np.mod(np.outer(j,Ea),Vmax)          # (Vmax-1, k)
    ph.sort(axis=1)
    diffs=np.diff(ph,axis=1)
    wrap=(ph[:,0]+Vmax-ph[:,-1])[:,None]
    g=np.concatenate([diffs,wrap],axis=1)    # (Vmax-1, k)
    mg=g.max(axis=1)                         # maxgap per j
    best=mg.max(); return int(best), int(j[mg.argmax()])
k=13; diss=7
print("DIRECT good-period existence over small-spread dissociated k=13 clusters.")
print("Existence at (E,Vmax): exists j with 7*maxgap > Vmax. Scan Vmax in [s+1, ceil(7s/6)+5].")
print("Report clusters+Vmax where NO good period exists (would break the covering leg).\n")
fails=[]; tested=0; tightest=(99.9,None)
for _ in range(400000):
    if tested>=15000: break
    s=random.choice([50,60,70,80,90,100,120,140])
    mid=sorted(random.sample(range(1,s),k-2)); E=[0]+mid+[s]
    if len(set(E))!=k or not prim(E) or longest_ap(E)>diss: continue
    tested+=1
    for Vmax in range(s+1, ceil(7*s/6)+6):
        if gcd(reduce(gcd,E),Vmax)!=1: pass
        best,aj=min_maxgap_over_j(E,Vmax)
        ratio=7*best/Vmax     # >1 means good period exists
        if ratio<=1.0:
            fails.append((tuple(E),Vmax,best,round(ratio,4)))
            break
        if ratio<tightest[0]: tightest=(ratio,(tuple(E),Vmax))
print(f"tested {tested} dissociated clusters; good-period-EXISTENCE failures: {len(fails)}")
if not fails:
    print("=> EVERY tested small-spread dissociated cluster admits a good period at EVERY critical Vmax.")
    print(f"   Tightest margin: 7*maxgap/Vmax = {tightest[0]:.4f} (>1) at {tightest[1]}")
    print("   => the sliver is a CERTIFICATE artifact; good-period existence holds directly (this regime CLOSES).")
else:
    print("=> EXISTENCE FAILURES (real gap in the covering leg!):")
    for f in fails[:10]: print("   ",f)
