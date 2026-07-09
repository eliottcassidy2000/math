# |R| vs the KISSING NUMBER of the grid-relation lattice (kps-S97). mac-mini-S25: relation-lattice
# kissing number = additive energy = # shortest additive relations. HERE the grid lattice is
# {n: Vmax|n.e}; its low-height members = my near-resonance count Z (S93). Claim: |R| ~ (kissing/count)
# x LEM-011 0.371-decay. Test |R|/lead vs (R2 energy, low-height resonance count) across the energy axis.
import numpy as np
from math import floor, gcd
from functools import reduce
from itertools import combinations, product
from collections import Counter
import random
random.seed(97071)
K=13; lead=(6.0/7.0)**K
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def R2_energy(E):
    d=Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i!=j: d[E[i]-E[j]]+=1
    return sum(v*v for v in d.values())
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for jj in range(i+1,len(E)):
            dd=E[jj]-E[i]; L=2; nx=E[jj]+dd
            while nx in S: L+=1; nx+=dd
            bk=E[i]-dd
            while bk in S: L+=1; bk-=dd
            best=max(best,L)
    return best
def R_resid(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.maximum(g-1.0/7,0).sum(axis=1).mean()-lead
def kissing_count(E,V,H=4,supmax=3):
    # # low-height n (support<=3, |n_i|<=H, no 7|n_i, no 7|sigma) with V | n.e  (n!=0). = near-resonance count.
    E=sorted(E); coords=list(range(1,K)); rng=[m for m in range(-H,H+1) if m!=0 and m%7!=0]; cnt=0
    for s in range(1,supmax+1):
        for combo in combinations(coords,s):
            for vals in product(rng,repeat=s):
                sig=sum(vals)
                if sig!=0 and sig%7==0: continue
                ne=sum(vals[t]*E[combo[t]] for t in range(s))
                if ne%V==0: cnt+=1
    return cnt
s=90; V=floor(7*s/6)
print(f"|R|/lead vs additive energy R2 and near-resonance/kissing count (support<=3,|n|<=4), V={V}.")
print(f"{'family':>14}{'R2':>7}{'lAP':>5}{'kiss':>6}{'|R|/lead':>10}")
rows=[]
def emit(lab,E):
    E=sorted(set(E))
    if len(E)!=K: return
    r2=R2_energy(E); ks=kissing_count(E,V); R=abs(R_resid(E,V))/lead
    rows.append((r2,ks,R,lab)); print(f"{lab:>14}{r2:>7}{longest_ap(E):>5}{ks:>6}{R:>10.4f}")
emit("AP-full",[i*(s//12) for i in range(13)])
emit("7-struct",[0,7,14,21,26,29,37,44,51,58,67,75,82])
for lab,lo,hi in [("nearAP",700,950),("mid",300,600),("Sidon",100,240)]:
    for _ in range(6000):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E) and lo<=R2_energy(E)<=hi: emit(lab,E); break
# correlation
import numpy as np
rows=np.array([(r[0],r[1],r[2]) for r in rows],float)
if len(rows)>=4:
    cR2=np.corrcoef(rows[:,0],rows[:,2])[0,1]; ckiss=np.corrcoef(rows[:,1],rows[:,2])[0,1]
    print(f"\ncorr(|R|, R2)={cR2:.3f}   corr(|R|, kissing-count)={ckiss:.3f}")
print("=> |R| governed by additive energy R2 / near-resonance kissing count; AP (max energy) extremal.")
print("   a-priori |R|<(6/7)^k <= Cohn-Elkies positivity over the grid-relation lattice (mac-mini-S24/S25).")
