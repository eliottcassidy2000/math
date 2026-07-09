# Is the grid resonance residual |R| controlled by the ADDITIVE ENERGY R2? (kps-S97)
# THM-660: Var(W) ~ 5.67e-5*R2, R2=sum_{d!=0} r_d^2 (r_d=#{(i,j):e_i-e_j=d}); AP maximizes R2, Sidon minimizes.
# If |R|/(6/7)^k <= g(R2) cleanly, the a-priori |R|<(6/7)^k reduces to an R2 (energy) bound (THM-660 supplies).
# Sweep clusters from Sidon (low R2) to AP (high R2) at critical V; plot |R|/(6/7)^k vs R2.
import numpy as np
from math import floor, gcd
from functools import reduce
from collections import Counter
import random
random.seed(97013)
K=13; lead=(6.0/7.0)**K
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def R2_energy(E):
    diffs=Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i!=j: diffs[E[i]-E[j]]+=1
    return sum(v*v for v in diffs.values())   # sum_{d!=0} r_d^2
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
def R_resid(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    W=np.maximum(g-1.0/7,0).sum(axis=1)
    return W.mean()-lead     # = R (grid residual)
print("R2 (additive energy) vs |R|/(6/7)^k at critical V=floor(7s/6). AP=max R2, Sidon=min R2.")
print(f"R2(AP_13)=k(k-1)(2k-1)/3={13*12*25//3}. (6/7)^13={lead:.5f}.")
print(f"{'family':>16}{'R2':>7}{'longAP':>7}{'|R|/lead':>10}")
s=90; V=floor(7*s/6)
# span the energy axis
fams=[]
fams.append(("AP-full",[i*(s//12) for i in range(13)]))
fams.append(("7-struct",[0,7,14,21,26,29,37,44,51,58,67,75,82]))
for lab,lo,hi in [("near-AP",700,900),("mid",300,600),("low-energy",150,300),("Sidon-ish",100,220)]:
    best=None
    for _ in range(4000):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))!=K or not prim(E): continue
        r2=R2_energy(E)
        if lo<=r2<=hi:
            best=E; break
    if best: fams.append((lab,best))
for lab,E in fams:
    E=sorted(set(E))
    if len(E)!=K: continue
    r2=R2_energy(E); R=abs(R_resid(E,V))/lead
    print(f"{lab:>16}{r2:>7}{longest_ap(E):>7}{R:>10.4f}")
print()
print("=> if |R|/lead rises monotonically with R2 and stays <1 for dissociated (low R2),")
print("   the a-priori |R|<(6/7)^k reduces to an R2 (additive-energy) threshold -- THM-660 territory.")
