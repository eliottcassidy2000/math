#!/usr/bin/env python3
"""mac-mini-S83: TEST my 'joint rigidity across facets' idea. Claim: covering breaks the facets
(Schur/E3, tournament Q(2)=E[2^X]), and JOINTLY the breaks lower-bound L. RED FLAG (S79): the
smallest-L family has the LARGEST Schur deficit. If facet-breaks and L are ANTI-correlated, no
joint of them can bound L. Measure L vs the facet-deficits over covering families; compute the
correlation. If L does NOT increase with the joint break (or anti-correlates), joint rigidity FAILS."""
import numpy as np
from math import comb
c=1.0/14
def facets_and_L(S,res=300000):
    S=sorted(set(S)); k=len(S)
    ts=(np.arange(res)+0.5)/res
    X=np.zeros(res,dtype=int)
    for v in S:
        r=(v*ts)%1.0; d=np.minimum(r,1-r); X+=(d<c).astype(int)
    L=(X==0).mean()
    Q2=(2.0**X).mean()      # tournament partition E[2^X]
    Aset=set(S); T=sum(1 for a in S for b in S if a+b in Aset)  # Schur
    return L,T,Q2
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    from math import gcd
    g=0
    for v in S: g=gcd(g,v)
    return g==1
# AP baseline
_,T_ap,Q_ap=facets_and_L(list(range(1,14)))
print(f"AP: Schur={T_ap}=C(13,2), Q(2)={Q_ap:.1f}\n")
# gather covering families: near-AP perturbations + deep-well variants + random covering
fams=[]
import random; random.seed(1)
fams+=[[*range(1,12),13,84],[*range(1,13),182],[*range(1,13),364],[*range(1,11),22,13,84],
       [*range(1,12),13,168],list(range(2,15)),[*range(1,12),13,84*2],[*range(1,10),20,11,13,84]]
tries=0
while len([f for f in fams])<28 and tries<4000:
    tries+=1
    S=sorted(set([1]+random.sample([v for v in range(2,120)],12)))
    if len(S)==13 and prim(S) and is_cov(S): fams.append(S)
rows=[]
for S in fams:
    S=sorted(set(S))
    if len(S)!=13 or not(prim(S) and is_cov(S)): continue
    L,T,Q2=facets_and_L(S)
    rows.append((L, T_ap-T, Q_ap-Q2, S))  # L, Schur-deficit, Q2-deficit
rows.sort()
print(f"{'L (small=hard)':>13} | Schur-def | Q2-def | family")
for L,dS,dQ,S in rows[:14]:
    print(f"{L:13.5f} | {dS:9d} | {dQ:6.1f} | {S[:5]}...{S[-1]}")
Ls=np.array([r[0] for r in rows]); dS=np.array([r[1] for r in rows]); dQ=np.array([r[2] for r in rows])
print(f"\ncorr(L, Schur-deficit) = {np.corrcoef(Ls,dS)[0,1]:+.3f}")
print(f"corr(L, Q2-deficit)    = {np.corrcoef(Ls,dQ)[0,1]:+.3f}")
print("\nVERDICT: if corr(L, facet-deficit) is NEGATIVE or near 0, the facet-breaks do NOT predict")
print("L (bigger structural break can mean SMALLER L). Then NO joint of the facets lower-bounds L")
print("=> the joint rigidity FAILS: the facets measure structural distance-from-AP, L measures")
print("safe-set openness, and they DISAGREE. The residue is irreducibly metric, not even jointly")
print("reducible to the structural facets. Honest end of the 'joint rigidity' idea.")
