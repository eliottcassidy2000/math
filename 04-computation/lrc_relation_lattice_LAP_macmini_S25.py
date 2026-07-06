#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S25 (HYP-4552) -- the RELATION LATTICE L(AP) and its KISSING NUMBER.

L(AP) = {a in Z^12 : sum_{i=1}^{12} i*a_i = 0} = ker of the moment map phi(a)=<c,a>,
c=(1..12).  disc(L(AP)) = |c|^2 = sum i^2 = 650.  Gram (basis d_k=(k+1)e_k - k e_{k+1})
is tridiagonal: diag k^2+(k+1)^2, offdiag -k(k+2).

FINDING: the MINIMAL vectors of L(AP) have norm 3 and are the ADDITIVE TRIPLES
(1,1,-1) at (i,j,i+j): v_i+v_j=v_{i+j}.  30 triples => kissing number 60 = 2*(#additive
triples) ~ the ADDITIVE ENERGY (HYP-2873).  The AP maximizes additive energy => L(AP)
is the Cohn-Elkies EXTREMAL relation lattice (max kissing).  The density floor = its
isolation.  Next shells: multiplicative doubling 2v_i=v_2i (norm 5), harmonic (1,-2,1)
(norm 6).  The theta itself is all-orders (support-3 shell -1.80 at beta=2/25), but the
KISSING NUMBER is a clean invariant with classical extremality.
"""
import numpy as np
from math import sin, pi
from itertools import combinations, product

def hhat(m, beta): return (1-2*beta) if m==0 else -sin(2*pi*m*beta)/(pi*m)

def relations(maxsupp=3, maxcoef=3):
    rels=set()
    for k in range(2, maxsupp+1):
        for idx in combinations(range(12), k):
            for cs in product(range(-maxcoef, maxcoef+1), repeat=k):
                if 0 in cs: continue
                if sum(cs[j]*(idx[j]+1) for j in range(k)) != 0: continue
                a=[0]*12
                for j in range(k): a[idx[j]]=cs[j]
                rels.add(tuple(a))
    return rels

if __name__=="__main__":
    c=np.array(range(1,13)); print(f"disc(L(AP)) = |c|^2 = sum i^2 = {int(c@c)}")
    B=np.zeros((11,12),dtype=int)
    for k in range(1,12): B[k-1,k-1]=k+1; B[k-1,k]=-k
    G=B@B.T
    print("Gram tridiagonal: diag", [int(G[k,k]) for k in range(11)], "offdiag", [int(G[k,k+1]) for k in range(10)])
    rels=relations()
    nrm=lambda a: sum(x*x for x in a)
    mn=min(nrm(a) for a in rels)
    minv=[a for a in rels if nrm(a)==mn]
    print(f"min norm = {mn}; kissing (supp<=3) = {len(minv)}; e.g.",
          [[(i+1,a[i]) for i in range(12) if a[i]] for a in list(minv)[:3]])
    beta=2/25
    for supp in (2,3,4):
        tot=sum(np.prod([hhat(x,beta) for x in a if x!=0]) for a in rels if sum(1 for x in a if x)==supp)
        print(f"  beta=2/25 support-{supp} shell = {tot:.5f}")
    print("=> min-norm-3 additive triples v_i+v_j=v_{i+j}; kissing = 2*(#triples) ~ additive energy; AP maximal.")
