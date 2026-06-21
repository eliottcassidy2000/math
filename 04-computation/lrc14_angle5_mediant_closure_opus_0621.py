#!/usr/bin/env python3
"""Is the binder characterized by MEDIANT/additive closure (Stern-Brocot / Frobenius
coin structure)? Test: among all P of given size, does the binder maximize the number
of additive relations a+b=c (the 'most Frobenius-degenerate' set)?"""
from fractions import Fraction as F
from itertools import combinations
def measGP(P):
    bad=[]
    for p in P:
        p=int(p)
        for kk in range(0,p+1):
            lo=max((F(kk)-F(1,14))/p,F(0)); hi=min((F(kk)+F(1,14))/p,F(1))
            if lo<hi: bad.append((lo,hi))
    bad.sort(); tot=F(0); cur=cb=None
    for lo,hi in bad:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return F(1)-tot
def nrel(P):
    P=sorted(P)
    return sum(1 for i in range(len(P)) for j in range(i+1,len(P)) for l in range(len(P)) if P[i]+P[j]==P[l])
for sz in [3,4,5]:
    rows=sorted((measGP(list(P)),P) for P in combinations(range(1,14),sz))
    m,binder=rows[0]
    maxrel=max(nrel(P) for P in combinations(range(1,14),sz))
    print(f"|P|={sz}: binder={sorted(binder)} nrel={nrel(binder)} (max nrel at size={maxrel}); binder maximizes additive-relations: {nrel(binder)==maxrel}")
    # also: does binder consist of {1} + a near-AP top cluster?
    print(f"   top-cluster check: binder = {sorted(binder)}")
