#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Task 2: is H_1^-(G_5) (dim 7) a Z_7/Fano-line structure under a Singer Z_7?  (klein-S6)

b_1^-(5)=7 is confirmed. Test whether the 7 R-odd cycles carry a genuine Singer Z_7 / Fano-plane
(octonion) structure, or whether 7 is coincidental. Three checks:
 (A) does Aut(G_5) contain an order-7 element (a Singer Z_7 acting on the 12 vertices: 7-orbit+5 fixed)?
 (B) if so, does it act on H_1^- as a single 7-cycle (regular rep) with Fano-line incidence?
 (C) where do the 7 R-odd cycles live -- are they tied to the 2 NS complement-pairs (the 'ribs')?
"""
from __future__ import annotations
import itertools
from math import comb
from collections import deque
import numpy as np

def pairs(n): return [(i, j) for i in range(n) for j in range(i+1, n)]
def perm_tables(n):
    P = pairs(n); idx = {p: k for k, p in enumerate(P)}; T = []
    for perm in itertools.permutations(range(n)):
        row = [( (idx[(perm[i],perm[j])],False) if perm[i]<perm[j] else (idx[(perm[j],perm[i])],True))
               for (i,j) in P]
        T.append(row)
    return T
def make_canon(T):
    def canon(bits):
        best=None
        for row in T:
            v=0
            for q,(t,inv) in enumerate(row):
                b=(bits>>t)&1
                if inv:b^=1
                v|=b<<q
            if best is None or v<best:best=v
        return best
    return canon

def build(n):
    d=comb(n,2);full=(1<<d)-1;canon=make_canon(perm_tables(n))
    s=canon(0);seen={s};q=deque([s]);edges=set();sig={}
    while q:
        c=q.popleft();sig[c]=canon(c^full)
        for a in range(d):
            nb=canon(c^(1<<a))
            if nb!=c:edges.add((min(c,nb),max(c,nb)))
            if nb not in seen:seen.add(nb);q.append(nb)
    cl=sorted(seen);ci={c:i for i,c in enumerate(cl)}
    E=[(ci[u],ci[v]) for (u,v) in edges]
    sigv=[ci[sig[c]] for c in cl]
    return cl,ci,sigv,E

def adj(V,E):
    A=[[0]*V for _ in range(V)]
    for u,v in E:A[u][v]=A[v][u]=1
    return A

def automorphisms(V,E):
    """all graph automorphisms (V=12 small): backtracking with degree refinement."""
    A=adj(V,E);deg=[sum(r) for r in A]
    nbr=[set(j for j in range(V) if A[i][j]) for i in range(V)]
    # refine by (deg, sorted neighbor degs)
    sig=[(deg[i],tuple(sorted(deg[j] for j in nbr[i]))) for i in range(V)]
    order=sorted(range(V),key=lambda i:(-deg[i]))
    autos=[]
    perm=[-1]*V;used=[False]*V
    def bt(k):
        if k==V:
            autos.append(perm[:]);return
        i=order[k]
        for j in range(V):
            if used[j] or sig[j]!=sig[i]:continue
            ok=True
            for ii in order[:k]:
                if A[i][ii]!=A[j][perm[ii]]:ok=False;break
            if ok:
                perm[i]=j;used[j]=True;bt(k+1);used[j]=False;perm[i]=-1
    bt(0)
    return autos

def order_of(p):
    V=len(p);o=1;q=p[:]
    while q!=list(range(V)):
        q=[p[q[i]] for i in range(V)];o+=1
        if o>10000:return -1
    return o

def cycle_space_basis(V,E):
    B=np.zeros((V,len(E)))
    for k,(u,v) in enumerate(E):B[u,k]-=1;B[v,k]+=1
    _,s,vt=np.linalg.svd(B);rank=int(np.sum(s>1e-9))
    return vt[rank:].T   # E x beta1

def signed_edge_perm(perm,E):
    """matrix of a vertex automorphism on edge space (oriented u<v), with orientation sign."""
    eidx={e:k for k,e in enumerate(E)};M=np.zeros((len(E),len(E)))
    for k,(u,v) in enumerate(E):
        a,b=perm[u],perm[v];t=(a,b) if a<b else (b,a);M[eidx[t],k]=1.0 if a<b else -1.0
    return M

if __name__=="__main__":
    n=5;cl,ci,sigv,E=build(n);V=len(cl)
    print("="*78);print(f" Task 2: Singer Z_7 / Fano test on H_1^-(G_{n}), dim 7  (klein-S6)");print("="*78)
    print(f" G_{n}: V={V}, E={len(E)}, SC(fixed by R)={sum(1 for i in range(V) if sigv[i]==i)}")
    autos=automorphisms(V,E)
    orders=sorted(set(order_of(p) for p in autos))
    print(f" |Aut(G_{n})| = {len(autos)};  element orders present = {orders}")
    has7=any(order_of(p)==7 for p in autos)
    print(f" (A) order-7 automorphism (Singer Z_7 on the 12 vertices) exists: {has7}")

    # R_* on H_1^-: confirm dim and structure
    ker=cycle_space_basis(V,E);beta1=ker.shape[1]
    R=signed_edge_perm(sigv,E)
    MR=ker.T@R@ker
    evR=np.linalg.eigvals(MR).real
    bminus=int(round(np.sum(evR<0)))
    print(f" H_1 dim beta1={beta1};  b_1^-={bminus}  (R_* -1-eigenspace)")

    # (B) if some automorphism of order 7 exists, act on H_1^- and test 7-cycle/Fano
    if has7:
        tau=next(p for p in autos if order_of(p)==7)
        T7=ker.T@signed_edge_perm(tau,E)@ker
        evt=np.linalg.eigvals(T7)
        print(f" (B) order-7 tau acts on H_1 with eigenvalues (|.|):近 7th-roots? "
              f"{sorted(round(abs(x),3) for x in evt)[:8]}")
    else:
        print(" (B) no order-7 automorphism -> no graph-Singer-Z_7; any Fano structure is NOT from graph symmetry.")

    # (C) localize the 7 R-odd cycles: NS complement-pairs (swapped vertices)
    swapped=[(i,sigv[i]) for i in range(V) if sigv[i]!=i]
    print(f" (C) R-swapped vertex pairs (NS complement-pairs): {len(swapped)//2} pairs "
          f"-> {sorted(set(tuple(sorted(p)) for p in swapped))}")
    # how much of H_1^- is supported away from SC backbone? (report b1- vs #NS pairs)
    print(f"     #NS pairs = {len(swapped)//2}; b_1^-={bminus}.")
    print("\n VERDICT: see whether (A)/(B) reveal a Singer Z_7; else 7 = dim coincidence (b_1^- = 0,1,7,119,1772,")
    print(" not a 7-pattern; 7|b_1^- holds n=5,6 but fails n=7).")
