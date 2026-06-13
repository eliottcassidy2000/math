#!/usr/bin/env python3
"""
paley_starstar_topological_factor_monad.py
monad-explorer-2026-06-07 (deep-research, 7th session)

THE TOPOLOGICAL FACTORIZATION of (star-star).  Builds on THM-438 ADDENDUM-4.

Every even-series pattern sigma has e = #series-classes (= #distinct edge-flow-lines
= #edges E_H of its TOPOLOGICAL reduction H, suppressing degree-2 chain vertices).
The even chain-length assignment (each series-class has length 2c, sum 2c = 2k, i.e.
sum c = k over the e classes) is INDEPENDENT of the mu-weight prod_v(|B_v|-1)! and of
the cycle rank m.  So group by e:

   G(k,e) := sum_{even-series sigma, #series-classes = e} (-1)^{m} prod_v(|B_v|-1)!

CLAIM (topological factorization):   G(k,e) = g_e * C(k-1, e-1),   g_e k-INDEPENDENT.

and then  S_k = sum_e G(k,e) = sum_e g_e C(k-1,e-1),  with (from binomial inversion of
S_k=(-1)^k C_k):   g_e = -1, 3, -10, 36, -137, 543, ...  = (-1)^e * A002212(e).

This script:
 (A) verifies the factorization G(k,e) = g_e C(k-1,e-1) directly (NON-tautological:
     it tests that the grouped sums really factor as binomial * constant), and
 (B) confirms g_e = (-1)^e A002212(e).

If (A)+(B) hold, then (star-star) reduces to the single sequence identity
g_e = (-1)^e A002212(e) PLUS the (classical, GF-provable) binomial-transform identity
sum_e (-1)^e A002212(e) C(k-1,e-1) = (-1)^k C_k.  No primes, no characters.
"""
import sys, math
from collections import defaultdict
import numpy as np

def catalan(k): return math.comb(2*k,k)//(k+1)

# A002212: 1,1,3,10,36,137,543,2219,9285,...  (a(0)=1). g.f. A: xA^2-(1-3x)A+1? we just hardcode.
A002212 = [1,1,3,10,36,137,543,2219,9285,39587]

def edge_flow_groups(edges, nb):
    """Return (#series-classes e, cycle rank m) or None if not even-series."""
    E=len(edges); Bm=np.zeros((nb,E))
    for ei,(u,v) in enumerate(edges):
        Bm[v,ei]+=1.0; Bm[u,ei]-=1.0
    u,s,vh=np.linalg.svd(Bm,full_matrices=True)
    tol=1e-9; rank=int((s>tol).sum()); m=E-rank
    if m==0: return None
    ns=vh[rank:]
    groups=defaultdict(int)
    for e in range(E):
        v=ns[:,e]
        if np.max(np.abs(v))<1e-7: return None  # bridge
        v=v/np.max(np.abs(v))
        for x in v:
            if abs(x)>1e-7:
                if x<0: v=-v
                break
        key=tuple(round(float(x),5) for x in v)
        groups[key]+=1
    if not all(c%2==0 for c in groups.values()): return None
    return len(groups), m   # e = #series classes, m = cycle rank

def gen_rgs(L):
    n=L+1; a=[0]*n; mx=[0]*n
    def rec(i):
        if i==n:
            yield a[:]; return
        for v in range(mx[i-1]+2):
            if v==a[i-1]: continue
            a[i]=v; mx[i]=max(mx[i-1],v); yield from rec(i+1)
    a[0]=0; mx[0]=0
    if n==1: yield [0]; return
    yield from rec(1)

KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 5
print("="*72)
print("TOPOLOGICAL FACTORIZATION  G(k,e) = g_e * C(k-1,e-1),  g_e=(-1)^e A002212(e)")
print("="*72)
all_ok=True
for k in range(1,KMAX+1):
    L=2*k
    G=defaultdict(int)
    for a in gen_rgs(L):
        nb=max(a)+1
        edges=[(a[i],a[i+1]) for i in range(L)]
        # connectivity
        adj=defaultdict(set)
        for (u,v) in edges: adj[u].add(v); adj[v].add(u)
        seen={0}; stk=[0]
        while stk:
            x=stk.pop()
            for w in adj[x]:
                if w not in seen: seen.add(w); stk.append(w)
        if len(seen)!=nb: continue
        res=edge_flow_groups(edges,nb)
        if res is None: continue
        e,m=res
        sz=defaultdict(int)
        for x in a: sz[x]+=1
        w=1
        for b in sz.values(): w*=math.factorial(b-1)
        G[(e)] += ((-1)**m)*w
    line=[]
    for e in sorted(G):
        g_pred=((-1)**e)*A002212[e]
        binom=math.comb(k-1,e-1)
        expect=g_pred*binom
        ok = (G[e]==expect)
        all_ok = all_ok and ok
        line.append(f"e={e}: G={G[e]} = g_{e}({g_pred})*C({k-1},{e-1})={expect} {'OK' if ok else 'XX'}")
    Stot=sum(G.values()); tgt=(-1)**k*catalan(k)
    print(f"k={k}  S_k={Stot} (tgt {tgt})")
    for s in line: print("      ", s)
print("="*72)
print(f"Factorization + A002212 identification holds at all tested k: {all_ok}")
print("=> (star-star) REDUCED to:  g_e := sum_{topo maps, E_H=e}(-1)^m prod(b-1)! = (-1)^e A002212(e)")
print("   (a finite-per-e topological-map identity), plus the classical binomial transform")
print("   sum_e (-1)^e A002212(e) C(k-1,e-1) = (-1)^k C_k  (verified separately by GF).")
