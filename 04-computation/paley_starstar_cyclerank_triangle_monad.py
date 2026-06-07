#!/usr/bin/env python3
"""
paley_starstar_cyclerank_triangle_monad.py
monad-explorer-2026-06-07 (deep-research, 7th session)

Extend the CYCLE-RANK triangle  t(k,m) = sum_{even-series sigma, rank m} prod_v (|B_v|-1)!
(THM-438 ADDENDUM-3) to as many k as feasible, to expose a recursion whose alternating
row sum  sum_m (-1)^m t(k,m) = (-1)^k C_k  (i.e. (star-star)).

Optimised enumerator: restricted-growth strings (RGS) for partitions of {0..2k} with
EARLY PRUNING  (no two consecutive positions equal => a[i]!=a[i-1]).  Then cheap filters
(connected, min vertex multiplicity, no bridge proxy) before the SVD even-series test.
"""
import sys, math
from collections import defaultdict
import numpy as np

def catalan(k): return math.comb(2*k,k)//(k+1)

def edge_flow_lines(edges, nb):
    E=len(edges); Bm=np.zeros((nb,E))
    for ei,(u,v) in enumerate(edges):
        Bm[v,ei]+=1.0; Bm[u,ei]-=1.0
    u,s,vh=np.linalg.svd(Bm,full_matrices=True)
    tol=1e-9; rank=int((s>tol).sum()); m=E-rank
    if m==0: return None,0
    ns=vh[rank:]; lines=[]
    for e in range(E):
        v=ns[:,e]
        if np.max(np.abs(v))<1e-7: lines.append(("Z",)); continue
        v=v/np.max(np.abs(v))
        for x in v:
            if abs(x)>1e-7:
                if x<0: v=-v
                break
        lines.append(tuple(round(float(x),5) for x in v))
    return lines,m

def even_series_rank(a, nb, L):
    # a: RGS assignment length L+1 (block id per position). edges:
    edges=[(a[i],a[i+1]) for i in range(L)]
    # connected
    adj=defaultdict(set)
    for (u,v) in edges: adj[u].add(v); adj[v].add(u)
    seen={0}; stk=[0]
    while stk:
        x=stk.pop()
        for w in adj[x]:
            if w not in seen: seen.add(w); stk.append(w)
    if len(seen)!=nb: return None
    lines,m=edge_flow_lines(edges,nb)
    if m==0: return None
    if any(l==("Z",) for l in lines): return None
    groups=defaultdict(int)
    for l in lines: groups[l]+=1
    if not all(c%2==0 for c in groups.values()): return None
    return m

def gen_rgs(L):
    """Yield restricted-growth strings of length L+1 with a[i]!=a[i-1], a[0]=0."""
    n=L+1
    a=[0]*n
    mx=[0]*n  # mx[i]=max used in a[0..i]
    def rec(i):
        if i==n:
            yield a[:]; return
        hi=mx[i-1]+1
        for v in range(hi+1):
            if v==a[i-1]: continue   # no self-loop
            a[i]=v; mx[i]=max(mx[i-1],v)
            yield from rec(i+1)
    a[0]=0; mx[0]=0
    if n==1:
        yield [0]; return
    yield from rec(1)

KMAX=int(sys.argv[1]) if len(sys.argv)>1 else 5
print("CYCLE-RANK TRIANGLE  t(k,m)=sum_{even-series rank m} prod(|B|-1)!")
for k in range(1,KMAX+1):
    L=2*k
    tri=defaultdict(int); cnt=0; S=0
    for a in gen_rgs(L):
        nb=max(a)+1
        m=even_series_rank(a,nb,L)
        if m is None: continue
        # weight = prod (|B|-1)!  ; block sizes:
        sz=defaultdict(int)
        for x in a: sz[x]+=1
        w=1
        for b in sz.values(): w*=math.factorial(b-1)
        tri[m]+=w; cnt+=1
        S += ((-1)**m)*w
    row=" ".join(str(tri[m]) for m in range(1,max(tri)+1))
    tgt=(-1)**k*catalan(k)
    print(f"k={k}: [{row}]   altsum={S}  (-1)^kC_k={tgt}  match={S==tgt}  #ev-ser={cnt}")
