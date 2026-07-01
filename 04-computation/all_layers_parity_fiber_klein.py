#!/usr/bin/env python3
"""
all_layers_parity_fiber_klein.py  --  klein-2026-07-01-S76

SYNTHESIS: unify HYP-1772 (d=1 wiggly edge-balance) + klein-S75 (d=m blue/black lines) into an ALL-LAYERS
parity law on the merged metagraph, and study the FIBER (which tilings share a node) via F(C)=H(C)/|Aut(C)|.

ALL-LAYERS PARITY LAW (conjecture): for each Hamming layer d (flip d tiles) and merged node u,
    tau(u,d) := #{(t,t'): t in u, dist(t,t')=d, node(t')!=u}  satisfies  2*lambda + tau = C(m,d)*M_u,
so  tau(u,d) ≡ C(m,d) * M_u  (mod 2).  Since M_u ≡ [u is SC] (HYP-1772) and C(m,d) mod 2 = Lucas(m,d),
    tau(u,d) is ODD  <=>  (u is SC)  AND  (C(m,d) is odd).
This makes HYP-1772 (d=1: C(m,1)=m) and klein-S75 (d=m: C(m,m)=1 => all SC have odd cross-degree) two
special cases of ONE Lucas-graded law.

FIBER: F(C) = #tilings mapping to class C = H(C)/|Aut(C)| (HYP-1772), ODD. Verify + study which tilings
share a node (the Ham-path/Aut orbit structure).
"""
import itertools, math
from collections import Counter, defaultdict

def build(n):
    pairs=[(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    pidx={p:k for k,p in enumerate(pairs)}
    tiles=[(x,y) for x in range(1,n+1) for y in range(1,x) if x-y>=2]
    perms=list(itertools.permutations(range(1,n+1)))
    return pairs,pidx,tiles,perms
def tourney_mask(tv,n,tiles,pidx):
    A=[[0]*(n+1) for _ in range(n+1)]
    for a in range(2,n+1): A[a][a-1]=1
    for b,(x,y) in enumerate(tiles):
        if (tv>>b)&1: A[x][y]=1
        else: A[y][x]=1
    mask=0
    for (i,j) in [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]:
        if A[i][j]: mask|=1<<pidx[(i,j)]
    return mask
def canon(mask,n,pairs,pidx,perms):
    best=None
    for pi in perms:
        v=0
        for k,(i,j) in enumerate(pairs):
            u,w=(i,j) if ((mask>>k)&1) else (j,i)
            a,b=pi[u-1],pi[w-1]
            if a<b: v|=1<<pidx[(a,b)]
        if best is None or v<best: best=v
    return best
def opp(mask,pairs,pidx):
    v=0
    for k,(i,j) in enumerate(pairs):
        if not((mask>>k)&1): v|=1<<pidx[(i,j)]
    return v
def ham_paths(mask,n,pairs,pidx):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k,(i,j) in enumerate(pairs):
        if (mask>>k)&1: A[i][j]=1
        else: A[j][i]=1
    cnt=0
    for perm in itertools.permutations(range(1,n+1)):
        if all(A[perm[t]][perm[t+1]] for t in range(n-1)): cnt+=1
    return cnt
def aut(mask,n,pairs,pidx,perms):
    return sum(1 for pi in perms if _apply(mask,pi,pairs,pidx)==mask)
def _apply(mask,pi,pairs,pidx):
    v=0
    for k,(i,j) in enumerate(pairs):
        u,w=(i,j) if ((mask>>k)&1) else (j,i)
        a,b=pi[u-1],pi[w-1]
        if a<b: v|=1<<pidx[(a,b)]
    return v

if __name__=="__main__":
    for n in [4,5,6]:
        pairs,pidx,tiles,perms=build(n); m=len(tiles)
        node=[0]*(1<<m); is_sc={}
        cmask={}
        for tv in range(1<<m):
            mk=tourney_mask(tv,n,tiles,pidx); c=canon(mk,n,pairs,pidx,perms); co=canon(opp(mk,pairs,pidx),n,pairs,pidx,perms)
            key=min(c,co); node[tv]=key; is_sc[key]=(c==co); cmask[key]=c
        M=Counter(node)                      # merged fiber (mass)
        keys=list(M)
        print(f"\n===== n={n}: m={m}, merged nodes {len(keys)}, #SC={sum(is_sc.values())} =====")
        # ALL-LAYERS parity law
        # tau(u,d): sum over t in u of #{t': dist=d, node!=u}
        tau=defaultdict(lambda: [0]*(m+1))
        for t in range(1<<m):
            u=node[t]
            for tp in range(1<<m):
                if node[tp]!=u:
                    d=bin(t^tp).count("1"); tau[u][d]+=1
        viol=0
        for u in keys:
            for d in range(1,m+1):
                pred=(math.comb(m,d)%2)*(M[u]%2)
                if tau[u][d]%2 != pred: viol+=1
        # which layers have any odd-cross-incidence, and are they exactly the SC nodes?
        lucas_odd_d=[d for d in range(1,m+1) if math.comb(m,d)%2==1]
        print(f"  ALL-LAYERS parity: tau(u,d)≡C(m,d)*M_u mod2 -> violations={viol} (0 => LAW HOLDS)")
        print(f"    C(m,d) odd at layers d in {lucas_odd_d} (Lucas); at those d, odd-cross = SC nodes; else all even")
        # check d=1 (HYP-1772) and d=m (klein-S75) specialize
        odd1=[u for u in keys if tau[u][1]%2==1]; oddm=[u for u in keys if tau[u][m]%2==1]
        print(f"    d=1 odd-cross nodes: {len(odd1)} (=SC iff m odd, m={m} {'odd' if m%2 else 'even'}); d=m odd-cross: {len(oddm)} (=#SC={sum(is_sc.values())})")
        # FIBER = H/|Aut|, odd
        fib_ok=True; fodd=True
        for u in keys[:60]:
            c=cmask[u]; H=ham_paths(c,n,pairs,pidx); A=aut(c,n,pairs,pidx,perms)
            Fclass=sum(1 for t in range(1<<m) if canon(tourney_mask(t,n,tiles,pidx),n,pairs,pidx,perms)==c)
            if Fclass!=H//A or H%A!=0: fib_ok=False
            if Fclass%2==0: fodd=False
        print(f"  FIBER F(C)=H/|Aut| holds: {fib_ok}; F(C) odd: {fodd} (=> merged mass parity = [SC], HYP-1772 proved via |Aut| odd)")
