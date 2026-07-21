#!/usr/bin/env python3
"""
TOURNAMENTS AS RECURSIVELY COMPOSED OF SEEDS (opus-S444).
Owner: tournaments are composed of recursive smaller subtournaments; what iso-class SEEDS
correspond to larger tournaments? Connect to the three recursion modes + the octonion (n=9) wall.

MODULAR/SUBSTITUTION decomposition. A MODULE M of T: a subset with every outside vertex uniform to
all of M (beats all, or loses to all). PRIME (simple/indecomposable) = only trivial modules
(singletons and V). SUBSTITUTION T[S_1..S_k]: blow up vertex i of T into a copy of S_i; between
blocks arcs follow T. The PRIMES are the SEEDS; every tournament is a substitution tree over seeds
(Gallai/modular decomposition). The transitive tournament is the fully-decomposable (linear) one;
the 3-cycle C3 is the smallest nontrivial seed.

Compute: (1) the seed census n=3..6; (2) char_S / SC4 / H composition under substitution;
(3) the octonion object C3[C3] (n=9); (4) patterns.
"""
import itertools, numpy as np
from math import comb

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj

def canon(adj,n):
    best=None
    for p in itertools.permutations(range(n)):
        key=0; b=0
        for i in range(n):
            for j in range(n):
                if i!=j: key|=adj[p[i]][p[j]]<<b; b+=1
        if best is None or key<best: best=key
    return best

def is_module(adj,n,M):
    Mset=set(M); out=[v for v in range(n) if v not in Mset]
    for v in out:
        # v must be uniform to M: all M beaten by v, or all beat v
        beats=[adj[v][m] for m in M]
        if len(set(beats))>1: return False
    return True

def is_prime(adj,n):
    for r in range(2,n):
        for M in itertools.combinations(range(n),r):
            if is_module(adj,n,M): return False
    return True

def substitute(T,k,blocks):
    """T on k vertices; blocks[i]=(adjS_i, m_i). Return big adj on sum m_i vertices."""
    sizes=[b[1] for b in blocks]; off=[0];
    for s in sizes: off.append(off[-1]+s)
    N=off[-1]; big=[[0]*N for _ in range(N)]
    for i in range(k):
        Si,mi=blocks[i]
        for a in range(mi):
            for c in range(mi):
                if a!=c: big[off[i]+a][off[i]+c]=Si[a][c]
        for j in range(k):
            if i!=j and T[i][j]:
                for a in range(sizes[i]):
                    for c in range(sizes[j]):
                        big[off[i]+a][off[j]+c]=1
    return big,N

def charS(adj,n):
    S=np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    return tuple(int(round(c)) for c in np.poly(S).real)
def sc4(adj,n):
    cnt=0
    for q in itertools.combinations(range(n),4):
        c=sum(1 for i,j,k in itertools.combinations(q,3)
              if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]))
        if c==2: cnt+=1
    return cnt
def Hcount(adj,n):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                for u in range(n):
                    if not (mask>>u)&1 and adj[v][u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[size-1])

# ---- (1) seed census ----
print("(1) SEED (prime/indecomposable) tournament census")
print("="*60)
seedreps={}
for n in range(1,7):
    seen=set(); primes=[]; total=0
    for adj in edges_iter(n) if n>=2 else [[[0]]]:
        c=canon(adj,n) if n>=2 else 0
        if c in seen: continue
        seen.add(c); total+=1
        if n==1 or is_prime(adj,n): primes.append(adj)
    seedreps[n]=primes
    print(f"  n={n}: {total} iso classes, {len(primes)} PRIME (seeds)")
print("  seed counts n=1..6:", [len(seedreps[n]) for n in range(1,7)], " (OEIS: simple tournaments)")

C3=[[0,1,0],[0,0,1],[1,0,0]]                       # the 3-cycle seed
print(f"\n  C3 prime? {is_prime(C3,3)};  transitive prime? {is_prime([[0,1,1],[0,0,1],[0,0,0]],3)}")

# ---- (2) substitution composition laws ----
print("\n(2) SUBSTITUTION composition: C3[S,S,S] (blow up the 3-cycle by S)")
print("="*60)
for name,S,m in [("C1",[[0]],1),("C3",C3,3)]:
    big,N=substitute(C3,3,[(S,m)]*3)
    print(f"  C3[{name}^3]: N={N}, regular? {len(set(sum(big[i]) for i in range(N)))==1} (out-deg {sum(big[0])}),"
          f" charS={charS(big,N)}, SC4={sc4(big,N)}")

# ---- (3) the octonion object C3[C3] at n=9 ----
print("\n(3) THE OCTONION OBJECT C3[C3] (n=9)")
print("="*60)
big,N=substitute(C3,3,[(C3,3)]*3)
sc=sum(big[0]); reg=len(set(sum(big[i]) for i in range(N)))==1
print(f"  C3[C3]: N=9, regular={reg} (out-deg {sc}); charS={charS(big,9)}")
print(f"  SC4(C3[C3]) = {sc4(big,9)}; H(C3[C3]) = {Hcount(big,9)}")
# compare var(lambda^2)
S=np.array([[0.0 if i==j else (1.0 if big[i][j] else -1.0) for j in range(9)] for i in range(9)])
l2=sorted(abs(e)**2 for e in np.linalg.eigvals(S))
print(f"  lambda^2 multiset = {[round(v,3) for v in l2]};  var(lambda^2) = {round(float(np.var(l2)),4)}")
print(f"  (transitive n=9 var = 2*C(9,3) = {2*comb(9,3)}; this measures its GIT-instability)")

# ---- (4) how SC4 composes: SC4(C3[S^3]) vs SC4 within/spanning ----
print("\n(4) PATTERN: SC4 under substitution (the var/octonion-relevant census)")
print("="*60)
for name,S,m in [("C1",[[0]],1),("C3",C3,3)]:
    big,N=substitute(C3,3,[(S,m)]*3)
    within=3*sc4(S,m) if m>=4 else 0
    print(f"  C3[{name}^3] (N={N}): SC4={sc4(big,N)}  (per-block SC4={sc4(S,m) if m>=4 else 0})")
print("\n READING: the PRIME seeds are C3 and the larger simple tournaments; C3[C3] (n=9) is the")
print(" substitution-SQUARE of the smallest seed, a regular tournament at the octonion level.")
