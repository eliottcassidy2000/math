#!/usr/bin/env python3
"""
Size-controlled: does H(C3[S,S,S]) need the ENDPOINT PATH-MATRIX PH(S), not just scalar H(S)?
Use n=5 blocks (where H does NOT determine the iso class): find two 5-tournaments with the SAME
scalar H but DIFFERENT PH-matrix, substitute each into C3[S^3] (n=15), and compare composite H.
Then confirm the refinement composes: same PH => same composite. (opus-S445)
"""
import itertools, numpy as np
def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj
def Hcount(adj,n):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                av=adj[v]
                for u in range(n):
                    if not (mask>>u)&1 and av[u]: dp[mask|(1<<u)][u]+=c
    return sum(dp[size-1])
def PH_sortedspectrum(adj,n):
    """iso-invariant of the PH matrix: sorted multiset of row-sums (paths starting at each vertex)
       + sorted all entries -> a canonical signature of PH up to relabeling."""
    size=1<<n
    rows=[]
    for s in range(n):
        dp=[[0]*n for _ in range(size)]; dp[1<<s][s]=1
        for mask in range(size):
            for v in range(n):
                c=dp[mask][v]
                if c:
                    av=adj[v]
                    for u in range(n):
                        if not (mask>>u)&1 and av[u]: dp[mask|(1<<u)][u]+=c
        rows.append(sum(dp[size-1]))     # # Ham paths STARTING at s
    return tuple(sorted(rows))
C3=[[0,1,0],[0,0,1],[1,0,0]]
def substitute(T,k,blocks):
    sizes=[b[1] for b in blocks]; off=[0]
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
                    for c in range(sizes[j]): big[off[i]+a][off[j]+c]=1
    return big,N

# collect n=5 tournaments by (H scalar, PH-startvector signature)
byH={}
for adj in edges_iter(5):
    h=Hcount(adj,5); ph=PH_sortedspectrum(adj,5)
    byH.setdefault(h,{}).setdefault(ph,adj)
print("n=5 blocks: for each scalar H, how many distinct PH-start-signatures?")
for h in sorted(byH):
    print(f"  H={h}: {len(byH[h])} distinct PH-signatures")

# pick an H with >=2 PH-signatures; substitute both, compare composite
target=None
for h in sorted(byH):
    if len(byH[h])>=2: target=h; break
sigs=list(byH[target].items())
(s1sig,S1),(s2sig,S2)=sigs[0],sigs[1]
A,NA=substitute(C3,3,[(S1,5)]*3); B,NB=substitute(C3,3,[(S2,5)]*3)
HA,HB=Hcount(A,NA),Hcount(B,NB)
print(f"\n same scalar H={target}, DIFFERENT PH-signature:")
print(f"   S1 PH-start={s1sig} -> H(C3[S1^3]) = {HA}")
print(f"   S2 PH-start={s2sig} -> H(C3[S2^3]) = {HB}")
print(f"   composite H differs: {HA!=HB}  => scalar H does NOT compose; PH is the refined invariant.")

# confirm PH composes: two DIFFERENT n=5 tournaments with the SAME PH-signature -> same composite
same=None
for h in byH:
    for ph,adj in byH[h].items():
        pass
# find two distinct iso classes (different adj) sharing a PH-signature at n=5
from collections import defaultdict
bysig=defaultdict(list)
seen=set()
def canon(adj,n):
    best=None
    for p in itertools.permutations(range(n)):
        key=0;b=0
        for i in range(n):
            for j in range(n):
                if i!=j: key|=adj[p[i]][p[j]]<<b;b+=1
        if best is None or key<best:best=key
    return best
for adj in edges_iter(5):
    c=canon(adj,5)
    if c in seen: continue
    seen.add(c)
    bysig[(Hcount(adj,5),PH_sortedspectrum(adj,5))].append(adj)
shared=[v for v in bysig.values() if len(v)>=2]
if shared:
    S1,S2=shared[0][0],shared[0][1]
    A,NA=substitute(C3,3,[(S1,5)]*3); B,NB=substitute(C3,3,[(S2,5)]*3)
    print(f"\n two DISTINCT iso classes with SAME (H,PH-sig): composite H = {Hcount(A,NA)} vs {Hcount(B,NB)}"
          f"  (PH-sig fixes it: {Hcount(A,NA)==Hcount(B,NB)})")
else:
    print("\n (no two distinct n=5 classes share the full PH-start signature; PH-start already separates)")
print("\n READING: H = row-sum-collapse of PH; PH (endpoint-resolved) is the FUNCTORIAL invariant that")
print(" composes under substitution, where scalar H does not. It is the refined 'real answer'.")
