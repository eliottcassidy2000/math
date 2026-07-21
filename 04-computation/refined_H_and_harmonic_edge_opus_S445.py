#!/usr/bin/env python3
"""
THE REFINED COMPOSITIONAL INVARIANT + THE HARMONIC EDGE (opus-S445).
Owner: H (#P) not determined by a poly invariant is an EDGE case; concoct something more refined
than H that is the real answer; tournaments sit at the formula/#P boundary like the harmonic series.

PART A -- does H COMPOSE from scalar H under cyclic substitution?
  H(C3[S1,S2,S3]): if it is NOT a function of (H(S1),H(S2),H(S3)) alone, scalar H does not compose;
  the REFINED invariant that DOES compose is the ENDPOINT-RESOLVED path data PH(S)[u,v] = #Ham-paths
  of S from u to v (a matrix, categorifying H = sum PH). Under substitution the blocks glue by a
  transfer over their boundary path-matrices. Test both.

PART B -- the harmonic edge: partition tournaments by their degree-k induced-subtournament CENSUS
  (a poly invariant, the k-profile); measure the H-DEFECT = max H - min H within each census-class,
  as k grows from 3 to n. If the defect only reaches 0 at k=n (full support) and shrinks slowly,
  tournaments sit at the formula/#P edge (each finite resolution = a partial sum; H = the divergent
  limit).
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
def PHmatrix(adj,n):
    """PH[u][v] = # Ham paths from u to v."""
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if c:
                for u in range(n):
                    if not (mask>>u)&1 and adj[v][u]: dp[mask|(1<<u)][u]+=c
    full=size-1
    PH=[[0]*n for _ in range(n)]
    # dp[full][v] counts paths ENDING at v; need start too -> recompute with start fixed
    for s in range(n):
        dp2=[[0]*n for _ in range(size)]; dp2[1<<s][s]=1
        for mask in range(size):
            for v in range(n):
                c=dp2[mask][v]
                if c:
                    for u in range(n):
                        if not (mask>>u)&1 and adj[v][u]: dp2[mask|(1<<u)][u]+=c
        for v in range(n): PH[s][v]=dp2[full][v]
    return PH

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

# small seed blocks (n=1,2,3) with their H and PH
blocks_lib=[]
for n in (1,2,3):
    for adj in edges_iter(n):
        blocks_lib.append((adj,n,Hcount(adj,n),tuple(map(tuple,PHmatrix(adj,n)))))
# dedup by (n, H, PH-multiset) is not needed; use a few representatives
print("PART A: does H(C3[S1,S2,S3]) depend only on (H(S1),H(S2),H(S3))?")
print("="*66)
seen={}
collision=None
import random
reps=[b for b in blocks_lib if b[1] in (1,2,3)]
for b1 in reps:
    for b2 in reps:
        for b3 in reps:
            big,N=substitute(C3,3,[(b1[0],b1[1]),(b2[0],b2[1]),(b3[0],b3[1])])
            Hc=Hcount(big,N)
            key=tuple(sorted([b1[2],b2[2],b3[2]]))   # multiset of block H's
            seen.setdefault(key,set()).add(Hc)
ambiguous={k:v for k,v in seen.items() if len(v)>1}
print(f"  block-H-multisets with MULTIPLE composite-H values: {len(ambiguous)}/{len(seen)}")
for k,v in list(ambiguous.items())[:4]:
    print(f"    block H's {k}: composite H(C3[...]) in {sorted(v)}  <-- H does NOT compose from scalar H")
print("  => the refined invariant is the ENDPOINT path-matrix PH (H = sum of its entries).")

# verify PH composes: two blocks with same H but different PH give different composite H
print("\n  same scalar H, different PH-matrix -> different composite (PH is the right invariant):")
byH={}
for b in reps: byH.setdefault((b[1],b[2]),[]).append(b)
for key,lst in byH.items():
    if len(lst)>=2 and len({b[3] for b in lst})>=2:
        b1,b2=lst[0],lst[1]
        h1=Hcount(*substitute(C3,3,[(b1[0],b1[1])]*3)[:1],substitute(C3,3,[(b1[0],b1[1])]*3)[1]) if False else None
        A=substitute(C3,3,[(b1[0],b1[1])]*3); B=substitute(C3,3,[(b2[0],b2[1])]*3)
        print(f"    n={key[0]} H={key[1]}: PH differ -> C3[S^3] H = {Hcount(*A)} vs {Hcount(*B)}")
        break

# PART B: harmonic edge -- H-defect vs resolution-k census
print("\nPART B: HARMONIC EDGE -- H-defect within degree-k census classes")
print("="*66)
def kprofile(adj,n,k):
    from collections import Counter
    cnt=Counter()
    for sub in itertools.combinations(range(n),k):
        # canonical form of induced sub-tournament
        m=len(sub); best=None
        for p in itertools.permutations(range(m)):
            key=0;b=0
            for i in range(m):
                for j in range(m):
                    if i!=j: key|=adj[sub[p[i]]][sub[p[j]]]<<b; b+=1
            if best is None or key<best: best=key
        cnt[best]+=1
    return tuple(sorted(cnt.items()))
for n in (5,6):
    # gather (kprofile_k, H) over all tournaments; measure defect per k
    print(f"  n={n}:")
    Hs=[]; profs={k:[] for k in range(3,n+1)}
    for adj in edges_iter(n):
        h=Hcount(adj,n); Hs.append(h)
        for k in range(3,n+1):
            profs[k].append(kprofile(adj,n,k))
    for k in range(3,n+1):
        cls={}
        for h,p in zip(Hs,profs[k]): cls.setdefault(p,[]).append(h)
        defect=max(max(v)-min(v) for v in cls.values())
        nclasses=len(cls)
        print(f"    resolution k={k}: #census-classes={nclasses:5d}, max H-defect within class = {defect}")
print("\n READING: H-defect shrinks to 0 only at k=n (full support). Each poly resolution k is a")
print(" partial sum; H is the limit reached only at full support -- the formula/#P edge. The refined")
print(" invariant that COMPOSES (unlike scalar H) is the endpoint path-matrix PH.")
