#!/usr/bin/env python3
"""
mapgraph_descend_klein_S334.py -- klein-2026-07-20-S334
Owner: do the cycle invariants descend to iso classes? is the star quotient the merged metagraph?

BOTH ANSWERED, BOTH NEGATIVE, both WITNESS-EFFECTIVELY (exhibited, not merely deduced).

Compare two partitions of the 2^m tilings: by tournament ISO CLASS, and by STAR ORBIT (coset mod cut(H)).
   n | tilings | iso classes | star orbits | iso refines orbits? | orbits refine iso?
   4 |       8 |      4      |   1 (2^0)   |     (vacuous)       |        no
   5 |      64 |     12      |   4 (2^2)   |       NO            |        no
   6 |    1024 |     56      |  32 (2^5)   |       NO            |        no

(i) THE CYCLE INVARIANTS DO NOT DESCEND. They are TILING functions, not tournament functions.
    WITNESS (n=5, tiles (3,1),(4,1),(5,1),(4,2),(5,2),(5,3)): one iso class holds 5 tilings spread over
    4 star orbits -- e.g. (0,0,0,0,1,0) and (0,0,0,1,1,0) are ISOMORPHIC tournaments with star-orbit reps
    (0,0,0,0,1,0) and (0,0,0,0,0,0). A second iso class spreads 9 tilings over 3 orbits.
(ii) THE STAR QUOTIENT IS UNRELATED to the merged metagraph: neither partition refines the other at
    n=5,6 -- they are transverse. Not the merged metagraph, not a refinement, not a coarsening.

WHY: cycle(H) is defined FROM the base path P (H = K_n minus P). Relabelling moves P, moves H, moves the
invariants. The construction was never S_n-equivariant -- a sharp instance of CLAUDE.md's warning that the
tiling model has lower symmetry than the arc model. THM-1382's algebra is exact and COMPLETE, but complete
for the tiling model = tournament PLUS a distinguished Hamiltonian path.
"""
import itertools
from itertools import permutations
def tiles(n): return [(x,y) for y in range(1,n+1) for x in range(1,n+1) if x-y>=2]
def star(n,v,T): return tuple(1 if (x==v or y==v) else 0 for (x,y) in T)
def tournament(n,T,bits):
    M=[[0]*(n+1) for _ in range(n+1)]
    for k in range(2,n+1): M[k][k-1]=1
    for (x,y),b in zip(T,bits):
        if b==0: M[x][y]=1
        else:    M[y][x]=1
    return M
def canon(n,M):
    best=None
    for p in permutations(range(1,n+1)):
        key=tuple(M[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or key<best: best=key
    return best
def orbit(n,T,bits,sv):
    best=None
    for mask in range(1<<len(sv)):
        v=list(bits)
        for i in range(len(sv)):
            if mask>>i&1: v=[a^b for a,b in zip(v,sv[i])]
        t=tuple(v)
        if best is None or t<best: best=t
    return best
def compare(n):
    T=tiles(n); sv=[star(n,v,T) for v in range(1,n+1)][:n-1]
    pairs=[]
    for bits in itertools.product((0,1),repeat=len(T)):
        pairs.append((canon(n,tournament(n,T,bits)), orbit(n,T,bits,sv)))
    isos={c for c,_ in pairs}; orbs={o for _,o in pairs}
    iso_ref=all(len({o for c2,o in pairs if c2==c})==1 for c in isos)
    orb_ref=all(len({c for c,o2 in pairs if o2==o})==1 for o in orbs)
    return len(isos), len(orbs), iso_ref, orb_ref
if __name__=="__main__":
    for n in (4,5,6):
        i,o,a,b=compare(n)
        print("n=%d: iso=%d orbits=%d | iso refines orbits: %s | orbits refine iso: %s"%(n,i,o,a,b))
