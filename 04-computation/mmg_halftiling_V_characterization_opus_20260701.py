"""Characterize V (vertical refl of half-region, x->(n+3)-x): (1) is it WELL-DEFINED on SC classes (respects
tournament-iso)? (2) V-fixed half-tilings + classes. (3) confirm V in (s,d) swaps span<->anti-diagonal."""
from itertools import permutations
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES); ti={t:i for i,t in enumerate(TILES)}
    TRANS=[ti[(n-y+1,n-x+1)] for (x,y) in TILES]; perms=list(permutations(range(n)))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL; A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
        return A
    def canon(A):
        b=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if b is None or s<b: b=s
        return b
    return n,m,TILES,ti,TRANS,adj,canon
for n in [5,6,7]:
    n,m,TILES,ti,TRANS,adj,canon=build(n)
    half=[i for i,(x,y) in enumerate(TILES) if x+y<=n+1]
    Vmap=[ti.get((n+3-x,y)) for (x,y) in TILES]
    def unfold(hb): return [hb[i if TRANS[i]==i else min(i,TRANS[i])] for i in range(m)]
    def applyV(bits):
        nb=[0]*m
        for i in range(m):
            if Vmap[i] is not None: nb[Vmap[i]]=bits[i]
        return nb
    # group grid-sym tilings by class; check if V maps a class consistently
    from collections import defaultdict
    cls_of={}; members=defaultdict(list); Vfix=0
    tilings=[]
    for mask in range(1<<len(half)):
        hb={half[k]:(mask>>k)&1 for k in range(len(half))}; b=unfold(hb)
        if not all(b[i]==b[TRANS[i]] for i in range(m) if TRANS[i]!=i): continue
        A=adj(b); c=canon(A); cls_of[tuple(b)]=c; members[c].append(tuple(b)); tilings.append(tuple(b))
        vb=tuple(applyV(list(b)))
        if vb==tuple(b): Vfix+=1
    # well-defined on classes? for each class, do all members map to the SAME target class?
    welldef=True
    for c,mem in members.items():
        tgt=set(canon(adj(list(applyV(list(b))))) for b in mem)
        if len(tgt)>1: welldef=False; break
    # V-fixed classes (classes with a V-fixed tiling)
    vfix_classes=set(canon(adj(list(b))) for b in tilings if tuple(applyV(list(b)))==tuple(b))
    print(f"n={n}: #SC classes={len(members)}, grid-sym tilings={len(tilings)}, V-fixed tilings={Vfix}")
    print(f"   V well-defined on SC classes (respects tournament-iso)? {welldef}")
    print(f"   => V is {'a METAGRAPH symmetry' if welldef else 'ONLY a GEOMETRIC symmetry of the tile-cube (does NOT respect iso)'}")
    print(f"   #classes touched by a V-fixed tiling: {len(vfix_classes)}")
print("\n(s,d)=(x+y, x-y): R:(s,d)->(2n+2-s, d) [fixes span d]; V:(s,d)->(n+3-d, n+3-s) [SWAPS s<->d].")
print("Tile lattice = perpendicular families: constant-span d (slope +1) x constant-antidiag s (slope -1), sublattice s=d mod2.")
