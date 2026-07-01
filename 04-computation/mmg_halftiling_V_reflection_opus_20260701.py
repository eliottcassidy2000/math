"""Symmetries of the half-tiling model. V: x->(n+3)-x (vertical refl of half-region). Is V class-preserving
(a relabeling), the complement, or a NEW fold? Also enumerate all reflection symmetries + the (span,anti-diag) grid."""
from itertools import permutations
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    ti={t:i for i,t in enumerate(TILES)}; TRANS=[ti[(n-y+1,n-x+1)] for (x,y) in TILES]
    perms=list(permutations(range(n)))
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
    half=[i for i,(x,y) in enumerate(TILES) if x+y<=n+1]     # transversal reps + fixed
    Vmap={}   # V on ALL tiles restricted meaning: (x,y)->(n+3-x,y)
    okV=True
    for i,(x,y) in enumerate(TILES):
        img=(n+3-x,y)
        if img in ti: Vmap[i]=ti[img]
        else: Vmap[i]=None
    # V must preserve half-tile set
    Vhalf_ok=all(Vmap[i] is not None and Vmap[i] in half for i in half)
    # apply V to a grid-sym tiling: bits over half determine full grid-sym tiling
    def unfold(hbits):  # hbits: dict tile-> bit for half; produce full bits (grid-sym)
        b=[0]*m
        for i in range(m):
            j = i if TRANS[i]==i else min(i,TRANS[i])
            b[i]=hbits[j]
        return b
    def applyV_full(bits):  # V permutes tiles; new bits: b'[V(i)]=b[i]
        nb=[0]*m
        for i in range(m):
            if Vmap[i] is not None: nb[Vmap[i]]=bits[i]
        return nb
    # test V on all grid-sym tilings
    import itertools
    Hbits=[i for i in half]
    same=comp=other=0
    for mask in range(1<<len(half)):
        hb={half[k]:(mask>>k)&1 for k in range(len(half))}
        b=unfold(hb); 
        # grid-sym check
        if not all(b[i]==b[TRANS[i]] for i in range(m) if TRANS[i]!=i): continue
        A=adj(b); c=canon(A)
        bv=applyV_full(b)
        # bv should still be grid-sym (V preserves half)
        Av=adj(bv); cv=canon(Av); cop=canon([[Av[j][i] for j in range(n)] for i in range(n)])
        c_cop=canon([[A[j][i] for j in range(n)] for i in range(n)])
        if cv==c: same+=1
        elif cv==c_cop: comp+=1
        else: other+=1
    fixedV=sum(1 for i in half if Vmap[i]==i)
    print(f"n={n}: |half|={len(half)}, V preserves half-tiles? {Vhalf_ok}; V-fixed half-tiles={fixedV}")
    print(f"   V action on grid-sym classes: same(relabel)={same}, complement={comp}, other={other}")
