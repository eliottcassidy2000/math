"""TASK 2: SC-spine blue-degree vs blue-graph spectrum. Blue subgraph = SC nodes (grid-sym=SC classes) + blue
flip-lines. Compute blue-degree (odd, T-join) and the adjacency SPECTRUM; compare (Sum lam^2 = 2|E|+diag, radius
vs max deg, integrality, gap)."""
import numpy as np
from itertools import permutations
def build_blue(n):
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
    # grid-sym tilings (SC); nodes=SC classes; blue edges = flip-pairs
    gs=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        if all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i): gs.append((mask,bits))
    cls={}; 
    def cid(bits):
        c=canon(adj(bits))
        if c not in cls: cls[c]=len(cls)
        return cls[c]
    node={}; 
    for mask,bits in gs: node[mask]=cid(bits)
    N=len(cls); B=np.zeros((N,N))
    seen=set()
    for mask,bits in gs:
        fm=mask^((1<<m)-1)
        key=(min(mask,fm),max(mask,fm))
        if key in seen: continue
        seen.add(key)
        a=node[mask]; b=node[fm]
        if a==b: B[a,a]+=1      # self-loop
        else: B[a,b]+=1; B[b,a]+=1
    return B,N
for n in [4,5,6,7]:
    B,N=build_blue(n)
    deg=(B.sum(1)-np.diag(B)).astype(int)      # edge-degree (exclude self-loop from off-diag sum already; diag=self-loops)
    offdeg=(B - np.diag(np.diag(B))).sum(1).astype(int)  # pure edge degree
    ev=np.sort(np.linalg.eigvalsh(B))
    evr=np.round(ev,4)
    E=int((B-np.diag(np.diag(B))).sum()/2); SL=int(np.diag(B).sum())
    integral=np.allclose(evr, np.round(evr))
    print(f"n={n}: {N} SC nodes, {E} blue edges, {SL} blue self-loops")
    print(f"  blue edge-degrees (odd T-join): {sorted(offdeg.tolist())}")
    print(f"  blue spectrum: {evr.tolist()}")
    print(f"  radius {round(ev[-1],3)} vs maxdeg {int(offdeg.max())}; sum lam^2={round((ev**2).sum(),2)} vs 2E+sum(self^2)={2*E+int((np.diag(B)**2).sum())}; integer spectrum? {integral}")
