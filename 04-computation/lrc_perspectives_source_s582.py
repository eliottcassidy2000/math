#!/usr/bin/env python3
"""The perspectives curiosity, nailed + the LRC source slice. For tournaments on n
vertices: #vertex-orbits (PERSPECTIVES) per iso class, and #SOURCE-orbits (observer can
sit there & be a source = LRC-relevant). Confirm perspectives sum (3+1=4, 4+4+2+2=12)
and source-perspectives = A000568(n-1) (THM-381). opus-2026-06-03-S582."""
from itertools import permutations, combinations
def tournaments(n):
    edges=list(combinations(range(n),2)); E=len(edges)
    for mask in range(2**E):
        adj=[[0]*n for _ in range(n)]
        for b,(i,j) in enumerate(edges):
            if mask>>b&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def canon(adj,n):
    best=None
    for p in permutations(range(n)):
        key=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or key<best: best=key
    return best
def vertex_orbits(adj,n):
    # orbits of the automorphism group on vertices
    autos=[]
    for p in permutations(range(n)):
        if all(adj[i][j]==adj[p[i]][p[j]] for i in range(n) for j in range(n)): autos.append(p)
    seen=set(); orbits=[]
    for v in range(n):
        if v in seen: continue
        orb=set(p[v] for p in autos); orbits.append(sorted(orb)); seen|=orb
    return orbits
def outdeg(adj,v,n): return sum(adj[v][j] for j in range(n))
def analyze(n):
    classes={}
    for adj in tournaments(n):
        c=canon(adj,n)
        if c not in classes: classes[c]=adj
    persp=0; srcpersp=0; detail=[]
    for c,adj in classes.items():
        orbs=vertex_orbits(adj,n); p=len(orbs)
        s=sum(1 for o in orbs if outdeg(adj,o[0],n)==n-1)  # source-orbits
        persp+=p; srcpersp+=s; detail.append((p,s))
    return len(classes),persp,srcpersp,sorted(detail,reverse=True)
A000568={1:1,2:1,3:2,4:4,5:12,6:56,7:456}
print("n: #iso-classes | PERSPECTIVES(sum vertex-orbits) | SOURCE-perspectives | per-class (persp,src)")
for n in [3,4,5,6]:
    cl,per,src,det=analyze(n)
    print(f"  n={n}: classes={cl} (A000568={A000568[n]}); perspectives={per} (=A000568(n+1)? {A000568.get(n+1)}); "
          f"source-persp={src} (=A000568(n-1)={A000568[n-1]}? {src==A000568[n-1]})")
    print(f"        per-class (persp,src): {det}")
