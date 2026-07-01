#!/usr/bin/env python3
"""mac-mini-2026-07-01-S90b -- the flip-rank EXCESS lower bound via super-symmetric SC classes.
Owner's hypothesis: excess(n) >= #{SC classes with |Aut| > n} (each forces a covering dimension), and the
super-symmetric (|Aut|>n) classes are SELF-COMPLEMENTARY (T-join parity x rarity).
Compute per iso class (via the tiling model): |Aut|, SC(self-complementary) status; tabulate #{SC & |Aut|>n},
#{|Aut|>n}, and whether every |Aut|>n class is SC. Compare to excess = 0,0,0,1,3 (n=3..7)."""
import itertools
from math import comb
from collections import Counter
def analyze(n):
    VERTS=list(range(n,0,-1)); TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); vpos={v:i for i,v in enumerate(VERTS)}; perms=list(itertools.permutations(range(n)))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi,yi=vpos[xL],vpos[yL]
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    def canon(A): return min(tuple(A[p[i]][p[j]] for i in range(n) for j in range(n)) for p in perms)
    seen={}
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]; A=adj(bits); c=canon(A)
        if c not in seen: seen[c]=A
    # per class: |Aut|, SC (canon(complement)==canon)
    def aut(A): return sum(1 for p in perms if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j))
    def comp(A): return [[0 if i==j else 1-A[i][j] for j in range(n)] for i in range(n)]
    rows=[]
    for c,A in seen.items():
        a=aut(A); sc=(canon(comp(A))==c); rows.append((a,sc))
    nclass=len(rows)
    super_sym=[r for r in rows if r[0]>n]
    sc_super=[r for r in super_sym if r[1]]
    all_super_sc=all(r[1] for r in super_sym)
    return dict(n=n, nclass=nclass, autdist=dict(Counter(r[0] for r in rows)),
                n_super=len(super_sym), n_sc_super=len(sc_super), all_super_sc=all_super_sc,
                nsc=sum(1 for r in rows if r[1]))
excess={3:0,4:0,5:0,6:1,7:3}
print("n | #class | #SC | |Aut|>n count | #{SC & |Aut|>n} | all |Aut|>n are SC? | excess")
for n in range(4,8):
    d=analyze(n)
    print(f"{n} | {d['nclass']:4d} | {d['nsc']:3d} | {d['n_super']:2d} | {d['n_sc_super']:2d} | {d['all_super_sc']} | {excess.get(n,'?')}")
    print(f"    |Aut| distribution: {dict(sorted(d['autdist'].items()))}")
