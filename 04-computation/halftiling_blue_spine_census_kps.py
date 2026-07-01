#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
BLUE-SPINE via the HALF-TILING: the SC/blue structure lives ONLY on the 2^D grid-sym tilings.

kind-pasteur-2026-07-01-S16. New mental model: to study the SC spine we do NOT need the full 2^m tiling cube --
only the 2^D grid-sym (=half-tiling) tilings, D=floor((n-1)^2/4).  This makes n=7 (512 half-tilings) and even
n=8 cheap.  Verify BLUE-SPINE CONNECTIVITY at n=7; extend the census (#SC, #pure-blue, #mixed, g_C odd) to n=7.
 - grid-sym tiling t => class is SC (transpose-self).  g_C = #grid-sym tilings in class C.  total tiling count
   = H/|Aut|.  PURE-BLUE iff g_C == H/|Aut| (no non-grid-sym tilings), else MIXED.
 - BLUE graph: edge class(t) -- class(phi t) (phi = flip-all = antipode of the half-tiling cube).
"""
import sys, itertools
from collections import defaultdict
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def analyze(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    # sigma-orbits of tiles
    seen=set(); orbits=[]
    for i in range(m):
        if i in seen: continue
        o=[i] if TRANS[i]==i else [i,TRANS[i]]
        for x in o: seen.add(x)
        orbits.append(o)
    D=len(orbits)
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=VERTS.index(xL); yi=VERTS.index(yL)
            if bits[i]==0: A[xi][yi]=1
            else: A[yi][xi]=1
        return A
    P=list(itertools.permutations(range(n)))
    def canon(A):
        best=None
        for p in P:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    def H(A): return sum(1 for p in P if all(A[p[t]][p[t+1]] for t in range(n-1)))
    def autc(A): return sum(1 for p in P if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)))
    # enumerate grid-sym tilings
    gs=[]
    for a in range(1<<D):
        bits=[0]*m
        for oi,orb in enumerate(orbits):
            v=(a>>oi)&1
            for t in orb: bits[t]=v
        gs.append((a,bits))
    # map each grid-sym tiling (by its assignment-int a over orbits) to class
    amask={}
    for a,bits in gs:
        amask[a]=bits
    def antipode_a(a): return a ^ ((1<<D)-1)   # flip all orbit-bits = flip all tiles = phi
    canon_of={}; gcount=defaultdict(int); repbits={}
    for a,bits in gs:
        c=canon(adj(bits)); canon_of[a]=c; gcount[c]+=1
        if c not in repbits: repbits[c]=bits
    SCclasses=list(gcount.keys())
    # H,|Aut|,total per class
    info={}
    for c in SCclasses:
        A=adj(repbits[c]); Hc=H(A); ac=autc(A); info[c]=(Hc,ac,Hc//ac)
    # pure-blue vs mixed
    pure=[c for c in SCclasses if gcount[c]==info[c][2]]
    mixed=[c for c in SCclasses if gcount[c]< info[c][2]]
    allodd=all(g%2==1 for g in gcount.values())
    # blue graph
    badj=defaultdict(set)
    for a,bits in gs:
        c1=canon_of[a]; c2=canon_of[antipode_a(a)]
        if c1!=c2: badj[c1].add(c2); badj[c2].add(c1)
    # connectivity
    comp=set(); stack=[SCclasses[0]]
    while stack:
        u=stack.pop()
        if u in comp: continue
        comp.add(u); stack+=[v for v in badj[u] if v not in comp]
    conn = comp>=set(SCclasses)
    # self-loops (antipode in same class)
    selfloops=sum(1 for a,bits in gs if canon_of[a]==canon_of[antipode_a(a)])//2
    return dict(n=n,m=m,D=D,nGS=len(gs),nSC=len(SCclasses),npure=len(pure),nmix=len(mixed),
                allodd=allodd,conn=conn,selfloops=selfloops,gcounts=sorted(gcount.values()))

print("="*96); print(" SC/BLUE structure via the HALF-TILING (only 2^D grid-sym tilings needed)"); print("="*96)
print(f"  {'n':>2} {'D':>3} {'2^D=#half':>10} {'#SC':>5} {'#pure-blue':>11} {'#mixed':>7} "
      f"{'g_C all odd':>11} {'BLUE-SPINE conn':>15} {'blue self-loops':>15}")
for n in [4,5,6,7]:
    r=analyze(n)
    print(f"  {n:>2} {r['D']:>3} {r['nGS']:>10} {r['nSC']:>5} {r['npure']:>11} {r['nmix']:>7} "
          f"{str(r['allodd']):>11} {str(r['conn']):>15} {r['selfloops']:>15}")
print("\n  => extends the census to n=7 CHEAPLY (512 half-tilings vs 32768 full).")
print("     #pure-blue = 1,3,2,? ; blue-spine connectivity tested at n=7.")
print("DONE.")
