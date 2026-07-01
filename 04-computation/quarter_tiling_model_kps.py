#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THE QUARTER-TILING MODEL = Q_m / <sigma, phi>  (Klein 4-group: sigma=grid-reflection=tournament-complement,
phi=complement TILING=flip-all).  One fold below the half-tiling (Q_m/<sigma>).

kind-pasteur-2026-07-01-S15. The half-tiling folds by sigma (grid reflection = complement phi(T^op)); the
QUARTER folds again by phi (flip-all, fixed-point-free = the blue-line pairing).  Orbit structure under the
Klein 4-group {id, sigma, phi, sigma.phi}:
  * grid-sym tilings (sigma-fixed): orbit {t, phi t} has SIZE 2  = a BLUE LINE.
  * non-grid-sym tilings: orbit {t, sigma t, phi t, sigma.phi t} SIZE 4 = a 'black quad' (two sigma-linked
    complement pairs).  (phi, sigma.phi are fixed-point-free for n>=4, f>=1.)
Burnside: |quarter| = (2^m + Fix(sigma))/4 = (2^m + 2^{(m+f)/2})/4.
Also checks BLUE-SPINE CONNECTIVITY (is the SC/blue graph connected?).
"""
import sys, itertools
from collections import defaultdict
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def struct(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    f=sum(1 for i in range(m) if TRANS[i]==i)
    return TILES,m,TRANS,f

print("="*100); print(" QUARTER-TILING orbit counts (Burnside over Klein 4-group <sigma,phi>)"); print("="*100)
print(f"  {'n':>2} {'m':>3} {'f':>3} {'D=(m+f)/2':>9} {'tilings 2^m':>12} {'HALF Q/<s>':>12} "
      f"{'QUARTER Q/<s,phi>':>16} {'blue-lines 2^(D-1)':>17} {'black-quads':>11}")
for n in range(4,9):
    TILES,m,TRANS,f=struct(n); D=(m+f)//2
    half=(2**m + 2**D)//2         # orbits under <sigma>
    quarter=(2**m + 2**D)//4      # orbits under <sigma,phi>
    blue=2**(D-1)                 # size-2 orbits (grid-sym / 2)
    quads=(2**m - 2**D)//4        # size-4 orbits
    print(f"  {n:>2} {m:>3} {f:>3} {D:>9} {2**m:>12} {half:>12} {quarter:>16} {blue:>17} {quads:>11}"
          + ("  (blue+quads=quarter: %s)"%(blue+quads==quarter) if n<=6 else ""))

# verify by enumeration n<=6 + blue-spine connectivity
print("\n"+"="*100); print(" VERIFY orbit sizes by enumeration + BLUE-SPINE CONNECTIVITY (n=4,5,6)"); print("="*100)
def build_meta(n):
    VERTS=[n-i for i in range(n)]
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m=len(TILES); tileIdx={t:i for i,t in enumerate(TILES)}
    TRANS=[tileIdx[(n-y+1,n-x+1)] for (x,y) in TILES]
    def gsym(bits): return all(TRANS[i]==i or bits[i]==bits[TRANS[i]] for i in range(m))
    def sig(mask):
        bits=[(mask>>k)&1 for k in range(m)]; tb=[0]*m
        for i in range(m): tb[TRANS[i]]=bits[i]
        return sum(b<<k for k,b in enumerate(tb))
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
    FULL=(1<<m)-1
    # Klein-4 orbit sizes
    sizes=defaultdict(int); seen=set()
    for mask in range(1<<m):
        if mask in seen: continue
        orb={mask, mask^FULL, sig(mask), sig(mask)^FULL}
        for x in orb: seen.add(x)
        sizes[len(orb)]+=1
    # metagraph for blue-spine connectivity
    T=[]
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        T.append(dict(mask=mask,canon=canon(adj(bits)),g=gsym(bits),fl=mask^FULL))
    sigs=sorted(set(t['canon'] for t in T)); ci={s:i for i,s in enumerate(sigs)}
    for t in T: t['ci']=ci[t['canon']]
    bym={t['mask']:t for t in T}
    for t in T: t['tci']=bym[sig(t['mask'])]['ci']
    tgt={c:[t for t in T if t['ci']==c][0]['tci'] for c in range(len(sigs))}
    par=list(range(len(sigs)))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for c in range(len(sigs)):
        a,b=find(c),find(tgt[c])
        if a!=b: par[max(a,b)]=min(a,b)
    noc=[find(c) for c in range(len(sigs))]
    isSC={nd:tgt[nd]==nd for nd in set(noc)}
    # blue graph among SC nodes
    blue_adj=defaultdict(set); seen2=set()
    for t in T:
        pr=frozenset((t['mask'],t['fl']))
        if pr in seen2: continue
        seen2.add(pr)
        if not t['g']: continue
        a=noc[t['ci']]; b=noc[bym[t['fl']]['ci']]
        if a!=b: blue_adj[a].add(b); blue_adj[b].add(a)
    SC=[nd for nd in set(noc) if isSC[nd]]
    # connectivity of blue graph on SC nodes
    if SC:
        comp=set(); stack=[SC[0]]
        while stack:
            u=stack.pop()
            if u in comp: continue
            comp.add(u); stack+=[v for v in blue_adj[u] if v not in comp]
        conn = comp>=set(SC)
    else: conn=True
    return sizes, len(SC), conn, blue_adj
for n in [4,5,6]:
    sizes,nSC,conn,ba=build_meta(n)
    print(f"  n={n}: Klein-4 orbit sizes {dict(sizes)}  | SC nodes={nSC}; BLUE graph on SC connected? {conn}")
print("\nSYNTHESIS: half=fold by sigma (tournament complement); QUARTER=fold again by phi (tiling complement,")
print("  flip-all, fpf). Blue lines = the SIZE-2 quarter-orbits; black quads = size-4. The blue subgraph =")
print("  the quarter-tiling's own complement structure (recursion). Blue-spine connectivity checked above.")
print("DONE.")
