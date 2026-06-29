#!/usr/bin/env python3
"""
affine-compression-search.py

Face search (axis-aligned) gave minimum covering dimensions:
  n=4 -> 2 (info bound 2, TIGHT)
  n=5 -> 4 (info bound 4, TIGHT)
  n=6 -> 7 (info bound 6, NOT tight: +1 overhead)

Question: does a general GF(2)-AFFINE slice (a coset of a linear subspace of the
tiling space, i.e. fixing linear combinations of tiles rather than single tiles)
recover information-theoretic tightness? This reframes "compression" as coding
theory: tilings = points of GF(2)^m, the iso-class map = a colouring, and we want
the smallest-dimension affine subspace meeting every colour (a covering coset).

We reuse the exact iso-class colouring from the face script (recomputed here) and
do a bounded randomized search for a dim-d affine cover, d = info bound.
"""
import itertools, sys
from math import comb

# ---- minimal copy of the exact iso-class colouring machinery (n<=6) ----
def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def base_path_pairs(n): return [(k,k+1) for k in range(n-1)]
def tile_pairs(n):
    bp=set(base_path_pairs(n)); return [p for p in pairs(n) if p not in bp]
def tiling_bits(n,flips):
    idx={p:k for k,p in enumerate(pairs(n))}; b=0
    for p in flips: b|=1<<idx[p]
    return b
def perm_tables(n):
    P=pairs(n); idx={p:k for k,p in enumerate(P)}; T=[]
    for perm in itertools.permutations(range(n)):
        row=[]
        for (i,j) in P:
            a,b=perm[i],perm[j]
            row.append((idx[(a,b)],False) if a<b else (idx[(b,a)],True))
        T.append(row)
    return T
def canonical(bits,T):
    best=None
    for row in T:
        v=0
        for q,(tgt,inv) in enumerate(row):
            bit=(bits>>tgt)&1
            if inv: bit^=1
            v|=bit<<q
        if best is None or v<best: best=v
    return best

def colouring(n):
    tiles=tile_pairs(n); m=len(tiles); T=perm_tables(n)
    col={}
    for mask in range(2**m):
        flips=[tiles[k] for k in range(m) if (mask>>k)&1]
        col[mask]=canonical(tiling_bits(n,flips),T)
    return m,col,len(set(col.values()))

# ---- deterministic LCG (Math.random is unavailable; keep it reproducible) ----
class LCG:
    def __init__(self,seed=12345): self.s=seed
    def next(self):
        self.s=(self.s*6364136223846793005+1442695040888963407)&((1<<64)-1)
        return self.s
    def randint(self,n): return self.next()%n

def affine_cover_search(n, m, col, nclasses, dim, tries=40000, seed=1):
    """Random dim-dimensional GF(2)-affine subspaces; report first full cover."""
    full=set(col.values()); rng=LCG(seed)
    for t in range(tries):
        # random offset + dim basis vectors in GF(2)^m
        v0=rng.randint(1<<m)
        basis=[]
        # ensure independence by rejection
        span={0}
        for _ in range(dim):
            for _try in range(64):
                g=rng.randint(1<<m)
                if g not in span:
                    break
            basis.append(g)
            span={s^ (g if (i>>j)&1 else 0) for s in list(span) for i in [1] for j in [0]}  # placeholder
            # rebuild span properly:
            newspan=set()
            for s in span: newspan.add(s)
            cur=set(newspan)
            for s in list(cur): newspan.add(s^g)
            span=newspan
        if len(span) < (1<<dim):
            continue  # not independent; skip
        # enumerate the affine subspace
        seen=set()
        pts=list(span)
        ok=True
        cset=set()
        for s in pts:
            cset.add(col[v0^s])
        if cset==full:
            return (t, v0, basis)
    return None

if __name__=="__main__":
    for n in (5,6):
        m,col,nc=colouring(n)
        ib=0
        while (1<<ib)<nc: ib+=1
        print(f"n={n}: m={m}, classes={nc}, info bound={ib}")
        res=affine_cover_search(n,m,col,nc,ib,tries=60000,seed=7)
        if res:
            t,v0,basis=res
            print(f"  FOUND affine {ib}-dim cover after {t} tries  "
                  f"(offset={v0}, basis={basis})  -> linear compression is TIGHT")
        else:
            print(f"  no affine {ib}-dim cover found in budget; "
                  f"(axis-aligned face needed {'4' if n==5 else '7'})")
