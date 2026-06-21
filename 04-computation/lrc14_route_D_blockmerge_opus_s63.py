#!/usr/bin/env python3
"""
lrc14_route_D_blockmerge_opus_s63.py   (opus-2026-06-20-S63)

ROUTE D, part 3.  Part 2 found the OBSTRUCTION: the Ly-metagraph has a few
LOCAL MAXIMA besides consec, and they are all "UNIONS OF CONSECUTIVE BLOCKS"
(e.g. k=8: {0,1,2,3} U {7,8,9,10} = two coherent blocks; k=9: one block + a far
singleton).  Local single-offset moves get trapped at these multi-block configs.

This reframes Route D in the RIGHT category.  The metagraph splits into TWO
scales:
   - WITHIN a block-pattern: single-offset moves close gaps INSIDE a block -- and
     these DO climb (each block wants to be solid).  -> reduces any shape to a
     "BLOCK SHAPE": a disjoint union of solid consecutive blocks.
   - BETWEEN block-patterns: the coarse metagraph whose VERTICES are block
     PARTITIONS (compositions of k into block sizes with gaps), and whose EDGE
     is a BLOCK-MERGE (slide two blocks together, or merge their sizes).
     consec = the single-block partition (k).  CLAIM to test: block-merge is
     Ly-MONOTONE on the block metagraph, with the single block as unique sink.

So the two-scale structure mirrors the project's "single coherent block is the
extremum" (THM-557).  We test, EXACTLY:

  (1) Are ALL local maxima of the full single-offset Ly-metagraph BLOCK SHAPES
      (disjoint solid consecutive blocks)?  Enumerate every local max.
  (2) On the BLOCK-SHAPE sub-metagraph (vertices = block compositions, edges =
      merge/slide two adjacent blocks closer), is Ly MONOTONE up to the single
      block?  i.e. does merging two blocks ALWAYS raise Ly?
  (3) The decisive reduction: if (1)+(2) hold, then
         max over all shapes  =  max over block shapes  =  single block = consec.
      Report whether this two-step gives a clean monotone ascent.

  (4) MERGE GEOMETRY: when two solid blocks of sizes a,b sit at gap d, slide them
      to gap d-1.  Track Delta Ly as a function of (a,b,d).  Is Delta Ly>=0 for
      every gap-closing slide between two solid blocks?  (the atomic merge lemma)

python3 stdlib only; exact Fractions.
"""
import sys, itertools, functools
from math import comb, gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

@functools.lru_cache(maxsize=None)
def pdist(Etuple):
    E=sorted(set(Etuple))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            for t in (F(j,7),F(j+1,7)):
                for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1)
    P=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        occ=[False]*7
        for e in E:
            occ[int(((e*xm)%1)*7)]=True
        N=sum(1 for j in range(1,7) if not occ[j])
        P[N]+=x1-x0
    return tuple(P)

def yvec(k):
    if k==8: g=lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40)
    elif k in (9,10): g=lambda t:F(-(t-2)*(t-3)*(t-6),36)
    else: g=lambda t:F((t-3)*(t-4),12)
    return [g(t) for t in range(7)]

@functools.lru_cache(maxsize=None)
def Ly(Etuple,k):
    g=yvec(k); P=pdist(Etuple)
    return sum(g[t]*P[t] for t in range(7))

def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=reduce(gcd,E,0)
    if g>1: E=[e//g for e in E]
    return tuple(E)
def span(E): return max(E)-min(E)

def all_compress_moves(E):
    E=sorted(set(E)); out=set()
    for i,e in enumerate(E):
        rest=[x for x in E if x!=e]
        for nv in range(min(rest),max(rest)+1):
            if nv in rest: continue
            cand=norm(rest+[nv])
            if len(cand)!=len(E) or span(cand)>span(E): continue
            out.add(cand)
    return list(out)

def all_shapes(k, max_span):
    seen=set()
    for rest in itertools.combinations(range(1,max_span+1),k-1):
        E=norm([0]+list(rest))
        if span(E)<=max_span and E not in seen:
            seen.add(E); yield E

def is_block_shape(E):
    """True iff E is a disjoint union of solid consecutive blocks (any gaps>=1
       between blocks).  Equivalently: every maximal run is consecutive -- which
       is automatic; a block shape is ANY set.  We mean: the BLOCK STRUCTURE,
       i.e. partition into maximal runs of consecutive integers.  Return the
       tuple of (block_size, gap_after) describing the composition."""
    E=sorted(set(E)); blocks=[]; i=0
    while i<len(E):
        j=i
        while j+1<len(E) and E[j+1]==E[j]+1: j+=1
        size=j-i+1
        gap = (E[j+1]-E[j]) if j+1<len(E) else 0
        blocks.append((size,gap)); i=j+1
    return blocks

def block_signature(E):
    """coarse signature: (block sizes, inter-block gaps)."""
    bl=is_block_shape(E)
    return tuple(b[0] for b in bl), tuple(b[1] for b in bl[:-1])

# build a block shape from sizes and gaps
def make_block_shape(sizes, gaps):
    """sizes: list of block sizes; gaps: list of gaps BETWEEN blocks (len=len(sizes)-1).
       gap g means g empty slots between blocks (so next block starts g+1 after prev end)."""
    E=[]; pos=0
    for i,s in enumerate(sizes):
        for _ in range(s): E.append(pos); pos+=1
        if i<len(sizes)-1: pos += gaps[i]   # pos jumps by gap (already +1 from last block end implicit)
    return norm(E)

if __name__=="__main__":
    print("="*72)
    print("ROUTE D part 3: local maxima are BLOCK shapes; block-merge monotone?")
    print("="*72)
    for k in (8,9,10):
        ms = {8:11,9:12,10:13}[k]
        consec=tuple(range(k)); lyc=Ly(consec,k)
        shapes=list(all_shapes(k,ms))
        print(f"\n### k={k}, span<={ms}, {len(shapes)} shapes.  Ly(consec)={float(lyc):.6f} ###")

        # (1) enumerate all local maxima; check each is a block shape
        locmax=[]
        for E in shapes:
            ly=Ly(E,k)
            if all(Ly(nb,k)<=ly for nb in all_compress_moves(E)):
                locmax.append(E)
        print(f"(1) #local maxima = {len(locmax)}")
        for E in locmax:
            sizes,gaps=block_signature(E)
            allgap1 = all(g>=1 for g in gaps) if gaps else True
            print(f"    {E}  Ly={float(Ly(E,k)):.6f}  blocks={sizes} gaps={gaps}"
                  f"  {'= CONSEC (single block)' if E==consec else ''}")

        # (2) BLOCK-MERGE monotonicity: vertices = block compositions of k.
        #     Use canonical (no-gap-bigger-than-needed) block shapes: gaps all =1
        #     between blocks (touching shapes differ by gap; we test minimal gap=1).
        #     Edge: merge two adjacent blocks (slide gap 1 -> 0, fusing sizes).
        print(f"(2) block-merge metagraph (blocks separated by gap=1):")
        # compositions of k into parts (ordered) -- block size sequences
        def compositions(n):
            if n==0: yield (); return
            for first in range(1,n+1):
                for rest in compositions(n-first):
                    yield (first,)+rest
        comps=list(compositions(k))
        # for each composition with gap-1 separators, compute Ly; an edge fuses
        # two adjacent blocks (sizes a,b -> a+b), which is "merge".  Monotone?
        viol=0; total=0; details=[]
        for comp in comps:
            if len(comp)==1: continue
            gaps=[1]*(len(comp)-1)
            E=make_block_shape(list(comp),gaps)
            lyE=Ly(E,k)
            # merge each adjacent pair
            for i in range(len(comp)-1):
                merged=list(comp[:i])+[comp[i]+comp[i+1]]+list(comp[i+2:])
                Em=make_block_shape(merged,[1]*(len(merged)-1))
                lyM=Ly(Em,k)
                total+=1
                if lyM<lyE:
                    viol+=1
                    details.append((comp,merged,float(lyE),float(lyM)))
        print(f"    adjacent-merge moves: {total}, Ly-DECREASING (violations): {viol}")
        for d in details[:6]:
            print(f"      {d[0]} -> {d[1]}: Ly {d[2]:.6f} -> {d[3]:.6f}  DROP")
        # also: GAP-CLOSE (slide blocks closer, gap d->d-1, same sizes) monotone?
        print(f"(2b) gap-close slides (two blocks, sizes a,b, gap d -> d-1):")
        gviol=0; gtot=0; gdet=[]
        for a in range(1,k):
            b=k-a
            for d in range(1, ms):  # gap from 1 up
                E1=make_block_shape([a,b],[d])
                E0=make_block_shape([a,b],[d-1]) if d-1>=0 else None
                if E0 is None or span(E1)>ms: continue
                ly1=Ly(E1,k); ly0=Ly(E0,k)
                gtot+=1
                if ly0<ly1:
                    gviol+=1; gdet.append((a,b,d,float(ly1),float(ly0)))
        print(f"    two-block gap-close: {gtot} slides, Ly-DECREASING: {gviol}")
        for d in gdet[:6]:
            print(f"      sizes({d[0]},{d[1]}) gap {d[2]}->{d[2]-1}: Ly {d[3]:.6f}->{d[4]:.6f} DROP")
    print("\nDONE part 3.")
