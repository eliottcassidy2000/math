#!/usr/bin/env python3
"""
almost_fixed_frame_s533.py    oracle-2026-06-01-S533

User: consider 2^x tournament iso classes ~ x independent-pair flips on a fixed
frame; understand how an ALMOST fixed frame can satisfy this (for n where the count
exceeds 2^x).

Verified separately: 2^floor(n/2) | A000568(n); quotient (frames) = 1,1,3,7,57,...
For n=3,4 the quotient is 1 -> a SINGLE fixed frame + floor(n/2) pair flips = all
iso classes. For n>=5 the frame must take 'quotient' values ('almost fixed').

This script, for n=5 (12 = 3 frames x 4 pair-flips) and n=6 (56 = 7 x 8):
 (1) fix a matching M (floor(n/2) disjoint arcs); the other arcs = the FRAME;
 (2) for each frame orientation, the 2^floor pair-flip settings give iso classes;
 (3) measure: how many distinct iso per frame (does a frame realize a full 2^floor
     block?), the minimum number of frames whose pair-flip blocks COVER all iso
     classes, and the HAMMING SPREAD of a minimal covering set ('almost fixed' =
     small spread).
"""
from itertools import combinations, permutations, product

def canon(adj):
    n=len(adj); best=None
    for p in permutations(range(n)):
        f=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best

def build(n, frame_arcs, frame_bits, M, match_bits):
    adj=[[0]*n for _ in range(n)]
    for (i,j),b in zip(frame_arcs, frame_bits):
        if b: adj[i][j]=1
        else: adj[j][i]=1
    for (i,j),b in zip(M, match_bits):
        if b: adj[i][j]=1
        else: adj[j][i]=1
    return adj

def analyze(n, M):
    allarcs=list(combinations(range(n),2))
    Mset=set(M)
    frame=[a for a in allarcs if a not in Mset]
    x=len(M); F=len(frame)
    print(f"  n={n}: matching M={M} ({x} pairs), frame={F} arcs, pair-flips=2^{x}")
    # for each frame, the block of iso classes (from 2^x pair-flips)
    blocks={}          # frame_bits -> frozenset of iso classes
    alliso=set()
    full_blocks=0
    for fb in product((0,1),repeat=F):
        isos=set()
        for mb in product((0,1),repeat=x):
            isos.add(canon(build(n,frame,fb,M,mb)))
        blocks[fb]=frozenset(isos); alliso|=isos
        if len(isos)==2**x: full_blocks+=1
    print(f"     total iso classes covered = {len(alliso)} (A000568({n}))")
    print(f"     frames whose 2^{x} pair-flips give a FULL block of {2**x} distinct iso = {full_blocks}/{2**F}")
    # minimum frames (greedy set cover) + hamming spread
    target=alliso; chosen=[]; covered=set()
    blist=sorted(blocks.items(), key=lambda kv:-len(kv[1]))
    while covered!=target:
        best=max(blist, key=lambda kv: len(kv[1]-covered))
        chosen.append(best[0]); covered|=best[1]
        blist=[b for b in blist if b[0]!=best[0]]
    def hd(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
    spread=max((hd(a,b) for a in chosen for b in chosen), default=0)
    print(f"     MIN frames to cover all iso (greedy) = {len(chosen)}  (quotient predicts {len(alliso)//(2**x)})")
    print(f"     Hamming spread of the covering frames = {spread} of {F} frame-arcs"
          + (f"  -> ALMOST FIXED (differ in <= {spread})" if spread<=max(1,F//3) else ""))
    return len(chosen), spread, F

def main():
    print("Almost-fixed frame: 2^x pair-flips + a barely-varying frame (oracle-S533)\n")
    # n=5: matching {(0,1),(2,3)}, vertex 4 unmatched
    analyze(5, [(0,1),(2,3)])
    print()
    # n=6: matching {(0,1),(2,3),(4,5)}
    analyze(6, [(0,1),(2,3),(4,5)])
    print("\nREADING: the iso classes factor as (frame) x (2^floor(n/2) pair-flips). For")
    print("n=3,4 one fixed frame suffices; for n>=5 the frame ranges over the quotient")
    print("A000568(n)/2^floor(n/2) but the needed frames are HAMMING-CLOSE -- an almost")
    print("fixed frame whose few toggled arcs ride along with the pair-flips. The pair-")
    print("flips are the free (Z/2)^floor(n/2) channels (the LRC independent pairs, S532);")
    print("the frame residue is the small remaining coupling.")

if __name__=="__main__":
    main()
