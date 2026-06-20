#!/usr/bin/env python3
"""
lrc_halftiling_foundations_macmini_0620s4.py  (mac-mini-2026-06-20-S4)

Verify the HALF-TILING framework (user's S4 prompt) from first principles.

TILING MODEL: tournament on {1..n}, fixed base Ham path n->n-1->...->1 (arc k->k-1).
 - base-path arcs (k,k-1): NOT tiles.  TILES = arcs (x,y), x>y, x-y>=2 (m=C(n-1,2) of them).
 - a TILING = orientation of every tile (fwd x->y 'with' descending path, or bwd y->x).
 - reflection sigma over y=x : tile (x,y) -> (n+1-y, n+1-x)  [the grid map, CLAUDE.md].

CLAIMS TO VERIFY:
 (1) sigma is an involution on tiles; fixed tiles = {x+y=n+1} ; #fixed = floor((n-1)/2).
 (2) half-region size = (m + #fixed)/2 = floor((n-1)^2/4).
 (3) sigma on the labeled tournament = relabel phi(i)=n+1-i COMPOSED with reverse-all-arcs (T^op),
     i.e. sigma(T) = phi(T^op); this preserves the base path. (The 'mirror = reverse all arcs'.)
 (4) even/odd half-size recursions: even h(n)=2h(n-1)-h(n-2); odd h(n)=2h(n-1)-2h(n-3)+h(n-4).
"""
import itertools, sys
sys.stdout.reconfigure(line_buffering=True)

def tiles(n):  # tile list (x,y), x>y, x-y>=2
    return [(x,y) for y in range(1,n-1) for x in range(n, y+1, -1) if x-y>=2]

def sigma(t,n):  # reflection over y=x
    x,y=t; return (n+1-y, n+1-x)

print(f"{'n':>3}{'m=C(n-1,2)':>11}{'#fixed':>8}{'half=(m+f)/2':>13}{'floor((n-1)^2/4)':>18}{'sigma invol?':>13}")
for n in range(2,13):
    T=tiles(n); m=len(T)
    fixed=[t for t in T if sigma(t,n)==t]
    # check sigma maps tiles to tiles and is an involution
    ok_inv = all(sigma(t,n) in T for t in T) and all(sigma(sigma(t,n),n)==t for t in T)
    half=(m+len(fixed))//2
    qs=(n-1)**2//4
    print(f"{n:>3}{m:>11}{len(fixed):>8}{half:>13}{qs:>18}{str(ok_inv and half==qs and len(fixed)==(n-1)//2):>13}")

# (3) sigma on tournament = phi . T^op.  Represent tournament by arc-direction dict on all pairs.
def base_arc(i,j):  # base path orientation for consecutive; for tiles default fwd (x->y)
    # returns the 'forward' head for pair (i,j) i<j: forward means j->i (descending). head=i.
    return (j,i)  # j->i
def tournament_from_tiling(n, tiling):
    """tiling: dict tile-> +1 (fwd j->i) or -1 (bwd i->j). base path always fwd. return arc set."""
    arcs=set()
    for j in range(2,n+1):
        for i in range(1,j):
            if j-i==1: arcs.add((j,i))           # base path k->k-1 (fwd, fixed)
            else:
                o=tiling[(j,i)]
                arcs.add((j,i) if o==1 else (i,j))
    return arcs
def reverse_all(arcs): return set((b,a) for (a,b) in arcs)
def relabel(arcs,n): return set((n+1-a, n+1-b) for (a,b) in arcs)
def tiling_reflect(n, tiling):
    """apply sigma to the tiling: new orientation of tile t = orientation of sigma(t), but
       reflection also flips fwd/bwd? Test both and see which gives phi.T^op."""
    out={}
    for t in tiling:
        out[t]=tiling[sigma(t,n)]   # plain relabel of orientation
    return out

import random
random.seed(1)
for n in [3,4,5,6]:
    T=tiles(n); good_plain=0; good_flip=0; total=0
    # sample tilings
    samples = list(itertools.product([1,-1],repeat=len(T))) if len(T)<=10 else [tuple(random.choice([1,-1]) for _ in T) for _ in range(200)]
    for vals in samples:
        total+=1
        tiling=dict(zip(T,vals))
        arcs=tournament_from_tiling(n,tiling)
        target=relabel(reverse_all(arcs),n)   # phi(T^op)
        # candidate A: reflected tiling, plain orientation copy
        refl_plain=dict((t, tiling[sigma(t,n)]) for t in T)
        arcs_plain=tournament_from_tiling(n,refl_plain)
        # candidate B: reflected tiling, flipped orientation
        refl_flip=dict((t, -tiling[sigma(t,n)]) for t in T)
        arcs_flip=tournament_from_tiling(n,refl_flip)
        if arcs_plain==target: good_plain+=1
        if arcs_flip==target: good_flip+=1
    print(f"n={n}: of {total} tilings, reflect(plain)==phi(T^op): {good_plain}; reflect(flip)==phi(T^op): {good_flip}")
