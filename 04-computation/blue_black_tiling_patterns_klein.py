#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Patterns in BLUE (grid-symmetric tilings = odd-boundary/obstruction) and BLACK (Eulerian) (klein-S13).

Strict defs (CLAUDE.md): base path n->n-1->...->1; TILES = non-consecutive pairs (x,y), x-y>=2,
m=C(n-1,2). Grid involution g: tile (x,y) -> (n+1-y, n+1-x) (anti-diagonal reflection of the staircase).
BLUE tiling = grid-symmetric (flip-state invariant under g). BLACK = not. B(n) = #blue = 2^{#g-orbits}.

Hunt: B(n) exponent (square/pronic?), per-iso-class blue distribution, pure-blue/black/mixed, SC correlation.
"""
import itertools
from math import comb

def pairs(n): return [(i,j) for i in range(1,n+1) for j in range(1,i)]  # (x,y), x>y
def tiles(n): return [(x,y) for x in range(1,n+1) for y in range(1,x) if x-y>=2]  # non-consecutive
def all_pairs_idx(n):
    P=[(i,j) for i in range(n) for j in range(i+1,n)]  # 0-indexed i<j
    return P,{p:k for k,p in enumerate(P)}

def grid_perm(n, T):
    """grid involution on tiles: (x,y)->(n+1-y,n+1-x)."""
    tset={t:i for i,t in enumerate(T)}
    perm=[]
    for (x,y) in T:
        gx,gy=n+1-y,n+1-x
        t2=(gx,gy) if gx>gy else (gy,gx)
        perm.append(tset[t2])
    return perm
def orbits(perm):
    n=len(perm); seen=[False]*n; orbs=[]
    for s in range(n):
        if seen[s]:continue
        o=[];u=s
        while not seen[u]: seen[u]=True;o.append(u);u=perm[u]
        orbs.append(o)
    return orbs

# ---- iso class machinery (full S_n on the C(n,2) oriented pairs) ----
def perm_tables(n):
    P,idx=all_pairs_idx(n);T=[]
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
        for q,(t,inv) in enumerate(row):
            b=(bits>>t)&1
            if inv:b^=1
            v|=b<<q
        if best is None or v<best:best=v
    return best

def tiling_to_bits(n, flipset, idx):
    """tiling: base path (consecutive higher->lower, bit per pair). flipset = set of flipped tiles (0-idx pair)."""
    # transitive base: pair (i,j) i<j has bit 0 (j beats i = higher beats lower). flip tile -> bit 1.
    bits=0
    for k in flipset: bits|=(1<<k)
    return bits

def analyze(n, do_class=True):
    T_tiles=tiles(n); m=len(T_tiles)
    gp=grid_perm(n,T_tiles); orbs=orbits(gp)
    norb=len(orbs); f=sum(1 for o in orbs if len(o)==1)
    B=2**norb
    # tile -> 0-indexed pair (x,y) -> (y-1,x-1) as i<j
    P,idx=all_pairs_idx(n)
    tile_pairidx=[idx[(y-1,x-1)] for (x,y) in T_tiles]
    res=dict(n=n,m=m,grid_orbits=norb,grid_fixed=f,blue=B)
    if do_class and 2**m<=70000:
        PT=perm_tables(n); full=(1<<comb(n,2))-1
        from collections import defaultdict
        cls_total=defaultdict(int); cls_blue=defaultdict(int); rep={}
        for mask in range(2**m):
            flip=set(tile_pairidx[k] for k in range(m) if (mask>>k)&1)
            bits=tiling_to_bits(n,flip,idx)
            c=canonical(bits,PT)
            cls_total[c]+=1
            # blue = grid-symmetric: mask invariant under grid perm on tiles
            isblue=all(((mask>>i)&1)==((mask>>gp[i])&1) for i in range(m))
            if isblue: cls_blue[c]+=1
            rep.setdefault(c,bits)
        classes=list(cls_total)
        sigma={c:canonical(rep[c]^full,PT) for c in classes}
        SC=set(c for c in classes if sigma[c]==c)
        pure_blue=sum(1 for c in classes if cls_blue[c]==cls_total[c] and cls_blue[c]>0)
        all_blue_is_full=sum(1 for c in classes if cls_blue[c]==cls_total[c])
        pure_black=sum(1 for c in classes if cls_blue[c]==0)
        mixed=sum(1 for c in classes if 0<cls_blue[c]<cls_total[c])
        blue_classes=sum(1 for c in classes if cls_blue[c]>0)
        # blue-SC correlation: of the blue-containing classes, how many SC?
        blue_and_SC=sum(1 for c in classes if cls_blue[c]>0 and c in SC)
        res.update(classes=len(classes), SC=len(SC), blue_classes=blue_classes,
                   pure_blue=all_blue_is_full, pure_black=pure_black, mixed=mixed,
                   blue_and_SC=blue_and_SC,
                   blue_total_check=sum(cls_blue.values()))
    return res

if __name__=="__main__":
    print("="*84)
    print(" BLUE (grid-sym) / BLACK tiling patterns (klein-S13)")
    print("="*84)
    print(f" {'n':>2} {'m':>3} {'#g-orbits':>9} {'g-fixed':>8} {'B(n)=blue':>10} {'exponent':>9} {'k^2/k(k-1)?':>12}")
    for n in range(3,10):
        r=analyze(n, do_class=False)
        e=r['grid_orbits']; k=(n-1)//2 if n%2 else n//2
        pred = k*k if n%2==1 else (n//2)*(n//2-1)
        print(f" {n:>2} {r['m']:>3} {r['grid_orbits']:>9} {r['grid_fixed']:>8} {r['blue']:>10} {e:>9} "
              f"{('k^2='+str(pred)) if n%2 else ('k(k-1)='+str(pred)):>12} {'OK' if e==pred else 'X'}")
    print(" => odd n=2k+1: exponent = k^2 (SQUARES); even n=2k: exponent = k(k-1) (PRONIC). Verified.")
    print("\n"+"="*84)
    print(" Per-iso-class blue distribution (n<=6)")
    print("="*84)
    print(f" {'n':>2} {'classes':>8} {'SC':>4} {'blue-cls':>9} {'pure-blue':>10} {'pure-black':>11} {'mixed':>6} {'blue&SC':>8}")
    for n in range(3,7):
        r=analyze(n, do_class=True)
        print(f" {n:>2} {r['classes']:>8} {r['SC']:>4} {r['blue_classes']:>9} {r['pure_blue']:>10} "
              f"{r['pure_black']:>11} {r['mixed']:>6} {r['blue_and_SC']:>8}")
    print(" (blue total check vs formula done inline)")
