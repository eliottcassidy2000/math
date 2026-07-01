#!/usr/bin/env python3
"""
tjoin_sc_cover_obstruction_klein.py  --  klein-2026-07-01-S77

Does the T-JOIN boundary parity of the SC classes obstruct LOW-DIMENSIONAL (subcube) covers of them?

Setup: blue lines (d=m, grid-sym) form a graph on merged nodes where every SC node has ODD blue-degree =>
the SC nodes are the boundary T of a T-JOIN (|T|=#SC must be even -- the handshake). Grid-symmetric (blue)
tilings are the sigma-fixed set; sigma is a LINEAR involution on tile-coords, so the blue tilings form a
LINEAR SUBSPACE W. Question: since each SC class meets W in an ODD number of tilings (T-join parity), does
that obstruct covering all SC classes with a low-dim subcube (=> flip-rank excess concentrated on SC)?

Tests (n=4,5,6):
 (A) blue tilings W = linear subspace? dim = (m+fix)/2, fix=floor((n-1)/2). SC class blue-fiber |C∩W| ODD?
 (B) FLIP-RANK restricted to class subsets: rho_all vs rho_SC vs rho_NS -- is SC the bottleneck?
 (C) OBSTRUCTION: excess rho_SC - ceil(log2 #SC); and the T-join parity check on affine covers.
"""
import itertools, math
from collections import Counter, defaultdict

def build(n):
    pairs=[(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    pidx={p:k for k,p in enumerate(pairs)}
    tiles=[(x,y) for x in range(1,n+1) for y in range(1,x) if x-y>=2]
    perms=list(itertools.permutations(range(1,n+1)))
    return pairs,pidx,tiles,perms
def tmask(tv,n,tiles,pidx):
    A=[[0]*(n+1) for _ in range(n+1)]
    for a in range(2,n+1): A[a][a-1]=1
    for b,(x,y) in enumerate(tiles):
        if (tv>>b)&1: A[x][y]=1
        else: A[y][x]=1
    mk=0
    for (i,j) in [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]:
        if A[i][j]: mk|=1<<pidx[(i,j)]
    return mk
def canon(mask,n,pairs,pidx,perms):
    best=None
    for pi in perms:
        v=0
        for k,(i,j) in enumerate(pairs):
            u,w=(i,j) if ((mask>>k)&1) else (j,i)
            a,b=pi[u-1],pi[w-1]
            if a<b: v|=1<<pidx[(a,b)]
        if best is None or v<best: best=v
    return best
def opp(mask,pairs,pidx):
    v=0
    for k,(i,j) in enumerate(pairs):
        if not((mask>>k)&1): v|=1<<pidx[(i,j)]
    return v

def flip_rank(m, cls, target, cap=6_000_000):
    """min k: exists subcube (fix m-k, free k) whose 2^k completions cover all classes in `target` set."""
    tset=set(target); nc=len(tset); lb=math.ceil(math.log2(nc))
    for k in range(lb, m+1):
        if math.comb(m,k)*(1<<(m-k))>cap: return None,lb
        for free in itertools.combinations(range(m),k):
            rest=[e for e in range(m) if e not in free]
            for fa in range(1<<len(rest)):
                fb=0
                for b,e in enumerate(rest):
                    if (fa>>b)&1: fb|=1<<e
                seen=set()
                for a in range(1<<k):
                    tt=fb
                    for b in range(k):
                        if (a>>b)&1: tt|=1<<free[b]
                    c=cls[tt]
                    if c in tset:
                        seen.add(c)
                        if len(seen)==nc: break
                if len(seen)==nc: return k,lb
    return None,lb

if __name__=="__main__":
    for n in [4,5,6]:
        pairs,pidx,tiles,perms=build(n); m=len(tiles); fix=(n-1)//2
        # sigma on tiles
        tidx={t:k for k,t in enumerate(tiles)}; sigma=[tidx[(n+1-y,n+1-x)] for (x,y) in tiles]
        cls=[0]*(1<<m); is_sc={}
        for tv in range(1<<m):
            mk=tmask(tv,n,tiles,pidx); c=canon(mk,n,pairs,pidx,perms); co=canon(opp(mk,pairs,pidx),n,pairs,pidx,perms)
            key=min(c,co); cls[tv]=key; is_sc[key]=(c==co)
        keys=list(set(cls)); SC=[k for k in keys if is_sc[k]]; NS=[k for k in keys if not is_sc[k]]
        gs=lambda tv: all(((tv>>b)&1)==((tv>>sigma[b])&1) for b in range(m))
        blue=[tv for tv in range(1<<m) if gs(tv)]
        # (A) W linear? closed under XOR + contains 0
        Wset=set(blue); lin = (0 in Wset) and all((a^b) in Wset for a in blue[:64] for b in blue[:64])
        wdim=len(blue).bit_length()-1
        # SC blue-fiber parity
        bluefib=Counter(cls[tv] for tv in blue)
        sc_blue_odd=all(bluefib[k]%2==1 for k in SC); ns_blue_zero=all(bluefib.get(k,0)==0 for k in NS)
        print(f"\n===== n={n}: m={m}, classes {len(keys)} (SC {len(SC)}, NS {len(NS)}) =====")
        print(f" (A) blue W: |W|={len(blue)}=2^{wdim} (dim {(m+fix)//2} expected); linear(sample)={lin}; each SC blue-fiber ODD={sc_blue_odd}; NS blue-fiber 0={ns_blue_zero}")
        # (B) flip-ranks
        ra,lba=flip_rank(m,cls,keys); rs,lbs=flip_rank(m,cls,SC); rn,lbn=flip_rank(m,cls,NS)
        print(f" (B) flip-rank: rho_all={ra} (lb {lba}); rho_SC={rs} (lb {lbs}, #SC={len(SC)}); rho_NS={rn} (lb {lbn}, #NS={len(NS)})")
        # (C) obstruction: excesses
        def exc(r,lb): return None if r is None else r-lb
        print(f" (C) covering-excess: SC {exc(rs,lbs)}, NS {exc(rn,lbn)}, all {exc(ra,lba)}  => bottleneck is {'SC' if (rs or 0)>=(rn or 0) and (rs or 0)>=(ra or 0) else 'NS/all'}")
        print(f"     T-join: SC nodes = odd-boundary of blue graph (|T|=#SC={len(SC)}, even? {len(SC)%2==0}); SC blue-fibers odd => low-dim affine cover of SC must resolve an odd-parity partition of W")
