#!/usr/bin/env python3
"""
quarter_tiling_complement_H_klein.py  --  klein-2026-07-01-S78

THREE targets:
 (Q) Does the tiling cube fold to a QUARTER? The half = sigma-fixed (grid-sym=complement-mirror) blue
     subspace W. A SECOND commuting involution flip=+1 (flip all tiles) folds W (and the whole cube) again.
     Verify: sigma = complement (T(sigma t) ~= T(t)^op), sigma^2=flip^2=id, they COMMUTE => <sigma,flip> =
     (Z/2)^2, and tilings/<sigma,flip> = the QUARTER. The blue LINES = W/flip = the sigma-fixed part of the
     quarter.
 (H) COMPLEMENT-H-PAIRING: H(T)=#Hamiltonian paths satisfies H(T)=H(T^op); complement pairs classes of
     equal H; each merged node carries a single H. Report H over categories.
 (R) NEXT TARGET: the BLUE FLIP-RANK rho_W = min subcube in the half-address (W-basis, one bit per
     sigma-orbit) coordinates hitting all SC classes; compare to rho_SC (full cube) and ceil(log2 #SC) --
     a covering bound intrinsic to the odd-cluster partition of the linear blue subspace.
"""
import itertools, math
from collections import Counter, defaultdict

def build(n):
    pairs=[(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]; pidx={p:k for k,p in enumerate(pairs)}
    tiles=[(x,y) for x in range(1,n+1) for y in range(1,x) if x-y>=2]; tidx={t:k for k,t in enumerate(tiles)}
    perms=list(itertools.permutations(range(1,n+1)))
    sigma=[tidx[(n+1-y,n+1-x)] for (x,y) in tiles]
    return pairs,pidx,tiles,tidx,perms,sigma
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
            u,w=(i,j) if ((mask>>k)&1) else (j,i); a,b=pi[u-1],pi[w-1]
            if a<b: v|=1<<pidx[(a,b)]
        if best is None or v<best: best=v
    return best
def opp(mask,pairs,pidx):
    v=0
    for k,(i,j) in enumerate(pairs):
        if not((mask>>k)&1): v|=1<<pidx[(i,j)]
    return v
def apply_sigma(tv,sigma,m):
    v=0
    for b in range(m):
        if (tv>>sigma[b])&1: v|=1<<b
    return v
def ham(mask,n,pairs):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k,(i,j) in enumerate(pairs):
        if (mask>>k)&1: A[i][j]=1
        else: A[j][i]=1
    return sum(1 for p in itertools.permutations(range(1,n+1)) if all(A[p[t]][p[t+1]] for t in range(n-1)))

def flip_rank_bits(points, cls, w, target):
    """min k: subcube in w-bit coords (fix w-k, free k) covering all classes in target (points: list of (bits->class))."""
    tset=set(target); nc=len(tset); lb=math.ceil(math.log2(nc)) if nc>1 else 0
    val={}
    for b,c in points: val[b]=c
    for k in range(lb,w+1):
        for free in itertools.combinations(range(w),k):
            rest=[e for e in range(w) if e not in free]
            for fa in range(1<<len(rest)):
                fb=0
                for bb,e in enumerate(rest):
                    if (fa>>bb)&1: fb|=1<<e
                seen=set(); good=True
                for a in range(1<<k):
                    bits=fb
                    for bb in range(k):
                        if (a>>bb)&1: bits|=1<<free[bb]
                    c=val.get(bits)
                    if c in tset:
                        seen.add(c)
                        if len(seen)==nc: break
                if len(seen)==nc: return k,lb
    return None,lb

if __name__=="__main__":
    for n in [4,5,6]:
        pairs,pidx,tiles,tidx,perms,sigma=build(n); m=len(tiles); full=(1<<m)-1
        cls=[0]*(1<<m); is_sc={}; cmask={}
        for tv in range(1<<m):
            mk=tmask(tv,n,tiles,pidx); c=canon(mk,n,pairs,pidx,perms); co=canon(opp(mk,pairs,pidx),n,pairs,pidx,perms)
            key=min(c,co); cls[tv]=key; is_sc[key]=(c==co); cmask[key]=c
        # (Q) sigma = complement?
        sig_is_comp=all(canon(tmask(apply_sigma(tv,sigma,m),n,tiles,pidx),n,pairs,pidx,perms)==canon(opp(tmask(tv,n,tiles,pidx),pairs,pidx),n,pairs,pidx,perms) for tv in range(min(1<<m,512)))
        # commute + involutions
        comm=all(apply_sigma(tv^full,sigma,m)==(apply_sigma(tv,sigma,m)^full) for tv in range(min(1<<m,512)))
        # orbits under <sigma,flip>
        seen=[False]*(1<<m); orbits=Counter()
        for tv in range(1<<m):
            if seen[tv]: continue
            orb={tv, apply_sigma(tv,sigma,m), tv^full, apply_sigma(tv,sigma,m)^full}
            for x in orb: seen[x]=True
            orbits[len(orb)]+=1
        half=sum((1 for _ in set(min(tv,apply_sigma(tv,sigma,m)) for tv in range(1<<m))))
        quarter=sum(orbits.values())
        print(f"\n===== n={n}: m={m}, 2^m={1<<m} tilings =====")
        print(f" (Q) sigma=complement: {sig_is_comp}; sigma,flip commute: {comm}; <sigma,flip>=(Z/2)^2 => FULL {1<<m} -> HALF(/sigma) {half} -> QUARTER(/<sigma,flip>) {quarter}; orbit sizes {dict(orbits)}")
        w=(m+(n-1)//2)//2
        print(f"     blue subspace W dim w={w}; blue LINES = W/flip = 2^(w-1)={1<<(w-1)} (= sigma-fixed part of the quarter)")
        # (H) complement-H-pairing
        Hpair=all(ham(cmask[k],n,pairs)==ham(opp(cmask[k],pairs,pidx),n,pairs) for k in list(is_sc)[:40])
        Hvals={k:ham(cmask[k],n,pairs) for k in list(is_sc)[:40]}
        print(f" (H) complement-H-pairing H(T)=H(T^op): {Hpair}; H values (sample) {sorted(set(Hvals.values()))} (all ODD, Redei)")
        # (R) blue flip-rank in half-address coords
        # half-address: pick representative tile per sigma-orbit
        orb_reps=[]; usedt=set()
        for b in range(m):
            o=min(b,sigma[b])
            if o not in usedt: usedt.add(o); orb_reps.append(o)
        # blue tiling <-> half-address: for a blue (sigma-fixed) tiling tv, bits = tv restricted to orb_reps
        blue=[tv for tv in range(1<<m) if apply_sigma(tv,sigma,m)==tv]
        pts=[]
        for tv in blue:
            bits=0
            for idx,b in enumerate(orb_reps):
                if (tv>>b)&1: bits|=1<<idx
            pts.append((bits, cls[tv]))
        SC=[k for k in is_sc if is_sc[k]]
        rW,lbW=flip_rank_bits(pts, cls, len(orb_reps), SC)
        rho_SC_full={4:1,5:4,6:6}.get(n,'?')
        print(f" (R) BLUE flip-rank rho_W (SC cover within W, w={len(orb_reps)} half-addr bits) = {rW} (log2 #SC={lbW}, #SC={len(SC)}); vs full-cube rho_SC={rho_SC_full}; excess rho_W-lb={None if rW is None else rW-lbW}")
