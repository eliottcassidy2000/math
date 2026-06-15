#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FKN on the tiling cube: the score vector is the level-<=1 (affine) Fourier part of a
tournament; the genuine cyclic content (c3, OCF) lives at level >=2. The Friedgut-
Kalai-Naor "dictator vs spread" hierarchy = the project's "transitive baseline vs
cyclic deviation" = Arrow's "rational vs Condorcet-cyclic".  kind-pasteur-2026-06-15-S6.

Setup (matches CLAUDE.md / tournament-tiling-explorer EXACTLY):
  vertices 1..n; base Hamiltonian path n -> n-1 -> ... -> 1 (higher beats next-lower).
  base-path arcs (k,k-1), k=2..n: FIXED (not tiles).
  TILES = non-consecutive arcs (x,y), x >= y+2; m = C(n-1,2) of them.
  zero tiling = transitive tournament (i beats j iff i>j), H=1.
  tile bit b=0: x beats y (transitive orientation).  b=1: y beats x (reversed).
We treat each tournament invariant as f:{0,1}^m -> R, take the Walsh-Hadamard
transform on {+-1}^m, and report the Fourier weight by level (degree).
"""
import sys, itertools
from math import comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
import numpy as np

# ---------- tiling model ----------
def tiles(n):
    """explorer enumeration order: for y=1..n-2: for x=n down to y+2: tile (x,y)."""
    T=[]
    for y in range(1,n-1):
        for x in range(n, y+1, -1):
            T.append((x,y))
    return T

def adj_from_bits(n, T, bits):
    """A[i][j]=1 iff i beats j (0-indexed vertices 0..n-1 == labels 1..n)."""
    A=np.zeros((n,n),dtype=int)
    # base path arcs: label k beats k-1  (k=2..n) -> index (k-1) beats (k-2)
    for k in range(2,n+1):
        A[k-1][k-2]=1
    # transitive closure of NON-tile, NON-base pairs is fixed by transitivity:
    # for any pair (x,y) x>y that is NOT a tile and NOT base it's still x beats y.
    # Simpler: start fully transitive (i>j -> i beats j) then apply tiles+base override.
    A[:]=0
    for i in range(n):
        for j in range(n):
            if i>j: A[i][j]=1     # transitive default: higher index beats lower
    # apply tile flips
    for (x,y),b in zip(T,bits):
        if b:  # reverse tile arc x,y  (x>y): now y beats x
            A[x-1][y-1]=0; A[y-1][x-1]=1
    return A

def scores(A): return A.sum(axis=1)           # out-degrees
def count_c3(A):
    n=len(A); c=0
    for i,j,k in itertools.combinations(range(n),3):
        # cyclic 3-set iff the 3 arcs form a directed cycle (either orientation)
        e=A[i][j]+A[j][k]+A[k][i]
        f=A[j][i]+A[k][j]+A[i][k]
        if e==3 or f==3: c+=1
    return c
def count_H(A):
    """# directed Hamiltonian paths (any start/end)."""
    n=len(A); cnt=0
    for perm in itertools.permutations(range(n)):
        if all(A[perm[t]][perm[t+1]] for t in range(n-1)): cnt+=1
    return cnt
def iso_key(A):
    """canonical form over relabelings (small n only)."""
    n=len(A); best=None
    for perm in itertools.permutations(range(n)):
        bits=tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or bits<best: best=bits
    return best

# ---------- Fourier (Walsh-Hadamard) ----------
def fwht(vec):
    """in-place-style WHT; returns Fourier coeffs hat f(S) with x in {+-1}, normalized
       so that f(x)=sum_S hat f(S) prod_{i in S} x_i (i.e. divide by 2^m)."""
    a=vec.astype(float).copy(); m=int(round(np.log2(len(a)))); h=1
    while h<len(a):
        for i in range(0,len(a),h*2):
            for j in range(i,i+h):
                x=a[j]; y=a[j+h]; a[j]=x+y; a[j+h]=x-y
        h*=2
    return a/len(a)

def level_weights(coeffs, m):
    """sum of squared Fourier coeffs grouped by |S| (Hamming weight of index)."""
    w=np.zeros(m+1)
    for idx,c in enumerate(coeffs):
        w[bin(idx).count('1')] += c*c
    return w

def analyze(n, verbose=True):
    T=tiles(n); m=len(T); N=1<<m
    assert m==comb(n-1,2)
    # bit ordering: index idx's bit t (LSB? ) — use idx bit t = (idx>>t)&1 for tile t.
    # fwht uses standard ordering where coordinate t corresponds to bit t. Build arrays.
    sc=np.zeros((N,n)); c3=np.zeros(N); Hv=np.zeros(N); istrans=np.zeros(N)
    isckey={}; keys=[]
    for idx in range(N):
        bits=[(idx>>t)&1 for t in range(m)]
        A=adj_from_bits(n,T,bits)
        s=scores(A); sc[idx]=s; c3[idx]=count_c3(A); Hv[idx]=count_H(A)
        istrans[idx]=1.0 if sorted(s)==list(range(n)) else 0.0
        k=iso_key(A); keys.append(k); isckey.setdefault(k,len(isckey))
    if verbose:
        print(f"\n===== n={n}, m={m} tiles, {N} tilings =====")
        print("zero-tiling check: H =",int(Hv[0]),"(want 1); scores=",sc[0].astype(int).tolist())
        print("#iso classes among tilings:",len(isckey),"(A000568:",[1,1,2,4,12,56,244,2704][n],")" if n<8 else ")")
    # --- score Fourier levels (each coordinate score, then aggregate) ---
    maxlvl_score=0
    for v in range(n):
        cf=fwht(sc[:,v]); lw=level_weights(cf,m)
        nz=[L for L in range(m+1) if lw[L]>1e-9]
        maxlvl_score=max(maxlvl_score, max(nz) if nz else 0)
    print(f"  SCORES: every coordinate's Fourier support max level = {maxlvl_score}  "
          f"(claim: <=1, i.e. AFFINE in tile bits)")
    # --- c3 ---
    cf=fwht(c3); lw=level_weights(cf,m); tot=lw.sum()
    print(f"  c3  level weights: "+", ".join(f"L{L}={lw[L]:.3f}" for L in range(min(m,5)+1) if lw[L]>1e-9))
    print(f"      c3 mean (L0 coeff) = {cf[0]:.4f}  (random-tiling E[c3]=C(n,3)/4={comb(n,3)/4:.4f})")
    print(f"      c3 FKN-defect (mass above L1)/(mass above L0) = {lw[2:].sum()/(tot-lw[0]):.4f}")
    # --- H ---
    cf=fwht(Hv); lw=level_weights(cf,m); tot=lw.sum()
    aboveL0=tot-lw[0]
    print(f"  H   level weights: "+", ".join(f"L{L}={lw[L]:.3f}" for L in range(min(m,6)+1) if lw[L]>1e-9))
    print(f"      H mean = {cf[0]:.4f} ; H level-1 weight / (var) = {lw[1]/aboveL0:.4f} "
          f"(FKN: low => H not score-determined)")
    print(f"      H FKN-defect (mass >=L2)/(var) = {lw[2:].sum()/aboveL0:.4f}")
    # --- transitive indicator ---
    cf=fwht(istrans); lw=level_weights(cf,m); tot=lw.sum()
    print(f"  [transitive] mass: L0={lw[0]:.4f} L1={lw[1]:.4f} >=L2={lw[2:].sum():.4f}  "
          f"(point-like => spread across ALL levels)")
    return T,m

if __name__=="__main__":
    for n in (4,5,6):
        analyze(n)
    print("\n--- n=7 (scores & c3 only; H/iso too slow) ---")
    # quick level check at n=7
    n=7; T=tiles(n); m=len(T); N=1<<m
    c3=np.zeros(N); sc0=np.zeros(N)
    for idx in range(N):
        bits=[(idx>>t)&1 for t in range(m)]
        A=adj_from_bits(n,T,bits); c3[idx]=count_c3(A); sc0[idx]=scores(A)[0]
    cf=fwht(sc0); lw=level_weights(cf,m)
    print("n=7 score[0] max Fourier level:", max([L for L in range(m+1) if lw[L]>1e-9]))
    cf=fwht(c3); lw=level_weights(cf,m)
    print("n=7 c3 level weights:", ", ".join(f"L{L}={lw[L]:.2f}" for L in range(6) if lw[L]>1e-9))
