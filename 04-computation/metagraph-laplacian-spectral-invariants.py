#!/usr/bin/env python3
"""
metagraph-laplacian-spectral-invariants.py   (klein-2026-06-29-S3)

New small invariants of the arc-flip (dominance-reversal) tournament metagraph,
all derived from THM-587's signed cycle index mult(k) -- so NO 2^{C(n,2)}
enumeration, reaching large n cheaply. Plus a verification (n<=6) that the
slow/Fiedler mode is the CYCLICITY (3-cycle count).

Background (THM-584/587): the arc-flip metagraph A (single dominance reversal =
one Q_d edge; d=C(n,2)) is d-regular (weighted) with eigenvalues d-2k of
multiplicity mult(k) = [x^k] P_n(x), P_n = signed cycle index.

NEW DEFINITIONS (this script):
  * Laplacian L = dI - A. Spectrum = {2k : k=0..d} with multiplicity mult(k).
  * Algebraic connectivity a(n) = smallest nonzero Laplacian eigenvalue
        = 2 * min{k>=1 : mult(k)>0}.  CLAIM: = 4 for all n>=3 (mult(1)=0, mult(2)=1).
  * Dominance-reversal walk P = A/d. Spectral gap = 1 - lambda_2 = (smallest
        Laplacian eigenvalue)/d = 4/d. CLAIM exact (sharpens opus-S268's c/n conj
        which was for a DIFFERENT, merged-wiggly normalization).
  * Cyclicity Fiedler mode: eigenvector at eigenvalue d-4 (level 2) ∝ 3-cycle count.
  * Weighted complexity (spanning trees) kappa(n) = (1/N) prod_{k>=1} (2k)^{mult(k)}.
  * Neutral-flip trace nu(n) = tr(A) = sum_k mult(k)(d-2k) = total self-loop
        (silent-mutation) weight summed over classes.
  * Heat trace Theta(n;b) = sum_k mult(k) e^{-2k b}; spectral zeta zeta(n;s)
        = sum_{k>=1} mult(k)/(2k)^s.
"""
import itertools
from math import comb, factorial
from fractions import Fraction
import numpy as np

# ---------- signed cycle index (THM-587) ----------
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] += ai*bj
    return r

def signed_pair_cycles(perm, n):
    pid = {}; pairs = []
    for i in range(n):
        for j in range(i+1, n):
            pid[(i, j)] = len(pairs); pairs.append((i, j))
    def step(p):
        i, j = p; a, b = perm[i], perm[j]
        return ((a, b), +1) if a < b else ((b, a), -1)
    seen = [False]*len(pairs); cyc = []
    for s0 in range(len(pairs)):
        if seen[s0]: continue
        p = pairs[s0]; s = 1; L = 0
        while not seen[pid[p]]:
            seen[pid[p]] = True; p, st = step(p); s *= st; L += 1
        cyc.append((L, s))
    return cyc

def signed_cycle_index(n):
    d = comb(n, 2); total = [0]*(d+1)
    for perm in itertools.permutations(range(n)):
        poly = [1]
        for (L, s) in signed_pair_cycles(perm, n):
            fac = [0]*(L+1); fac[0] = 1; fac[L] = s
            poly = pmul(poly, fac)
        for i, c in enumerate(poly):
            total[i] += c
    nf = factorial(n)
    return [c//nf for c in total]

# ---------- invariants from mult(k) ----------
def invariants(n):
    d = comb(n, 2); mult = signed_cycle_index(n)
    N = sum(mult)                                   # #classes = A000568
    lap = []                                        # Laplacian eigenvalues 2k with mult
    for k, m in enumerate(mult):
        lap += [2*k]*m
    lap.sort()
    nonzero = [x for x in lap if x > 0]
    alg_conn = min(nonzero)                         # algebraic connectivity
    gap_walk = Fraction(alg_conn, d)                # P=A/d spectral gap
    # weighted spanning-tree complexity = (1/N) prod nonzero Laplacian eigenvalues
    kappa = Fraction(1, N)
    for x in nonzero:
        kappa *= x
    nu = sum(m*(d-2*k) for k, m in enumerate(mult))  # tr(A) self-loop weight
    k1 = mult[1] if d >= 1 else 0
    k2 = mult[2] if d >= 2 else 0
    return dict(d=d, N=N, mult=mult, alg_conn=alg_conn, gap=gap_walk,
                kappa=kappa, nu=nu, mult1=k1, mult2=k2,
                heat_b1=sum(m*np.exp(-2*k*1.0) for k,m in enumerate(mult)),
                spectral_zeta_2=float(sum(Fraction(m,(2*k)**2) for k,m in enumerate(mult) if k>=1)))

# ---------- verify Fiedler mode = cyclicity (n<=6, build A) ----------
def pairs(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
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
def scores(n,bits):
    P=pairs(n); idx={p:k for k,p in enumerate(P)}; s=[0]*n
    for (i,j) in P:
        if (bits>>idx[(i,j)])&1: s[i]+=1
        else: s[j]+=1
    return s
def cyclicity(n,bits):
    s=scores(n,bits)
    return comb(n,3) - sum(comb(si,2) for si in s)   # # cyclic triangles

def verify_cyclicity_fiedler(n):
    d=comb(n,2); T=perm_tables(n)
    classes=[]; class_of={}
    for bits in range(2**d):
        c=canonical(bits,T); class_of[bits]=c
        if c==bits: classes.append(c)
    cidx={c:i for i,c in enumerate(classes)}; N=len(classes)
    A=np.zeros((N,N))
    for c in classes:
        for a in range(d):
            A[cidx[c]][cidx[class_of[c^(1<<a)]]]+=1
    w,V=np.linalg.eig(A)
    w=w.real; V=V.real
    # eigenvector(s) at eigenvalue d-4 (level 2)
    target=d-4
    idxs=[i for i in range(N) if abs(w[i]-target)<1e-6]
    cyc=np.array([cyclicity(n,c) for c in classes],dtype=float)
    cyc=cyc-cyc.mean()
    best=0.0
    for i in idxs:
        v=V[:,i]
        if np.linalg.norm(v)>0 and np.linalg.norm(cyc)>0:
            r=abs(np.dot(v,cyc)/(np.linalg.norm(v)*np.linalg.norm(cyc)))
            best=max(best,r)
    return target, len(idxs), best

if __name__=="__main__":
    print("="*78)
    print(" Arc-flip (dominance-reversal) metagraph: new spectral invariants (from THM-587)")
    print("="*78)
    print(f" {'n':>2} {'d':>3} {'N=A000568':>10} {'mult(1)':>7} {'mult(2)':>7} "
          f"{'alg.conn':>8} {'walk gap 4/d':>14} {'tr(A)=nu':>9} {'kappa(span.trees)':>20}")
    for n in range(3,10):
        I=invariants(n)
        print(f" {n:>2} {I['d']:>3} {I['N']:>10} {I['mult1']:>7} {I['mult2']:>7} "
              f"{I['alg_conn']:>8} {str(I['gap']):>14} {I['nu']:>9} {str(I['kappa']):>20}")
    print("\n CLAIM: mult(1)=0 and mult(2)=1 for all n -> algebraic connectivity = 4 for all n>=3.")
    print(" (Laplacian spectrum is exactly {2k : mult(k)>0}; gap of walk P=A/d is exactly 4/d.)")

    print("\n" + "="*78)
    print(" Verify the level-2 (eigenvalue d-4) mode IS the cyclicity (3-cycle count), n=4,5,6")
    print("="*78)
    for n in (4,5,6):
        tgt,multg,corr=verify_cyclicity_fiedler(n)
        print(f"   n={n}: eigenvalue d-4={tgt}, geom mult={multg}, "
              f"|corr(eigvec, 3-cycle count)| = {corr:.6f}")

    print("\n" + "="*78)
    print(" Heat trace Theta(b=1) and spectral zeta(s=2) (from mult(k)):")
    print("="*78)
    for n in range(3,10):
        I=invariants(n)
        print(f"   n={n}: Theta(1)={I['heat_b1']:.4f}   zeta(2)={I['spectral_zeta_2']:.5f}")
