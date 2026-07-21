#!/usr/bin/env python3
"""
LOCATE THE EDGE: the first n where poly-time invariants FAIL to determine H, and measure the
2-adic depth of H's poly-determinability there.
kind-pasteur-2026-07-21-S128c143.

At n<=6 the joint poly-tuple P=(score,specA,specS,disc,arb) DETERMINES the iso class (56 cells =
56 classes at n=6), so H is trivially poly-determined. The #P-ness of H is asymptotic: it needs
n>=7, where non-isomorphic tournaments can share P. We SAMPLE random labeled T_n, bucket by P, and
find buckets with >1 distinct H (edge witnesses). The gcd of within-bucket H-differences = M;
v2(M) = number of poly-determined LOW BITS of H at the edge. Redei forces v2(M) >= 1 (all H odd).

Everything EXACT integer (Faddeev-LeVerrier char polys, Bareiss determinant) so hashing is safe.
"""
import sys, math
from functools import reduce
from collections import defaultdict
import random

def rand_A(n, rng):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.getrandbits(1): A[i][j] = 1
            else: A[j][i] = 1
    return A

def matmul(A, B, n):
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def charpoly_int(M, n):
    """Faddeev-LeVerrier, exact python ints. Returns (c0=1,...,cn)."""
    c = [1]; Mk = [[1 if i==j else 0 for j in range(n)] for i in range(n)]
    for k in range(1, n+1):
        Mk = matmul(M, Mk, n)
        tr = sum(Mk[i][i] for i in range(n))
        ck = -tr // k            # exact: k | tr at each step for integer matrices
        c.append(ck)
        for i in range(n): Mk[i][i] += ck
    return tuple(c)

def bareiss_det(M, n):
    """Exact integer determinant (Bareiss)."""
    M = [row[:] for row in M]; sign = 1; prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            sw = next((r for r in range(k+1, n) if M[r][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def score(A, n): return tuple(sorted(sum(A[i]) for i in range(n)))
def specA(A, n): return charpoly_int(A, n)
def Kmat(A, n): return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
def specS(A, n): return charpoly_int(Kmat(A, n), n)
def disc(A, n):
    cK = charpoly_int(Kmat(A, n), n)           # det(xI-K)
    val = sum(cK[k]*((-1)**(n-k)) for k in range(n+1))   # charpoly_K(-1) = det(-I-K)
    detIK = ((-1)**n) * val                    # det(I+K)
    from fractions import Fraction
    return Fraction(abs(detIK), 2**(n-1))
def arb(A, n):
    # in-degree Laplacian for out-trees rooted at 0: L = diag(indeg) - A^T ; delete row/col 0
    indeg = [sum(A[r][c] for r in range(n)) for c in range(n)]
    Lm = [[ (indeg[i] if i==j else 0) - A[j][i] for j in range(n)] for i in range(n)]
    minor = [row[1:] for row in Lm[1:]]
    return bareiss_det(minor, n-1)
def H_count(A, n):
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1<<n):
        for v in range(n):
            if not (mask>>v)&1 or dp[mask][v]==0: continue
            for w in range(n):
                if (mask>>w)&1 or not A[v][w]: continue
                dp[mask|(1<<w)][w] += dp[mask][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))

def gcd_list(xs):
    xs = [x for x in xs if x];
    return reduce(math.gcd, xs) if xs else 0
def v2(m):
    if m == 0: return math.inf
    k = 0
    while m % 2 == 0: m//=2; k+=1
    return k

def run(n, N, seed=1):
    rng = random.Random(seed)
    buckets = defaultdict(set)      # P -> set of H
    for _ in range(N):
        A = rand_A(n, rng)
        P = (score(A,n), specA(A,n), specS(A,n), disc(A,n), arb(A,n))
        buckets[P].add(H_count(A,n))
    split = {P:Hs for P,Hs in buckets.items() if len(Hs) > 1}
    diffs = []
    for Hs in split.values():
        s = sorted(Hs)
        for i in range(len(s)):
            for j in range(i+1,len(s)):
                diffs.append(s[j]-s[i])
    M = gcd_list(diffs)
    print(f"n={n}: samples={N}  distinct P-buckets={len(buckets)}  "
          f"H-SPLIT buckets (edge witnesses)={len(split)}")
    if split:
        M2 = 2**v2(M)
        print(f"   M = gcd(H-differences) = {M}   v2(M) = {v2(M)}   odd-part={M>>v2(M) if M else 0}")
        print(f"   => at the edge, H mod {M2} is poly-time (FORMULA); H mod {M2*2} is NOT poly-determined")
        ex = sorted(split.items(), key=lambda kv: -len(kv[1]))[:3]
        for P,Hs in ex:
            print(f"     witness bucket: H-values {sorted(Hs)}   (all odd? {all(h%2==1 for h in Hs)})")
    else:
        print(f"   NO edge witness found in {N} samples — poly-tuple still determines H at this sampling depth")
        print(f"   (detection floor: cospectral+coscore+codisc+coarb non-iso pairs rarer than ~1/{N})")
    return M, len(split)

if __name__ == "__main__":
    print("="*76); print("EDGE SEARCH — where do poly invariants first fail to determine H?"); print("="*76)
    run(6, 40000, seed=1)     # control: should find 0 (poly-tuple determines iso at n<=6)
    run(7, 300000, seed=2)
    run(7, 300000, seed=3)    # second seed to accumulate witnesses
    run(8, 120000, seed=4)
    print("DONE.")
