#!/usr/bin/env python3
"""
transfer_matrix_F_connection.py — Can F(T,x) be expressed via the transfer matrix?

The transfer matrix W(x) for tournament T has:
  W[u][v] = x if u->v (forward), 1 if v->u (backward)
  W[u][u] = 0

Key properties from S44:
  W[u][v] * W[v][u] = x always
  tr(W) = 0
  tr(W^2) = n(n-1)x (universal)
  tr(W^3) depends on t3

F(T,x) = sum over all permutations P of prod_{i=0}^{n-2} W[P[i]][P[i+1]]
        = sum over Hamiltonian paths of weight product
        = permanent-like sum (but without the circular closure)

QUESTION: Is F(T,x) = per_path(W) where per_path is a "path permanent"?

The PERMANENT of W: per(W) = sum_{sigma in S_n} prod W[i][sigma(i)]
This counts weighted cycle covers, NOT Hamiltonian paths.

The "Hamiltonian path sum" is:
  H(W) = sum_{(v_0,...,v_{n-1}) permutation} prod_{i=0}^{n-2} W[v_i][v_{i+1}]

This is related to the permanent of the (n-1) x n matrix obtained by
removing one row of W, but it's not exactly a standard permanent.

Actually: H(W) = sum_sigma W[sigma(0)][sigma(1)] * W[sigma(1)][sigma(2)] * ... * W[sigma(n-2)][sigma(n-1)]

This is the sum of all open walk products of length n-1 visiting all vertices.

FORMULA: H(W) = tr(adj(I - zW) * W^{n-1}) / (n-1)! evaluated at... no, this isn't right.

Let me just compute and see what H(W) gives.

Author: opus-2026-03-07-S45
"""
from itertools import permutations
from fractions import Fraction
import math

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    """F(T,x) coefficients."""
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def transfer_matrix_poly(adj, n):
    """W(x)[u][v] as polynomial in x: W[u][v] = x if adj[u][v], 1 if adj[v][u], 0 if u=v.
    Represent as (constant, x_coefficient)."""
    W = [[(0, 0)]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v:
                W[u][v] = (0, 0)  # 0
            elif adj[u][v]:
                W[u][v] = (0, 1)  # x
            else:
                W[u][v] = (1, 0)  # 1
    return W

def poly_from_pair(pair):
    """Convert (constant, x_coeff) to polynomial list."""
    c, x = pair
    return [c, x] if x != 0 or c != 0 else [0]

def poly_mul(a, b):
    """Multiply two polynomials (lists of coefficients)."""
    if not a or not b:
        return [0]
    n = len(a) + len(b) - 1
    c = [0]*n
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i+j] += ai*bj
    return c

def poly_add(a, b):
    """Add two polynomials."""
    n = max(len(a), len(b))
    c = [0]*n
    for i in range(len(a)):
        c[i] += a[i]
    for i in range(len(b)):
        c[i] += b[i]
    return c

# ============================================================
# COMPUTE H(W,x) = sum over Ham paths of product of W entries
# ============================================================
print("=" * 60)
print("H(W,x) = F(T,x) verification")
print("=" * 60)

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    for bits in [0, 1, 10, (1 << m) - 1]:
        if bits >= (1 << m):
            continue
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)

        # Compute H(W,x) by summing over permutations
        H_poly = [0]*n
        for P in permutations(range(n)):
            # Product of W[P[i]][P[i+1]] for i=0..n-2
            # Each W entry is either x (if forward) or 1 (if backward)
            fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
            # Product = x^fwd * 1^(n-1-fwd) = x^fwd
            H_poly[fwd] += 1

        match = (H_poly == F)
        print(f"  bits={bits}: F={F}, H(W)={H_poly}, match={match}")

# ============================================================
# CHARACTERISTIC POLYNOMIAL of W(x) — does it encode F?
# ============================================================
print("\n" + "=" * 60)
print("CHARACTERISTIC POLYNOMIAL of W(x)")
print("=" * 60)

# For small n, compute det(lambda*I - W(x)) as a polynomial in lambda and x.
# Since W entries are linear in x, det is a polynomial in x of degree <= n.

# Actually, let's compute det(W) and per(W) as polynomials in x.

def compute_det_poly(adj, n):
    """det(W(x)) where W[u][v] = x*adj[u][v] + (1-adj[u][v])*(1-delta_{u,v})."""
    # Sum over permutations with sign
    det = [0]*(n+1)
    for P in permutations(range(n)):
        # Sign of permutation
        sign = 1
        inv = sum(1 for i in range(n) for j in range(i+1, n) if P[i] > P[j])
        sign = (-1)**inv

        # Product of W[i][P[i]]
        fwd = 0
        valid = True
        for i in range(n):
            if i == P[i]:
                valid = False  # W[i][i] = 0
                break
            if adj[i][P[i]]:
                fwd += 1
        if not valid:
            continue
        det[fwd] += sign

    return det

def compute_per_poly(adj, n):
    """per(W(x)): unsigned version."""
    per = [0]*(n+1)
    for P in permutations(range(n)):
        fwd = 0
        valid = True
        for i in range(n):
            if i == P[i]:
                valid = False
                break
            if adj[i][P[i]]:
                fwd += 1
        if not valid:
            continue
        per[fwd] += 1
    return per

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    seen = set()
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        det = compute_det_poly(adj, n)
        per = compute_per_poly(adj, n)

        print(f"  F={F}, det(W)={det}, per(W)={per}")

# ============================================================
# det(I - zW) — Irving-Omar generating function approach
# ============================================================
print("\n" + "=" * 60)
print("GENERATING FUNCTION: sum_k tr(W^k) * z^k")
print("=" * 60)

# tr(W^k) as polynomial in x gives the weighted closed walk count
# by length k. This is related to det(I-zW)^{-1}.
# The coefficient of z^{n-1} in [SOMETHING] should give F(T,x).

# Let's compute tr(W^k) for k=1,...,n and look for patterns.

def matrix_pow_poly(adj, n, k):
    """Compute W^k where W[u][v] is a polynomial in x.
    Entries of W^k are polynomials of degree <= k.
    Return as list of lists of polynomial lists."""
    # Initialize W
    W = [[[0]*(k+1) for _ in range(n)] for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u == v:
                pass  # already 0
            elif adj[u][v]:
                W[u][v][1] = 1  # x
            else:
                W[u][v][0] = 1  # 1

    # Compute W^k by repeated squaring (just use naive multiplication for small k)
    result = [[[0]*(k+1) for _ in range(n)] for _ in range(n)]
    # Identity matrix
    for i in range(n):
        result[i][i][0] = 1

    Wk = result
    Wcur = W

    for _ in range(k):
        new = [[[0]*(k+1) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for m in range(n):
                    # Wk[i][m] * Wcur[m][j]
                    for a in range(k+1):
                        for b in range(k+1):
                            if a + b <= k:
                                new[i][j][a+b] += Wk[i][m][a] * Wcur[m][j][b]
        Wk = new

    return Wk

for n in [4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    for bits in [0, (1 << m)//3, (1 << m) - 1]:
        if bits >= (1 << m):
            continue
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)

        print(f"  bits={bits}, F={F}")
        for k in range(1, n+1):
            Wk = matrix_pow_poly(adj, n, k)
            # tr(W^k)
            tr = [0]*(k+1)
            for i in range(n):
                for j in range(k+1):
                    tr[j] += Wk[i][i][j]
            print(f"    tr(W^{k}) = {[t for t in tr[:k+1]]}")

# ============================================================
# OPEN WALK SUM: does sum_v row_sum(W^{n-1})[v] = F(T,x)?
# ============================================================
print("\n" + "=" * 60)
print("ROW SUM OF W^{n-1} vs F(T,x)")
print("=" * 60)

# F(T,x) = sum over permutations of prod W[P[i]][P[i+1]]
# This is NOT the same as row/column sums of W^{n-1} because
# W^{n-1}[u][v] sums over ALL walks of length n-1 from u to v,
# not just Hamiltonian paths.

# But maybe there's a correction factor?

for n in [4]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    for bits in [0, 5, 10]:
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)

        Wk = matrix_pow_poly(adj, n, n-1)

        # Sum all entries of W^{n-1}
        total = [0]*n
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    total[k] += Wk[i][j][k]

        # Row sums
        row_sum = [0]*n
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    row_sum[k] += Wk[i][j][k]

        print(f"  bits={bits}: F={F}")
        print(f"    sum(W^{n-1}) = {total}")
        print(f"    ratio: {[total[k]/F[k] if F[k] != 0 else 'inf' for k in range(n)]}")

# ============================================================
# CAYLEY-HAMILTON: Can we express F via char poly of W?
# ============================================================
print("\n" + "=" * 60)
print("EIGENVALUES OF W(x) AT SPECIFIC x VALUES")
print("=" * 60)

import numpy as np

for n in [5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")

    seen = set()
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        # Evaluate W at x=1: W(1) has all off-diagonal entries = 1 or 1.
        # Actually W(1)[u][v] = 1 for u!=v, 0 for u=v. That's J-I.
        # This is the same for ALL tournaments!
        # eigenvalues of J-I: n-1 (once), -1 (n-1 times)

        # At x=2:
        W2 = np.zeros((n, n))
        for u in range(n):
            for v in range(n):
                if u == v:
                    continue
                W2[u][v] = 2 if adj[u][v] else 1

        eigs = sorted(np.linalg.eigvals(W2).real, reverse=True)

        # F(T,2) = sum F[k]*2^k
        F2 = sum(F[k]*2**k for k in range(n))

        if len(seen) <= 8:
            print(f"  F={F}, F(2)={F2}, eigs(W(2))={[f'{e:.2f}' for e in eigs]}")
            # Product of eigenvalues = det(W(2))
            det_W2 = np.linalg.det(W2)
            print(f"    det(W(2)) = {det_W2:.1f}")
