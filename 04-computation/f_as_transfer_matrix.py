#!/usr/bin/env python3
"""
f_as_transfer_matrix.py — Express F(T,x) via the transfer matrix method.

F(T,x) = sum over all n! permutations P of prod_{i=0}^{n-2} w(P_i, P_{i+1})
where w(u,v) = x if A[u][v] = 1, else 1 (= x^{A[u][v]}).

This is a SUM OVER HAMILTONIAN PATHS of a weighted tournament.
The weight matrix is W(x)[u][v] = x * A[u][v] + (1 - A[u][v]).

For a tournament, A[u][v] + A[v][u] = 1, so:
  W[u][v] = x * A[u][v] + A[v][u]
  W[v][u] = x * A[v][u] + A[u][v]
  W[u][v] * W[v][u] = (x*a + (1-a)) * (x*(1-a) + a)
    where a = A[u][v] in {0,1}
  = (x*1 + 0)*(0 + 1) = x if a=1 (so W[u][v]*W[v][u] = x always!)
  = (0 + 1)*(x + 0) = x if a=0

So W[u][v] * W[v][u] = x for ALL pairs! The matrix W satisfies
a MULTIPLICATIVE CONSTRAINT: opposite entries multiply to x.

This means W(x) = (x^{1/2}) * M where M[u][v] * M[v][u] = 1.
So M is "skew-unitary" in a multiplicative sense.

For x > 0: W = sqrt(x) * M where M[u][v] = sqrt(x)^{2A[u][v]-1}
  If A[u][v]=1: M[u][v] = sqrt(x), M[v][u] = 1/sqrt(x)
  If A[u][v]=0: M[u][v] = 1/sqrt(x), M[v][u] = sqrt(x)

F(T,x) = sum over HP of prod W = x^{(n-1)/2} * sum over HP of prod M

Actually: prod W = prod x^{A[P_i][P_{i+1}]} = x^{fwd(P)}.
And prod M = x^{fwd(P) - (n-1)/2} = x^{fwd - (n-1)/2}.

The palindrome means: E[x^{fwd}] = E[x^{n-1-fwd}] (path reversal).

NEW IDEA: Can we compute F(T,x) as det(I - x*M) or similar algebraic expression?

For Hamiltonian path counting, the PERMANENT of the adjacency matrix counts
cycle covers, not paths. But the IMMANANT or the inclusion-exclusion permanent
might work.

F(T,1) = n! = permanent of J-I (all-ones minus identity).

Let me try: Does det(I - t*W(x)) relate to F(T,x)?

Author: opus-2026-03-07-S44
"""
import numpy as np
from itertools import permutations
import math
import random

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def weighted_matrix(A, n, x):
    """W[i][j] = x * A[i][j] + (1 - A[i][j]) for i != j, 0 on diagonal."""
    W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                W[i][j] = x * A[i][j] + (1 - A[i][j])
    return W

# Test: does tr(adj(I - t*W)) give anything useful?
print("=== Transfer matrix explorations ===")

for n in [3, 4, 5]:
    A = tournament_from_bits(0, n)  # transitive
    F = compute_F(A, n)
    print(f"\nn={n}: F = {F}")

    # Compute F(T,x) via matrix power sum
    # F(T,x) = sum_v [row v of W^{n-1}]_{sum} ... no, paths not walks.

    # Inclusion-exclusion for Hamiltonian paths:
    # H(T, weights) = sum_{S subset V} (-1)^{n-|S|} * (products of walks in S)
    # This is the permanent-like formula.

    # Actually: sum over HP of weight products =
    # sum_{sigma in S_n as path} prod W[sigma(i)][sigma(i+1)]
    # This is related to immanant computations.

    # Let me try: characteristic polynomial of W(x)
    for x_val in [0.5, 1.0, 2.0]:
        W = weighted_matrix(A, n, x_val)
        det_val = np.linalg.det(W)
        eigenvals = np.linalg.eigvals(W)

        # F(T, x_val) by direct computation
        F_val = sum(F[k] * x_val**k for k in range(n))

        print(f"  x={x_val}: F(T,x)={F_val:.1f}, det(W)={det_val:.4f}, "
              f"eigs={[f'{e:.3f}' for e in sorted(eigenvals, key=lambda z: -abs(z))]}")

    # Check: permanent of W(x) vs F(T,x)
    # Permanent by definition: perm(M) = sum_{sigma in S_n} prod M[i][sigma(i)]
    # This counts CYCLE COVERS, not paths.
    # F(T,x) counts HAMILTONIAN PATHS with weights.
    # They are different objects!

    # But for n=3, let's compare
    if n <= 5:
        for x_val in [2.0]:
            W = weighted_matrix(A, n, x_val)
            # Permanent
            perm = 0
            for sigma in permutations(range(n)):
                prod_val = 1
                for i in range(n):
                    prod_val *= W[i][sigma[i]]
                perm += prod_val

            F_val = sum(F[k] * x_val**k for k in range(n))
            print(f"  x={x_val}: perm(W)={perm:.1f}, F(T,x)={F_val:.1f}")

# ============================================================
# NOVEL: F(T,x) via INCLUSION-EXCLUSION on CYCLE STRUCTURE
# ============================================================
print("\n" + "=" * 60)
print("F(T,x) VIA CYCLE STRUCTURE DECOMPOSITION")
print("=" * 60)

# A permutation of [n] decomposes into disjoint cycles.
# A Hamiltonian PATH is a permutation with a SINGLE cycle of length n
# ... no, a path is not a cycle.

# Actually: a permutation sigma can be viewed as sigma(1), sigma(2), ..., sigma(n).
# This IS a sequence (path) visiting each vertex once.
# The "forward edge weight" is prod x^{A[sigma(i)][sigma(i+1)]} for i=1..n-1.

# The cycle decomposition of sigma (as a bijection) is different from the path.
# A sigma with cycle type (n) [single n-cycle] corresponds to n distinct paths.

# Hmm, this is getting complicated. Let me try a completely different approach.

# ============================================================
# NOVEL: F(T,x) and the IHARA ZETA FUNCTION
# ============================================================
print("\n" + "=" * 60)
print("IHARA ZETA FUNCTION CONNECTION")
print("=" * 60)

# The Ihara zeta function of a graph G:
# Z_G(u) = prod_{[C] prime} (1 - u^|C|)^{-1} = 1/det(I - uA + u^2(D-I))
# For directed graphs: Z_D(u) = 1/det(I - u*A)  (simpler!)

# For a tournament: A is the adjacency matrix (0-1).
# det(I - u*A) = characteristic polynomial of A evaluated at 1/u.

# Connection to F(T,x): F(T,x) involves PATHS not cycles.
# But the generating function for WALKS is (I - x*A)^{-1}.
# And Hamiltonian paths can be extracted via inclusion-exclusion
# from the walk generating function.

# Let me check: is det(I - u*A) related to cycle counts?
n = 5
A_mat = tournament_from_bits(0, n)
A_np = np.array(A_mat, dtype=float)

char_poly = np.poly(A_np)  # characteristic polynomial
print(f"n=5 transitive: char_poly = {[f'{c:.1f}' for c in char_poly]}")
print(f"  This is det(lambda*I - A)")

# det(I - u*A) = u^n * det(1/u * I - A) = u^n * char_poly(1/u)
# For u=1: det(I-A). For u=-1: det(I+A).

det_IA = np.linalg.det(np.eye(n) - A_np)
det_IpA = np.linalg.det(np.eye(n) + A_np)
print(f"  det(I-A) = {det_IA:.1f}")
print(f"  det(I+A) = {det_IpA:.1f}")
print(f"  H(T) = {compute_F(A_mat, n)[-1]}")

# Check for several tournaments
print("\nComparing det(I+A), det(I-A) with H(T) at n=5:")
m = n*(n-1)//2
seen = set()
for bits in range(1 << m):
    A_mat = tournament_from_bits(bits, n)
    F = compute_F(A_mat, n)
    H = F[n-1]
    key = H
    if key in seen:
        continue
    seen.add(key)

    A_np = np.array(A_mat, dtype=float)
    det_IpA = np.linalg.det(np.eye(n) + A_np)
    det_ImA = np.linalg.det(np.eye(n) - A_np)

    # Also try det(I + 2A) and det(I - 2A)
    det_I2A = np.linalg.det(np.eye(n) + 2*A_np)

    # B = 2A - J + I (signed matrix)
    B = 2*A_np - np.ones((n,n)) + np.eye(n)
    det_B = np.linalg.det(B)

    print(f"  H={H:3d}: det(I+A)={det_IpA:8.1f}, det(I-A)={det_ImA:8.1f}, "
          f"det(I+2A)={det_I2A:8.1f}, det(B)={det_B:8.1f}")
