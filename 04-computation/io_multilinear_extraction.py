#!/usr/bin/env python3
"""
Irving-Omar Walk GF: Multilinear Extraction and Transfer Matrix Connection.

IO Walk GF: W(z) = det(I + zA^T) / det(I - zA)
where A is the adjacency matrix of tournament T.

At the commutative level, this is a product over eigenvalues:
  W(z) = prod_i (1 + z*lambda_i^*) / (1 - z*lambda_i)

For the MULTILINEAR version with vertex-tracking variables x_1,...,x_n:
  W(z,x) = det(I + zXA^T) / det(I - zXA)  where X = diag(x_1,...,x_n)

The multilinear part (coefficient of x_1*x_2*...*x_n) in the expansion
relates to Hamiltonian structures.

KEY QUESTION: Can we extract M[a,b] from the multilinear part of W?

The transfer matrix M involves SIGNED path decompositions.
The IO GF involves cycle-cover expansions with SIGNED weights.
Both have (-1)^k type signs. Can these be identified?

APPROACH: Compute W(z) explicitly for small tournaments and
extract coefficients, comparing with M entries.
"""

import numpy as np
from itertools import permutations
from collections import defaultdict

def adjacency_matrix(n, tournament_dict):
    A = np.zeros((n, n))
    for (i,j), v in tournament_dict.items():
        A[i][j] = v
    return A

def make_tournament(n, bits):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    T = {}
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T

def ham_paths(A_mat):
    n = A_mat.shape[0]
    paths = []
    for perm in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if A_mat[perm[k]][perm[k+1]] < 0.5: ok = False; break
        if ok: paths.append(perm)
    return paths

def transfer_matrix_consec(A_mat):
    n = A_mat.shape[0]
    M = np.zeros((n, n), dtype=int)
    for p in ham_paths(A_mat):
        for a in range(n):
            M[a][a] += (-1)**(list(p).index(a))
        for j in range(n-1):
            a, b = p[j], p[j+1]
            M[a][b] += (-1)**j
            M[b][a] += (-1)**j
    return M

# =====================================================================
# W(z) = det(I + zA^T) / det(I - zA) as a rational function
# =====================================================================
print("=" * 70)
print("IRVING-OMAR WALK GF: W(z) COMPUTATION")
print("=" * 70)

from numpy.polynomial import polynomial as P

# For n=3: 3-cycle
n = 3
A3 = np.array([[0,1,0],[0,0,1],[1,0,0]], dtype=float)
print(f"\nn=3, 3-cycle:")
print(f"  A eigenvalues: {sorted(np.linalg.eigvals(A3).real, reverse=True)}")

# det(I - zA) = 1 - z^3 (for 3-cycle permutation matrix)
# det(I + zA^T) = 1 + z^3 (A^T = A^2 = A^{-1})
# W(z) = (1+z^3)/(1-z^3) = 1 + 2z^3/(1-z^3) = 1 + 2z^3 + 2z^6 + ...
print(f"  W(z) = (1+z^3)/(1-z^3)")
print(f"  Coefficients: W_0=1, W_3=2, W_6=2, ...")

# Taylor expansion via matrix power series
# (I - zA)^{-1} = I + zA + z^2 A^2 + z^3 A^3 + ...
# At z, det(I-zA) has roots at eigenvalues^{-1}

# Compute W(z) numerically at several points
for z in [0.1, 0.2, 0.5]:
    num = np.linalg.det(np.eye(n) + z * A3.T)
    den = np.linalg.det(np.eye(n) - z * A3)
    W = num / den
    W_exact = (1 + z**3) / (1 - z**3)
    print(f"  W({z}) = {W:.6f} (exact: {W_exact:.6f})")

# M for 3-cycle
M3 = transfer_matrix_consec(A3)
H3 = len(ham_paths(A3))
print(f"  H = {H3}, M = {M3.tolist()}")

# =====================================================================
# n=5: Paley tournament
# =====================================================================
print(f"\nn=5, Paley:")
n = 5
A5 = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i != j and (j-i)%n in [1, 2]:
            A5[i][j] = 1

eigs = np.linalg.eigvals(A5)
print(f"  A eigenvalues: {sorted(eigs.real, reverse=True)}")

# Compute W(z) at z=1/2
z = 0.5
num = np.linalg.det(np.eye(n) + z * A5.T)
den = np.linalg.det(np.eye(n) - z * A5)
W = num / den
print(f"  W(0.5) = {W:.6f}")
print(f"  num = {num:.6f}, den = {den:.6f}")

# Check IO reciprocity: W(-z) at commutative level
# W_T(-z) = det(I - zA^T) / det(I + zA)
num_neg = np.linalg.det(np.eye(n) - z * A5.T)
den_neg = np.linalg.det(np.eye(n) + z * A5)
W_neg = num_neg / den_neg
print(f"  W(-0.5) = {W_neg:.6f}")
print(f"  W(0.5) * W(-0.5) = {W * W_neg:.6f}")

# Actually, IO reciprocity in even parity:
# W(-z,-r) = W(z,r) where r tracks arc weight
# At commutative level with r=1: W(-z) should relate to W(z)
# det(I - zA^T) / det(I + zA) vs det(I + zA^T) / det(I - zA)
# = [det(I-zA^T) * det(I-zA)] / [det(I+zA) * det(I+zA^T)] * [det(I+zA^T)/det(I-zA)] * ... complex

# Simpler: W(-z) * W(z) = ?
# W(z) * W(-z) = [det(I+zA^T)det(I-zA^T)] / [det(I-zA)det(I+zA)]
#              = det(I - z^2 (A^T)^2) / det(I - z^2 A^2)
# For tournament: A + A^T = J - I (all 1s minus identity)
# A^2 and (A^T)^2 have same eigenvalues (since A^T has conjugate eigs)
# So det(I - z^2 A^2) = det(I - z^2 (A^T)^2) and W(z)*W(-z) = 1!

print(f"\n  IDENTITY: W(z) * W(-z) = {W * W_neg:.6f} (should = 1)")

# =====================================================================
# W(z) * W(-z) = 1 for all tournaments?
# =====================================================================
print()
print("=" * 70)
print("IDENTITY: W(z) * W(-z) = 1?")
print("=" * 70)

for n in [3, 4, 5]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    all_pass = True
    for bits in range(min(2**len(edges), 100)):
        T = make_tournament(n, bits)
        A = adjacency_matrix(n, T)
        z = 0.3
        Wp = np.linalg.det(np.eye(n) + z * A.T) / np.linalg.det(np.eye(n) - z * A)
        Wm = np.linalg.det(np.eye(n) - z * A.T) / np.linalg.det(np.eye(n) + z * A)
        if abs(Wp * Wm - 1) > 1e-10:
            all_pass = False
            break
    print(f"  n={n}: W(z)*W(-z) = 1? {all_pass}")

# PROOF: W(z)*W(-z) = det(I+zA^T)*det(I-zA^T) / [det(I-zA)*det(I+zA)]
# Numerator: det((I+zA^T)(I-zA^T)) = det(I - z^2 A^T A^T) = det(I - z^2 (A^T)^2)
# Denominator: det((I-zA)(I+zA)) = det(I - z^2 A^2)
# For tournaments: (A^T)^2 and A^2 have same eigenvalues
# (A^T has eigenvalues = complex conjugates of A's eigenvalues)
# So det(I - z^2 (A^T)^2) = conj(det(I - z^2 A^2)) ... hmm
# Actually for REAL matrices: (A^T)^2 = (A^2)^T, so same determinant.
# det(I - z^2 (A^2)^T) = det((I - z^2 A^2)^T) = det(I - z^2 A^2). ✓
print("  PROVED: W(z)*W(-z) = 1 for all tournaments (since det(M^T) = det(M))")

# =====================================================================
# W(z) COEFFICIENT EXTRACTION: W(z) = sum c_k z^k
# =====================================================================
print()
print("=" * 70)
print("W(z) TAYLOR COEFFICIENTS AND TRANSFER MATRIX")
print("=" * 70)

# W(z) = det(I+zA^T)/det(I-zA)
# = [1 + z*tr(A^T) + z^2*(tr^2-tr(A^T^2))/2 + ...]
#   / [1 - z*tr(A) + z^2*(tr^2-tr(A^2))/2 - ...]

# For tournament: tr(A) = 0 (no self-loops).
# W(z) near z=0: W(0) = 1.
# W'(0) = tr(A^T) * 1 + tr(A) * 1 = 0 + 0 = 0? Let me check.

# Actually: W(z) = det(I+zA^T) * det(I-zA)^{-1}
# d/dz [det(I+zA^T)] at z=0 = tr(A^T) * det(I) = tr(A^T) = 0
# d/dz [det(I-zA)^{-1}] at z=0 = tr(A) * 1 = 0
# So W'(0) = 0. And W(0) = 1.

# Numerical coefficient extraction
n = 5
A5_paley = A5.copy()
zvals = np.linspace(-0.4, 0.4, 1000)
Wvals = []
for z in zvals:
    num = np.linalg.det(np.eye(n) + z * A5_paley.T)
    den = np.linalg.det(np.eye(n) - z * A5_paley)
    Wvals.append(num / den)
Wvals = np.array(Wvals)

# Fit polynomial: W(z) ≈ sum c_k z^k
from numpy.polynomial.polynomial import polyfit
coeffs = polyfit(zvals, Wvals, 10)
print(f"  Paley T_5 W(z) Taylor coefficients:")
for k in range(11):
    if abs(coeffs[k]) > 1e-6:
        print(f"    c_{k} = {coeffs[k]:.4f}")

# For Paley T_5: eigenvalues of A are related to Gauss sums
# A has eigenvalues: 2 (once), (-1±sqrt(5))/2 (each twice)
# The golden ratio appears!
eigs5 = np.linalg.eigvals(A5_paley)
print(f"\n  Paley T_5 eigenvalues of A: {sorted(eigs5.real, reverse=True)}")

# W(z) = prod_i (1+z*conj(lambda_i))/(1-z*lambda_i)
# For real eigenvalues: (1+z*lambda_i)/(1-z*lambda_i)
# For complex pair (a+bi, a-bi):
#   [(1+z(a-bi))(1+z(a+bi))] / [(1-z(a+bi))(1-z(a-bi))]
#   = (1+2az+z^2(a^2+b^2)) / (1-2az+z^2(a^2+b^2))

print(f"\n  W(z) as eigenvalue product:")
for e in sorted(eigs5, key=lambda x: -x.real):
    print(f"    lambda = {e:.4f}: factor = (1+{e.real:.4f}*z)/(1-{e.real:.4f}*z)")

# =====================================================================
# KEY CONNECTION: W evaluated at z=1 gives...?
# =====================================================================
print()
print("=" * 70)
print("W(z) AT SPECIAL VALUES")
print("=" * 70)

for n in [3, 5]:
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    for bits in range(2**len(edges)):
        A_arr = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1: A_arr[i][j] = 1
            else: A_arr[j][i] = 1
        key = tuple(tuple(row) for row in A_arr)
        min_key = key
        for perm in permutations(range(n)):
            pkey = tuple(tuple(A_arr[perm[i]][perm[j]] for j in range(n)) for i in range(n))
            if pkey < min_key: min_key = pkey
        if min_key in seen: continue
        seen.add(min_key)

        A_np = np.array(A_arr, dtype=float)
        H = len(ham_paths(A_np))

        # W at z = 1/(n-1): ?
        z1 = 1.0 / (n-1)
        num = np.linalg.det(np.eye(n) + z1 * A_np.T)
        den = np.linalg.det(np.eye(n) - z1 * A_np)
        if abs(den) > 1e-12:
            W1 = num / den
        else:
            W1 = float('inf')

        # W at z = 1/sqrt(n):
        z2 = 1.0 / np.sqrt(n)
        num2 = np.linalg.det(np.eye(n) + z2 * A_np.T)
        den2 = np.linalg.det(np.eye(n) - z2 * A_np)
        W2 = num2 / den2 if abs(den2) > 1e-12 else float('inf')

        scores = sorted([sum(A_arr[i]) for i in range(n)], reverse=True)
        print(f"  n={n}, H={H:3d}, scores={scores}: "
              f"W(1/{n-1})={W1:8.3f}, W(1/sqrt({n}))={W2:8.3f}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
KEY RESULTS:

1. W(z)*W(-z) = 1 for ALL tournaments.
   This is the commutative IO reciprocity identity.
   PROOF: det(M^T) = det(M) for any matrix.

2. W(z) Taylor coefficients: c_0 = 1, c_1 = 0, c_2 = 0
   (since tr(A) = 0 for tournaments).
   First nonzero coefficient is c_3 = 2*(#3-cycles) at n≥3.

3. The IO walk GF W(z) encodes CYCLE COVERS of T.
   The transfer matrix M encodes SIGNED PATH statistics.
   These are DIFFERENT combinatorial structures.

4. W(z) at rational z values does NOT directly give H or M entries.
   The connection to M must go through the MULTILINEAR extraction
   from the noncommuting version W(z, x_1,...,x_n).

5. The identity W(z)*W(-z)=1 is the IO analog of THM-030 symmetry.
   Both encode "even parity" structure:
   - THM-030: M symmetric despite asymmetric tournament
   - IO: W(z)*W(-z) = 1 despite W(z) ≠ W(-z)
""")
