#!/usr/bin/env python3
"""
irving_omar_detper.py — Irving-Omar det+per formula for H(T).

FROM THEIR PAPER (arXiv:2412.10572):
  ham(D) = sum_{S ⊆ [n]} det(Ā[S]) * per(A[S^c])

where A is the adjacency matrix, Ā = J-I-A is the complement.

For tournaments: Ā[i][j] = 1 - A[i][j] for i≠j, so Ā = transpose(A)
(since A[i][j] + A[j][i] = 1 for i≠j, and Ā[i][j] = A[j][i]).

So for tournaments:
  H(T) = sum_{S ⊆ [n]} det(A^T[S]) * per(A[S^c])
        = sum_{S ⊆ [n]} det(A[S]) * per(A[S^c])   (since det(A^T) = det(A))

QUESTION 1: Verify this formula.
QUESTION 2: Can we put x-weights in to get F(T,x)?
QUESTION 3: What is the "walk generating function" W_D(z)?

Also their formula:
  W_D(z) = det(I + z*X*Ā) / det(I - z*X*A)

where X = diag(x_1,...,x_n). At X=I:
  W_D(z) = det(I + z*A^T) / det(I - z*A)

The coefficient of z^k in W_D gives the k-th moment of walks.

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
import numpy as np

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

def compute_H(adj, n):
    count = 0
    for P in permutations(range(n)):
        if all(adj[P[i]][P[i+1]] for i in range(n-1)):
            count += 1
    return count

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def det_submatrix(adj, S):
    """det of adj restricted to rows and columns in S."""
    k = len(S)
    if k == 0:
        return 1
    mat = [[adj[S[i]][S[j]] for j in range(k)] for i in range(k)]
    # Compute det by Leibniz formula (small matrices)
    if k == 1:
        return mat[0][0]
    total = 0
    for perm in permutations(range(k)):
        sign = 1
        inv = sum(1 for i in range(k) for j in range(i+1, k) if perm[i] > perm[j])
        sign = (-1)**inv
        prod = 1
        for i in range(k):
            prod *= mat[i][perm[i]]
        total += sign * prod
    return total

def per_submatrix(adj, S):
    """Permanent of adj restricted to rows and columns in S."""
    k = len(S)
    if k == 0:
        return 1
    mat = [[adj[S[i]][S[j]] for j in range(k)] for i in range(k)]
    total = 0
    for perm in permutations(range(k)):
        prod = 1
        for i in range(k):
            prod *= mat[i][perm[i]]
        total += sign * prod  # wait, permanent has no sign
    return total

def per_submatrix(adj, S):
    """Permanent of adj restricted to rows and columns in S."""
    k = len(S)
    if k == 0:
        return 1
    total = 0
    for perm in permutations(range(k)):
        prod = 1
        for i in range(k):
            prod *= adj[S[i]][S[perm[i]]]
        total += prod
    return total

# ============================================================
# VERIFY: H(T) = sum_S det(A^T[S]) * per(A[S^c])
# ============================================================
print("=" * 60)
print("IRVING-OMAR DET+PER FORMULA VERIFICATION")
print("=" * 60)

for n in [3, 4, 5]:
    m = n*(n-1)//2
    ok = 0
    total = 0
    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        H = compute_H(adj, n)

        # Compute sum_S det(A^T[S]) * per(A[S^c])
        # A^T[i][j] = adj[j][i]
        adj_T = [[adj[j][i] for j in range(n)] for i in range(n)]

        formula = 0
        for r in range(n+1):
            for S in combinations(range(n), r):
                Sc = [x for x in range(n) if x not in S]
                d = det_submatrix(adj_T, list(S))
                p = per_submatrix(adj, Sc)
                formula += d * p

        total += 1
        if formula == H:
            ok += 1
        else:
            if total - ok <= 3:
                print(f"  FAIL at n={n} bits={bits}: H={H}, formula={formula}")

    print(f"n={n}: {ok}/{total} pass")

# ============================================================
# WEIGHTED VERSION: can we get F(T,x)?
# ============================================================
print("\n" + "=" * 60)
print("WEIGHTED DET+PER FOR F(T,x)?")
print("=" * 60)

# The natural weight: replace adj[i][j] by x if adj[i][j]=1, 1 if adj[j][i]=1
# i.e., use W(x) matrix. Then per(W[S^c]) counts cycle covers of S^c weighted by x.
# And det(W^T[S]) counts signed cycle covers of S weighted by x.
# But this gives CYCLE COVER weights, not Hamiltonian path weights.

# Alternative: track forward edges differently.
# F(T,x) = sum_P x^{fwd(P)} = sum_P prod_{i} (x if adj[P[i]][P[i+1]] else 1)
# This is the "path permanent" of W(x), NOT the det+per sum.

# But maybe the det+per sum with W(x) gives something interesting?

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

        # Compute det+per sum with x-weights
        # W(x)[i][j] = x if adj[i][j]=1, 1 if adj[j][i]=1, 0 if i=j
        # W^T(x)[i][j] = W(x)[j][i] = x if adj[j][i]=1, 1 if adj[i][j]=1
        # But W^T(x)[i][j] = x if adj_T[i][j]=1, 1 otherwise (for i≠j)
        # So W^T(x) has x on BACKWARD arcs and 1 on FORWARD arcs
        # = W(1/x) * x (up to scaling)... hmm

        # Let's just compute the polynomial sum_S det(W^T[S]) * per(W[S^c])
        # Each entry is either x or 1 (or 0 on diagonal)
        # det and per are polynomials in x

        # For small n, compute numerically at several x values
        results = {}
        for x_val in [1, 2, 3, -1]:
            W = [[0]*n for _ in range(n)]
            WT = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    W[i][j] = x_val if adj[i][j] else 1
                    WT[i][j] = x_val if adj[j][i] else 1

            formula = 0
            for r in range(n+1):
                for S in combinations(range(n), r):
                    Sc = [v for v in range(n) if v not in S]
                    d = det_submatrix(WT, list(S))
                    p = per_submatrix(W, Sc)
                    formula += d * p
            results[x_val] = formula

        F_vals = {x: sum(F[k]*x**k for k in range(n)) for x in [1, 2, 3, -1]}

        if len(seen) <= 5:
            print(f"  F={F}")
            print(f"  F(x) at 1,2,3,-1: {[F_vals[x] for x in [1,2,3,-1]]}")
            print(f"  det+per at 1,2,3,-1: {[results[x] for x in [1,2,3,-1]]}")
            # Check if they're proportional
            ratios = []
            for x in [1, 2, 3, -1]:
                if F_vals[x] != 0:
                    ratios.append(results[x] / F_vals[x])
                else:
                    ratios.append('inf')
            print(f"  ratio: {ratios}")

# ============================================================
# WALK GENERATING FUNCTION: W_D(z) = det(I+zA^T)/det(I-zA)
# ============================================================
print("\n" + "=" * 60)
print("WALK GF: W_D(z) = det(I+zA^T)/det(I-zA)")
print("=" * 60)

# Expand as power series in z. Coefficient of z^k gives tr(A^k) type info.
# For tournaments, A^T = J-I-A (complement).

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

        A = np.array(adj, dtype=float)
        AT = A.T

        # Characteristic polynomials
        # det(I - zA) = char poly of A evaluated at 1/z... no.
        # det(I - zA) as polynomial in z
        # = sum_k (-1)^k e_k(eigenvalues of A) z^k
        # For n=4, degree 4 in z.

        eigs_A = np.linalg.eigvals(A)
        eigs_AT = np.linalg.eigvals(AT)  # same eigenvalues

        # det(I+zA^T) = prod_i (1 + z*lambda_i)
        # det(I-zA) = prod_i (1 - z*lambda_i)
        # W(z) = prod_i (1+z*lambda_i)/(1-z*lambda_i)

        # At z=1: W(1) = det(I+A^T)/det(I-A)
        detIpAT = np.linalg.det(np.eye(n) + AT)
        detImA = np.linalg.det(np.eye(n) - A)

        if len(seen) <= 8:
            H = sum(F)
            W1 = detIpAT / detImA if abs(detImA) > 1e-10 else float('inf')
            print(f"  F={F}, H={H}")
            print(f"  det(I+A^T)={detIpAT:.0f}, det(I-A)={detImA:.0f}, W(1)={W1:.4f}")
            print(f"  eigenvalues of A: {sorted([f'{e.real:.2f}+{e.imag:.2f}i' for e in eigs_A])}")

# ============================================================
# KEY INSIGHT: det(I-zA) FOR TOURNAMENTS
# ============================================================
print("\n" + "=" * 60)
print("CHARACTERISTIC POLYNOMIAL OF TOURNAMENT ADJACENCY MATRIX")
print("=" * 60)

# For a tournament, A + A^T = J - I (all-ones minus identity)
# So eigenvalues of A satisfy: if lambda is eigenvalue of A,
# then lambda + lambda_bar is eigenvalue of J-I = n-1 (once) and -1 (n-1 times)
# Actually A + A^T has eigenvalues n-1 (once) and -1 (n-1 times)
# But A is not symmetric, so this doesn't directly constrain eigs of A.

# However, A + A^T = J-I implies:
# tr(A) = 0 (diagonal is 0)
# tr(A^2) = sum_{i,j} A[i][j]*A[j][i] = sum_{i<j} (A[i][j]*A[j][i] + A[j][i]*A[i][j])
# But A[i][j] + A[j][i] = 1, and A[i][j]*A[j][i] = 0 always (exactly one is 1).
# So tr(A^2) = 0. Wait, that means sum of squares of eigenvalues = 0!

# tr(A^k) counts closed walks of length k.
# tr(A^2) = # arcs (i,j) with A[i][j]=1 AND (j,i) with A[j][i]=1 for some path i->j->i
# No, tr(A^2) = sum_i (A^2)[i][i] = sum_i sum_j A[i][j]*A[j][i] = 0
# since A[i][j]*A[j][i] = 0 for all i,j (tournament: exactly one direction).

for n in [4, 5]:
    m = n*(n-1)//2
    seen = set()
    char_polys = set()

    for bits in range(1 << m):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        A = np.array(adj, dtype=float)
        # Characteristic polynomial coefficients
        coeffs = np.round(np.poly(A)).astype(int).tolist()
        char_polys.add(tuple(coeffs))

        if len(seen) <= 5:
            H = sum(F)
            t3 = sum(1 for i in range(n) for j in range(n) for k in range(n)
                     if i<j<k and adj[i][j] and adj[j][k] and adj[k][i]) + \
                 sum(1 for i in range(n) for j in range(n) for k in range(n)
                     if i<j<k and adj[i][k] and adj[k][j] and adj[j][i])
            print(f"  F={F}, H={H}, t3={t3}")
            print(f"  char poly of A: {coeffs}")
            print(f"  tr(A^2)=0? {int(np.trace(A@A))}")
            print(f"  tr(A^3) = {int(np.round(np.trace(A@A@A)))} (should be 6*t3_directed? = 3*t3)")

    print(f"\nn={n}: {len(char_polys)} distinct characteristic polynomials, {len(seen)} F-vectors")
