#!/usr/bin/env python3
"""
HYPOTHESIS AA: M satisfies a functional equation involving A.

Since M is not a polynomial in A, maybe it satisfies:
  M * A = A * M (commutation)
  or M * A + A^T * M = something simple
  or det(M) = some function of H

Let me check what algebraic relations M satisfies.
"""
from itertools import permutations
import numpy as np
import random

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_M(A, n):
    M = [[0]*n for _ in range(n)]
    V = set(range(n))
    for a in range(n):
        for b in range(n):
            if a == b:
                total = 0
                for P in permutations(range(n)):
                    if all(A[P[i]][P[i+1]] for i in range(n-1)):
                        total += (-1)**P.index(a)
                M[a][a] = total
            else:
                others = sorted(V - {a, b})
                total = 0
                for mask in range(1 << len(others)):
                    S = set()
                    for idx in range(len(others)):
                        if (mask >> idx) & 1:
                            S.add(others[idx])
                    R = set(others) - S
                    verts_E = sorted(S | {a})
                    E_a = sum(1 for P in permutations(verts_E)
                              if P[-1] == a and all(A[P[i]][P[i+1]] for i in range(len(P)-1)))
                    verts_B = sorted({b} | R)
                    B_b = sum(1 for P in permutations(verts_B)
                              if P[0] == b and all(A[P[i]][P[i+1]] for i in range(len(P)-1)))
                    total += (-1)**len(S) * E_a * B_b
                M[a][b] = total
    return M

random.seed(42)
n = 5

print("=== Algebraic relations of M ===")
for trial in range(5):
    A = random_tournament(n)
    M = np.array(compute_M(A, n), dtype=float)
    A_np = np.array(A, dtype=float)
    B = 2*A_np - np.ones((n,n)) + np.eye(n)  # Signed adjacency

    H = int(np.trace(M))
    det_M = np.linalg.det(M)
    eigs_M = np.sort(np.linalg.eigvals(M))[::-1]

    print(f"\n  trial {trial}: H={H}")

    # Does M commute with A?
    comm = M @ A_np - A_np @ M
    print(f"  [M,A] max entry: {np.max(np.abs(comm)):.2f}")

    # Does M commute with B?
    comm_B = M @ B - B @ M
    print(f"  [M,B] max entry: {np.max(np.abs(comm_B)):.2f}")

    # Does M anticommute with B?
    anti_B = M @ B + B @ M
    print(f"  {{M,B}} max entry: {np.max(np.abs(anti_B)):.2f}")

    # M^2?
    M2 = M @ M
    print(f"  tr(M^2) = {np.trace(M2):.0f}, sum(M_ij^2) = {np.sum(M**2):.0f}")

    # Is M^2 a polynomial in A?
    # Try M^2 = a*I + b*M
    # If so: eigenvalues satisfy lambda^2 = a + b*lambda
    # All eigenvalues are roots of x^2 - b*x - a = 0
    # So M has at most 2 distinct eigenvalues!
    print(f"  eigenvalues of M: {np.round(eigs_M, 3)}")
    print(f"  det(M) = {det_M:.2f}")

    # Check if M has ≤ 2 distinct eigenvalues
    unique_eigs = np.unique(np.round(eigs_M, 2))
    print(f"  distinct eigenvalues: {unique_eigs}")

    # Is M^2 = aI + bM + cJ?
    J = np.ones((n,n))
    X = np.column_stack([np.eye(n).flatten(), M.flatten(), J.flatten()])
    beta = np.linalg.lstsq(X, M2.flatten(), rcond=None)[0]
    resid = M2.flatten() - X @ beta
    if np.max(np.abs(resid)) < 0.01:
        print(f"  M^2 = {beta[0]:.2f}*I + {beta[1]:.2f}*M + {beta[2]:.2f}*J (EXACT!)")
    else:
        print(f"  M^2 ≠ aI + bM + cJ (resid={np.max(np.abs(resid)):.2f})")

# ========== DEEP: M(r) parameterized version ==========
print("\n\n=== M(r) parameterized transfer matrix ===")
# M(r)[a,b] should reduce to M[a,b] at r=1/2.
# Let's compute M(r) at several r values and see if it has a clean
# polynomial structure in r.

def compute_M_r(A, n, r):
    """Transfer matrix at parameter r: each path weighted by prod(r + s_e)"""
    M = [[0.0]*n for _ in range(n)]
    for P in permutations(range(n)):
        weight = 1.0
        for i in range(n-1):
            s = 0.5 if A[P[i]][P[i+1]] else -0.5
            weight *= (r + s)

        # How to attribute weight to matrix entry?
        # From the definition: M(r)[a,b] relates to paths FROM a TO b
        # weighted by the product of (r + s_e)
        # The diagonal: M(r)[a,a] = sum_P (-1)^{pos(a,P)} * prod(r+s_e)?
        # Wait, this isn't right. The r-parameterized version comes from
        # the W-polynomial being tr(M(r)).

        # Actually M(r) is defined so that sum_a M(r)[a,a] = W(T,r).
        # The off-diagonals contribute 0 to the trace.
        # But how are the off-diagonals defined?

        # From the formula:
        # M(r)[a,b] = sum_S (-1)^|S| * E_a(S,r) * B_b(R,r)
        # where E_a(S,r) and B_b(R,r) are r-weighted path counts.
        pass

    # Let me just compute it the brute force way for diagonal
    diag = [0.0] * n
    for a in range(n):
        for P in permutations(range(n)):
            weight = 1.0
            for i in range(n-1):
                s = 0.5 if A[P[i]][P[i+1]] else -0.5
                weight *= (r + s)
            pos_a = list(P).index(a)
            diag[a] += (-1)**pos_a * weight
    return diag

A = random_tournament(n)
print(f"  Diagonal of M(r) at several r values:")
for r in [0.0, 0.25, 0.5, 1.0, 1.5]:
    diag = compute_M_r(A, n, r)
    print(f"  r={r:.2f}: diag = {[f'{d:.2f}' for d in diag]}, sum={sum(diag):.2f}")

# At r=1/2: sum of diag = H(T)
# At r=0: sum of diag = c_0 * ... hmm
# W(r) = tr(M(r)) and W has only even powers of r
# So each M(r)[a,a] should also have only even powers of r?

# Let's check: compute M(r)[a,a] at r and -r
print(f"\n  M(r)[a,a] at r vs -r:")
for r in [0.3, 0.7, 1.5]:
    diag_pos = compute_M_r(A, n, r)
    diag_neg = compute_M_r(A, n, -r)
    print(f"  r={r}: diag(r) = {[f'{d:.2f}' for d in diag_pos]}")
    print(f"  r={-r}: diag(-r)= {[f'{d:.2f}' for d in diag_neg]}")
    print(f"   diff = {[f'{a-b:.2f}' for a,b in zip(diag_pos, diag_neg)]}")
    print(f"   sum  = {[f'{a+b:.2f}' for a,b in zip(diag_pos, diag_neg)]}")
