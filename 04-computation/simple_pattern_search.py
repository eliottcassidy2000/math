#!/usr/bin/env python3
"""
STEPPING BACK: What is the SIMPLE underlying pattern?

The user says "the underlying patterns here can't be too complex."
Let me look for the simplest possible characterization of M[a,b].

Key known facts:
1. M is symmetric (proved computationally)
2. tr(M) = H(T) for odd n
3. M[a,a] = sum_P (-1)^{pos(a,P)}
4. For vertex-transitive T: M = (H/n)*I (scalar)
5. W(T,r) = W(T^op,r) for odd n

The SIMPLEST hypothesis: M is a function of the "local neighborhood structure"
around vertices a and b.

HYPOTHESIS V: M[a,b] depends only on the INDUCED TOURNAMENT on N[a] ∩ N[b]
(common neighbors) and the edge direction between a and b.

HYPOTHESIS W: M[a,b] = (-1)^{I(a→b)} * f(something simple)
where f counts some kind of signed walk from a to b.

HYPOTHESIS X: M is related to the COFACTOR MATRIX of (I - A) or similar.
For an n×n matrix, the cofactor matrix C satisfies A * C^T = det(A) * I.
If det(I-A) ≈ H in some sense, maybe M ≈ cofactor(I-A)?

Let me check this numerically.
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

print("=== HYPOTHESIS X: M vs cofactor(I-A) ===")
n = 5
for trial in range(5):
    A = random_tournament(n)
    M = np.array(compute_M(A, n), dtype=float)
    A_np = np.array(A, dtype=float)

    # Various matrix constructions
    IminusA = np.eye(n) - A_np
    detIA = np.linalg.det(IminusA)

    # Cofactor matrix = det * inverse^T (when invertible)
    if abs(detIA) > 0.01:
        cofactor = detIA * np.linalg.inv(IminusA).T
    else:
        cofactor = np.zeros((n,n))

    # Adjugate of A
    adj_A = np.linalg.det(A_np) * np.linalg.inv(A_np).T if np.linalg.det(A_np) != 0 else np.zeros((n,n))

    # (I-A)^{-1} type
    if abs(detIA) > 0.01:
        inv_IA = np.linalg.inv(IminusA)
    else:
        inv_IA = np.zeros((n,n))

    H = int(np.trace(M))

    print(f"\n  trial {trial}: H={H}, det(I-A)={detIA:.2f}")
    print(f"  M = \n{M.astype(int)}")

    # Try various linear combinations
    for name, candidate in [
        ("cofactor(I-A)", cofactor),
        ("(I-A)^{-1}", inv_IA),
        ("A^2 - (n-1)/2 * A", A_np@A_np - (n-1)/2 * A_np),
    ]:
        if np.max(np.abs(candidate)) < 0.01:
            print(f"  {name}: all zeros, skip")
            continue
        # Is M proportional to candidate?
        ratios = []
        for i in range(n):
            for j in range(n):
                if abs(candidate[i,j]) > 0.01:
                    ratios.append(M[i,j] / candidate[i,j])
        if ratios:
            if max(ratios) - min(ratios) < 0.01:
                print(f"  {name}: M = {ratios[0]:.4f} * candidate (PROPORTIONAL!)")
            else:
                print(f"  {name}: ratio range = [{min(ratios):.2f}, {max(ratios):.2f}] (not proportional)")

    # Try: M = a*I + b*A + c*A^T + d*A^2 + ...
    # Flatten M and fit as linear combination of matrix powers of A
    terms = [np.eye(n), A_np, A_np.T, A_np@A_np, A_np.T@A_np.T,
             A_np@A_np.T, A_np.T@A_np]
    X = np.column_stack([t.flatten() for t in terms])
    beta = np.linalg.lstsq(X, M.flatten(), rcond=None)[0]
    resid = M.flatten() - X @ beta
    max_resid = np.max(np.abs(resid))
    if max_resid < 0.01:
        print(f"  M = {beta[0]:.2f}*I + {beta[1]:.2f}*A + {beta[2]:.2f}*A^T + "
              f"{beta[3]:.2f}*A^2 + {beta[4]:.2f}*(A^T)^2 + "
              f"{beta[5]:.2f}*A*A^T + {beta[6]:.2f}*A^T*A (EXACT)")
    else:
        print(f"  Linear combination of A powers: max_resid={max_resid:.2f} (not exact)")

# ========== HYPOTHESIS Y: M[a,b] via signed walks ==========
print("\n\n=== HYPOTHESIS Y: M[a,b] as signed walk count ===")
# M[a,a] = sum_P (-1)^{pos(a,P)}
# For off-diagonal: is there a signed walk interpretation?
#
# M[a,b] = sum_S (-1)^|S| E_a(S) B_b(R)
# The alternating sign creates an inclusion-exclusion.
# What if M[a,b] counts "signed paths from a to b" in some signed graph?

# At n=3: let's check ALL tournaments
n = 3
edges = [(0,1), (0,2), (1,2)]
print(f"\n  ALL tournaments at n=3:")
for bits in range(8):
    A = [[0]*3 for _ in range(3)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    M = compute_M(A, n)
    H = sum(1 for P in permutations(range(n))
            if all(A[P[i]][P[i+1]] for i in range(n-1)))

    edge_str = "".join(f"{i}→{j}" if A[i][j] else f"{j}→{i}" for i,j in edges)
    print(f"  {edge_str}: H={H}, M={M}")

# ========== HYPOTHESIS Z: M as Möbius function of tournament poset ==========
print("\n\n=== HYPOTHESIS Z: Simplest characterization ===")
# At n=3, there are only 2 tournament types:
# Transitive: H=1, and cyclic: H=3
# For transitive (0→1→2→0... wait, 0→1, 1→2, 0→2):
# M should be easy to compute and understand.

# Actually, the simplest characterization might be:
# M[a,b] = #{HPs where a appears BEFORE b, weighted by (-1)^{pos(a)}}
#        - #{HPs where a appears AFTER b, weighted by (-1)^{pos(a)}}
#
# Or more precisely:
# M[a,b] for a≠b: this is some inclusion-exclusion count.
#
# Let me check: is M[a,b] related to the POSITION CORRELATION between a and b?
print("  Checking position correlation interpretation...")
n = 5
for trial in range(3):
    A = random_tournament(n)
    M = np.array(compute_M(A, n), dtype=float)

    # For each HP P, compute (-1)^{pos(a,P)} * (-1)^{pos(b,P)}
    C = [[0]*n for _ in range(n)]
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        for a in range(n):
            for b in range(n):
                pos_a = P.index(a)
                pos_b = P.index(b)
                C[a][b] += (-1)**pos_a * (-1)**pos_b

    C = np.array(C, dtype=float)
    print(f"\n  trial {trial}:")
    print(f"  M = \n{M.astype(int)}")
    print(f"  C (signed position product) = \n{C.astype(int)}")
    print(f"  M == C? {np.allclose(M, C)}")
