#!/usr/bin/env python3
"""
HYPOTHESIS CC: Transfer matrix symmetry from deletion-contraction

The other agent verified: H(T) = H(T\e) + H(T/e)
where e is an edge, \e deletes it, /e contracts it.

For the TRANSFER MATRIX version:
Does M(T) = M_\e(T\e) + M_/e(T/e) in some sense?

If so, and if M(T\e) and M(T/e) are both symmetric (by induction),
then M(T) would be symmetric.

The issue: T\e is NOT a tournament (it has a non-edge).
T/e IS a tournament on n-1 vertices.

So the induction would be:
- M(T/e) is symmetric by induction hypothesis
- M(T\e) needs a separate argument

Actually, the deletion-contraction H(T) = H(T\e) + H(T/e) gives us
a recursion on H. Can we lift it to M?

Let me first verify the deletion-contraction formula and understand
the contraction operation precisely.

Convention (from the other agent):
"contraction merges tail/head, w inherits IN from tail, OUT from head"

So for edge u→v, contracting:
- Remove u and v, add new vertex w
- For any other vertex x:
  - w→x iff v→x (inherit OUT from head v)
  - x→w iff x→u (inherit IN from tail u)
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

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def delete_edge(A, u, v):
    """Delete edge u→v: set A[u][v]=0, A[v][u]=0 (no longer a tournament!)"""
    n = len(A)
    B = [row[:] for row in A]
    B[u][v] = 0
    # Don't set B[v][u] = 1! The edge is just removed.
    return B

def contract_edge(A, u, v):
    """Contract edge u→v: merge u,v into w.
    w inherits IN from tail u, OUT from head v.
    """
    n = len(A)
    others = [x for x in range(n) if x != u and x != v]
    k = len(others)
    B = [[0]*(k+1) for _ in range(k+1)]
    # w is at index 0
    # others are at indices 1..k
    for i, x in enumerate(others):
        for j, y in enumerate(others):
            B[i+1][j+1] = A[x][y]
        # w→x iff v→x (head's outgoing)
        B[0][i+1] = A[v][x]
        # x→w iff x→u (tail's incoming)
        B[i+1][0] = A[x][u]
    return B, [None] + others  # vertex map

def count_HP_general(A, n):
    """Count Hamiltonian paths in a general digraph (not necessarily tournament)"""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

random.seed(42)

# Verify H(T) = H(T\e) + H(T/e)
print("=== Verify deletion-contraction H(T) = H(T\\e) + H(T/e) ===")
for n in [4, 5, 6]:
    violations = 0
    tests = 0
    for trial in range(20 if n < 6 else 5):
        A = random_tournament(n)
        H_T = ham_path_count_dp(A, n)

        for u in range(n):
            for v in range(n):
                if u == v or not A[u][v]: continue
                tests += 1

                # Deletion: remove edge u→v
                B_del = delete_edge(A, u, v)
                H_del = count_HP_general(B_del, n)

                # Contraction: merge u,v
                B_con, vmap = contract_edge(A, u, v)
                H_con = ham_path_count_dp(B_con, n-1)

                if H_T != H_del + H_con:
                    violations += 1
                    if violations <= 3:
                        print(f"  VIOLATION at n={n}: H={H_T}, H_del={H_del}, H_con={H_con}, "
                              f"edge {u}→{v}")

    print(f"  n={n}: {violations}/{tests} violations")

# ========== Now test: does M satisfy deletion-contraction? ==========
print("\n\n=== Does M satisfy deletion-contraction? ===")

def compute_M_small(A, n):
    """Transfer matrix for tournament or general digraph"""
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
    return np.array(M, dtype=float)

n = 4
for trial in range(3):
    A = random_tournament(n)
    M_T = compute_M_small(A, n)

    # Pick an edge
    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]: continue

            # M for T\e (deletion)
            B_del = delete_edge(A, u, v)
            M_del = compute_M_small(B_del, n)

            # M for T/e (contraction)
            B_con, vmap = contract_edge(A, u, v)
            M_con = compute_M_small(B_con, n-1)

            # The contraction matrix is (n-1)×(n-1).
            # To compare with M_T (n×n), we need to embed M_con.
            # The contracted vertex w represents {u,v}.
            # How should M_con contribute to M_T?

            # For diagonal: M_T[a,a] = M_del[a,a] + ??
            # Let's just check if tr(M_T) = tr(M_del) + tr(M_con)
            tr_T = np.trace(M_T)
            tr_del = np.trace(M_del)
            tr_con = np.trace(M_con)

            if abs(tr_T - tr_del - tr_con) < 0.01:
                pass  # expected since H = H_del + H_con and tr = H (odd n) or 0 (even n)
            break
        break

    # At n=4 (even): tr(M_T) = 0, H = H_del + H_con
    # tr(M_del) for the non-tournament may not be 0
    # tr(M_con) for the (n-1=3)-tournament = H_con (odd)

    H_T = ham_path_count_dp(A, n)
    u, v = 0, [x for x in range(n) if A[0][x]][0]
    B_del = delete_edge(A, u, v)
    B_con, vmap = contract_edge(A, u, v)
    M_del = compute_M_small(B_del, n)
    M_con = compute_M_small(B_con, n-1)

    print(f"\n  trial {trial}: n={n}, edge {u}→{v}")
    print(f"  H_T={H_T}, H_del={count_HP_general(B_del,n)}, H_con={ham_path_count_dp(B_con,n-1)}")
    print(f"  tr(M_T)={np.trace(M_T):.0f}, tr(M_del)={np.trace(M_del):.0f}, tr(M_con)={np.trace(M_con):.0f}")
    print(f"  M_T symmetric: {np.allclose(M_T, M_T.T)}")
    print(f"  M_del symmetric: {np.allclose(M_del, M_del.T)}")
    print(f"  M_con symmetric: {np.allclose(M_con, M_con.T)}")
