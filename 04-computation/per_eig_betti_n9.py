#!/usr/bin/env python3
"""
per_eig_betti_n9.py — opus-2026-03-13-S70

Verify eigenspace Betti uniformity at n=9 for the Interval tournament.
Q_k values range from 0.28 to 8.29 — very non-uniform.
Do all eigenspaces still have identical Betti?
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths_from_0(A, m):
    """Get regular m-paths starting from vertex 0."""
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    dfs([0], 0, -1)
    return paths

n = 9
S = {1, 2, 3, 4}
A = circulant_tournament(n, S)
omega_root = np.exp(2j * np.pi / n)

print(f"Interval n={n}, S={sorted(S)}")
print(f"Q_k values:")
for k in range(n):
    if k == 0:
        print(f"  Q_0 = {len(S)**2} (trivial)")
    else:
        S_hat = sum(omega_root**(k*s) for s in S)
        print(f"  Q_{k} = {abs(S_hat)**2:.4f}")

# Compute per-eigenspace Betti for select eigenspaces
per_eig_omega = []
for m in range(n):
    paths = get_regular_paths_from_0(A, m)
    per_eig_omega.append(len(paths))
    print(f"Omega_{m}: total={len(paths)*n}, per-eig={len(paths)}")

# Build and analyze boundary matrices for k=0 and k=1
for k_test in [0, 1, 2]:
    print(f"\nEigenspace k={k_test}:")
    ranks = {}
    for m in range(1, n):
        paths_m = get_regular_paths_from_0(A, m)
        paths_m1 = get_regular_paths_from_0(A, m-1)

        if not paths_m or not paths_m1:
            ranks[m] = 0
            continue

        path_to_idx = {p: i for i, p in enumerate(paths_m1)}
        B_k = np.zeros((len(paths_m1), len(paths_m)), dtype=complex)
        path_len = m

        for j, path in enumerate(paths_m):
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                shift = face[0]
                canonical_face = tuple((v - shift) % n for v in face)
                phase = omega_root ** (k_test * shift)
                if canonical_face in path_to_idx:
                    B_k[path_to_idx[canonical_face], j] += sign * phase

        ranks[m] = np.linalg.matrix_rank(B_k, tol=1e-8)

    betti_k = []
    for m in range(n):
        omega_m = per_eig_omega[m]
        rk_dm = ranks.get(m, 0)
        rk_dm1 = ranks.get(m+1, 0)
        beta = omega_m - rk_dm - rk_dm1
        betti_k.append(beta)

    print(f"  ranks: {[ranks.get(m,0) for m in range(1,n)]}")
    print(f"  beta:  {betti_k}")

print("\nDONE.")
