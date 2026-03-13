#!/usr/bin/env python3
"""
per_eigenspace_betti.py — opus-2026-03-13-S70

THM-154: Betti of circulant on Z_p is p * (per-eigenspace Betti).
THM-145: per-eigenspace Omega profile is determined by Q_k.

Now: compute per-eigenspace Betti DIRECTLY and see how they relate to Q_k.

For circulant on Z_p with connection set S:
- The chain complex decomposes over eigenspaces k=0,...,p-1
- Each eigenspace has its own Omega and Betti
- The k≠0 eigenspaces depend on Q_k = |Ŝ(k)|²
- For Paley: all Q_k equal => all eigenspaces identical
- For Interval: Q_k vary => eigenspaces differ

We can compute the per-eigenspace chain complex directly using the
DFT-block-diagonalization. For the k-th eigenspace:
- Ω_m^(k) = regular m-paths where sum of vertex labels ≡ 0 mod p
  (projected onto the k-th Fourier mode)
- Actually: the eigenspace decomposition of the boundary map
"""

import numpy as np
from itertools import permutations
from collections import defaultdict

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
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
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def boundary_matrix(paths_m, paths_m_minus_1):
    if not paths_m or not paths_m_minus_1:
        return np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)
    path_to_idx = {p: i for i, p in enumerate(paths_m_minus_1)}
    m = len(paths_m[0]) - 1
    matrix = np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in path_to_idx:
                matrix[path_to_idx[face], j] += sign
    return matrix

def per_eigenspace_decomposition(n, S):
    """Decompose the chain complex into eigenspaces of the Z_n rotation.

    For each eigenspace k, project the boundary matrix onto the k-th
    Fourier mode and compute the rank.
    """
    A = circulant_tournament(n, S)

    # Get all regular paths
    all_paths = {}
    for m in range(n):
        all_paths[m] = get_regular_paths(A, m)

    # For each path, compute its "total label mod n" for the rotation action.
    # Under rotation by 1: path (v_0,...,v_m) -> (v_0+1,...,v_m+1).
    # The eigenvalue for eigenspace k is ω^{k*(m+1)} where ω = e^{2πi/n}
    # since all vertices shift by 1.

    # Actually, the Z_n action on Ω_m is: σ(v_0,...,v_m) = (v_0+1,...,v_m+1) mod n.
    # An orbit of this action has size n (since p is prime, no fixed points for non-trivial rotation).
    # The projection onto eigenspace k is: P_k(f) = (1/n) Σ_{j=0}^{n-1} ω^{-jk} σ^j(f)

    # For the basis of orbits: each orbit {path, σ(path), ..., σ^{n-1}(path)}
    # contributes 1 to each eigenspace (each orbit gives a 1-dim representation space for each k).

    # Since the action is free and transitive on orbits of size n,
    # dim(Ω_m^(k)) = dim(Ω_m)/n for all k.

    # The boundary map in eigenspace k is: ∂_m^(k) is the projection of ∂_m.
    # Since ∂_m commutes with σ, the block-diagonalization works.

    omega = np.exp(2j * np.pi / n)
    results = {}

    for k in range(n):
        betti_k = []
        for m in range(n):
            # Group paths into orbits
            if not all_paths[m]:
                betti_k.append(0)
                continue

            # Each path (v_0,...,v_m) belongs to the orbit of its canonical form
            # (the one starting from 0, i.e., (v_0-v_0, v_1-v_0,...,v_m-v_0) mod n)
            seen = set()
            orbit_reps = []
            for path in all_paths[m]:
                # Canonical form: rotate so first element is 0
                shift = path[0]
                canonical = tuple((v - shift) % n for v in path)
                if canonical not in seen:
                    seen.add(canonical)
                    orbit_reps.append(canonical)

            # dim of eigenspace k = number of orbits
            # (since action is free and transitive on each orbit)

        # For the boundary map in eigenspace k, we need the Fourier transform
        # of the boundary matrix restricted to orbit representatives.

        # Actually, let's compute it properly by building the projected matrix.
        # For eigenspace k, the boundary matrix becomes:
        # ∂_m^(k)[orbit_β, orbit_α] = Σ_{j=0}^{n-1} ω^{jk} * ∂_m[σ^j(β), α]

        # Simpler approach: use orbit representatives starting from 0.
        # The boundary of a path starting from 0 might not start from 0,
        # but we can rotate it back and multiply by ω^k.

        omegas_k = []
        for m_dim in range(n):
            paths_m = [p for p in all_paths[m_dim] if p[0] == 0]
            omegas_k.append(len(paths_m))

        results[k] = {'omega': omegas_k}

    # Report
    return results, all_paths

# ============================================================
# Interval n=7: per-eigenspace analysis
# ============================================================
print("="*70)
print("PER-EIGENSPACE ANALYSIS: INTERVAL n=7")
print("="*70)

n = 7
S = {1, 2, 3}
A = circulant_tournament(n, S)

# Get orbits (paths starting from 0)
print("\n  Orbit counts (= per-eigenspace Omega):")
for m in range(n):
    all_p = get_regular_paths(A, m)
    from_0 = [p for p in all_p if p[0] == 0]
    print(f"    Omega_{m}: total={len(all_p)}, from-0={len(from_0)}, total/n={len(all_p)//n}")

# Compute per-eigenspace boundary matrices using paths from vertex 0
# For eigenspace k, the boundary operator on orbit representatives is:
# ∂^(k)(path) = Σ_{interior i} (-1)^i * ω^{k * face_shift} * face
# where face_shift rotates the face back to start from 0.

omega_root = np.exp(2j * np.pi / n)

print(f"\n  Per-eigenspace boundary ranks:")
per_eig_omega = []
per_eig_betti = []

for m in range(n):
    paths_from_0 = [p for p in get_regular_paths(A, m) if p[0] == 0]
    per_eig_omega.append(len(paths_from_0))

for k in range(n):
    ranks = {}
    for m in range(1, n):
        paths_m = [p for p in get_regular_paths(A, m) if p[0] == 0]
        paths_m1 = [p for p in get_regular_paths(A, m-1) if p[0] == 0]

        if not paths_m or not paths_m1:
            ranks[m] = 0
            continue

        # Build the projected boundary matrix for eigenspace k
        path_to_idx = {p: i for i, p in enumerate(paths_m1)}
        B_k = np.zeros((len(paths_m1), len(paths_m)), dtype=complex)

        path_len = len(paths_m[0]) - 1  # = m

        for j, path in enumerate(paths_m):
            for i in range(1, path_len):  # interior vertices
                face = path[:i] + path[i+1:]
                sign = (-1)**i

                # Rotate face so it starts from 0
                shift = face[0]
                canonical_face = tuple((v - shift) % n for v in face)
                phase = omega_root ** (k * shift)

                if canonical_face in path_to_idx:
                    B_k[path_to_idx[canonical_face], j] += sign * phase

        ranks[m] = np.linalg.matrix_rank(B_k, tol=1e-8)

    # Compute Betti for this eigenspace
    betti_k = []
    for m in range(n):
        omega_m = per_eig_omega[m]
        rk_dm = ranks.get(m, 0)
        rk_dm1 = ranks.get(m+1, 0)
        beta = omega_m - rk_dm - rk_dm1
        betti_k.append(beta)

    print(f"  k={k}: ranks={[ranks.get(m,0) for m in range(1,n)]}, "
          f"beta={betti_k}")

    per_eig_betti.append(betti_k)

# Check: sum over k should give total Betti
print(f"\n  Total β (sum over k): {[sum(per_eig_betti[k][m] for k in range(n)) for m in range(n)]}")

# Compute Q_k
print(f"\n  Q_k values:")
for k in range(n):
    if k == 0:
        Q_k = len(S)**2  # |Ŝ(0)|² = |S|²
        print(f"    Q_0 = {Q_k} (trivial)")
    else:
        S_hat = sum(omega_root**(k*s) for s in S)
        Q_k = abs(S_hat)**2
        print(f"    Q_{k} = {Q_k:.4f}")

# ============================================================
# Paley p=7: same analysis
# ============================================================
print(f"\n{'='*70}")
print("PER-EIGENSPACE ANALYSIS: PALEY p=7")
print("="*70)

S_pal = {1, 2, 4}
A_pal = circulant_tournament(n, S_pal)

per_eig_omega_pal = []
for m in range(n):
    paths_from_0 = [p for p in get_regular_paths(A_pal, m) if p[0] == 0]
    per_eig_omega_pal.append(len(paths_from_0))

print(f"  Per-eigenspace Omega: {per_eig_omega_pal}")

per_eig_betti_pal = []
for k in range(n):
    ranks = {}
    for m in range(1, n):
        paths_m = [p for p in get_regular_paths(A_pal, m) if p[0] == 0]
        paths_m1 = [p for p in get_regular_paths(A_pal, m-1) if p[0] == 0]

        if not paths_m or not paths_m1:
            ranks[m] = 0
            continue

        path_to_idx = {p: i for i, p in enumerate(paths_m1)}
        B_k = np.zeros((len(paths_m1), len(paths_m)), dtype=complex)
        path_len = len(paths_m[0]) - 1

        for j, path in enumerate(paths_m):
            for i in range(1, path_len):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                shift = face[0]
                canonical_face = tuple((v - shift) % n for v in face)
                phase = omega_root ** (k * shift)
                if canonical_face in path_to_idx:
                    B_k[path_to_idx[canonical_face], j] += sign * phase

        ranks[m] = np.linalg.matrix_rank(B_k, tol=1e-8)

    betti_k = []
    for m in range(n):
        omega_m = per_eig_omega_pal[m]
        rk_dm = ranks.get(m, 0)
        rk_dm1 = ranks.get(m+1, 0)
        beta = omega_m - rk_dm - rk_dm1
        betti_k.append(beta)

    print(f"  k={k}: beta={betti_k}")
    per_eig_betti_pal.append(betti_k)

# Verify: all k≠0 eigenspaces should be identical for Paley
print(f"\n  All k≠0 eigenspaces identical?")
ref = per_eig_betti_pal[1]
all_same = all(per_eig_betti_pal[k] == ref for k in range(1, n))
print(f"    {all_same}")
if all_same:
    print(f"    Per-eigenspace Betti (k≠0): {ref}")
    print(f"    k=0 eigenspace Betti: {per_eig_betti_pal[0]}")

print("\nDONE.")
