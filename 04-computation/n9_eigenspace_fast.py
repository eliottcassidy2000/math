#!/usr/bin/env python3
"""
n9_eigenspace_fast.py - Eigenspace decomposition of n=9 circulant maximizer b5=10

Uses GLMY boundary: d(v0,...,vp) = sum_{i=0}^{p} (-1)^i (v0,...,vi_hat,...,vp)
Omega_p = { u in A_p : du in A_{p-1} }

For the circulant maximizer on Z_9 with S={1,5,6,7}, we have b5=10.
Question: how does Z/9Z decompose this across eigenspaces?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def build_circulant(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def enumerate_paths(A, n, length):
    """Enumerate all allowed p-paths (sequences with directed edges)."""
    if length == 0:
        return [(v,) for v in range(n)]
    paths = []
    def dfs(path, used):
        if len(path) == length + 1:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v not in used and A[last][v]:
                dfs(path + [v], used | {v})
    for s in range(n):
        dfs([s], {s})
    return paths

def boundary_faces(path):
    """All faces of a path with signs: [(sign, face), ...]"""
    p = len(path) - 1
    return [((-1)**i, path[:i] + path[i+1:]) for i in range(p + 1)]

def compute_omega(A, n, all_paths, p):
    """Compute Omega_p = { u in A_p : du in A_{p-1} }.
    For p=0: Omega_0 = A_0.
    For p>=1: need all faces of u (including terminal) to be in A_{p-1}."""
    if p == 0:
        return all_paths[0]

    allowed_pm1 = set(all_paths[p-1])
    regular = []
    for path in all_paths[p]:
        is_omega = True
        for sign, face in boundary_faces(path):
            if len(set(face)) == len(face) and face not in allowed_pm1:
                is_omega = False
                break
        if is_omega:
            regular.append(path)
    return regular

def build_boundary(omega_p, omega_pm1):
    """Build GLMY boundary matrix d_p: Omega_p -> Omega_{p-1}.
    Uses FULL boundary: d = sum_{i=0}^{p} (-1)^i face_i."""
    pm1_idx = {tuple(p): i for i, p in enumerate(omega_pm1)}
    m_p = len(omega_p)
    m_pm1 = len(omega_pm1)

    bd = np.zeros((m_pm1, m_p), dtype=float)
    for j, path in enumerate(omega_p):
        for sign, face in boundary_faces(path):
            if face in pm1_idx:
                bd[pm1_idx[face], j] += sign
    return bd

# ===== Main =====
print("=" * 70)
print("EIGENSPACE DECOMPOSITION: n=9 CIRCULANT MAXIMIZER")
print("=" * 70)

n = 9
S = {1, 5, 6, 7}
A = build_circulant(n, S)
H = H_tournament(A, n)
print(f"Z_{n}, S={sorted(S)}, H={H}")
assert H == 3357

t0 = time.time()

# Enumerate allowed paths
print("\n--- Allowed paths ---")
all_paths = {}
for d in range(n):
    t1 = time.time()
    paths = enumerate_paths(A, n, d)
    all_paths[d] = paths
    if not paths:
        print(f"  dim {d}: 0 paths")
        break
    print(f"  dim {d}: {len(paths)} ({time.time()-t1:.1f}s)")

# Compute Omega
print("\n--- Omega (regular) paths ---")
omega = {}
for d in range(n):
    if d not in all_paths or not all_paths[d]:
        break
    om = compute_omega(A, n, all_paths, d)
    omega[d] = om
    print(f"  Om_{d}: {len(om)} paths")

# Build boundary matrices
print("\n--- Boundary matrices ---")
boundaries = {}
for d in range(1, max(omega.keys()) + 1):
    if d not in omega or d-1 not in omega:
        continue
    bd = build_boundary(omega[d], omega[d-1])
    rank = np.linalg.matrix_rank(bd, tol=1e-8)
    boundaries[d] = bd
    print(f"  d_{d}: {bd.shape[0]}x{bd.shape[1]}, rank={rank}")

# Compute global Betti
print("\n--- Global Betti ---")
max_d = max(omega.keys())
betti = []
for d in range(max_d + 1):
    m_d = len(omega[d])
    # ker(d_d)
    if d not in boundaries:
        ker = m_d  # d_0 = 0
    else:
        ker = m_d - np.linalg.matrix_rank(boundaries[d], tol=1e-8)
    # im(d_{d+1})
    if d+1 not in boundaries:
        im = 0
    else:
        im = np.linalg.matrix_rank(boundaries[d+1], tol=1e-8)
    betti.append(ker - im)
print(f"  b = {betti}")

# ===== Eigenspace decomposition =====
print("\n" + "=" * 70)
print("EIGENSPACE BETTI NUMBERS (Z/9Z)")
print("=" * 70)

omega_9 = np.exp(2j * np.pi / 9)

for k in range(n):
    betti_k = []
    for d in range(max_d + 1):
        m_d = len(omega[d])

        # Build shift permutation on Omega_d
        p_idx = {tuple(p): i for i, p in enumerate(omega[d])}
        sigma = np.zeros((m_d, m_d), dtype=complex)
        for j, p in enumerate(omega[d]):
            shifted = tuple((v + 1) % n for v in p)
            if shifted in p_idx:
                sigma[p_idx[shifted], j] = 1.0

        # Projector onto eigenspace k
        P = np.zeros((m_d, m_d), dtype=complex)
        sigma_j = np.eye(m_d, dtype=complex)
        for jj in range(n):
            P += (omega_9**(-k*jj)) * sigma_j
            sigma_j = sigma_j @ sigma
        P /= n
        proj_dim = round(np.real(np.trace(P)))

        if proj_dim == 0:
            betti_k.append(0)
            continue

        # ker(d_d) in eigenspace k
        if d not in boundaries:
            ker_k = proj_dim
        else:
            bd_proj = boundaries[d].astype(complex) @ P
            ker_k = proj_dim - np.linalg.matrix_rank(bd_proj, tol=1e-6)

        # im(d_{d+1}) in eigenspace k
        if d+1 not in boundaries:
            im_k = 0
        else:
            # Need projector on Omega_{d+1}
            m_dp1 = len(omega[d+1])
            p_idx_dp1 = {tuple(p): i for i, p in enumerate(omega[d+1])}
            sigma_dp1 = np.zeros((m_dp1, m_dp1), dtype=complex)
            for j, p in enumerate(omega[d+1]):
                shifted = tuple((v + 1) % n for v in p)
                if shifted in p_idx_dp1:
                    sigma_dp1[p_idx_dp1[shifted], j] = 1.0
            P_dp1 = np.zeros((m_dp1, m_dp1), dtype=complex)
            sigma_j = np.eye(m_dp1, dtype=complex)
            for jj in range(n):
                P_dp1 += (omega_9**(-k*jj)) * sigma_j
                sigma_j = sigma_j @ sigma_dp1
            P_dp1 /= n

            bd_next_proj = boundaries[d+1].astype(complex) @ P_dp1
            im_k = np.linalg.matrix_rank(bd_next_proj, tol=1e-6)

        betti_k.append(max(0, int(round(ker_k - im_k))))

    if any(b > 0 for b in betti_k):
        # Get eigenspace dimensions
        dims_k = []
        for d in range(max_d + 1):
            m_d = len(omega[d])
            p_idx_d = {tuple(p): i for i, p in enumerate(omega[d])}
            sigma_d = np.zeros((m_d, m_d), dtype=complex)
            for j, p in enumerate(omega[d]):
                shifted = tuple((v + 1) % n for v in p)
                if shifted in p_idx_d:
                    sigma_d[p_idx_d[shifted], j] = 1.0
            P_d = np.zeros((m_d, m_d), dtype=complex)
            sigma_j = np.eye(m_d, dtype=complex)
            for jj in range(n):
                P_d += (omega_9**(-k*jj)) * sigma_j
                sigma_j = sigma_j @ sigma_d
            P_d /= n
            dims_k.append(round(np.real(np.trace(P_d))))
        print(f"  k={k}: Om_dims={dims_k}, b={betti_k}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("Expected: b = [1,0,0,0,0,10,0,0,0]")
print(f"Got:      b = {betti}")
print(f"\nFor P_7: 6 non-trivial eigenspaces each contribute b4=1.")
print(f"Question: How does Z_9 decompose b5=10?")
print(f"\nTotal time: {time.time()-t0:.1f}s")
