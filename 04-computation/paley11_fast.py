#!/usr/bin/env python3
"""
paley11_fast.py — opus-2026-03-13-S71

Fast GLMY Ω dim computation for Paley P_11 using circulant structure.

Key insight: For a circulant tournament on Z_p, vertex-transitivity means
we can compute Ω_m by:
1. Enumerate paths from vertex 0 only (factor of p speedup)
2. Use the circulant structure of the junk map

For directed m-paths from vertex 0, there are |A_m|/p of them.
The junk map J restricted to these paths gives rank(J_0).
Then dim(Ω_m) = p * (|A_m|/p - rank(J_0)).

Wait — this is NOT quite right because the junk map involves faces
from ALL paths, not just those starting at 0. Need to be more careful.

Actually: the Z_p action on paths gives a representation.
The junk map J commutes with this action.
By Schur's lemma, J block-diagonalizes into p eigenspaces.
Each eigenspace has the same rank (for a prime p and transitive action).
So rank(J) = p * rank(J restricted to eigenspace k=0).

The k=0 eigenspace corresponds to Z_p-invariant combinations.
For paths starting at 0, the k=0 projection is (1/p) * sum over shifts.
But we can work with paths starting at 0 directly!

Let me just enumerate paths from 0 and count them first.
"""

import numpy as np
import time

def circulant_adj(p, S):
    """Returns adjacency list for circulant tournament."""
    adj = [[] for _ in range(p)]
    for i in range(p):
        for s in S:
            adj[i].append((i + s) % p)
    return adj

# P_11: QR = {1, 3, 4, 5, 9}
p = 11
QR = {a*a % p for a in range(1, p)}
print(f"P_{p}: QR = {sorted(QR)}")
adj = circulant_adj(p, QR)

# Enumerate paths from vertex 0
def paths_from_0(adj, n, m):
    """Enumerate directed m-paths starting from vertex 0."""
    if m == 0: return [(0,)]
    paths = []
    def dfs(path, depth, visited):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in adj[last]:
            if v not in visited:
                path.append(v)
                visited.add(v)
                dfs(path, depth+1, visited)
                visited.remove(v)
                path.pop()
    dfs([0], 0, {0})
    return paths

print(f"\nDirected path counts from vertex 0:")
path_counts = {}
all_paths_from0 = {}
for m in range(p):
    t0 = time.time()
    paths = paths_from_0(adj, p, m)
    t1 = time.time()
    path_counts[m] = len(paths)
    all_paths_from0[m] = paths
    total = len(paths) * p
    print(f"  m={m:2d}: from_0={len(paths):8d}, total={total:10d}, time={t1-t0:.2f}s")
    if t1 - t0 > 60:
        print(f"  (stopping enumeration — too slow)")
        break

# ============================================================
# Now compute Ω dims using the junk map on paths from 0
# ============================================================
print(f"\n{'='*70}")
print("OMEGA DIM COMPUTATION")
print("="*70)

# For paths starting from 0: faces can start from any vertex.
# A face of (0, v_1, ..., v_m) by deleting v_i gives a (m-1)-path.
# The face starts at 0 if i > 0, or at v_1 if i = 0.
#
# For the junk map, a face is "junk" if it's not a directed (m-1)-path.
# Due to vertex transitivity: a (m-1)-tuple starting at 0 is a directed path
# iff the corresponding shifted tuple starting at any vertex is.
#
# So we can classify ALL (m-1)-tuples as junk or not by checking if they're
# directed paths in the tournament.
#
# For the FULL junk map J (all paths), J has a block-circulant structure.
# rank(J) = sum over eigenspaces of rank(J_k)
# By the transitive Z_p action on a prime, J_k has the same rank for all k≠0
# (or possibly differs for k=0 vs k≠0).
#
# SIMPLER approach: just compute for small m where we have all paths.

def is_directed_path(A_adj, path):
    """Check if a tuple is a directed path (all edges present, all distinct)."""
    if len(set(path)) != len(path): return False
    for i in range(len(path) - 1):
        if path[i+1] not in A_adj[path[i]]:
            return False
    return True

# Build adjacency matrix for fast lookup
A = np.zeros((p,p), dtype=int)
for i in range(p):
    for s in QR:
        A[i][(i+s)%p] = 1

# Compute Ω for small m using ALL paths
max_m_full = min(5, max(m for m in all_paths_from0))

for m in range(max_m_full + 1):
    t0 = time.time()

    # All m-paths (using vertex transitivity)
    paths_0 = all_paths_from0[m]
    # All m-paths total
    all_m_paths = []
    for start in range(p):
        for path_0 in paths_0:
            shifted = tuple((v + start) % p for v in path_0)
            all_m_paths.append(shifted)

    # All (m-1)-paths
    if m >= 1:
        paths_m1_0 = all_paths_from0[m-1]
        all_m1_paths = set()
        for start in range(p):
            for path_0 in paths_m1_0:
                shifted = tuple((v + start) % p for v in path_0)
                all_m1_paths.add(shifted)
    else:
        all_m1_paths = set()

    if m <= 1:
        omega_m = len(all_m_paths)
        print(f"  m={m}: Ω_m = {omega_m}, Ω_m/{p} = {omega_m//p}")
        continue

    # Build junk map
    junk = {}
    junk_count = 0
    for path in all_m_paths:
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face not in all_m1_paths and face not in junk:
                junk[face] = junk_count
                junk_count += 1

    if junk_count == 0:
        omega_m = len(all_m_paths)
        print(f"  m={m}: no junk, Ω_m = {omega_m}, Ω_m/{p} = {omega_m//p}")
        continue

    # Build sparse J matrix
    from scipy.sparse import lil_matrix
    J = lil_matrix((junk_count, len(all_m_paths)))
    for j, path in enumerate(all_m_paths):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in junk:
                J[junk[face], j] += (-1)**i

    J_csr = J.tocsr()

    # Compute rank via J^H J (smaller if one dim is small)
    if junk_count < len(all_m_paths):
        # J^H J is junk_count × junk_count  — use JJ^T
        M = (J_csr @ J_csr.T).toarray()
        rank_J = np.linalg.matrix_rank(M, tol=1e-8)
    else:
        # J^T J is |A_m| × |A_m|
        M = (J_csr.T @ J_csr).toarray()
        rank_J = np.linalg.matrix_rank(M, tol=1e-8)

    omega_m = len(all_m_paths) - rank_J
    t1 = time.time()
    print(f"  m={m}: |A_m|={len(all_m_paths)}, junk={junk_count}, rank(J)={rank_J}, "
          f"Ω_m={omega_m}, Ω_m/{p}={omega_m//p if omega_m%p==0 else f'{omega_m}/{p}'}, "
          f"time={t1-t0:.1f}s")

print("\nDONE.")
