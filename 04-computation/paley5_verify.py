#!/usr/bin/env python3
"""
paley5_verify.py — opus-2026-03-13-S71

Double-check Paley P_5 = C5_{1,4} GLMY Betti.
β = [1, 6, 0, 0, 0] seems surprisingly large (β_1=6 from Ω_1=10).
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

# P_5 = C5_{1,4}: i→j iff j-i ∈ {1,4} mod 5
A = circulant_tournament(5, {1,4})
print("P_5 adjacency matrix:")
print(A)

print("\nOut-degrees:", [sum(A[i]) for i in range(5)])
# Score should be (2,2,2,2,2) since |{1,4}|=2

# Enumerate directed 1-paths (just edges)
edges = []
for i in range(5):
    for j in range(5):
        if A[i][j]:
            edges.append((i,j))
print(f"\nEdges ({len(edges)}):")
for e in edges:
    print(f"  {e[0]}→{e[1]}")

# Enumerate directed 2-paths
paths2 = []
for i in range(5):
    for j in range(5):
        if i == j or not A[i][j]: continue
        for k in range(5):
            if k == i or k == j or not A[j][k]: continue
            paths2.append((i,j,k))

print(f"\n2-paths: {len(paths2)}")

# Check which 2-path faces are NOT edges
print(f"\nFaces of 2-paths that are NOT directed 1-paths:")
edge_set = set(edges)
for p in paths2[:10]:
    faces = []
    for i in range(3):
        face = p[:i] + p[i+1:]
        is_edge = face in edge_set
        faces.append(f"{face}{'✓' if is_edge else '✗'}")
    print(f"  {p}: faces = {faces}")

# Build junk map for 2-paths
junk_faces = {}
junk_count = 0
for path in paths2:
    for i in range(3):
        face = path[:i] + path[i+1:]
        if face not in edge_set and face not in junk_faces:
            junk_faces[face] = junk_count
            junk_count += 1

print(f"\nJunk faces (non-edges appearing as faces of 2-paths): {junk_count}")
for face, idx in sorted(junk_faces.items(), key=lambda x: x[1]):
    print(f"  {face}")

# Build J matrix
J = np.zeros((junk_count, len(paths2)))
for j, path in enumerate(paths2):
    for i in range(3):
        face = path[:i] + path[i+1:]
        if face in junk_faces:
            J[junk_faces[face], j] += (-1)**i

print(f"\nJ matrix shape: {J.shape}")
U, s, Vh = np.linalg.svd(J, full_matrices=True)
rank_J = np.sum(s > 1e-10)
print(f"rank(J) = {rank_J}")
print(f"dim(Ω_2) = {len(paths2)} - {rank_J} = {len(paths2) - rank_J}")

# So Ω_2 = null(J)
omega2_dim = len(paths2) - rank_J

# Now build boundary ∂_1: Ω_1 → Ω_0
# Ω_0 = 5 vertices, Ω_1 = edges (all 10 of them)
# ∂_1(i,j) = j - i (in GLMY full boundary)

B1 = np.zeros((5, len(edges)))
for j, (a,b) in enumerate(edges):
    B1[b, j] += 1   # face from deleting index 0
    B1[a, j] -= 1   # face from deleting index 1, sign (-1)^1

print(f"\n∂_1 matrix shape: {B1.shape}")
rk_B1 = np.linalg.matrix_rank(B1)
print(f"rank(∂_1) = {rk_B1}")
print(f"ker(∂_1) = {len(edges) - rk_B1}")

# Build boundary ∂_2: Ω_2 → Ω_1
# First, ∂_2 as a map from A_2 to A_1
B2_full = np.zeros((len(edges), len(paths2)))
edge_idx = {e: i for i, e in enumerate(edges)}
for j, path in enumerate(paths2):
    for i in range(3):
        face = path[:i] + path[i+1:]
        if face in edge_idx:
            B2_full[edge_idx[face], j] += (-1)**i

# Restrict to Ω_2
omega2_basis = Vh[rank_J:].T if omega2_dim > 0 else np.zeros((len(paths2), 0))
B2_omega = B2_full @ omega2_basis

print(f"\n∂_2 restricted to Ω_2:")
print(f"  B2_omega shape: {B2_omega.shape}")
rk_B2 = np.linalg.matrix_rank(B2_omega, tol=1e-8)
print(f"  rank(∂_2|Ω_2) = {rk_B2}")

# β_1 = dim(ker ∂_1) - dim(im ∂_2|Ω_2)
beta_1 = (len(edges) - rk_B1) - rk_B2
print(f"\nβ_1 = ker(∂_1) - im(∂_2|Ω_2) = {len(edges) - rk_B1} - {rk_B2} = {beta_1}")

# Also check β_0
# β_0 = dim(Ω_0) - 0 - rk(∂_1) = 5 - 4 = 1
print(f"β_0 = {5} - 0 - {rk_B1} = {5 - rk_B1}")

# Summary
print(f"\nSummary for P_5 = C5_{{1,4}}:")
print(f"  A_m: [5, {len(edges)}, {len(paths2)}]")
print(f"  Ω_m: [5, {len(edges)}, {omega2_dim}]")
print(f"  β: [{5 - rk_B1}, {beta_1}, ...]")

# Compare with C5_{1,2}
print(f"\n{'='*50}")
print("Comparison: C5_{1,2} (regular/interval)")
A2 = circulant_tournament(5, {1,2})
print(f"C5_{{1,2}} adjacency:")
print(A2)
print(f"Out-degrees: {[sum(A2[i]) for i in range(5)]}")

# Are they isomorphic?
from itertools import permutations
def canon(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

print(f"\nIsomorphic? {canon(A) == canon(A2)}")

# P_5 edge pattern: i→i+1 and i→i+4=i-1
# This means BOTH i→i+1 AND i→i-1 (bidirectional nearest neighbor)
# Wait, that can't be right for a tournament...
# {1,4} mod 5: 4 ≡ -1. So i→i+1 and i→i-1.
# This creates edges BOTH ways between adjacent vertices? No!
# Tournament: exactly one of (i,j) or (j,i) has an edge.
# For |i-j|=1: if j=i+1, j-i=1∈S, so i→j. Also if j=i-1, then i-j=1, j-i=-1≡4∈S, so j→i... wait
# j-i = -1 ≡ 4 mod 5, so i→j when j-i ∈ {1,4}
# For j=i+1: j-i=1 ∈ S, so i→i+1 ✓
# For j=i-1: j-i=-1≡4 ∈ S, so i→i-1 ✓
# But then 0→1 AND 0→4, 1→2 AND 1→0...
# So 0→1 and 1→0 both? That can't be a tournament!

# Wait, tournament means A[i][j] + A[j][i] = 1 for i≠j
# Check:
for i in range(5):
    for j in range(i+1, 5):
        if A[i][j] + A[j][i] != 1:
            print(f"  NOT tournament at ({i},{j}): A[{i}][{j}]={A[i][j]}, A[{j}][{i}]={A[j][i]}")

print("\nDONE.")
