#!/usr/bin/env python3
"""
betti_divisibility.py — opus-2026-03-13-S70

Key observation: Betti numbers for circulant tournaments on Z_p appear
to be divisible by p!

Interval n=7: β = (7, 7, 0, 14, 14, 7, 0) — all divisible by 7
Interval n=9: β = (9, 9, 0, 27, 27, 18, 27, 27, 27) — all divisible by 9
Paley p=7: β = (7, 0, 0, 21, 21, 21, 21) — all divisible by 7

This is likely because the Z_p symmetry group acts on the chain complex,
and the Betti numbers decompose over p eigenspaces.

Verify and explore this pattern for multiple circulant tournaments.
Also: check if β_m/p are related to per-eigenspace Betti numbers.
"""

import numpy as np
from itertools import permutations
from collections import Counter

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

def get_regular_paths(A, m):
    """Get all regular m-paths as tuples."""
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
    """Build boundary matrix ∂_m: Ω_m → Ω_{m-1}.
    ∂_m(v_0,...,v_m) = Σ_{i=1}^{m-1} (-1)^i (v_0,...,v̂_i,...,v_m)
    Only delete interior vertices (indices 1 to m-1).
    """
    if not paths_m or not paths_m_minus_1:
        return np.zeros((len(paths_m_minus_1), len(paths_m)))

    path_to_idx = {p: i for i, p in enumerate(paths_m_minus_1)}
    m = len(paths_m[0]) - 1  # path length

    matrix = np.zeros((len(paths_m_minus_1), len(paths_m)), dtype=int)

    for j, path in enumerate(paths_m):
        for i in range(1, m):  # interior vertices only
            face = path[:i] + path[i+1:]
            sign = (-1)**i
            if face in path_to_idx:
                matrix[path_to_idx[face], j] += sign

    return matrix

def compute_betti(A, max_dim=None):
    """Compute all Betti numbers."""
    n = A.shape[0]
    if max_dim is None:
        max_dim = n - 1

    # Get regular paths for each dimension
    all_paths = {}
    for m in range(max_dim + 1):
        all_paths[m] = get_regular_paths(A, m)

    # Compute boundary matrices
    ranks = {}
    for m in range(1, max_dim + 1):
        if all_paths[m] and all_paths[m-1]:
            B = boundary_matrix(all_paths[m], all_paths[m-1])
            ranks[m] = np.linalg.matrix_rank(B)
        else:
            ranks[m] = 0

    # Betti numbers
    betti = []
    for m in range(max_dim + 1):
        omega_m = len(all_paths[m])
        rank_dm = ranks.get(m, 0)  # rank of ∂_m
        rank_dm_plus_1 = ranks.get(m+1, 0)  # rank of ∂_{m+1}
        ker_dm = omega_m - rank_dm
        im_dm_plus_1 = rank_dm_plus_1
        beta_m = ker_dm - im_dm_plus_1
        betti.append(beta_m)

    omegas = [len(all_paths[m]) for m in range(max_dim + 1)]
    return betti, omegas

# ============================================================
# Check divisibility by n for various circulant tournaments
# ============================================================
print("="*70)
print("BETTI DIVISIBILITY BY n FOR CIRCULANT TOURNAMENTS")
print("="*70)

def legendre(a, p):
    if a % p == 0: return 0
    return pow(a, (p-1)//2, p)

tests = []

# Interval tournaments
for n in [5, 7, 9, 11]:
    m = (n-1)//2
    S = set(range(1, m+1))
    tests.append((f"Interval n={n}", n, S))

# Paley tournaments (p ≡ 3 mod 4)
for p in [3, 7, 11]:
    if p % 4 != 3: continue
    QR = {a % p for a in range(1, p) if legendre(a, p) == 1}
    tests.append((f"Paley p={p}", p, QR))

# Other circulants at n=7
for S_list in [[1,2,4], [1,3,5]]:
    S = set(S_list)
    if S | {7-s for s in S} == set(range(1,7)):
        tests.append((f"C7 S={S_list}", 7, S))

for name, n, S in tests:
    print(f"\n  {name}:")
    A = circulant_tournament(n, S)
    betti, omegas = compute_betti(A)

    print(f"    Omega: {omegas}")
    print(f"    Betti: {betti}")

    # Check divisibility
    all_div = all(b % n == 0 for b in betti)
    betti_over_n = [b // n for b in betti]
    print(f"    β/n:   {betti_over_n}")
    print(f"    All β divisible by n={n}: {all_div}")

    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    print(f"    chi = {chi}, chi/n = {chi/n}")

# ============================================================
# For n=5: check all 12 types — which have β divisible by 5?
# ============================================================
print(f"\n{'='*70}")
print("n=5: DIVISIBILITY BY 5 FOR ALL TYPES")
print("="*70)

def adj_matrix_5(bits):
    A = np.zeros((5,5), dtype=int)
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = []
        for i in range(n):
            for j in range(i+1, n):
                enc.append(A[perm[i]][perm[j]])
        enc = tuple(enc)
        if best is None or enc < best:
            best = enc
    return best

types = {}
for bits in range(2**10):
    A = adj_matrix_5(bits)
    cf = canon_form(A)
    if cf not in types:
        types[cf] = bits

for cf, bits in sorted(types.items()):
    A = adj_matrix_5(bits)
    betti, omegas = compute_betti(A)
    scores = tuple(sorted(sum(A[v][w] for w in range(5) if w != v) for v in range(5)))
    all_div = all(b % 5 == 0 for b in betti)
    betti_5 = [b // 5 if b % 5 == 0 else f"{b}" for b in betti]

    # Check if it's a circulant tournament
    is_circulant = False
    for S_candidate in [[1,2], [2,4], [1,4], [2,3], [1,3], [3,4]]:
        S = set(S_candidate)
        if S | {5-s for s in S} == set(range(1,5)):
            A_circ = circulant_tournament(5, S)
            if canon_form(A_circ) == cf:
                is_circulant = True
                break

    circ_tag = " [circulant]" if is_circulant else ""
    print(f"  {str(scores):20s} β={betti} div5={all_div} β/5={betti_5}{circ_tag}")

# ============================================================
# Per-eigenspace Betti for Paley p=7
# ============================================================
print(f"\n{'='*70}")
print("PALEY p=7: PER-EIGENSPACE BETTI")
print("="*70)

p = 7
QR = {1, 2, 4}
A = circulant_tournament(p, QR)
betti, omegas = compute_betti(A)
print(f"  Total Betti: {betti}")
print(f"  Total Omega: {omegas}")

# Since all Q_k are equal, all eigenspaces should contribute equally
# β_m / p should give the per-eigenspace Betti
print(f"\n  Per-eigenspace (β/p):")
per_eig = [b // p for b in betti]
print(f"    {per_eig}")
print(f"  This is the Betti of the 'eigenspace chain complex'")

# Note: β_0/p = 1, β_1/p = 0, β_2/p = 0, β_3/p = 3, β_4/p = 3, β_5/p = 3, β_6/p = 3
# What does this mean?

# For Interval n=7:
S_int = {1, 2, 3}
A_int = circulant_tournament(p, S_int)
betti_int, _ = compute_betti(A_int)
per_eig_int = [b // p for b in betti_int]
print(f"\n  Interval n=7: β/p = {per_eig_int}")
print(f"  Interval total β: {betti_int}")

print("\nDONE.")
