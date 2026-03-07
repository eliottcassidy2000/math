#!/usr/bin/env python3
"""
Spectral analysis of the blue line skeleton at n=5.

The skeleton is a bipartite graph. Questions:
1. What are its eigenvalues?
2. Do they relate to tournament invariants?
3. Is the adjacency matrix related to the transfer matrix?

Also: explore the "magnetic" interpretation.
The skeleton is an antiferromagnetic Ising model on tournaments:
- Each SC class has a "spin" (t3 parity)
- GS flip edges always connect opposite spins
- The skeleton IS the interaction graph of this magnetic system

kind-pasteur-2026-03-06-S25h
"""
from itertools import permutations, combinations
from collections import defaultdict

def tournament_from_tiling(n, tiling_bits):
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (tiling_bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def num_tiling_bits(n):
    return n*(n-1)//2 - (n-1)

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def count_3cycles(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_ham_paths(A, n):
    """Count Hamiltonian paths using DP."""
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def dp(mask, last):
        if mask == (1 << n) - 1:
            return 1
        total = 0
        for v in range(n):
            if not (mask & (1 << v)) and A[last][v]:
                total += dp(mask | (1 << v), v)
        return total
    return sum(dp(1 << v, v) for v in range(n))

def tiling_transpose_pairs(n):
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}
    pairs = []
    fixed = []
    seen = set()
    for idx, (i, j) in enumerate(edges):
        if idx in seen:
            continue
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if (ti, tj) == (i, j):
            fixed.append(idx)
            seen.add(idx)
        elif (ti, tj) in edge_to_idx:
            tidx = edge_to_idx[(ti, tj)]
            pairs.append((idx, tidx))
            seen.add(idx)
            seen.add(tidx)
    return pairs, fixed

def gen_gs_tilings(n, pairs, fixed):
    gs_dof = len(pairs) + len(fixed)
    result = []
    for free_val in range(2**gs_dof):
        bits = 0
        for k, (idx1, idx2) in enumerate(pairs):
            if (free_val >> k) & 1:
                bits |= (1 << idx1) | (1 << idx2)
        for k, fidx in enumerate(fixed):
            if (free_val >> (len(pairs) + k)) & 1:
                bits |= (1 << fidx)
        result.append(bits)
    return result

def flip_tiling(tiling_bits, m):
    return tiling_bits ^ ((1 << m) - 1)

n = 5
m = num_tiling_bits(n)
print(f"SPECTRAL ANALYSIS at n={n}")
print(f"{'='*60}")

pairs, fixed = tiling_transpose_pairs(n)
gs_tilings = gen_gs_tilings(n, pairs, fixed)

# Build class database
canon_db = {}
class_list = []
bits_to_class = {}

for bits in range(2**m):
    A = tournament_from_tiling(n, bits)
    c = canonical(A, n)
    if c not in canon_db:
        canon_db[c] = len(class_list)
        class_list.append({
            'rep': A, 'sc': is_SC(A, n),
            'scores': score_seq(A, n),
            't3': count_3cycles(A, n),
            'H': count_ham_paths(A, n),
            'gs_tilings': set(), 'tilings': set()
        })
    idx = canon_db[c]
    class_list[idx]['tilings'].add(bits)
    bits_to_class[bits] = idx

for bits in gs_tilings:
    idx = bits_to_class[bits]
    class_list[idx]['gs_tilings'].add(bits)

sc_indices = [i for i, c in enumerate(class_list) if c['sc']]
print(f"  {len(class_list)} classes, {len(sc_indices)} SC")

# Build adjacency matrix for the SC skeleton
sc_to_idx = {v: i for i, v in enumerate(sc_indices)}
nsc = len(sc_indices)
adj_matrix = [[0]*nsc for _ in range(nsc)]

gs_edges = defaultdict(int)
for bits in gs_tilings:
    c_from = bits_to_class[bits]
    flipped = flip_tiling(bits, m)
    c_to = bits_to_class[flipped]
    if c_from != c_to and c_from in sc_to_idx and c_to in sc_to_idx:
        edge = (min(c_from, c_to), max(c_from, c_to))
        gs_edges[edge] += 1

for (c_from, c_to), weight in gs_edges.items():
    i, j = sc_to_idx[c_from], sc_to_idx[c_to]
    adj_matrix[i][j] = weight
    adj_matrix[j][i] = weight

# Print adjacency matrix
print(f"\n  Adjacency matrix (weighted):")
header = "   " + " ".join(f"{sc_indices[j]:3d}" for j in range(nsc))
print(header)
for i in range(nsc):
    row = f"{sc_indices[i]:3d}" + " ".join(f"{adj_matrix[i][j]:3d}" for j in range(nsc))
    print(row)

# Eigenvalues (using numpy if available, else manual for small matrix)
try:
    import numpy as np
    A_np = np.array(adj_matrix, dtype=float)
    eigenvalues = sorted(np.linalg.eigvalsh(A_np))[::-1]
    print(f"\n  Eigenvalues: {[f'{e:.4f}' for e in eigenvalues]}")
    print(f"  Spectral radius: {max(abs(e) for e in eigenvalues):.4f}")

    # For bipartite graph, eigenvalues are symmetric around 0
    print(f"  Symmetric around 0: {all(abs(eigenvalues[i] + eigenvalues[-(i+1)]) < 0.001 for i in range(nsc))}")

    # Eigenvectors
    eigenvalues_full, eigenvectors = np.linalg.eigh(A_np)
    idx_sorted = np.argsort(-eigenvalues_full)

    print(f"\n  Principal eigenvector (Perron):")
    v1 = eigenvectors[:, idx_sorted[0]]
    for i, sc_class in enumerate(sc_indices):
        c = class_list[sc_class]
        print(f"    Class {sc_class} (t3={c['t3']}, H={c['H']}, scores={c['scores']}): v1={v1[i]:.4f}")

except ImportError:
    print("\n  numpy not available, computing eigenvalues manually for 8x8...")
    # For 8x8 bipartite, eigenvalues come in +/- pairs
    # Just print the adjacency matrix structure
    pass

# Analyze: does H(T) correlate with graph-theoretic position?
print(f"\n  CLASS PROPERTIES:")
for i, sc_class in enumerate(sc_indices):
    c = class_list[sc_class]
    degree = sum(adj_matrix[i])
    t3_side = "A" if c['t3'] % 2 == 1 else "B"
    print(f"    Class {sc_class}: t3={c['t3']}, H={c['H']}, degree={degree}, "
          f"side={t3_side}, #GS={len(c['gs_tilings'])}, #total={len(c['tilings'])}")

# GS fraction: what fraction of each class's tilings are GS?
print(f"\n  GS FRACTION per class:")
for i, sc_class in enumerate(sc_indices):
    c = class_list[sc_class]
    gs_frac = len(c['gs_tilings']) / len(c['tilings'])
    print(f"    Class {sc_class}: {len(c['gs_tilings'])}/{len(c['tilings'])} = {gs_frac:.4f}")

# H(T) on each side
side_A = [i for i in sc_indices if class_list[i]['t3'] % 2 == 1]
side_B = [i for i in sc_indices if class_list[i]['t3'] % 2 == 0]
H_A = [class_list[i]['H'] for i in side_A]
H_B = [class_list[i]['H'] for i in side_B]
print(f"\n  H(T) on side A (odd t3): {sorted(H_A)}")
print(f"  H(T) on side B (even t3): {sorted(H_B)}")
print(f"  Mean H: side A = {sum(H_A)/len(H_A):.1f}, side B = {sum(H_B)/len(H_B):.1f}")
