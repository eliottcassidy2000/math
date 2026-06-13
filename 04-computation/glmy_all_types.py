#!/usr/bin/env python3
"""
glmy_all_types.py — opus-2026-03-13-S71

Compute CORRECT GLMY path homology Betti numbers for all tournament
isomorphism types at n=5 and n=6. Compare with TRH (interior-only on
regular paths).

Also: characterize which tournaments have TRH d²=0.
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict

def enumerate_allowed_paths(A, n, p):
    """Directed p-paths: (v0,...,vp) with edges v_i→v_{i+1}."""
    if p == 0: return [(v,) for v in range(n)]
    paths = []
    def dfs(path, depth):
        if depth == p:
            paths.append(tuple(path))
            return
        last = path[-1]
        visited = set(path)
        for v in range(n):
            if v not in visited and A[last][v]:
                path.append(v)
                dfs(path, depth+1)
                path.pop()
    for s in range(n):
        dfs([s], 0)
    return paths

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0: return np.zeros((0,0))
    if p == 0: return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed = {}
    na_count = 0
    for path in allowed_p:
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face not in allowed_pm1_set and face not in non_allowed:
                non_allowed[face] = na_count
                na_count += 1
    if na_count == 0: return np.eye(dim_Ap)
    J = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for i in range(p+1):
            face = path[:i] + path[i+1:]
            if face in non_allowed:
                J[non_allowed[face], j] += (-1)**i
    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    if rank == dim_Ap: return np.zeros((dim_Ap, 0))
    return Vh[rank:].T

def glmy_betti(A):
    n = A.shape[0]
    max_dim = n - 1
    allowed = {}
    omega = {}
    for p in range(max_dim+2):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        omega[p] = compute_omega_basis(A, n, p, allowed[p],
                                        allowed.get(p-1, []))
    ranks = {}
    for p in range(1, max_dim+1):
        dim_om = omega[p].shape[1]
        dim_om_prev = omega[p-1].shape[1]
        if dim_om == 0 or dim_om_prev == 0:
            ranks[p] = 0
            continue
        idx = {path: i for i, path in enumerate(allowed[p-1])}
        B = np.zeros((len(allowed[p-1]), len(allowed[p])))
        for j, path in enumerate(allowed[p]):
            for i in range(p+1):
                face = path[:i] + path[i+1:]
                if face in idx:
                    B[idx[face], j] += (-1)**i
        B_omega = omega[p-1].T @ B @ omega[p]
        ranks[p] = np.linalg.matrix_rank(B_omega, tol=1e-8)
    betti = []
    for p in range(max_dim+1):
        dim_p = omega[p].shape[1]
        rk_p = ranks.get(p, 0)
        rk_p1 = ranks.get(p+1, 0)
        betti.append(dim_p - rk_p - rk_p1)
    omega_dims = [omega[p].shape[1] for p in range(max_dim+1)]
    return betti, omega_dims

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    return best

def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]+A[j][k]+A[k][i]==3 or A[j][i]+A[i][k]+A[k][j]==3:
                    count += 1
    return count

# ============================================================
print("="*70)
print("GLMY PATH HOMOLOGY FOR ALL n=5 TOURNAMENT TYPES")
print("="*70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
types = {}
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    cf = canon_form(A)
    if cf not in types:
        types[cf] = A

for idx, (cf, A) in enumerate(sorted(types.items())):
    betti, omega_dims = glmy_betti(A)
    t3 = count_3cycles(A, n)
    scores = tuple(sorted(int(sum(A[v][w] for w in range(n) if w != v))
                          for v in range(n)))
    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    print(f"  Type {idx+1:2d}: score={scores}, t3={t3}, "
          f"Ω={omega_dims}, β={betti}, chi={chi}")

# β_2 check
all_b2_zero = all(glmy_betti(A)[0][2] == 0 for A in types.values())
print(f"\n  GLMY β_2 = 0 for ALL? {all_b2_zero}")

# ============================================================
print(f"\n{'='*70}")
print("GLMY PATH HOMOLOGY FOR ALL n=6 TOURNAMENT TYPES")
print("="*70)

n = 6
pairs6 = [(i,j) for i in range(n) for j in range(i+1,n)]
types6 = {}
for bits in range(2**len(pairs6)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs6):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    cf = canon_form(A)
    if cf not in types6:
        types6[cf] = A

print(f"  {len(types6)} isomorphism types")

glmy_dist = Counter()
glmy_b2 = Counter()
for cf, A in types6.items():
    betti, _ = glmy_betti(A)
    glmy_dist[tuple(betti)] += 1
    glmy_b2[betti[2]] += 1

print(f"\n  β_2 distribution: {dict(sorted(glmy_b2.items()))}")
all_b2_zero_6 = all(k == 0 for k in glmy_b2.keys())
print(f"  GLMY β_2 = 0 for ALL? {all_b2_zero_6}")

print(f"\n  Distinct GLMY Betti vectors:")
for bv, count in sorted(glmy_dist.items()):
    chi = sum((-1)**m * b for m, b in enumerate(bv))
    print(f"    β={list(bv)}, chi={chi}, count={count}")

print("\nDONE.")
