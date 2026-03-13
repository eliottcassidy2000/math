#!/usr/bin/env python3
"""
glmy_circulant_compare.py — opus-2026-03-13-S71

Compare GLMY (correct) vs TRH vs circulant_homology Betti numbers
for key circulant tournaments. Also compute GLMY for n=7 circulants.

The three chain complexes:
1. GLMY: directed paths + full boundary + Ω subspace
2. TRH: regular paths + interior-only boundary + A_m directly
3. circulant_homology: regular paths + full boundary + (some Ω?)
"""

import numpy as np
from itertools import permutations

# ============================================================
# GLMY correct implementation (directed paths + Ω subspace)
# ============================================================

def enumerate_directed_paths(A, n, p):
    """Directed p-paths: (v0,...,vp) with edges v_i→v_{i+1}, all distinct."""
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
    if p <= 1: return np.eye(dim_Ap)
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
        allowed[p] = enumerate_directed_paths(A, n, p)
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
    a_dims = [len(allowed[p]) for p in range(max_dim+1)]
    return betti, omega_dims, a_dims

# ============================================================
# TRH implementation (regular paths + interior boundary)
# ============================================================

def enumerate_regular_paths(A, n, m):
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

def trh_betti(A):
    n = A.shape[0]
    all_paths = {}
    for m in range(n):
        all_paths[m] = enumerate_regular_paths(A, n, m)
    ranks = {}
    for m in range(1, n):
        if not all_paths[m] or not all_paths[m-1]:
            ranks[m] = 0
            continue
        path_to_idx = {p: i for i, p in enumerate(all_paths[m-1])}
        path_m = len(all_paths[m][0]) - 1
        B = np.zeros((len(all_paths[m-1]), len(all_paths[m])), dtype=int)
        for j, path in enumerate(all_paths[m]):
            for i in range(1, path_m):
                face = path[:i] + path[i+1:]
                sign = (-1)**i
                if face in path_to_idx:
                    B[path_to_idx[face], j] += sign
        ranks[m] = np.linalg.matrix_rank(B)
    betti = []
    for m in range(n):
        omega_m = len(all_paths[m])
        rank_dm = ranks.get(m, 0)
        rank_dm_plus_1 = ranks.get(m+1, 0)
        betti.append(omega_m - rank_dm - rank_dm_plus_1)
    omegas = [len(all_paths[m]) for m in range(n)]
    return betti, omegas

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

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
print("GLMY vs TRH FOR KEY CIRCULANT TOURNAMENTS")
print("="*70)

test_cases = [
    ("C3 cyclic", 3, {1}),
    ("C5 regular", 5, {1,2}),
    ("C5 interval", 5, {1,2}),  # same as regular at n=5
    ("C7 Paley", 7, {1,2,4}),
    ("C7 interval {1,2,3}", 7, {1,2,3}),
    ("C7 {1,3,5}", 7, {1,3,5}),
]

for name, n, S in test_cases:
    A = circulant_tournament(n, S)
    t3 = count_3cycles(A, n)

    # GLMY
    bg, og, ag = glmy_betti(A)
    chi_g = sum((-1)**m * bg[m] for m in range(len(bg)))

    # TRH
    bt, ot = trh_betti(A)
    chi_t = sum((-1)**m * bt[m] for m in range(len(bt)))

    print(f"\n  {name} (n={n}, t3={t3}):")
    print(f"    GLMY directed: A={ag}")
    print(f"    GLMY Ω:        {og}")
    print(f"    GLMY β:        {bg}, chi={chi_g}")
    print(f"    TRH  regular:  {ot}")
    print(f"    TRH  β:        {bt}, chi={chi_t}")

    # Check divisibility by n
    if all(b % n == 0 for b in bg):
        print(f"    GLMY β/n: {[b//n for b in bg]}")
    if all(b % n == 0 for b in bt):
        print(f"    TRH  β/n: {[b//n for b in bt]}")

# ============================================================
print(f"\n{'='*70}")
print("DIRECTED vs REGULAR PATH COUNTS AT n=7")
print("="*70)

for name, n, S in [("C7 Paley", 7, {1,2,4}), ("C7 interval", 7, {1,2,3})]:
    A = circulant_tournament(n, S)
    print(f"\n  {name}:")
    for m in range(n):
        dp = enumerate_directed_paths(A, n, m)
        rp = enumerate_regular_paths(A, n, m)
        print(f"    m={m}: directed={len(dp):6d}, regular={len(rp):5d}, ratio={len(dp)/max(len(rp),1):.2f}")

# ============================================================
print(f"\n{'='*70}")
print("GLMY Ω STRUCTURE FOR n=7 CIRCULANTS — DETAILED")
print("="*70)

for name, n, S in [("C7 Paley", 7, {1,2,4}), ("C7 interval", 7, {1,2,3}),
                    ("C7 {1,3,5}", 7, {1,3,5})]:
    A = circulant_tournament(n, S)
    bg, og, ag = glmy_betti(A)
    print(f"\n  {name}:")
    print(f"    A_m:  {ag}")
    print(f"    Ω_m:  {og}")
    print(f"    β_m:  {bg}")
    gap = [ag[m] - og[m] for m in range(len(ag))]
    print(f"    Gap:  {gap}")

    # Check β divisibility
    for d in [n, n-1, n+1]:
        if all(b % d == 0 for b in bg if b > 0):
            print(f"    ALL nonzero β divisible by {d}")

    # Ω divisibility
    for d in [n]:
        if all(o % d == 0 for o in og):
            print(f"    ALL Ω divisible by {d}: {[o//d for o in og]}")

# ============================================================
print(f"\n{'='*70}")
print("GLMY FOR TRANSITIVE TOURNAMENTS T_n")
print("="*70)

for n in [3, 4, 5, 6, 7]:
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    bg, og, ag = glmy_betti(A)
    from math import comb
    expected = [comb(n, m+1) for m in range(n)]
    print(f"  T_{n}: β={bg}, Ω={og}, A={ag}")
    if og == expected:
        print(f"        Ω = C(n,m+1) ✓")

# ============================================================
print(f"\n{'='*70}")
print("GLMY vs TRH — WHERE DO THEY AGREE/DIFFER?")
print("="*70)

# Check all n=5 types
from itertools import permutations as perms

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
types = {}
for bits in range(2**len(pairs)):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    best = None
    for perm in perms(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best: best = enc
    if best not in types:
        types[best] = A

agree_count = 0
differ_count = 0
for cf, A in sorted(types.items()):
    bg, og, ag = glmy_betti(A)
    bt, ot = trh_betti(A)

    # Compare (note: TRH β_0 = n, GLMY β_0 = 1, so skip β_0)
    glmy_tail = bg[1:]
    trh_tail = bt[1:]

    if glmy_tail == trh_tail:
        agree_count += 1
    else:
        differ_count += 1
        scores = tuple(sorted(sum(A[v][w] for w in range(n) if w != v) for v in range(n)))
        print(f"  DIFFER: score={scores}")
        print(f"    GLMY β={bg}")
        print(f"    TRH  β={bt}")

print(f"\n  n=5: {agree_count} agree (β_1..β_{n-1}), {differ_count} differ")

print("\nDONE.")
