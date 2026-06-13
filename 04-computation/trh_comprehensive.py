#!/usr/bin/env python3
"""
trh_comprehensive.py — opus-2026-03-13-S71

Comprehensive analysis of the Tournament Regular Homology (TRH) chain complex.
TRH uses:
  - Allowed paths = regular paths (v_i→v_{i+1} AND v_{i-1}→v_{i+1})
  - Boundary = interior-only deletion (indices 1 to m-1)
  - d² = 0 verified

Key questions:
1. Is TRH a known construction? (regular path complex + interior boundary)
2. What topological properties does it capture?
3. Why β_0 = n always? (Because ∂_1 = 0 in TRH)
4. The eigenspace Betti uniformity and divisibility — why?
5. Comparison with GLMY at n=5, n=6

Also: explore whether TRH = "magnitude homology" or something related.
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict

def get_all_tournaments(n):
    """Generate all tournaments on n vertices."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    num_pairs = len(pairs)
    tournaments = []
    for bits in range(2**num_pairs):
        A = np.zeros((n,n), dtype=int)
        for k, (i,j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        tournaments.append(A)
    return tournaments

def canon_form(A):
    n = A.shape[0]
    best = None
    for perm in permutations(range(n)):
        enc = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(i+1,n))
        if best is None or enc < best:
            best = enc
    return best

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

def trh_betti(A):
    """TRH Betti numbers: regular paths + interior boundary."""
    n = A.shape[0]
    all_paths = {}
    for m in range(n):
        all_paths[m] = get_regular_paths(A, m)

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

def count_3cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] + A[j][k] + A[k][i] == 3 or
                    A[j][i] + A[i][k] + A[k][j] == 3):
                    count += 1
    return count

# ============================================================
# TRH for all tournament types at n=5
# ============================================================
print("="*70)
print("TRH BETTI FOR ALL n=5 TOURNAMENT ISOMORPHISM TYPES")
print("="*70)

n = 5
types = {}
for A in get_all_tournaments(n):
    cf = canon_form(A)
    if cf not in types:
        types[cf] = A

for idx, (cf, A) in enumerate(sorted(types.items())):
    betti, omegas = trh_betti(A)
    t3 = count_3cycles(A, n)
    scores = tuple(sorted(sum(A[v][w] for w in range(n) if w != v) for v in range(n)))
    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    print(f"  Type {idx+1:2d}: score={scores}, t3={t3}, "
          f"Omega={omegas}, β={betti}, chi={chi}")

# ============================================================
# Verify: TRH β_2 = 0 for all n=5
# ============================================================
print(f"\n  TRH β_2 = 0 for ALL n=5 types? "
      f"{all(trh_betti(A)[0][2] == 0 for A in types.values())}")

# ============================================================
# TRH for all tournament types at n=6
# ============================================================
print(f"\n{'='*70}")
print("TRH BETTI FOR ALL n=6 TOURNAMENT ISOMORPHISM TYPES")
print("="*70)

n = 6
types6 = {}
count = 0
for A in get_all_tournaments(n):
    cf = canon_form(A)
    if cf not in types6:
        types6[cf] = A
        count += 1

print(f"  Found {count} isomorphism types at n=6")
print()

# Collect statistics
beta_dist = Counter()
beta2_vals = Counter()
chi_vals = Counter()

for idx, (cf, A) in enumerate(sorted(types6.items())):
    betti, omegas = trh_betti(A)
    chi = sum((-1)**m * betti[m] for m in range(len(betti)))
    beta_dist[tuple(betti)] += 1
    beta2_vals[betti[2]] += 1
    chi_vals[chi] += 1

print(f"  β_2 distribution: {dict(sorted(beta2_vals.items()))}")
print(f"  chi distribution: {dict(sorted(chi_vals.items()))}")
print(f"  β_2 = 0 for ALL? {all(k == 0 for k in beta2_vals.keys())}")

# Show all distinct Betti vectors
print(f"\n  Distinct TRH Betti vectors (n=6):")
for betti_vec, count in sorted(beta_dist.items()):
    chi = sum((-1)**m * b for m, b in enumerate(betti_vec))
    print(f"    β={list(betti_vec)}, chi={chi}, count={count}")

# ============================================================
# Key structural analysis: why β_0 = n?
# ============================================================
print(f"\n{'='*70}")
print("WHY β_0 = n IN TRH?")
print("="*70)

print("""
  In TRH, ∂_1 = 0 because:
  - A regular 1-path is (v_0, v_1) with v_0 → v_1.
  - Interior vertices = indices 1 to m-1 = 1 to 0 = empty.
  - So ∂_1(v_0, v_1) = 0 (empty sum).
  Therefore ker(∂_0) = Ω_0 = all vertices, and im(∂_1) = {0}.
  β_0 = dim(ker ∂_0) - dim(im ∂_1) = n - 0 = n.

  In GLMY, ∂_1(v_0, v_1) = (v_1) - (v_0), giving rank n-1, so β_0 = 1.

  The TRH β_0 = n is a "multiplicity" reflecting that each vertex
  is an isolated 0-cycle. This is analogous to homology with fixed endpoints.
""")

# ============================================================
# Check: is TRH chi always related to n?
# ============================================================
print("="*70)
print("TRH CHI ANALYSIS")
print("="*70)

for n_test in [3, 4, 5]:
    types_test = {}
    for A in get_all_tournaments(n_test):
        cf = canon_form(A)
        if cf not in types_test:
            types_test[cf] = A
    chi_dist = Counter()
    for A in types_test.values():
        betti, _ = trh_betti(A)
        chi = sum((-1)**m * betti[m] for m in range(len(betti)))
        chi_dist[chi] += 1
    print(f"  n={n_test}: chi distribution = {dict(sorted(chi_dist.items()))}")

# Same for n=6 (already computed)
print(f"  n=6: chi distribution = {dict(sorted(chi_vals.items()))}")

print("\nDONE.")
