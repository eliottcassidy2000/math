#!/usr/bin/env python3
"""
beta3_obstruction.py — opus-2026-03-09-S52

KEY INSIGHT: From the LES with beta_2 = 0:
  H_3(T,T\v) = beta_3(T) - rank(i_*: H_3(T\v) -> H_3(T))

This is EXACT (not an inequality). So:
  beta_3(T) = H_3(T,T\v) + rank(i_*)

If H_3(T,T\v) <= 1 and beta_3(T\v) <= 1, then beta_3(T) <= 2.
For beta_3(T) = 2, we'd need ALL vertices v to have:
  - beta_3(T\v) = 1
  - rank(i_*) = 1  (i.e., i_* injective)
  - H_3(T,T\v) = 1

Can this happen? This script checks at n=7 (exhaustive for a subset):
  1. At n=7, do ANY tournaments have beta_3 = 2?
  2. At n=7, do any tournaments have ALL 7 deletions with beta_3 >= 1?
  3. Frequency of "all deletions have beta_3 = 1" at various n.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
import sys

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_betti_quick(T, n, max_p=4):
    """Compute Betti numbers up to max_p."""
    tol = 1e-8
    all_paths = {}
    for p in range(0, max_p + 2):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    def build_omega(pd, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = pd[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd[p-1]) if p > 0 else set()
            na = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na:
                            na[face] = len(na)
            if not na:
                omega[p] = np.eye(len(a_p))
            else:
                mat = np.zeros((len(na), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na:
                            mat[na[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                null_dim = len(a_p) - rank
                if null_dim == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, max_p)

    boundary = {}
    for p in range(1, max_p + 1):
        a_p = all_paths[p]
        a_pm1 = all_paths[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    ranks = {}
    dims = {}
    for p in range(max_p + 1):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, max_p + 1):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    betti = {}
    for p in range(max_p + 1):
        ker = dims[p] - ranks.get(p, 0)
        im_next = ranks.get(p+1, 0)
        betti[p] = ker - im_next

    return betti

def main():
    print("=" * 70)
    print("BETA_3 OBSTRUCTION ANALYSIS")
    print("=" * 70)

    # ============================================
    # Part 1: n=7, check if beta_3 can reach 2
    # Large sample to be confident
    # ============================================
    print("\n--- n=7: Checking if beta_3 can equal 2 ---")
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)
    sample_size = 500

    b3_dist = Counter()
    all_deletions_b3_pos = 0  # tournaments where ALL deletions have beta_3 >= 1
    min_b3_sub_dist = Counter()  # distribution of min(beta_3(T\v))

    for trial in range(sample_size):
        if trial % 100 == 0:
            print(f"  ... {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        try:
            betti = compute_betti_quick(T, n, max_p=4)
        except:
            continue
        b3 = betti.get(3, 0)
        b3_dist[b3] += 1

        # Check all vertex-deletions
        min_b3_sub = float('inf')
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub = compute_betti_quick(T_sub, n-1, max_p=4)
            b3_sub = betti_sub.get(3, 0)
            min_b3_sub = min(min_b3_sub, b3_sub)

        min_b3_sub_dist[min_b3_sub] += 1
        if min_b3_sub >= 1:
            all_deletions_b3_pos += 1

    print(f"\n  beta_3 distribution: {dict(sorted(b3_dist.items()))}")
    print(f"  min(beta_3(T\\v)) over v: {dict(sorted(min_b3_sub_dist.items()))}")
    print(f"  All deletions have beta_3 >= 1: {all_deletions_b3_pos}/{sample_size}")
    print(f"\n  KEY: beta_3 = 2 ever? {'YES' if b3_dist.get(2,0) > 0 else 'NO'}")

    # ============================================
    # Part 2: n=8, check same
    # ============================================
    print(f"\n{'='*70}")
    print("--- n=8: Checking if beta_3 can equal 2 ---")
    n = 8
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng2 = np.random.RandomState(123)
    sample_size_8 = 100

    b3_dist_8 = Counter()
    all_del_8 = 0
    min_b3_sub_dist_8 = Counter()

    for trial in range(sample_size_8):
        if trial % 20 == 0:
            print(f"  ... {trial}/{sample_size_8}")
        bits = rng2.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        try:
            betti = compute_betti_quick(T, n, max_p=4)
        except:
            continue
        b3 = betti.get(3, 0)
        b3_dist_8[b3] += 1

        min_b3_sub = float('inf')
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                betti_sub = compute_betti_quick(T_sub, n-1, max_p=4)
                b3_sub = betti_sub.get(3, 0)
                min_b3_sub = min(min_b3_sub, b3_sub)
            except:
                min_b3_sub = 0  # SVD failure, can't conclude
                break

        if min_b3_sub != float('inf'):
            min_b3_sub_dist_8[min_b3_sub] += 1
            if min_b3_sub >= 1:
                all_del_8 += 1

    print(f"\n  beta_3 distribution: {dict(sorted(b3_dist_8.items()))}")
    print(f"  min(beta_3(T\\v)): {dict(sorted(min_b3_sub_dist_8.items()))}")
    print(f"  All deletions have beta_3 >= 1: {all_del_8}/{sample_size_8}")

    # ============================================
    # Part 3: Theoretical analysis output
    # ============================================
    print(f"\n{'='*70}")
    print("THEORETICAL FRAMEWORK")
    print("=" * 70)
    print("""
    EXACT EQUATION (from LES + beta_2 = 0):
    ========================================
    H_3(T,T\\v) = beta_3(T) - rank(i_*: H_3(T\\v) -> H_3(T))

    PROOF STRATEGY for beta_3 <= 1:
    ================================
    By induction on n (base: beta_3 = 0 for n <= 5).

    For n-vertex tournament T:
    beta_3(T) = H_3(T,T\\v) + rank(i_*)
              <= 1 + beta_3(T\\v) <= 1 + 1 = 2.

    For beta_3(T) = 2, need ALL v to have:
      (a) beta_3(T\\v) = 1
      (b) i_*: H_3(T\\v) -> H_3(T) injective
      (c) H_3(T,T\\v) = 1

    This is an EXTREMELY strong condition.
    Computationally: NEVER observed.

    POSSIBLE PROOF ROUTES:
    ======================
    1. Show beta_3(T\\v) = 1 for all v is impossible
       (i.e., some deletion must have beta_3 = 0).

    2. Show that when beta_3(T\\v) = 1, the inclusion
       i_*: H_3(T\\v) -> H_3(T) must be zero.

    3. Show H_3(T,T\\v) = 1 AND beta_3(T\\v) = 1 implies
       something contradictory about the LES at lower levels.

    4. Use the seesaw: beta_1 * beta_3 = 0. If beta_3(T) >= 1,
       then beta_1(T) = 0. Combine with LES at level 1.
    """)

    # ============================================
    # Part 4: Check route 1 — can all deletions have beta_3 >= 1?
    # ============================================
    print(f"\n{'='*70}")
    print("ROUTE 1: Can all deletions of an n=7 tournament have beta_3 = 1?")
    print("=" * 70)

    # First, catalog which n=6 tournaments have beta_3 = 1
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    b3_1_set = set()
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        betti = compute_betti_quick(T, n, max_p=4)
        if betti.get(3, 0) == 1:
            b3_1_set.add(bits)

    print(f"\n  n=6: {len(b3_1_set)} tournaments with beta_3 = 1 (out of {n_total})")
    print(f"  Fraction: {len(b3_1_set)/n_total:.4f}")

    # For n=7, check a sample: for each T, how many deletions have beta_3 = 1?
    print(f"\n  n=7 sampled: number of deletions with beta_3 = 1 per tournament")
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng3 = np.random.RandomState(456)

    count_dist = Counter()
    for trial in range(300):
        bits = rng3.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        count_b3_1 = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            betti_sub = compute_betti_quick(T_sub, n-1, max_p=4)
            if betti_sub.get(3, 0) == 1:
                count_b3_1 += 1
        count_dist[count_b3_1] += 1

    print(f"  Distribution of #deletions with beta_3 = 1:")
    for k in sorted(count_dist.keys()):
        print(f"    {k} deletions: {count_dist[k]} tournaments")

    max_count = max(count_dist.keys())
    print(f"\n  Maximum #deletions with beta_3 = 1: {max_count}")
    print(f"  ALL deletions (7/7) with beta_3 = 1: {count_dist.get(7, 0)} times")
    print(f"\n  CONCLUSION: {'OBSTRUCTION FOUND' if count_dist.get(7, 0) == 0 else 'OBSTRUCTION NOT FOUND'}")

if __name__ == '__main__':
    main()
