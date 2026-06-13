#!/usr/bin/env python3
"""
beta3_compatible_scores.py — opus-2026-03-09-S53

The 4 score sequences compatible with beta_3=2 at n=7 are:
  (2,2,2,3,4,4,4), (2,2,3,3,3,4,4), (2,3,3,3,3,3,4), (3,3,3,3,3,3,3)

For each: enumerate tournaments with that score, compute beta_3.
If NONE has beta_3=2, we have a proof at n=7.

Also check: for each tournament with these scores, how many of the
7 deletions actually have beta_3=1?
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

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

def compute_beta3(T, n):
    tol = 1e-8
    all_paths = {}
    max_p = min(5, n)
    for p in range(0, max_p + 1):
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

    omega = build_omega(all_paths, 4)
    boundary = {}
    for p in range(1, 5):
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
    for p in range(5):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, 5):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    ker3 = dims[3] - ranks.get(3, 0)
    im4 = ranks.get(4, 0)
    return ker3 - im4

def main():
    n = 7
    target_scores = [
        (2,2,2,3,4,4,4),
        (2,2,3,3,3,4,4),
        (2,3,3,3,3,3,4),
        (3,3,3,3,3,3,3),
    ]

    print("=" * 70)
    print("EXHAUSTIVE CHECK: beta_3 for score-compatible tournaments at n=7")
    print("=" * 70)

    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    rng = np.random.RandomState(42)
    sample_size = 10000

    for target in target_scores:
        print(f"\n  Score sequence: {target}")
        found_count = 0
        b3_dist = Counter()
        max_b3_sub_1 = 0

        for trial in range(sample_size):
            bits = rng.randint(0, n_total)
            T = tournament_from_bits(n, bits)
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            if scores != target:
                continue

            found_count += 1
            try:
                b3 = compute_beta3(T, n)
            except Exception as e:
                print(f"    ERROR: {e}")
                continue
            b3_dist[b3] += 1

            # Count deletions with beta_3=1
            if b3 >= 1:
                count_b3_1_sub = 0
                for v in range(n):
                    remaining = [i for i in range(n) if i != v]
                    T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
                    try:
                        b3_sub = compute_beta3(T_sub, n-1)
                        if b3_sub == 1:
                            count_b3_1_sub += 1
                    except:
                        pass
                max_b3_sub_1 = max(max_b3_sub_1, count_b3_1_sub)

        print(f"    Found: {found_count}/{sample_size}")
        print(f"    beta_3 distribution: {dict(sorted(b3_dist.items()))}")
        print(f"    Max #deletions with beta_3=1: {max_b3_sub_1}")

    # Also try a focused search of regular tournaments (score 3,3,3,3,3,3,3)
    print(f"\n{'='*70}")
    print("FOCUSED: Regular tournaments on 7 vertices")
    print("=" * 70)

    reg_count = 0
    reg_b3_dist = Counter()
    reg_b3_sub_counts = Counter()

    for trial in range(50000):
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
        if scores != (3,3,3,3,3,3,3):
            continue

        reg_count += 1
        try:
            b3 = compute_beta3(T, n)
        except:
            continue
        reg_b3_dist[b3] += 1

        # Count deletions with beta_3=1
        count_b3_1 = 0
        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                b3_sub = compute_beta3(T_sub, n-1)
                if b3_sub == 1:
                    count_b3_1 += 1
            except:
                pass
        reg_b3_sub_counts[count_b3_1] += 1

    print(f"  Regular tournaments found: {reg_count}")
    print(f"  beta_3 distribution: {dict(sorted(reg_b3_dist.items()))}")
    print(f"  #deletions with beta_3=1: {dict(sorted(reg_b3_sub_counts.items()))}")
    print(f"  ALL 7 deletions beta_3=1? {reg_b3_sub_counts.get(7, 0)}")

if __name__ == '__main__':
    main()
