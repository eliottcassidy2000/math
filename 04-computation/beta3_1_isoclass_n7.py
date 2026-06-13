#!/usr/bin/env python3
"""
beta3_1_isoclass_n7.py — opus-2026-03-09-S53

At n=7, how many isomorphism classes have beta_3=1?
For each, what's the deletion pattern?

Key: at n=6, beta_3=1 has exactly 2 iso classes (80+240=320 tours).
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
    for p in range(0, 6):
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
                if len(a_p) - rank == 0:
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

def count_3cycles(T, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    count += 1
                if T[i][k] and T[k][j] and T[j][i]:
                    count += 1
    return count

def main():
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)

    print("=" * 70)
    print("beta_3=1 ISOMORPHISM CLASSES at n=7")
    print("=" * 70)

    b3_1_data = []
    sample_size = 3000

    print(f"\n  Sampling {sample_size} tournaments at n={n}...")

    for trial in range(sample_size):
        if trial % 500 == 0:
            print(f"    ... {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        try:
            b3 = compute_beta3(T, n)
        except:
            continue
        if b3 == 1:
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            c3 = count_3cycles(T, n)

            # Deletion beta_3 pattern
            del_b3 = []
            del_scores = []
            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
                try:
                    b3_sub = compute_beta3(T_sub, n-1)
                except:
                    b3_sub = -1
                del_b3.append(b3_sub)
                sub_scores = tuple(sorted([sum(T_sub[a][b] for b in range(n-1) if b != a) for a in range(n-1)]))
                del_scores.append(sub_scores)

            b3_1_data.append({
                'bits': bits,
                'scores': scores,
                'c3': c3,
                'del_b3': tuple(sorted(del_b3)),
                'del_scores_when_b3_1': [del_scores[i] for i in range(7) if del_b3[i] == 1]
            })

    print(f"\n  Found {len(b3_1_data)} beta_3=1 tournaments")

    # Group by score sequence
    by_score = {}
    for d in b3_1_data:
        by_score.setdefault(d['scores'], []).append(d)

    print(f"\n  Score sequence distribution:")
    for s in sorted(by_score.keys()):
        items = by_score[s]
        c3_counts = Counter(d['c3'] for d in items)
        del_patterns = Counter(d['del_b3'] for d in items)
        print(f"    {s}: {len(items)} tournaments")
        print(f"      3-cycle counts: {dict(sorted(c3_counts.items()))}")
        print(f"      Deletion beta_3 profiles: {dict(sorted(del_patterns.items()))}")

        # What are the deletion score sequences when beta_3(T\v)=1?
        del_score_when_1 = Counter()
        for d in items:
            for ds in d['del_scores_when_b3_1']:
                del_score_when_1[ds] += 1
        if del_score_when_1:
            print(f"      Deletion scores when beta_3=1:")
            for ds, cnt in sorted(del_score_when_1.items(), key=lambda x: -x[1]):
                print(f"        {ds}: {cnt}")

    # Summary: probability of beta_3=1 at n=7
    total_sample = sample_size  # approximate
    print(f"\n  Estimated P(beta_3=1) at n=7: {len(b3_1_data)}/{total_sample} = {len(b3_1_data)/total_sample:.4f}")

if __name__ == '__main__':
    main()
