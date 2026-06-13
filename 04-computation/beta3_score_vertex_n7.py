#!/usr/bin/env python3
"""
beta3_score_vertex_n7.py — opus-2026-03-09-S53

At n=7: which vertex scores guarantee beta_3(T\v) = 0?
If extreme scores (0,1 or 5,6) guarantee this, and every tournament
has at least one such vertex, we're done.

Sampling approach at n=7.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

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

def main():
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(789)

    print("=" * 70)
    print("SCORE-VERTEX RELATION AT n=7")
    print("=" * 70)

    score_b3_relation = defaultdict(Counter)
    score_b3_given_T_b3_1 = defaultdict(Counter)
    sample_size = 1000

    print(f"\n  Sampling {sample_size} tournaments at n={n}...")

    for trial in range(sample_size):
        if trial % 200 == 0:
            print(f"    ... {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]

        try:
            b3_T = compute_beta3(T, n)
        except:
            continue

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                b3_sub = compute_beta3(T_sub, n-1)
            except:
                continue
            score_b3_relation[scores[v]][b3_sub] += 1
            if b3_T == 1:
                score_b3_given_T_b3_1[scores[v]][b3_sub] += 1

    print(f"\n  Score of v → beta_3(T\\v) distribution at n=7:")
    print(f"  {'Score':>6} | beta_3=0     | beta_3=1     | Total      | P(beta_3=1)")
    print(f"  {'-'*65}")
    for score in sorted(score_b3_relation.keys()):
        d = score_b3_relation[score]
        total = sum(d.values())
        b0 = d.get(0, 0)
        b1 = d.get(1, 0)
        p1 = b1/total if total > 0 else 0
        print(f"  {score:>6} | {b0:>12} | {b1:>12} | {total:>10} | {p1:.6f}")

    # Which scores guarantee beta_3(T\v) = 0?
    guaranteed_good = []
    for score in sorted(score_b3_relation.keys()):
        if score_b3_relation[score].get(1, 0) == 0:
            guaranteed_good.append(score)

    print(f"\n  Scores GUARANTEEING beta_3(T\\v)=0 (in sample): {guaranteed_good}")

    if score_b3_given_T_b3_1:
        print(f"\n  Same, for beta_3(T)=1 tournaments only:")
        print(f"  {'Score':>6} | beta_3=0     | beta_3=1     | P(beta_3=1)")
        print(f"  {'-'*55}")
        for score in sorted(score_b3_given_T_b3_1.keys()):
            d = score_b3_given_T_b3_1[score]
            total = sum(d.values())
            b0 = d.get(0, 0)
            b1 = d.get(1, 0)
            p1 = b1/total if total > 0 else 0
            print(f"  {score:>6} | {b0:>12} | {b1:>12} | {p1:.6f}")

    # Check: does deleting a SOURCE (score=6) or SINK (score=0) always give beta_3=0?
    # And: what's the score sequence of T\v when beta_3(T\v)=1?
    print(f"\n{'='*70}")
    print("SCORE SEQUENCE OF T\\v WHEN beta_3(T\\v)=1")
    print("=" * 70)

    b3_1_del_scores = Counter()
    for trial in range(2000):
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            try:
                b3_sub = compute_beta3(T_sub, n-1)
            except:
                continue
            if b3_sub == 1:
                sub_scores = tuple(sorted([sum(T_sub[a][b] for b in range(n-1) if b != a) for a in range(n-1)]))
                parent_score = sum(T[v][j] for j in range(n) if j != v)
                b3_1_del_scores[(parent_score, sub_scores)] += 1

    print(f"  (score of v, score seq of T\\v) when beta_3(T\\v)=1:")
    for key, cnt in sorted(b3_1_del_scores.items(), key=lambda x: -x[1]):
        print(f"    v_score={key[0]}, T\\v scores={key[1]}: {cnt}")

if __name__ == '__main__':
    main()
