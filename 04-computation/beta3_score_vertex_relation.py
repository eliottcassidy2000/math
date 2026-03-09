#!/usr/bin/env python3
"""
beta3_score_vertex_relation.py — opus-2026-03-09-S53

KEY QUESTION: Is beta_3(T\v) related to the score (out-degree) of v?

If we can show that for extreme-score vertices (sources/sinks or near),
beta_3(T\v) = 0, then the good vertex property follows since every
tournament has vertices with different scores.

At n=6 exhaustive: check all 32768 tournaments.
For each vertex v, record (score of v, beta_3(T\v)).
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
    print("=" * 70)
    print("SCORE-VERTEX RELATION: does score of v predict beta_3(T\\v)?")
    print("=" * 70)

    # n=6 exhaustive
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    # For each (score_of_v, beta_3_of_T, beta_3_of_T_minus_v), count occurrences
    score_b3_relation = defaultdict(Counter)  # (score_v) -> Counter of beta_3(T\v)
    score_b3_when_b3_1 = defaultdict(Counter)  # same but only when beta_3(T)=1

    print(f"\n  Exhaustive enumeration at n={n} ({n_total} tournaments)...")

    for bits in range(n_total):
        if bits % 5000 == 0:
            print(f"    ... {bits}/{n_total}")
        T = tournament_from_bits(n, bits)
        scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]

        # beta_3(T)
        b3_T = compute_beta3(T, n)

        for v in range(n):
            remaining = [i for i in range(n) if i != v]
            T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            b3_sub = compute_beta3(T_sub, n-1)

            score_b3_relation[scores[v]][b3_sub] += 1
            if b3_T == 1:
                score_b3_when_b3_1[scores[v]][b3_sub] += 1

    print(f"\n  Score of v → beta_3(T\\v) distribution:")
    print(f"  {'Score':>6} | beta_3=0     | beta_3=1     | P(beta_3=1)")
    print(f"  {'-'*55}")
    for score in sorted(score_b3_relation.keys()):
        d = score_b3_relation[score]
        total = sum(d.values())
        b0 = d.get(0, 0)
        b1 = d.get(1, 0)
        p1 = b1/total if total > 0 else 0
        print(f"  {score:>6} | {b0:>12} | {b1:>12} | {p1:.6f}")

    print(f"\n  Same, restricted to beta_3(T)=1 tournaments:")
    print(f"  {'Score':>6} | beta_3=0     | beta_3=1     | P(beta_3=1)")
    print(f"  {'-'*55}")
    for score in sorted(score_b3_when_b3_1.keys()):
        d = score_b3_when_b3_1[score]
        total = sum(d.values())
        b0 = d.get(0, 0)
        b1 = d.get(1, 0)
        p1 = b1/total if total > 0 else 0
        print(f"  {score:>6} | {b0:>12} | {b1:>12} | {p1:.6f}")

    # Check: which scores GUARANTEE beta_3(T\v) = 0?
    guaranteed_good = []
    for score in sorted(score_b3_relation.keys()):
        if score_b3_relation[score].get(1, 0) == 0:
            guaranteed_good.append(score)

    print(f"\n  Scores that GUARANTEE beta_3(T\\v) = 0: {guaranteed_good}")
    print(f"  (These are 'automatically good' vertices)")

    if guaranteed_good:
        # Check: does every tournament have at least one vertex with a guaranteed-good score?
        violations = 0
        for bits in range(n_total):
            T = tournament_from_bits(n, bits)
            scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
            if not any(s in guaranteed_good for s in scores):
                violations += 1
        print(f"  Tournaments with NO guaranteed-good vertex: {violations}/{n_total}")

    # Also: for beta_3(T)=1 tournaments, what are the SCORES?
    print(f"\n  Score sequences of beta_3=1 tournaments at n=6:")
    b3_1_scores = Counter()
    for bits in range(n_total):
        T = tournament_from_bits(n, bits)
        b3 = compute_beta3(T, n)
        if b3 == 1:
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            b3_1_scores[scores] += 1
    for s, c in sorted(b3_1_scores.items()):
        print(f"    {s}: {c}")

if __name__ == '__main__':
    main()
