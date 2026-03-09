#!/usr/bin/env python3
"""
beta3_1_structure_n6.py — opus-2026-03-09-S53

Exhaustive characterization of beta_3=1 tournaments at n=6.

We know there are exactly 320:
  - 80 with score (1,1,1,4,4,4)
  - 240 with score (2,2,2,3,3,3)

Questions:
1. How many isomorphism classes?
2. What is the 3-cycle structure?
3. What is the allowed path structure at level 3?
4. Is there a simple graph-theoretic characterization?

The (1,1,1,4,4,4) tournaments: 3 "losers" (score 1) and 3 "winners" (score 4).
Each loser beats exactly 1 other vertex, each winner loses to exactly 1.
This is: take a bipartite tournament between {losers} and {winners} where
winners beat all losers, plus internal tournaments on each group.

The (2,2,2,3,3,3) tournaments: balanced bipartition structure?
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
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    count += 1
                if T[i][k] and T[k][j] and T[j][i]:
                    count += 1
    return count

def is_isomorphic(T1, T2, n):
    """Check if T1 and T2 are isomorphic."""
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if T1[i][j] != T2[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def canonical_form(T, n):
    """Get a canonical form for isomorphism class (first lex adjacency matrix under permutation)."""
    best = None
    for perm in permutations(range(n)):
        # Build permuted adjacency
        adj = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or adj < best:
            best = adj
    return best

def main():
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    print("=" * 70)
    print("STRUCTURAL CHARACTERIZATION of beta_3=1 tournaments at n=6")
    print("=" * 70)

    b3_1_tournaments = []
    b3_1_by_score = {}

    for bits in range(n_total):
        if bits % 5000 == 0:
            print(f"  Scanning {bits}/{n_total}...", flush=True)
        T = tournament_from_bits(n, bits)
        try:
            b3 = compute_beta3(T, n)
        except:
            continue
        if b3 == 1:
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            b3_1_tournaments.append((bits, T, scores))
            b3_1_by_score.setdefault(scores, []).append((bits, T))

    print(f"\n  Total beta_3=1 tournaments: {len(b3_1_tournaments)}")

    # 3-cycle counts
    cycle_counts = Counter()
    for bits, T, scores in b3_1_tournaments:
        c3 = count_3cycles(T, n)
        cycle_counts[(scores, c3)] += 1

    print(f"\n  3-cycle counts by score:")
    for (s, c), cnt in sorted(cycle_counts.items()):
        print(f"    score={s}, 3-cycles={c}: {cnt} tournaments")

    # Isomorphism classes
    print(f"\n  Computing isomorphism classes...")
    iso_classes = {}
    for scores in sorted(b3_1_by_score.keys()):
        print(f"\n  Score {scores}:")
        tournaments = b3_1_by_score[scores]
        classes = []
        for bits, T in tournaments:
            canon = canonical_form(T, n)
            found = False
            for cls_canon, cls_bits in classes:
                if canon == cls_canon:
                    found = True
                    break
            if not found:
                classes.append((canon, bits))
        iso_classes[scores] = classes
        print(f"    {len(tournaments)} tournaments in {len(classes)} isomorphism classes")

        # For each class, show a representative
        for idx, (canon, bits) in enumerate(classes):
            T = tournament_from_bits(n, bits)
            scores_v = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]
            c3 = count_3cycles(T, n)
            # Count allowed paths by level
            path_counts = []
            for p in range(6):
                count = 0
                if p + 1 <= n:
                    for verts in combinations(range(n), p+1):
                        for perm in permutations(verts):
                            if is_allowed_path(T, perm):
                                count += 1
                path_counts.append(count)

            print(f"    Class {idx+1}: 3-cycles={c3}, vertex scores={scores_v}")
            print(f"      Allowed paths by level: {path_counts}")
            print(f"      Adjacency matrix:")
            for row in T:
                print(f"        {row}")

    # Check: are the (1,1,1,4,4,4) tournaments "near-transitive"?
    # A transitive tournament has score sequence (0,1,2,3,4,5).
    # The (1,1,1,4,4,4) ones have 3 winners and 3 losers.
    print(f"\n{'='*70}")
    print("STRUCTURAL ANALYSIS")
    print("=" * 70)

    # For (1,1,1,4,4,4): partition into low={score 1} and high={score 4}
    for scores_key in [(1,1,1,4,4,4), (2,2,2,3,3,3)]:
        print(f"\n  Score {scores_key}:")
        if scores_key not in b3_1_by_score:
            continue
        bits0, T0 = b3_1_by_score[scores_key][0]
        scores_v = [sum(T0[i][j] for j in range(n) if j != i) for i in range(n)]

        if scores_key == (1,1,1,4,4,4):
            low = [i for i in range(n) if scores_v[i] == 1]
            high = [i for i in range(n) if scores_v[i] == 4]
            print(f"    Low (score 1): {low}")
            print(f"    High (score 4): {high}")
            # Arcs between groups
            for i in low:
                for j in high:
                    print(f"      T[{i}][{j}] = {T0[i][j]} (low→high)")
            # Internal arcs
            print(f"    Internal arcs among low:")
            for i in low:
                for j in low:
                    if i != j:
                        print(f"      T[{i}][{j}] = {T0[i][j]}")
            print(f"    Internal arcs among high:")
            for i in high:
                for j in high:
                    if i != j:
                        print(f"      T[{i}][{j}] = {T0[i][j]}")

        elif scores_key == (2,2,2,3,3,3):
            low = [i for i in range(n) if scores_v[i] == 2]
            high = [i for i in range(n) if scores_v[i] == 3]
            print(f"    Low (score 2): {low}")
            print(f"    High (score 3): {high}")
            for i in low:
                for j in high:
                    print(f"      T[{i}][{j}] = {T0[i][j]} (low→high)")

if __name__ == '__main__':
    main()
