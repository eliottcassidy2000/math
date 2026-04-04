"""
seesaw_n8_analysis.py — opus-2026-04-04-S1

The seesaw β₁·β₃=0 holds exhaustively at n≤7 but FAILS at n=8.
Question: What is the correct generalization?

Hypotheses to test:
1. β₁ + β₃ ≤ 2 (weaker bound)
2. β₁·β₃ = 0 OR (β₁=1 AND β₃=1) only (exact characterization)
3. β₁·β₃ > 0 implies score sequence is self-complementary
4. β₁·β₃ > 0 requires specific cycle counts (c₃ thresholds)
5. The adjacent-odd seesaw β_{2k-1}·β_{2k+1}=0 fails ONLY for (k=1)

Also compute: β₁·β₃ rate as function of c₃ and score sequence.
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
import random

def make_tournament(n, seed=None):
    """Random tournament on n vertices."""
    if seed is not None:
        rng = np.random.RandomState(seed)
    else:
        rng = np.random
    T = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.rand() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def scores(T, n):
    return tuple(sorted(int(T[i].sum()) for i in range(n)))

def count_3cycles(T, n):
    c3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if ((T[a][b] and T[b][c] and T[c][a]) or
                    (T[a][c] and T[c][b] and T[b][a])):
                    c3 += 1
    return c3

def is_sc_score(scores_tuple):
    """Check if score sequence is self-complementary: d_i + d_{n-1-i} = n-1."""
    n = len(scores_tuple)
    s = sorted(scores_tuple)
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def boundary_matrix_d1(T, n):
    """Compute d₁ boundary matrix for GLMY path homology.
    Basis for Ω₁: directed edges (i→j) where T[i][j]=1.
    Basis for Ω₀: vertices {0,...,n-1}.
    d₁(i→j) = j - i."""
    edges = [(i,j) for i in range(n) for j in range(n) if i != j and T[i][j]]
    ne = len(edges)
    d1 = np.zeros((n, ne), dtype=np.int8)
    for col, (i,j) in enumerate(edges):
        d1[j][col] += 1
        d1[i][col] -= 1
    return d1, edges

def allowed_2paths(T, n):
    """Compute basis for Ω₂: directed 2-paths (i→j→k) where
    T[i][j]=1, T[j][k]=1, T[i][k]=1 (allowed = not a non-allowed face)."""
    paths = []
    for i in range(n):
        for j in range(n):
            if j == i or not T[i][j]: continue
            for k in range(n):
                if k == i or k == j or not T[j][k]: continue
                if T[i][k]:  # allowed: i→k exists
                    paths.append((i,j,k))
    return paths

def boundary_matrix_d2(T, n, edges, paths2):
    """d₂(i→j→k) = (j→k) - (i→k) + (i→j)."""
    edge_idx = {e: idx for idx, e in enumerate(edges)}
    ne = len(edges)
    np2 = len(paths2)
    d2 = np.zeros((ne, np2), dtype=np.int8)
    for col, (i,j,k) in enumerate(paths2):
        d2[edge_idx[(j,k)]][col] += 1
        d2[edge_idx[(i,k)]][col] -= 1
        d2[edge_idx[(i,j)]][col] += 1
    return d2

def allowed_3paths(T, n):
    """Compute basis for Ω₃."""
    paths = []
    for i in range(n):
        for j in range(n):
            if j == i or not T[i][j]: continue
            for k in range(n):
                if k == i or k == j or not T[j][k]: continue
                if not T[i][k]: continue  # need i→k
                for l in range(n):
                    if l == i or l == j or l == k or not T[k][l]: continue
                    if T[i][l] and T[j][l]:  # all forward edges exist
                        paths.append((i,j,k,l))
    return paths

def boundary_matrix_d3(T, n, paths2, paths3):
    """d₃(i→j→k→l) = (j→k→l) - (i→k→l) + (i→j→l) - (i→j→k)."""
    p2_idx = {p: idx for idx, p in enumerate(paths2)}
    np2 = len(paths2)
    np3 = len(paths3)
    d3 = np.zeros((np2, np3), dtype=np.int8)
    for col, (i,j,k,l) in enumerate(paths3):
        if (j,k,l) in p2_idx: d3[p2_idx[(j,k,l)]][col] += 1
        if (i,k,l) in p2_idx: d3[p2_idx[(i,k,l)]][col] -= 1
        if (i,j,l) in p2_idx: d3[p2_idx[(i,j,l)]][col] += 1
        if (i,j,k) in p2_idx: d3[p2_idx[(i,j,k)]][col] -= 1
    return d3

def compute_betti_1_3(T, n):
    """Compute β₁ and β₃ using rank computations mod p (fast)."""
    # β₁ = dim ker(d₁) - rank(d₂)
    # β₃ = dim ker(d₃) - rank(d₄)
    # We use: β₁ = C(n,2) - (n-1) - rank(d₂) [since rank(d₁) = n-1 for tournaments]

    d1, edges = boundary_matrix_d1(T, n)
    paths2 = allowed_2paths(T, n)

    if not paths2:
        beta1 = len(edges) - (n-1)  # no d₂ range
        return beta1, 0

    d2 = boundary_matrix_d2(T, n, edges, paths2)

    # Rank via mod-p Gaussian elimination (p=251 for speed)
    rank_d2 = mod_rank(d2, 251)
    beta1 = len(edges) - (n-1) - rank_d2

    # For β₃, we need d₃ and d₄
    paths3 = allowed_3paths(T, n)
    if not paths3:
        return beta1, 0

    d3 = boundary_matrix_d3(T, n, paths2, paths3)
    rank_d3 = mod_rank(d3, 251)
    ker_d3 = len(paths3) - rank_d3

    # d₄ requires 4-paths... expensive. For β₃, we need rank(d₄).
    # β₃ = ker(d₃) - rank(d₄)
    # We can compute rank(d₄) or use: β₃ = ker(d₃) - im(d₄)

    # For now, compute d₄ if feasible
    paths4 = allowed_4paths(T, n)
    if not paths4:
        return beta1, ker_d3

    d4 = boundary_matrix_d4(T, n, paths3, paths4)
    rank_d4 = mod_rank(d4, 251)
    beta3 = ker_d3 - rank_d4

    return beta1, beta3

def allowed_4paths(T, n):
    """Compute basis for Ω₄: directed 4-paths with all forward edges."""
    paths = []
    for i in range(n):
        for j in range(n):
            if j == i or not T[i][j]: continue
            for k in range(n):
                if k in (i,j) or not T[j][k] or not T[i][k]: continue
                for l in range(n):
                    if l in (i,j,k) or not T[k][l] or not T[i][l] or not T[j][l]: continue
                    for m in range(n):
                        if m in (i,j,k,l): continue
                        if T[l][m] and T[i][m] and T[j][m] and T[k][m]:
                            paths.append((i,j,k,l,m))
    return paths

def boundary_matrix_d4(T, n, paths3, paths4):
    """d₄ boundary."""
    p3_idx = {p: idx for idx, p in enumerate(paths3)}
    np3 = len(paths3)
    np4 = len(paths4)
    d4 = np.zeros((np3, np4), dtype=np.int8)
    for col, (i,j,k,l,m) in enumerate(paths4):
        faces = [(j,k,l,m), (i,k,l,m), (i,j,l,m), (i,j,k,m), (i,j,k,l)]
        signs = [1, -1, 1, -1, 1]
        for face, sign in zip(faces, signs):
            if face in p3_idx:
                d4[p3_idx[face]][col] += sign
    return d4

def mod_rank(M, p):
    """Gaussian elimination over Z/pZ."""
    M = M.astype(np.int64) % p
    rows, cols = M.shape
    rank = 0
    for col in range(cols):
        # Find pivot
        pivot = -1
        for row in range(rank, rows):
            if M[row][col] % p != 0:
                pivot = row
                break
        if pivot == -1:
            continue
        # Swap
        M[[rank, pivot]] = M[[pivot, rank]]
        # Scale
        inv = pow(int(M[rank][col] % p), p-2, p)
        M[rank] = (M[rank] * inv) % p
        # Eliminate
        for row in range(rows):
            if row != rank and M[row][col] % p != 0:
                M[row] = (M[row] - M[row][col] * M[rank]) % p
        rank += 1
    return rank

def main():
    n = 8
    num_samples = 3000
    random.seed(42)

    seesaw_violations = []  # β₁>0 AND β₃>0
    all_data = []

    print(f"Sampling {num_samples} tournaments at n={n}...")
    for trial in range(num_samples):
        if trial % 500 == 0:
            print(f"  trial {trial}/{num_samples}...", flush=True)

        T = make_tournament(n, seed=trial)
        sc = scores(T, n)
        c3 = count_3cycles(T, n)

        beta1, beta3 = compute_betti_1_3(T, n)

        record = {
            'trial': trial,
            'scores': sc,
            'c3': c3,
            'beta1': beta1,
            'beta3': beta3,
            'sc_score': is_sc_score(sc),
            'violation': beta1 > 0 and beta3 > 0,
        }
        all_data.append(record)

        if beta1 > 0 and beta3 > 0:
            seesaw_violations.append(record)

    print(f"\n{'='*60}")
    print(f"SEESAW ANALYSIS at n={n} ({num_samples} samples)")
    print(f"{'='*60}")

    n_b1 = sum(1 for d in all_data if d['beta1'] > 0)
    n_b3 = sum(1 for d in all_data if d['beta3'] > 0)
    n_both = len(seesaw_violations)

    print(f"β₁ > 0: {n_b1}/{num_samples} ({100*n_b1/num_samples:.1f}%)")
    print(f"β₃ > 0: {n_b3}/{num_samples} ({100*n_b3/num_samples:.1f}%)")
    print(f"β₁·β₃ > 0 (violations): {n_both}/{num_samples} ({100*n_both/num_samples:.2f}%)")

    if seesaw_violations:
        print(f"\nViolation details:")
        for v in seesaw_violations[:20]:
            print(f"  trial={v['trial']}: β₁={v['beta1']}, β₃={v['beta3']}, "
                  f"scores={v['scores']}, c₃={v['c3']}, SC_score={v['sc_score']}")

        # Check if all violations have SC score
        all_sc = all(v['sc_score'] for v in seesaw_violations)
        print(f"\nAll violations have SC score: {all_sc}")

        # β₁ and β₃ values in violations
        b1_vals = set(v['beta1'] for v in seesaw_violations)
        b3_vals = set(v['beta3'] for v in seesaw_violations)
        print(f"β₁ values in violations: {b1_vals}")
        print(f"β₃ values in violations: {b3_vals}")

        # c₃ distribution in violations vs general
        viol_c3 = [v['c3'] for v in seesaw_violations]
        all_c3 = [d['c3'] for d in all_data]
        print(f"c₃ in violations: {sorted(set(viol_c3))}")
        print(f"c₃ overall: range [{min(all_c3)}, {max(all_c3)}], mean={np.mean(all_c3):.1f}")

    # β₁ + β₃ bound
    sums = [d['beta1'] + d['beta3'] for d in all_data]
    print(f"\nβ₁ + β₃ max: {max(sums)}")
    print(f"β₁ + β₃ distribution: {dict(sorted(defaultdict(int, {s: sums.count(s) for s in set(sums)}).items()))}")

    # Check what predicts violations
    print(f"\n{'='*60}")
    print(f"VIOLATION PREDICTORS")
    print(f"{'='*60}")

    # Score sequence analysis
    sc_counts = defaultdict(lambda: {'total': 0, 'b1': 0, 'b3': 0, 'both': 0})
    for d in all_data:
        key = d['scores']
        sc_counts[key]['total'] += 1
        if d['beta1'] > 0: sc_counts[key]['b1'] += 1
        if d['beta3'] > 0: sc_counts[key]['b3'] += 1
        if d['violation']: sc_counts[key]['both'] += 1

    # Show score sequences with violations
    viol_scores = [(k, v) for k, v in sc_counts.items() if v['both'] > 0]
    if viol_scores:
        print(f"\nScore sequences producing violations:")
        for sc_seq, counts in sorted(viol_scores, key=lambda x: -x[1]['both']):
            print(f"  {sc_seq}: {counts['both']}/{counts['total']} violations "
                  f"({100*counts['both']/counts['total']:.1f}%), "
                  f"SC={is_sc_score(sc_seq)}")

    # Euler characteristic check
    print(f"\n{'='*60}")
    print(f"EULER CHARACTERISTIC CHECK")
    print(f"{'='*60}")
    # χ = β₀ - β₁ + β₂ - β₃ + ... should be 1 for tournaments
    # Since β₀=1, β₂=0: χ = 1 - β₁ - β₃ + higher
    # For violations: χ = 1 - β₁ - β₃ + β₄ - β₅ + ...
    # We don't compute β₄+ here, but note that χ=1 constrains β₁+β₃
    for v in seesaw_violations[:5]:
        print(f"  trial={v['trial']}: β₁={v['beta1']}, β₃={v['beta3']}, "
              f"β₁+β₃={v['beta1']+v['beta3']} "
              f"(needs β₄-β₅+...= {v['beta1']+v['beta3']} for χ=1)")

if __name__ == "__main__":
    main()
