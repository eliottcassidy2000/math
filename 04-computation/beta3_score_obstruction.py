#!/usr/bin/env python3
"""
beta3_score_obstruction.py — opus-2026-03-09-S53

GOAL: Prove beta_3(T)=2 is impossible by SCORE SEQUENCE analysis alone.

From S52: if beta_3(T)=2, ALL n vertex-deletions must have beta_3(T\v)=1.
At n=6, beta_3=1 occurs ONLY for score sequences:
  (1,1,1,4,4,4) — 80 tournaments
  (2,2,2,3,3,3) — 240 tournaments

QUESTION: Can a 7-vertex tournament T have ALL 7 deletions with score
sequences in {(1,1,1,4,4,4), (2,2,2,3,3,3)}?

This is a purely combinatorial question about score sequences.
If the answer is NO, then beta_3=2 is impossible at n=7.

By induction: at n=8, if beta_3=2, all deletions have beta_3=1.
Their deletions (at n=7) must also have beta_3 ≤ 1, so the n=7
constraint propagates up.

APPROACH: Enumerate all possible 7-vertex score sequences and check
if any is compatible with all deletions being in the target set.
"""

from itertools import combinations, product
from collections import Counter
import numpy as np

def check_score_compatibility_n7():
    """
    For a tournament T on 7 vertices with scores (d_0,...,d_6):
    - Sum = 21
    - 0 ≤ d_i ≤ 6
    - Removing vertex i: remaining scores are d_j - T[j][i] for j ≠ i
    - T[j][i] ∈ {0,1}, and exactly d_i of the 6 edges from i go OUT (T[i][j]=1)
    - So removing i: d_i vertices keep their score, 6-d_i vertices lose 1

    For the deletion of i to have sorted scores in {(1,1,1,4,4,4), (2,2,2,3,3,3)}:
    we need a specific partition of the remaining vertices into "kept" and "decreased".
    """
    target_scores = {(1,1,1,4,4,4), (2,2,2,3,3,3)}

    print("=" * 70)
    print("SCORE SEQUENCE OBSTRUCTION AT n=7")
    print("=" * 70)

    # Enumerate all 7-vertex tournament score sequences
    # Score sequence is a non-decreasing tuple summing to 21
    # with each entry in [0,6]
    score_seqs = []
    def gen_scores(pos, remaining, min_val, current):
        if pos == 7:
            if remaining == 0:
                score_seqs.append(tuple(current))
            return
        max_val = min(6, remaining - (6 - pos) * 0)  # upper bound
        min_bound = max(min_val, 0)
        # remaining must be achievable: each remaining position gets ≥ min_val
        # and ≤ 6
        remaining_positions = 7 - pos
        if remaining < min_bound * remaining_positions:
            return
        if remaining > 6 * remaining_positions:
            return
        for v in range(min_bound, max_val + 1):
            gen_scores(pos + 1, remaining - v, v, current + [v])

    gen_scores(0, 21, 0, [])
    print(f"\n  Total 7-vertex tournament score sequences: {len(score_seqs)}")

    # For each score sequence, check if ALL 7 deletions COULD have
    # sorted scores in target_scores.
    #
    # When we delete vertex with score d_i:
    # - d_i vertices among the remaining 6 are "kept" (they beat vertex i)
    # - 6-d_i vertices are "decreased" (they lost to vertex i)
    # - The kept vertices have scores d_j, decreased have d_j - 1
    # - The sorted result must be in target_scores
    #
    # Note: we don't know the EXACT assignment (which specific vertices
    # are kept vs decreased). We need to check if ANY valid assignment exists.
    # But we can check necessary conditions from the multiset of remaining scores.

    compatible = []

    for scores in score_seqs:
        all_deletions_ok = True

        for i in range(7):
            d_i = scores[i]
            remaining = list(scores[:i]) + list(scores[i+1:])

            # When deleting vertex i with degree d_i:
            # d_i of the 6 remaining vertices beat i → score unchanged
            # 6-d_i remaining vertices lost to i → score decreased by 1

            # We need to choose which d_i vertices are "kept" and which 6-d_i are "decreased"
            # The result sorted must be in target_scores

            # Try all ways to partition remaining 6 vertices into
            # d_i "kept" and (6-d_i) "decreased"
            found_valid = False

            for kept_indices in combinations(range(6), d_i):
                result = list(remaining)
                for j in range(6):
                    if j not in kept_indices:
                        result[j] -= 1
                result_sorted = tuple(sorted(result))
                if result_sorted in target_scores:
                    found_valid = True
                    break

            if not found_valid:
                all_deletions_ok = False
                break

        if all_deletions_ok:
            compatible.append(scores)

    print(f"\n  Score sequences compatible with ALL deletions having beta_3=1:")
    print(f"  Count: {len(compatible)}")
    if compatible:
        for s in compatible:
            print(f"    {s}")
    else:
        print(f"    NONE — beta_3=2 is IMPOSSIBLE at n=7 by score constraint!")

    return compatible

def check_score_compatibility_n8():
    """Same analysis at n=8."""
    # At n=7, the beta_3=1 score sequences need to be determined.
    # From S52 data: sampled, but not exhaustive.
    # However, we can check a necessary condition: if beta_3(T)=2 at n=8,
    # then all 8 deletions have beta_3=1 at n=7.
    # From the n=7 analysis above, if no score sequence is compatible at n=7,
    # then beta_3=2 is impossible at n=7.
    # But beta_3=1 at n=7 CAN occur (38/500 in sampling).
    # We'd need the n=7 beta_3=1 score sequences to do this analysis at n=8.

    print(f"\n{'='*70}")
    print("SCORE SEQUENCE ANALYSIS AT n=7 (beta_3=1 characterization)")
    print("=" * 70)

    # Let me check: what score sequences give beta_3=1 at n=7?
    # This requires actual computation (not just scores).
    # From S52 sampling (500 tournaments): 38 with beta_3=1.
    # Let me compute the score sequences of those.

    from itertools import permutations

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
        """Quick beta_3 computation."""
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

    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)

    # Collect score sequences of beta_3=1 tournaments at n=7
    b3_1_scores = Counter()
    total_checked = 0
    sample_size = 300

    print(f"\n  Sampling n=7 tournaments to find beta_3=1 score sequences...")
    for trial in range(sample_size):
        if trial % 50 == 0:
            print(f"    ... {trial}/{sample_size}")
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)
        total_checked += 1
        try:
            b3 = compute_beta3(T, n)
        except:
            continue
        if b3 == 1:
            scores = tuple(sorted([sum(T[i][j] for j in range(n) if j != i) for i in range(n)]))
            b3_1_scores[scores] += 1

    print(f"\n  Total checked: {total_checked}")
    print(f"  beta_3=1 score sequences at n=7:")
    for s, c in sorted(b3_1_scores.items(), key=lambda x: -x[1]):
        print(f"    {s}: {c} tournaments")

    # Now check: at n=8, can all 8 deletions have beta_3=1?
    # The deletions must have score sequences from the above set.
    if b3_1_scores:
        n7_target = set(b3_1_scores.keys())
        print(f"\n  n=7 beta_3=1 score sequences: {n7_target}")

        # Generate 8-vertex score sequences
        score_seqs_8 = []
        def gen_scores_8(pos, remaining, min_val, current):
            if pos == 8:
                if remaining == 0:
                    score_seqs_8.append(tuple(current))
                return
            max_val = min(7, remaining)
            remaining_positions = 8 - pos
            if remaining < min_val * remaining_positions or remaining > 7 * remaining_positions:
                return
            for v in range(min_val, max_val + 1):
                gen_scores_8(pos + 1, remaining - v, v, current + [v])

        gen_scores_8(0, 28, 0, [])
        print(f"  Total 8-vertex score sequences: {len(score_seqs_8)}")

        compatible_8 = []
        for scores in score_seqs_8:
            all_ok = True
            for i in range(8):
                d_i = scores[i]
                remaining = list(scores[:i]) + list(scores[i+1:])
                found = False
                for kept_indices in combinations(range(7), d_i):
                    result = list(remaining)
                    for j in range(7):
                        if j not in kept_indices:
                            result[j] -= 1
                    result_sorted = tuple(sorted(result))
                    if result_sorted in n7_target:
                        found = True
                        break
                if not found:
                    all_ok = False
                    break
            if all_ok:
                compatible_8.append(scores)

        print(f"\n  n=8 score sequences compatible with all deletions having beta_3=1:")
        print(f"  Count: {len(compatible_8)}")
        for s in compatible_8[:10]:
            print(f"    {s}")

def main():
    compatible_n7 = check_score_compatibility_n7()
    check_score_compatibility_n8()

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("=" * 70)
    if not compatible_n7:
        print("""
  THEOREM: beta_3(T) ≠ 2 for any tournament T on n=7 vertices.

  PROOF (by score sequence obstruction):
  If beta_3(T)=2, then by the exact LES equation (HYP-359) with beta_2=0:
    H_3(T,T\\v) = beta_3(T) - rank(i_*) for all v.
  Since H_3(T,T\\v) <= 1 (HYP-351) and beta_3(T\\v) <= 1 (induction):
    2 = H_3 + rank(i_*) requires H_3=1 and rank(i_*)=1,
    hence beta_3(T\\v)=1 for ALL v.

  But at n=6, beta_3=1 only for score sequences (1,1,1,4,4,4) and (2,2,2,3,3,3).
  No 7-vertex score sequence produces ALL 7 deletions with these scores.
  QED.

  This proof extends by induction: if beta_3 <= 1 at n <= k (proved),
  then at n=k+1, beta_3=2 requires all deletions beta_3=1, which restricts
  the score sequences. If no compatible score exists, beta_3 <= 1 at n=k+1.
        """)

if __name__ == '__main__':
    main()
