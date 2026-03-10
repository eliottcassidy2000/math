"""
beta3_n9_analysis.py — Analyze the beta_3=2 tournament found at n=9

The search (seed 54321, trial 1946) found exactly 1 beta_3=2 tournament
in 3000 random n=9 samples. Score: [2,3,3,3,4,5,5,5,6].

Key questions:
1. How many good/bad vertices? (Is it all-good like the n=8 SC type?)
2. What is the deletion pattern? (b3(T\v) for each v)
3. What is H_3(T,T\v) for each vertex?
4. Compare with n=8 findings (Claim II failure)

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import gc
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, adj_to_bits
from beta3_lean import fast_beta3_lean


def main():
    print("=" * 70)
    print("BETA_3=2 AT n=9: STRUCTURAL ANALYSIS")
    print("=" * 70)

    n = 9

    # Reproduce the exact tournament from the search
    rng = np.random.RandomState(54321)
    target_trial = 1946

    print(f"\nReproducing trial {target_trial}...")
    t0 = time.time()
    A = None
    for trial in range(target_trial + 1):
        A = random_tournament(n, rng)
        if trial == target_trial:
            break

    b3 = fast_beta3_lean(A, n)
    scores = sorted([int(sum(A[i])) for i in range(n)])
    bits = adj_to_bits(A, n)
    print(f"  Trial {target_trial}: beta_3 = {b3}, scores = {scores}")
    print(f"  Bits encoding: {bits}")
    print(f"  Reproduced in {time.time()-t0:.1f}s")

    if b3 != 2:
        print("ERROR: Did not reproduce beta_3=2!")
        return

    # Vertex deletion analysis
    print(f"\n{'='*60}")
    print("VERTEX DELETION ANALYSIS")
    print(f"{'='*60}")

    good_vertices = []
    bad_vertices = []
    b3_deletions = []

    for v in range(n):
        remaining = [i for i in range(n) if i != v]
        A_sub = np.array([[A[remaining[i]][remaining[j]]
                          for j in range(n-1)] for i in range(n-1)], dtype=np.int8)
        gc.collect()
        b3_v = fast_beta3_lean(A_sub, n-1)
        out_deg = int(sum(A[v]))
        in_deg = n - 1 - out_deg

        tag = "GOOD" if b3_v == 0 else "BAD"
        if b3_v == 0:
            good_vertices.append(v)
        else:
            bad_vertices.append(v)
        b3_deletions.append(b3_v)

        print(f"  v={v} (out={out_deg}, in={in_deg}): beta_3(T\\v) = {b3_v} [{tag}]")

    print(f"\n  Good vertices: {good_vertices} ({len(good_vertices)}/{n})")
    print(f"  Bad vertices:  {bad_vertices} ({len(bad_vertices)}/{n})")
    print(f"  Deletion beta_3 values: {b3_deletions}")

    # For good vertices: H_3(T,T\v) = beta_3(T) = 2 (since i_* = 0)
    # For bad vertices: H_3(T,T\v) = beta_3(T) - rank(i_*)
    # If ALL vertices are good, Claim II (H_3^rel <= 1) fails maximally
    n_good = len(good_vertices)
    if n_good == n:
        print(f"\n  *** ALL {n} VERTICES ARE GOOD ***")
        print(f"  => H_3(T,T\\v) = 2 for ALL vertices")
        print(f"  => Claim II FAILS at n={n} (same as all-good n=8 case)")
    elif n_good > 0:
        print(f"\n  Mixed case: {n_good} good, {n - n_good} bad")
        print(f"  For good vertices: H_3(T,T\\v) = beta_3(T) = 2")
        print(f"  For bad vertices: need rank(i_*) computation")

    # Double-deletion analysis (for good vertices)
    if n_good <= 4:  # Only do this for a few vertices to save time
        print(f"\n{'='*60}")
        print("DOUBLE DELETION: b3(T\\v\\w) for good vertices")
        print(f"{'='*60}")

        for v in good_vertices[:3]:
            remaining_v = [i for i in range(n) if i != v]
            A_v = np.array([[A[remaining_v[i]][remaining_v[j]]
                            for j in range(n-1)] for i in range(n-1)], dtype=np.int8)

            print(f"\n  Vertex v={v} (out={int(sum(A[v]))}): beta_3(T\\v) = 0")
            for w_idx in range(n-1):
                w = remaining_v[w_idx]
                remaining_vw = [i for i in range(n-1) if i != w_idx]
                A_vw = np.array([[A_v[remaining_vw[i]][remaining_vw[j]]
                                 for j in range(n-2)] for i in range(n-2)], dtype=np.int8)
                gc.collect()
                b3_vw = fast_beta3_lean(A_vw, n-2)
                out_w = int(sum(A[w])) - (1 if A[w][v] else 0)
                if b3_vw > 0:
                    print(f"    T\\{v}\\{w}: beta_3 = {b3_vw} (w out_deg_in_T\\v = {out_w})")

    # Score sequence analysis
    print(f"\n{'='*60}")
    print("SCORE SEQUENCE ANALYSIS")
    print(f"{'='*60}")

    score_seq = sorted([int(sum(A[i])) for i in range(n)])
    print(f"  Score sequence: {score_seq}")
    print(f"  Sum: {sum(score_seq)} (expected: {n*(n-1)//2})")

    # Check self-complementary
    is_sc = all(score_seq[i] + score_seq[n-1-i] == n-1 for i in range(n))
    print(f"  Self-complementary score: {is_sc}")

    score_var = np.var(score_seq)
    print(f"  Score variance: {score_var:.3f}")

    # 3-cycle count
    c3 = 0
    for triple in __import__('itertools').combinations(range(n), 3):
        a, b, c = triple
        if (A[a][b] + A[b][c] + A[c][a] == 3) or (A[a][c] + A[c][b] + A[b][a] == 3):
            c3 += 1
    print(f"  3-cycle count: c3 = {c3}")

    # Compare with n=8 beta_3=2 types
    print(f"\n{'='*60}")
    print("COMPARISON WITH n=8 FINDINGS")
    print(f"{'='*60}")
    print(f"  n=8 all-good type: score (3,3,3,3,4,4,4,4), SC, all 8 good")
    print(f"  n=8 mixed type: score (2,3,3,3,4,4,4,5), not SC, 6 good + 2 bad")
    print(f"  n=9 this case: score {score_seq}, SC={is_sc}, {n_good} good + {n-n_good} bad")

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  beta_3(T) = {b3} at n={n}")
    print(f"  Score: {score_seq} (SC: {is_sc})")
    print(f"  Good/Bad: {n_good}/{n-n_good}")
    if n_good == n:
        print(f"  => Claim II fails: H_3(T,T\\v) = 2 for ALL v")
    else:
        print(f"  => Mixed good/bad pattern")
    print(f"  Frequency: 1/3000 = 0.033% (vs 4/5000 = 0.08% at n=8)")
    print(f"{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
