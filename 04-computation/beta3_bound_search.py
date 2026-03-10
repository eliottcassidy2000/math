"""
beta3_bound_search.py — Search for max beta_3 at n=8,9,10

Key question: what replaces beta_3 <= 1 at n >= 8?
- Is it beta_3 <= 2?
- Does it grow with n?
- Is beta_3 <= floor(n/4)?

Uses large random sampling with mod-p exact arithmetic.

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp


def search_max_beta3(n, n_trials, seed=42):
    """Search for maximum beta_3 across random tournaments."""
    print(f"\n--- n={n}: sampling {n_trials} tournaments (seed={seed}) ---")
    rng = np.random.RandomState(seed)

    b3_dist = {}
    max_b3 = 0
    max_b3_scores = []
    total_nonzero = 0

    t0 = time.time()
    for trial in range(n_trials):
        A = random_tournament(n, rng)
        res = full_chain_complex_modp(A, n, max_p=5)
        b3 = res['bettis'].get(3, 0)

        b3_dist[b3] = b3_dist.get(b3, 0) + 1
        if b3 > 0:
            total_nonzero += 1

        if b3 > max_b3:
            max_b3 = b3
            max_b3_scores = [tuple(sorted([int(sum(A[i])) for i in range(n)]))]
            # Also get full profile
            if b3 > 1:
                res_full = full_chain_complex_modp(A, n, max_p=7)
                profile = tuple(res_full['bettis'].get(p, 0) for p in range(n))
                print(f"  NEW MAX beta_3={b3} at trial {trial}! "
                      f"score={max_b3_scores[-1]}, profile={profile}")
        elif b3 == max_b3 and b3 > 1:
            score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            max_b3_scores.append(score)

        if (trial + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (trial + 1) / elapsed
            print(f"  {trial+1}/{n_trials}, {elapsed:.1f}s ({rate:.1f}/s), "
                  f"max_b3={max_b3}, b3>0: {total_nonzero}", flush=True)

    elapsed = time.time() - t0
    print(f"\n  n={n} RESULTS ({n_trials} tours, {elapsed:.1f}s):")
    print(f"    max beta_3 = {max_b3}")
    print(f"    b3 distribution: {dict(sorted(b3_dist.items()))}")
    print(f"    b3 > 0: {total_nonzero}/{n_trials} ({100*total_nonzero/n_trials:.1f}%)")
    if max_b3 > 1:
        print(f"    b3={max_b3} scores: {max_b3_scores[:5]}")

    return max_b3, b3_dist


def main():
    print("=" * 70)
    print("BETA_3 BOUND SEARCH — What replaces beta_3 <= 1 at n >= 8?")
    print("=" * 70)

    # n=8: large sample to get better statistics
    search_max_beta3(8, 10000, seed=99999)

    # n=9: moderate sample
    search_max_beta3(9, 2000, seed=99999)

    # n=10: small sample (slow)
    search_max_beta3(10, 200, seed=99999)

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
