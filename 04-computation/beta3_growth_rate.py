"""
beta3_growth_rate.py — Investigate beta_3 growth rate across n=6-9

Key question: does max(beta_3) stay bounded at 2, or grow with n?

Known data:
  n=6: max beta_3 = 1 (exhaustive, 32768 tournaments)
  n=7: max beta_3 = 1 (exhaustive, 2097152 tournaments)
  n=8: max beta_3 = 2 (sampled, 4/5000 = 0.08%)
  n=9: max beta_3 = 2 (sampled, 1/3000 = 0.033%)

This script computes:
1. beta_3 distribution statistics at each n
2. Rate of beta_3 > 0 (nonzero rate)
3. Rate of beta_3 = 2 (max rate)
4. Mean and variance of beta_3

Author: kind-pasteur-S49 (2026-03-09)
"""
import sys
import time
import gc
import numpy as np
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament
from beta3_lean import fast_beta3_lean


def sample_beta3(n, n_trials, seed):
    """Sample beta_3 distribution at given n."""
    rng = np.random.RandomState(seed)
    dist = {}
    vals = []

    t0 = time.time()
    for trial in range(n_trials):
        A = random_tournament(n, rng)
        gc.collect()
        b3 = fast_beta3_lean(A, n)
        dist[b3] = dist.get(b3, 0) + 1
        vals.append(b3)

        if (trial + 1) % 500 == 0:
            elapsed = time.time() - t0
            print(f"    n={n}: {trial+1}/{n_trials}, {elapsed:.1f}s, "
                  f"dist={dict(sorted(dist.items()))}", flush=True)

    elapsed = time.time() - t0
    return dist, vals, elapsed


def main():
    print("=" * 70)
    print("BETA_3 GROWTH RATE ANALYSIS")
    print("=" * 70)

    seed = 77777
    results = {}

    # n=6: fast, can do many
    for n, n_trials in [(6, 5000), (7, 3000), (8, 3000), (9, 1500)]:
        print(f"\n--- n={n}: {n_trials} random tournaments ---")
        dist, vals, elapsed = sample_beta3(n, n_trials, seed)

        mean_b3 = np.mean(vals)
        var_b3 = np.var(vals)
        max_b3 = max(vals)
        nonzero_rate = sum(1 for v in vals if v > 0) / len(vals)
        gt1_rate = sum(1 for v in vals if v > 1) / len(vals) if max_b3 > 1 else 0

        results[n] = {
            'dist': dict(sorted(dist.items())),
            'mean': mean_b3,
            'var': var_b3,
            'max': max_b3,
            'nonzero_rate': nonzero_rate,
            'gt1_rate': gt1_rate,
            'elapsed': elapsed,
        }

        print(f"  n={n}: dist={results[n]['dist']}")
        print(f"    mean={mean_b3:.4f}, var={var_b3:.6f}, max={max_b3}")
        print(f"    P(b3>0)={nonzero_rate:.4f}, P(b3>1)={gt1_rate:.6f}")
        print(f"    Time: {elapsed:.1f}s ({n_trials/elapsed:.1f}/s)")

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY TABLE")
    print(f"{'='*70}")
    print(f"{'n':>3} | {'max b3':>6} | {'P(b3>0)':>8} | {'P(b3>1)':>8} | {'mean':>8} | {'var':>10}")
    print("-" * 60)
    for n in sorted(results.keys()):
        r = results[n]
        print(f"{n:>3} | {r['max']:>6} | {r['nonzero_rate']:>8.4f} | {r['gt1_rate']:>8.6f} | "
              f"{r['mean']:>8.4f} | {r['var']:>10.6f}")

    # Growth analysis
    print(f"\n{'='*70}")
    print("GROWTH ANALYSIS")
    print(f"{'='*70}")
    print("Known exact results:")
    print("  n<=5: max beta_3 = 0 (trivially)")
    print("  n=6:  max beta_3 = 1 (exhaustive)")
    print("  n=7:  max beta_3 = 1 (exhaustive)")
    print("  n=8:  max beta_3 = 2 (5000 samples, 0.08%)")
    print("  n=9:  max beta_3 = 2 (3000 samples, 0.033%)")
    print()
    print("Question: beta_3 <= floor(n/4)?")
    for n in range(5, 13):
        print(f"  n={n}: floor(n/4) = {n//4}")
    print()
    print("beta_3 <= 2 holds through n=9 (sampling).")
    print("floor(n/4) bound: 1,1,1,2,2,2,2,3 for n=5..12")
    print("beta_3 first exceeds 1 at n=8 (floor(8/4)=2). Consistent!")
    print("Prediction: beta_3 <= 2 should hold through n=11.")
    print("beta_3 = 3 first possible at n=12 if floor(n/4) is correct.")

    print(f"\n{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
