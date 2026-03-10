# ⚠️ WARNING: This script assumes beta_3 ≤ 1 for all tournaments, which is
# FALSE at n=8 (beta_3=2 found in 0.08% of tournaments). See MISTAKE-018.
# Results are valid only for n ≤ 7.

"""
beta3_large_n_hybrid.py — Verify beta_3 <= 1 at n=8 and n=9 using hybrid method

Uses compute_betti_hybrid from tournament_utils.py (numpy int64 + mod-p Gauss)
to compute beta_3 for random tournaments at n=8 (1000 samples) and n=9 (200 samples).

Author: kind-pasteur (2026-03-09)
"""
import sys
import time
import numpy as np
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, compute_betti_hybrid, bits_to_adj, adj_to_bits

def run_experiment(n, num_samples, seed, progress_interval):
    """Run beta_3 verification for num_samples random tournaments of size n."""
    print(f"\n{'='*60}")
    print(f"  n={n}: Testing {num_samples} random tournaments (seed={seed})")
    print(f"{'='*60}")

    rng = np.random.RandomState(seed)
    beta3_values = []
    timings = []
    violations = []  # tournaments with beta_3 > 1

    t_total_start = time.time()

    for trial in range(num_samples):
        A = random_tournament(n, rng)
        bits = adj_to_bits(A, n)

        t_start = time.time()
        b3 = compute_betti_hybrid(A, n, 3, max_p=5)
        t_elapsed = time.time() - t_start

        beta3_values.append(b3)
        timings.append(t_elapsed)

        if b3 > 1:
            violations.append((trial, bits, b3))
            print(f"  *** VIOLATION at trial {trial}: beta_3 = {b3}, bits = {bits}")

        if (trial + 1) % progress_interval == 0:
            avg_t = sum(timings[-progress_interval:]) / progress_interval
            print(f"  Progress: {trial+1}/{num_samples} done | "
                  f"last {progress_interval} avg {avg_t:.3f}s/tour | "
                  f"max beta_3 so far = {max(beta3_values)}")

    t_total = time.time() - t_total_start
    dist = Counter(beta3_values)

    print(f"\n--- Results for n={n} ---")
    print(f"  Samples:       {num_samples}")
    print(f"  Total time:    {t_total:.1f}s ({t_total/num_samples:.3f}s avg per tournament)")
    print(f"  Max beta_3:    {max(beta3_values)}")
    print(f"  Min beta_3:    {min(beta3_values)}")
    print(f"  Distribution:  {dict(sorted(dist.items()))}")
    if violations:
        print(f"  VIOLATIONS (beta_3 > 1): {len(violations)}")
        for trial, bits, b3 in violations:
            print(f"    trial={trial}, bits={bits}, beta_3={b3}")
    else:
        print(f"  No violations found (beta_3 <= 1 for all samples)")
    print(f"  Timing range:  [{min(timings):.3f}s, {max(timings):.3f}s]")
    print()

    return beta3_values, violations


def main():
    print("beta3_large_n_hybrid.py")
    print(f"Using compute_betti_hybrid (numpy int64 + mod-p Gauss)")
    print(f"Verifying beta_3 <= 1 at n=8 and n=9")
    print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # n=8: 1000 random tournaments
    vals8, viol8 = run_experiment(n=8, num_samples=1000, seed=42, progress_interval=100)

    # n=9: 200 random tournaments
    vals9, viol9 = run_experiment(n=9, num_samples=200, seed=42, progress_interval=20)

    # Summary
    print("=" * 60)
    print("  SUMMARY")
    print("=" * 60)
    print(f"  n=8: max beta_3 = {max(vals8)}, violations = {len(viol8)}")
    print(f"  n=9: max beta_3 = {max(vals9)}, violations = {len(viol9)}")
    if not viol8 and not viol9:
        print(f"  CONCLUSION: beta_3 <= 1 verified for all tested tournaments.")
    else:
        print(f"  WARNING: Violations found! See details above.")
    print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
