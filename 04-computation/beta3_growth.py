"""
beta3_growth.py — How does max(beta_3) grow with n?

Known: beta_3 <= 1 at n <= 7 (exhaustive n<=7).
       beta_3 = 2 at n=8 (0.08%), n=9 (0.05%).

Question: Does beta_3 grow with n? If so, how?
- beta_3 <= floor(n/4)?
- beta_3 <= n - 6?
- Something else?

Also: what Betti profiles occur at each n?

Author: opus-2026-03-09-S56
"""
import sys
import time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    full_chain_complex_modp,
    RANK_PRIME
)


def survey_bettis(n, num_samples, seed=42):
    """Survey Betti number distributions at given n."""
    rng = np.random.RandomState(seed)
    t0 = time.time()

    betti_profiles = Counter()
    max_b3 = 0
    b3_dist = Counter()

    for i in range(num_samples):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        bettis = tuple(cc['bettis'].get(p, 0) for p in range(n))
        betti_profiles[bettis] += 1
        b3 = bettis[3] if len(bettis) > 3 else 0
        b3_dist[b3] += 1
        if b3 > max_b3:
            max_b3 = b3
            elapsed = time.time() - t0
            print(f"    n={n}: NEW max beta_3 = {b3} at sample {i+1}, {elapsed:.1f}s")
            scores = sorted(sum(A[i]) for i in range(n))
            print(f"      bettis={bettis}, scores={tuple(scores)}")

        if (i + 1) % 1000 == 0:
            elapsed = time.time() - t0
            print(f"    n={n}: {i+1}/{num_samples}, max_b3={max_b3}, {elapsed:.1f}s")

    elapsed = time.time() - t0
    return {
        'n': n,
        'samples': num_samples,
        'max_b3': max_b3,
        'b3_dist': dict(sorted(b3_dist.items())),
        'num_profiles': len(betti_profiles),
        'top_profiles': betti_profiles.most_common(10),
        'elapsed': elapsed,
    }


def main():
    print("=" * 70)
    print("BETA_3 GROWTH SURVEY")
    print("=" * 70)

    # Survey each n with appropriate sample sizes
    configs = [
        (7, 5000),    # exhaustive is 2M, but beta_3<=1 is known
        (8, 5000),    # beta_3=2 at 0.08%
        (9, 3000),    # beta_3=2 at 0.05%
        (10, 2000),   # new territory
        (11, 1000),   # slower, fewer samples
    ]

    all_results = []

    for n, num_samples in configs:
        print(f"\n--- n={n}, {num_samples} samples ---")
        r = survey_bettis(n, num_samples)
        all_results.append(r)

        print(f"  max(beta_3) = {r['max_b3']}")
        print(f"  beta_3 distribution: {r['b3_dist']}")
        print(f"  {r['num_profiles']} distinct Betti profiles")
        print(f"  Top 5 profiles:")
        for prof, count in r['top_profiles'][:5]:
            print(f"    {prof}: {count} ({100*count/r['samples']:.1f}%)")
        print(f"  Time: {r['elapsed']:.1f}s")

    print(f"\n{'='*70}")
    print("SUMMARY: max(beta_3) by n")
    print("=" * 70)
    for r in all_results:
        print(f"  n={r['n']}: max(beta_3) = {r['max_b3']} "
              f"(from {r['samples']} samples)")

    # Check growth hypotheses
    print(f"\nGrowth hypotheses:")
    for r in all_results:
        n = r['n']
        m = r['max_b3']
        print(f"  n={n}: max_b3={m}, floor(n/4)={n//4}, n-6={n-6}, "
              f"floor((n-5)/2)={(n-5)//2}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
