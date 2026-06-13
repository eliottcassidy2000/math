"""
betti_profile_n8.py — Deep investigation of n=8 Betti profile distribution

Focus on understanding:
1. What score sequences give beta_3>0 vs beta_4>0?
2. What is the exact rate of each profile at n=8?
3. What tournament structure correlates with beta_4>0?
4. Verify the chi constraint: chi = 1 or higher based on profiles

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp


def score_sequence(A, n):
    """Compute sorted score sequence."""
    return tuple(sorted(sum(A[i]) for i in range(n)))


def c3_count(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if {i,j,k} form a 3-cycle
                a, b, c = A[i][j], A[j][k], A[i][k]
                if a*b*(1-c) + (1-a)*(1-b)*c:
                    count += 1
    return count


def main():
    print("n=8 BETTI PROFILE DEEP INVESTIGATION")
    print("=" * 60)

    n = 8
    rng = np.random.RandomState(42)
    samples = 2000
    max_p = n - 1  # = 7 (full complex)

    profiles = Counter()  # betti profile -> count
    score_by_profile = {}  # betti profile -> list of score sequences
    c3_by_profile = {}  # betti profile -> list of c3 counts

    t0 = time.time()
    for trial in range(samples):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, max_p=max_p)

        betti = tuple(cc['bettis'].get(k, 0) for k in range(n))
        chi = sum((-1)**k * betti[k] for k in range(n))

        # Only store interesting profiles
        profile = (betti, chi)
        profiles[profile] += 1

        if profile not in score_by_profile:
            score_by_profile[profile] = Counter()
            c3_by_profile[profile] = Counter()

        sc = score_sequence(A, n)
        score_by_profile[profile][sc] += 1
        c3_by_profile[profile][c3_count(A, n)] += 1

        if (trial + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/{samples} ({elapsed:.0f}s) ...", flush=True)

    elapsed = time.time() - t0
    print(f"\nTotal: {samples} samples, {elapsed:.1f}s")
    print(f"\n{'='*60}")
    print("BETTI PROFILES (sorted by frequency):")
    print(f"{'='*60}")

    for (betti, chi), count in sorted(profiles.items(), key=lambda x: -x[1]):
        rate = count / samples * 100
        b_str = ", ".join(f"b{k}={v}" for k, v in enumerate(betti) if v > 0 or k == 0)
        print(f"\n  Profile: chi={chi}, [{b_str}]")
        print(f"  Count: {count}/{samples} ({rate:.2f}%)")

        # Score sequences
        scores = score_by_profile[(betti, chi)]
        top_scores = sorted(scores.items(), key=lambda x: -x[1])[:5]
        print(f"  Top score sequences: {top_scores}")

        # c3 distribution
        c3s = c3_by_profile[(betti, chi)]
        print(f"  c3 range: {min(c3s.keys())}-{max(c3s.keys())}, most common: {c3s.most_common(3)}")

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    # Count by chi value
    chi_dist = Counter()
    for (betti, chi), count in profiles.items():
        chi_dist[chi] += count
    print("chi distribution:")
    for chi, count in sorted(chi_dist.items()):
        print(f"  chi={chi}: {count}/{samples} ({count/samples*100:.2f}%)")

    # Count by which betas are nonzero
    nonzero = Counter()
    for (betti, chi), count in profiles.items():
        nz = tuple(k for k in range(1, n) if betti[k] > 0)
        nonzero[nz] += count
    print("\nNon-zero beta indices:")
    for nz, count in sorted(nonzero.items(), key=lambda x: -x[1]):
        print(f"  {nz}: {count}/{samples} ({count/samples*100:.2f}%)")


if __name__ == '__main__':
    main()
    print("\nDONE.")
