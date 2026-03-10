"""
beta3_n9_lean.py — Search for max beta_3 at n=9 using lean implementation

The key question: does beta_3 grow beyond 2 at n=9?

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


def main():
    print("=" * 70)
    print("BETA_3 SEARCH AT n=9 (LEAN)")
    print("=" * 70)

    n = 9
    rng = np.random.RandomState(54321)

    b3_dist = {}
    max_b3 = 0
    total = 0
    b3_gt1_examples = []

    t0 = time.time()
    N_TRIALS = 3000

    for trial in range(N_TRIALS):
        A = random_tournament(n, rng)
        gc.collect()
        b3 = fast_beta3_lean(A, n)

        b3_dist[b3] = b3_dist.get(b3, 0) + 1
        total += 1

        if b3 > max_b3:
            max_b3 = b3
            if b3 > 1:
                scores = sorted([int(sum(A[i])) for i in range(n)])
                print(f"  NEW MAX beta_3={b3} at trial {trial}, scores={scores}")

        if b3 > 1 and len(b3_gt1_examples) < 5:
            scores = sorted([int(sum(A[i])) for i in range(n)])
            b3_gt1_examples.append((trial, b3, scores))

        if (trial + 1) % 200 == 0:
            elapsed = time.time() - t0
            rate = (trial + 1) / elapsed
            gc.collect()
            print(f"  {trial+1}/{N_TRIALS}, {elapsed:.1f}s ({rate:.1f}/s), "
                  f"max_b3={max_b3}, dist={dict(sorted(b3_dist.items()))}", flush=True)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"RESULTS: n={n}, {total} tournaments, {elapsed:.1f}s")
    print(f"  max beta_3 = {max_b3}")
    print(f"  distribution: {dict(sorted(b3_dist.items()))}")
    if b3_gt1_examples:
        print(f"  beta_3 > 1 examples:")
        for trial, b3, scores in b3_gt1_examples:
            print(f"    trial={trial}: b3={b3}, scores={scores}")
    print(f"{'='*70}")
    print("DONE.")


if __name__ == '__main__':
    main()
