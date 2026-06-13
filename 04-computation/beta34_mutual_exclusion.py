"""
beta34_mutual_exclusion.py — Does beta_3 > 0 force beta_4 = 0?

If yes, chi(T) = 0 for all beta_3=1 tournaments, which gives the
chi_rel dichotomy algebraically.

Key test: check n=8 where beta_4 > 0 is possible (rare).
At n=7, Paley T_7 has beta_4=6 but beta_3=0.

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp, compute_betti_hybrid
)


def main():
    print("=" * 70)
    print("BETA_3 vs BETA_4 MUTUAL EXCLUSION TEST")
    print("=" * 70)

    # Part 1: n=8 — search for beta_3=1 tournaments and check beta_4
    print("\n--- Part 1: n=8 targeted search ---")
    n = 8
    rng = np.random.RandomState(42)

    b3_count = 0
    b4_count = 0
    b3_and_b4 = 0
    b3_chi = Counter()
    b3_betti = Counter()
    total = 0

    t0 = time.time()
    for trial in range(2000):
        A = random_tournament(n, rng)
        # First compute beta_3 (fast)
        b3 = compute_betti_hybrid(A, n, 3, max_p=5)
        total += 1

        if b3 > 0:
            b3_count += 1
            # Now compute full chain complex including beta_4
            res = full_chain_complex_modp(A, n, max_p=7)
            b = res['bettis']
            b4 = b.get(4, 0)
            if b4 > 0:
                b3_and_b4 += 1
                print(f"  COEXISTENCE: trial={trial}, bettis={tuple(b.get(p,0) for p in range(8))}")

            chi = sum((-1)**p * b.get(p, 0) for p in range(8))
            b3_chi[chi] += 1
            bv = tuple(b.get(p, 0) for p in range(8))
            b3_betti[bv] += 1

        # Also search for beta_4 > 0 separately
        if b3 == 0:
            b4 = compute_betti_hybrid(A, n, 4, max_p=6)
            if b4 > 0:
                b4_count += 1

        if (trial+1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  {trial+1}/2000, {elapsed:.1f}s, b3={b3_count}, b4={b4_count}, coexist={b3_and_b4}",
                  flush=True)

    t1 = time.time()
    print(f"\n  n=8: {total} tournaments in {t1-t0:.1f}s")
    print(f"  beta_3 > 0: {b3_count}")
    print(f"  beta_4 > 0 (with b3=0): {b4_count}")
    print(f"  beta_3 > 0 AND beta_4 > 0: {b3_and_b4}")
    print(f"\n  chi(T) for beta_3=1: {dict(sorted(b3_chi.items()))}")
    print(f"\n  Full Betti vectors for beta_3=1:")
    for bv, cnt in sorted(b3_betti.items()):
        print(f"    {bv}: {cnt}")

    # Part 2: n=8 targeted search for beta_4 > 0 (Paley-like)
    print("\n--- Part 2: n=8 beta_4 search (score-biased) ---")
    rng2 = np.random.RandomState(123)
    b4_found = 0
    b4_betti = []

    for trial in range(500):
        A = random_tournament(n, rng2)
        b4 = compute_betti_hybrid(A, n, 4, max_p=6)
        if b4 > 0:
            b4_found += 1
            res = full_chain_complex_modp(A, n, max_p=7)
            b = res['bettis']
            bv = tuple(b.get(p, 0) for p in range(8))
            score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
            b4_betti.append((bv, score))
            if b4_found <= 10:
                print(f"  beta_4={b4} at trial={trial}, bettis={bv}, score={score}")

    print(f"\n  beta_4 > 0 found: {b4_found}/500")
    if b4_betti:
        b4_b3_coexist = sum(1 for bv, _ in b4_betti if bv[3] > 0)
        print(f"  Among these, beta_3 > 0: {b4_b3_coexist}")

    # Part 3: Summary
    print("\n--- Part 3: Summary ---")
    if b3_and_b4 == 0:
        print("  beta_3 and beta_4 appear MUTUALLY EXCLUSIVE at n=8 (0 coexistence in 2000+500 samples)")
        if all(chi == 0 for chi in b3_chi.keys()):
            print("  chi(T) = 0 for ALL beta_3=1 tournaments at n=8")
            print("  => chi_rel dichotomy follows:")
            print("     Good vertex: chi(T\\v) = 1, chi_rel = -1, H_3^rel = 1")
            print("     Bad vertex:  chi(T\\v) = 0, chi_rel = 0, H_3^rel = 0")
        else:
            print(f"  BUT chi(T) not always 0: {dict(b3_chi)}")
    else:
        print(f"  COEXISTENCE FOUND: {b3_and_b4} cases at n=8")

    print("\nDONE.")


if __name__ == '__main__':
    main()
