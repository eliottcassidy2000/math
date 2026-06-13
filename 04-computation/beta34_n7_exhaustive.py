"""
beta34_n7_exhaustive.py — Exhaustive check: do beta_3>0 and beta_4>0 coexist at n=7?

Uses targeted Betti computation (not full chain complex) for speed.
Only computes beta_4 when beta_3 > 0.

Author: kind-pasteur-S48 (2026-03-09)
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '.')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    bits_to_adj, compute_betti_hybrid
)


def main():
    print("=" * 70)
    print("EXHAUSTIVE n=7: beta_3 vs beta_4 COEXISTENCE CHECK")
    print("=" * 70)

    n = 7
    total = 2 ** (n*(n-1)//2)
    print(f"  Total tournaments: {total}")

    b3_count = 0
    b4_count = 0
    b3_and_b4 = 0
    b4_vals_when_b3 = Counter()
    b3_vals_when_b4 = Counter()

    t0 = time.time()
    for bits in range(total):
        A = bits_to_adj(bits, n)
        b3 = compute_betti_hybrid(A, n, 3, max_p=5)

        if b3 > 0:
            b3_count += 1
            b4 = compute_betti_hybrid(A, n, 4, max_p=6)
            b4_vals_when_b3[b4] += 1
            if b4 > 0:
                b3_and_b4 += 1
                score = tuple(sorted([int(sum(A[i])) for i in range(n)]))
                print(f"  COEXISTENCE: bits={bits}, b3={b3}, b4={b4}, score={score}")
        else:
            # Only check beta_4 for a sample to save time
            if bits % 64 == 0:
                b4 = compute_betti_hybrid(A, n, 4, max_p=6)
                if b4 > 0:
                    b4_count += 1
                    b3_vals_when_b4[b3] += 1

        if (bits+1) % 100000 == 0:
            elapsed = time.time() - t0
            rate = (bits+1) / elapsed
            remaining = (total - bits - 1) / rate
            print(f"  {bits+1}/{total} ({(bits+1)*100/total:.1f}%), {rate:.0f}/s, "
                  f"ETA {remaining:.0f}s, b3={b3_count}, coexist={b3_and_b4}",
                  flush=True)

    t1 = time.time()
    print(f"\n  n=7: {total} tournaments in {t1-t0:.1f}s")
    print(f"  beta_3 > 0: {b3_count}")
    print(f"  beta_4 > 0 (sampled 1/64 of b3=0): {b4_count}")
    print(f"  beta_3 AND beta_4 > 0: {b3_and_b4}")
    print(f"\n  beta_4 distribution when beta_3 > 0: {dict(sorted(b4_vals_when_b3.items()))}")

    if b3_and_b4 == 0:
        print("\n  CONFIRMED: beta_3 and beta_4 are MUTUALLY EXCLUSIVE at n=7")
        print("  => chi(T) = 0 for ALL beta_3=1 tournaments at n=7")
    else:
        print(f"\n  COEXISTENCE FOUND: {b3_and_b4} cases")

    print("\nDONE.")


if __name__ == '__main__':
    main()
