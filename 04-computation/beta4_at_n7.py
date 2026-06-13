"""
beta4_at_n7.py — Check β_4 for b3=1 tournaments at n=7

Quick check: is β_4(T) always 0 for tournaments on 7 vertices with β_3=1?
If so, the LES simplifies: δ is injective, H_4^rel = im(δ) = ker(i_*).

Author: opus-2026-03-10-S59
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament, full_chain_complex_modp
)


def main():
    for n in [5, 6, 7, 8]:
        print(f"\n{'='*50}")
        print(f"BETTI NUMBERS AT n={n}")
        print(f"{'='*50}")

        rng = np.random.RandomState(42)
        betti_data = []
        t0 = time.time()
        target = 1000 if n <= 7 else 500

        for trial in range(200000):
            if len(betti_data) >= target:
                break
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, min(n - 1, 6))
            bettis = cc['bettis']
            betti_data.append(bettis)

        elapsed = time.time() - t0
        print(f"  {len(betti_data)} tournaments, {elapsed:.1f}s")

        # Filter b3=1
        b3_1 = [b for b in betti_data if b.get(3, 0) == 1]
        print(f"  β_3=1: {len(b3_1)}/{len(betti_data)}")

        if b3_1:
            for k in range(2, 7):
                vals = Counter(b.get(k, 0) for b in b3_1)
                print(f"  β_{k} (among β_3=1): {dict(sorted(vals.items()))}")

        # All tournaments
        for k in range(2, 7):
            vals = Counter(b.get(k, 0) for b in betti_data)
            nonzero = sum(v for key, v in vals.items() if key > 0)
            print(f"  β_{k} (all): nonzero in {nonzero}/{len(betti_data)}, dist={dict(sorted(vals.items()))}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
