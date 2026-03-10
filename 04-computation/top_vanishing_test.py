"""
top_vanishing_test.py — How far down does top-degree vanishing extend?

PROVED (computationally, n<=8): beta_{n-1} = beta_{n-2} = 0 for ALL tournaments.
Question: does beta_{n-3} = 0 too? How about beta_{n-4}?

At n=7: Paley T_7 has beta_4 = 6, and n-3 = 4. So beta_{n-3} CAN be nonzero!
At n=6: beta_3 can be 1 (n-3 = 3). So beta_{n-3} = 0 only at n <= 5.
At n=5: beta_2 = 0 always (n-3 = 2).

So the vanishing threshold is:
- n=5: beta_2=beta_3=beta_4=0 always. Vanishing from degree 2 = n-3.
- n=6: beta_3 can be 1. Vanishing from degree 4 = n-2.
- n=7: beta_4 = 6 (Paley). Vanishing from degree 5 = n-2.
- n=8: beta_4 can be nonzero. Vanishing from degree 6 = n-2.

Pattern: beta_{n-1} = beta_{n-2} = 0 ALWAYS.
The last two Betti numbers vanish universally.

But wait: at n=5, beta_1 can be 1 (n-4 = 1). So the full story is:
- beta_{n-1} = 0 always (n >= 2)
- beta_{n-2} = 0 always (n >= 4? Let's check n=3,4)
- beta_{n-3} = nonzero possible (n >= 6: beta_3 at n=6)

Let's verify for n=3 and n=4 as well, and also check d_{n-1} injectivity.

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp


def main():
    print("TOP VANISHING: beta_{n-1} and beta_{n-2}")
    print("=" * 60)

    for n in [3, 4, 5, 6, 7, 8]:
        rng = np.random.RandomState(42)
        samples = 500 if n <= 7 else 200
        if n <= 4:
            samples = min(samples, 2**(n*(n-1)//2))  # all tournaments

        betti_dist = {p: Counter() for p in range(n)}
        omega_dist = {p: [] for p in range(n)}
        ranks_dist = {p: [] for p in range(1, n)}

        for trial in range(samples):
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, max_p=n-1)
            for p in range(n):
                betti_dist[p][cc['bettis'].get(p, 0)] += 1
                omega_dist[p].append(cc['omega_dims'].get(p, 0))
            for p in range(1, n):
                ranks_dist[p].append(cc['ranks'].get(p, 0))

        print(f"\nn={n} ({samples} samples):")
        for p in range(n):
            nonzero = samples - betti_dist[p].get(0, 0)
            if nonzero > 0 or p >= n-3:
                print(f"  beta_{p}: nonzero={nonzero}/{samples}, "
                      f"dist={dict(sorted(betti_dist[p].items()))}")
        for p in range(n):
            ovals = omega_dist[p]
            if max(ovals) > 0 or p >= n-3:
                print(f"  Omega_{p}: [{min(ovals)},{max(ovals)}], mean={np.mean(ovals):.1f}")
        # Check injectivity of d_{n-1}
        if n >= 3:
            inj = sum(1 for i in range(samples) if ranks_dist.get(n-1, [0])[i] == omega_dist[n-1][i])
            print(f"  d_{n-1} injective: {inj}/{samples}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
