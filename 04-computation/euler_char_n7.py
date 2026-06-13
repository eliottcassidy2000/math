"""
euler_char_n7.py — Is chi = 0 universal at n=7?

chi = sum(-1)^p * beta_p = sum(-1)^p * dim(Omega_p) (alternating Omega dims)

If chi = 0 at n=7 and we know b1*b3=0, b5=0, b6=0:
- When b1=1, b3=0: chi = 1 - 1 + b4 = b4 = 0, so b4 = 0
- When b1=0, b3=1: chi = 1 - 1 + b4 = b4 = 0, so b4 = 0
- When b1=0, b3=0: chi = 1 + b4 = 0, so b4 = -1??? That's impossible!

Wait, I need to be more careful. Let me recompute.

Author: opus-2026-03-09-S56
"""
import sys
import time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    full_chain_complex_modp,
    RANK_PRIME
)


def main():
    print("=" * 70)
    print("EULER CHARACTERISTIC UNIVERSALITY")
    print("=" * 70)

    for n in [5, 6, 7, 8]:
        rng = np.random.RandomState(42)
        chi_dist = Counter()
        betti_profiles = Counter()
        num = 2000 if n <= 7 else 1000

        for trial in range(num):
            A = random_tournament(n, rng)
            cc = full_chain_complex_modp(A, n, n - 1)
            bettis = tuple(cc['bettis'].get(p, 0) for p in range(n))
            chi = sum((-1)**p * bettis[p] for p in range(n))
            chi_dist[chi] += 1
            betti_profiles[bettis] += 1

        print(f"\n  n={n}: chi distribution: {dict(sorted(chi_dist.items()))}")
        for prof, count in betti_profiles.most_common(5):
            chi = sum((-1)**p * prof[p] for p in range(n))
            print(f"    {prof}: {count} times, chi={chi}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
