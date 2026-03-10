"""chi distribution for ALL tournaments at n=7."""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp

for n in [5, 6, 7, 8]:
    rng = np.random.RandomState(42)
    chi_dist = Counter()
    chi_betti = Counter()
    num = 5000 if n <= 7 else 2000

    for _ in range(num):
        A = random_tournament(n, rng)
        cc = full_chain_complex_modp(A, n, n - 1)
        bettis = [cc['bettis'].get(p, 0) for p in range(n)]
        chi = sum((-1)**p * bettis[p] for p in range(n))
        chi_dist[chi] += 1

    print(f"n={n}: chi distribution: {dict(sorted(chi_dist.items()))}")
    print(f"  chi <= 1: {sum(v for k, v in chi_dist.items() if k <= 1)}/{num}")
