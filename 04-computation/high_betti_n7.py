"""Quick check: b5, b6 for b3=1 tournaments at n=7."""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp

n = 7
rng = np.random.RandomState(42)
found = 0
target = 500
b4_dist = Counter()
b5_dist = Counter()
b6_dist = Counter()
chi_dist = Counter()

for trial in range(50000):
    A = random_tournament(n, rng)
    cc = full_chain_complex_modp(A, n, n - 1)
    if cc['bettis'].get(3, 0) != 1:
        continue
    found += 1
    bettis = [cc['bettis'].get(p, 0) for p in range(n)]
    b4_dist[bettis[4]] += 1
    b5_dist[bettis[5]] += 1
    b6_dist[bettis[6]] += 1
    chi = sum((-1)**p * bettis[p] for p in range(n))
    chi_dist[chi] += 1
    if found >= target:
        break

print(f"n={n}, {found} b3=1 tournaments:")
print(f"  b4: {dict(sorted(b4_dist.items()))}")
print(f"  b5: {dict(sorted(b5_dist.items()))}")
print(f"  b6: {dict(sorted(b6_dist.items()))}")
print(f"  chi: {dict(sorted(chi_dist.items()))}")
print(f"\n  => b4=b5=b6=0 always? {all(b4_dist.get(k, 0) == 0 for k in b4_dist if k != 0)}")
