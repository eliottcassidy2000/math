from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import random

def gaps(speeds, t):
    pts = {F(0)}
    for v in speeds:
        pts.add((v * t) % 1)
    p = sorted(pts)
    n = len(p)
    g = []
    for i in range(n):
        nxt = p[(i+1) % n] if i+1 < n else p[0] + 1
        g.append(nxt - p[i])
    return g

def distinct(speeds, t):
    return len(set(gaps(speeds, t)))

print("="*70)
print("PART 2: GENERAL speed sets - number of distinct gap values")
print("="*70)
random.seed(12345)

for m in [4,5,6]:
    overall = Counter()
    maxd = 0
    max_example = None
    n_sets = 600
    for _ in range(n_sets):
        # random distinct positive speeds
        speeds = random.sample(range(1, 40), m)
        # sample exact rational t over various denominators
        for q in [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]:
            for a in range(1, q):
                t = F(a, q)
                d = distinct(speeds, t)
                overall[d] += 1
                if d > maxd:
                    maxd = d
                    max_example = (sorted(speeds), str(t), d)
    npts = m + 1  # at most m+1 points -> at most m+1 gaps
    print(f"  m={m}  (max possible points={npts}, so max possible distinct gaps={npts})")
    print(f"    distinct-gap-count distribution: {dict(sorted(overall.items()))}")
    print(f"    MAX distinct observed: {maxd}   example: {max_example}")
