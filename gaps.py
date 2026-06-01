from fractions import Fraction as F
from itertools import combinations
from collections import Counter
import random

def points(speeds, t):
    """Return sorted set of {0} ∪ {v_i t mod 1} as fractions in [0,1)."""
    pts = {F(0)}
    for v in speeds:
        pts.add((v * t) % 1)
    return sorted(pts)

def gaps(speeds, t):
    """Cyclic consecutive arc-lengths."""
    p = points(speeds, t)
    n = len(p)
    g = []
    for i in range(n):
        nxt = p[(i+1) % n] if i+1 < n else p[0] + 1
        g.append(nxt - p[i])
    return g  # sums to 1

def observer_gaps(speeds, t):
    """The two gaps adjacent to point 0."""
    p = points(speeds, t)
    # 0 is p[0]. right gap = p[1]-p[0]; left gap = (p[0]+1)-p[-1]
    if len(p) == 1:
        return F(0), F(0)
    right = p[1] - p[0]
    left = (p[0] + 1) - p[-1]
    return left, right

def distinct_gap_count(speeds, t):
    return len(set(gaps(speeds, t)))

# ---------- PART 1: AP speeds, Three-Distance Theorem ----------
def part1():
    print("="*70)
    print("PART 1: AP speeds v_i = i, Three-Distance (Steinhaus) Theorem")
    print("="*70)
    for m in [4,5,6,7]:
        speeds = list(range(1, m+1))  # 1..m  -> points k t mod 1 for k=0..m
        counts = Counter()
        violations_3val = 0
        samples = 0
        # exact rational t = a/q
        qs = range(2, 60)
        for q in qs:
            for a in range(1, q):
                t = F(a, q)
                g = gaps(speeds, t)
                vals = sorted(set(g))
                counts[len(vals)] += 1
                samples += 1
                if len(vals) == 3:
                    # largest gap == sum of other two?
                    if vals[2] != vals[0] + vals[1]:
                        violations_3val += 1
        # also random floats rounded
        print(f"  m={m} (points k t, k=0..{m}): {samples} rational t sampled")
        print(f"    distinct-gap-count distribution: {dict(sorted(counts.items()))}")
        print(f"    3-value cases violating (largest = sum of other two): {violations_3val}")

part1()
