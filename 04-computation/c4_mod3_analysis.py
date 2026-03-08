"""
c4_mod3_analysis.py
kind-pasteur-2026-03-07-S37

Analyze c_4(T) = sum_P C(asc(P), 4) to understand why 3|c_4 for n>=7.

Strategy: decompose c_4 by overlap pattern of the 4 ascent positions.
For c_2: c_2 = A_non + (n-2)! * dp(T), both terms 0 mod 3.
For c_4: more complex overlap patterns, but same principle.

c_4 = (1/24) * sum_P asc(P) * (asc(P)-1) * (asc(P)-2) * (asc(P)-3)
    = (1/24) * F^{(4)}(T, 1)
    = sum_P C(asc(P), 4)
    = #{(P, {i1<i2<i3<i4}) : all ij are ascent positions of P in T}

Group by overlap pattern of {i1,i2,i3,i4}.
"""

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
import random
from math import comb, factorial
from itertools import combinations


def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj


def compute_F_dp(adj, n):
    dp = [[[0] * n for _ in range(n)] for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v][0] = 1
    for mask in range(1, 1 << n):
        for last in range(n):
            if not (mask & (1 << last)):
                continue
            for fwd in range(n):
                if dp[mask][last][fwd] == 0:
                    continue
                for nxt in range(n):
                    if mask & (1 << nxt):
                        continue
                    new_mask = mask | (1 << nxt)
                    if adj[last][nxt]:
                        dp[new_mask][nxt][fwd + 1] += dp[mask][last][fwd]
                    else:
                        dp[new_mask][nxt][fwd] += dp[mask][last][fwd]
    full = (1 << n) - 1
    F = [0] * n
    for last in range(n):
        for fwd in range(n):
            F[fwd] += dp[full][last][fwd]
    return F


def v3(n):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % 3 == 0:
        v += 1
        n //= 3
    return v


def overlap_type(positions):
    """Classify 4 positions by their overlap pattern.
    Two positions i,j overlap if |i-j| = 1.
    Returns a canonical descriptor."""
    positions = sorted(positions)
    gaps = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
    overlaps = sum(1 for g in gaps if g == 1)
    return tuple(sorted(gaps))


# For c_4: the 4 positions (i1 < i2 < i3 < i4) in [0, n-2] involve
# vertices P[i1], P[i1+1], P[i2], P[i2+1], P[i3], P[i3+1], P[i4], P[i4+1].
# Depending on overlaps, some of these are the same vertex.

# Key: an ascent at position i means P[i] -> P[i+1] in T.
# For 4 ascent positions, the contribution depends on how many DISTINCT
# vertices are involved and what tournament structure they have.

# At n=7, positions range in [0,5] (6 positions).
# 4 positions from [0,5]: C(6,4) = 15 possible sets.
# But many have overlaps.

# Enumerate position sets and their vertex counts
print("=" * 70)
print("Position set analysis for c_4 at n=7")
print("=" * 70)

n = 7
num_positions = n - 1  # positions 0 to n-2

for pos_set in combinations(range(num_positions), 4):
    vertices = set()
    for p in pos_set:
        vertices.add(p)
        vertices.add(p + 1)
    gaps = tuple(pos_set[i+1] - pos_set[i] for i in range(3))
    print(f"  positions={pos_set}, #vertices={len(vertices)}, gaps={gaps}")


# Now compute c_4 decomposed by position-set type
print("\n" + "=" * 70)
print("c_4 decomposition by position overlap type (n=6, exhaustive)")
print("=" * 70)

n = 6
num_positions = n - 1
m = n * (n - 1) // 2
num_T = 1 << m

# Classify position sets
pos_sets_by_type = {}
for ps in combinations(range(num_positions), 4):
    vertices = set()
    for p in ps:
        vertices.add(p)
        vertices.add(p + 1)
    nv = len(vertices)
    gaps = tuple(ps[i+1] - ps[i] for i in range(3))
    key = (nv, gaps)
    if key not in pos_sets_by_type:
        pos_sets_by_type[key] = []
    pos_sets_by_type[key].append(ps)

print(f"\n  Position set types (n={n}):")
for key in sorted(pos_sets_by_type.keys()):
    nv, gaps = key
    count = len(pos_sets_by_type[key])
    print(f"    #verts={nv}, gaps={gaps}: {count} position set(s)")

# For each tournament, compute the contribution of each type
type_sums = {key: [] for key in pos_sets_by_type}
c4_values = []

for bits in range(num_T):
    adj = tournament_from_bits(n, bits)
    F = compute_F_dp(adj, n)
    c4 = sum(comb(k, 4) * F[k] for k in range(n))
    c4_values.append(c4)

print(f"\n  c_4 at n={n}: min={min(c4_values)}, max={max(c4_values)}")
print(f"  c_4 mod 3: all 0? {all(c % 3 == 0 for c in c4_values)}")
print(f"  v_3(c_4): min={min(v3(c) for c in c4_values if c != 0)}")

# Check c_4 at n=7 (sampled) and decompose
print("\n" + "=" * 70)
print("c_4 structure at n=7 (sampled)")
print("=" * 70)

n = 7
m = n * (n - 1) // 2
num_T = 1 << m
random.seed(42)
num_samples = 5000

c4_mod3 = [0, 0, 0]
c4_mod9 = [0] * 9

for _ in range(num_samples):
    bits = random.randint(0, num_T - 1)
    adj = tournament_from_bits(n, bits)
    F = compute_F_dp(adj, n)
    c4 = sum(comb(k, 4) * F[k] for k in range(n))
    c4_mod3[c4 % 3] += 1
    c4_mod9[c4 % 9] += 1

print(f"  c_4 mod 3 distribution: {c4_mod3}")
print(f"  c_4 mod 9 distribution: {c4_mod9}")
print(f"  Always 0 mod 3: {c4_mod3[1] == 0 and c4_mod3[2] == 0}")


# Key test: is c_4(T) tournament-independent mod 3?
# If so, c_4(T) mod 3 = c_4(transitive) mod 3 for all T.
print("\n" + "=" * 70)
print("Is c_4 tournament-independent (mod 3 or mod 9)?")
print("=" * 70)

for n in [5, 6]:
    m = n * (n - 1) // 2
    num_T = 1 << m
    c4_set = set()
    for bits in range(num_T):
        adj = tournament_from_bits(n, bits)
        F = compute_F_dp(adj, n)
        c4 = sum(comb(k, 4) * F[k] for k in range(n))
        c4_set.add(c4)

    c4_mod3_set = set(c % 3 for c in c4_set)
    c4_mod9_set = set(c % 9 for c in c4_set)
    print(f"  n={n}: {len(c4_set)} distinct c_4 values, "
          f"mod 3 values: {sorted(c4_mod3_set)}, "
          f"mod 9 values: {sorted(c4_mod9_set)}")


# Explore: what does c_4 depend on?
# At n=6: c_4 should NOT be tournament-independent (it's the first free coefficient).
# What tournament invariant determines c_4 mod 3?
print("\n" + "=" * 70)
print("What determines c_4? (n=6)")
print("=" * 70)

n = 6
m = n * (n - 1) // 2
num_T = 1 << m

# Group by score sequence
from collections import defaultdict
score_groups = defaultdict(list)

for bits in range(num_T):
    adj = tournament_from_bits(n, bits)
    scores = tuple(sorted(sum(adj[i]) for i in range(n)))
    F = compute_F_dp(adj, n)
    c4 = sum(comb(k, 4) * F[k] for k in range(n))
    score_groups[scores].append(c4)

print(f"\n  c_4 by score sequence:")
for scores in sorted(score_groups.keys()):
    vals = score_groups[scores]
    mod3_vals = sorted(set(v % 3 for v in vals))
    print(f"    scores={scores}: {len(vals)} T's, c_4 mod 3 in {mod3_vals}, "
          f"c_4 range=[{min(vals)},{max(vals)}]")


print("\nDONE")
