#!/usr/bin/env python3
"""
ip_hierarchy.py — I(Ω, k) for k = 2, 3, 5 (Cartan determinants).

H = I(Ω, 2). What are I(Ω, 3) and I(Ω, 5)?
Do they have gaps? What values are impossible?

The hierarchy KEY₁=2, KEY₂=3, KEY₃=5:
  I(Ω, 2) = H(T) — Hamiltonian paths
  I(Ω, 3) = ? — "cuboid count"
  I(Ω, 5) = ? — "dodecahedral count"

For a single cycle: I({v}, x) = 1+x.
  I(2) = 3 = KEY₂
  I(3) = 4 = KEY₃-1
  I(5) = 6

For the empty graph: I(∅, x) = 1.
  I(2) = 1 = transitive tournament H

For K_m (clique): I(K_m, x) = 1 + mx.
  I(K_3, 2) = 7
  I(K_3, 3) = 10
  I(K_3, 5) = 16

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def fast_hp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp.get((mask, v), 0)
            if not c or not (mask & (1 << v)):
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get((full, v), 0) for v in range(n))

def get_cycle_data(A, n):
    groups = defaultdict(int)
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            v0 = verts[0]
            for perm in permutations(verts[1:]):
                cycle = (v0,) + perm
                ok = True
                for i in range(length):
                    if A[cycle[i]][cycle[(i+1) % length]] != 1:
                        ok = False
                        break
                if ok:
                    groups[frozenset(verts)] += 1
    return groups

def compute_ip_at_x(groups, x):
    """Compute I(Ω, x) using vertex-set decomposition."""
    vs_list = list(groups.items())
    n_vs = len(vs_list)

    # Build adjacency
    adj = [[False]*n_vs for _ in range(n_vs)]
    for i in range(n_vs):
        for j in range(i+1, n_vs):
            if vs_list[i][0] & vs_list[j][0]:
                adj[i][j] = adj[j][i] = True

    # I(Ω, x) = sum over independent sets of vertex sets,
    # for each: product of d(S_i) * x^k
    total = 0
    for mask in range(2**n_vs):
        verts_in = [i for i in range(n_vs) if mask & (1 << i)]
        independent = True
        for i in range(len(verts_in)):
            for j in range(i+1, len(verts_in)):
                if adj[verts_in[i]][verts_in[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            k = len(verts_in)
            mult = 1
            for idx in verts_in:
                mult *= vs_list[idx][1]
            total += mult * (x**k)

    return total

print("=" * 70)
print("I(Ω, k) HIERARCHY AT k = 2, 3, 5")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# At n=5: compute I(Ω, 2), I(Ω, 3), I(Ω, 5)
# ═══════════════════════════════════════════════════════════════════
print("\n--- n=5: I(Ω, k) for k=2,3,5 ---")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

ip2_vals = Counter()
ip3_vals = Counter()
ip5_vals = Counter()
ip_pairs = defaultdict(set)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_cycle_data(A, n)
    H = fast_hp(A, n)
    i2 = compute_ip_at_x(groups, 2)
    i3 = compute_ip_at_x(groups, 3)
    i5 = compute_ip_at_x(groups, 5)

    assert i2 == H, f"OCF mismatch: H={H}, I(Ω,2)={i2}"

    ip2_vals[i2] += 1
    ip3_vals[i3] += 1
    ip5_vals[i5] += 1
    ip_pairs[(i2, i3)].add(bits)

print(f"  I(Ω,2) values (= H): {sorted(ip2_vals.keys())}")
print(f"  I(Ω,3) values: {sorted(ip3_vals.keys())}")
print(f"  I(Ω,5) values: {sorted(ip5_vals.keys())}")

# Gaps
max_i3 = max(ip3_vals.keys())
i3_gaps = [v for v in range(1, max_i3+1) if v not in ip3_vals and v % 1 == 0]
print(f"\n  I(Ω,3) gaps (missing integers): {i3_gaps[:20]}")

max_i5 = max(ip5_vals.keys())
i5_gaps = [v for v in range(1, max_i5+1) if v not in ip5_vals]
print(f"  I(Ω,5) gaps (first 20): {i5_gaps[:20]}")

# Does I(Ω,3)=10 exist? (analogous to H=7 gap at I(Ω,2))
# I(K₃, 3) = 1+9 = 10
print(f"\n  I(Ω,3) = 10: {'EXISTS' if 10 in ip3_vals else 'MISSING'}")
print(f"  I(Ω,3) = 28: {'EXISTS' if 28 in ip3_vals else 'MISSING'}")  # 28 = 1+3*9 = analogue of 21?

# Show (H, I(Ω,3)) mapping
print(f"\n  (H, I(Ω,3)) pairs:")
for h in sorted(ip2_vals.keys()):
    i3_for_h = set()
    for (h2, i3), bits_set in ip_pairs.items():
        if h2 == h:
            i3_for_h.add(i3)
    if i3_for_h:
        print(f"    H={h:3d}: I(Ω,3) ∈ {sorted(i3_for_h)}")

# ═══════════════════════════════════════════════════════════════════
# At n=6: same analysis (slower)
# ═══════════════════════════════════════════════════════════════════
print("\n--- n=6: I(Ω, k) for k=2,3 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

ip2_vals_6 = set()
ip3_vals_6 = set()

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_cycle_data(A, n)
    i2 = compute_ip_at_x(groups, 2)
    i3 = compute_ip_at_x(groups, 3)

    ip2_vals_6.add(i2)
    ip3_vals_6.add(i3)

print(f"  I(Ω,2) achievable: {sorted(ip2_vals_6)}")
print(f"  I(Ω,3) achievable: {sorted(ip3_vals_6)}")

# Check analogy: I(Ω,3)=10 is the "K₃ evaluation" at x=3
# I(K₃,3) = 1+3·3 = 10
# I(K₁₀,3) = 1+10·3 = 31
# Do these have gaps analogous to H=7 and H=21?

i3_gaps_6 = [v for v in range(1, max(ip3_vals_6)+1) if v not in ip3_vals_6]
print(f"\n  I(Ω,3) gaps at n=6: {i3_gaps_6[:20]}")

# Does I(Ω,3) have permanent gaps?
print(f"\n  I(Ω,3) = 10 at n=6: {'YES' if 10 in ip3_vals_6 else 'NO'}")
print(f"  I(Ω,3) = 31 at n=6: {'YES' if 31 in ip3_vals_6 else 'NO'}")
