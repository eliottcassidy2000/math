#!/usr/bin/env python3
"""
delta_h_gap.py — opus-2026-03-14-S71f

ΔH distribution under single arc flip at n=5 skips ±10.
Is this connected to the H=21 sum gap (α₁+2α₂=10)?

If ΔH = 2·Δα₁ + 4·Δα₂, and ΔH=10 never occurs, then
  Δα₁ + 2·Δα₂ = 5 never occurs under single arc flip.

Check: what Δα₁ values occur under single arc flip?
"""

import numpy as np
from itertools import permutations, combinations
from collections import Counter

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def hp_dp(A, n):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if A[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

# n=5: exhaustive
n = 5
num_arcs = n*(n-1)//2

H_val = {}
for mask in range(1 << num_arcs):
    A = make_tournament(mask, n)
    H_val[mask] = hp_dp(A, n)

delta_dist = Counter()
for mask in range(1 << num_arcs):
    for arc in range(num_arcs):
        d = H_val[mask ^ (1 << arc)] - H_val[mask]
        delta_dist[d] += 1

print("n=5 ΔH distribution:")
for d in sorted(delta_dist.keys()):
    print(f"  ΔH = {d:4d}: {delta_dist[d]:6d}")

# Now n=4
print()
n = 4
num_arcs = 6
H4 = {}
for mask in range(1 << num_arcs):
    A = make_tournament(mask, n)
    H4[mask] = hp_dp(A, n)

dd4 = Counter()
for mask in range(1 << num_arcs):
    for arc in range(num_arcs):
        d = H4[mask ^ (1 << arc)] - H4[mask]
        dd4[d] += 1

print("n=4 ΔH distribution:")
for d in sorted(dd4.keys()):
    print(f"  ΔH = {d:4d}: {dd4[d]:6d}")

# n=6: exhaustive
print()
n = 6
num_arcs = 15
H6 = {}
for mask in range(1 << num_arcs):
    A = make_tournament(mask, n)
    H6[mask] = hp_dp(A, n)

dd6 = Counter()
for mask in range(1 << num_arcs):
    for arc in range(num_arcs):
        d = H6[mask ^ (1 << arc)] - H6[mask]
        dd6[d] += 1

print("n=6 ΔH distribution:")
for d in sorted(dd6.keys()):
    print(f"  ΔH = {d:4d}: {dd6[d]:6d}")

# Analysis: what values are skipped?
print(f"\n{'='*60}")
print("Gap Analysis")
print(f"{'='*60}")

for nn, dd in [(4, dd4), (5, delta_dist), (6, dd6)]:
    vals = sorted(dd.keys())
    max_val = max(abs(v) for v in vals)
    expected = set(range(-max_val, max_val+1, 2))
    actual = set(vals)
    gaps = sorted(expected - actual)
    print(f"\n  n={nn}: ΔH range [{min(vals)}, {max(vals)}]")
    print(f"       All even? {all(v % 2 == 0 for v in vals)}")
    print(f"       Skipped even values: {gaps if gaps else 'NONE'}")
    print(f"       |ΔH|/2 values: {sorted(set(abs(v)//2 for v in vals))}")

# Key question: is |ΔH|/2 = |Δα₁| always?
# At n=5 (α₂=0): ΔH/2 = Δα₁, so |ΔH|/2 ∈ {0,1,2,3,4,5,6}
# But |ΔH|/2 values are {0,1,2,3,4,6} — skipping 5!

# Is |Δα₁| = 5 impossible under single arc flip?
print(f"\n{'='*60}")
print("Is |Δα₁| = 5 impossible under single arc flip at n=5?")
print(f"{'='*60}")

# At n=5, α₁ = number of directed 3-cycles + directed 5-cycles
# Flipping one arc (i,j)→(j,i):
# - 3-cycles affected: those containing the arc (i,j)
#   Each such cycle has vertices {i, j, k} for some k
#   The arc (i,j) participates in n-2 = 3 vertex triples
# - Changes at most 3 directed 3-cycles (each triple either gains or loses one)
# - Plus changes in 5-cycles

# But the directed cycle count:
# At n=5, a 3-cycle on {i,j,k}: contributes 1 directed cycle
# Flipping arc (i,j) can create/destroy cycles on {i,j,k} for k ∈ {0,..,4}\{i,j}
# That's 3 possible affected triples
# Each triple either:
#   (a) had a cycle through (i,j) → loses 1 cycle, may gain 1 (opposite direction)
#   (b) had a cycle through (j,i) → loses 1 (opposite), may gain 1
#   (c) had no cycle through either → may gain 1

# Actually for 3-cycles: flipping (i,j) creates or destroys EXACTLY the cycle
# through (i,j) on each triple {i,j,k}. The net change in t₃ is:
# Δt₃ ∈ {-3, -1, 1, 3} (bounded by 3)

# For 5-cycles at n=5: the arc (i,j) appears in some 5-cycles
# At n=5, the only possible 5-cycle set is {0,1,2,3,4} (all vertices)
# So flipping (i,j) can change the directed 5-cycle count on {0,1,2,3,4}

# Directed 5-cycles change: at most ±3 (multiplicity range is [0,3])
# So |Δα₁| ≤ |Δt₃| + |Δd₅| ≤ 3 + 3 = 6

# But |Δα₁| = 5 requires |Δt₃| + |Δd₅| = 5 with same sign
# This requires specific coordination between 3-cycle and 5-cycle changes

# The fact that |ΔH|/2 = 5 never occurs means this coordination is impossible!

print("\nFlipping arc (i,j) at n=5 affects:")
print("  - 3-cycles on triples {i,j,k}: at most 3 triples")
print("  - 5-cycles on {0,1,2,3,4}: the single full set")
print("  |Δ(3-cycles)| ≤ 3, |Δ(5-cycles)| ≤ 3")
print("  |Δα₁| = |Δ(3-cyc) + Δ(5-cyc)| ≤ 6")
print("  But |Δα₁| = 5 never occurs (ΔH = ±10 gap)")
print()
print("This suggests a PARITY CONSTRAINT on Δ(3-cyc) ± Δ(5-cyc).")

# Verify: compute Δ(3-cycles) and Δ(5-cycles) separately
from itertools import permutations

def count_3cyc_5cyc(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1
    d5 = 0
    for verts in combinations(range(n), 5):
        for p in permutations(verts[1:]):
            order = [verts[0]] + list(p)
            ok = True
            for idx in range(5):
                if A[order[idx]][order[(idx+1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                d5 += 1
                # Each VERTEX SET can have 1-3 directed cycles
                # but we want total count, so don't break
    return t3, d5

n = 5
num_arcs = 10
delta_pairs = Counter()  # (Δt3, Δd5)
for mask in range(1 << num_arcs):
    A0 = make_tournament(mask, n)
    t3_0, d5_0 = count_3cyc_5cyc(A0, n)
    for arc in range(num_arcs):
        A1 = make_tournament(mask ^ (1 << arc), n)
        t3_1, d5_1 = count_3cyc_5cyc(A1, n)
        dt3 = t3_1 - t3_0
        dd5 = d5_1 - d5_0
        delta_pairs[(dt3, dd5)] += 1

print("\nJoint (Δt₃, Δd₅) distribution at n=5:")
print(f"{'Δt₃':>5s} {'Δd₅':>5s} {'Count':>6s} {'Δα₁':>5s} {'ΔH':>5s}")
for (dt3, dd5), count in sorted(delta_pairs.items()):
    da1 = dt3 + dd5  # since α₁ = t₃ + d₅ (vertex sets)
    # Wait, at n=5 there's only 1 possible 5-vertex set, so d5 = # directed 5-cycles on ALL vertices
    # These are SEPARATE vertices in Ω
    dh = 2 * (dt3 + dd5)  # when α₂=0
    print(f"{dt3:>5d} {dd5:>5d} {count:>6d} {da1:>5d} {dh:>5d}")

# Check: is Δt₃ + Δd₅ = 5 ever achieved?
has_5 = any(abs(dt3 + dd5) == 5 for (dt3, dd5) in delta_pairs)
print(f"\n|Δt₃ + Δd₅| = 5 achieved? {has_5}")

# What about Δt₃ + Δd₅ parity?
print("\nParity of Δα₁ = Δt₃ + Δd₅:")
parity = Counter()
for (dt3, dd5), count in delta_pairs.items():
    parity[(dt3 + dd5) % 2] += count
print(f"  Even: {parity[0]:6d}")
print(f"  Odd:  {parity[1]:6d}")
