#!/usr/bin/env python3
"""
alpha2_gaps.py — opus-2026-03-14-S71d

Investigating GAPS in the (α₁, α₂) lattice.

At n=6, dc3=8: α₂ ∈ {1,2,4} but NOT 3!
This is a "gap" in the α₂ values. Why?

Also: connect to the forbidden H values discovered by S63.
H=21 is impossible — is this because the (α₁,α₂) decomposition is impossible?
"""

import sys
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                pm = mask ^ (1 << v)
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def count_ham_cycles(A, n):
    if n < 3: return 0
    full_mask = (1 << n) - 1
    dp = {(1 << 0, 0): 1}
    for ms in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != ms: continue
            if not (mask & 1): continue
            for v in range(n):
                if not (mask & (1 << v)): continue
                if v == 0 and ms < n: continue
                pm = mask ^ (1 << v)
                if not (pm & 1): continue
                t = 0
                for u in range(n):
                    if (pm & (1 << u)) and A[u][v]:
                        t += dp.get((pm, u), 0)
                if t:
                    dp[(mask, v)] = t
    total = 0
    for v in range(1, n):
        if A[v][0] and (full_mask, v) in dp:
            total += dp[(full_mask, v)]
    return total

def count_directed_k_cycles(A, n, k):
    if k > n: return 0
    total = 0
    for combo in combinations(range(n), k):
        verts = list(combo)
        sub = np.zeros((k, k), dtype=int)
        for i in range(k):
            for j in range(k):
                sub[i][j] = A[verts[i]][verts[j]]
        total += count_ham_cycles(sub, k)
    return total

# ======================================================================
# PART 1: The full (α₁, α₂) lattice at n=6 (exhaustive)
# ======================================================================
print("=" * 70)
print("PART 1: FULL (α₁, α₂) LATTICE AT n=6")
print("=" * 70)

n = 6
tb = n*(n-1)//2

lattice = defaultdict(int)
dc_data = defaultdict(lambda: defaultdict(int))

for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    lattice[(a1, a2)] += 1
    dc_data[(dc3, dc5)][a2] += 1

print(f"\n  Achievable (α₁, α₂) points:")
for a1 in sorted(set(k[0] for k in lattice)):
    a2_vals = sorted(set(k[1] for k in lattice if k[0] == a1))
    counts = [lattice[(a1, a2)] for a2 in a2_vals]
    print(f"    α₁={a1:3d}: α₂ ∈ {a2_vals} (counts: {counts})")

# Find gaps in α₂ for each α₁
print(f"\n  GAPS in α₂ for each α₁:")
for a1 in sorted(set(k[0] for k in lattice)):
    a2_vals = sorted(set(k[1] for k in lattice if k[0] == a1))
    if len(a2_vals) > 1:
        full_range = list(range(min(a2_vals), max(a2_vals)+1))
        missing = [x for x in full_range if x not in a2_vals]
        if missing:
            print(f"    α₁={a1:3d}: α₂ range [{min(a2_vals)},{max(a2_vals)}], MISSING: {missing}")

# Full grid
a1_vals = sorted(set(k[0] for k in lattice))
a2_vals = sorted(set(k[1] for k in lattice))
print(f"\n  Achievable lattice grid (count):")
header = "  α₂\\α₁  " + "".join(f"{a1:5d}" for a1 in a1_vals)
print(header)
for a2 in a2_vals:
    row = f"  {a2:3d}     " + "".join(
        f"{lattice.get((a1,a2),0):5d}" if lattice.get((a1,a2),0) > 0 else "    ."
        for a1 in a1_vals)
    print(row)

# ======================================================================
# PART 2: dc3=8 at n=6 — why α₂=3 is missing
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: dc3=8 AT n=6 — THE MISSING α₂=3")
print("=" * 70)

# At n=6: dc3=8 means 8 out of C(6,3)=20 triples are cyclic.
# The complement: 12 transitive.
# α₂ = #{complementary pairs (T, T^c) where both cyclic}

# There are 10 complementary pairs. With 8 cyclic triples:
# Each pair (T, T^c): T has 1 complement. Need both cyclic.
# 8 cyclic triples. Each cyclic triple T has complement T^c.
# If T is cyclic, T^c might or might not be.

# With 8 cyclic: 12 transitive. So 12 triples are NOT cyclic.
# Among 10 complementary pairs, some have T cyclic and T^c transitive.
# We need pairs where BOTH are cyclic.

# How many transitive triples can be "complementary to cyclic"?
# If t triples are transitive, and they account for p complementary pairs
# where one is transitive, then α₂ = 10 - #{pairs with at least one trans.}

print(f"  With dc3=8 at n=6:")
print(f"    8 cyclic, 12 transitive out of C(6,3)=20 triples")
print(f"    10 complementary pairs")

# Let's enumerate all tournaments with dc3=8 at n=6
dc8_data = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    if dc3 != 8:
        continue
    H = count_ham_paths(A, n)
    dc5 = count_directed_k_cycles(A, n, 5)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4

    # Check complementary pairs
    all_triples = list(combinations(range(6), 3))
    cyclic_set = set()
    for combo in all_triples:
        a, b, c = combo
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            cyclic_set.add(frozenset(combo))

    comp_pairs = []
    all_v = frozenset(range(6))
    for combo in combinations(range(6), 3):
        t1 = frozenset(combo)
        t2 = all_v - t1
        if t1 < t2:
            c1 = t1 in cyclic_set
            c2 = t2 in cyclic_set
            comp_pairs.append((c1, c2))

    both_cyclic = sum(1 for c1, c2 in comp_pairs if c1 and c2)
    one_trans = sum(1 for c1, c2 in comp_pairs if not (c1 and c2))
    neither = sum(1 for c1, c2 in comp_pairs if not c1 and not c2)

    dc8_data.append({
        'bits': bits, 'a2': a2, 'both_cyclic': both_cyclic,
        'one_trans': one_trans, 'neither': neither,
        'dc5': dc5, 'H': H
    })

# Group by α₂
a2_groups = defaultdict(list)
for d in dc8_data:
    a2_groups[d['a2']].append(d)

print(f"\n  dc3=8 tournaments: {len(dc8_data)} total")
for a2 in sorted(a2_groups.keys()):
    g = a2_groups[a2]
    ex = g[0]
    print(f"    α₂={a2}: {len(g)} tournaments. Example: both_cyclic={ex['both_cyclic']}, "
          f"one_trans={ex['one_trans']}, neither={ex['neither']}, dc5={ex['dc5']}")

print(f"\n  So α₂=3 is missing because... let's check the complementary pair constraint.")
print(f"  With 8 cyclic triples among 20:")
print("    10 complementary pairs. Let X = number of (cyclic T with cyclic T^c)")
print(f"    Then 2X of the 8 cyclic triples have cyclic complements.")
print(f"    So 8-2X cyclic triples have transitive complements.")
print(f"    And X pairs have both cyclic, 8-2X pairs have one cyclic+one trans,")
print(f"    10-X-(8-2X) = 10-8+X = 2+X pairs have both transitive.")
print(f"    But we need 2+X ≥ 0, so X ≥ -2 (always true).")
print(f"    Also need 8-2X ≥ 0, so X ≤ 4.")
print(f"    And need 2+X ≤ 12/2=6 (since 12 transitive triples form at most 6 pairs),")
print(f"    so X ≤ 4.")
print(f"    X = α₂ can be 0,1,2,3,4 a priori.")
print(f"    But the actual achievable values are {sorted(a2_groups.keys())}!")
print(f"    α₂=0 missing: with 8 cyclic, it's impossible for NONE to be complementary.")
print(f"    α₂=3 missing: some geometric constraint prevents exactly 3 pairs.")

# ======================================================================
# PART 3: Score (2,2,2,3,3,3) analysis
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: REGULAR-LIKE SCORES AT n=6")
print("=" * 70)

# Score (2,2,2,3,3,3) has dc3=8 always. What's the score structure?
score_233 = []
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    if scores != (2, 2, 2, 3, 3, 3):
        continue
    H = count_ham_paths(A, n)
    dc3 = count_directed_k_cycles(A, n, 3)
    dc5 = count_directed_k_cycles(A, n, 5)
    a1 = dc3 + dc5
    a2 = (H - 1 - 2*a1) // 4
    score_233.append({'bits': bits, 'dc3': dc3, 'dc5': dc5, 'a1': a1, 'a2': a2, 'H': H})

print(f"\n  Score (2,2,2,3,3,3): {len(score_233)} tournaments")
a2_dist = Counter(d['a2'] for d in score_233)
H_dist = Counter(d['H'] for d in score_233)
print(f"    α₂ distribution: {dict(sorted(a2_dist.items()))}")
print(f"    H distribution: {dict(sorted(H_dist.items()))}")

# ======================================================================
# PART 4: Forbidden H values and the lattice
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: FORBIDDEN H VALUES")
print("=" * 70)

# Collect all achievable H at n=6
all_H = set()
for bits in range(1 << tb):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)
    all_H.add(H)

max_H = max(all_H)
print(f"\n  n=6: achievable H = {sorted(all_H)}")
print(f"  Missing odd H in [1,{max_H}]: {sorted(set(range(1, max_H+1, 2)) - all_H)}")

# At n=5:
n5 = 5
all_H5 = set()
for bits in range(1 << (5*4//2)):
    A = bits_to_adj(bits, n5)
    H = count_ham_paths(A, n5)
    all_H5.add(H)
print(f"  n=5: achievable H = {sorted(all_H5)}")
print(f"  Missing odd H in [1,{max(all_H5)}]: {sorted(set(range(1, max(all_H5)+1, 2)) - all_H5)}")

# H = 1 + 2α₁ + 4α₂ with α₁, α₂ ≥ 0
# H is always ≡ 1 (mod 2)
# H ≡ 1 or 3 (mod 4) depending on α₁ parity
# H = 7: need α₁=3, α₂=0 (7=1+6). Is α₁=3 achievable at n=5? n=6?
n = 6
# Check: can α₁=3 at n=6?
count_a1_3 = sum(1 for bits in range(1 << tb)
    if True  # placeholder
)
# Use the lattice we already computed
a1_3_exists = any(k[0] == 3 for k in lattice)
print(f"\n  α₁=3 achievable at n=6? {a1_3_exists}")
if not a1_3_exists:
    print(f"  → H=7 (=1+2*3+4*0) is IMPOSSIBLE at n=6!")
    # What about decompositions with α₂>0?
    # H=7: 1+2α₁+4α₂=7, so 2α₁+4α₂=6, α₁+2α₂=3.
    # Possible: (α₁,α₂)=(3,0), (1,1).
    # Is (1,1) achievable?
    a1_1_a2_1 = (1,1) in lattice
    print(f"  → (α₁,α₂)=(1,1) achievable at n=6? {a1_1_a2_1}")
    if not a1_1_a2_1:
        print(f"  → ALL decompositions of H=7 are impossible at n=6!")

# ======================================================================
# PART 5: The achievable (α₁, α₂) as a polytope
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: (α₁, α₂) AS A POLYTOPE")
print("=" * 70)

# What are the EXACT constraints on (α₁, α₂) at n=6?
# α₁ = dc3 + dc5, α₂ = dp33
# dc3 = C(6,3) - Σ C(s_i,2) ∈ {0,1,...,8} at n=6
# dc5 ∈ {0,...,many}
# dp33 ∈ {0,...,4} at n=6 (constrained by complementary structure)

# Let's map (dc3, dc5) → achievable α₂
print(f"\n  (dc3, dc5) → achievable α₂ at n=6:")
for (dc3, dc5) in sorted(dc_data.keys()):
    a2_vals = sorted(dc_data[(dc3, dc5)].keys())
    counts = [dc_data[(dc3, dc5)][a2] for a2 in a2_vals]
    print(f"    dc3={dc3:2d}, dc5={dc5:3d}: α₂ ∈ {a2_vals}")

# The boundary: max α₂ as function of α₁
print(f"\n  max(α₂) as function of α₁:")
for a1 in sorted(set(k[0] for k in lattice)):
    max_a2 = max(k[1] for k in lattice if k[0] == a1)
    count = sum(lattice[(a1, a2)] for a2 in range(max_a2+1) if (a1, a2) in lattice)
    print(f"    α₁={a1:3d}: max(α₂)={max_a2}, H_max = {1+2*a1+4*max_a2}")

# The MISSING α₁ values
all_a1 = set(k[0] for k in lattice)
max_a1 = max(all_a1)
missing_a1 = sorted(set(range(max_a1+1)) - all_a1)
print(f"\n  Missing α₁ values at n=6: {missing_a1}")
print(f"  (These correspond to α₁=3 being impossible — HYP-926!)")

# ======================================================================
# PART 6: Cross-check with HYP-927 (H=21 impossible)
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: H=21 IMPOSSIBILITY MECHANISM")
print("=" * 70)

# H=21: 2α₁+4α₂ = 20, so α₁+2α₂ = 10.
# Possible decompositions:
# (α₁,α₂) = (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)
print(f"\n  H=21: α₁+2α₂ = 10")
print(f"  Possible decompositions:")
for a2 in range(6):
    a1 = 10 - 2*a2
    achievable = (a1, a2) in lattice
    print(f"    (α₁={a1:2d}, α₂={a2}): achievable? {achievable}")

# H=7: α₁+2α₂ = 3
print(f"\n  H=7: α₁+2α₂ = 3")
for a2 in range(2):
    a1 = 3 - 2*a2
    achievable = (a1, a2) in lattice
    print(f"    (α₁={a1:2d}, α₂={a2}): achievable? {achievable}")

# H=17: α₁+2α₂ = 8
print(f"\n  H=17: α₁+2α₂ = 8")
for a2 in range(5):
    a1 = 8 - 2*a2
    achievable = (a1, a2) in lattice
    print(f"    (α₁={a1:2d}, α₂={a2}): achievable? {achievable}")

# List all impossible H at n=6
impossible_H = set()
for h in range(1, max_H+2, 2):  # odd only
    needed = (h - 1) // 2  # α₁ + 2α₂ = needed
    found = False
    for a2 in range(needed // 2 + 1):
        a1 = needed - 2*a2
        if (a1, a2) in lattice:
            found = True
            break
    if not found:
        impossible_H.add(h)

print(f"\n  Impossible odd H values at n=6 (≤{max_H}):")
print(f"    {sorted(impossible_H)}")
print(f"  Compare with actual missing H: {sorted(set(range(1, max_H+1, 2)) - all_H)}")

# Check if they match
print(f"  Match: {sorted(impossible_H) == sorted(set(range(1, max_H+1, 2)) - all_H)}")

print("\nDone.")
