#!/usr/bin/env python3
"""
h21_modular_constraint.py — Find the modular constraint that blocks H=21.

β₁+2β₂ skips the value 10 (goes 9→11 at n=6). WHY?

Key idea: β₁ = t₃ + d₅ + d₇ + ... where:
  t₃ = number of 3-cycles (each vertex set has exactly 1 directed 3-cycle)
  d₅ = total directed 5-cycles across all 5-element vertex sets
  d₇ = total directed 7-cycles across all 7-element vertex sets

Is there a modular constraint on t₃+d₅ that prevents certain values?

At n=5: β₁ ∈ {0,1,2,4,5,6,7} — note β₁=3 is MISSING!
This is exactly α₁=3 missing (HYP-1022: dc3=3 forces dc5≥1).

So β₁=3 impossible at n=5 ⟹ H=7 impossible at n=5.
Does β₁+2β₂=10 have a similar structural obstruction?

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def get_cycle_data_detailed(A, n):
    """Get t₃ (3-cycles), d₅ (directed 5-cycles), and their vertex sets."""
    t3 = 0
    t3_sets = []
    d5_sets = defaultdict(int)

    # 3-cycles
    for verts in combinations(range(n), 3):
        v0, v1, v2 = verts
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            t3 += 1
            t3_sets.append(frozenset(verts))
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            t3 += 1
            t3_sets.append(frozenset(verts))

    # 5-cycles (directed Hamiltonian cycles of 5-subtournaments)
    d5_total = 0
    for verts in combinations(range(n), 5):
        v0 = verts[0]
        count = 0
        for perm in permutations(verts[1:]):
            cycle = (v0,) + perm
            ok = True
            for i in range(5):
                if A[cycle[i]][cycle[(i+1) % 5]] != 1:
                    ok = False
                    break
            if ok:
                count += 1
        if count > 0:
            d5_sets[frozenset(verts)] = count
            d5_total += count

    return t3, t3_sets, d5_total, d5_sets

def count_disjoint_pairs_weighted(t3_sets, d5_sets):
    """Count β₂ = Σ d(S)·d(T) over disjoint pairs."""
    all_vs = []
    for s in t3_sets:
        all_vs.append((s, 1))
    for s, d in d5_sets.items():
        all_vs.append((s, d))

    # Merge: same vertex set → add multiplicities
    merged = defaultdict(int)
    for s, d in all_vs:
        merged[s] += d

    vs_list = list(merged.items())
    beta2 = 0
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            if not (vs_list[i][0] & vs_list[j][0]):
                beta2 += vs_list[i][1] * vs_list[j][1]

    return beta2

print("=" * 70)
print("MODULAR CONSTRAINT BLOCKING H=21")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: Decompose β₁ = t₃ + d₅ at n=5
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: β₁ = t₃ + d₅ decomposition at n=5 ---")

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

t3d5_counts = Counter()
t3d5_to_beta1 = defaultdict(set)
beta1_from_t3d5 = defaultdict(set)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3, t3_sets, d5, d5_sets = get_cycle_data_detailed(A, n)
    beta1 = t3 + d5
    t3d5_counts[(t3, d5)] += 1
    t3d5_to_beta1[(t3, d5)].add(beta1)
    beta1_from_t3d5[beta1].add((t3, d5))

print(f"  (t₃, d₅) → β₁ at n=5:")
for key in sorted(t3d5_counts.keys()):
    t3, d5 = key
    beta1 = t3 + d5
    print(f"    t₃={t3}, d₅={d5}: β₁={beta1}, count={t3d5_counts[key]}")

print(f"\n  β₁ decompositions:")
for b1 in sorted(beta1_from_t3d5.keys()):
    decomps = sorted(beta1_from_t3d5[b1])
    print(f"    β₁={b1}: from {decomps}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: At n=5, verify the GAP at β₁=3
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: The β₁=3 gap at n=5 ---")
print(f"  β₁=3 would need (t₃,d₅) with t₃+d₅=3:")
print(f"    (3,0): t₃=3 forces d₅≥1 (HYP-1022 — 3 directed 3-cycles force a 5-cycle)")
print(f"    (2,1): d₅=1 means one 5-set has exactly 1 directed Ham cycle")
print(f"    (1,2): d₅=2 means a 5-set has 2 directed Ham cycles")
print(f"    (0,3): d₅=3 means a 5-set has 3 directed Ham cycles")
print()
print(f"  Do (2,1), (1,2), (0,3) exist?")
for key in [(2,1), (1,2), (0,3), (3,0)]:
    count = t3d5_counts.get(key, 0)
    print(f"    ({key[0]},{key[1]}): {'EXISTS' if count else 'DOES NOT EXIST'} ({count} tournaments)")

# ═══════════════════════════════════════════════════════════════════
# Part 3: At n=6, decompose and find why β₁+2β₂=10 is impossible
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: β₁+2β₂ decomposition at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

target_counter = Counter()
target_decomp = defaultdict(list)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3, t3_sets, d5, d5_sets = get_cycle_data_detailed(A, n)
    beta1 = t3 + d5
    beta2 = count_disjoint_pairs_weighted(t3_sets, d5_sets)
    target = beta1 + 2*beta2

    target_counter[target] += 1
    if target in [9, 10, 11]:
        target_decomp[target].append((t3, d5, beta2))

print(f"  β₁+2β₂ distribution at n=6 (near 10):")
for t in range(max(0, 7), min(14, max(target_counter.keys())+1)):
    count = target_counter.get(t, 0)
    marker = " ← TARGET (H=21)" if t == 10 else ""
    print(f"    β₁+2β₂={t:2d}: {count:5d} tournaments{marker}")

# Detailed decomposition for target=9 and target=11
for t_val in [9, 10, 11]:
    entries = target_decomp.get(t_val, [])
    decomp_counter = Counter(entries)
    if entries:
        print(f"\n  β₁+2β₂={t_val} decompositions (t₃, d₅, β₂):")
        for key, cnt in sorted(decomp_counter.items()):
            t3, d5, b2 = key
            b1 = t3 + d5
            print(f"    t₃={t3}, d₅={d5}, β₂={b2} → β₁+2β₂={b1+2*b2}: {cnt} tournaments")
    else:
        print(f"\n  β₁+2β₂={t_val}: NO TOURNAMENTS")

# ═══════════════════════════════════════════════════════════════════
# Part 4: What t₃ values are possible at n=6?
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: t₃ distribution at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

t3_dist = Counter()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3 = 0
    for v0, v1, v2 in combinations(range(n), 3):
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            t3 += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            t3 += 1

    t3_dist[t3] += 1

print(f"  t₃ values at n=6: {sorted(t3_dist.keys())}")
for t3 in sorted(t3_dist.keys()):
    print(f"    t₃={t3:2d}: {t3_dist[t3]:5d} tournaments")

print("\n  Key: C(6,3)=20, each triple contributes exactly 1 directed 3-cycle")
print("  so t₃ = number of 'cyclic triples' in the tournament")
print("  t₃ ranges from 0 (transitive) to 10 (maximum for n=6)")

# ═══════════════════════════════════════════════════════════════════
# Part 5: The CRITICAL gap analysis
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: Why β₁+2β₂ skips 10 ---")
print("  At n=6 with no 7-cycles (n<7):")
print("  β₁ = t₃ + d₅")
print("  β₂ = Σ_{disjoint (S,T)} d(S)·d(T)")
print()
print("  For β₁+2β₂=10 we need:")
print("  Case β₂=0: β₁=10, so t₃+d₅=10. Max t₃ at n=6 is 10.")
print("  Case β₂=1: β₁=8, so t₃+d₅=8, plus exactly 1 disjoint weighted pair")
print("  etc.")
print()
print("  But t₃+d₅=10 with β₂=0 means all cycle vertex sets intersect.")
print("  At t₃=10 (max), does d₅=0 and β₂=0 hold?")

# Check tournaments with t₃=10 at n=6
n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3 = 0
    for v0, v1, v2 in combinations(range(n), 3):
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            t3 += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            t3 += 1

    if t3 != 10:
        continue

    t3_full, t3_sets, d5, d5_sets = get_cycle_data_detailed(A, n)
    beta1 = t3_full + d5
    beta2 = count_disjoint_pairs_weighted(t3_sets, d5_sets)

    print(f"  t₃=10: bits={bits}, β₁={beta1} (t₃={t3_full}, d₅={d5}), β₂={beta2}")
    print(f"         β₁+2β₂={beta1+2*beta2}")
    break  # Just show first one

# Check range where β₁+2β₂ could be 10
print("\n  Scanning all n=6 tournaments with β₁+2β₂ in [8,12]:")
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    t3, t3_sets, d5, d5_sets = get_cycle_data_detailed(A, n)
    beta1 = t3 + d5
    beta2 = count_disjoint_pairs_weighted(t3_sets, d5_sets)
    target = beta1 + 2*beta2

    if target == 10:
        print(f"  FOUND β₁+2β₂=10! bits={bits}, t₃={t3}, d₅={d5}, β₂={beta2}")

print("  (If no FOUND message, β₁+2β₂=10 is impossible at n=6)")
