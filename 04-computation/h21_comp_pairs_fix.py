#!/usr/bin/env python3
"""
h21_comp_pairs_fix.py — Fixed complementary pair analysis at n=6.

Bug in h21_cases_81_62.py Part 4: `vs < comp` comparison fails
for equal-sized frozenset vs set. Fix: use canonical ordering.

Key question: For which t₃ values can complementary 3-cycle
vertex sets BOTH be cyclic?

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def get_directed_cycles(A, n):
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

def compute_alpha(groups):
    vs_list = list(groups.items())
    alpha1 = sum(d for _, d in vs_list)
    alpha2 = 0
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            if not (vs_list[i][0] & vs_list[j][0]):
                alpha2 += vs_list[i][1] * vs_list[j][1]
    return alpha1, alpha2

print("=" * 70)
print("COMPLEMENTARY PAIR ANALYSIS (FIXED)")
print("=" * 70)

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

# Generate all 10 complementary pairs
all_v = frozenset(range(6))
comp_pairs_list = []
for verts in combinations(range(6), 3):
    vs = frozenset(verts)
    comp = all_v - vs
    if sorted(vs)[0] < sorted(comp)[0]:  # canonical ordering
        comp_pairs_list.append((vs, comp))

print(f"  {len(comp_pairs_list)} complementary pairs")

# For each tournament, count how many complementary pairs are both cyclic
t3_to_cp = defaultdict(Counter)  # t₃ → {#comp_pairs: count}
a1_to_cp = defaultdict(Counter)  # alpha_1 → {#comp_pairs: count}

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    # Count 3-cycles
    t3 = sum(d for vs, d in groups.items() if len(vs) == 3)

    # Count complementary pairs where both are cyclic
    cp = 0
    for vs, comp in comp_pairs_list:
        if vs in groups and comp in groups:
            cp += 1

    t3_to_cp[t3][cp] += 1
    a1_to_cp[a1][cp] += 1

print(f"\n--- t₃ → complementary pair count ---")
for t3 in sorted(t3_to_cp.keys()):
    cp_dist = t3_to_cp[t3]
    cp_vals = sorted(cp_dist.keys())
    total = sum(cp_dist.values())
    with_cp = sum(v for k, v in cp_dist.items() if k > 0)
    print(f"  t₃={t3:2d}: comp_pairs ∈ {cp_vals}, {with_cp}/{total} have ≥1")

print(f"\n--- alpha_1 → complementary pair count ---")
for a1 in sorted(a1_to_cp.keys()):
    cp_dist = a1_to_cp[a1]
    cp_vals = sorted(cp_dist.keys())
    total = sum(cp_dist.values())
    with_cp = sum(v for k, v in cp_dist.items() if k > 0)
    print(f"  α₁={a1:2d}: comp_pairs ∈ {cp_vals}, {with_cp}/{total} have ≥1")

# KEY: For which alpha_1 values does alpha_2 > 0?
# alpha_2 counts disjoint pairs of cycle vertex sets.
# At n=6, the ONLY disjoint pairs are complementary 3-cycle pairs.
# WAIT: also need to check 3-cycle vs 5-cycle disjointness.

print(f"\n--- Disjoint pair sources at n=6 ---")
print(f"  A 3-cycle uses 3 of 6 vertices.")
print(f"  A 5-cycle uses 5 of 6 vertices.")
print(f"  3+3=6: disjoint iff complementary.")
print(f"  3+5=8>6: ALWAYS share ≥2 vertices.")
print(f"  5+5=10>6: ALWAYS share ≥4 vertices.")
print(f"  Therefore alpha_2 at n=6 = #complementary both-cyclic pairs (weighted).")

# Verify this claim
print(f"\n--- Verify: alpha_2 = weighted comp pair count ---")
mismatch = 0
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    # Compute alpha_2 from complementary pairs only
    a2_from_cp = 0
    for vs, comp in comp_pairs_list:
        if vs in groups and comp in groups:
            a2_from_cp += groups[vs] * groups[comp]

    if a2 != a2_from_cp:
        mismatch += 1
        if mismatch <= 3:
            print(f"  MISMATCH at bits={bits}: alpha_2={a2}, from_cp={a2_from_cp}")

if mismatch == 0:
    print(f"  CONFIRMED: alpha_2 = weighted complementary pair count at n=6")
else:
    print(f"  MISMATCHES: {mismatch}")

# Now the real question: for which t₃ can complementary pairs be both cyclic?
print(f"\n--- KEY: When are complementary 3-cycles both cyclic? ---")

# Score analysis: if {a,b,c} is cyclic, each vertex has out-degree 1 internally.
# If {d,e,f} is also cyclic, each has out-degree 1 internally.
# Total score s_v = internal_out + external_out.
# Each vertex has internal_out = 1 (in its cyclic triple).
# Sum of all scores = 15.
# Sum of internal_out = 6 (each of 6 vertices has internal_out=1).
# Sum of external_out = 15-6 = 9 (matching cross-edges).

# For score sequence: each s_v = 1 + external_v.
# external_v ∈ {0,1,2,3} (beats 0-3 of the 3 cross-vertices).
# Score range: s_v ∈ {1,2,3,4}.

# So tournaments with a both-cyclic complementary pair have all scores in [1,4].
# The score sum constraint: Σ s_v = 15, all s_v ∈ {1,4}.
# Score sequences: need 6 values in [1,4] summing to 15.

# Min sum: 6·1 = 6. Max sum: 6·4 = 24. 15 is in range.
# Average: 2.5.

# For t₃ = C(6,3) - Σ C(s_i,2) = 20 - Σ C(s_i,2).
# With scores in [1,4]: Σ C(s_i,2) ranges from min 6·C(2,2)=6 (all scores=2.5,
# nearest: 3 twos and 3 threes) to max: e.g., (1,1,4,4,4,1) → Σ=0+0+6+6+6+0=18.

# Hmm, but the score must come from a tournament with BOTH complementary triples cyclic.
# The cross-edges form a BIPARTITE tournament between {a,b,c} and {d,e,f}.
# This bipartite tournament has 9 edges, each directed.

# Let external_a = #(cross-edges from a to {d,e,f}).
# Total: Σ external_a = Σ external_d = 9... no.
# From {a,b,c} to {d,e,f}: some number k of edges go from abc side to def side.
# From {d,e,f} to {a,b,c}: 9-k edges.
# Σ external_a = k (edges from abc to def).
# Σ external_d = 9-k (edges from def to abc).
# Sum = k + (9-k) = 9. ✓

# Score constraint: Σ s_v = (Σ external) + 6 = 9 + 6 = 15. ✓

# Now: t₃ depends on ALL triples, not just the two complementary ones.
# The cross-triples (one from each group) can also be cyclic.
# Cross-triples: triples with 1 or 2 vertices from each group.

# 1 from {a,b,c} + 2 from {d,e,f}: C(3,1)*C(3,2) = 9 triples
# 2 from {a,b,c} + 1 from {d,e,f}: C(3,2)*C(3,1) = 9 triples
# Total cross-triples: 18
# Pure triples: 2 (the complementary pair itself)
# Total: 18+2 = 20 = C(6,3). ✓

# So t₃ = 2 (from the two cyclic complementary triples) + cross_cycles.
# For t₃ to be small, few cross-triples should be cyclic.

# Let's find the min t₃ that allows a both-cyclic complementary pair
min_t3_with_cp = float('inf')
max_t3_with_cp = 0
t3_dist_with_cp = Counter()

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    t3 = sum(d for vs, d in groups.items() if len(vs) == 3)

    for vs, comp in comp_pairs_list:
        if vs in groups and comp in groups:
            min_t3_with_cp = min(min_t3_with_cp, t3)
            max_t3_with_cp = max(max_t3_with_cp, t3)
            t3_dist_with_cp[t3] += 1
            break

print(f"  Min t₃ with both-cyclic comp pair: {min_t3_with_cp}")
print(f"  Max t₃ with both-cyclic comp pair: {max_t3_with_cp}")
print(f"  t₃ distribution (tournaments with ≥1 both-cyclic comp pair):")
for t3 in sorted(t3_dist_with_cp.keys()):
    print(f"    t₃={t3}: {t3_dist_with_cp[t3]} tournaments")

# CRITICAL: does this mean alpha_1=8 CANNOT have a both-cyclic comp pair?
# Because max t₃ at n=6 is 8, and we need to check if t₃=8 can have one.

print(f"\n--- CONCLUSION ---")
if 8 not in t3_dist_with_cp:
    print(f"  t₃=8 NEVER has a both-cyclic complementary pair!")
    print(f"  This is because t₃=8 requires score (2,2,2,3,3,3),")
    print(f"  and this score FORCES all complementary pairs to have")
    print(f"  one cyclic and one transitive triple.")
else:
    print(f"  t₃=8 CAN have both-cyclic complementary pairs!")

# Also show: what alpha_1 values have alpha_2>0?
print(f"\n--- alpha_1 values with alpha_2>0 at n=6 ---")
for a1 in sorted(a1_to_cp.keys()):
    cp_dist = a1_to_cp[a1]
    a2_vals_for_a1 = set()
    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
        groups = get_directed_cycles(A, n)
        a1_v, a2_v = compute_alpha(groups)
        if a1_v == a1:
            a2_vals_for_a1.add(a2_v)

    print(f"  α₁={a1:2d}: alpha_2 ∈ {sorted(a2_vals_for_a1)}")

print(f"\n  For H=21: need T=10 = alpha_1 + 2*alpha_2.")
print(f"  Check each alpha_1:")
for a1 in sorted(a1_to_cp.keys()):
    if a1 > 10:
        break
    if (10 - a1) % 2 == 0:
        needed_a2 = (10 - a1) // 2
        a2_vals = set()
        for bits in range(2**ne):
            A = [[0]*n for _ in range(n)]
            for idx, (i,j) in enumerate(edges):
                if bits & (1 << idx):
                    A[i][j] = 1
                else:
                    A[j][i] = 1
            groups = get_directed_cycles(A, n)
            a1_v, a2_v = compute_alpha(groups)
            if a1_v == a1:
                a2_vals.add(a2_v)

        status = "ACHIEVABLE" if needed_a2 in a2_vals else "IN GAP"
        print(f"    α₁={a1:2d}: need α₂={needed_a2}, achievable={sorted(a2_vals)}: {status}")
