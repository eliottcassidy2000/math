#!/usr/bin/env python3
"""
h21_cases_81_62.py — Investigate WHY (8,1) and (6,2) are impossible.

From the Phase Transition Table:
  alpha_1=8: alpha_2 ∈ {0} at n≤7, {0,7} at n=8+
  alpha_1=6: alpha_2 ∈ {0,1} at n≤7, {0,1,5} at n=8+

The required values (1 and 2 respectively) always fall in the gap.

opus-2026-03-14-S71e
"""

import sys
from itertools import combinations, permutations
from collections import defaultdict, Counter

sys.stdout.reconfigure(line_buffering=True)

def get_directed_cycles_3only(A, n):
    """Get 3-cycles only (fast for large n)."""
    groups = defaultdict(int)
    for verts in combinations(range(n), 3):
        v0, v1, v2 = verts
        if A[v0][v1] and A[v1][v2] and A[v2][v0]:
            groups[frozenset(verts)] += 1
        elif A[v0][v2] and A[v2][v1] and A[v1][v0]:
            groups[frozenset(verts)] += 1
    return groups

def get_directed_cycles(A, n, max_len=None):
    """Get all directed odd cycles."""
    if max_len is None:
        max_len = n
    groups = defaultdict(int)
    for length in range(3, min(n, max_len)+1, 2):
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
print("CASES (8,1) AND (6,2): STRUCTURAL ANALYSIS")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# Part 1: At n=6 exhaustive, alpha_1=8 structure
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 1: alpha_1=8 at n=6 ---")

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
ne = len(edges)

a8_examples = []
a6_examples = defaultdict(list)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    a1, a2 = compute_alpha(groups)

    if a1 == 8:
        a8_examples.append((bits, dict(groups), a2))
    if a1 == 6:
        a6_examples[a2].append((bits, dict(groups)))

print(f"  alpha_1=8: {len(a8_examples)} tournaments, alpha_2 = {Counter(x[2] for x in a8_examples)}")

# Analyze structure: how many 3-cycles, 5-cycles?
for bits, groups, a2 in a8_examples[:3]:
    t3 = sum(d for vs, d in groups.items() if len(vs) == 3)
    d5 = sum(d for vs, d in groups.items() if len(vs) == 5)
    print(f"\n  bits={bits}: t₃={t3}, d₅={d5}, alpha_2={a2}")
    for vs in sorted(groups.keys(), key=lambda x: (len(x), sorted(x))):
        d = groups[vs]
        print(f"    {sorted(vs)} d={d}")

    # Check which pairs are disjoint
    vs_list = list(groups.keys())
    disj_pairs = []
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            if not (vs_list[i] & vs_list[j]):
                disj_pairs.append((sorted(vs_list[i]), sorted(vs_list[j])))
    print(f"    Disjoint pairs: {len(disj_pairs)}")
    for a, b in disj_pairs:
        print(f"      {a} ⊥ {b}")

# ═══════════════════════════════════════════════════════════════════
# Part 2: alpha_1=6 at n=6 with alpha_2 values
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 2: alpha_1=6 at n=6 ---")

for a2_val in sorted(a6_examples.keys()):
    cases = a6_examples[a2_val]
    print(f"\n  alpha_2={a2_val}: {len(cases)} tournaments")
    bits, groups = cases[0]
    t3 = sum(d for vs, d in groups.items() if len(vs) == 3)
    d5 = sum(d for vs, d in groups.items() if len(vs) == 5)
    print(f"    Example (bits={bits}): t₃={t3}, d₅={d5}")
    for vs in sorted(groups.keys(), key=lambda x: (len(x), sorted(x))):
        d = groups[vs]
        print(f"    {sorted(vs)} d={d}")

    vs_list = list(groups.keys())
    disj_pairs = []
    for i in range(len(vs_list)):
        for j in range(i+1, len(vs_list)):
            if not (vs_list[i] & vs_list[j]):
                disj_pairs.append((sorted(vs_list[i]), sorted(vs_list[j])))
    print(f"    Disjoint pairs: {len(disj_pairs)}")

# ═══════════════════════════════════════════════════════════════════
# Part 3: The t₃ → alpha_2 forcing at n=6
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 3: t₃ → alpha_2 forcing ---")
print("  At n=6: how does t₃ (number of cyclic triples) constrain alpha_2?")

t3_to_a2 = defaultdict(Counter)

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles(A, n)
    t3 = sum(d for vs, d in groups.items() if len(vs) == 3)
    a1, a2 = compute_alpha(groups)

    t3_to_a2[t3][a2] += 1

print(f"  {'t₃':>4} | achievable alpha_2")
for t3 in sorted(t3_to_a2.keys()):
    a2_vals = sorted(t3_to_a2[t3].keys())
    total = sum(t3_to_a2[t3].values())
    print(f"  {t3:4d} | {a2_vals}  (n={total})")

# ═══════════════════════════════════════════════════════════════════
# Part 4: The number of complementary disjoint 3-cycle pairs
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 4: Disjoint 3-cycle pair count at n=6 ---")

# At n=6, two 3-cycles are disjoint iff they partition {0,...,5}.
# There are C(6,3)/2 = 10 complementary pairs.
# For a given tournament, how many complementary pairs are both cyclic?

comp_pair_count = Counter()

for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1

    groups = get_directed_cycles_3only(A, n)
    t3 = sum(groups.values())

    # Count complementary pairs where both are cyclic
    all_v = set(range(6))
    cp = 0
    for verts in combinations(range(6), 3):
        vs = frozenset(verts)
        comp = all_v - vs
        if vs < comp:  # avoid double-counting
            if vs in groups and comp in groups:
                cp += 1

    comp_pair_count[(t3, cp)] += 1

print(f"  (t₃, #complementary_both_cyclic): count")
for (t3, cp) in sorted(comp_pair_count.keys()):
    print(f"    t₃={t3:2d}, comp_pairs={cp}: {comp_pair_count[(t3,cp)]:5d}")

# The KEY insight: alpha_2 for 3-cycles-only at n=6 equals comp_pairs
# because the only way to have disjoint 3-cycles at n=6 is complementary pairs.

# ═══════════════════════════════════════════════════════════════════
# Part 5: The t₃ vs complementary pair relationship
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 5: t₃ vs disjoint pair relationship ---")
print("  At n=6: can we have t₃=8 with any complementary pair both cyclic?")

t3_8_with_cp = sum(v for (t3, cp), v in comp_pair_count.items() if t3 == 8 and cp > 0)
t3_8_total = sum(v for (t3, cp), v in comp_pair_count.items() if t3 == 8)
print(f"  t₃=8: {t3_8_with_cp}/{t3_8_total} have complementary pairs")

if t3_8_with_cp == 0:
    print("  PROVED: t₃=8 at n=6 → no complementary pair both cyclic!")
    print("  This means alpha_1=8 (3-cycles only) → alpha_2=0 at n=6. ✓")

# Check for other relevant t₃ values
for t3_target in [6, 7, 8, 9, 10]:
    total = sum(v for (t3, cp), v in comp_pair_count.items() if t3 == t3_target)
    with_cp = sum(v for (t3, cp), v in comp_pair_count.items() if t3 == t3_target and cp > 0)
    cp_vals = set(cp for (t3, cp), v in comp_pair_count.items() if t3 == t3_target and v > 0)
    print(f"  t₃={t3_target:2d}: {with_cp}/{total} have comp pairs, achievable cp = {sorted(cp_vals)}")

# ═══════════════════════════════════════════════════════════════════
# Part 6: STRUCTURAL REASON for the gap
# ═══════════════════════════════════════════════════════════════════
print("\n--- Part 6: Why does t₃=8 block complementary pairs? ---")

# At n=6, C(6,3)=20 triples. t₃=8 means 8 are cyclic, 12 are transitive.
# A complementary pair {a,b,c} and {d,e,f} being both cyclic means
# BOTH triples are cyclic. This uses 2 of the 8 cyclic triples.

# For t₃=8 with Σ C(s_i,2) = 20-8 = 12:
# Score sequences with Σ C(s_i,2) = 12:
# Sum of scores = C(6,2) = 15.
# Need 6 scores summing to 15 with Σ C(s_i,2) = 12.

scores_for_t3_8 = set()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    groups = get_directed_cycles_3only(A, n)
    t3 = sum(groups.values())
    if t3 == 8:
        scores_for_t3_8.add(scores)

print(f"  Score sequences giving t₃=8: {sorted(scores_for_t3_8)}")

# Check the complementary forcing:
# If {a,b,c} is cyclic (a→b→c→a), what does this mean for {d,e,f}?
# The 6 edges between {a,b,c} and {d,e,f} are independent of the
# internal structure. So a complementary pair can be both cyclic
# INDEPENDENTLY... unless the score sequence constrains it.

# Score of vertex a: out-degree to other 5 vertices.
# Within {a,b,c}: a beats exactly 1 vertex (in a 3-cycle).
# So a's score from internal edges = 1.
# External: a beats some of {d,e,f}.
# s_a = 1 + external_a.

# If complementary is also cyclic: each of d,e,f beats exactly 1 of {d,e,f}.
# So d beats 1 internally. Plus d beats some of {a,b,c} externally.
# s_d = 1 + external_d.

# Edge balance: Σ external_a (a∈{a,b,c}) + Σ external_d (d∈{d,e,f}) = 9
# (there are 9 cross-edges).
# Σ external_a = Σ (s_a - 1) for a∈{a,b,c} = Σ s_a - 3
# Σ external_d = Σ (s_d - 1) for d∈{d,e,f} = Σ s_d - 3
# Total: (Σ s_a - 3) + (Σ s_d - 3) = 15 - 6 = 9. ✓ (consistent)

# For both to be cyclic: each vertex has internal out-degree = 1.
# So s_i = 1 + external_i for all i. Minimum s_i = 1 (external=0).
# Maximum s_i = 4 (external=3, beats all cross-vertices).
# So all scores ∈ [1, 4]. Sum of scores = 15.

# For t₃=8: need Σ C(s_i,2) = 12.
# With scores in [1,4]: C(1,2)=0, C(2,2)=1, C(3,2)=3, C(4,2)=6.
# Need 6 scores summing to 15, each ∈ [1,4], with Σ C(s_i,2) = 12.

# Let's enumerate:
# Scores (a, b, c) for the two groups, each summing to some value,
# with Σ C(s_i,2) total = 12.

# If BOTH triples are cyclic, scores are in [1,4].
# If score (2,2,2,3,3,3): Σ C = 3+9 = 12, sum = 15. ✓
# But is this compatible with both complementary triples being cyclic?
# If {a,b,c} has scores {2,2,2}: each beats 1 internal + 1 external.
#   Internal out-degree always 1 (3-cycle), external = s-1 = 1.
#   So each of a,b,c beats exactly 1 of {d,e,f}.
# If {d,e,f} has scores {3,3,3}: each beats 1 internal + 2 external.
#   d beats 2 of {a,b,c}.
# Cross-edges from {a,b,c} to {d,e,f}: Σ external_a = 3.
# Cross-edges from {d,e,f} to {a,b,c}: Σ external_d = 6.
# Total cross-edges = 3+6 = 9. ✓

print(f"\n  Can score (2,2,2,3,3,3) have a complementary cyclic pair?")

# Check all tournaments with this score sequence
found = False
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    if scores == (2, 2, 2, 3, 3, 3):
        groups = get_directed_cycles_3only(A, n)
        t3 = sum(groups.values())
        # Check complementary pairs
        all_v = set(range(6))
        for verts in combinations(range(6), 3):
            vs = frozenset(verts)
            comp = all_v - vs
            if vs in groups and comp in groups:
                print(f"    YES! bits={bits}, {sorted(vs)} and {sorted(comp)} both cyclic, t₃={t3}")
                found = True
                break
    if found:
        break

if not found:
    print("    NO — score (2,2,2,3,3,3) never has both-cyclic complementary pair")

# Check ALL score sequences for both-cyclic complementary pairs
print(f"\n  ALL score sequences that allow complementary cyclic pairs:")
scores_with_cp = set()
for bits in range(2**ne):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            A[i][j] = 1
        else:
            A[j][i] = 1
    scores = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
    groups = get_directed_cycles_3only(A, n)
    all_v = set(range(6))
    for verts in combinations(range(6), 3):
        vs = frozenset(verts)
        comp = all_v - vs
        if vs in groups and comp in groups:
            scores_with_cp.add(scores)
            break

print(f"    {sorted(scores_with_cp)}")

# For these score sequences, what is t₃?
print(f"\n  Score sequences with comp pairs and their t₃ range:")
for ss in sorted(scores_with_cp):
    t3_vals = set()
    for bits in range(2**ne):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
        s = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        if s == ss:
            groups = get_directed_cycles_3only(A, n)
            t3 = sum(groups.values())
            all_v = set(range(6))
            for verts in combinations(range(6), 3):
                vs = frozenset(verts)
                comp = all_v - vs
                if vs in groups and comp in groups:
                    t3_vals.add(t3)
                    break

    sigma_c = sum(si*(si-1)//2 for si in ss)
    t3_formula = 20 - sigma_c
    print(f"    {ss}: t₃={t3_formula}, has comp pair: t₃∈{sorted(t3_vals) if t3_vals else 'none'}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
