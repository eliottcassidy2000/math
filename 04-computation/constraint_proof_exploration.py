#!/usr/bin/env python3
"""
constraint_proof_exploration.py -- kind-pasteur-2026-03-13-S61

Exploring WHY c5_dir + 2*disj_33 = 56 for regular n=7 tournaments.

The constraint means: within the regular n=7 class, 5-cycle count and
3-cycle disjointness trade off EXACTLY. Why?

Approach 1: Count the same quantity two ways.
  What does c5_dir + 2*disj_33 count?
  - c5_dir = sum over 5-subsets of (directed Ham cycles from first vertex)
  - disj_33 = pairs of disjoint 3-cycle vertex sets

Approach 2: Relate to a global invariant.
  For regular tournaments, many invariants take fixed values.
  c3 = 14 vertex sets, c5 varies. But maybe some COMBINATION is fixed?

Approach 3: Double counting via directed paths.
  Count ordered triples (C3_a, C3_b, v) where C3_a, C3_b are disjoint 3-cycles
  and v is the remaining vertex. This gives 7*disj_33 (each pair leaves 1 vertex,
  there are 7 possible free vertices).

Actually: each disjoint 3-3 pair uses 6 of 7 vertices, leaving 1.
So 7*disj_33 / 7 = disj_33 disjoint pairs.

Let me try to connect c5_dir to disj_33 via a "complement" argument.

At n=7: a 5-subset S has complement {a,b} (2 vertices).
A 3-cycle vertex set T intersects S in |T cap S| vertices.
Since |T| = 3 and |S| = 5 and |S| + |T| = 8 > 7:
  |T cap S| >= 3 + 5 - 7 = 1. So every 3-cycle intersects every 5-subset.

But TWO disjoint 3-cycles T1, T2 cover 6 vertices, with complement {v}.
The 5-subset S containing both T1 and T2 minus v doesn't exist since |T1 cup T2| = 6 > 5.
Actually T1 subset S requires S to contain all 3 of T1's vertices,
and T2 subset S requires S to contain all 3 of T2's vertices.
T1 union T2 has 6 vertices, so no 5-subset contains both.

Let me think differently. For each pair of vertices (a,b), the complementary
5-subset is V \ {a,b}... no, the complement of a 5-subset is a 2-subset.

Actually, for each 5-subset of {0,...,6}, the complement is a pair of vertices.
And each 3-cycle vertex set either is fully contained in the 5-subset (3 of 5),
or has exactly 1 or 2 vertices outside.

Let me count: for a fixed 5-subset S, how many 3-cycle vertex sets are
fully contained in S? This is the number of 3-cycles in the subtournament
on S.

Author: kind-pasteur-2026-03-13-S61
"""

import math
from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_directed_ham_cycles_subset(A, verts):
    k = len(verts)
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


n = 7
m = n * (n - 1) // 2
total = 1 << m

print("=" * 70)
print("EXPLORING THE CONSTRAINT c5_dir + 2*disj_33 = 56")
print("=" * 70)

# First, collect data for regular tournaments
regular = []
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue

    H = count_ham_paths(A, n)

    # 3-cycle vertex sets
    c3_vsets = []
    for a, b, c in combinations(range(n), 3):
        cnt = count_directed_ham_cycles_subset(A, [a, b, c])
        if cnt > 0:
            c3_vsets.append(frozenset([a, b, c]))

    # Directed 5-cycles per 5-subset
    c5_per_subset = {}
    for subset in combinations(range(n), 5):
        cnt = count_directed_ham_cycles_subset(A, list(subset))
        c5_per_subset[frozenset(subset)] = cnt

    c5_dir = sum(c5_per_subset.values())

    # 3-cycles within each 5-subset
    c3_in_5sub = {}
    for subset in combinations(range(n), 5):
        fs = frozenset(subset)
        count = sum(1 for vs in c3_vsets if vs.issubset(fs))
        c3_in_5sub[fs] = count

    # Disjoint 3-3 pairs
    disj_33 = 0
    disj_33_by_free = defaultdict(int)  # free vertex -> count
    for i in range(len(c3_vsets)):
        for j in range(i+1, len(c3_vsets)):
            if not (c3_vsets[i] & c3_vsets[j]):
                disj_33 += 1
                free = set(range(n)) - c3_vsets[i] - c3_vsets[j]
                assert len(free) == 1
                disj_33_by_free[list(free)[0]] += 1

    regular.append({
        'bits': bits, 'H': H,
        'c3_vsets': len(c3_vsets), 'c5_dir': c5_dir,
        'disj_33': disj_33, 'disj_33_by_free': dict(disj_33_by_free),
        'c5_per_subset': c5_per_subset, 'c3_in_5sub': c3_in_5sub
    })
    break  # Just analyze one representative from each class first

# Actually, get one representative per H class
print("\nGathering one representative per H class...")
representatives = {}
for bits in range(total):
    A = binary_to_tournament(bits, n)
    scores = [sum(A[v]) for v in range(n)]
    if any(s != 3 for s in scores):
        continue
    H = count_ham_paths(A, n)
    if H in representatives:
        continue

    c3_vsets = []
    for a, b, c in combinations(range(n), 3):
        if count_directed_ham_cycles_subset(A, [a, b, c]) > 0:
            c3_vsets.append(frozenset([a, b, c]))

    c5_per_subset = {}
    c3_in_5sub = {}
    for subset in combinations(range(n), 5):
        fs = frozenset(subset)
        c5_per_subset[fs] = count_directed_ham_cycles_subset(A, list(subset))
        c3_in_5sub[fs] = sum(1 for vs in c3_vsets if vs.issubset(fs))

    c5_dir = sum(c5_per_subset.values())

    disj_33 = 0
    disj_33_by_free = defaultdict(int)
    for i in range(len(c3_vsets)):
        for j in range(i+1, len(c3_vsets)):
            if not (c3_vsets[i] & c3_vsets[j]):
                disj_33 += 1
                free = set(range(n)) - c3_vsets[i] - c3_vsets[j]
                disj_33_by_free[list(free)[0]] += 1

    representatives[H] = {
        'bits': bits, 'H': H,
        'c3_vsets': c3_vsets, 'c5_dir': c5_dir,
        'disj_33': disj_33, 'disj_33_by_free': dict(disj_33_by_free),
        'c5_per_subset': c5_per_subset, 'c3_in_5sub': c3_in_5sub
    }

    if len(representatives) == 3:
        break


# ========================================================================
# ANALYSIS 1: 5-SUBSET DECOMPOSITION
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 1: WHAT c5_dir COUNTS PER 5-SUBSET")
print("=" * 70)

for H in sorted(representatives.keys()):
    d = representatives[H]
    print(f"\n--- H={H} ---")
    print(f"  c5_dir = {d['c5_dir']}, disj_33 = {d['disj_33']}")
    print(f"  c5_dir + 2*disj_33 = {d['c5_dir'] + 2*d['disj_33']}")

    # For each 5-subset, show c5 and c3-within
    print(f"\n  5-subset analysis (complement = 2 vertices):")
    print(f"  {'complement':>12s} | {'c5_dir':>7s} | {'c3_in':>6s} | {'c3_in*(c3_in-1)/2':>18s}")
    print(f"  {'':->12s}-+-{'':->7s}-+-{'':->6s}-+-{'':->18s}")

    total_c5 = 0
    total_c3_in = 0
    total_c3_pairs_in = 0
    for subset in combinations(range(n), 5):
        fs = frozenset(subset)
        comp = set(range(n)) - set(subset)
        c5 = d['c5_per_subset'][fs]
        c3_in = d['c3_in_5sub'][fs]
        c3_pairs = c3_in * (c3_in - 1) // 2
        total_c5 += c5
        total_c3_in += c3_in
        total_c3_pairs_in += c3_pairs
        print(f"  {str(sorted(comp)):>12s} | {c5:>7d} | {c3_in:>6d} | {c3_pairs:>18d}")

    print(f"  {'TOTAL':>12s} | {total_c5:>7d} | {total_c3_in:>6d} | {total_c3_pairs_in:>18d}")

    # Key: each 3-cycle vertex set is in C(7-3, 5-3) = C(4,2) = 6 five-subsets
    # So sum of c3_in over all 5-subsets = 14 * 6 = 84
    print(f"\n  Check: each 3-cycle is in C(4,2)={math.comb(4,2)} five-subsets")
    print(f"  Expected sum(c3_in) = 14 * 6 = 84, actual = {total_c3_in}")

    # Each pair of 3-cycles (whether disjoint or not) is counted by
    # the number of 5-subsets containing both.
    # If they overlap in k vertices, they jointly use 6-k vertices.
    # A 5-subset containing both must contain all 6-k vertices,
    # choosing the remaining 5-(6-k) = k-1 from the 7-(6-k) = k+1 remaining.
    # So: C(k+1, k-1) five-subsets contain both.
    # k=0 (disjoint): C(1, -1) = 0. No 5-subset contains both!
    # k=1: C(2, 0) = 1.
    # k=2: C(3, 1) = 3.

    # So:
    # Sum of c3_pairs_in = sum over 5-subsets of C(c3_in, 2)
    # = (# pairs with overlap 1)*1 + (# pairs with overlap 2)*3
    # Disjoint pairs contribute 0 to this sum!

    # Total 3-3 pairs: C(14,2) = 91
    # = overlap-0 (disj_33) + overlap-1 + overlap-2
    # Since each 3-vertex set has 3 pairs of vertices, and each pair determines
    # the overlap with another 3-vertex set...

    # Let me compute the overlap distribution directly
    c3_vsets = d['c3_vsets']
    ov_dist = defaultdict(int)
    for i in range(len(c3_vsets)):
        for j in range(i+1, len(c3_vsets)):
            ov = len(c3_vsets[i] & c3_vsets[j])
            ov_dist[ov] += 1

    print(f"\n  3-cycle overlap distribution:")
    for ov in sorted(ov_dist.keys()):
        print(f"    overlap={ov}: {ov_dist[ov]} pairs")

    # Verify: sum_5sub C(c3_in, 2) = ov1*1 + ov2*3
    expected_pairs_in = ov_dist.get(1, 0) * 1 + ov_dist.get(2, 0) * 3
    print(f"  sum(C(c3_in,2)) = {total_c3_pairs_in}")
    print(f"  ov1*1 + ov2*3 = {expected_pairs_in}")
    print(f"  Match? {total_c3_pairs_in == expected_pairs_in}")

    # Now: can we relate c5_dir to the overlap distribution?
    # For each 5-subset S, c5_dir(S) is the number of directed Ham cycles on S.
    # The subtournament on S is a 5-vertex tournament with c3_in(S) three-cycles.
    # At n=5: if score = (2,2,2,2,2), c3=5, c5_dir=2.
    # If score = (1,2,2,2,3), c3=4, c5_dir in {1,2,3}.
    # So c5_dir(S) is NOT determined by c3_in(S) alone at n=5!

    # But maybe there's a FORMULA relating sum of c5_dir over all 5-subsets
    # to some global invariant?

    # At this point, let's just check the formula on the three classes.
    print(f"\n  c5_dir per 5-subset stats:")
    c5_vals = list(d['c5_per_subset'].values())
    c5_dist = defaultdict(int)
    for v in c5_vals:
        c5_dist[v] += 1
    print(f"    Distribution: {dict(sorted(c5_dist.items()))}")
    print(f"    Sum = {sum(c5_vals)}")


# ========================================================================
# ANALYSIS 2: LOOKING FOR THE COMBINATORIAL IDENTITY
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: COMBINATORIAL IDENTITY SEARCH")
print("=" * 70)

# c5_dir + 2*disj_33 = 56 for all regular n=7.
# 56 = C(7,2) * ... no, 56 = C(8,3) = 56. Hmm.
# 56 = 7 * 8 = 56. Yes! 7*8 = 56.
# 56 = 7 * C(4,2) = 7 * 6 = 42. No.
# Actually 56 = C(8,5) = C(8,3) = 56. And 56 = 7 * 8 = 56.

# Could this be related to the number of ordered vertex pairs times something?
# C(7,2) = 21, 2*21 = 42, 3*21 = 63. Not 56.
# C(7,3) = 35, 2*35 = 70. No.
# C(7,5) = 21, 2*21 + 14 = 56? 2*21 = 42 + 14 = 56. Hmm!
# 2*C(7,5) + C(7,3) = 42 + 35 = 77. No.
# Wait: 56 = 2*C(7,5) + c3_vsets? 2*21 + 14 = 56. YES!
# So c5_dir + 2*disj_33 = 2*C(n,5) + c3_vsets?
# At n=7: 2*21 + 14 = 56. YES!

# Is this a coincidence? The formula would be:
# c5_dir + 2*disj_33 = 2*C(n, n-2) + c3_vsets
# = 2*C(7,5) + 14
# But c3_vsets = 14 = n(n-1)(n+1)/24 for regular = 7*6*8/24 = 14.

# This is specific to regular n=7. Let's check if a similar formula holds elsewhere.

print(f"56 = 2*C(7,5) + c3_vsets = 2*21 + 14 = 56")
print(f"This would mean: c5_dir + 2*disj_33 = 2*C(n,n-2) + c3 for regular tournaments")

# What about non-regular n=6?
# At n=6, score (2,2,2,3,3,3):
# c3_vsets = 8, disj_33 varies (1-4), c5_dir varies (6-12)
# c5+2*disj = {12, 13, 14} (from the output). Not constant!
# So the constraint is specific to regular n=7 (or more generally, to some special class)

# Let me verify for the score (2,2,2,3,3,3) at n=6
n6 = 6
m6 = n6 * (n6 - 1) // 2
total6 = 1 << m6

combos_6 = set()
for bits in range(total6):
    A = binary_to_tournament(bits, n6)
    scores = tuple(sorted([sum(A[v]) for v in range(n6)]))
    if scores != (2, 2, 2, 3, 3, 3):
        continue

    c3_vsets = []
    for a, b, c in combinations(range(n6), 3):
        if count_directed_ham_cycles_subset(A, [a, b, c]) > 0:
            c3_vsets.append(frozenset([a, b, c]))

    c5_dir = sum(count_directed_ham_cycles_subset(A, list(subset))
                 for subset in combinations(range(n6), 5))

    disj_33 = 0
    for i in range(len(c3_vsets)):
        for j in range(i+1, len(c3_vsets)):
            if not (c3_vsets[i] & c3_vsets[j]):
                disj_33 += 1

    combos_6.add((c5_dir, disj_33, c5_dir + 2*disj_33))

print(f"\nn=6, score (2,2,2,3,3,3): c5_dir + 2*disj_33 = {sorted(set(x[2] for x in combos_6))}")
print(f"  (c5_dir, disj_33) pairs: {sorted(combos_6)}")

# Check if the constraint holds for any OTHER regular-like class at n=6
for sc_target in [(2, 2, 2, 2, 3, 4), (1, 2, 3, 3, 3, 3)]:
    combos = set()
    for bits in range(total6):
        A = binary_to_tournament(bits, n6)
        scores = tuple(sorted([sum(A[v]) for v in range(n6)]))
        if scores != sc_target:
            continue

        c3_vsets = []
        for a, b, c in combinations(range(n6), 3):
            if count_directed_ham_cycles_subset(A, [a, b, c]) > 0:
                c3_vsets.append(frozenset([a, b, c]))

        c5_dir = sum(count_directed_ham_cycles_subset(A, list(subset))
                     for subset in combinations(range(n6), 5))

        disj_33 = 0
        for i in range(len(c3_vsets)):
            for j in range(i+1, len(c3_vsets)):
                if not (c3_vsets[i] & c3_vsets[j]):
                    disj_33 += 1

        combos.add((c5_dir, disj_33, c5_dir + 2*disj_33))

    vals = sorted(set(x[2] for x in combos))
    print(f"\nn=6, score {sc_target}: c5_dir + 2*disj_33 = {vals}")
    if len(vals) == 1:
        print(f"  CONSTANT!")
    else:
        print(f"  NOT constant (varies)")


# ========================================================================
# ANALYSIS 3: WHAT MAKES REGULAR n=7 SPECIAL?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: REGULARITY AND THE CONSTRAINT")
print("=" * 70)

# For regular n=7: each vertex has score 3 (beats exactly 3 others).
# c3_vsets = 14 (always), c5_vsets = 21 (always = C(7,5)).
# This means EVERY 5-vertex subset supports at least one directed 5-cycle.
# At regular n=7, no 5-subset is "cycle-free" at the 5-cycle level.

# The constraint c5_dir + 2*disj_33 = 56 might follow from:
# 1. Every vertex is in exactly the same number of 3-cycles (by regularity)
# 2. The complementary structure is uniform

# Let me check: for each vertex v, how many 3-cycles contain v?
print(f"\n3-cycles per vertex (for regular n=7):")
for H in sorted(representatives.keys()):
    d = representatives[H]
    c3_per_v = defaultdict(int)
    for vs in d['c3_vsets']:
        for v in vs:
            c3_per_v[v] += 1
    print(f"  H={H}: {dict(sorted(c3_per_v.items()))}")

# For regular: each vertex is in C(6,2)/2 * (something)... actually:
# Each vertex v beats 3 others and loses to 3. A 3-cycle through v is
# either: v->a->b->v (a in out-nbrs, b in in-nbrs of v AND out-nbr of a)
# or similar. By regularity, each vertex is in exactly n(n-1)/6 = 7*6/6 = 7
# Wait: total 3-cycle vertex-sets = 14. Each vertex is in 14*3/7 = 6 of them.

# Check disj_33 per free vertex
print(f"\nDisjoint 3-3 pairs per free vertex:")
for H in sorted(representatives.keys()):
    d = representatives[H]
    print(f"  H={H}: disj_33={d['disj_33']}, by_free={d['disj_33_by_free']}")

# And c5_dir per excluded 2-set (which is the complement of a 5-subset)
print(f"\nc5_dir per 5-subset (grouped by complement pair):")
for H in sorted(representatives.keys()):
    d = representatives[H]
    by_comp = {}
    for subset in combinations(range(n), 5):
        comp = frozenset(set(range(n)) - set(subset))
        by_comp[comp] = d['c5_per_subset'][frozenset(subset)]
    print(f"\n  H={H}: c5_dir = {d['c5_dir']}")
    for comp in sorted(by_comp.keys()):
        print(f"    complement {sorted(comp)}: c5 = {by_comp[comp]}")


print("\n" + "=" * 70)
print("DONE.")
