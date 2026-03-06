#!/usr/bin/env python3
"""
Prove alpha_1=3 is impossible for ALL n.

Strategy: Show that c3=3 forces the 3 cyclic triples to span at most 5 vertices
and share a common vertex, for ANY n.

Proof attempt:
- Let T be a tournament on n vertices with exactly 3 cyclic triples.
- Consider the "triple graph" G: vertices = {a,b,c,...,n}, edges = union of triples.
- Each cyclic triple uses 3 vertices and 3 arcs forming a cycle.
- With 3 triples and <= 3*3 = 9 vertex slots, the triples can use at most 9 vertices.
- But actually they use at most 3*3 - (shared vertices) = ...

Key insight: Moon's formula c3 = C(n,3) - sum C(s_i, 2).
For c3 = 3: sum C(s_i, 2) = C(n,3) - 3.

Let's check this constraint and see what score sequences are possible.

Also: prove that the "triple graph" (union of 3 cyclic triples) has structure
that forces a common vertex.

kind-pasteur-2026-03-06-S22
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations
from collections import Counter

def get_cyclic_triples(T, n):
    result = []
    for combo in combinations(range(n), 3):
        a, b, c = combo
        if T[a][b] and T[b][c] and T[c][a]:
            result.append(frozenset(combo))
        elif T[a][c] and T[c][b] and T[b][a]:
            result.append(frozenset(combo))
    return result

# Verify at n=7: check all tournaments with c3=3
print("=" * 60)
print("n=7: VERIFYING alpha_1=3 impossibility")
print("=" * 60)

# n=7 has 2^21 = 2M tournaments - sample instead
import random
random.seed(42)
n = 7
m = n*(n-1)//2

found_c3_3 = 0
all_share_vertex = True
all_span_5 = True
score_dist = Counter()

# Systematic search: try to construct c3=3 tournaments
# Moon: c3 = C(7,3) - sum C(s_i, 2) = 35 - sum C(s_i, 2) = 3
# So sum C(s_i, 2) = 32
# Possible score sequences (sum = 7*6/2 = 21):
# For each partition of 21 into 7 non-negative integers <= 6...
# sum C(s_i, 2) = sum s_i*(s_i-1)/2 = 32
# => sum s_i^2 - sum s_i = 64
# => sum s_i^2 = 64 + 21 = 85

# Score sequences summing to 21, sum of squares = 85
from itertools import combinations_with_replacement

print("\nScore sequences with c3=3 (sum=21, sum_sq=85):")
valid_scores = []
for combo in combinations_with_replacement(range(7), 7):
    if sum(combo) == 21 and sum(s*s for s in combo) == 85:
        valid_scores.append(combo)
        print(f"  {combo}")

# Now sample tournaments and check those with c3=3
for trial in range(500000):
    bits = random.randint(0, (1 << m) - 1)
    T = tournament_from_bits(n, bits)

    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    if sum(s*s for s in scores) != 85:
        continue

    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue

    found_c3_3 += 1

    # Check common vertex
    common = triples[0] & triples[1] & triples[2]
    if not common:
        all_share_vertex = False
        print(f"  NO COMMON VERTEX! bits={bits}, triples={[sorted(t) for t in triples]}")

    # Check span
    all_verts = set()
    for t in triples:
        all_verts.update(t)
    if len(all_verts) > 5:
        all_span_5 = False
        print(f"  SPAN > 5! span={len(all_verts)}, bits={bits}")

    score_dist[scores] += 1

print(f"\nFound {found_c3_3} tournaments with c3=3 in 500000 samples")
print(f"All share common vertex? {all_share_vertex}")
print(f"All span <= 5 vertices? {all_span_5}")
print(f"Score distribution: {dict(score_dist)}")

# Now check n=8 similarly
print("\n" + "=" * 60)
print("n=8: VERIFYING")
print("=" * 60)

n = 8
m = n*(n-1)//2
# Moon: c3 = C(8,3) - sum C(s_i,2) = 56 - sum C(s_i,2) = 3
# sum C(s_i,2) = 53, so sum s_i^2 = 106 + 28 = 134 (sum s_i = 28)

print(f"Need sum s_i = {n*(n-1)//2}, sum s_i^2 = {2*(n*(n-1)*(n-2)//6 - 3) + n*(n-1)//2}")
target_sum = n*(n-1)//2
target_sum_sq = 2*(n*(n-1)*(n-2)//6 - 3) + target_sum

print(f"Score sequences with c3=3 at n=8 (sum={target_sum}, sum_sq={target_sum_sq}):")
valid_scores_8 = []
for combo in combinations_with_replacement(range(n), n):
    if sum(combo) == target_sum and sum(s*s for s in combo) == target_sum_sq:
        valid_scores_8.append(combo)
if len(valid_scores_8) <= 20:
    for s in valid_scores_8:
        print(f"  {s}")
else:
    print(f"  {len(valid_scores_8)} valid score sequences")

found_c3_3 = 0
all_share_vertex = True
all_span_5 = True

random.seed(123)
for trial in range(500000):
    bits = random.randint(0, (1 << m) - 1)
    T = tournament_from_bits(n, bits)

    scores = tuple(sorted(sum(T[i]) for i in range(n)))
    if sum(s*s for s in scores) != target_sum_sq:
        continue

    triples = get_cyclic_triples(T, n)
    if len(triples) != 3:
        continue

    found_c3_3 += 1

    common = triples[0] & triples[1] & triples[2]
    if not common:
        all_share_vertex = False
        print(f"  NO COMMON VERTEX! bits={bits}")

    all_verts = set()
    for t in triples:
        all_verts.update(t)
    if len(all_verts) > 5:
        all_span_5 = False
        print(f"  SPAN > 5! span={len(all_verts)}, bits={bits}")

print(f"\nFound {found_c3_3} tournaments with c3=3 in 500000 samples")
print(f"All share common vertex? {all_share_vertex}")
print(f"All span <= 5 vertices? {all_span_5}")

# THEORETICAL ARGUMENT
print("\n" + "=" * 60)
print("THEORETICAL ANALYSIS")
print("=" * 60)

# Three triples T1, T2, T3, each a set of 3 vertices.
# How can 3 sets of 3 elements intersect?
#
# Cases by union size:
# |T1 ∪ T2 ∪ T3| = 3: all three are the same set. But only 1 undirected
#   triple => at most 1 cyclic triple. Need 3 distinct triples.
# |union| = 4: 3 triples on 4 vertices. Each pair shares >= 2 vertices.
#   Possible: {a,b,c}, {a,b,d}, {a,c,d} (common vertex a, union = {a,b,c,d})
#   Or: {a,b,c}, {a,b,d}, {b,c,d} (pairwise share >= 1)
# |union| = 5: 3 triples on 5 vertices.
#   e.g. {a,b,c}, {a,d,e}, {b,c,d} - no common vertex!
# |union| = 6,7,...: triples can be more spread out.
#
# Key question: with c3=3 in a TOURNAMENT, are the triples forced to share
# a common vertex?

# Let's check: can we have 3 triples spanning 6+ vertices with c3=3?
# If triples are {a,b,c}, {d,e,f}, and one more:
# Two disjoint triples use 6 vertices. The third triple picks 3 from the
# remaining vertices. On n vertices, there are C(n,3) - 3 other triples.
# For these to ALL be transitive... that's very restrictive.
# But with n large enough, you might have room.

# Actually the constraint is: EXACTLY 3 cyclic triples total.
# If two triples are disjoint ({a,b,c} and {d,e,f}), the third triple
# must include vertices from both to stay at exactly c3=3.
# But adding a third triple on, say, {a,d,g} means we need:
# - {a,d,g} is cyclic
# - No other triple is cyclic
# This is possible if all other triples are transitive.

# BUT: does having {a,b,c} cyclic AND {d,e,f} cyclic AND {a,d,g} cyclic
# force any other triples to be cyclic?
# Not necessarily - it depends on the arc orientations.

# Let me construct a specific counterexample attempt at n=7:
# Try to build a tournament with c3=3 where triples don't share a common vertex.

print("\nAttempting to construct counterexample at n=7...")
print("Triples: {0,1,2}, {3,4,5}, {0,3,6} - no common vertex to all 3")

# We need: exactly these 3 cyclic, all C(7,3)-3 = 32 others transitive.
# {0,1,2} cyclic: 0->1->2->0 (or reverse)
# {3,4,5} cyclic: 3->4->5->3
# {0,3,6} cyclic: 0->3->6->0
# All other triples must be transitive.

# Let's try to build this and check.
# Fix: 0->1, 1->2, 2->0 (triple {0,1,2} cyclic)
# Fix: 3->4, 4->5, 5->3 (triple {3,4,5} cyclic)
# Fix: 0->3, 3->6, 6->0 (triple {0,3,6} cyclic)

# Now we need all other triples transitive.
# Already fixed arcs: 0->1, 1->2, 2->0, 3->4, 4->5, 5->3, 0->3, 3->6, 6->0
# Missing arcs to determine: all pairs not yet specified.
# Missing pairs: (0,4),(0,5),(1,3),(1,4),(1,5),(1,6),(2,3),(2,4),(2,5),(2,6),(3,5) already set,
# (4,6),(5,6)
# Wait, let me be more careful.

# Pairs and their arcs:
# (0,1): 0->1
# (0,2): 2->0
# (0,3): 0->3
# (0,4): ?
# (0,5): ?
# (0,6): 6->0
# (1,2): 1->2
# (1,3): ?
# (1,4): ?
# (1,5): ?
# (1,6): ?
# (2,3): ?
# (2,4): ?
# (2,5): ?
# (2,6): ?
# (3,4): 3->4
# (3,5): 5->3
# (3,6): 3->6
# (4,5): 4->5
# (4,6): ?
# (5,6): ?

# 12 undetermined arcs. We need all remaining triples to be transitive.
# Let's try a brute-force search over the 2^12 = 4096 possibilities.

fixed = {}
# 0->1, 1->2, 2->0, 3->4, 4->5, 5->3, 0->3, 3->6, 6->0
fixed[(0,1)] = 1; fixed[(1,0)] = 0
fixed[(1,2)] = 1; fixed[(2,1)] = 0
fixed[(2,0)] = 1; fixed[(0,2)] = 0
fixed[(3,4)] = 1; fixed[(4,3)] = 0
fixed[(4,5)] = 1; fixed[(5,4)] = 0
fixed[(5,3)] = 1; fixed[(3,5)] = 0
fixed[(0,3)] = 1; fixed[(3,0)] = 0
fixed[(3,6)] = 1; fixed[(6,3)] = 0
fixed[(6,0)] = 1; fixed[(0,6)] = 0

free_pairs = [(0,4),(0,5),(1,3),(1,4),(1,5),(1,6),(2,3),(2,4),(2,5),(2,6),(4,6),(5,6)]
assert len(free_pairs) == 12

found_counterexample = False
for bits in range(1 << 12):
    T = [[0]*7 for _ in range(7)]
    # Set fixed arcs
    for (i,j), val in fixed.items():
        T[i][j] = val
    # Set free arcs
    for k, (i,j) in enumerate(free_pairs):
        if bits & (1 << k):
            T[i][j] = 1
            T[j][i] = 0
        else:
            T[j][i] = 1
            T[i][j] = 0

    # Count cyclic triples
    triples = get_cyclic_triples(T, 7)
    if len(triples) == 3:
        # Check if exactly our 3 triples
        expected = {frozenset({0,1,2}), frozenset({3,4,5}), frozenset({0,3,6})}
        if set(triples) == expected:
            common = triples[0] & triples[1] & triples[2]
            all_verts = set()
            for t in triples:
                all_verts.update(t)
            print(f"  Found! bits={bits}, span={len(all_verts)}, common={common}")
            if not common:
                found_counterexample = True
                print(f"  *** COUNTEREXAMPLE: No common vertex!")
                # Check c5
                # Count directed 5-cycles
                from itertools import combinations as combos
                c5 = 0
                for combo in combos(range(7), 5):
                    verts = list(combo)
                    dp = {}
                    dp[(1, 0)] = 1
                    for mask in range(1, 1 << 5):
                        if not (mask & 1): continue
                        for vi in range(5):
                            if not (mask & (1 << vi)): continue
                            c = dp.get((mask, vi), 0)
                            if c == 0: continue
                            for ui in range(5):
                                if mask & (1 << ui): continue
                                if T[verts[vi]][verts[ui]]:
                                    key = (mask | (1 << ui), ui)
                                    dp[key] = dp.get(key, 0) + c
                    full = (1 << 5) - 1
                    for vi in range(1, 5):
                        c = dp.get((full, vi), 0)
                        if c > 0 and T[verts[vi]][verts[0]]:
                            c5 += c
                print(f"  c5={c5}, alpha_1={3+c5}")
                break

if not found_counterexample:
    print("  No tournament found with these 3 specific triples and c3=3")
    print("  (Might have other triples becoming cyclic)")

# Let's also check: can ANY 3 triples without common vertex give c3=3 at n=7?
print("\nSearching for ANY c3=3 tournament at n=7 with triples not sharing common vertex...")
random.seed(999)
found_no_common = False
checked = 0
for trial in range(1000000):
    bits = random.randint(0, (1 << 21) - 1)
    T = tournament_from_bits(7, bits)
    scores = tuple(sorted(sum(T[i]) for i in range(7)))
    if sum(s*s for s in scores) != 85:
        continue
    triples = get_cyclic_triples(T, 7)
    if len(triples) != 3:
        continue
    checked += 1
    common = triples[0] & triples[1] & triples[2]
    if not common:
        found_no_common = True
        print(f"  FOUND: bits={bits}, triples={[sorted(t) for t in triples]}")
        break

print(f"Checked {checked} c3=3 tournaments at n=7")
print(f"Found without common vertex? {found_no_common}")

print("\nDone.")
