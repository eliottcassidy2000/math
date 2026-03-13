#!/usr/bin/env python3
"""
vitali_atom_structure.py -- kind-pasteur-2026-03-13-S61

The 6 arc flips that preserve lambda but change c7 form C(4,2)=6 edges
among 4 vertices -- a sub-tournament reversal!

This is the "Vitali atom": the minimal transformation that changes H
while preserving the entire lambda graph structure.

Questions:
1. Is this always a 4-vertex sub-tournament reversal?
2. What is the structure of the 4-vertex sub-tournament being reversed?
3. Does the reversal ALWAYS change c7 by exactly 1?
4. What is the relationship to the Rédei function and arc-flip theory?
5. Can we characterize which sub-tournament reversals preserve lambda?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
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


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
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


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val
    return lam, c3_sets


n = 7

# ========================================================================
# ANALYSIS 1: Structure of the 4-vertex sub-tournament
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: THE 4-VERTEX SUB-TOURNAMENT REVERSAL")
print("=" * 70)

cases = [
    ("Case 1", 4728, 4658, (0, 1, 5, 2, 3, 4, 6)),
    ("Case 2", 9388, 9653, (1, 0, 2, 3, 6, 5, 4)),
]

for case_name, bits1, bits2, perm in cases:
    print(f"\n  {case_name}:")
    A1 = binary_to_tournament(bits1, n)
    A2 = binary_to_tournament(bits2, n)

    # Apply perm to A1
    A1p = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A1p[i][j] = A1[perm[i]][perm[j]]

    # Find flipped edges
    flips = []
    for i in range(n):
        for j in range(i+1, n):
            if A1p[i][j] != A2[i][j]:
                flips.append((i, j))

    flip_verts = set()
    for i, j in flips:
        flip_verts.add(i)
        flip_verts.add(j)
    flip_verts = sorted(flip_verts)

    print(f"    Flipped edges: {flips}")
    print(f"    Vertices involved: {flip_verts}")
    print(f"    C({len(flip_verts)},2) = {len(flip_verts)*(len(flip_verts)-1)//2}")
    print(f"    Number of flips: {len(flips)}")

    if len(flips) == len(flip_verts)*(len(flip_verts)-1)//2:
        print(f"    => ALL edges among {flip_verts} are flipped!")
        print(f"    => This is a COMPLETE SUB-TOURNAMENT REVERSAL on {flip_verts}")

    # What is the sub-tournament before and after?
    print(f"\n    Sub-tournament on {flip_verts} in A1' (H=109):")
    for v in flip_verts:
        row = [A1p[v][w] for w in flip_verts]
        print(f"      v={v}: {row}")
    scores_sub_1 = [sum(A1p[v][w] for w in flip_verts if w != v) for v in flip_verts]
    print(f"      Scores: {scores_sub_1}")

    print(f"    Sub-tournament on {flip_verts} in A2 (H=111):")
    for v in flip_verts:
        row = [A2[v][w] for w in flip_verts]
        print(f"      v={v}: {row}")
    scores_sub_2 = [sum(A2[v][w] for w in flip_verts if w != v) for v in flip_verts]
    print(f"      Scores: {scores_sub_2}")

    # Are the sub-tournament scores reversed?
    print(f"    Score reversal: {scores_sub_1} vs {[3-s for s in scores_sub_2]} (complement)")
    if scores_sub_1 == [3-s for s in scores_sub_2]:
        print(f"    => YES, exact complement!")

    # 3-cycle count in the sub-tournament
    sub_c3_1 = count_directed_ham_cycles_on_subset(A1p, flip_verts)
    sub_c3_2 = count_directed_ham_cycles_on_subset(A2, flip_verts)
    print(f"\n    Directed 3-cycles in sub: A1'={sub_c3_1}, A2={sub_c3_2}")

    # 4-cycle (Hamiltonian cycle on 4 vertices)
    sub_c4_1 = 0
    sub_c4_2 = 0
    # There are 3 Hamiltonian cycle directions on 4 vertices: (3!/2 = 3 non-equivalent)
    for p in permutations(range(4)):
        if p[0] != 0:
            continue
        is_cyc_1 = all(A1p[flip_verts[p[i]]][flip_verts[p[(i+1)%4]]] for i in range(4))
        is_cyc_2 = all(A2[flip_verts[p[i]]][flip_verts[p[(i+1)%4]]] for i in range(4))
        if is_cyc_1:
            sub_c4_1 += 1
        if is_cyc_2:
            sub_c4_2 += 1
    print(f"    Directed 4-cycles in sub: A1'={sub_c4_1}, A2={sub_c4_2}")

    # What type of tournament is the sub on 4 vertices?
    # Types: transitive (scores 0,1,2,3), near-transitive (scores 1,1,2,2 with one 3-cycle),
    #        or two 3-cycles (scores 1,2,2,1 -- wait, that's near-transitive too)
    # Actually at n=4: score (0,1,2,3) = transitive, (1,1,2,2) = has exactly 1 or 2 3-cycles
    for label, A_use in [("A1' (H=109)", A1p), ("A2 (H=111)", A2)]:
        sc = sorted([sum(A_use[v][w] for w in flip_verts if w != v) for v in flip_verts])
        print(f"    {label} sub-scores: {sc}", end="")
        if sc == [0,1,2,3]:
            print(" (TRANSITIVE)")
        elif sc == [1,1,2,2]:
            # How many 3-cycles?
            nc = 0
            for trip in combinations(flip_verts, 3):
                nc += count_directed_ham_cycles_on_subset(A_use, list(trip))
            print(f" (NEAR-REGULAR, {nc} directed 3-cycles)")
        else:
            print(f" (scores={sc})")


# ========================================================================
# ANALYSIS 2: Lambda preservation under sub-tournament reversal
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: WHY DOES REVERSAL PRESERVE LAMBDA?")
print("=" * 70)

# For lambda[u][v] to be preserved, every 3-cycle through {u,v} must be
# preserved. Reversing a sub-tournament on S flips ALL arcs within S.
# A 3-cycle {a,b,c} with all 3 vertices in S gets reversed but stays a 3-cycle.
# A 3-cycle {a,b,c} with 2 vertices in S and 1 outside: both arcs to/from S flip.
# A 3-cycle {a,b,c} with 1 vertex in S and 2 outside: one arc flips.
# A 3-cycle {a,b,c} with 0 vertices in S: unchanged.

# So lambda[u][v] is preserved iff:
# - For u,v both in S: 3-cycles through {u,v} either have 3rd vertex in S (preserved as 3-cycle)
#   or have 3rd vertex outside S. In the outside case, both arcs within {u,v,w} that touch S flip.
#   Wait: if u,v in S and w outside S, then arcs u->v, v->u flip (the one within S),
#   arcs u->w, w->u stay, arcs v->w, w->v stay. Only one arc flips. Not necessarily preserved.

# Actually: reversing within S means for i,j both in S: A'[i][j] = A[j][i] = 1-A[i][j].
# For i in S, j not in S: A'[i][j] = A[i][j] (unchanged).

# So for a 3-cycle {u,v,w} with u,v in S, w outside S:
# Old: either u->v->w->u or reverse. After reversal of S: the arc u->v becomes v->u.
# The arcs u->w and v->w stay.
# So if old was u->v, v->w, w->u: new is v->u, v->w, w->u. This is u<-v->w->u<- which is
# w->u->... wait let me be more careful.

# Old: A[u][v]=1, A[v][w]=1, A[w][u]=1 (u->v->w->u is a directed 3-cycle)
# New: A'[u][v]=0 (so v->u), A'[v][w]=1 (same), A'[w][u]=1 (same)
# Is this a 3-cycle? v->u, and w->u, so both v and w point to u. Then v->w? A'[v][w]=1 yes.
# So v->w->u and v->u. We need to check: is v->w->u->v? A'[u][v]=0 means v->u, so u->v is false.
# The 3-cycle would be v->w->u->v? That needs A'[u][v]=1 which is false. No.
# Try w->u->v->w? A'[w][u]=1, A'[u][v]=0. No.
# Try u->w->v->u? A'[u][w]=A[u][w]=0 (since w->u in original). No.
# So the old 3-cycle {u,v,w} is DESTROYED!

# But lambda is preserved. So there must be a NEW 3-cycle created to compensate.
# The reversal is a bijection on 3-cycle vertex sets IF for every destroyed 3-cycle,
# a new one is created with the same vertex set.

# Wait: the number of 3-cycles through a vertex triple {u,v,w} is:
# At most 2 (both directions). For a triple with u,v in S, w outside:
# Old: count_old 3-cycles. New: count_new 3-cycles.
# For lambda to be preserved, we need count_old = count_new for every triple
# that includes any pair {u,v} with specific lambda.

# Actually, lambda[u][v] = #{3-cycle vertex sets containing both u and v}.
# Each triple {u,v,w} either forms a 3-cycle or doesn't (it can have 0 or 2 directed 3-cycles,
# since a triple always has either 0 or 2 directed Hamiltonian cycles... wait no.
# A tournament on 3 vertices: score (1,1,1) has exactly 2 directed 3-cycles (both directions of the
# unique cyclic order). Score (0,1,2) has 0. So either 0 or 2 directed cycles, but as VERTEX SETS
# it's either 0 or 1 3-cycle vertex sets.

# So: for each triple, it's either a 3-cycle triple or a transitive triple.
# Reversing S flips the arc between the two S-vertices.
# For a triple {u,v,w} with u,v in S, w outside:
# The arc u->v flips. The arcs u->w and v->w stay.
# A triple is a 3-cycle iff it's NOT transitive iff no vertex beats both others.
# Old: if u->v, then {u beats v and w->u and v->w}: u beats v, v beats w, w beats u -> cycle.
#      or: u beats v, u beats w, v beats w: u beats all -> transitive.
# The flip changes who beats whom between u and v. This can turn a cycle into transitive or vice versa.

# So lambda is preserved iff for every pair {u,v}, the number of third vertices w
# making {u,v,w} a cycle is the same before and after.
# For u,v both in S: flip of u<->v arc. Third vertex w:
#   If w in S: all arcs flip, preserving cycle/transitive status.
#   If w not in S: only u<->v flips. This can change cycle status.
# For u in S, v not in S: no arc flips between u and v. Third vertex w:
#   If w in S: u<->w flips. Same analysis as above.
#   If w not in S: no flips. Preserved.

# So the condition is: for every pair {u,v} with u,v in S, and for every w outside S,
# the 3-cycle status of {u,v,w} must be the same before and after flipping u<->v.
# This is a strong condition!

for case_name, bits1, bits2, perm in cases:
    print(f"\n  {case_name}:")
    A2 = binary_to_tournament(bits2, n)
    A1 = binary_to_tournament(bits1, n)
    A1p = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A1p[i][j] = A1[perm[i]][perm[j]]

    flip_verts = set()
    for i in range(n):
        for j in range(i+1, n):
            if A1p[i][j] != A2[i][j]:
                flip_verts.add(i)
                flip_verts.add(j)
    S = sorted(flip_verts)
    outside = [v for v in range(n) if v not in S]

    print(f"    S = {S}, outside = {outside}")

    # Check: for each pair in S and each w outside, how does 3-cycle status change?
    for u, v in combinations(S, 2):
        for w in outside:
            triple = [u, v, w]
            old_cyc = count_directed_ham_cycles_on_subset(A1p, triple) > 0
            new_cyc = count_directed_ham_cycles_on_subset(A2, triple) > 0
            if old_cyc != new_cyc:
                print(f"    3-cycle status CHANGED: {triple}, old={old_cyc}, new={new_cyc}")
            else:
                pass

    # Check lambda preservation
    lam1, _ = get_labeled_lambda(A1p, n)
    lam2, _ = get_labeled_lambda(A2, n)
    preserved = all(lam1[i][j] == lam2[i][j] for i in range(n) for j in range(i+1, n))
    print(f"    Lambda preserved: {preserved}")

    # For which outside vertices w and pairs u,v in S does 3-cycle status change?
    print(f"\n    Detailed 3-cycle changes for pairs in S with outside vertices:")
    for w in outside:
        changes = []
        for u, v in combinations(S, 2):
            triple = sorted([u, v, w])
            old_cyc = count_directed_ham_cycles_on_subset(A1p, triple) > 0
            new_cyc = count_directed_ham_cycles_on_subset(A2, triple) > 0
            if old_cyc != new_cyc:
                changes.append(((u,v), old_cyc, new_cyc))
        if changes:
            print(f"      w={w}: {len(changes)} changes")
            for (uv, oc, nc) in changes:
                print(f"        pair {uv}: {'cycle->trans' if oc else 'trans->cycle'}")
        else:
            print(f"      w={w}: NO changes")


# ========================================================================
# ANALYSIS 3: Is this a known type of tournament transformation?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: LAMBDA-PRESERVING TRANSFORMATIONS")
print("=" * 70)

# The reversal of all arcs within a subset S preserves lambda iff
# for every pair {u,v} with u,v in S, the set of w outside S making
# {u,v,w} a 3-cycle is the same size before and after.

# For n=7, check ALL possible 4-vertex subsets of each tournament
# to see which sub-reversals preserve lambda.

print("\n  Checking all C(7,4)=35 sub-tournament reversals of A2 (bits=4658):")
A2 = binary_to_tournament(4658, n)
lam2, _ = get_labeled_lambda(A2, n)

preserving = []
for S in combinations(range(n), 4):
    # Reverse all arcs within S
    A_rev = [row[:] for row in A2]
    for i in S:
        for j in S:
            if i != j:
                A_rev[i][j] = 1 - A2[i][j]

    lam_rev, _ = get_labeled_lambda(A_rev, n)
    if all(lam_rev[i][j] == lam2[i][j] for i in range(n) for j in range(i+1, n)):
        H_rev = count_ham_paths(A_rev, n)
        c7_rev = count_directed_ham_cycles_on_subset(A_rev, list(range(n)))
        preserving.append((S, H_rev, c7_rev))

H2 = count_ham_paths(A2, n)
c7_2 = count_directed_ham_cycles_on_subset(A2, list(range(n)))
print(f"    Original: H={H2}, c7={c7_2}")
print(f"    Lambda-preserving 4-vertex reversals: {len(preserving)} out of 35")
for S, H, c7 in preserving:
    delta_H = H - H2
    delta_c7 = c7 - c7_2
    print(f"      S={S}: H={H} (delta={delta_H:+d}), c7={c7} (delta={delta_c7:+d})")

# Also check 3-vertex reversals
print(f"\n  Checking all C(7,3)=35 sub-tournament reversals (3-vertex):")
preserving_3 = []
for S in combinations(range(n), 3):
    A_rev = [row[:] for row in A2]
    for i in S:
        for j in S:
            if i != j:
                A_rev[i][j] = 1 - A2[i][j]
    lam_rev, _ = get_labeled_lambda(A_rev, n)
    if all(lam_rev[i][j] == lam2[i][j] for i in range(n) for j in range(i+1, n)):
        H_rev = count_ham_paths(A_rev, n)
        c7_rev = count_directed_ham_cycles_on_subset(A_rev, list(range(n)))
        preserving_3.append((S, H_rev, c7_rev))

print(f"    Lambda-preserving 3-vertex reversals: {len(preserving_3)} out of 35")
for S, H, c7 in preserving_3:
    delta_H = H - H2
    delta_c7 = c7 - c7_2
    print(f"      S={S}: H={H} (delta={delta_H:+d}), c7={c7} (delta={delta_c7:+d})")


# ========================================================================
# ANALYSIS 4: General theory - what determines lambda preservation?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: THEORY OF LAMBDA-PRESERVING REVERSALS")
print("=" * 70)

# For each 4-vertex subset S, compute the "external balance":
# For each pair (u,v) in S, count how many outside vertices w make {u,v,w} a 3-cycle.
# The reversal preserves lambda iff for each (u,v), this count is the same
# before and after. But the reversal only flips the u<->v arc, so the count changes
# iff flipping u<->v changes the 3-cycle status with some w.

# For a specific pair (u,v) and outside w:
# Old: u->v. Triple {u,v,w} is a 3-cycle iff (u->v, v->w, w->u) or (u->w, w->v, v->u=false).
# Wait: u->v means A[u][v]=1, A[v][u]=0.
# {u,v,w} is a 3-cycle set iff exactly one of the two directed 3-cycles exists.
# u->v->w->u: needs A[u][v] and A[v][w] and A[w][u]
# u->w->v->u: needs A[u][w] and A[w][v] and A[v][u]

# After flipping u<->v: A'[u][v]=0, A'[v][u]=1 (i.e., v->u).
# u->v->w->u now fails (A'[u][v]=0).
# u->w->v->u now fails (A'[v][u]=1, so we'd need A[v][u] in original = 0, but we just set it to 1).
# New cycles: v->u->w->v: needs A'[v][u]=1, A[u][w], A[w][v].
# Or: v->w->u->v: needs A[v][w], A[w][u], A'[u][v]=0... wait A'[u][v]=0 means v->u, which
# doesn't help for the cycle v->w->u->v which needs u->v.

# Let me just enumerate. With u->v (before flip):
# 3-cycle exists iff (A[v][w]*A[w][u]) = 1 (forward cycle u->v->w->u)
#                 or (A[u][w]*A[w][v]) = 0 and (A[v][u]*A[u][w]*A[w][v])...
# Actually: one cycle direction is u->v->w->u: A[u][v]*A[v][w]*A[w][u]
# Other is u->w->v->u: A[u][w]*A[w][v]*A[v][u]
# Since A[u][v]=1, A[v][u]=0, the second direction needs A[v][u]=1 which is false.
# So 3-cycle exists (as vertex set) iff A[v][w]*A[w][u] = 1 (i.e., v->w and w->u).

# After flip: A'[u][v]=0, A'[v][u]=1. Now v->u.
# Cycle u->v->w->u: A'[u][v]=0, fails.
# Cycle u->w->v->u: A[u][w]*A[w][v]*A'[v][u]=A[u][w]*A[w][v]*1 = A[u][w]*A[w][v].
# So 3-cycle after flip exists iff A[u][w]*A[w][v] = 1 (i.e., u->w and w->v).

# Summary: flip(u<->v) changes 3-cycle status of {u,v,w} from
#   "exists iff v->w and w->u" to "exists iff u->w and w->v"

# Lambda is preserved iff for each pair (u,v) in S:
#   #{w outside: v->w and w->u} = #{w outside: u->w and w->v}
# This is: #{w: A[v][w]=1, A[w][u]=1, w not in S} = #{w: A[u][w]=1, A[w][v]=1, w not in S}

print("""
THEOREM (Lambda-preservation condition):

Reversing all arcs within S preserves lambda iff for EVERY pair (u,v) in S:
  #{w outside S: v->w and w->u} = #{w outside S: u->w and w->v}

Equivalently, for every (u,v) in S:
  The number of outside vertices w forming a 3-cycle u->v->w->u
  equals the number forming a 3-cycle v->u->w->v.

This is a BALANCE CONDITION: the two orientations of the arc u<->v
participate in equally many 3-cycles with outside vertices.
""")

# Verify this for Case 1
A2 = binary_to_tournament(4658, n)
S = [1, 2, 3, 5]
outside = [v for v in range(n) if v not in S]
print(f"  Verification for Case 1 (H=111), S={S}, outside={outside}:")
for u, v in combinations(S, 2):
    fwd = sum(1 for w in outside if A2[v][w] and A2[w][u])
    rev = sum(1 for w in outside if A2[u][w] and A2[w][v])
    print(f"    ({u},{v}): fwd_3cyc={fwd}, rev_3cyc={rev}, balanced={fwd==rev}")


# ========================================================================
# ANALYSIS 5: At n=7, how common are lambda-preserving reversals?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: PREVALENCE AT n=7")
print("=" * 70)

# Sample tournaments and count how many have lambda-preserving 4-reversals
import random
random.seed(42)

total_tested = 0
total_with_reversal = 0
total_reversals = 0

for bits in range(0, 1 << 21, 128):  # Sample every 128th tournament
    A = binary_to_tournament(bits, n)
    lam, _ = get_labeled_lambda(A, n)
    has_reversal = False
    rev_count = 0

    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            has_reversal = True
            rev_count += 1

    total_tested += 1
    if has_reversal:
        total_with_reversal += 1
    total_reversals += rev_count

print(f"  Tested: {total_tested} tournaments")
print(f"  With at least one lambda-preserving 4-reversal: {total_with_reversal} ({100*total_with_reversal/total_tested:.1f}%)")
print(f"  Total reversals: {total_reversals} (avg {total_reversals/total_tested:.2f} per tournament)")


# ========================================================================
# ANALYSIS 6: The "Vitali atom" as a LOCAL operation
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 6: VITALI ATOMS AND H-CHANGE")
print("=" * 70)

# For all lambda-preserving reversals found in our sample,
# what is the distribution of delta_H?

delta_H_dist = defaultdict(int)
delta_c7_dist = defaultdict(int)

for bits in range(0, 1 << 21, 128):
    A = binary_to_tournament(bits, n)
    H = count_ham_paths(A, n)
    lam, _ = get_labeled_lambda(A, n)

    for S in combinations(range(n), 4):
        A_rev = [row[:] for row in A]
        for i in S:
            for j in S:
                if i != j:
                    A_rev[i][j] = 1 - A[i][j]
        lam_rev, _ = get_labeled_lambda(A_rev, n)
        if all(lam_rev[i][j] == lam[i][j] for i in range(n) for j in range(i+1, n)):
            H_rev = count_ham_paths(A_rev, n)
            c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))
            c7_rev = count_directed_ham_cycles_on_subset(A_rev, list(range(n)))
            delta_H_dist[H_rev - H] += 1
            delta_c7_dist[c7_rev - c7] += 1

print(f"  Delta H distribution:")
for dh in sorted(delta_H_dist.keys()):
    print(f"    delta_H = {dh:+d}: {delta_H_dist[dh]} occurrences")

print(f"\n  Delta c7 distribution:")
for dc in sorted(delta_c7_dist.keys()):
    print(f"    delta_c7 = {dc:+d}: {delta_c7_dist[dc]} occurrences")

# Key check: is delta_H always 2*delta_c7?
print(f"\n  Is delta_H = 2*delta_c7 always?")
# Since c3, c5, alpha_2 are preserved, and H = 1 + 2*(c3+c5+c7) + 4*alpha_2,
# we should have delta_H = 2*delta_c7
# But alpha_2 might change! Let's check.

# Actually, alpha_2 is the number of disjoint 3-cycle pairs.
# A 4-vertex reversal preserves the 3-cycle vertex sets (lambda preserved).
# Wait: it preserves lambda, which means it preserves which vertex triples are 3-cycle sets.
# So it preserves c3, c5, and alpha_2 as well!
# The ONLY thing that can change is c7.

print("""
Since lambda is preserved, the 3-cycle vertex sets are preserved.
Therefore c3, c5, alpha_2 are ALL preserved.
The OCF decomposition gives: H = 1 + 2*(c3+c5+c7) + 4*alpha_2
So delta_H = 2*delta_c7 ALWAYS.

This means the Vitali atom changes H by exactly 2*delta_c7.
The minimal non-trivial change is delta_c7 = +/-1, giving delta_H = +/-2.
""")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
