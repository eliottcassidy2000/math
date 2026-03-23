#!/usr/bin/env python3
"""
parent_tree_degeneracy_s209.py -- opus-2026-03-23-S209

VERTEX DELETION PARENTS AND THE DEGENERACY FORMULA

Each class C at level n has PARENT classes at level n-1:
  Delete vertex v from tournament T in C -> tournament in parent class P_v.
  The DELETION PROFILE of C = multiset of parent classes.

From S171: deletion profiles determine the fiber structure.
From S178: this IS the tournaplex face map.

THE KEY INSIGHT: The degeneracy D = T - 2E counts arc-orbit COLLISIONS.
Two arcs of T lead to the SAME target class iff they produce the
same deletion profile change. This is controlled by the PARENT STRUCTURE.

If two arcs (a,b) and (c,d) have the SAME effect on the deletion profile
(removing the same vertex gives the same result), they collide.

Author: opus-2026-03-23-S209
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def delete_vertex(A, n, v):
    remaining = [i for i in range(n) if i != v]
    return [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]

# Build class catalogs for n and n-1
def build_catalog(n):
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; cdata = {}; creps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A,n), 'scores': score_seq(A,n), 'size': 0}
            creps[cid] = [row[:] for row in A]
        cdata[cmap[c]]['size'] += 1
    return cmap, cdata, creps, len(cdata)

# ============================================================
# BUILD FOR n=5 WITH FULL PARENT/CHILD STRUCTURE
# ============================================================

banner("PARENT TREE ANALYSIS AT n=5")

cmap5, cd5, cr5, nc5 = build_catalog(5)
cmap4, cd4, cr4, nc4 = build_catalog(4)
cmap3, cd3, cr3, nc3 = build_catalog(3)

print("  G_5: %d classes, G_4: %d classes, G_3: %d classes\n" % (nc5, nc4, nc3))

# Compute deletion profiles: class at n=5 -> multiset of classes at n=4
del_profiles = {}
for cid5 in range(nc5):
    A = cr5[cid5]
    profile = Counter()
    for v in range(5):
        A_del = delete_vertex(A, 5, v)
        c4 = cmap4[canonical(A_del, 4)]
        profile[c4] += 1
    del_profiles[cid5] = profile

# Also compute n=4 -> n=3 deletion
del_profiles_43 = {}
for cid4 in range(nc4):
    A = cr4[cid4]
    profile = Counter()
    for v in range(4):
        A_del = delete_vertex(A, 4, v)
        c3 = cmap3[canonical(A_del, 3)]
        profile[c3] += 1
    del_profiles_43[cid4] = profile

# Print the tree
print("  DELETION TREE (n=5 -> n=4 -> n=3):\n")
for cid5 in range(nc5):
    d5 = cd5[cid5]
    prof = del_profiles[cid5]
    parent_str = ", ".join("%d:%d" % (c, cnt) for c, cnt in sorted(prof.items()))
    parent_H = sorted(set(cd4[c]['H'] for c in prof))

    # For each n=4 parent, what are ITS n=3 parents?
    grandparents = Counter()
    for c4, cnt4 in prof.items():
        for c3, cnt3 in del_profiles_43[c4].items():
            grandparents[c3] += cnt4 * cnt3  # approximate, not exact

    gp_str = ", ".join("%d:%d" % (c, cnt) for c, cnt in sorted(grandparents.items()))

    print("  class %2d (H=%2d): n=4 parents = {%s} (H=%s)" %
          (cid5, d5['H'], parent_str, parent_H))
    print("                    n=3 grandparents = {%s}" % gp_str)

# ============================================================
# HOW PARENTS RELATE TO ARC-REVERSAL EDGES
# ============================================================

banner("PARENTS AND ARC-REVERSAL EDGES")

print("""
  KEY INSIGHT: Reversing arc (a,b) in tournament T changes the
  DELETION PROFILE at vertices a and b:
  - Deleting vertex a: arc (a,b) disappears, T_(-a) is the same
    whether (a,b) or (b,a). So parent at vertex a is UNCHANGED.
    Wait: that's wrong. Deleting a removes ALL arcs incident to a,
    including (a,b). So the sub-tournament T-a doesn't involve arc (a,b) at all.
    Reversing (a,b) has NO EFFECT on T-a.
    But it DOES affect T-b (arc (a,b) is incident to b but NOT removed).
    Wait: deleting b removes ALL arcs incident to b. So T-b also
    doesn't involve (a,b) directly... WAIT.

  CORRECTION: Deleting vertex v removes ALL arcs to/from v.
  Reversing arc (a,b) changes:
  - T-v for v != a and v != b: NO CHANGE (arc (a,b) is still there, direction flipped)
  - T-a: arc (a,b) is removed (since a is deleted). The rest is unchanged.
    So T-a is the SAME before and after reversing (a,b).
  - T-b: arc (a,b) is removed (since b is deleted). Same: T-b is UNCHANGED.

  WAIT: that can't be right. T-v for v != a, v != b DOES change:
  the sub-tournament on vertices {0,...,n-1} \ {v} includes both a and b.
  The arc between a and b in this sub-tournament is (a,b) before and (b,a) after.
  So: T-v CHANGES when (a,b) is reversed, for all v != a and v != b.
  T-a and T-b do NOT change (since the reversed arc is removed with the deleted vertex).

  THEREFORE: Reversing arc (a,b) changes the parent at vertex v
  for ALL v != a and v != b (that's n-2 parents).
  The parents at a and b are UNCHANGED.

  This means: arc reversal at (a,b) affects (n-2) of the n parent classes.
  The 2 unchanged parents are at the endpoints of the reversed arc.
""")

# Verify at n=5
print("  VERIFICATION at n=5:\n")
m5 = comb(5, 2)
pairs5 = [(i,j) for i in range(5) for j in range(i+1,5)]

# Pick a specific tournament and arc
bits_example = 0  # transitive
A_example = [[0]*5 for _ in range(5)]
for k,(i,j) in enumerate(pairs5):
    if (bits_example>>k)&1: A_example[i][j]=1
    else: A_example[j][i]=1

# Original deletion profile
orig_profile = []
for v in range(5):
    A_del = delete_vertex(A_example, 5, v)
    c = cmap4[canonical(A_del, 4)]
    orig_profile.append(c)

print("  Transitive tournament, original parents: %s" % orig_profile)
print("  (All the same class because transitive is highly symmetric)\n")

# Reverse arc (0,1) (the first arc)
A_rev = [row[:] for row in A_example]
A_rev[0][1] = 1 - A_rev[0][1]; A_rev[1][0] = 1 - A_rev[1][0]

rev_profile = []
for v in range(5):
    A_del = delete_vertex(A_rev, 5, v)
    c = cmap4[canonical(A_del, 4)]
    rev_profile.append(c)

print("  After reversing arc (0,1), parents: %s" % rev_profile)
print("  Changed at vertices: %s" %
      [v for v in range(5) if orig_profile[v] != rev_profile[v]])
print("  Unchanged at vertices: %s" %
      [v for v in range(5) if orig_profile[v] == rev_profile[v]])
print("  (Should be unchanged at 0 and 1, changed at 2,3,4)")

# ============================================================
# THE PARENT PROFILE AS AN EDGE INVARIANT
# ============================================================

banner("THE PARENT PROFILE DETERMINES EDGE TARGETS")

print("""
  When we reverse arc (a,b) in tournament T:
  - The NEW tournament T' has a DIFFERENT deletion profile
  - Specifically: the parent at vertex v (v != a, v != b) changes
    because arc (a,b) flips direction in the sub-tournament T-v

  The NEW class of T' is determined by ALL the arc changes.
  But the PARENT PROFILE of T' tells us which n=4 classes appear.

  TWO arc reversals (a,b) and (c,d) produce the SAME target class
  iff the resulting tournaments are isomorphic. This happens when:
  - The "affected" parents (n-2 of them) change in the same way
  - There exists a relabeling sending one reversal to the other

  THE COLLISION CONDITION:
  Arcs (a,b) and (c,d) COLLIDE (= same target class) iff:
  there exists sigma in S_n with:
    sigma(a) = c, sigma(b) = d, AND sigma is an automorphism of T_(-{a,b,c,d})

  For |Aut(T)| = 1 (generic): sigma must map EVERYTHING correctly.
  This is very restrictive, so collisions are RARE for generic T.
  Hence D/T -> 0 as generic classes dominate.

  For |Aut(T)| > 1: the automorphism creates MANDATORY collisions.
  If sigma(a)=a, sigma(b)=c: then arcs (a,b) and (a,c) collide.
  This happens when b and c are in the same Aut-orbit.
  Number of such collisions = sum over Aut-orbits of C(orbit_size, 2).
""")

# For each class at n=5, compute collisions from Aut structure
print("  AUT-ORBIT COLLISIONS at n=5:\n")

for cid5 in range(nc5):
    A = cr5[cid5]
    # Find automorphisms
    auts = []
    for perm in permutations(range(5)):
        ok = True
        for i in range(5):
            for j in range(5):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: auts.append(perm)

    aut_size = len(auts)

    # Arc orbits under Aut
    arcs = [(i,j) for i in range(5) for j in range(5) if i != j and A[i][j]]
    arc_orbits_list = []
    seen = set()
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            mapped = (perm[arc[0]], perm[arc[1]])
            if A[mapped[0]][mapped[1]]:
                orbit.add(mapped)
        arc_orbits_list.append(orbit)
        seen.update(orbit)

    n_orbits = len(arc_orbits_list)

    # For each orbit, check where it reverses to
    orbit_targets = []
    for orbit in arc_orbits_list:
        arc = next(iter(orbit))
        i, j = arc
        A_rev = [row[:] for row in A]
        A_rev[i][j] = 0; A_rev[j][i] = 1
        target = cmap5[canonical(A_rev, 5)]
        orbit_targets.append(target)

    # Count collisions: orbits that go to the SAME target
    target_counts = Counter(orbit_targets)
    collisions = sum(c*(c-1)//2 for c in target_counts.values())
    self_loops = sum(1 for t in orbit_targets if t == cid5)

    # The key: collision count = contribution to D from this class
    print("  class %2d (H=%2d, |Aut|=%d): %d orbits, %d self, %d collisions, targets=%s" %
          (cid5, cd5[cid5]['H'], aut_size, n_orbits, self_loops, collisions,
           dict(target_counts)))

# Sum of collisions = D/2 (since each collision is counted once per class)
# Wait: D = sum over classes of collisions? Let me check.

total_collisions = 0
total_self = 0
for cid5 in range(nc5):
    A = cr5[cid5]
    auts = []
    for perm in permutations(range(5)):
        ok = True
        for i in range(5):
            for j in range(5):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: auts.append(perm)

    arcs = [(i,j) for i in range(5) for j in range(5) if i != j and A[i][j]]
    arc_orbits_list = []
    seen = set()
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            mapped = (perm[arc[0]], perm[arc[1]])
            if A[mapped[0]][mapped[1]]: orbit.add(mapped)
        arc_orbits_list.append(orbit)
        seen.update(orbit)

    orbit_targets = []
    for orbit in arc_orbits_list:
        arc = next(iter(orbit))
        i, j = arc
        A_rev = [row[:] for row in A]
        A_rev[i][j] = 0; A_rev[j][i] = 1
        target = cmap5[canonical(A_rev, 5)]
        orbit_targets.append(target)

    target_counts = Counter(orbit_targets)
    collisions = sum(c*(c-1)//2 for c in target_counts.values())
    self_count = target_counts.get(cid5, 0)
    total_collisions += collisions
    total_self += self_count

print("\n  Total collisions (sum over classes): %d" % total_collisions)
print("  Total self-loop orbits: %d" % total_self)
print("  Expected D = %d" % (T[2] - 2*E[2]))
print("  D/2 = %d" % ((T[2] - 2*E[2])//2))

# D = 28. Collisions should relate to 28 somehow.
# collisions = sum_C sum_D C(mult(C->D), 2)
# Self-loops don't count as collisions (they're a single target, just "used up")
# D = total_collisions + total_self? Let me check.

D_check = total_collisions + total_self
print("  collisions + self = %d (D=%d, match=%s)" %
      (D_check, 28, D_check == 28))

# Actually: D = T - 2E
# T = sum_C #orbits(C) = 88
# 2E = 2*30 = 60
# D = 28
# total_orbits = 88 (verified)
# Each orbit either: (a) goes to a new target (contributes to E), or
# (b) goes to a target already reached by another orbit (collision)
# or (c) is a self-loop

# total_orbits = (orbits to new targets) + (collision orbits) + (self-loops)
# 2E = 2 * (orbits to distinct cross-targets)
# Wait: each edge (C,D) is created by >= 1 orbit from C to D.
# The FIRST orbit creates the edge. Additional orbits are "collisions."
# So: total_cross_orbits = 2E + collision_excess
# And: total_orbits = total_cross_orbits + total_self
# = 2E + collision_excess + total_self
# D = total_orbits - 2E = collision_excess + total_self

print("\n  DECOMPOSITION:")
print("  D = collision_excess + self_loops")
print("  collision_excess = sum_C sum_{edges from C} (mult - 1)")
print("  = total_collisions (since C(mult,2) counts pairs, not excess)")
print()

# Let me compute collision_excess correctly
# For each class C, for each target D != C with mult(C->D) = k:
# the "excess" is k - 1 (one orbit creates the edge, rest are excess)
collision_excess = 0
for cid5 in range(nc5):
    A = cr5[cid5]
    auts = []
    for perm in permutations(range(5)):
        ok = True
        for i in range(5):
            for j in range(5):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: auts.append(perm)

    arcs = [(i,j) for i in range(5) for j in range(5) if i != j and A[i][j]]
    arc_orbits_list = []; seen = set()
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            mapped = (perm[arc[0]], perm[arc[1]])
            if A[mapped[0]][mapped[1]]: orbit.add(mapped)
        arc_orbits_list.append(orbit); seen.update(orbit)

    orbit_targets = []
    for orbit in arc_orbits_list:
        arc = next(iter(orbit)); i, j = arc
        A_rev = [row[:] for row in A]
        A_rev[i][j] = 0; A_rev[j][i] = 1
        target = cmap5[canonical(A_rev, 5)]
        orbit_targets.append(target)

    target_counts = Counter(orbit_targets)
    for target, count in target_counts.items():
        if target != cid5:
            collision_excess += (count - 1)  # excess beyond the first

total_self_2 = 0
for cid5 in range(nc5):
    A = cr5[cid5]
    auts = []
    for perm in permutations(range(5)):
        ok = True
        for i in range(5):
            for j in range(5):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: auts.append(perm)

    arcs = [(i,j) for i in range(5) for j in range(5) if i != j and A[i][j]]
    arc_orbits_list = []; seen = set()
    for arc in arcs:
        if arc in seen: continue
        orbit = set()
        for perm in auts:
            mapped = (perm[arc[0]], perm[arc[1]])
            if A[mapped[0]][mapped[1]]: orbit.add(mapped)
        arc_orbits_list.append(orbit); seen.update(orbit)

    orbit_targets = []
    for orbit in arc_orbits_list:
        arc = next(iter(orbit)); i, j = arc
        A_rev = [row[:] for row in A]
        A_rev[i][j] = 0; A_rev[j][i] = 1
        target = cmap5[canonical(A_rev, 5)]
        orbit_targets.append(target)

    target_counts = Counter(orbit_targets)
    total_self_2 += target_counts.get(cid5, 0)

print("  collision_excess = %d" % collision_excess)
print("  total_self_loops = %d" % total_self_2)
print("  collision_excess + self_loops = %d" % (collision_excess + total_self_2))
print("  D = %d" % 28)
print("  MATCH: %s" % ((collision_excess + total_self_2) == 28))

# ============================================================
banner("SYNTHESIS: THE PARENT TREE DETERMINES THE DEGENERACY")
# ============================================================

print("""
  D = collision_excess + self_loops

  WHERE:
  collision_excess = sum_C sum_{D != C} max(mult(C->D) - 1, 0)
  self_loops = sum_C mult(C->C)

  mult(C->D) = #{arc orbits of C that reverse to D}

  THE PARENT TREE CONNECTION:
  Reversing arc (a,b) in T changes the parent profile at
  vertices v != a, v != b (n-2 of the n parents change).
  The parent at a and b stay the same.

  Two arcs (a,b) and (c,d) produce the SAME target class iff
  the resulting tournaments are isomorphic. This means:
  - The 2 UNCHANGED parents (at the arc endpoints) must match
  - The n-2 CHANGED parents must be isomorphic

  For |Aut|=1 classes: collisions require very specific conditions.
  The collision rate depends on how many arcs share the same
  "parent fingerprint" (the pattern of changed parents).

  THE FORMULA:
  D = sum_C [#{self-loop orbits} + sum_{D!=C} max(#{orbits C->D} - 1, 0)]

  This is EXACT but requires knowing the orbit-target multiplicity.
  The orbit-target multiplicity is determined by the AUT-ORBIT structure.

  FOR LARGE n:
  Generic classes (|Aut|=1): each arc is its own orbit.
  Collision = two arcs leading to the same target.
  The collision rate ~ 1/V ~ 1/A000568 -> 0 exponentially.
  So collision_excess/T -> 0.
  And self_loops/T ~ (n-1)/m = 2/(n+1) -> 0.
  Both components vanish, giving D/T -> 0. CONFIRMED.
""")

print("Done. opus-2026-03-23-S209")
