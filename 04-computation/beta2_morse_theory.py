#!/usr/bin/env python3
"""
beta2_morse_theory.py - Discrete Morse theory approach to beta_2 = 0

Idea: Find an acyclic matching on Omega_2(T) that pairs generators
with either Omega_1 or Omega_3 generators, leaving no unpaired 2-cells.
If successful, this proves beta_2 = 0 constructively.

For the GLMY path complex:
- Omega_0: vertices [n of them]
- Omega_1: edges a->b [n(n-1)/2 for tournaments, i.e. all edges]
- Omega_2: allowed 2-paths in Omega
- Omega_3: allowed 3-paths in Omega

An acyclic matching on the chain complex is a partial matching of
generators in adjacent dimensions with d-compatibility.

Alternative approach: Find the KERNEL of d_2|Omega_2 explicitly and
show it equals the IMAGE of d_3|Omega_3 by constructing explicit
preimages.

Actually, let me try the simplest possible approach first:
Is there a TOPOLOGICAL ORDERING of the tournament vertices such that
the resulting "decreasing 2-paths" form a basis for ker(d_2) that
is visibly in im(d_3)?

At n=5, the transitive tournament has score (0,1,2,3,4) and beta_1=0.
In this case, ker(d_2) has dimension 4 and im(d_3) has dimension 4.
The Omega complex for the transitive tournament IS the simplex complex
(all faces are transitive), so this is just simplicial homology of
the 4-simplex, which is trivially acyclic.

For non-transitive tournaments, the Omega complex is more complex.
But maybe we can find a "shelling order" or "collapsing sequence."

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


# ============================================================
# PART 1: What is the structure of ker(d_2) for different tournaments?
# ============================================================
print("=" * 70)
print("KERNEL STRUCTURE OF d_2|Omega_2")
print("=" * 70)

n = 5

# For each tournament, find ker(d_2|Omega_2) and express generators
# in terms of "elementary swap cycles" or other basic elements.

# A key question: does ker(d_2) have a basis of elements supported
# on a SINGLE vertex triple {i,j,k}?

print(f"\nn=5: Support analysis of ker(d_2) generators")

single_triple = 0
multi_triple = 0
total_gens = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        continue

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2
    U, S, Vt = np.linalg.svd(D2_omega, full_matrices=True)
    rank = sum(s > 1e-8 for s in S)
    ker_dim = dim_O2 - rank

    if ker_dim == 0:
        continue

    ker_basis = Vt[rank:]  # rows = ker vectors in Omega_2 coords
    ker_A2 = omega2 @ ker_basis.T  # columns in A_2 coords

    for g in range(ker_A2.shape[1]):
        v = ker_A2[:, g]
        total_gens += 1

        # Find which triples are involved
        triples = set()
        for j, path in enumerate(paths2):
            if abs(v[j]) > 1e-10:
                triples.add(tuple(sorted(path)))

        if len(triples) == 1:
            single_triple += 1
        else:
            multi_triple += 1

print(f"  Total ker(d_2) generators: {total_gens}")
print(f"  Supported on single triple: {single_triple}")
print(f"  Supported on multiple triples: {multi_triple}")


# ============================================================
# PART 2: Can we find a NICE basis for ker(d_2)?
# For 3-cycle triples {i,j,k} with i->j->k->i:
#   The 3 allowed 2-paths are: (i,j,k), (j,k,i), (k,i,j)
#   with coefficients in Omega_2 constrained by the non-allowed faces.
#   If all 3 are unconstrained, the "rotation element"
#   r_{ijk} = (i,j,k) + (j,k,i) + (k,i,j) might be in ker(d_2).
# ============================================================
print(f"\n{'='*70}")
print("3-CYCLE ROTATION ELEMENTS IN ker(d_2)")
print("=" * 70)

n = 5
rotation_in_ker = 0
rotation_total = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    path2_idx = {p: i for i, p in enumerate(paths2)}

    D2 = build_full_boundary_matrix(paths2, paths1)

    # Find 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 3 directed orientations
                if A[i][j] and A[j][k] and A[k][i]:
                    # 3-cycle i->j->k->i
                    p1, p2, p3 = (i,j,k), (j,k,i), (k,i,j)
                elif A[j][i] and A[i][k] and A[k][j]:
                    # 3-cycle j->i->k->j
                    p1, p2, p3 = (j,i,k), (i,k,j), (k,j,i)
                else:
                    continue

                # These 3 paths should all be allowed
                if p1 not in path2_idx or p2 not in path2_idx or p3 not in path2_idx:
                    continue

                rotation_total += 1

                # Build the rotation element
                r = np.zeros(len(paths2))
                r[path2_idx[p1]] = 1
                r[path2_idx[p2]] = 1
                r[path2_idx[p3]] = 1

                # Check if d_2(r) = 0
                dr = D2 @ r
                if np.max(np.abs(dr)) < 1e-8:
                    rotation_in_ker += 1
                else:
                    # How about with alternating signs?
                    r_alt = np.zeros(len(paths2))
                    r_alt[path2_idx[p1]] = 1
                    r_alt[path2_idx[p2]] = -1
                    r_alt[path2_idx[p3]] = 1
                    dr_alt = D2 @ r_alt
                    if np.max(np.abs(dr_alt)) < 1e-8:
                        rotation_in_ker += 1

print(f"n=5: {rotation_total} 3-cycles, {rotation_in_ker} have rotation/alternating in ker(d_2)")


# ============================================================
# PART 3: The d_2 map at the level of individual 2-paths
# ============================================================
print(f"\n{'='*70}")
print("BOUNDARY MAP d_2 STRUCTURE")
print("=" * 70)

# For an allowed 2-path (a,b,c):
# d_2(a,b,c) = (b,c) - (a,c) + (a,b)
#
# In the tournament:
# (a,b): ALWAYS allowed (either a->b or b->a in T, and (a,b) means a->b since it's a 1-path)
# (b,c): ALWAYS allowed
# (a,c): Allowed iff a->c (transitive triple). Non-allowed iff c->a (3-cycle).
#
# For transitive triple (TT) (a,b,c) with a->b->c AND a->c:
#   d_2(a,b,c) = (b,c) - (a,c) + (a,b) -- all three faces are allowed
#
# For 3-cycle triple (3C) (a,b,c) with a->b->c AND c->a:
#   d_2(a,b,c) = (b,c) - [non-allowed (a,c)] + (a,b)
#   The non-allowed face (a,c) doesn't exist in Omega_1,
#   so in the Omega complex: d_2(a,b,c) = (b,c) + (a,b)
#   Wait, that's wrong. The boundary map goes A_2 -> A_1,
#   and THEN we project to Omega_1. But Omega_1 = A_1 for tournaments!
#   (All 1-paths are allowed in a tournament since every pair has an edge.)
#
# Actually, (a,c) as a 1-path means "the directed path from a to c."
# If c->a in T, then (a,c) is NOT an allowed 1-path.
# So d_2(a,b,c) in A_1 has (a,c) which is NOT in Omega_1.
# But... the Omega boundary map goes Omega_2 -> Omega_1,
# and we need d_2(Omega_2) subset Omega_1. This is automatic for
# elements of Omega_2 because the Omega constraint kills the
# non-allowed parts.

# Actually, I need to be more careful. Let me check: for a single
# 3-cycle path (a,b,c) with c->a: is (a,b,c) in Omega_2?

# The Omega_2 constraint: for each non-allowed 1-path (x,y),
# the sum of coefficients of all 2-paths with (x,y) as the
# deleted-middle face must vanish.
#
# For (a,b,c) with c->a: the middle face is (a,c) which IS non-allowed.
# Constraint: sum over all 2-paths with middle face (a,c) must be 0.
# These are paths (a,w,c) for w with a->w->c.
# If (a,b,c) is the ONLY such path, then its coefficient must be 0.
# If there are other mediators w1,w2,..., the sum must be 0.
#
# For a single (a,b,c): it's in Omega_2 ONLY IF there exists another
# 2-path (a,w,c) with a->w->c AND the sum constraint is satisfiable.

print(f"\nn=5: Analysis of 3-cycle vs TT paths in Omega_2")

for bits in [76]:  # regular n=5
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    print(f"  T#{bits} scores={scores}")

    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    print(f"  dim(A_2)={len(paths2)}, dim(Omega_2)={dim_O2}")

    # Group by triple
    by_triple = defaultdict(list)
    for j, (a,b,c) in enumerate(paths2):
        triple = tuple(sorted([a,b,c]))
        face_type = "TT" if A[a][c] else "3C"
        by_triple[triple].append(((a,b,c), face_type, j))

    print(f"  Triples ({len(by_triple)}):")
    for triple in sorted(by_triple.keys()):
        entries = by_triple[triple]
        types = [e[1] for e in entries]
        tt_count = types.count("TT")
        c3_count = types.count("3C")
        print(f"    {triple}: {len(entries)} paths ({tt_count} TT, {c3_count} 3C)")

    # For each non-allowed pair (a,c) with c->a, how many mediators?
    print(f"\n  Non-allowed pairs and mediators:")
    for a in range(n):
        for c in range(n):
            if a != c and A[a][c] == 0:  # c->a
                meds = [b for b in range(n) if b != a and b != c
                        and A[a][b] == 1 and A[b][c] == 1]
                print(f"    ({a},{c}) [c->a]: {len(meds)} mediators: {meds}")


# ============================================================
# PART 4: The Omega_2 constraint matrix
# ============================================================
print(f"\n{'='*70}")
print("OMEGA_2 CONSTRAINT MATRIX STRUCTURE")
print("=" * 70)

# The constraint matrix C maps A_2 -> R^{#non-allowed-pairs}
# where C[non-allowed (a,c), path (a,b,c)] = 1 if b is a mediator.
# C has exactly one row per non-allowed pair (a,c).
# C has exactly one column per 2-path, but only 3-cycle-type paths
# have nonzero entries (TT paths have all faces allowed).

# Actually: the constraint is for EACH non-allowed 1-path (x,y):
# sum over all 2-paths in A_2 whose "middle face" is (x,y) = 0.
# The middle face of (a,b,c) is (a,c) (delete middle vertex b).
# So the constraint matrix has:
# Row for non-allowed (a,c): column (a,b,c) gets +1 for each b mediator.

# For TOURNAMENT: non-allowed (a,c) means c->a.
# Mediators: b with a->b AND b->c. (a->b->c with c->a = 3-cycle on {a,b,c})

# KEY INSIGHT: Each 3-cycle triple {a,b,c} with i->j->k->i
# contributes constraints at THREE non-allowed pairs:
# (j,i) [i->j, constraint: sum of paths (j,?,i) = 0]
# (k,j) [j->k, constraint: sum of paths (k,?,j) = 0]
# (i,k) [k->i, constraint: sum of paths (i,?,k) = 0]
#
# But EACH non-allowed pair can have mediators from MULTIPLE triples!
# E.g., (a,c) with c->a: mediators can be in {a,b1,c}, {a,b2,c}, etc.

# The constraint structure is a BIPARTITE GRAPH:
# - One side: non-allowed pairs (a,c)
# - Other side: 3-cycle-type 2-paths (a,b,c)
# - Edge: (a,c) -- (a,b,c) if b is a mediator

# For the constraint matrix to have rank = #{non-allowed pairs with mediators},
# we need each row to have at least one nonzero entry (obvious for pairs with mediators),
# AND the rows to be linearly independent.
# This is guaranteed because different non-allowed pairs (a,c) involve
# DISJOINT sets of 2-paths (different starting vertex a or ending vertex c).

# Wait, is that true? (a,c) involves paths (a,*,c) and (a',c') involves (a',*,c').
# These are disjoint iff (a,c) != (a',c'), which is true for different pairs.
# YES! Each 2-path (a,b,c) appears in EXACTLY ONE constraint row (for the pair (a,c)).

# THIS IS KEY: the constraint matrix has orthogonal rows!
# Therefore: dim(Omega_2) = dim(A_2) - #{non-allowed pairs with at least 1 mediator}
# EXACTLY.

print(f"Constraint matrix has orthogonal rows (each 2-path in exactly 1 constraint).")
print(f"This is because each 2-path (a,b,c) has a UNIQUE middle face (a,c).")
print(f"Therefore:")
print(f"  dim(Omega_2) = dim(A_2) - #{'{'}non-allowed pairs with mediators{'}'}")
print(f"  = C(n,3) + 2*c3 - #{'{'}non-allowed pairs with mediators{'}'}")

# Now: what IS #{non-allowed pairs with mediators}?
# Non-allowed (a,c) means c->a. Mediator b means a->b AND b->c.
# So: #{a,c: c->a AND exists b with a->b->c}
# = #{a,c: c->a AND a,c are NOT the only path from a to c in T}

# For TOURNAMENTS: a->c does NOT hold (c->a), so is there a path a ~> c?
# a ~> c of length 2: exists b with a->b->c. This is what we need.
# If a->b->c for some b, then (a,c) has a mediator.
# If no such b exists: every vertex b has either b->a or c->b.
# This means: N+(a) subset N+(c) union {c}? No...
# a->b but NOT b->c means c->b. And b->a means b is in N-(a).
# So: for every b != a,c: either b->a OR c->b (or both).
# In other words: N+(a) \ {c} = emptyset (everything that a beats... nothing except c, but c->a).
# Wait, a doesn't beat c. So N+(a) \ {c} = N+(a).
# We need: for all b in N+(a), c->b. This means N+(a) subset N-(c).
# Since |N+(a)| = out-degree of a, and |N-(c)| = in-degree of c:
# N+(a) subset N-(c) means out(a) <= in(c) AND all of a's successors are c's predecessors.

# This is VERY restrictive. It means a is "dominated by" c in a strong sense.
# For regular tournaments (all degrees (n-1)/2), this is impossible unless n is very small.

# So: #{non-allowed pairs WITHOUT mediators} should be small.
# Let's count.

print(f"\nn=5: Non-allowed pairs without mediators:")
count_no_med = 0
count_with_med = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    no_med = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[a][c] == 0:  # c->a
                has_med = any(A[a][b] == 1 and A[b][c] == 1
                             for b in range(n) if b != a and b != c)
                if not has_med:
                    no_med += 1

    if no_med > 0:
        count_no_med += 1
    else:
        count_with_med += 1

print(f"  Tournaments with ALL non-allowed pairs having mediators: {count_with_med}/1024")
print(f"  Tournaments with SOME no-mediator pairs: {count_no_med}/1024")

# Check at n=6
n = 6
count_no_med_6 = 0
count_with_med_6 = 0
max_no_med = 0

for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    no_med = 0
    for a in range(n):
        for c in range(n):
            if a != c and A[a][c] == 0:
                has_med = any(A[a][b] == 1 and A[b][c] == 1
                             for b in range(n) if b != a and b != c)
                if not has_med:
                    no_med += 1

    if no_med > 0:
        count_no_med_6 += 1
        max_no_med = max(max_no_med, no_med)
    else:
        count_with_med_6 += 1

print(f"\nn=6:")
print(f"  ALL pairs have mediators: {count_with_med_6}/32768")
print(f"  SOME no-mediator pairs: {count_no_med_6}/32768")
print(f"  Max no-mediator count: {max_no_med}")


# ============================================================
# PART 5: Is there a topological sort - based collapsing?
# ============================================================
print(f"\n{'='*70}")
print("TOPOLOGICAL SORT BASED ANALYSIS")
print("=" * 70)

# For a transitive tournament, the topological sort is unique.
# For a general tournament, we can still define a "topological sort"
# using the score sequence: order vertices by out-degree.
# This gives a linear order v_1 < v_2 < ... < v_n where
# "most" edges go from high to low (like a transitive tournament).
# The non-transitive edges are exactly the 3-cycles.

# Idea: Filter Omega_p by the "topological defect" = number of
# 3-cycle-type faces. Transitive paths have defect 0.
# Show that defect > 0 paths can be collapsed.

# First: what's the defect distribution?
n = 5
for bits in [0, 76, 148]:  # transitive, regular, near-regular
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    paths2 = enumerate_allowed_paths(A, n, 2)
    defect_counts = defaultdict(int)

    for (a,b,c) in paths2:
        defect = 0
        if A[a][c] == 0: defect += 1  # middle face (a,c) is non-allowed
        defect_counts[defect] += 1

    print(f"\n  T#{bits} scores={scores}: defect distribution: {dict(defect_counts)}")
    betti = path_betti_numbers(A, n, max_dim=3)
    print(f"    beta={list(betti)}")


print("\n\nDone.")
