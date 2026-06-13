#!/usr/bin/env python3
"""Approach β_2=0 via the junk matrix rank structure.

The chain complex is:
  ... → A_3 →^{∂_3} A_2 →^{∂_2} A_1 → ...

The GLMY path complex truncates to:
  ... → Ω_3 →^{∂_3} Ω_2 →^{∂_2} Ω_1 → ...

where Ω_p = ker(J_p) and J_p is the "junk matrix" that maps
p-chains to their non-allowed faces.

KEY STRUCTURAL PROPERTY FOR TOURNAMENTS:
- J_2: each 2-path (a,b,c) has AT MOST one bad face (a,c), iff c→a
- Bad faces in DISJOINT groups (each bad pair indexes a group)
- J_2 is block-diagonal with all-ones blocks

FOR β_2=0: we need rank(∂_3|_{Ω_3}) = dim(ker(∂_2|_{Ω_2})).

Alternative approach via the EXACT SEQUENCE of the junk quotient:
  0 → Ω_2 →^{inc} A_2 →^{J_2} J_1
where J_1 = "junk 1-path space" = span of non-allowed faces.

The boundary ∂_2 on A_2 sends to R_1 (all 1-chains), not just A_1.
But ∂_2(Ω_2) ⊂ A_1 by definition.

Can we use the SNAKE LEMMA or a diagram chase to relate
the homology of the full complex (A_*) to that of the truncated (Ω_*)?
"""
import numpy as np
from itertools import combinations
import sys
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Approach 1: Full A_* complex homology =====
print("=" * 70)
print("FULL A_* COMPLEX HOMOLOGY")
print("=" * 70)

# In the full complex A_*, the boundary ∂_p: A_p → R_{p-1}
# where R_{p-1} = ALL (p-1)-paths (not just allowed ones).
# But we can restrict to ∂_p: A_p → A_{p-1} ⊕ J_{p-1}
# where J_{p-1} is the junk space.

# Actually, ∂_p maps A_p to the span of all (p-1)-paths.
# Each face of an allowed p-path is either:
#   - An allowed (p-1)-path (in A_{p-1}), or
#   - A non-allowed (p-1)-path (a "junk" face)

# Define: J_p = set of non-allowed p-paths that appear as faces of A_{p+1} paths.
# Then ∂_{p+1}: A_{p+1} → A_p ⊕ J_p

# The GLMY path complex Ω_{p+1} is precisely the kernel of the projection
# A_{p+1} → J_p.

# So we have the short exact sequence:
# 0 → Ω_{p+1} → A_{p+1} →^{π} J_p
# where π sends each (p+1)-path to the "junk part" of its boundary.

# The boundary ∂_{p+1} = ∂_{p+1}^A + ∂_{p+1}^J where
#   ∂_{p+1}^A: A_{p+1} → A_p (allowed part of boundary)
#   ∂_{p+1}^J: A_{p+1} → J_p (junk part of boundary)
# And Ω_{p+1} = ker(∂_{p+1}^J).

# Now, ∂ ∘ ∂ = 0 in the FULL chain complex R_*.
# For u ∈ A_{p+1}: ∂_{full}(∂_{full}(u)) = 0, where ∂_{full} goes to R_{p-1}.
# The boundary of the junk faces ∂(J_p) goes to R_{p-1} as well.

n = 5
print(f"\nn={n}: Checking full complex structure")

for tidx, A in enumerate(all_tournaments_gen(n)):
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    if scores != (2,2,2,2,2):
        continue

    # Enumerate all paths at each level
    a = {}
    for p in range(5):
        a[p] = [tuple(x) for x in enumerate_allowed_paths(A, n, p)]

    # Junk faces at level 1: non-allowed 1-paths (non-edges)
    # For tournaments: ALL pairs have edges, so J_1 = ∅!
    j1 = [(i,j) for i in range(n) for j in range(n) if i != j and not A[i][j]]
    print(f"  J_1 (non-edges): {len(j1)}")
    # For tournaments: j1 = {(a,c): c→a}, which are the REVERSE edges.
    # These are NOT "non-edges" — they're edges in the other direction.
    # Actually, (a,c) is a 1-path. It's ALLOWED if a→c. If c→a, then (a,c) is NOT allowed.
    # So j1 = {(a,c): c→a} = reverse edges. |j1| = C(n,2).
    print(f"  Actually for tournaments: J_1 = reverse edges = {len(j1)}")

    # Junk faces at level 2: non-A_2 2-paths that appear as faces of A_3 paths
    all_2_paths = set()
    for p3 in a[3]:
        a_, b_, c_, d_ = p3
        faces = [(b_,c_,d_), (a_,c_,d_), (a_,b_,d_), (a_,b_,c_)]
        for f in faces:
            if f not in set(a[2]):
                all_2_paths.add(f)
    j2 = list(all_2_paths)
    print(f"  J_2 (non-A_2 faces of A_3): {len(j2)}")

    # What makes a 2-path non-allowed?
    # (x,y,z) is allowed iff x→y→z and x,y,z distinct
    # Non-allowed from A_3 faces:
    for f in sorted(j2)[:5]:
        x,y,z = f
        xy = A[x][y]
        yz = A[y][z]
        distinct = len({x,y,z}) == 3
        print(f"    {f}: x→y={xy}, y→z={yz}, distinct={distinct}")

    # For tournaments: (x,y,z) with x,y,z distinct is in A_2 iff x→y AND y→z.
    # Non-allowed: x→y but z→y (wrong direction at y→z), or y→x (wrong at x→y).
    # From faces of A_3 path (a,b,c,d):
    #   (a,c,d): allowed iff a→c. If c→a: not allowed.
    #   (a,b,d): allowed iff b→d. If d→b: not allowed.
    # (b,c,d) and (a,b,c) are always allowed (edges a→b→c→d given).

    # So J_2 faces come in two types:
    # Type A: (a,c,d) with c→a (from face at position 1)
    # Type B: (a,b,d) with d→b (from face at position 2)
    type_a = [(a_,c_,d_) for (a_,b_,c_,d_) in a[3] if A[c_][a_]]
    type_b = [(a_,b_,d_) for (a_,b_,c_,d_) in a[3] if A[d_][b_]]
    type_a_set = set(type_a)
    type_b_set = set(type_b)
    print(f"\n  J_2 type A (c→a at pos 1): {len(type_a_set)}")
    print(f"  J_2 type B (d→b at pos 2): {len(type_b_set)}")
    print(f"  Overlap: {len(type_a_set & type_b_set)}")

    # For each junk face, how many A_3 paths produce it?
    junk_mult_a = Counter()
    for (a_,b_,c_,d_) in a[3]:
        if A[c_][a_]:
            junk_mult_a[(a_,c_,d_)] += 1
    junk_mult_b = Counter()
    for (a_,b_,c_,d_) in a[3]:
        if A[d_][b_]:
            junk_mult_b[(a_,b_,d_)] += 1

    print(f"\n  Type A multiplicities: {dict(Counter(junk_mult_a.values()))}")
    print(f"  Type B multiplicities: {dict(Counter(junk_mult_b.values()))}")

    # KEY INSIGHT: each Type A face (a,c,d) comes from deleting b in (a,b,c,d).
    # The intermediary b satisfies: a→b→c (and c→a), b→d... wait, no.
    # b appears in (a,b,c,d) so a→b→c→d and all distinct.
    # For face (a,c,d) to be junk: c→a.
    # The intermediary b: a→b, b→c, b ∉ {a,c,d}, and c→d.
    # So b ∈ N^+(a) ∩ N^-(c) \ {a,c,d}.

    # Each Type A junk face (a,c,d) has mult = |N^+(a) ∩ N^-(c) \ {a,c,d}|.
    # Similarly for Type B.

    break

# ===== Approach 2: The "junk-face acyclicity" condition =====
print(f"\n\n{'='*70}")
print("JUNK FACE BOUNDARY STRUCTURE")
print("=" * 70)

# For β_2=0, we need: ∂_3(Ω_3) = ker(∂_2|_{Ω_2}).
# Equivalently: there's no 2-cycle in Ω_2 that isn't a 3-boundary.
#
# Think of it as: the "extended boundary" ∂_3: A_3 → A_2 ⊕ J_2
# has the property that Ω_3 = ker(π_{J_2} ∘ ∂_3), where π_{J_2} is the
# projection onto J_2.
#
# The key diagram is:
#            ∂_3
#   A_3 ————————→ A_2 ⊕ J_2
#    ↑              ↑
#   Ω_3 ————————→ A_2   (∂_3 restricted, image in A_2 since Ω_3 = ker of J_2 part)
#                   ↓ ∂_2
#                  A_1 ⊕ J_1

# The composition ∂_2 ∘ ∂_3 = 0 in R_1 (all 1-paths), but
# ∂_2(A_2) may have J_1 components.
# For Ω_2: ∂_2(Ω_2) ⊂ A_1 (by definition).

# So the chain complex Ω_* has ∂∂=0 because the full R_* complex does.

# Let's check: is the "A_* complex" (with boundary restricted to A_*) also exact?
print(f"\nn={n}: H_2 of A_* complex (boundary mod J)")

for tidx, A in enumerate(all_tournaments_gen(n)):
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    if tidx > 500:
        break

    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    # ∂_2^{A→A}: A_2 → A_1 (only the allowed part of boundary)
    bd2_full = build_full_boundary_matrix(a2, a1)
    # This already only has A_1 columns, so it IS the A→A part.

    # ∂_3^{A→A}: A_3 → A_2 (only the allowed part of boundary)
    bd3_full = build_full_boundary_matrix(a3, a2)

    # H_2(A) = ker(∂_2^{A→A}) / im(∂_3^{A→A})
    rank_bd2 = np.linalg.matrix_rank(bd2_full, tol=1e-8)
    ker_bd2 = len(a2) - rank_bd2

    rank_bd3 = np.linalg.matrix_rank(bd3_full, tol=1e-8)

    h2_A = ker_bd2 - rank_bd3

    if h2_A != 0:
        print(f"  Tournament #{tidx}, scores={scores}: H_2(A)={h2_A}")
        break
else:
    print(f"  H_2(A) = 0 for all checked (first 500)")

# Check ALL
h2_A_dist = Counter()
for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    bd2_full = build_full_boundary_matrix(a2, a1)
    bd3_full = build_full_boundary_matrix(a3, a2)

    rank_bd2 = np.linalg.matrix_rank(bd2_full, tol=1e-8)
    ker_bd2 = len(a2) - rank_bd2
    rank_bd3 = np.linalg.matrix_rank(bd3_full, tol=1e-8)
    h2_A = ker_bd2 - rank_bd3

    h2_A_dist[h2_A] += 1

print(f"\n  n={n}: H_2(A_*) distribution: {dict(sorted(h2_A_dist.items()))}")
print(f"  H_2(A_*) = 0 for all? {'YES' if set(h2_A_dist.keys()) == {0} else 'NO'}")

# ===== CHECK n=4 =====
n = 4
h2_A_dist_4 = Counter()
for A in all_tournaments_gen(n):
    a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
    a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
    a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

    bd2_full = build_full_boundary_matrix(a2, a1)
    bd3_full = build_full_boundary_matrix(a3, a2)

    rank_bd2 = np.linalg.matrix_rank(bd2_full, tol=1e-8)
    ker_bd2 = len(a2) - rank_bd2
    rank_bd3 = np.linalg.matrix_rank(bd3_full, tol=1e-8)
    h2_A = ker_bd2 - rank_bd3

    h2_A_dist_4[h2_A] += 1

print(f"  n={n}: H_2(A_*) distribution: {dict(sorted(h2_A_dist_4.items()))}")

# ===== KEY QUESTION: Is H_2(A_*) = 0 stronger than H_2(Ω_*) = 0? =====
# If H_2(A_*) = 0, then for any z ∈ ker(∂_2|_{A_2}), there exists w ∈ A_3
# with ∂_3(w) = z (in A_2). But w may not be in Ω_3!
# We need w ∈ Ω_3, i.e., ∂_3(w) has no J_2 component.
#
# But if z ∈ Ω_2 ⊂ A_2 and ∂_3(w) = z for some w ∈ A_3,
# then does w ∈ Ω_3? Not necessarily — the J_2 component of ∂_3(w)
# might cancel in A_2 but w itself might have nontrivial J_2 projection.
#
# Wait, if ∂_3(w) = z ∈ A_2, that means the J_2 part of ∂_3(w) is zero.
# So w IS in Ω_3! Because Ω_3 = {w ∈ A_3 : J_2-part of ∂_3(w) = 0}.

# WAIT — this might be the proof!
# If H_2(A_*) = 0, then for any z ∈ ker(∂_2|_{A_2}), ∃w ∈ A_3: ∂_3^{A→A}(w) = z.
# This means the A_2 part of ∂_3(w) equals z and the J_2 part of ∂_3(w) = 0.
# So w ∈ Ω_3 and ∂_3(w) = z.
# In particular, if z ∈ Ω_2 ⊂ ker(∂_2|_{A_2}), then z = ∂_3(w) for some w ∈ Ω_3.
# This gives β_2(Ω) = 0!

# BUT WAIT: build_full_boundary_matrix(a3, a2) only includes the A_2 components
# of ∂_3, i.e., it projects out the J_2 part. Let me verify this.
# If ∂_3(w) ∈ A_2 (i.e., the full boundary has no non-A_2 components),
# then w ∈ Ω_3. And conversely.

# So: H_2(A_*) = 0 ⟹ β_2(Ω_*) = 0, where A_* has boundary projected to A_{*-1}.

# But IS H_2(A_*) = 0? Let me verify more carefully.

print(f"\n\n{'='*70}")
print("VERIFICATION: H_2(A_*) = 0 IMPLIES β_2(Ω_*) = 0")
print("=" * 70)

# Careful: the boundary matrix bd3 = build_full_boundary_matrix(a3, a2)
# gives the PROJECTION of ∂_3 onto A_2. This means:
# bd3[i,j] = coefficient of a2[i] in ∂_3(a3[j]) IF a2[i] appears as a face.
# But there might be ADDITIONAL non-A_2 faces that are dropped.

# H_2(A_*) = ker(bd2: A_2 → A_1) / im(bd3: A_3 → A_2)
# where bd3 is the PROJECTION of ∂_3 onto A_2.

# If H_2(A_*) = 0: for any z ∈ A_2 with bd2·z = 0 (in A_1),
# there exists w ∈ A_3 with bd3·w = z.
# This means: the A_2-projection of ∂_3(w) equals z.
# But ∂_3(w) = (A_2 part) + (J_2 part), and bd3·w = z means A_2 part = z.
# So the J_2 part of ∂_3(w) is determined by w but NOT constrained to be 0.

# CORRECTION: bd3·w = z does NOT mean w ∈ Ω_3!
# The J_2 part of ∂_3(w) could be nonzero.

# So H_2(A_*) = 0 does NOT directly imply β_2 = 0.
# We would need w to ALSO satisfy J_2-part = 0.

# But maybe we can MODIFY w: w' = w + correction, where correction ∈ ker(bd3)
# but the J_2 part of ∂_3(correction) cancels the J_2 part of ∂_3(w)?

# This is possible if the J_2 part of ∂_3(ker(bd3 projection)) spans
# the J_2 part of ∂_3(A_3). But this is hard to verify in general.

# However, if H_2(A_*) = 0, then bd3 is surjective onto ker(bd2).
# And Ω_3 = {w ∈ A_3 : J_2(∂_3(w)) = 0}.
# The question is: can we find w ∈ Ω_3 with bd3·w = z?
# i.e., can we find w in the intersection of Ω_3 and the preimage of z under bd3?

# This is a system: bd3·w = z AND J·w = 0 (where J is the Ω_3 constraint).

print("The relationship between H_2(A_*) and β_2(Ω_*) is subtle.")
print("H_2(A_*)=0 is NECESSARY but not SUFFICIENT for β_2=0.")
print()

# Actually wait - let me re-examine. The boundary matrix bd3 only has
# the ALLOWED face coefficients. But build_full_boundary_matrix
# computes ALL boundary coefficients for faces that are in the reference list.
# Let me check what build_full_boundary_matrix actually does.

print("Checking build_full_boundary_matrix behavior...")

n = 4
A = [[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,0,0,0]]  # transitive
a3 = [tuple(x) for x in enumerate_allowed_paths(A, n, 3)]
a2 = [tuple(x) for x in enumerate_allowed_paths(A, n, 2)]
a1 = [tuple(x) for x in enumerate_allowed_paths(A, n, 1)]

print(f"  A_3 = {a3}")
print(f"  A_2 = {a2}")

bd3 = build_full_boundary_matrix(a3, a2)
print(f"  bd3 shape: {bd3.shape}")
print(f"  bd3 = {bd3}")

# ∂(0,1,2,3) = (1,2,3) - (0,2,3) + (0,1,3) - (0,1,2)
# All faces should be in A_2 for transitive tournament
print(f"  Expected: (1,2,3)-(0,2,3)+(0,1,3)-(0,1,2)")

# For non-transitive: check what happens
A2 = [[0,1,0,1],[0,0,1,0],[1,0,0,1],[0,1,0,0]]  # has 3-cycles
a3_2 = [tuple(x) for x in enumerate_allowed_paths(A2, n, 3)]
a2_2 = [tuple(x) for x in enumerate_allowed_paths(A2, n, 2)]
print(f"\n  Non-transitive: A_3 = {a3_2}")
print(f"  A_2 = {a2_2}")
if a3_2:
    bd3_2 = build_full_boundary_matrix(a3_2, a2_2)
    print(f"  bd3 shape: {bd3_2.shape}")
    for j, p3 in enumerate(a3_2):
        faces = [(a2_2[i], bd3_2[i,j]) for i in range(len(a2_2)) if abs(bd3_2[i,j]) > 1e-10]
        # What are the ACTUAL faces?
        a_,b_,c_,d_ = p3
        all_faces = [(b_,c_,d_), (a_,c_,d_), (a_,b_,d_), (a_,b_,c_)]
        signs = [1, -1, 1, -1]
        print(f"  ∂({p3}): in A_2: {faces}")
        for f, s in zip(all_faces, signs):
            in_a2 = f in a2_2
            print(f"    {'+' if s>0 else '-'}{f}: in_A2={in_a2}")

print("\nDone.")
