#!/usr/bin/env python3
"""
Combinatorial investigation: WHY is β_m^{(0)} = m(m-3)/2?

m(m-3)/2 appears in several combinatorial contexts:
1. C(m-1,2) - 1 = (m-1)(m-2)/2 - 1 — nope, that's (m²-3m)/2 = m(m-3)/2 ✓
2. C(m,2) - m = m(m-1)/2 - m = m(m-3)/2 ✓
3. Number of NON-ADJACENT pairs in a cycle graph C_m: m(m-3)/2 ✓
4. Dimension of the space of diagonals in a convex m-gon: m(m-3)/2 ✓
5. Number of edges in the complement of C_m: C(m,2) - m = m(m-3)/2 ✓

The match with "diagonals of an m-gon" is suggestive!

Z_m acts on QR by multiplication. QR has a natural cyclic structure.
Are the β_m generators indexed by "diagonals" (non-adjacent pairs) of QR?

For m=3: m(m-3)/2 = 0. A triangle has 0 diagonals. β_m = 0. ✓
For m=5: m(m-3)/2 = 5. A pentagon has 5 diagonals. β_m = 5. ✓
For m=7: m(m-3)/2 = 14. A heptagon has 14 diagonals. β_m = 14. ✓ (predicted)
For m=9: m(m-3)/2 = 27. A nonagon has 27 diagonals. β_m = 27. ✓ (predicted)

The ORBIT version: β_m^{orb} = (m-3)/2.
In the m-gon, Z_m acts by rotation. Each diagonal orbit has exactly m diagonals
(for prime m, the action is free on diagonals since any diagonal has a unique "gap size").
Number of orbits = m(m-3)/2 / m = (m-3)/2. ✓

Each diagonal corresponds to a "gap" g ∈ {2, 3, ..., m-2} (where gap 1 and m-1 are edges).
Number of non-edge gaps = m - 3. Orbits under Z_m: (m-3) gaps, but each gap g and m-g
give the same diagonal. So (m-3)/2 orbits for odd m. ✓

THIS IS THE KEY INSIGHT: β_m^{orb} counts diagonal TYPES in a regular m-gon!

Can we make this precise? Is there a bijection between:
  - Independent m-cycles in the orbit complex (non-boundary cycles in Ω_m/Z_m)
  - Diagonal types {g, m-g} for g = 2, ..., (m-1)/2 in the regular m-gon?

For m=5: (m-3)/2 = 1 diagonal type (g=2, which equals m-2=3 by symmetry).
For m=9: (m-3)/2 = 3 diagonal types (g=2, g=3, g=4, i.e., {2,7}, {3,6}, {4,5}).

opus-2026-03-13-S71b
"""

# Let's verify this connection more carefully
for m in range(3, 20, 2):
    p = 2*m + 1
    beta_m = m*(m-3)//2
    beta_orb = (m-3)//2
    diag = m*(m-3)//2  # diagonals of m-gon
    diag_types = (m-3)//2  # diagonal types (gap orbits under rotation)
    print(f"m={m}: β_m={beta_m}, β_orb={beta_orb}, "
          f"diag(m-gon)={diag}, diag_types={diag_types}")

print()

# For the m-gon: vertices are {0, 1, ..., m-1} = Z_m.
# Edges connect i to i+1 (and i to i-1). So edge gaps are {1, m-1}.
# Diagonals: gaps {2, 3, ..., m-2}.
# Gap g and gap m-g give the same "unoriented diagonal type".
# For odd m: all gaps are distinct under g ↔ m-g (since m is odd, g ≠ m-g for any g).
# So diagonal types = (m-3)/2 (the number of gap pairs {g, m-g} with 2 ≤ g ≤ (m-1)/2).

print("Diagonal types for small m:")
for m in range(3, 14, 2):
    gaps = list(range(2, (m+1)//2))
    print(f"  m={m}: gaps {gaps}, types {len(gaps)}, (m-3)/2={len(gaps)}")

print()

# Now: in the Paley tournament context, QR = {quadratic residues mod p}.
# Z_m = QR acts on QR by multiplication. The "vertices" of our "polygon" are
# elements of QR: {q_0, q_1, ..., q_{m-1}}.
#
# Under the Z_m action, the "graph" structure on QR is:
#   q_i → q_{i+1} if q_{i+1}/q_i ∈ QR (multiplication by a QR element).
#   This gives a COMPLETE graph on QR (since QR·QR = QR).
#
# Wait, that means every pair of QR elements is connected by a Z_m-element.
# So the "adjacency" is trivial — every pair is equivalent under Z_m.
# There's no natural "cycle graph" structure.

# BUT: the DIFF-SEQ complex has a natural notion of "adjacency":
# Two QR elements s, t are "adjacent" in the diff-seq sense if they can
# appear as consecutive entries in an allowed diff-seq.
# (s, t) ∈ A_2 iff partial sum s+t mod p is not 0 AND (s+t) not ∈ visited.
# But for just degree 2: (s, t) ∈ A_2 iff s+t ≢ 0 mod p.
# Since s, t ∈ QR, s+t = 0 iff t = -s. For p ≡ 3 mod 4, -1 ∉ QR,
# so -s ∉ QR for s ∈ QR. So s + t ≢ 0 for s, t ∈ QR.
# Thus ALL pairs (s, t) ∈ QR × QR are in A_2!
# |A_2| = m² (confirmed).
#
# But Ω_2 < A_2 because of the constraint: the "junk face"
# face_1(s, t) = (s+t) must be in QR for the chain to be in Omega.
# If s + t ∉ QR, then face_1 is "junk" and contributes to the constraint.
#
# An element (s, t) ∈ Ω_2 iff s + t ∈ QR (so that all faces are allowed).
# Wait, the constraint says chains whose boundary has junk components are zero.
# But (s, t) has faces (t), (s+t), (s). Face 0 = (t) ∈ A_1 always.
# Face 2 = (s) ∈ A_1 always. Face 1 = (s+t).
# If s+t ∉ QR, then face_1 is NOT in A_1, so it's "junk".
# ∂(s,t) = (t) - (s+t) + (s), where (s+t) might not be in A_1.
# The constraint: chains in Ω_2 have no junk face contribution.
# For individual diff-seqs (s,t), the junk face is (s+t) when s+t ∉ QR.
# The constraint row for junk face (s+t): coefficient of (s,t) is (-1)^1 = -1.
# So the constraint says: for each non-QR value v=s+t, Σ_{(s,t): s+t=v} (-1) · c_{(s,t)} = 0.
# i.e., c_{(s,t)} = 0 for each (s,t) with s+t ∉ QR (since each pair maps to a unique v).
# Wait, but multiple (s,t) could give the same v! s+t=v means t=v-s.
# For fixed v ∉ QR: t = v-s. We need s ∈ QR and t = v-s ∈ QR.
# Number of such pairs = #{s ∈ QR : v-s ∈ QR} = N_QQ(v).
# For v ∉ QR (nonzero): N_QQ(v) = (m+1)/2 by Jacobi sums.

# So the constraint for junk face v says: Σ_{s: v-s ∈ QR, s ∈ QR} c_{(s,v-s)} = 0.
# This relates (m+1)/2 different diff-seqs, not just one!
# So Ω_2 ≠ {(s,t) : s+t ∈ QR}. It's more subtle — linear combinations can cancel.

# For Ω_2^orb = m-1: the constraint kills 1 orbit dimension at degree 2.
# There are m orbits (since |A_2|/m = m), and rank(C) = 1 orbit.
# So Ω_2 = m(m-1), Ω_2^orb = m-1.

# The key question: is there a "polygon-like" structure in the Paley orbit
# complex that explains β_m^{orb} = "diagonals"?

# Perhaps the orbit complex is homotopy equivalent to the cycle graph C_m
# with some suspension? The cycle graph has β_1 = 1. An (m-1)-fold suspension
# would give β_m = 1. But β_m^{orb} = (m-3)/2, not 1.

# Alternatively: the orbit complex might be homotopy equivalent to a
# wedge of (m-3)/2 copies of S^m ∨ S^{m+1}.

print("\nSUMMARY:")
print("β_m = m(m-3)/2 = diagonals of a regular m-gon")
print("β_m^{orb} = (m-3)/2 = diagonal TYPES (gap orbits) of m-gon")
print()
print("This suggests a bijection between:")
print("  - Orbit m-cycles in H_m of the k=0 chain complex")
print("  - Diagonal types {g, m-g} for g = 2,...,(m-1)/2")
print()
print("Each diagonal type corresponds to one independent cycle,")
print("and the Z_m orbit of that cycle gives m copies (one per QR element).")
print()
print("OPEN: Is there an explicit construction mapping gap g to an m-cycle?")
