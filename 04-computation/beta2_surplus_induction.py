#!/usr/bin/env python3
"""
beta2_surplus_induction.py — Induction proof via surplus connectivity

THEOREM (Computational, n≤6):
  The arc-flip graph restricted to {T : surplus(T) ≥ 0} is connected.
  The transitive tournament has surplus = C(n-1,4) ≥ 0.
  Therefore surplus ≥ 0 for all tournaments, hence β₂ = 0.

PROOF STRATEGY:
  1. Start from transitive tournament (all i→j for i<j)
  2. Any tournament can be reached by a sequence of arc flips
  3. Show surplus stays ≥ 0 along some path

We already verified exhaustive connectivity at n=5 (1024/1024) and n=6 (32768/32768).

The MISSING PIECE: an algebraic proof that the restricted graph is connected.
Equivalently: for any tournament T with surplus > 0, there exists an arc
u→v such that flipping it doesn't make surplus negative, AND such that
the flip brings T closer to ANY target.

SIMPLER APPROACH: Show surplus ≥ 0 DIRECTLY (not via connectivity).

Key identities we know:
  dim(Ω₂) = |A₂| - |BF₂| where BF₂ = {(a,c) : c→a, ∃b with a→b→c}
  rk(∂₂) ≤ dim(Ω₂)
  surplus = dim(Ω₃) - (dim(Ω₂) - rk(∂₂))

Can we show dim(Ω₃) ≥ dim(Ω₂) - rk(∂₂)?

Actually, WAIT. β₂ = 0 means Z₂ = im(∂₃|Ω₃).
This is equivalent to: rk(∂₃|Ω₃) = dim(Z₂).
And surplus = dim(Ω₃) - dim(Z₂) = dim(ker ∂₃|Ω₃) = β₃ ≥ 0.

So SURPLUS = β₃ ≥ 0 IFF β₂ = 0!

This means proving surplus ≥ 0 IS proving β₂ = 0. It's not circular —
we need to show that the kernel of ∂₃ on Ω₃ is large enough.

Actually, I realize surplus = dim(Ω₃) - dim(Z₂) ≥ 0 is NOT the same as β₂=0.
surplus = dim(Ω₃) - dim(Z₂)
β₂ = dim(Z₂) - rk(∂₃) = dim(Z₂) - (dim(Ω₃) - dim(ker ∂₃)) = dim(Z₂) - dim(Ω₃) + β₃
So β₂ = -surplus + β₃.
Since β₂ ≥ 0 and β₃ ≥ 0: surplus = β₃ - β₂.
surplus ≥ 0 iff β₃ ≥ β₂.
surplus = 0 iff β₃ = β₂.

For β₂ = 0: need surplus = β₃ ≥ 0 (automatically true!).
WAIT: if β₂ = 0 then surplus = β₃ ≥ 0. And conversely if surplus ≥ 0 then β₃ ≥ β₂.
But surplus ≥ 0 doesn't FORCE β₂ = 0.

Hmm, but computationally surplus is ALWAYS ≥ 0 and β₂ is always 0.
The surplus being ≥ 0 is a NECESSARY condition for β₂=0 but not sufficient.

Let me reconsider. The ACTUAL proof needs rk(∂₃) ≥ dim(Z₂).
rk(∂₃) ≤ dim(Ω₃). So dim(Ω₃) ≥ dim(Z₂) is necessary.
But rk(∂₃) could be less than dim(Ω₃) (when β₃ > 0).

So β₂ = 0 iff rk(∂₃) = dim(Z₂) iff im(∂₃) = Z₂.

NEW APPROACH: Instead of surplus, study rk(∂₃) directly.
rk(∂₃) changes under arc flips. Can we track it?

Actually, let me think about what we KNOW:
- dim(Ω₃) ≥ dim(Z₂) (= surplus ≥ 0, verified n≤8)
- rk(∂₃) = dim(Z₂) (= β₂=0, verified n≤8)
- These together give: dim(ker ∂₃) = dim(Ω₃) - dim(Z₂) = surplus = β₃

So β₃ = surplus exactly. The surplus IS β₃.

Strategy: Prove β₂ = 0 by showing im(∂₃) ⊇ Z₂.
This requires an algebraic argument about the boundary map.

Let me try a COMPLETELY DIFFERENT approach: vertex induction.

H₂(T, T\\v) = 0 for all T, v (verified n≤9).
Long exact sequence gives:
  H₂(T\\v) → H₂(T) → H₂(T, T\\v) → H₁(T\\v) → H₁(T)
  0 → H₂(T) → 0 → H₁(T\\v) → H₁(T)

So H₂(T) = 0! IF we know H₂(T\\v) = 0 (induction) and H₂(T, T\\v) = 0.

The base case: n=3, β₂=0 trivially.
Inductive step: assume β₂(T')=0 for all tournaments on n-1 vertices.
For tournament T on n, pick any vertex v.
T\\v is a tournament on n-1, so β₂(T\\v) = 0 by induction.
H₂(T,T\\v) = 0 is what we need to prove.

This reduces the problem to: PROVE H₂(T,T\\v) = 0.

Let me study what H₂(T,T\\v) = 0 means structurally.

Author: opus-2026-03-08-S44
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def delete_vertex(A, n, v):
    """Delete vertex v from tournament, returning (n-1) × (n-1) matrix."""
    B = []
    for i in range(n):
        if i == v:
            continue
        row = []
        for j in range(n):
            if j == v:
                continue
            row.append(A[i][j])
        B.append(row)
    return B

print("=" * 70)
print("RELATIVE HOMOLOGY H₂(T, T\\v) STRUCTURE")
print("=" * 70)

# The relative chain complex:
# R_p = Ω_p(T) / Ω_p(T\\v)
# The relative homology H_p(T,T\\v) = ker(∂_p^rel) / im(∂_{p+1}^rel)
#
# R_p consists of Ω_p elements that use vertex v.
# More precisely: R_p = span of Ω_p(T) elements that have nonzero
# projection onto A_p paths through v.
#
# KEY: A p-path through v = (v₀,...,v_p) where v ∈ {v₀,...,v_p}.
# These paths exist in T but NOT in T\\v.
#
# For H₂(T,T\\v) = 0, we need: every relative 2-cycle is a relative boundary.

print("\n--- Relative complex dimensions for n=5 ---")
n = 5

rel_dims = defaultdict(lambda: Counter())
rel_details = []

for T_idx, A in enumerate(all_tournaments(n)):
    for v in range(n):
        B = delete_vertex(A, n, v)

        # Compute Ω_p(T) and Ω_p(T\\v) for p=1,2,3
        dims_T = {}
        dims_Tv = {}

        for p in range(4):
            ap_T = enumerate_allowed_paths(A, n, p)
            ap_prev_T = enumerate_allowed_paths(A, n, p-1) if p > 0 else []

            if ap_T and (p == 0 or ap_prev_T):
                if p == 0:
                    dims_T[p] = n
                else:
                    om = compute_omega_basis(A, n, p, ap_T, ap_prev_T)
                    dims_T[p] = om.shape[1] if om.ndim == 2 else 0
            else:
                dims_T[p] = n if p == 0 else 0

            ap_Tv = enumerate_allowed_paths(B, n-1, p)
            ap_prev_Tv = enumerate_allowed_paths(B, n-1, p-1) if p > 0 else []

            if ap_Tv and (p == 0 or ap_prev_Tv):
                if p == 0:
                    dims_Tv[p] = n-1
                else:
                    om = compute_omega_basis(B, n-1, p, ap_Tv, ap_prev_Tv)
                    dims_Tv[p] = om.shape[1] if om.ndim == 2 else 0
            else:
                dims_Tv[p] = n-1 if p == 0 else 0

        # Relative dimensions
        for p in range(4):
            rel_dim = dims_T.get(p, 0) - dims_Tv.get(p, 0)
            rel_dims[p][rel_dim] += 1

        # Detailed tracking for p=2,3
        r2 = dims_T.get(2, 0) - dims_Tv.get(2, 0)
        r3 = dims_T.get(3, 0) - dims_Tv.get(3, 0)
        rel_details.append((r2, r3))

for p in range(4):
    print(f"\n  R_{p} = Ω_{p}(T) - Ω_{p}(T\\v) dimension distribution:")
    for dim, count in sorted(rel_dims[p].items()):
        print(f"    dim={dim}: {count}")

r2r3_dist = Counter(rel_details)
print(f"\n  Joint (R₂, R₃) distribution:")
for (r2, r3), count in sorted(r2r3_dist.items()):
    print(f"    R₂={r2}, R₃={r3}: {count} ({100*count/len(rel_details):.1f}%)")
    # For H₂^rel = 0: need ker(∂₂^rel) ⊆ im(∂₃^rel)
    # This means rk(∂₃^rel) ≥ dim(ker ∂₂^rel)
    # dim(ker ∂₂^rel) ≤ R₂
    # rk(∂₃^rel) ≤ min(R₂, R₃)
    # So need R₃ ≥ dim(ker ∂₂^rel) ≥ R₂ - rk(∂₂^rel)

# ===== PATHS THROUGH v: COMBINATORIAL STRUCTURE =====
print(f"\n{'='*70}")
print("PATHS THROUGH v: WHAT MAKES THE RELATIVE COMPLEX")
print("=" * 70)

n = 5
# For a given tournament and vertex v:
# 2-paths through v: (v,b,c), (a,v,c), (a,b,v)
# 3-paths through v: (v,b,c,d), (a,v,c,d), (a,b,v,d), (a,b,c,v)

# Count how many of these are TT/DT
tt_through_v = Counter()
dt_through_v = Counter()

for A in all_tournaments(n):
    for v in range(n):
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        # 2-paths through v
        tt_v2 = sum(1 for p in a2 if v in p and A[p[0]][p[2]] == 1)
        nt_v2 = sum(1 for p in a2 if v in p and A[p[0]][p[2]] == 0)
        total_v2 = tt_v2 + nt_v2

        # 3-paths through v
        dt_v3 = sum(1 for p in a3 if v in p and A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1)
        nondt_v3 = sum(1 for p in a3 if v in p and not (A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1))
        total_v3 = dt_v3 + nondt_v3

        tt_through_v[(tt_v2, nt_v2)] += 1
        dt_through_v[(dt_v3, nondt_v3)] += 1

print(f"\n  2-paths through v: (TT, NT) distribution:")
for key, count in sorted(tt_through_v.items()):
    pct = 100*count/(1024*5)
    if pct > 0.5:
        print(f"    TT={key[0]:>2}, NT={key[1]:>2}: {count} ({pct:.1f}%)")

print(f"\n  3-paths through v: (DT, nonDT) distribution:")
for key, count in sorted(dt_through_v.items()):
    pct = 100*count/(1024*5)
    if pct > 0.5:
        print(f"    DT={key[0]:>2}, nonDT={key[1]:>2}: {count} ({pct:.1f}%)")

# ===== POSITION OF v IN PATHS =====
print(f"\n{'='*70}")
print("POSITION OF v IN RELATIVE 2-CYCLES AND 3-CHAINS")
print("=" * 70)

# For the relative complex, the key is:
# A relative 2-cycle z ∈ ker(∂₂^rel) is a 2-chain using v whose boundary
# (projected to relative complex) is zero.
#
# A relative 3-chain w using v has ∂₃^rel(w) which is a relative 2-cycle.
# Need im(∂₃^rel) = ker(∂₂^rel).
#
# v appears at positions 0,1,2 in a 2-path, and 0,1,2,3 in a 3-path.
# The position of v determines which faces use v and which don't.

# Position analysis: for 2-paths (a,b,c) with v at position i:
# pos 0: v=a, faces are (b,c), (v,c), (v,b). Face (b,c) doesn't use v.
# pos 1: v=b, faces are (v,c), (a,c), (a,v). Face (a,c) doesn't use v.
# pos 2: v=c, faces are (b,v), (a,v), (a,b). Face (a,b) doesn't use v.

# In each case, exactly 1 face doesn't use v (goes to Ω₁(T\\v)) and 2 faces use v.

# For 3-paths (a,b,c,d) with v at position i:
# pos 0: v=a, faces (b,c,d), (v,c,d), (v,b,d), (v,b,c). Face (b,c,d) doesn't use v.
# pos 1: v=b, faces (v,c,d), (a,c,d), (a,v,d), (a,v,c). Face (a,c,d) doesn't use v.
# pos 2: v=c, faces (b,v,d), (a,v,d), (a,b,d), (a,b,v). Face (a,b,d) doesn't use v.
# pos 3: v=d, faces (b,c,v), (a,c,v), (a,b,v), (a,b,c). Face (a,b,c) doesn't use v.

# In each case, exactly 1 face doesn't use v and 3 faces use v.

# This is the key structural fact about relative complexes of vertex deletions.

print(f"\n  For vertex deletion:")
print(f"  - Each 2-path through v has exactly 1 face NOT through v")
print(f"  - Each 3-path through v has exactly 1 face NOT through v")
print(f"  - The 'non-v face' of position-i path is the face obtained by deleting v")

# ===== CAN WE PROVE H₂(T,T\\v)=0 ALGEBRAICALLY? =====
print(f"\n{'='*70}")
print("KEY: H₂(T,T\\v) = 0 ALGEBRAIC APPROACH")
print("=" * 70)

# The relative complex R_* at level 2:
# R₂ = Ω₂(T) / Ω₂(T\\v) ≅ span of Ω₂-basis elements using v
# R₃ = Ω₃(T) / Ω₃(T\\v) ≅ span of Ω₃-basis elements using v
#
# ∂₂^rel: R₂ → R₁  and  ∂₃^rel: R₃ → R₂
# H₂^rel = ker(∂₂^rel) / im(∂₃^rel)
#
# For the INDUCTIVE proof to work, we need H₂^rel = 0 for all v.
# This is verified computationally through n=9.
#
# Can we prove it? Let's think about what R₂ and R₃ look like.
#
# R₂: Elements of Ω₂ using v. These are:
#   TT triples through v: (v,b,c) with v→c, (a,v,c) with a→c, (a,b,v) with a→v
#   NT cancellation combos through v
#
# R₃: Elements of Ω₃ using v. These are:
#   DT 4-paths through v
#   Ω₃ cancellation elements through v
#
# The map ∂₃^rel sends a 3-path through v to its boundary, modulo paths NOT through v.
#
# For a DT 4-path (a,b,c,d) with v at position i:
#   ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
#   Exactly 3 faces use v, 1 doesn't. The non-v face is in Ω₂(T\\v).
#   So ∂₃^rel picks up the 3 faces using v.

# ===== COMPUTE: For each (n,v), what is dim(ker ∂₂^rel)? =====
print(f"\n  Computing ker(∂₂^rel) dimensions...")

n = 5
ker_rel2 = Counter()
surjectivity_check = 0
total_verts = 0

for A in all_tournaments(n):
    for v in range(n):
        total_verts += 1

        # Compute using explicit relative complex
        # R₂ elements: Ω₂-basis vectors with nonzero v-path components
        a0 = enumerate_allowed_paths(A, n, 0)
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        a3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, a1, a0)
        om2 = compute_omega_basis(A, n, 2, a2, a1)
        om3 = compute_omega_basis(A, n, 3, a3, a2) if a3 else np.zeros((0,0))

        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
        dim_om3 = om3.shape[1] if om3.ndim == 2 else 0

        if dim_om2 == 0:
            continue

        # Boundary ∂₂: Ω₂ → Ω₁
        bd2 = build_full_boundary_matrix(a2, a1)
        bd2_om = bd2 @ om2  # in A₁ coordinates
        # Project to Ω₁ coords
        coords2, _, _, _ = np.linalg.lstsq(om1, bd2_om, rcond=None)

        # Identify which A₁ paths use v
        v_mask_1 = np.array([v in p for p in a1])
        # Identify which A₂ paths use v
        v_mask_2 = np.array([v in p for p in a2])

        # Relative boundary: ∂₂^rel projects ∂₂ onto v-using components
        # Actually, the relative complex works in Ω coordinates.
        # R₂ = span of Ω₂ basis vectors that project onto v-using A₂ paths.
        # This is complex. Let me use a simpler approach.

        # Delete vertex v
        B = delete_vertex(A, n, v)
        b1 = enumerate_allowed_paths(B, n-1, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b3 = enumerate_allowed_paths(B, n-1, 3)

        om2_B = compute_omega_basis(B, n-1, 2, b2, b1) if b2 else np.zeros((0,0))
        om3_B = compute_omega_basis(B, n-1, 3, b3, b2) if b3 else np.zeros((0,0))

        dim_om2_B = om2_B.shape[1] if om2_B.ndim == 2 and om2_B.shape[0] > 0 else 0
        dim_om3_B = om3_B.shape[1] if om3_B.ndim == 2 and om3_B.shape[0] > 0 else 0

        r2 = dim_om2 - dim_om2_B
        r3 = dim_om3 - dim_om3_B

        # For H₂^rel=0: need im(∂₃^rel) = ker(∂₂^rel)
        # dim(ker ∂₂^rel) ≤ r2
        # For surjectivity: r3 ≥ dim(ker ∂₂^rel)

        ker_rel2[r2 - r3 if r3 <= r2 else 0] += 0  # placeholder

        # Simple bound: r3 ≥ r2 would be sufficient (if ∂₃^rel is injective)
        if r3 >= r2:
            surjectivity_check += 1

print(f"  R₃ ≥ R₂ (sufficient for H₂^rel=0 if ∂₃^rel injective): "
      f"{surjectivity_check}/{total_verts} ({100*surjectivity_check/total_verts:.1f}%)")

# ===== DETAILED: WHEN IS R₃ < R₂? =====
print(f"\n  Cases where R₃ < R₂:")
cases_r3_lt_r2 = 0
for A in all_tournaments(n):
    for v in range(n):
        B = delete_vertex(A, n, v)

        a2_T = enumerate_allowed_paths(A, n, 2)
        a3_T = enumerate_allowed_paths(A, n, 3)
        a1_T = enumerate_allowed_paths(A, n, 1)
        b2 = enumerate_allowed_paths(B, n-1, 2)
        b3 = enumerate_allowed_paths(B, n-1, 3)
        b1 = enumerate_allowed_paths(B, n-1, 1)

        om2_T = compute_omega_basis(A, n, 2, a2_T, a1_T)
        om3_T = compute_omega_basis(A, n, 3, a3_T, a2_T) if a3_T else np.zeros((0,0))
        om2_B = compute_omega_basis(B, n-1, 2, b2, b1) if b2 else np.zeros((0,0))
        om3_B = compute_omega_basis(B, n-1, 3, b3, b2) if b3 else np.zeros((0,0))

        d2T = om2_T.shape[1] if om2_T.ndim == 2 else 0
        d3T = om3_T.shape[1] if om3_T.ndim == 2 and om3_T.shape[0] > 0 else 0
        d2B = om2_B.shape[1] if om2_B.ndim == 2 and om2_B.shape[0] > 0 else 0
        d3B = om3_B.shape[1] if om3_B.ndim == 2 and om3_B.shape[0] > 0 else 0

        r2 = d2T - d2B
        r3 = d3T - d3B

        if r3 < r2:
            cases_r3_lt_r2 += 1
            if cases_r3_lt_r2 <= 5:
                scores = sorted(sum(A[i]) for i in range(n))
                dv = sum(A[v])
                print(f"    scores={scores}, v={v} (d_v={dv}): R₂={r2}, R₃={r3}")

print(f"\n  Total R₃ < R₂ cases: {cases_r3_lt_r2}/{total_verts}")
print(f"  (When R₃ < R₂, H₂^rel=0 requires ∂₃^rel to have nontrivial kernel)")

# ===== THE PROOF STRUCTURE =====
print(f"\n{'='*70}")
print("PROOF STRUCTURE SUMMARY")
print("=" * 70)
print("""
To prove β₂(T) = 0 for all tournaments T:

1. BASE CASE: n=3. Only 8 tournaments, all have β₂=0. ✓

2. INDUCTIVE STEP: Assume β₂(T')=0 for all tournaments on <n vertices.
   For tournament T on n vertices, pick any vertex v.
   Long exact sequence:
     H₂(T\\v) → H₂(T) → H₂(T,T\\v) → H₁(T\\v) → H₁(T)
       = 0         ?         ?

   By induction: H₂(T\\v) = 0.
   Need to show: H₂(T,T\\v) = 0.
   Then: 0 → H₂(T) → 0, so H₂(T) = 0.

3. RELATIVE HOMOLOGY: H₂(T,T\\v) = 0.
   The relative complex R_*(v) = Ω_*(T) / Ω_*(T\\v).
   H₂^rel = ker(∂₂^rel) / im(∂₃^rel).

   Need: every relative 2-cycle is a relative 3-boundary.
   Equivalently: rk(∂₃^rel) ≥ dim(ker ∂₂^rel).

4. VERIFIED: H₂(T,T\\v) = 0 exhaustive n≤6, sampled n=7-9.

   KEY QUESTION: Why is H₂(T,T\\v) always 0?

   Possible approaches:
   a) Show R₃ ≥ ker(∂₂^rel) dimension-wise
   b) Construct an explicit section Z₂^rel → R₃
   c) Use tournament completeness to show "enough" relative 3-chains exist
""")

print("Done.")
