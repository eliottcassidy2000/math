#!/usr/bin/env python3
"""
beta2_tetra_prism.py — Connect S₃×Z₂ "prism" symmetry to β₂=0

The paper's Grid(n) has S₃ acting on barycentric coords (r',c',k')
with r'+c'+k'=n-3. This is the symmetry of a 2-simplex (triangle).
Combined with Z₂ (complement), it gives S₃×Z₂ = symmetry of prism.

HYPOTHESIS: The Ω₂/Ω₃ complex has a similar structure:
- Each DT 4-path (a,b,c,d) creates a "tetrahedron" of 4 triangle faces
- The boundary ∂₃ = simplicial boundary of tetrahedron
- The S₃ acting on (b,c, a↔d) encodes which face is TT vs NT
- β₂=0 follows from these tetrahedra covering all 2-cycles

Key question: can we map the grid position (r',c',k') to the 
DT path structure and show the acyclicity follows?

Author: opus-2026-03-08-S45
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(pairs)

print("=" * 70)
print("TETRAHEDRAL / PRISM CONNECTION TO β₂=0")
print("=" * 70)

# The grid Grid(n) has positions (r,c) with 1≤r, 1≤c, r+c≤n-1.
# In barycentric: r' = r-1, c' = c-1, k' = n-1-r-c.
# Each grid position corresponds to a tile (non-backbone arc pair).
# 
# For n=5: Grid(5) has positions with r+c ≤ 4, r≥1, c≥1.
# That's a triangle with vertices at (1,1), (1,3), (3,1).
# 6 positions = C(4,2) = 6 tiles.
#
# The S₃ symmetry permutes (r',c',k'):
# σ: (r',c',k') → (c',r',k')  [swap r'↔c']
# τ: (r',c',k') → (r',k',c')  [swap c'↔k']
# φ (complement): (r',c',k') → (k',c',r') [swap r'↔k', but also flip)
#
# For path homology at level 2-3:
# Ω₂: transitive triples (a,b,c) with a→b→c, a→c
# Ω₃: DT 4-paths (a,b,c,d) with a→b→c→d, a→c, b→d
#
# A TT triple (a,b,c) uses 3 vertices from [n].
# A DT 4-path (a,b,c,d) uses 4 vertices from [n].
# The "position" of a triple within the tournament is determined by
# how a,b,c relate to the rest of the tournament.

# KEY STRUCTURAL PARALLEL:
# 
# Grid(n) position (r',c',k'):
#   r' = "row parameter" = position of tile in one dimension
#   c' = "column parameter" = position in another dimension  
#   k' = "complement" = remaining degrees of freedom
#
# DT 4-path (a,b,c,d):
#   The 3 "internal" edges are a→c, b→d, and a↔d.
#   These 3 relationships have a natural S₃ structure.
#   a→c is guaranteed (DT condition #1)
#   b→d is guaranteed (DT condition #2)
#   a↔d is the FREE edge
#
# If we think of (a→c, b→d, a↔d) as three binary variables,
# the S₃ permuting them gives different "types" of DT paths.

# Let me check: how many DT paths exist as a function of tournament type?
print("\nDT path analysis:")

for bits in [0, 2, 8, 15, 31]:
    if bits >= (1 << m): break
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    
    ap3 = enumerate_allowed_paths(A, n, 3)
    scores = sorted(sum(A[i]) for i in range(n))
    
    # Classify DT paths by the a↔d edge
    dt_ad_fwd = 0  # a→d (super-DT)
    dt_ad_bwd = 0  # d→a (DT but not super-DT)
    non_dt = 0
    
    for path in ap3:
        a, b, c, d = path
        if A[a][c] and A[b][d]:
            if A[a][d]:
                dt_ad_fwd += 1
            else:
                dt_ad_bwd += 1
        else:
            non_dt += 1
    
    print(f"  T#{bits} (scores={scores}): DT(a→d)={dt_ad_fwd}, DT(d→a)={dt_ad_bwd}, non-DT={non_dt}")

# The "super-DT" paths (a→d too) are actually the most transitive 4-paths.
# For the transitive tournament: ALL 4-paths are super-DT.

# SIMPLEX STRUCTURE:
# A simplex Δ³ on 4 vertices {a,b,c,d} has:
#   - 4 triangular faces: {b,c,d}, {a,c,d}, {a,b,d}, {a,b,c}
#   - 6 edges: all pairs
#   - 4 vertices
# 
# The boundary ∂₃(Δ³) = Σ (-1)^i [omit vertex i]
# = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# 
# For DT path (a,b,c,d): ∂₃ = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# This is EXACTLY the simplicial boundary!
# 
# Face status:
# - (a,b,c): a→b, b→c, a→c (TT). Always TT for DT.
# - (b,c,d): b→c, c→d, b→d (TT). Always TT for DT.
# - (a,c,d): a→c, c→d, a→d? TT iff a→d. NT iff d→a.
# - (a,b,d): a→b, b→d, a→d? TT iff a→d. NT iff d→a.
# 
# So the 4 faces decompose as:
#   2 always-TT faces: (a,b,c), (b,c,d)
#   2 conditional faces: (a,c,d), (a,b,d) — TT iff a→d
# 
# The Z₂ choice (a→d vs d→a) determines whether the tetrahedron
# has 4 TT faces or 2 TT + 2 NT faces.

print(f"\n{'='*70}")
print("SIMPLICIAL STRUCTURE OF DT BOUNDARIES")
print("=" * 70)

# For the regular tournament (c₃=5):
A_reg = [[0, 1, 1, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 0, 1, 1],
         [1, 0, 0, 0, 1],
         [1, 1, 0, 0, 0]]

ap3 = enumerate_allowed_paths(A_reg, n, 3)
ap2 = enumerate_allowed_paths(A_reg, n, 2)
ap1 = enumerate_allowed_paths(A_reg, n, 1)

om3 = compute_omega_basis(A_reg, n, 3, ap3, ap2)
d3 = om3.shape[1] if om3.ndim == 2 else 0

print(f"\nRegular T (c₃=5): {d3} Ω₃ elements")
for col in range(d3):
    v = om3[:, col]
    terms = [(tuple(ap3[i]), round(v[i],3)) for i in range(len(v)) if abs(v[i]) > 1e-8]
    print(f"  e{col}: {terms}")
    
    # For each DT term, show face decomposition
    for path, coeff in terms:
        a, b, c, d = path
        is_dt = A_reg[a][c] and A_reg[b][d]
        ad = "a→d" if A_reg[a][d] else "d→a"
        faces_tt = sum([A_reg[a][c], A_reg[b][d], 
                        A_reg[a][d] if is_dt else 0,
                        A_reg[a][d] if is_dt else 0])
        print(f"    ({a},{b},{c},{d}): {'DT' if is_dt else 'nDT'}, {ad}, "
              f"faces: (abc)[TT] (bcd)[TT] (acd)[{'TT' if A_reg[a][d] else 'NT'}] "
              f"(abd)[{'TT' if A_reg[a][d] else 'NT'}]")

# The PRISM connection:
# 
# Each DT 4-path contributes a "tetrahedral" boundary.
# The tetrahedra "tile" the 2-cycle space.
# The S₃ acts on which triple {a,b,c,d}\{x} we choose as a face.
# The Z₂ is the a↔d direction choice.
#
# The prism = triangle × interval.
# - Triangle: the 3 "middle" face choices (omit a, b, c, or d)
# - Interval: the a→d / d→a choice
#
# In the TILING MODEL:
# A tile at position (r,c) in Grid(n) represents an arc between
# two non-backbone vertices. The S₃ permutes the barycentric
# coordinates (r',c',k') of this tile.
#
# The DT path structure maps:
# - A DT path (a,b,c,d) involves vertices from the backbone path
# - The internal structure (which faces are TT vs NT) is controlled
#   by the a↔d arc
# - This a↔d arc IS a tile in Grid(n)!

# CRUCIAL OBSERVATION:
# A DT path (a,b,c,d) in a Hamiltonian path P = (v₁,...,vₙ):
# a,b,c,d are 4 consecutive-ish vertices in P.
# The "extra" edge a↔d is a NON-BACKBONE arc = a TILE.
# Its position in Grid(n) determines the barycentric coordinates.
#
# The S₃ action on (r',c',k') corresponds to permutations of
# (a,d) and the intermediate vertices.

print(f"\n{'='*70}")  
print("TILE POSITIONS OF DT PATHS")
print("=" * 70)

# For n=5 transitive tournament, HP = 0→1→2→3→4
# A DT path (a,b,c,d) with a→b→c→d, a→c, b→d:
# In the HP ordering, these are consecutive or near-consecutive.

# Grid(5) tiles correspond to arcs (i,j) with i < j in HP and j-i ≥ 2
# (non-backbone arcs). These are:
# (0,2), (0,3), (0,4), (1,3), (1,4), (2,4) — all 6 tiles.
# Grid position: (r,c) = (i+1, j-i-1) = (i+1, gap-1)

print("\nGrid(5) tiles for transitive T:")
for i in range(n):
    for j in range(i+2, n):
        r = i + 1
        c = j - i - 1
        rp = r - 1
        cp = c - 1
        kp = n - 1 - r - c
        print(f"  arc ({i},{j}): grid ({r},{c}), bary ({rp},{cp},{kp})")

# For each DT path, identify its "free edge" tile position
print("\nDT paths and their free-edge tiles:")
A_trans = [[0,1,1,1,1],[0,0,1,1,1],[0,0,0,1,1],[0,0,0,0,1],[0,0,0,0,0]]
ap3 = enumerate_allowed_paths(A_trans, n, 3)
for path in ap3:
    a, b, c, d = tuple(path)
    # DT conditions: a→c, b→d (both satisfied in transitive)
    # Free edge: a→d (always forward in transitive)
    r = a + 1
    gap = d - a
    c_pos = gap - 1
    rp, cp, kp = r-1, c_pos-1, n-1-r-c_pos
    print(f"  DT ({a},{b},{c},{d}): free edge ({a},{d}), "
          f"grid ({r},{c_pos}), bary ({rp},{cp},{kp})")

# The DT paths in the transitive tournament on 5 vertices are:
# (0,1,2,3), (0,1,2,4), (0,1,3,4), (0,2,3,4), (1,2,3,4)
# Free edges: (0,3), (0,4), (0,4), (0,4), (1,4)
# Wait, there's overlap — multiple DT paths can share the same free edge.

# Actually for transitive: ALL 4-paths are DT (and super-DT).
# There are C(5,4)×... no, let me just count.
print(f"\n|A₃| for transitive n=5: {len(ap3)}")
print(f"Paths: {[tuple(p) for p in ap3]}")

# These are all monotone 4-paths: (a,b,c,d) with a<b<c<d.
# There are C(5,4) = 5 such paths.

# The free edge for (a,b,c,d) is always (a,d).
# Grid positions of free edges:
# (0,1,2,3): free=(0,3), grid (1,2), bary (0,1,1)
# (0,1,2,4): free=(0,4), grid (1,3), bary (0,2,0)
# (0,1,3,4): free=(0,4), grid (1,3), bary (0,2,0) [same tile!]
# (0,2,3,4): free=(0,4), grid (1,3), bary (0,2,0) [same tile!]
# (1,2,3,4): free=(1,4), grid (2,2), bary (1,1,0)

# So 3 DT paths share the tile (0,4), and the others have unique tiles.
# This DEGENERACY is important — it means the boundary map has
# higher rank than the number of distinct tiles.

print(f"\n{'='*70}")
print("PROOF SKETCH VIA TETRAHEDRAL DECOMPOSITION")
print("=" * 70)

print("""
OUTLINE: β₂=0 via simplicial acyclicity

1. SIMPLICIAL STRUCTURE: Each DT 4-path (a,b,c,d) gives a
   "simplicial tetrahedron" with boundary
   ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
   
   This is exactly the boundary of a 3-simplex Δ³ = conv{a,b,c,d}.

2. FACE STRUCTURE: For a DT path:
   - 2 faces always TT: (a,b,c), (b,c,d)
   - 2 faces depend on a↔d: (a,c,d), (a,b,d) — TT iff a→d
   
   The Z₂ choice (a→d vs d→a) is the "complement" direction.

3. S₃ ACTION: The 3 barycentric coordinates (r',c',k') of the
   free edge (a,d) in Grid(n) encode:
   - r' = position of a
   - c' = gap from a to d minus 2
   - k' = remaining space after d
   
   S₃ permutes these, giving equivalent DT paths with permuted
   internal structure.

4. COVERING THEOREM (TO PROVE):
   The DT tetrahedra COVER all 2-cycles in Ω₂.
   Equivalently: the span of {∂₃(DT path)} = Z₂.
   
   This is where tournament completeness is essential:
   - Every pair has an edge → enough DT paths exist
   - The complete graph K_n has H₂=0 (contractible)
   - Tournament is "half" of K_n but preserves acyclicity at level 2

5. CONNECTION TO PRISM:
   The S₃×Z₂ symmetry of Grid(n) corresponds to:
   - S₃: permuting internal structure of DT paths
   - Z₂: choosing a→d vs d→a (complement direction)
   
   This is the symmetry group of the triangular PRISM:
   - Triangle face = 3 barycentric choices
   - Prism height = 2 complement choices (a→d, d→a)
   
   The PRISM encodes all possible DT paths through a given
   set of 4 vertices, and their boundary contributions.
""")

# Verify: do DT paths ALONE span Z₂?
print("Verification: DT coverage of Z₂")
dt_coverage = Counter()

for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap1 = enumerate_allowed_paths(A, n, 1)
    if not ap2 or not ap3: continue
    
    om2 = compute_omega_basis(A, n, 2, ap2, ap1)
    om1 = compute_omega_basis(A, n, 1, ap1, enumerate_allowed_paths(A, n, 0))
    d2 = om2.shape[1] if om2.ndim == 2 else 0
    if d2 == 0: continue
    
    # Z₂
    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2 = d2 - rk2
    
    # DT boundaries in A₂ coords
    ap3_tuples = [tuple(p) for p in ap3]
    dt_indices = [i for i, (a,b,c,d) in enumerate(ap3_tuples) if A[a][c] and A[b][d]]
    
    if dt_indices:
        bd3 = build_full_boundary_matrix(ap3, ap2)
        dt_boundaries = bd3[:, dt_indices]
        
        # Project into Ω₂
        dt_in_om2 = np.linalg.lstsq(om2, dt_boundaries, rcond=None)[0]
        rk_dt = np.linalg.matrix_rank(dt_in_om2, tol=1e-8)
    else:
        rk_dt = 0
    
    sufficient = rk_dt >= z2
    dt_coverage[sufficient] += 1
    
    if not sufficient:
        scores = sorted(sum(A[i]) for i in range(n))
        print(f"  DT insufficient! T#{bits} scores={scores}, Z₂={z2}, DT_rank={rk_dt}")

print(f"\nDT coverage sufficient: {dt_coverage[True]}/{sum(dt_coverage.values())}")
print(f"DT coverage insufficient: {dt_coverage.get(False, 0)}/{sum(dt_coverage.values())}")

print("\nDone.")
