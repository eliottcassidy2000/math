#!/usr/bin/env python3
"""
self_dual_dimensions_s134.py — opus-2026-03-21-S134

SELF-DUAL DIMENSIONS IN THE SIMPLEX AND TOURNAMENT THEORY

At n=5: C(5,2) = C(5,3) = 10.
Edges = triangles. The 1-skeleton and 2-skeleton have the SAME SIZE.

This happens when k = (n-1)/2, i.e., the MIDDLE dimension of the simplex.
For odd n: the middle dimension has C(n, (n+1)/2) faces.
For even n: no exact self-duality (C(n, n/2) != C(n, n/2+1) in general).

THE DEEP QUESTION: What does self-duality mean for tournaments?

In tournament terms:
  k=1 faces (edges) = arcs = the tournament data
  k=2 faces (triangles) = 3-cycles = the atomic structure
  Self-duality at n=5: every arc has exactly one "dual" triangle partner.

Actually: C(5,2) = C(5,3) = 10 means there are as many arcs as triangles.
But the CORRESPONDENCE is not 1-to-1 in the obvious way.
Each arc is in exactly 3 triangles (C(3,1) = 3 = n-2 for each edge).
Each triangle has exactly 3 arcs.
The incidence is 3:3, giving a 3-regular bipartite graph on 10+10 vertices.

This bipartite graph IS the Petersen graph on one side!

Author: opus-2026-03-21-S134
"""

import sys
import numpy as np
from math import comb, factorial, sqrt
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

tau = 1.8392867552
phi = (1 + sqrt(5)) / 2

print("=" * 72)
print("  SELF-DUAL DIMENSIONS")
print("  opus-2026-03-21-S134")
print("=" * 72)

# ========================================================================
# PART 1: When does self-duality occur?
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: SELF-DUALITY IN THE SIMPLEX")
print("=" * 72)

print(f"\n  C(n, k) = C(n, n-k) always. But C(n, k) = C(n, k+1) is special.")
print(f"\n  Face counts C(n, k+1) for the (n-1)-simplex:")
print(f"  {'n':>3s}", end="")
for k in range(6):
    print(f"  {'k='+str(k):>6s}", end="")
print()

for n in range(3, 10):
    print(f"  {n:3d}", end="")
    for k in range(min(6, n)):
        c = comb(n, k+1)
        # Mark self-dual pairs
        marker = ""
        if k < n-1 and comb(n, k+1) == comb(n, k+2):
            marker = "*"
        print(f"  {c:>5d}{marker}", end="")
    print()

print(f"""
  * = self-dual pair: C(n, k+1) = C(n, k+2)

  This happens when k+1 = n - (k+2), i.e., 2k + 3 = n.
  So: n = 2k + 3 (odd), and the self-dual dimension is k.

  n=3 (k=0): C(3,1) = C(3,2) = 3.  Vertices = edges.
  n=5 (k=1): C(5,2) = C(5,3) = 10. Edges = triangles.
  n=7 (k=2): C(7,3) = C(7,4) = 35. Triangles = tetrahedra.
  n=9 (k=3): C(9,4) = C(9,5) = 126. Tetrahedra = 4-simplices.

  FOR TOURNAMENT THEORY:
  n=3: vertices = edges. Self-duality at the simplest level.
  n=5: arcs = triangles. THE KEY SELF-DUALITY.
  n=7: triangles = tetrahedra. Higher-order self-duality.
  n=9: tetrahedra = 4-simplices.

  THE ODD PATTERN: self-duality occurs at ALL odd n.
  And the self-dual dimension INCREASES by 1 for each step of 2 in n.
""")

# ========================================================================
# PART 2: The edge-triangle incidence at n=5
# ========================================================================
print("=" * 72)
print("  PART 2: THE EDGE-TRIANGLE INCIDENCE AT n=5")
print("=" * 72)

n = 5
edges = list(combinations(range(n), 2))
tris = list(combinations(range(n), 3))

print(f"  n={n}: {len(edges)} edges, {len(tris)} triangles")

# Incidence matrix: I[e, t] = 1 if edge e is in triangle t
inc = np.zeros((len(edges), len(tris)), dtype=int)
for i, e in enumerate(edges):
    for j, t in enumerate(tris):
        if set(e).issubset(set(t)):
            inc[i][j] = 1

# Properties of the incidence matrix
print(f"\n  Incidence matrix I (edges x triangles): {inc.shape}")
print(f"  Row sums (edges per triangle): {inc.sum(axis=1).tolist()}")
print(f"  Col sums (triangles per edge... wait, reversed)")
print(f"  Each edge is in {inc[0].sum()} triangles")
print(f"  Each triangle contains {inc.T[0].sum()} edges")

# The incidence matrix is 10x10. What are its properties?
print(f"\n  Incidence matrix I is {inc.shape[0]}x{inc.shape[1]} = SQUARE!")
print(f"  det(I) = {int(round(np.linalg.det(inc.astype(float))))}")
print(f"  rank(I) = {np.linalg.matrix_rank(inc)}")

# Eigenvalues
eigs = sorted(np.linalg.eigvals(inc.astype(float)).real, reverse=True)
print(f"  Eigenvalues: {[round(e, 4) for e in eigs]}")

# Is I related to the Petersen graph?
# The Petersen graph K(5,2) adjacency: edge (i,j) ~ (k,l) iff disjoint.
pet_adj = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(10):
        if i != j and not set(edges[i]) & set(edges[j]):
            pet_adj[i][j] = 1

# I * I^T = ?
IIT = inc @ inc.T
print(f"\n  I * I^T (edge-edge via shared triangles):")
print(f"  Diagonal: {IIT[0,0]} (each edge is in {IIT[0,0]} triangles)")
print(f"  Off-diagonal patterns:")
off_diag = set()
for i in range(10):
    for j in range(i+1, 10):
        off_diag.add(IIT[i,j])
print(f"    Distinct off-diagonal values: {sorted(off_diag)}")

# Check: is I*I^T related to the Petersen adjacency or complement?
print(f"\n  I*I^T vs Petersen adjacency:")
print(f"  I*I^T = 3*Identity + complement(Petersen)?")
comp_pet = 1 - pet_adj - np.eye(10, dtype=int)  # complement of Petersen (excl diagonal)
check = 3 * np.eye(10, dtype=int) + comp_pet
print(f"  Match: {np.array_equal(IIT, check)}")

# Actually: I*I^T[i,j] = number of triangles containing BOTH edges i and j.
# Two edges share a triangle iff they share a vertex (non-disjoint).
# So I*I^T[i,j] = 1 if edges share a vertex (and the shared vertex + 2 endpoints make a triangle)
# and 0 if edges are disjoint.
# PLUS the diagonal is 3 (each edge in 3 triangles).

# So: I*I^T = 3*I_10 + J(5,2) where J(5,2) is the Johnson graph adjacency!
# J(5,2) has edges between pairs that share exactly 1 element.

john_adj = np.zeros((10, 10), dtype=int)
for i in range(10):
    for j in range(i+1, 10):
        if len(set(edges[i]) & set(edges[j])) == 1:
            john_adj[i][j] = john_adj[j][i] = 1

print(f"\n  I*I^T = 3*I + J(5,2) (Johnson graph)?")
check2 = 3 * np.eye(10, dtype=int) + john_adj
print(f"  Match: {np.array_equal(IIT, check2)}")

print(f"""
  BEAUTIFUL: I * I^T = 3*I + J(5,2)

  The edge-triangle incidence matrix, when multiplied by its transpose,
  gives 3 times the identity PLUS the Johnson graph adjacency.

  The Johnson graph J(5,2):
    Vertices: 2-subsets of [5] (same as edges of K_5)
    Edges: pairs sharing exactly 1 element
    J(5,2) is the COMPLEMENT of the Petersen graph K(5,2)!

  So: I * I^T = 3*I + complement(K(5,2))
             = 3*I + (J - I - K(5,2))
             = (3-1)*I + J - K(5,2)... not quite clean.

  More precisely: I*I^T = 3*I + J(5,2) where J(5,2) + K(5,2) + I = J (all-ones).
  So: I*I^T = 3*I + (J - I - K(5,2)) = 2*I + J - K(5,2).

  THE INCIDENCE IS CONTROLLED BY THE PETERSEN GRAPH.
  Self-duality at n=5 connects edges to triangles via the Johnson/Petersen structure.
""")

# ========================================================================
# PART 3: The incidence as a map between chain spaces
# ========================================================================
print("=" * 72)
print("  PART 3: INCIDENCE AS A CHAIN MAP")
print("=" * 72)
print(f"""
  The incidence matrix I is a map: C_1 → C_2
  (from the 1-chain space to the 2-chain space).

  In the simplicial chain complex:
    C_0 ←d₁← C_1 ←d₂← C_2 ←d₃← C_3

  The BOUNDARY MAP d₂: C_2 → C_1 sends each triangle to its boundary edges.
  In matrix form: d₂ is a 10×10 matrix (edges × triangles).

  Is d₂ = I^T? Let me check.

  d₂({{a,b,c}}) = {{b,c}} - {{a,c}} + {{a,b}} (oriented boundary).
  But our I is UNSIGNED: I[e,t] = 1 if edge e is a face of triangle t.

  So: d₂ is a SIGNED version of I^T.
  |d₂| = I^T (absolute values).
  The SIGN pattern encodes the orientation.

  For the TOURNAMENT: the orientation ε_ij already provides signs.
  The tournament + the boundary map d₂ together determine c₃ (3-cycle count):
    c₃ = number of triangles where the tournament orientation makes a cycle.

  A triangle {{a,b,c}} with tournament ε gives a 3-cycle iff
    ε_ab * ε_bc * ε_ca = +1 (all three arcs go "forward" in the cycle).

  Actually: for a DIRECTED 3-cycle a→b→c→a, we need A[a][b]=A[b][c]=A[c][a]=1.
  In terms of ε: ε_ab * ε_bc * ε_ca = ??? depends on convention.

  THE 3-CYCLE INDICATOR: for each triangle, the tournament is either
  a 3-cycle (one of 2 directions) or transitive (one of 6 orderings,
  but actually 3 transitive out of 8 orientations on 3 arcs...
  wait: for 3 edges, 2^3 = 8 orientations. Of these, 2 are 3-cycles
  and 6 are transitive triples.)

  THE SELF-DUALITY: at n=5, the tournament (living on edges) and the
  3-cycle structure (living on triangles) occupy EQUAL-SIZED spaces.
  The map between them (the incidence I) is SQUARE and INVERTIBLE
  (det ≠ 0 if the incidence has full rank).
""")

print(f"  det(I) = {int(round(np.linalg.det(inc.astype(float))))}")
print(f"  rank(I) = {np.linalg.matrix_rank(inc)}")
if np.linalg.matrix_rank(inc) == 10:
    print(f"  I is INVERTIBLE! The map edges → triangles is bijective.")
    I_inv = np.linalg.inv(inc.astype(float))
    print(f"  I^{{-1}} exists. Max entry: {np.max(np.abs(I_inv)):.4f}")
else:
    print(f"  I is NOT invertible. rank < 10.")

# ========================================================================
# PART 4: What self-duality means for the OCR
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: SELF-DUALITY AND THE OCR")
print("=" * 72)
print(f"""
  THE OCR AT SELF-DUAL n:
    n=3: OCR = 1.000 (self-dual: vertices = edges, C(3,1)=C(3,2)=3)
    n=5: OCR = 0.970 (self-dual: edges = triangles, C(5,2)=C(5,3)=10)
    n=7: OCR = 0.957 (self-dual: triangles = tetrahedra, C(7,3)=C(7,4)=35)
    n=9: OCR = 0.966 (self-dual: C(9,4)=C(9,5)=126)

  These are ALL ODD n (where self-duality occurs).
  But the OCR is NOT higher at self-dual n vs non-self-dual.
  n=6 (not self-dual): OCR = 0.959. Between n=5 and n=7.

  SO: self-duality alone doesn't determine the OCR.

  BUT: self-duality determines the STRUCTURE of the OCR:
  At self-dual n=5: the tournament (edges) and the cycle structure
  (triangles) have the same size. The incidence matrix I is SQUARE.
  This means: in principle, you can INVERT I to go from cycles to edges.

  THE INVERTIBILITY:
  rank(I) = {np.linalg.matrix_rank(inc)} at n=5. {'INVERTIBLE' if np.linalg.matrix_rank(inc)==10 else 'NOT invertible'}.
""")

if np.linalg.matrix_rank(inc) < 10:
    print(f"  The incidence is NOT invertible. Self-duality is broken by rank deficiency.")
    print(f"  This means: some edge-information is lost when projecting to triangles.")
    print(f"  Specifically: ker(I) has dimension {10 - np.linalg.matrix_rank(inc)}.")
    print(f"  These kernel directions are tournament configurations that produce")
    print(f"  NO 3-cycles — i.e., sub-transitive structures.")

# ========================================================================
# PART 5: The self-dual dimension sequence
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: THE SELF-DUAL DIMENSION SEQUENCE")
print("=" * 72)

print(f"\n  Self-dual face count: C(n, (n+1)/2) for odd n")
print(f"  {'n':>3s} {'self-dual dim':>13s} {'face count':>10s} {'OCR':>7s} {'is boundary?':>12s}")

ocr_data = {3: 1.0, 5: 0.970, 7: 0.957, 9: 0.966}
for n in range(3, 14, 2):
    k = (n - 1) // 2  # self-dual dimension
    face_count = comb(n, k + 1)
    ocr = ocr_data.get(n, '?')
    boundary = "YES" if n == 5 else ("n=3 trivial" if n == 3 else "")
    print(f"  {n:3d} {k:13d} {face_count:10d} {str(ocr):>7s} {boundary:>12s}")

print(f"""
  The self-dual face count: 3, 10, 35, 126, 462, 1716, ...
  These are the CENTRAL BINOMIAL-like coefficients C(n, (n+1)/2).

  For n=2m+1: C(2m+1, m+1) = C(2m+1, m).
  This is the number of ways to choose the "majority" of elements.

  THE TOURNAMENT INTERPRETATION:
  At self-dual n = 2m+1:
  - k-faces (edges at k=1, triangles at k=2, ...): TOURNAMENT DATA
  - (n-1-k)-faces: DUAL STRUCTURE

  The self-dual face count is the size of BOTH the data and the dual.
  The INCIDENCE between them is a SQUARE matrix.

  FIBONACCI CONNECTION:
  C(3,2) = 3 = F_4
  C(5,3) = 10 = F_5 × 2? No. But 10 = C(5,2) = T(5) - 2.
  C(7,4) = 35 = F_9? No, F_9 = 34. Close!
  C(9,5) = 126 = F_11? No, F_11 = 89. Not Fibonacci.

  Actually these are just the central binomial coefficients.
  C(2m+1, m+1) grows as 4^m / sqrt(pi*m).
""")

# ========================================================================
# PART 6: Self-duality and Poincaré duality
# ========================================================================
print("=" * 72)
print("  PART 6: POINCARE DUALITY AND TOURNAMENTS")
print("=" * 72)
print(f"""
  POINCARE DUALITY for the simplex:
  The simplex Delta_{{n-1}} is a manifold with boundary.
  Its homology: H_k(Delta) = Z for k=0, H_k = 0 for k > 0.
  (The simplex is contractible.)

  But the BOUNDARY of the simplex (the (n-2)-sphere) has:
  H_k(partial Delta) = Z for k=0 and k=n-2, 0 otherwise.
  Poincare duality: H_k ≅ H_{{n-2-k}} for the boundary sphere.

  FOR TOURNAMENTS: we're not looking at the simplex's topology
  but at the COMBINATORICS of face incidence.

  The self-duality C(n,k) = C(n,n-k) is the COMBINATORIAL analog
  of Poincare duality. At the middle dimension:
  C(n, floor(n/2)) = C(n, ceil(n/2)).

  THE PATH HOMOLOGY of the tournament is DIFFERENT:
  beta_0 = 1, beta_1 in {{0,1}}, beta_2 = 0 always, beta_3 varies.
  This is NOT the homology of the simplex.
  It's the homology of the DIRECTED graph (GLMY theory).

  But the SIMPLICIAL structure of the underlying complete graph K_n
  does connect to the tournament via the chain complex.
  The BOUNDARY MAP d_1^T: edges → vertices IS the score projection.

  THE SELF-DUAL DIMENSION IS WHERE EDGES AND TRIANGLES BALANCE.
  This is the dimension where the tournament data (edges) and the
  cycle structure (triangles) have equal "weight" in the simplex.

  AT n=5: this balance is EXACT (10 = 10).
  The tournament and its cycle structure occupy equal space.
  The OCR (97%) measures how tightly they're linked.

  AT n=7: the balance shifts to triangles vs tetrahedra (35 = 35).
  The tournament no longer has self-dual edges vs triangles
  (C(7,2)=21 != C(7,3)=35). But a HIGHER self-duality holds.

  THE PROGRESSION OF SELF-DUAL PHENOMENA:
  n=3: vertex-edge duality → trivial OCR = 1
  n=5: edge-triangle duality → OCR = 97% (the BOUNDARY)
  n=7: triangle-tetrahedron duality → OCR = 96%
  n=9: tetrahedron-4-simplex duality → OCR = 97%

  Each self-dual level involves a HIGHER-ORDER duality.
  The OCR dips at n=7 because that's where the SECOND level of
  self-duality kicks in, requiring HIGHER-ORDER data (c_5, alpha_2)
  that scores don't capture.
""")

# ========================================================================
# PART 7: The middle dimension and the OCR minimum
# ========================================================================
print("=" * 72)
print("  PART 7: WHY n=7 IS THE OCR MINIMUM — SELF-DUAL PERSPECTIVE")
print("=" * 72)
print(f"""
  THE OCR MINIMUM AT n=7 — A SELF-DUAL EXPLANATION:

  At n=5 (self-dual k=1):
    Edges = triangles = 10.
    The tournament DATA (edges) perfectly mirrors the cycle STRUCTURE (triangles).
    Scores (= edge boundary) capture 97% of triangle structure.
    The self-duality means: there's a 1-to-1 correspondence between
    arc configurations and cycle configurations (up to the rank of I).

  At n=7 (self-dual k=2):
    Triangles = tetrahedra = 35.
    The CYCLE structure (triangles, 35) now exceeds the tournament data (edges, 21).
    There are MORE triangles than edges: 35 > 21.
    The tournament can't individually "see" each triangle.
    Scores capture only 21/35 ≈ 60% of the triangle dimensions.

    But OCR is still 96% because c_3 (linearly related to S_2)
    captures the TOTAL triangle count even if not individual triangles.

    The 4% residual comes from the INTERNAL structure of the 35 triangles
    that 21 edges can't individually resolve.

  At n=9 (self-dual k=3):
    Tetrahedra = 4-simplices = 126.
    The cycle structure is even HIGHER-dimensional.
    But the OCR RECOVERS because Var(H) grows factorially,
    drowning the higher-dimensional residual.

  THE SELF-DUAL EXPLANATION OF THE n=7 MINIMUM:
  n=7 is where the TRIANGLE-TO-EDGE RATIO is at its critical value:
    C(7,3)/C(7,2) = 35/21 = 5/3 ≈ 1.67.

  This ratio controls how many cycle degrees of freedom per arc:
""")

print(f"  {'n':>3s} {'edges':>6s} {'tris':>5s} {'ratio':>8s} {'OCR':>7s}")
for n in range(3, 10):
    ed = comb(n, 2)
    tr = comb(n, 3)
    ratio = tr / ed
    ocr = ocr_data.get(n, '?')
    marker = " ← MIN" if n == 7 else ""
    print(f"  {n:3d} {ed:6d} {tr:5d} {ratio:8.3f} {str(ocr):>7s}{marker}")

print(f"""
  The ratio C(n,3)/C(n,2) = (n-2)/3 increases linearly.
  At n=5: ratio = 1.000 (SELF-DUAL!).
  At n=7: ratio = 1.667.
  At n=9: ratio = 2.333.

  The OCR minimum occurs NOT at the self-dual point (n=5, ratio=1)
  but at n=7 (ratio=5/3). This is because:

  At n=5: self-duality means edges and triangles balance perfectly.
  The incidence matrix is square → the boundary captures almost everything.

  At n=7: triangles outnumber edges by 5/3. The tournament data can't
  individually address each triangle. The residual is MAXIMAL.

  At n=9+: the factorial growth of Var(H) overwhelms the triangle excess.
  Even though the ratio keeps growing, the OCR recovers.

  THE OCR MINIMUM = THE POINT WHERE TRIANGLE EXCESS MAXIMALLY IMPACTS
  THE VARIANCE RATIO, BEFORE FACTORIAL GROWTH TAKES OVER.
""")

# ========================================================================
# PART 8: The self-dual Euler characteristic
# ========================================================================
print("=" * 72)
print("  PART 8: EULER CHARACTERISTIC AND SELF-DUALITY")
print("=" * 72)

for n in range(3, 10):
    # Euler characteristic of the boundary of Delta_{n-1} = (n-1)-sphere
    # chi(S^{n-2}) = 1 + (-1)^{n-2}
    chi_sphere = 1 + (-1)**(n-2)
    # chi of the simplex itself
    chi_simplex = 1  # contractible

    # The ALTERNATING SUM of face counts:
    alt_sum = sum((-1)**k * comb(n, k+1) for k in range(n))

    # For the self-dual face count at odd n:
    if n % 2 == 1:
        sd_dim = (n-1)//2
        sd_count = comb(n, sd_dim + 1)
    else:
        sd_count = 0

    print(f"  n={n}: chi(boundary)={chi_sphere}, self-dual face={sd_count}, "
          f"alt_sum={alt_sum}")

print(f"""
  SYNTHESIS: SELF-DUAL DIMENSIONS AND TOURNAMENT THEORY

  1. Self-duality occurs at ALL odd n, at the middle dimension k=(n-1)/2.
  2. The self-dual face count C(n, (n+1)/2) grows as 4^{{n/2}}/sqrt(n).
  3. At n=5: the self-dual incidence I is 10x10 with I*I^T = 3I + J(5,2).
  4. The incidence is controlled by the Johnson/Petersen duality.
  5. The OCR minimum at n=7 occurs where triangle/edge ratio = 5/3,
     maximizing the impact of cycle excess before factorial growth dominates.
  6. Self-duality is the SYMMETRY of the simplex that makes tournaments possible:
     it ensures that the tournament data (edges) has a natural DUAL (triangles)
     that encodes the cycle structure.

  WITHOUT self-duality: the tournament and its cycle structure would
  have incommensurable sizes, making the OCR much harder to analyze.
  WITH self-duality (at odd n): the tournament and cycles live in
  equal-sized spaces, connected by the square incidence matrix I.
  The QUALITY of this connection (how close I is to invertible)
  determines the OCR.
""")

print("\nDone. opus-2026-03-21-S134")
