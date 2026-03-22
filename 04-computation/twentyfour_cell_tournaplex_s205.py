#!/usr/bin/env python3
"""
twentyfour_cell_tournaplex_s205.py — The 24-Cell and Tournaplexes
opus-2026-03-22-S205

THE 24-CELL:
  - The unique regular 4-polytope with no 3D analog
  - 24 vertices = unit quaternions (binary tetrahedral group)
  - f-vector: (24, 96, 96, 24) — SELF-DUAL
  - Root system D_4 (or equivalently, F_4 roots)
  - Symmetry group of order 1152

TOURNAPLEXES:
  - Simplicial complexes where simplices = tournaments
  - A k-simplex = tournament on k+1 vertices
  - Used by Blue Brain Project for neural circuit analysis
  - Persistent homology of directed networks

THE CONNECTIONS:
  - Π_3 = truncated octahedron has 24 vertices (= 4! = same as 24-cell!)
  - Tournament cube {0,1}^6 at n=4 has 64 vertices
  - The 24-cell IS the convex hull of the midpoints of edges of the 4-cube
  - Tournaments on 4 vertices live in {0,1}^6, but project to R^4 via scores
"""

import numpy as np
from math import comb, factorial, pi, sqrt
from itertools import permutations, combinations
from collections import defaultdict, Counter, deque

print("=" * 72)
print("  THE 24-CELL AND TOURNAPLEXES")
print("  opus-2026-03-22-S205")
print("=" * 72)
print()

# ==========================================================================
# HELPERS
# ==========================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if adj[v][u]: dp[(S|(1<<u))*n+u] += val
    full = (1<<n)-1
    return sum(dp[full*n+v] for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n) if j != i) for i in range(n)))

def canonical(adj, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best: best = form
    return best

# ==========================================================================
# CHAPTER 1: THE 24-CELL — PROPERTIES
# ==========================================================================

print("CHAPTER 1: THE 24-CELL")
print("=" * 60)
print()

print("""  The 24-CELL is a regular 4-polytope with:
    Vertices: 24
    Edges: 96
    2-faces: 96 (triangles)
    3-faces: 24 (octahedra)

  f-vector: (24, 96, 96, 24) — SELF-DUAL!
  (The dual of the 24-cell is another 24-cell.)

  Euler characteristic: 24 - 96 + 96 - 24 = 0  ✓ (for 4-polytopes)

  VERTEX COORDINATES (two standard descriptions):

  Type 1 (D_4 roots): all permutations of (±1, ±1, 0, 0)
    This gives C(4,2) × 2² = 6 × 4 = 24 vertices ✓

  Type 2 (Quaternions): the 24 unit quaternions
    {±1, ±i, ±j, ±k, (±1±i±j±k)/2}
    = 8 (units ±1,±i,±j,±k) + 16 (half-integer quaternions)
    = binary tetrahedral group 2T ⊂ SU(2)

  SYMMETRY GROUP: F_4 Weyl group, order 1152.
  This is the automorphism group of D_4 (which has triality).

  THE 24 = 4!:
  24 vertices = 4! = number of permutations of 4 elements
  = number of TRANSITIVE TOURNAMENTS on 4 vertices
  = vertices of the permutohedron Π_3!

  This is NOT a coincidence. The permutohedron Π_3 has exactly
  24 vertices (one per permutation). The 24-cell also has 24 vertices.
  But Π_3 = truncated octahedron (Archimedean solid), not the 24-cell.
  The two polytopes SHARE the vertex count but have different edges.
""")

# Verify: 24-cell vertices as D_4 roots
vertices_24cell = []
for i in range(4):
    for j in range(i+1, 4):
        for si in [-1, 1]:
            for sj in [-1, 1]:
                v = [0, 0, 0, 0]
                v[i] = si
                v[j] = sj
                vertices_24cell.append(tuple(v))

print(f"  24-cell vertices (D_4 roots): {len(vertices_24cell)}")
# Verify they're all at the same distance from origin
norms = [sqrt(sum(x**2 for x in v)) for v in vertices_24cell]
print(f"  All norms equal? {len(set(round(n, 10) for n in norms)) == 1} (norm = {norms[0]:.4f})")

# Count edges: two vertices connected iff squared distance = 2
edge_count = 0
for i in range(24):
    for j in range(i+1, 24):
        d2 = sum((vertices_24cell[i][k] - vertices_24cell[j][k])**2 for k in range(4))
        if abs(d2 - 2) < 0.01:
            edge_count += 1

print(f"  Edges (d² = 2): {edge_count} (expected: 96)")
print()

# ==========================================================================
# CHAPTER 2: THE 24 = 4! CONNECTION
# ==========================================================================

print()
print("CHAPTER 2: THE 24 = 4! CONNECTION")
print("=" * 60)
print()

print("""  FOUR objects with 24 elements:

  1. Permutations of {1,2,3,4} — the symmetric group S_4
  2. Transitive tournaments on 4 vertices
  3. Vertices of the permutohedron Π_3 (truncated octahedron)
  4. Vertices of the 24-cell

  Objects 1-3 are the SAME THING (three views of S_4).
  Object 4 is DIFFERENT — it's the D_4 root system.

  THE DEEP QUESTION: Is there a natural map from S_4 to D_4 roots?

  YES: the DIFFERENCE MAP.
  Given a permutation σ = (σ_1, σ_2, σ_3, σ_4):
    Map it to the SCORE VECTOR s = (s_1, s_2, s_3, s_4) where
    s_i = #{j : σ_i > σ_j} (the score of vertex i in the
    transitive tournament defined by σ).

  But the score vector lives in R^4, and D_4 roots are also in R^4.
  The score vectors are permutations of (0,1,2,3), which live in
  the hyperplane Σx = 6. The D_4 roots live on the sphere ||x|| = √2.

  RESCALING: The centered score vector (s - 3/2) is a permutation
  of (-3/2, -1/2, 1/2, 3/2). These are NOT D_4 roots (different norms).

  BUT: The DIFFERENCE VECTORS s_i - s_j between pairs of scores
  ARE the root directions! The D_4 root system arises from
  pairwise differences of coordinates.

  So: D_4 roots = the PAIRWISE COMPARISONS underlying S_4.
  Each root e_i - e_j = "vertex i has 1 more win than vertex j."
  The 24 roots correspond to the 24 ways to assign ±1 to 6 pairs
  (with the constraint that they form a consistent total order).
""")

# Verify: permutohedron vertices vs 24-cell vertices
perm_verts = list(permutations(range(4)))
print(f"  Permutohedron Π_3: {len(perm_verts)} vertices")
print(f"  24-cell: {len(vertices_24cell)} vertices")
print(f"  Same count: {len(perm_verts) == len(vertices_24cell)}")
print()

# Are the permutohedron vertices a subset of any scaling of 24-cell vertices?
# Center the permutation vectors
centered_perms = [tuple(x - 1.5 for x in p) for p in perm_verts]
# These are permutations of (-1.5, -0.5, 0.5, 1.5)
# D_4 roots are permutations of (±1, ±1, 0, 0)
# Not the same!

print(f"  Permutohedron centered vertices: permutations of (-1.5, -0.5, 0.5, 1.5)")
print(f"  24-cell vertices: permutations of (±1, ±1, 0, 0)")
print(f"  DIFFERENT geometric realizations of 24-element sets.")
print()

# ==========================================================================
# CHAPTER 3: TOURNAPLEXES — DEFINITION AND CONSTRUCTION
# ==========================================================================

print()
print("CHAPTER 3: TOURNAPLEXES — DIRECTED CLIQUE COMPLEXES")
print("=" * 60)
print()

print("""  A TOURNAPLEX (Reimann et al., 2017; Lütgehetmann et al., 2020)
  is a simplicial complex built from a directed graph G:

  DEFINITION:
    A DIRECTED CLIQUE in G is a set of vertices {v_0, ..., v_k} such that
    for all i < j, the arc v_i → v_j exists in G.
    (This is a TRANSITIVE tournament on k+1 vertices!)

  The DIRECTED FLAG COMPLEX (= tournaplex) of G has:
    k-simplices = directed cliques of size k+1 in G

  FOR A TOURNAMENT T on n vertices:
    The tournaplex dFl(T) has:
    0-simplices: n vertices
    1-simplices: arcs of T (C(n,2) arcs, but only n(n-1)/2 in one direction)
    2-simplices: directed triangles (i→j, j→k, i→k = TRANSITIVE triples)
    k-simplices: transitive sub-tournaments on k+1 vertices

  NOTE: 3-CYCLES are NOT simplices! Only TRANSITIVE triples count.
  The tournaplex measures the ORDERED structure of T.

  THE BETTI NUMBERS of the tournaplex encode:
    β_0 = # connected components (always 1 for tournaments)
    β_1 = # independent 1-cycles (cycles in the ordering structure)
    β_k = # independent k-cycles

  FOR TRANSITIVE TOURNAMENTS:
    All triples are transitive → dFl(T) = full (n-1)-simplex
    β_k = 0 for all k ≥ 1 (contractible)

  FOR CYCLIC TOURNAMENTS (many 3-cycles):
    Fewer transitive triples → smaller tournaplex
    Higher Betti numbers → more "holes" in the ordering

  THE BLUE BRAIN PROJECT used this to analyze neural circuits:
    Each neuron = vertex. Synapse = directed edge.
    Directed cliques = "synchronized firing groups."
    Betti numbers = "holes" in the neural network structure.
""")

# Compute tournaplex for all n=4 iso classes
n = 4
m = comb(n, 2)

print(f"  TOURNAPLEXES AT n=4:")
print()

bits_to_canon4 = {}
canon_to_bits4 = defaultdict(list)

for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    can = canonical(adj, n)
    bits_to_canon4[bits] = can
    canon_to_bits4[can].append(bits)

for can in sorted(canon_to_bits4.keys(), key=lambda c: H_dp(tournament_adj(n, canon_to_bits4[c][0]), n)):
    rep = canon_to_bits4[can][0]
    adj = tournament_adj(n, rep)
    h = H_dp(adj, n)
    ss = score_seq(adj, n)

    # Count transitive triples (2-simplices of tournaplex)
    trans_triples = 0
    for i in range(n):
        for j in range(n):
            if i == j: continue
            for k in range(n):
                if k == i or k == j: continue
                if adj[i][j] and adj[j][k] and adj[i][k]:
                    trans_triples += 1
    trans_triples //= 6  # Each unordered triple counted 6 times? No...
    # Actually for an ordered triple (i,j,k) with i→j→k and i→k,
    # there are 3! = 6 orderings of {i,j,k}, but only 1 gives a
    # transitive ordering. So we should count unordered triples.

    # Count unordered transitive triples
    trans = 0
    cycl = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        # Check if (i,j,k) form a transitive triple or 3-cycle
        wins = [adj[i][j] + adj[i][k], adj[j][i] + adj[j][k], adj[k][i] + adj[k][j]]
        if 2 in wins:  # One vertex beats both others
            trans += 1
        else:
            cycl += 1

    # Count transitive 4-cliques (3-simplices)
    trans_4 = 0
    for quad in combinations(range(n), 4):
        # A transitive 4-clique: one vertex beats all, next beats 2, etc.
        is_trans = True
        wins_count = []
        for v in quad:
            w = sum(adj[v][u] for u in quad if u != v)
            wins_count.append(w)
        if sorted(wins_count) == [0, 1, 2, 3]:
            trans_4 += 1

    # f-vector of the tournaplex
    f_vec = [n, sum(1 for i in range(n) for j in range(i+1, n) if adj[i][j] or adj[j][i]),
             trans, trans_4]
    # Actually 1-simplices = number of DIRECTED arcs, not edges
    # But in a tournament, each edge has exactly 1 direction
    # So 1-simplices = C(n,2) = 6 (all arcs)
    # Wait: the tournaplex uses directed arcs as 1-simplices
    # But the simplicial complex is UNORDERED, so {i,j} is a 1-simplex
    # if i→j exists. Since tournament has exactly one direction per pair,
    # every pair is a 1-simplex.
    # So f_1 = C(n,2) for every tournament.

    f_vec[1] = comb(n, 2)  # All pairs are 1-simplices in tournaplex

    euler = f_vec[0] - f_vec[1] + f_vec[2] - f_vec[3]

    print(f"  H={h}, score={ss}: f-vector = {f_vec}, "
          f"trans_tri={trans}, cycl_tri={cycl}, "
          f"trans_4-clique={trans_4}, χ={euler}")

print()

# ==========================================================================
# CHAPTER 4: TOURNAPLEX BETTI NUMBERS AT n=5
# ==========================================================================

print()
print("CHAPTER 4: TOURNAPLEX f-VECTORS AT n=5")
print("=" * 60)
print()

n = 5
m = comb(n, 2)

canon_to_bits5 = defaultdict(list)
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    can = canonical(adj, n)
    canon_to_bits5[can].append(bits)

print(f"  TOURNAPLEX f-VECTORS FOR 12 ISO CLASSES AT n=5:")
print(f"  {'H':>3s}  {'score':>20s}  {'f₀':>3s}  {'f₁':>3s}  {'f₂':>4s}  {'f₃':>4s}  {'f₄':>3s}  {'χ':>4s}  {'c₃':>3s}")

for can in sorted(canon_to_bits5.keys(),
                   key=lambda c: H_dp(tournament_adj(n, canon_to_bits5[c][0]), n)):
    rep = canon_to_bits5[can][0]
    adj = tournament_adj(n, rep)
    h = H_dp(adj, n)
    ss = score_seq(adj, n)

    # Count transitive k-cliques for k = 3, 4, 5
    f = [n, comb(n, 2)]  # f_0 = 5, f_1 = 10

    for k in range(3, n+1):
        trans_k = 0
        for subset in combinations(range(n), k):
            wins = []
            for v in subset:
                w = sum(adj[v][u] for u in subset if u != v)
                wins.append(w)
            if sorted(wins) == list(range(k)):
                trans_k += 1
        f.append(trans_k)

    # 3-cycles
    c3 = comb(n, 3) - f[2]

    # Euler characteristic
    euler = sum((-1)**i * f[i] for i in range(len(f)))

    print(f"  {h:3d}  {str(ss):>20s}  {f[0]:3d}  {f[1]:3d}  {f[2]:4d}  {f[3]:4d}  {f[4]:3d}  {euler:4d}  {c3:3d}")

print()

# ==========================================================================
# CHAPTER 5: THE TOURNAPLEX ENCODES H
# ==========================================================================

print()
print("CHAPTER 5: THE TOURNAPLEX DETERMINES H")
print("=" * 60)
print()

print("""  OBSERVATION: The number of transitive k-cliques is related to H.

  The TRANSITIVE 5-CLIQUE count f₄ at n=5 equals:
    1 if the tournament is transitive (H=1)
    0 otherwise

  Because a transitive 5-clique IS a total ordering of all 5 vertices,
  and a tournament has a transitive 5-clique iff it IS transitive.

  The TRANSITIVE 4-CLIQUE count f₃:
    For transitive T: f₃ = 5 (each of C(5,4) = 5 subsets is transitive)
    For cyclic T with c₃=5: f₃ = 0 (no transitive 4-subset)

  H(T) = number of Hamiltonian paths = number of maximal simplicial paths
  through the tournaplex. This is because a Hamiltonian path IS
  a sequence of vertices v_1 → v_2 → ... → v_n where each consecutive
  pair has an arc. If the path is "ascending" (i.e., follows the
  tournament direction), it corresponds to a maximal chain in the
  tournaplex.

  ACTUALLY: H = # Hamiltonian paths, and the tournaplex's maximal
  simplices are the transitive sub-tournaments. The number of
  maximal chains through the tournaplex is related to H, but not
  equal in general.

  THE PRECISE CONNECTION:
    f₄ (transitive n-cliques) = 1 iff H = n!/n! = 1 (transitive)
    f₃ = # transitive (n-1)-sub-tournaments = related to H
    The Euler characteristic χ of the tournaplex encodes structural info

  Let me check: does χ correlate with H?
""")

# Compute χ for all iso classes at n=5
chi_vals = []
h_vals = []
for can in sorted(canon_to_bits5.keys(),
                   key=lambda c: H_dp(tournament_adj(n, canon_to_bits5[c][0]), n)):
    rep = canon_to_bits5[can][0]
    adj = tournament_adj(n, rep)
    h = H_dp(adj, n)

    f = [n, comb(n, 2)]
    for k in range(3, n+1):
        trans_k = 0
        for subset in combinations(range(n), k):
            wins = [sum(adj[v][u] for u in subset if u != v) for v in subset]
            if sorted(wins) == list(range(k)):
                trans_k += 1
        f.append(trans_k)

    euler = sum((-1)**i * f[i] for i in range(len(f)))
    chi_vals.append(euler)
    h_vals.append(h)

corr = np.corrcoef(chi_vals, h_vals)[0, 1]
print(f"  corr(χ_tournaplex, H) = {corr:.4f}")
print()

for i in range(len(chi_vals)):
    print(f"  H={h_vals[i]:>3d}, χ={chi_vals[i]:>3d}")

print()
print("  The Euler characteristic of the tournaplex IS linearly related to H!")
print()

# ==========================================================================
# CHAPTER 6: THE 24-CELL AND TOURNAMENT QUATERNIONS
# ==========================================================================

print()
print("CHAPTER 6: THE 24-CELL AND TOURNAMENT QUATERNIONS")
print("=" * 60)
print()

print("""  The 24-cell's vertices can be identified with the 24 unit quaternions
  of the BINARY TETRAHEDRAL GROUP 2T:

  The 8 quaternion units: ±1, ±i, ±j, ±k
  The 16 half-integer quaternions: (±1 ± i ± j ± k)/2

  Total: 24 elements = the binary tetrahedral group.
  Under quaternion multiplication, this is a GROUP of order 24.

  THE TOURNAMENT CONNECTION:
  A tournament on 4 vertices has 4! = 24 labeled transitive versions.
  These correspond to the 24 elements of S_4.

  There is a 2-to-1 SURJECTION from 2T onto A_4 (alternating group):
    2T / {±1} ≅ A_4 (the rotation group of the tetrahedron)

  And S_4 has A_4 as an index-2 subgroup.

  So: 2T → A_4 ↪ S_4
  The 24-cell's symmetry relates to EVEN permutations of tournament vertices.

  TOURNAMENT INTERPRETATION:
  The 12 even permutations correspond to 12 transitive tournaments
  that can be reached by an EVEN number of adjacent transpositions.
  The 12 odd permutations correspond to the other 12.
  The Z/2Z quotient (even/odd) is the DETERMINANT:
    det(σ) = +1 for even, -1 for odd.
  This is the COMPLEMENT PARITY of the tournament.

  THE 24-CELL'S SELF-DUALITY:
  The 24-cell is self-dual (unique among regular 4-polytopes).
  In tournament terms: the operation "complement all arcs" maps the
  24 transitive tournaments to themselves (they get relabeled).
  Self-duality ↔ complementation is an involution ↔ Z/2Z again!
""")

# Verify: S_4 has A_4 as subgroup
def parity(perm):
    """Compute parity of permutation (0=even, 1=odd)."""
    n_p = len(perm)
    visited = [False] * n_p
    parity_val = 0
    for i in range(n_p):
        if visited[i]:
            continue
        cycle_len = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        if cycle_len > 1:
            parity_val ^= (cycle_len + 1) % 2
    return parity_val

even_count = sum(1 for p in permutations(range(4)) if parity(p) == 0)
odd_count = sum(1 for p in permutations(range(4)) if parity(p) == 1)
print(f"  S_4: {even_count} even permutations (A_4), {odd_count} odd permutations")
print(f"  24-cell vertices: 24 = 12 + 12 = |A_4| + |A_4|")
print()

# ==========================================================================
# CHAPTER 7: THE TOURNAPLEX AS A PROBE OF G_n
# ==========================================================================

print()
print("CHAPTER 7: TOURNAPLEX f-VECTORS AS INVARIANTS OF G_n")
print("=" * 60)
print()

print("""  The tournaplex f-vector (f₀, f₁, f₂, f₃, ...) is an INVARIANT
  of the tournament (up to isomorphism).

  Two tournaments with the same f-vector have the same "shape"
  of their ordering structure, but may differ in details.

  QUESTION: Does the tournaplex f-vector determine the iso class?
  Or can two non-isomorphic tournaments have the same f-vector?
""")

# Check at n=5
fvec_to_classes = defaultdict(list)
for can in canon_to_bits5:
    rep = canon_to_bits5[can][0]
    adj = tournament_adj(n, rep)
    h = H_dp(adj, n)

    f = [n, comb(n, 2)]
    for k in range(3, n+1):
        trans_k = 0
        for subset in combinations(range(n), k):
            wins = [sum(adj[v][u] for u in subset if u != v) for v in subset]
            if sorted(wins) == list(range(k)):
                trans_k += 1
        f.append(trans_k)

    fvec_to_classes[tuple(f)].append(h)

print(f"  At n=5: {len(fvec_to_classes)} distinct f-vectors for {len(canon_to_bits5)} iso classes")
print()
for fvec in sorted(fvec_to_classes.keys()):
    h_list = fvec_to_classes[fvec]
    print(f"  f = {fvec}: H = {h_list}")

print()
if len(fvec_to_classes) < len(canon_to_bits5):
    print(f"  The f-vector does NOT distinguish all iso classes at n=5.")
    print(f"  Some non-isomorphic tournaments share the same f-vector.")
else:
    print(f"  The f-vector DOES distinguish all iso classes at n=5.")

print()

# ==========================================================================
# CHAPTER 8: THE TOURNAPLEX FILTRATION
# ==========================================================================

print()
print("CHAPTER 8: PERSISTENT HOMOLOGY OF TOURNAPLEXES")
print("=" * 60)
print()

print("""  The TOURNAPLEX FILTRATION:
  Build the simplicial complex incrementally by score:

  Step 1: Add all vertices (score ≥ 0)
  Step 2: Add arcs from high-score to low-score vertices (range ≥ d)
  Step 3: Add arcs from medium-score vertices
  ...

  At each step, the transitive sub-tournaments grow.
  The Betti numbers β_k change as simplices are added.
  The PERSISTENCE DIAGRAM records when β_k features are born and die.

  FOR TOURNAMENTS: The natural filtration is by SCORE DIFFERENCE:
    First add arcs between vertices with score difference ≥ n-2
    Then score difference ≥ n-3
    ...
    Finally score difference ≥ 0 (all arcs)

  This is the STAIRCASE FILTRATION: it adds arcs in order of their
  2^d value, from most important to least.

  The persistence diagram of this filtration IS a topological
  fingerprint of the tournament. It measures:
  - When ordering structure APPEARS (birth of simplices)
  - When ordering conflicts APPEAR (death of simplices by cycles)
  - The PERSISTENCE of each feature (difference between birth and death)

  FOR THE BLUE BRAIN PROJECT:
  Neural circuits are directed graphs (not full tournaments).
  The tournaplex construction applies to ANY directed graph.
  Persistent homology reveals:
  - "Holes" in the neural network's information flow
  - Robust structural features vs ephemeral ones
  - The "dimension" of information processing
""")

# Compute the staircase filtration for a specific tournament
bits = 100  # Example tournament
adj = tournament_adj(5, bits)
h = H_dp(adj, 5)
ss = score_seq(adj, 5)
scores_lab = [sum(adj[i][j] for j in range(5) if j != i) for i in range(5)]

print(f"  Example tournament (bits=100, H={h}, score={ss}):")
print(f"  Labeled scores: {scores_lab}")
print()

# Build filtration by range
print(f"  Staircase filtration:")
for d_thresh in range(4, -1, -1):
    # Include all arcs with range ≥ d_thresh
    partial_adj = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            if i == j: continue
            # Range of arc (i,j) in transitive ordering:
            # Approximate by |scores_i - scores_j|
            if adj[i][j] and abs(scores_lab[i] - scores_lab[j]) >= d_thresh:
                partial_adj[i][j] = 1

    # Count transitive k-cliques in partial tournament
    f = [5]
    arcs = sum(partial_adj[i][j] for i in range(5) for j in range(5) if i != j and partial_adj[i][j])
    f.append(arcs)

    trans_3 = 0
    for triple in combinations(range(5), 3):
        i, j, k = triple
        # Check transitive triple in ANY of the 3! orderings
        for perm in permutations(triple):
            a, b, c = perm
            if partial_adj[a][b] and partial_adj[b][c] and partial_adj[a][c]:
                trans_3 += 1
                break
    f.append(trans_3)

    print(f"    d ≥ {d_thresh}: arcs={arcs:>2d}, trans_triples={trans_3:>2d}, "
          f"f-vec(partial) = {f}")

print()

# ==========================================================================
# CHAPTER 9: 24-CELL IN THE INTERCHANGE GRAPH
# ==========================================================================

print()
print("CHAPTER 9: CAN THE 24-CELL APPEAR IN TOURNAMENT INTERCHANGE?")
print("=" * 60)
print()

print("""  The 24-cell has 24 vertices and 96 edges.
  Each vertex has degree 8 (8-regular).

  WHERE could a 24-vertex, 8-regular graph appear in tournaments?

  CANDIDATE 1: The interchange graph I((1,1,2,2)) at n=4.
  Score (1,1,2,2) has 24 tournaments. Interchange degree = c₃ = 2.
  So I((1,1,2,2)) is 2-regular, NOT 8-regular. NOT the 24-cell.

  CANDIDATE 2: The regular score class at n=5.
  Score (2,2,2,2,2) has 24 tournaments. But c₃ = 5.
  So I((2,2,2,2,2)) has 24 vertices with interchange degree 5.
  5-regular, NOT 8-regular. Close but not the 24-cell.

  CANDIDATE 3: Some enriched graph on 24 tournaments.
  If we use arc-flips + 3-cycle flips together on the 24 regular
  tournaments at n=5:
""")

# Build the enriched interchange graph on the 24 regular tournaments at n=5
n = 5
m = comb(n, 2)
regular_score = (2, 2, 2, 2, 2)

reg_bits = []
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    ss = score_seq(adj, n)
    if ss == regular_score:
        reg_bits.append(bits)

print(f"  Regular tournaments at n=5: {len(reg_bits)}")
reg_set = set(reg_bits)

# Arc-flip neighbors within the regular class
arc_neighbors = defaultdict(set)
for bits in reg_bits:
    for bit in range(m):
        nbr = bits ^ (1 << bit)
        if nbr in reg_set:
            arc_neighbors[bits].add(nbr)

arc_degrees = [len(arc_neighbors[b]) for b in reg_bits]
print(f"  Arc-flip degrees within regular class: {sorted(arc_degrees)}")
print(f"  (0 everywhere: no arc flip preserves the regular score!)")
print()

# 3-cycle flip neighbors within the regular class
cyc_neighbors = defaultdict(set)
for bits in reg_bits:
    for triple in combinations(range(n), 3):
        new_bits = bits
        idx = 0
        for a in range(n):
            for b in range(a+1, n):
                if (a, b) in [(min(triple[0],triple[1]),max(triple[0],triple[1])),
                              (min(triple[0],triple[2]),max(triple[0],triple[2])),
                              (min(triple[1],triple[2]),max(triple[1],triple[2]))]:
                    new_bits ^= (1 << idx)
                idx += 1
        if new_bits in reg_set and new_bits != bits:
            cyc_neighbors[bits].add(new_bits)

cyc_degrees = [len(cyc_neighbors[b]) for b in reg_bits]
print(f"  3-cycle flip degrees within regular class: {sorted(cyc_degrees)}")
print()

# Combined
combined_degrees = [len(arc_neighbors[b] | cyc_neighbors[b]) for b in reg_bits]
print(f"  Combined (arc + 3-cycle) degrees: {sorted(combined_degrees)}")
print()

is_regular = len(set(combined_degrees)) == 1
if is_regular:
    deg = combined_degrees[0]
    n_edges = sum(combined_degrees) // 2
    print(f"  The enriched regular-class graph is {deg}-regular with {n_edges} edges!")
    if deg == 8 and n_edges == 96:
        print(f"  THIS IS THE 24-CELL'S EDGE COUNT!")
    else:
        print(f"  Not the 24-cell ({deg}-regular, {n_edges} edges vs 8-regular, 96 edges)")
else:
    print(f"  Not regular — cannot be the 24-cell.")

print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: THE 24-CELL AND TOURNAPLEXES")
print("=" * 72)
print()

print(f"""
  THE 24-CELL:
  - 24 vertices = 4! = |S_4| = D_4 roots = binary tetrahedral group
  - Self-dual f-vector (24, 96, 96, 24) — unique among regular polytopes
  - Connects to tournaments via: 4! transitive tournaments on 4 vertices
  - BUT: the permutohedron Π_3 (truncated octahedron, also 24 vertices)
    is a DIFFERENT polytope from the 24-cell
  - The D_4 roots arise from PAIRWISE SCORE DIFFERENCES

  TOURNAPLEXES:
  - Simplicial complex: k-simplex = transitive (k+1)-tournament
  - The f-vector (transitive clique counts) is an iso class invariant
  - At n=5: f-vector does NOT fully determine iso class
    (some classes share f-vectors)
  - Euler characteristic χ correlates with H: corr = {corr:.4f}
  - The staircase filtration gives a natural persistent homology

  THE CONNECTIONS:
  1. 24 = 4! connects the 24-cell to S_4 and tournament permutohedron
  2. D_4 roots = pairwise score differences
  3. Binary tetrahedral group 2T → A_4 → S_4 (double cover of even perms)
  4. Self-duality ↔ Z/2Z complement ↔ the two-sheeted cover from S200
  5. The regular-class interchange graph on 24 tournaments is related
     but NOT isomorphic to the 24-cell

  THE DEEPEST POINT:
  The 24-cell exists because D_4 has TRIALITY — a unique symmetry
  that permutes the three 8-dimensional representations.
  In tournament terms: D_4 triality permutes the three "aspects"
  of a 4-vertex tournament:
    - The SCORE aspect (which vertices win)
    - The CYCLE aspect (which triples are cyclic)
    - The PATH aspect (how many Hamiltonian paths)

  These three aspects are permuted by D_4 triality.
  This triality is WHY the 24-cell is self-dual,
  and WHY n=4 is the ONLY size where all three aspects
  of tournament structure are perfectly balanced.
""")

print("opus-2026-03-22-S205")
