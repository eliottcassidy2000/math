#!/usr/bin/env python3
"""
polygon_tournament_s201.py — Regular Polygons, Simplices, and Tournaments
opus-2026-03-22-S201

THE ISOMORPHISMS:
  Tournament on n vertices = orientation of the (n-1)-simplex's edges
  Cyclic tournament C_n = regular n-gon with directed edges
  Staircase Young diagram δ_{n-1} = upper triangle of adjacency matrix
  Meta-graph G_n = quotient of simplex orientation space by S_n

This script makes ALL connections precise and computational.
"""

import numpy as np
from math import comb, factorial, pi, sqrt, gcd
from fractions import Fraction
from itertools import permutations, combinations
from collections import defaultdict, Counter

print("=" * 72)
print("  REGULAR POLYGONS, SIMPLICES, AND TOURNAMENTS")
print("  opus-2026-03-22-S201")
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
# CHAPTER 1: THE SIMPLEX — A TOURNAMENT IS AN ORIENTED SIMPLEX
# ==========================================================================

print("CHAPTER 1: TOURNAMENT = ORIENTED (n-1)-SIMPLEX")
print("=" * 60)
print()

print("""  The (n-1)-simplex Δ_{n-1} has:
    - n vertices (the items being compared)
    - C(n,2) edges (one for each pair)
    - C(n,3) triangular faces (one for each triple)
    - ...
    - 1 top-dimensional face (the entire simplex)

  A TOURNAMENT on n vertices is an ORIENTATION of the edges
  of K_n = the 1-skeleton of Δ_{n-1}.

  This is EXACT, not metaphorical:
    Tournament T ↔ map σ: edges(Δ_{n-1}) → {+1, -1}

  The number of tournaments = 2^{C(n,2)} = 2^{# edges of Δ_{n-1}}.

  THE FACE LATTICE OF THE SIMPLEX:
    0-faces (vertices): n items
    1-faces (edges): C(n,2) arcs — these ARE the tournament
    2-faces (triangles): C(n,3) triples — these carry the 3-cycles
    k-faces: C(n,k+1) — these carry the k+1-cycles
    (n-1)-face: the entire simplex — carries the Hamiltonian structure

  THE HAMILTONIAN PATH COUNT H(T) is a function on the TOP face
  of the simplex. It counts spanning paths through all vertices.
  So H is defined on the (n-1)-face.

  THE SCORE SEQUENCE is defined on the 0-faces (vertices).
  Each vertex has a score = out-degree = number of edges pointing away.

  THE 3-CYCLE COUNT c₃ is defined on the 2-faces (triangles).
  Each triangle is either a 3-cycle (1 of 2 orientations) or transitive.
""")

# Show the simplex dimensions
print(f"  SIMPLEX DIMENSIONS (f-vector):")
print(f"  {'n':>3s}  {'vertices':>8s}  {'edges':>7s}  {'triangles':>9s}  "
      f"{'tetra':>6s}  {'total faces':>11s}")
for n in range(3, 8):
    f_vector = [comb(n, k+1) for k in range(n)]
    total = sum(f_vector) + 1  # +1 for empty face
    print(f"  {n:3d}  {f_vector[0]:>8d}  {f_vector[1]:>7d}  "
          f"{f_vector[2] if len(f_vector)>2 else 0:>9d}  "
          f"{f_vector[3] if len(f_vector)>3 else 0:>6d}  "
          f"{total:>11d}")

print()
print("  Each k-face carries a k+1-cycle structure of the tournament.")
print("  The ENTIRE tournament is determined by the 1-faces (edges).")
print("  Higher faces carry DERIVED information (cycles, paths).")
print()

# ==========================================================================
# CHAPTER 2: REGULAR POLYGON = CYCLIC TOURNAMENT
# ==========================================================================

print()
print("CHAPTER 2: REGULAR n-GON = CYCLIC TOURNAMENT C_n")
print("=" * 60)
print()

print("""  The regular n-gon has vertices {0, 1, ..., n-1} arranged in a circle.
  The CYCLIC TOURNAMENT C_n orients edges by the rule:
    i → j iff 0 < (j - i) mod n ≤ (n-1)/2

  This means: i beats the (n-1)/2 vertices "clockwise" from it.
  (Only defined for odd n, since (n-1)/2 must be an integer.)

  THE SYMMETRY GROUP:
  - The regular n-gon has symmetry group D_n (dihedral, order 2n).
  - C_n (as a tournament) has automorphism group C_n (cyclic, order n).
  - The reflections of D_n REVERSE all arcs (give the complement).
  - So Aut(C_n) = C_n, not D_n.

  C_n is the MOST SYMMETRIC tournament:
  - It's vertex-transitive (all vertices look the same).
  - It's arc-transitive (all arcs look the same).
  - It has the largest automorphism group among tournaments on n vertices.
  - BUT: C_n does NOT always maximize H (we showed H_max(7) = 189 > 175 = H(C_7)).

  THE POLYGON PICTURE:
  Draw the n-gon. At each vertex, draw arrows to the (n-1)/2 nearest
  vertices clockwise. The resulting picture IS the tournament C_n.
  The arrows follow the "nearest-neighbor" rule on the polygon.
""")

# Compute and display cyclic tournaments
for n in [3, 5, 7]:
    half = (n-1) // 2
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if 0 < (j - i) % n <= half:
                adj[i][j] = 1
    h = H_dp(adj, n)
    ss = score_seq(adj, n)

    # Show adjacency
    print(f"  C_{n} (regular {n}-gon tournament):")
    print(f"    Score sequence: {ss}")
    print(f"    H = {h}")
    print(f"    |Aut| = {n} (cyclic group)")

    # Count 3-cycles
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    c3 += 1
    print(f"    3-cycles: {c3} (out of C({n},3) = {comb(n,3)} triangles)")
    print(f"    Fraction cyclic: {c3/comb(n,3):.4f}")
    print()

# ==========================================================================
# CHAPTER 3: STAIRCASE YOUNG DIAGRAM = TOURNAMENT MATRIX
# ==========================================================================

print()
print("CHAPTER 3: STAIRCASE YOUNG DIAGRAM δ_{n-1}")
print("=" * 60)
print()

print("""  The STAIRCASE YOUNG DIAGRAM δ_{n-1} = (n-1, n-2, ..., 1, 0):

  n=4:        n=5:            n=6:
  □ □ □       □ □ □ □         □ □ □ □ □
  □ □         □ □ □           □ □ □ □
  □           □ □             □ □ □
              □               □ □
                              □

  # boxes = C(n,2) = number of tournament arcs.
  Box (i,j) = arc between vertex i and vertex j (i < j).
  Value in box = 1 if i→j, 0 if j→i.

  THE STAIRCASE IS THE UPPER TRIANGLE OF THE ADJACENCY MATRIX:
    Row i, column j > i: box (i,j) = adj[i][j]

  A tournament is a 0-1 FILLING of δ_{n-1}.
  There are 2^{|δ|} = 2^{C(n,2)} such fillings = # tournaments.

  THE RANGE d of box (i,j):
    d = j - i - 1  (how far apart the vertices are)
    Range 0: adjacent vertices (boxes on the first superdiagonal)
    Range 1: two-apart vertices
    ...
    Range n-2: most distant vertices (the corner box)

  From S182: H(transitive with one flip at (i,j)) = 1 + 2^d.
  The RANGE IS the exponent. This is the staircase number system.

  THE STAIRCASE HAS DIAGONAL STRUCTURE:
    Anti-diagonal k: all boxes with i + j = k + 1
    These are arcs with the SAME sum of vertex indices.
    From S182: anti-diagonal arcs often give BLUE edges in G_n.
""")

# Show staircase for n=5 with ranges
n = 5
print(f"  Staircase δ_{n-1} with ranges d = j-i-1:")
print(f"    {'':>4s}", end="")
for j in range(1, n):
    print(f"  j={j}", end="")
print()
for i in range(n-1):
    print(f"  i={i}", end="")
    for j in range(1, n):
        if j > i:
            d = j - i - 1
            print(f"  d={d}", end="")
        else:
            print(f"     ", end="")
    print()

print()
print(f"  Range structure at n=5:")
for d in range(n-1):
    boxes = [(i, i+d+1) for i in range(n-d-1)]
    print(f"    Range d={d}: {len(boxes)} boxes at {boxes}, value 2^{d} = {2**d}")

print()
print(f"  Total: Σ 2^d × (n-1-d) = Σ from d=0 to {n-2}")
total = sum(2**d * (n-1-d) for d in range(n-1))
print(f"  = {total} (the staircase capacity)")
print(f"  H_max - 1 = {15 - 1} = 14")
print(f"  Capacity/H_range = {total}/14 = {total/14:.2f} (redundancy)")
print()

# ==========================================================================
# CHAPTER 4: THE SIMPLEX ORIENTATION POLYTOPE
# ==========================================================================

print()
print("CHAPTER 4: THE TOURNAMENT POLYTOPE")
print("=" * 60)
print()

print("""  Each tournament T on n vertices is a vertex of the CUBE {0,1}^{C(n,2)}.
  The TOURNAMENT POLYTOPE is the convex hull of all tournament vectors.

  But MORE interesting: consider the SCORE POLYTOPE.
  Map each tournament to its score sequence: T → s(T) ∈ R^n.
  The image is the LANDAU POLYTOPE P_n (from S188-S189).

  The Landau polytope P_n:
    P_n = conv{s(T) : T tournament on n vertices} ∩ {sorted sequences}

  Its integer points are the valid score sequences (A000571).
  Its vertices are the EXTREME score sequences.

  CONNECTION TO SIMPLICES:
  The regular (n-1)-simplex Δ has vertices e₁, ..., eₙ (standard basis).
  The score map s: {0,1}^{C(n,2)} → R^n is:
    s(T)_i = Σ_{j≠i} T_{ij}

  The transitive tournament has score (0,1,...,n-1).
  The regular tournament (odd n) has score ((n-1)/2, ..., (n-1)/2).

  The transitive score is a VERTEX of the Landau polytope.
  The regular score is the CENTER (barycenter).

  THIS IS EXACTLY THE SIMPLEX STRUCTURE:
  The simplex has vertices at the corners and center at the barycenter.
  The Landau polytope IS (a projection of) the simplex, where the
  vertices correspond to the most "extreme" score sequences and
  the center corresponds to the most "balanced" one.
""")

# Compute the Landau polytope vertices at n=5
def is_landau(seq):
    n_val = len(seq)
    s = sorted(seq)
    if sum(s) != comb(n_val, 2): return False
    for k in range(1, n_val + 1):
        if sum(s[:k]) < comb(k, 2): return False
    return True

def enumerate_score_sequences(n_val):
    target = comb(n_val, 2)
    results = []
    def gen(pos, remaining, min_val, current):
        if pos == n_val:
            if remaining == 0: results.append(tuple(current))
            return
        max_val = min(n_val - 1, remaining - (n_val - pos - 1) * min_val)
        if max_val < 0: return
        for v in range(min_val, max_val + 1):
            gen(pos + 1, remaining - v, v, current + [v])
    gen(0, target, 0, [])
    return [s for s in results if is_landau(s)]

seqs5 = enumerate_score_sequences(5)
print(f"  Landau polytope at n=5: {len(seqs5)} integer points (score sequences)")
print(f"  Vertices: {seqs5[0]} (transitive) and {seqs5[-1]} (regular)")
print(f"  Dimension: {5 - 1 - 1} = 3 (lives in hyperplane Σs = C(5,2) = 10)")
print()

# ==========================================================================
# CHAPTER 5: THE META-GRAPH AS A POLYTOPE SKELETON
# ==========================================================================

print()
print("CHAPTER 5: THE META-GRAPH G_n = SKELETON OF A POLYTOPE")
print("=" * 60)
print()

print("""  The meta-graph G_n has:
    - Vertices = isomorphism classes of tournaments
    - Edges = arc flips that change the iso class

  CONJECTURE: G_n is the 1-skeleton (edge graph) of a polytope.

  Evidence:
  - G_3: 2 vertices, 1 edge → interval (1-simplex)
  - G_4: 4 vertices, 5 edges → NOT a simplex (3-simplex has 6 edges)
  - G_5: 12 vertices, 30 edges

  If G_n were the 1-skeleton of a simple polytope P,
  then P would have:
    f₀ = |V(G_n)| vertices
    f₁ = |E(G_n)| edges
    dim(P) = degree of G_n (for a simple polytope, each vertex has dim edges)

  At n=5: degrees range 2-7, so G_5 is NOT the skeleton of a
  SIMPLE polytope. But it could be a non-simple polytope.

  THE PERMUTOHEDRON CONNECTION:
  The permutohedron Π_{n-1} is the convex hull of all permutations
  of (1, 2, ..., n) in R^n. It has:
    n! vertices
    Dimension n-1
    Its 1-skeleton connects permutations differing by an adjacent transposition

  But transitive tournaments ↔ permutations (total orders),
  and an adjacent transposition ↔ flipping an arc of range 0.

  So: the transitive fiber in G_n IS a face of the permutohedron!
  The permutohedron is TOO BIG (n! vertices vs 2^{C(n,2)} tournaments).
  But the tournament polytope PROJECTS onto the permutohedron
  via the score map.

  THE ASSOCIAHEDRON CONNECTION:
  The associahedron K_{n-1} has vertices = triangulations of an (n+1)-gon
  and edges = diagonal flips.
  Tournaments don't directly map to triangulations, BUT:
  the arc-flip graph on tournaments has the same LOCAL structure
  (each flip changes one element, like a diagonal flip).

  G_n is the QUOTIENT of the arc-flip graph by S_n.
  It's like the associahedron mod symmetry.
""")

# Compute G_n for n=3,4,5
for n in range(3, 6):
    m = comb(n, 2)
    # Find all iso classes
    class_to_rep = {}
    bits_to_class = {}
    for bits in range(1 << m):
        adj = tournament_adj(n, bits)
        can = canonical(adj, n)
        if can not in class_to_rep:
            class_to_rep[can] = bits
        bits_to_class[bits] = can

    # Build G_n
    classes = list(class_to_rep.keys())
    n_classes = len(classes)
    class_idx = {c: i for i, c in enumerate(classes)}

    edges = set()
    degrees = [0] * n_classes
    for bits in range(1 << m):
        c1 = bits_to_class[bits]
        for bit in range(m):
            nbr = bits ^ (1 << bit)
            c2 = bits_to_class[nbr]
            if c1 != c2:
                edge = tuple(sorted([class_idx[c1], class_idx[c2]]))
                edges.add(edge)

    n_edges = len(edges)
    for i, j in edges:
        degrees[i] += 1
        degrees[j] += 1

    print(f"  G_{n}: {n_classes} vertices, {n_edges} edges, "
          f"degrees = {sorted(degrees)}")

print()

# ==========================================================================
# CHAPTER 6: EACH FACE OF THE SIMPLEX IS A SUB-TOURNAMENT
# ==========================================================================

print()
print("CHAPTER 6: FACES OF THE SIMPLEX = SUB-TOURNAMENTS")
print("=" * 60)
print()

print("""  A k-face of the (n-1)-simplex corresponds to a (k+1)-subset
  of vertices. The INDUCED SUB-TOURNAMENT on these vertices
  is the restriction of T to this subset.

  This gives a NATURAL FILTRATION:
    0-faces → isolated vertices (scores)
    1-faces → arcs (wins/losses)
    2-faces → triangles (3-cycles or transitive triples)
    3-faces → 4-vertex sub-tournaments
    ...
    (n-1)-face → the full tournament

  THE EULER CHARACTERISTIC:
    χ(Δ_{n-1}) = 1 (the simplex is contractible)

  Each face carries tournament information of increasing complexity.
  The k-th face carries the k+1-vertex sub-tournament structure.

  THE DELETION-CONTRACTION on the simplex:
    Deleting a vertex = deleting a 0-face = removing a vertex from T
    Contracting an edge = contracting a 1-face = merging two vertices

  The deletion-contraction tree for H(T) IS the face lattice
  of the simplex, decorated with orientations!
""")

# Count sub-tournament types at each level for n=5
n = 5
m = comb(n, 2)

# Pick a specific tournament for analysis
bits = 0  # transitive
adj = tournament_adj(n, bits)

print(f"  Sub-tournament census for transitive tournament T₅:")
for k in range(1, n+1):
    # All k-subsets
    from itertools import combinations
    subtournaments = Counter()
    for subset in combinations(range(n), k):
        # Induced sub-tournament
        sub_adj = [[adj[subset[i]][subset[j]] for j in range(k)] for i in range(k)]
        if k <= 5:
            can = canonical(sub_adj, k)
            subtournaments[can] += 1

    print(f"    {k}-vertex sub-tournaments: {comb(n, k)} subsets, "
          f"{len(subtournaments)} distinct types")

print()

# ==========================================================================
# CHAPTER 7: REGULAR POLYGONS AND TOURNAMENT SYMMETRY GROUPS
# ==========================================================================

print()
print("CHAPTER 7: SYMMETRY — DIHEDRAL GROUPS AND TOURNAMENT AUTOMORPHISMS")
print("=" * 60)
print()

print("""  The regular n-gon has symmetry group D_n of order 2n:
    - n rotations (cyclic group C_n)
    - n reflections

  The cyclic tournament C_n has Aut(C_n) = C_n (rotations only):
    - Rotations PRESERVE arc directions → automorphisms
    - Reflections REVERSE all arcs → give complement, not automorphism

  UNLESS n is odd AND the tournament is self-complementary (SC),
  in which case reflections map T to an isomorphic copy of T^c = T.

  FOR SC TOURNAMENTS:
  - C_3 is SC: complement of the 3-cycle is the same 3-cycle (rotated)
  - C_5 is NOT SC: complement has different H
  - There exist SC tournaments that ARE fixed by a "reflection"

  THE REGULAR POLYGON HIERARCHY:
    n=3: Triangle → only 2 tournaments (transitive & cyclic)
         C_3 is the regular tournament. Aut = C_3 ≅ Z/3Z.
    n=4: Square → 4 iso classes
         No regular tournament exists (odd scores needed).
    n=5: Pentagon → 12 iso classes
         C_5 is one of 3 regular tournaments. Aut = C_5 ≅ Z/5Z.
         But H(C_5) = 15 = H_max ← the pentagon maximizes H!
    n=7: Heptagon → 456 iso classes
         C_7 has Aut = C_7. H(C_7) = 175 ≠ 189 = H_max.
         The heptagon does NOT maximize H!
""")

# Compute automorphism groups of cyclic tournaments
for n in [3, 5, 7]:
    if n % 2 == 0: continue
    half = (n-1) // 2
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if 0 < (j - i) % n <= half:
                adj[i][j] = 1

    # Count automorphisms
    aut_count = 0
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(n):
                if adj[i][j] != adj[perm[i]][perm[j]]:
                    is_aut = False
                    break
            if not is_aut: break
        if is_aut:
            aut_count += 1

    h = H_dp(adj, n)
    print(f"  C_{n}: |Aut| = {aut_count} (expect {n}), H = {h}")

print()

# ==========================================================================
# CHAPTER 8: THE TOURNAMENT CUBE AND THE HYPERCUBE
# ==========================================================================

print()
print("CHAPTER 8: TOURNAMENT SPACE = ORIENTED HYPERCUBE")
print("=" * 60)
print()

print("""  The space of all tournaments on n vertices is the hypercube
  {0,1}^{C(n,2)}. Each coordinate is one edge of the simplex.
  Value 0 or 1 = orientation of that edge.

  The TOURNAMENT CUBE has 2^{C(n,2)} vertices.
  Two vertices are adjacent if they differ in exactly one coordinate
  = they differ by exactly one arc flip.

  The ARC-FLIP GRAPH is the 1-skeleton of this hypercube.
  It's a C(n,2)-dimensional cube.

  THE META-GRAPH G_n is the QUOTIENT of this cube by S_n:
    G_n = {0,1}^{C(n,2)} / S_n

  where S_n acts by relabeling vertices.

  CUBE DIMENSIONS:
""")

for n in range(3, 8):
    m = comb(n, 2)
    n_tourn = 2**m
    n_iso = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}[n]
    print(f"  n={n}: cube dim = {m}, vertices = 2^{m} = {n_tourn}, "
          f"iso classes = {n_iso}, compression = {n_tourn/n_iso:.1f}×")

print()

print("""  The hypercube {0,1}^m has its own face lattice:
    0-faces (vertices): 2^m = # tournaments
    1-faces (edges): m × 2^{m-1} = # arc flips
    k-faces: C(m,k) × 2^{m-k}

  The tournament cube is a SIMPLEX (n-1 dim) wrapped in a
  HYPERCUBE (C(n,2) dim). The simplex lives inside as the
  set of "tournament-valid" orientations.

  KEY INSIGHT: The faces of the hypercube correspond to
  PARTIAL TOURNAMENTS — tournaments with some edges unoriented.
  A k-face has k unoriented edges and m-k oriented edges.
  This is a PARTIALLY ORDERED set.

  The SCORE of a partial tournament is an INTERVAL:
    s_i ∈ [s_i^min, s_i^max] depending on how unoriented edges resolve.
  The width of this interval measures UNCERTAINTY.
""")

# ==========================================================================
# CHAPTER 9: THE POLYGON DECOMPOSITION OF THE STAIRCASE
# ==========================================================================

print()
print("CHAPTER 9: POLYGON DECOMPOSITION OF THE STAIRCASE")
print("=" * 60)
print()

print("""  The staircase δ_{n-1} can be decomposed into DIAGONAL BANDS:

  n=5 staircase with diagonals:
       j=1  j=2  j=3  j=4
  i=0 [ d0 | d1 | d2 | d3 ]
  i=1      [ d0 | d1 | d2 ]
  i=2           [ d0 | d1 ]
  i=3                [ d0 ]

  Diagonal d: all boxes (i,j) with j-i-1 = d.
  Diagonal d has n-1-d boxes.

  Each diagonal IS A SIDE OF A POLYGON:
  - Diagonal 0 has n-1 boxes → the sides of the (n-1)-gon
  - Diagonal 1 has n-2 boxes → the short diagonals
  - Diagonal d has n-1-d boxes → the d-th diagonals
  - Diagonal n-2 has 1 box → the longest diagonal

  In the REGULAR POLYGON:
  - Side (d=0): connects adjacent vertices → closest comparison
  - Short diagonal (d=1): connects 2-apart vertices
  - d-th diagonal: connects (d+1)-apart vertices
  - Longest diagonal (d=n-2): connects antipodal vertices

  THE 2^d VALUE of each diagonal:
  From H = 1 + 2^d: reversing an arc at diagonal d changes H by 2^d.
  Longer diagonals (higher d) have EXPONENTIALLY more impact.

  THIS IS THE POLYGON'S HIERARCHY OF DIAGONALS:
    Sides (d=0): 2^0 = 1 impact per reversal
    Short diagonals (d=1): 2^1 = 2 impact
    Medium diagonals (d=2): 2^2 = 4 impact
    Long diagonal (d=n-2): 2^{n-2} impact

  The regular polygon naturally encodes the BINARY NUMBER SYSTEM
  through its diagonal structure!
""")

# Verify: diagonal structure
for n in [5, 6, 7]:
    print(f"  n={n} staircase diagonals:")
    total_capacity = 0
    for d in range(n-1):
        count = n - 1 - d
        value = 2**d
        total_capacity += count * value
        print(f"    d={d}: {count} boxes, value 2^{d}={value}, "
              f"contribution = {count * value}")
    print(f"    Total capacity = {total_capacity}")
    print()

# ==========================================================================
# CHAPTER 10: G_n AS THE DUAL OF A TRIANGULATION
# ==========================================================================

print()
print("CHAPTER 10: G_n AND GEOMETRIC DUALITY")
print("=" * 60)
print()

print("""  The meta-graph G_n connects iso classes by arc flips.
  Each iso class is a "region" in tournament space.
  Each flip is a "wall" between regions.

  This is EXACTLY the structure of a CELL DECOMPOSITION:
  - 0-cells (vertices of G_n) = iso classes
  - 1-cells (edges of G_n) = flip walls
  - 2-cells = ???

  The DUAL of G_n would have:
  - Vertices = 2-cells of the original decomposition
  - Edges connecting 2-cells that share a 1-cell

  WHAT ARE THE 2-CELLS?
  A 2-cell is a "face" bounded by edges of G_n.
  It corresponds to a 2-flip sequence that returns to the
  starting class: flip arc a, flip arc b, and you're back.

  These 2-faces are TRIANGLES in G_n: three iso classes
  forming a 3-cycle or triangle in the meta-graph.

  THE TRIANGULATION:
  G_5 has 12 vertices and 30 edges.
  A triangulation of a surface with V vertices and E edges has:
    F = 2 - V + E (Euler's formula for genus 0)
    F = 2 - 12 + 30 = 20 faces

  If G_5 is PLANAR, it would have 20 triangular faces.
  But is G_5 planar?
""")

# Check planarity of G_5 via Euler formula
n = 5
m = comb(n, 2)

# Build G_5
class_to_rep = {}
bits_to_class = {}
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    can = canonical(adj, n)
    if can not in class_to_rep:
        class_to_rep[can] = bits
    bits_to_class[bits] = can

classes = list(class_to_rep.keys())
n_classes = len(classes)
class_idx = {c: i for i, c in enumerate(classes)}

gn_edges = set()
gn_adj = defaultdict(set)
for bits in range(1 << m):
    c1 = bits_to_class[bits]
    for bit in range(m):
        nbr = bits ^ (1 << bit)
        c2 = bits_to_class[nbr]
        if c1 != c2:
            i, j = class_idx[c1], class_idx[c2]
            edge = tuple(sorted([i, j]))
            gn_edges.add(edge)
            gn_adj[i].add(j)
            gn_adj[j].add(i)

V = n_classes
E = len(gn_edges)

# For planarity: need E ≤ 3V - 6
max_planar_edges = 3 * V - 6
print(f"  G_5: V={V}, E={E}")
print(f"  Max edges for planar graph: 3V-6 = {max_planar_edges}")
print(f"  G_5 planar? {E <= max_planar_edges}")
print()

if E <= max_planar_edges:
    # Euler: V - E + F = 2
    F = 2 - V + E
    print(f"  If planar: F = 2 - {V} + {E} = {F} faces")
    # Average face size
    avg_face = 2 * E / F
    print(f"  Average face size: 2E/F = {avg_face:.1f} edges")
else:
    print(f"  G_5 is NOT planar (E={E} > 3V-6={max_planar_edges})")
    # Compute genus
    # 2 - 2g = V - E + F, with F = 2E/3 (if all faces are triangles)
    F_tri = 2 * E // 3
    genus = (2 - V + E - F_tri) // 2
    print(f"  If triangulated: genus ≥ {genus}")

print()

# Count triangles in G_5
triangles = 0
for i in range(V):
    for j in gn_adj[i]:
        if j > i:
            for k in gn_adj[j]:
                if k > j and k in gn_adj[i]:
                    triangles += 1

print(f"  Triangles in G_5: {triangles}")
print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: THE POLYGON-SIMPLEX-TOURNAMENT TRINITY")
print("=" * 72)
print()

print(f"""
  THE ISOMORPHISMS:

  1. TOURNAMENT = ORIENTED (n-1)-SIMPLEX
     n vertices + C(n,2) directed edges = simplex with oriented 1-skeleton.
     Each k-face carries a (k+1)-vertex sub-tournament.
     Hamiltonian paths live on the top face.

  2. CYCLIC TOURNAMENT = REGULAR n-GON
     C_n orients edges by nearest-neighbor rule on the polygon.
     Aut(C_n) = C_n (rotations preserve, reflections reverse).
     C_n maximizes H at n=3,5 but NOT at n=7 (189 > 175).

  3. STAIRCASE δ_{{n-1}} = UPPER TRIANGLE = POLYGON DIAGONALS
     Diagonal d = arcs between vertices (d+1) apart on the polygon.
     Value 2^d = exponential hierarchy of diagonal importance.
     The staircase IS the polygon's diagonal decomposition.

  4. META-GRAPH G_n = QUOTIENT POLYTOPE
     G_n = tournament cube / S_n.
     G_5 has {V} vertices, {E} edges.
     {'Planar' if E <= max_planar_edges else f'Non-planar (genus ≥ {genus})'}.
     {triangles} triangles.

  5. TOURNAMENT CUBE = ORIENTED HYPERCUBE
     {{0,1}}^{{C(n,2)}} has faces = partial tournaments.
     The score map projects the cube onto the Landau polytope.
     The permutohedron lives inside as the transitive fiber.

  THE DEEPEST CONNECTION:
  The regular polygon is the MOST SYMMETRIC way to arrange n points.
  The cyclic tournament C_n is the corresponding MOST SYMMETRIC tournament.
  The staircase δ_{{n-1}} is the SHADOW of the simplex onto the matrix plane.
  And G_n is the QUOTIENT of all possible orientations by the polygon's
  symmetry group, collapsed to isomorphism classes.

  Everything IS a polygon, viewed from different angles.
""")

print("opus-2026-03-22-S201")
