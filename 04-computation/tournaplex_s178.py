#!/usr/bin/env python3
"""
tournaplex_s178.py -- opus-2026-03-22-S178

TOURNAPLEXES AND TOURNAMENT THEORY

A TOURNAPLEX (Lupo, Reimann et al., 2020) is a semi-simplicial complex
whose simplices are tournaments:
  - An n-simplex = a tournament on (n+1) vertices
  - Face map d_i = delete vertex i (giving a tournament on n vertices)
  - The chain complex and boundary operator are the usual semi-simplicial ones

THIS IS EXACTLY OUR VERTEX DELETION!
  From S171: each class at n maps to classes at n-1 via vertex deletion.
  The deletion profile = the tournaplex face map applied to iso classes.

THE KEY INSIGHT: Our entire project has been studying the tournaplex
of ALL tournaments, denoted T_n. We just didn't call it that.

CONNECTIONS:
  1. GLMY path homology (S38-S41) = the path complex of a single tournament
  2. The tournaplex = the simplicial complex of ALL tournaments
  3. beta_2 = 0 (our theorem) = a statement about the tournaplex homology
  4. The vertex deletion profiles = face maps of the tournaplex
  5. G_n = the 1-skeleton of the tournaplex (after quotienting by S_n)

This script:
  1. Constructs the tournaplex at small n
  2. Computes its chain complex and homology
  3. Connects to our existing results
  4. Explores the directionality filtration
  5. Finds what the tournaplex Betti numbers mean

Author: opus-2026-03-22-S178
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

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

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def delete_vertex(A, n, v):
    """Face map d_v: delete vertex v from tournament A."""
    remaining = [i for i in range(n) if i != v]
    A_del = [[A[remaining[i]][remaining[j]] for j in range(n-1)] for i in range(n-1)]
    return A_del

# ============================================================
# PART 1: THE TOURNAPLEX DEFINITION
# ============================================================

banner("PART 1: WHAT IS A TOURNAPLEX?")

print("""
A TOURNAPLEX is a semi-simplicial complex X where:
  - X_n = set of (n+1)-vertex tournaments (the n-simplices)
  - Face maps d_i: X_n -> X_{n-1} delete vertex i
  - d_i d_j = d_{j-1} d_i for i < j (semi-simplicial identity)

The FULL TOURNAPLEX T contains ALL tournaments at ALL sizes:
  T_0 = {single vertex} (1 element)
  T_1 = {2-vertex tournaments} = {A->B, B->A} (2 elements)
  T_2 = {3-vertex tournaments} (8 elements)
  T_3 = {4-vertex tournaments} (64 elements)
  T_n = {(n+1)-vertex tournaments} (2^{C(n+1,2)} elements)

The CHAIN COMPLEX:
  C_n = free abelian group on T_n
  d = sum_{i=0}^n (-1)^i d_i : C_n -> C_{n-1}
  d^2 = 0 (from the semi-simplicial identity)

The HOMOLOGY:
  H_n(T) = ker(d_n) / im(d_{n+1})
""")

# ============================================================
# PART 2: COMPUTE THE TOURNAPLEX AT SMALL n
# ============================================================

banner("PART 2: THE TOURNAPLEX CHAIN COMPLEX")

# At each level, enumerate tournaments and compute the boundary map
for N in [3, 4]:  # N = number of vertices = (n+1) for n-simplices
    n = N - 1  # simplex dimension
    print("  T_%d: %d-simplices = tournaments on %d vertices" % (n, n, N))
    print("  |T_%d| = 2^C(%d,2) = %d\n" % (n, N, 2**comb(N,2)))

    # Enumerate all tournaments at this level
    tourns_N = list(all_tournaments(N))
    print("  Number of %d-vertex tournaments: %d" % (N, len(tourns_N)))

    # For each tournament, compute all face maps
    # d_i deletes vertex i, giving a tournament on N-1 vertices
    if N >= 3:
        tourns_Nm1 = list(all_tournaments(N-1))
        # Build index for (N-1)-vertex tournaments
        idx_Nm1 = {}
        for j, A in enumerate(tourns_Nm1):
            key = tuple(A[r][c] for r in range(N-1) for c in range(N-1))
            idx_Nm1[key] = j

        # Compute boundary matrix
        # d(sigma) = sum_{i=0}^{n} (-1)^i d_i(sigma)
        # The boundary matrix B has rows indexed by T_{n-1} and cols by T_n
        nrow = len(tourns_Nm1)
        ncol = len(tourns_N)
        B = np.zeros((nrow, ncol), dtype=int)

        for j, A in enumerate(tourns_N):
            for v in range(N):
                A_del = delete_vertex(A, N, v)
                key = tuple(A_del[r][c] for r in range(N-1) for c in range(N-1))
                idx = idx_Nm1[key]
                sign = (-1)**v
                B[idx][j] += sign

        # Rank of boundary matrix
        rank_B = np.linalg.matrix_rank(B)
        print("  Boundary matrix d_%d: %d x %d, rank = %d" %
              (n, nrow, ncol, rank_B))
        print("  ker(d_%d) dimension = %d - %d = %d" %
              (n, ncol, rank_B, ncol - rank_B))
        print()

# ============================================================
# PART 3: TOURNAPLEX ON ISO CLASSES
# ============================================================

banner("PART 3: TOURNAPLEX ON ISO CLASSES (THE QUOTIENT)")

print("""
  The S_n action on tournaments gives the QUOTIENT tournaplex T/S_n.
  Its simplices are isomorphism classes, and face maps are the
  DELETION PROFILES from S171.

  The quotient T/S_n is NOT a semi-simplicial set (face maps may
  not be well-defined on classes since deleting different vertices
  of the same iso class may give different iso classes at N-1).

  But we CAN define a WEIGHTED version:
  The face map d_i on a class C produces a MULTISET of classes at N-1.
  The multiplicity = number of vertices in C whose deletion gives
  the target class.
""")

# Build the class-level face maps at N=5
N = 5; n = N - 1
print("  N=%d (dimension %d):\n" % (N, n))

# Build class catalogs for N and N-1
cmap = {}; cdata = {}; creps = {}
for A in all_tournaments(N):
    c = canonical(A, N)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid
        cdata[cid] = {'H': count_hp(A, N), 'size': 0}
        creps[cid] = [row[:] for row in A]
    cdata[cmap[c]]['size'] += 1
nc_N = len(cdata)

cmap4 = {}; cdata4 = {}
for A in all_tournaments(N-1):
    c = canonical(A, N-1)
    if c not in cmap4:
        cid = len(cmap4); cmap4[c] = cid
        cdata4[cid] = {'H': count_hp(A, N-1), 'size': 0}
    cdata4[cmap4[c]]['size'] += 1
nc_Nm1 = len(cdata4)

print("  %d iso classes at N=%d, %d iso classes at N=%d\n" %
      (nc_N, N, nc_Nm1, N-1))

# For each class at N, compute the deletion profile
print("  Deletion profiles (class at N=5 -> multiset of classes at N=4):")
del_profiles = {}
for cid in range(nc_N):
    A = creps[cid]
    profile = Counter()
    for v in range(N):
        A_del = delete_vertex(A, N, v)
        c_del = cmap4[canonical(A_del, N-1)]
        profile[c_del] += 1
    del_profiles[cid] = profile
    h = cdata[cid]['H']
    prof_str = ", ".join("%d:%d" % (c, cnt) for c, cnt in sorted(profile.items()))
    child_H = [cdata4[c]['H'] for c in sorted(profile.keys())]
    print("    class %2d (H=%2d): {%s} -> child H=%s" %
          (cid, h, prof_str, child_H))

# ============================================================
# PART 4: TOURNAPLEX BOUNDARY ON ISO CLASSES
# ============================================================

banner("PART 4: ISO CLASS BOUNDARY MAP")

print("  The boundary d(class C) = sum_{v=0}^{N-1} (-1)^v d_v(C)")
print("  On iso classes, this gives a SIGNED sum of child classes.\n")

print("  BUT: the sign (-1)^v depends on the VERTEX LABELING.")
print("  Different representatives of the same class may give")
print("  different signs. So the boundary is NOT well-defined on classes.\n")

print("  HOWEVER: the UNSIGNED version (absolute values) IS well-defined.")
print("  The unsigned boundary d_abs(C) = sum of multiplicities.\n")

print("  The TOTAL multiplicity sum = N (one deletion per vertex).\n")

for cid in range(nc_N):
    h = cdata[cid]['H']
    total_mult = sum(del_profiles[cid].values())
    print("    class %d (H=%d): total multiplicity = %d (should be %d)" %
          (cid, h, total_mult, N))

# ============================================================
# PART 5: CONNECTION TO GLMY PATH HOMOLOGY
# ============================================================

banner("PART 5: TOURNAPLEX vs GLMY PATH HOMOLOGY")

print("""
  GLMY path homology (Grigoryan-Lin-Muranov-Yau) is defined for
  a SINGLE digraph G. It constructs a chain complex from ALLOWED
  PATHS in G and computes homology groups H_p(G).

  OUR RESULTS (S38-S41):
  - beta_2 = 0 for ALL tournaments (verified n=3-8)
  - beta_1 * beta_3 = 0 (seesaw)
  - beta_3 onset at n=6

  The TOURNAPLEX is DIFFERENT from GLMY:
  - GLMY: one fixed digraph G, paths of various lengths
  - Tournaplex: the SPACE of all tournaments, with face maps = deletion

  But they're RELATED:
  - A tournament T on n vertices is an (n-1)-simplex in the tournaplex
  - The GLMY homology of T measures the "internal" structure
  - The tournaplex homology measures the "external" structure
    (how tournaments relate to each other via vertex deletion)

  THE KEY CONNECTION:
  beta_2(T) = 0 for all tournaments means: every 2-cycle in the
  GLMY path complex of T is a boundary. This is a statement about
  the INTERNAL structure of each simplex.

  The tournaplex boundary d_i deletes a vertex, which changes the
  INTERNAL structure. The tournaplex homology measures when these
  changes "cancel" across the deleted vertex positions.

  CONJECTURE: The tournaplex homology H_n(T) is related to the
  AVERAGE GLMY homology across all tournaments:
  H_n(tournaplex) ~ <H_n(GLMY)> averaged over tournaments.
""")

# ============================================================
# PART 6: THE DIRECTIONALITY FILTRATION
# ============================================================

banner("PART 6: DIRECTIONALITY FILTRATION")

print("""
  The DIRECTIONALITY FILTRATION (Lupo-Reimann) filters tournaments
  by how "transitive" they are.

  For a tournament T on vertices {v_0, ..., v_n}, the DIRECTIONALITY
  is the number of arcs that go "upward" (from lower to higher index)
  in the vertex ordering.

  A TRANSITIVE tournament has all n(n-1)/2 arcs going upward.
  A tournament with d arcs going upward has directionality d.

  The filtration: T_d = tournaments with directionality >= d.
  T_0 = all tournaments (no constraint)
  T_{C(n,2)} = transitive tournament only

  PERSISTENT HOMOLOGY of this filtration captures how topological
  features (homology classes) appear and disappear as we restrict
  to more transitive tournaments.

  CONNECTION TO OUR WORK:
  Directionality = C(n,2) - (number of "upsets" from the natural order)
  = C(n,2) - (Kemeny distance from transitive)

  The Kemeny distance IS the minimal number of arc reversals to reach
  the transitive tournament. This connects directly to G_n!

  In G_n: the distance between a class and the transitive class
  = minimum number of arc reversals = Kemeny distance (up to iso).
""")

# Compute directionality distribution at N=5
N = 5
m = comb(N, 2)
dir_dist = Counter()
for A in all_tournaments(N):
    # Directionality = arcs going i -> j with i < j
    d = sum(A[i][j] for i in range(N) for j in range(i+1, N))
    dir_dist[d] += 1

print("  Directionality distribution at N=%d:" % N)
for d in sorted(dir_dist.keys()):
    print("    d=%2d: %d tournaments (%.1f%%)" %
          (d, dir_dist[d], 100*dir_dist[d]/2**m))

print()
print("  Note: d and C(n,2)-d have the same count (complement symmetry).")
print("  d = C(n,2)/2 = %d is the 'equator' (maximum count).\n" % (m//2))

# ============================================================
# PART 7: THE TOURNAPLEX AND THE STAIRCASE
# ============================================================

banner("PART 7: TOURNAPLEX MEETS THE STAIRCASE")

print("""
  From S164: the staircase delta_{n-2} encodes the tournament.
  The tournaplex face map d_v removes a vertex.

  In the staircase: deleting vertex v removes:
  - All arcs incident to v = one row + one column of the staircase
  - The remaining staircase is delta_{n-3} = the inner tournament

  This is EXACTLY the source-sink decomposition from S171!
  Deleting the source (highest score) or sink (lowest score)
  = removing the L-border of the staircase.

  THE TOURNAPLEX FACE MAP d_v IS THE STAIRCASE BORDER REMOVAL.

  The chain complex of the tournaplex, restricted to the staircase,
  gives the HOMOLOGY OF THE STAIRCASE GROWTH SEQUENCE:
    delta_1 -> delta_2 -> delta_3 -> ...
  where each step adds an L-border.

  The TOURNAPLEX HOMOLOGY measures when these border additions
  "cancel" in a chain-complex sense.
""")

# ============================================================
# PART 8: THE TOURNAPLEX AND THE HALF-SIMPLEX
# ============================================================

banner("PART 8: THE TOURNAPLEX AND THE HALF-SIMPLEX (S177)")

print("""
  From S177: the tournament has fractional simplicial dimension 1/2.
  The a-simplicial number S_{1/2}(k) = C(2k,k)/4^k is the fiber fraction.

  In a standard simplicial complex, an n-simplex has n+1 vertices
  and the face maps d_i reduce the dimension by 1.

  In a TOURNAPLEX, the "simplex" is a tournament on n+1 vertices.
  But the tournament is MORE than just a simplex — it has a
  DIRECTION on every edge. The directionality adds 1/2 dimension.

  A standard (n-1)-simplex has volume 1/n!.
  A tournament on n vertices has "fiber fraction" S_{1/2}(n-2)
  = C(2(n-2), n-2) / 4^{n-2} ~ 1/sqrt(pi*(n-2)).

  The RATIO:
  fiber_fraction / simplex_volume = S_{1/2}(n-2) * (n-1)!
  This is a LARGE number, growing rapidly.

  INTERPRETATION: Each simplex in the tournaplex carries MUCH MORE
  information than a standard simplex. The extra information is the
  directionality — the tournament structure on the edges.

  The tournaplex is a "FAT" simplicial complex: each simplex is
  inflated by the tournament's directional information.
  The inflation factor = 2^{C(n,2)} / (number of standard simplices)
  = the number of possible directions per simplex.
""")

for n in range(3, 8):
    N = n
    m = comb(N, 2)
    n_tournaments = 2**m
    n_orderings = factorial(N)
    inflation = n_tournaments // n_orderings if n_orderings > 0 else 0
    print("  N=%d: 2^C(N,2) = %d tournaments, N! = %d orderings, ratio = %.2f" %
          (N, n_tournaments, n_orderings, n_tournaments/n_orderings))

# ============================================================
# PART 9: TOURNAPLEX BETTI NUMBERS
# ============================================================

banner("PART 9: TOURNAPLEX BETTI NUMBERS (small computation)")

print("  Computing H_*(T) for the full tournaplex up to N=4.\n")

# N=2 (T_0 level): 2 tournaments on 2 vertices
# N=3 (T_1 level): 8 tournaments on 3 vertices
# Boundary d_1: T_1 -> T_0

# T_0: tournaments on 2 vertices = {A->B, B->A} = 2 elements
t0 = [(0,), (1,)]  # 0 = A->B, 1 = B->A

# T_1: tournaments on 3 vertices = 8 elements
t1 = []
for A in all_tournaments(3):
    key = tuple(A[i][j] for i in range(3) for j in range(3) if i != j)
    t1.append(key)

# Build boundary map d_1: C_1 -> C_0
# d_1(sigma) = d_0(sigma) - d_1(sigma) + d_2(sigma)
# d_i = delete vertex i

# Index T_0 (2-vertex tournaments)
t0_list = list(all_tournaments(2))
t0_idx = {}
for i, A in enumerate(t0_list):
    key = tuple(A[r][c] for r in range(2) for c in range(2))
    t0_idx[key] = i

# Build boundary
B1 = np.zeros((len(t0_list), len(t1)), dtype=int)
t1_list = list(all_tournaments(3))
for j, A in enumerate(t1_list):
    for v in range(3):
        A_del = delete_vertex(A, 3, v)
        key_del = tuple(A_del[r][c] for r in range(2) for c in range(2))
        idx = t0_idx[key_del]
        B1[idx][j] += (-1)**v

rank_B1 = np.linalg.matrix_rank(B1)
ker_B1 = B1.shape[1] - rank_B1
print("  d_1 (3-vert -> 2-vert): %d x %d, rank = %d, ker = %d" %
      (B1.shape[0], B1.shape[1], rank_B1, ker_B1))

# For H_1, we need im(d_2): C_2 -> C_1
# C_2 = tournaments on 4 vertices
t2_list = list(all_tournaments(4))
t1_idx = {}
for i, A in enumerate(t1_list):
    key = tuple(A[r][c] for r in range(3) for c in range(3))
    t1_idx[key] = i

B2 = np.zeros((len(t1_list), len(t2_list)), dtype=int)
for j, A in enumerate(t2_list):
    for v in range(4):
        A_del = delete_vertex(A, 4, v)
        key_del = tuple(A_del[r][c] for r in range(3) for c in range(3))
        idx = t1_idx[key_del]
        B2[idx][j] += (-1)**v

rank_B2 = np.linalg.matrix_rank(B2)
ker_B2 = B2.shape[1] - rank_B2
print("  d_2 (4-vert -> 3-vert): %d x %d, rank = %d, ker = %d" %
      (B2.shape[0], B2.shape[1], rank_B2, ker_B2))

# H_0 = ker(d_0) / im(d_1), but d_0 = 0 (no lower level)
# So H_0 = C_0 / im(d_1) = dim(C_0) - rank(d_1)
H0 = len(t0_list) - rank_B1
print("\n  beta_0(T) = %d (connected components)" % H0)

# H_1 = ker(d_1) / im(d_2)
# ker(d_1) = B1.shape[1] - rank_B1 = 8 - rank_B1
# im(d_2) = rank_B2
H1 = ker_B1 - rank_B2
print("  beta_1(T) = ker(d_1) - rank(d_2) = %d - %d = %d" %
      (ker_B1, rank_B2, H1))

# Verify d^2 = 0
print("\n  Checking d^2 = 0 (B1 @ B2 = 0):")
product = B1 @ B2
print("  max|B1 @ B2| = %d (should be 0)" % abs(product).max())

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE TOURNAPLEX IS OUR PROJECT")

print("""
THE TOURNAPLEX CONNECTS EVERYTHING:

1. DEFINITION: A tournaplex is a semi-simplicial complex of tournaments.
   Face maps = vertex deletion. Our entire project studies the
   full tournaplex T = {all tournaments at all sizes}.

2. CHAIN COMPLEX: d = sum (-1)^i d_i gives a boundary operator.
   The tournaplex homology H_*(T) measures "cancellations" in the
   vertex deletion structure. We computed beta_0 = %d, beta_1 = %d.

3. CONNECTION TO GLMY: GLMY path homology = internal structure
   of each tournament simplex. Tournaplex homology = external
   structure of how tournaments relate via deletion.
   Our beta_2(T) = 0 is an INTERNAL property of tournament simplices.

4. CONNECTION TO G_n: The iso class quotient T/S_n gives
   a "quotient tournaplex" whose 1-skeleton IS G_n.
   The deletion profiles are the face maps.
   The fiber bundle G_n -> G_{n-2} is the tournaplex d_source + d_sink.

5. CONNECTION TO STAIRCASE: Vertex deletion in the tournaplex
   = border removal in the staircase delta_{n-2} -> delta_{n-3}.
   The tournaplex chain complex on staircases gives the homology
   of the staircase growth sequence.

6. DIRECTIONALITY FILTRATION: Filtering by "how transitive"
   gives persistent homology. The Kemeny distance to the transitive
   tournament = the filtration parameter. This connects to G_n
   where distance to the transitive class = Kemeny distance.

7. THE HALF-SIMPLEX: Each tournaplex simplex carries 2^{C(n,2)}
   orientations. The "effective dimension" is 1/2 (from S177).
   The tournaplex is a FAT simplicial complex: more structure
   per simplex than standard simplicial complexes.

8. NEUROSCIENCE: The Blue Brain Project uses flag tournaplexes
   to analyze neural connectivity. Our tournament tools (H count,
   OCR, spectral gap, fiber fractions) are directly applicable
   to this analysis pipeline.

THE DEEP POINT:
The tournaplex is the natural HOME for all our results.
  - OCF: H = I(Omega, 2) is a property of individual simplices
  - G_n: the quotient tournaplex modulo relabeling
  - Staircase: the fundamental domain of S_n in the tournaplex
  - GLMY: the internal homology of each simplex
  - Fiber fraction: the coset size in the tournaplex tower
  - Half-simplex: the effective dimensionality per simplex

The tournaplex IS the "meta-tournament" we've been building.
""" % (H0, H1))

print("Done. opus-2026-03-22-S178")
