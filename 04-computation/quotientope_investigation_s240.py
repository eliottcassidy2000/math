#!/usr/bin/env python3
"""
quotientope_investigation_s240.py -- opus-2026-03-23-S240

THE QUOTIENTOPE CONNECTION

The arc-reversal graph on labeled tournaments is EXACTLY the m-cube Q_m
where m = C(n,2). Each tournament = a vertex of Q_m (binary string of
arc orientations), each arc reversal = flipping one bit = one cube edge.

The isomorphism class graph G_n = Q_m / S_n, the ORBIT GRAPH of the cube
under the symmetric group acting on pair coordinates.

The merged meta-graph G_n/Z_2 further quotients by complement (bit-flip all).

QUESTIONS:
1. Is G_n (or G_n/Z_2) the 1-skeleton of a polytope?
2. What is the clique complex? Its f-vector?
3. Is the H-ordering a lattice? What congruences does it have?
4. How does this relate to Pilaud-Santos quotientopes?

The permutohedron Π_n has n! vertices = transitive tournaments.
The m-cube has 2^m vertices = ALL tournaments.
The quotientope spectrum: Π_n (bottom) → associahedron → m-cube (top).
Our G_n lives "above" the cube quotient — it's the orbit graph.

KNOWN: Pilaud-Santos quotientopes come from lattice congruences of the
weak order on S_n. Our construction comes from the S_n action on Q_m.
These are DIFFERENT but may be related.

Author: opus-2026-03-23-S240
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

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

# ============================================================
banner("PART 1: THE CUBE FRAMEWORK")
# ============================================================

print("""
  THE KEY INSIGHT: G_n = Q_m / S_n

  The arc-reversal graph on LABELED tournaments is the m-cube Q_m:
  - Vertices: 2^m binary strings (tournament orientations)
  - Edges: pairs differing in exactly 1 bit (one arc reversal)
  - m = C(n,2) = number of edge pairs

  Q_m properties:
    Vertices: 2^m
    Edges: m * 2^{m-1}
    k-faces: C(m,k) * 2^{m-k}  (each k-face is a k-subcube)
    f-vector: (2^m, m*2^{m-1}, C(m,2)*2^{m-2}, ...)

  The isomorphism class graph G_n = Q_m / S_n:
    Vertices: 2^m / n! * (correction for non-trivial stabilizers) = A000568(n)
    Edges: E(n) = our sequence 1, 5, 30, 290, 4086, 91161

  S_n acts on Q_m by permuting the m coordinates (pair labels).
  The action preserves the cube structure.
""")

for n in range(3, 9):
    m = comb(n, 2)
    print("  n=%d: m=%d, |Q_m|=%d, m*2^{m-1}=%d, n!=%d, A000568=%s" %
          (n, m, 2**m, m*2**(m-1), factorial(n),
           {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}.get(n, '?')))

# ============================================================
banner("PART 2: CLIQUE COMPLEX OF G_n/Z_2")
# ============================================================

# For each n, build G_n/Z_2 and compute all cliques (complete subgraphs)
# This gives the f-vector of the clique complex (flag complex)

for n in [5, 6]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap = {}; sizes = {}; reps = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    # Info
    info = {}
    for cid in range(nc):
        A = reps[cid]
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c]}

    # Merged
    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # F matrix and edges
    F = np.zeros((nc, nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            F[ci][tc[bits^(1<<k)]] += 1

    # Merged adjacency
    adj = [[False]*nm for _ in range(nm)]
    for ci in range(nc):
        mi = mid_map[ci]
        for cj in range(nc):
            if F[ci][cj] > 0:
                mj = mid_map[cj]
                if mi != mj:
                    adj[mi][mj] = True
                    adj[mj][mi] = True

    ne = sum(sum(1 for j in range(i+1, nm) if adj[i][j]) for i in range(nm))

    # Find all cliques using Bron-Kerbosch
    all_cliques = []  # list of frozensets

    def bron_kerbosch(R, P, X):
        if not P and not X:
            all_cliques.append(frozenset(R))
            return
        # Choose pivot
        pivot = max(P | X, key=lambda v: sum(1 for w in P if adj[v][w]))
        for v in list(P - {w for w in range(nm) if adj[pivot][w]}):
            N_v = {w for w in range(nm) if adj[v][w]}
            bron_kerbosch(R | {v}, P & N_v, X & N_v)
            P = P - {v}
            X = X | {v}

    bron_kerbosch(set(), set(range(nm)), set())

    # f-vector of clique complex
    f_vec = Counter()
    for clique in all_cliques:
        for k in range(1, len(clique)+1):
            for sub in combinations(clique, k):
                f_vec[k] += 1

    # Actually, f_vec[k] = number of k-cliques. But we're overcounting.
    # Better: count maximal cliques by size, then derive f-vector.
    max_clique_sizes = Counter(len(c) for c in all_cliques)

    # f-vector: f_k = number of complete subgraphs of size k+1
    # f_0 = V, f_1 = E, f_2 = triangles, f_3 = K_4 subgraphs, etc.

    # Count directly
    f = [0] * 10
    f[0] = nm  # vertices
    f[1] = ne  # edges

    # Triangles
    tri = 0
    for i in range(nm):
        for j in range(i+1, nm):
            if not adj[i][j]: continue
            for k in range(j+1, nm):
                if adj[i][k] and adj[j][k]:
                    tri += 1
    f[2] = tri

    # K_4
    k4 = 0
    for i in range(nm):
        for j in range(i+1, nm):
            if not adj[i][j]: continue
            for k in range(j+1, nm):
                if not (adj[i][k] and adj[j][k]): continue
                for l in range(k+1, nm):
                    if adj[i][l] and adj[j][l] and adj[k][l]:
                        k4 += 1
    f[3] = k4

    # K_5 (only for n=5 which is small enough)
    k5 = 0
    if nm <= 50:
        for combo in combinations(range(nm), 5):
            if all(adj[combo[a]][combo[b]] for a in range(5) for b in range(a+1, 5)):
                k5 += 1
    f[4] = k5

    # Maximal clique sizes
    print("  n=%d: V_merged=%d, E_merged=%d" % (n, nm, ne))
    print("  Maximal clique sizes: %s" % dict(max_clique_sizes))
    print("  Number of maximal cliques: %d" % len(all_cliques))
    print("  Clique number (max clique size): %d" % max(len(c) for c in all_cliques))

    # f-vector
    fv = [f[k] for k in range(5) if f[k] > 0]
    print("  f-vector (faces): %s" % fv)
    print("    f_0 (vertices) = %d" % f[0])
    print("    f_1 (edges) = %d" % f[1])
    print("    f_2 (triangles) = %d" % f[2])
    print("    f_3 (K_4s) = %d" % f[3])
    print("    f_4 (K_5s) = %d" % f[4])

    # Euler characteristic
    chi = sum((-1)**k * f[k] for k in range(5))
    print("  Euler characteristic (truncated): %d" % chi)

    # Compare with cube f-vector
    # m-cube f_k = C(m, k) * 2^{m-k}
    print("\n  Compare with Q_%d f-vector:" % m)
    for k in range(min(5, m+1)):
        cube_fk = comb(m, k) * 2**(m-k)
        ratio = f[k] / cube_fk if cube_fk > 0 and f[k] > 0 else 0
        print("    f_%d: G_n/Z_2 = %d, Q_%d = %d, ratio = %.6f" %
              (k, f[k], m, cube_fk, ratio))

    # Independence number (max independent set)
    # This is related to the chromatic number
    # For now, just compute max degree
    degs = [sum(1 for j in range(nm) if adj[i][j]) for i in range(nm)]
    print("\n  Degree stats: min=%d, max=%d, avg=%.2f" %
          (min(degs), max(degs), sum(degs)/nm))

    # Is it a polytope skeleton? Key test: is it k-connected for the right k?
    # A d-polytope skeleton is d-connected (Balinski's theorem)
    # Check: what is the vertex connectivity?
    # Simple test: remove each vertex, check if graph stays connected

    def is_connected_without(excluded):
        """Check if graph minus excluded vertices is connected."""
        remaining = [v for v in range(nm) if v not in excluded]
        if len(remaining) == 0: return True
        visited = {remaining[0]}
        queue = [remaining[0]]
        while queue:
            v = queue.pop(0)
            for w in remaining:
                if w not in visited and adj[v][w]:
                    visited.add(w); queue.append(w)
        return len(visited) == len(remaining)

    # Vertex connectivity: min number of vertices to remove to disconnect
    # Test 1-connectivity
    min_cut = nm  # upper bound
    for v in range(nm):
        if not is_connected_without({v}):
            min_cut = 1; break

    if min_cut > 1:
        # Test 2-connectivity
        found_2cut = False
        for v in range(nm):
            for w in range(v+1, nm):
                if not is_connected_without({v, w}):
                    min_cut = 2; found_2cut = True; break
            if found_2cut: break

    if min_cut > 2:
        # Test 3-connectivity
        found_3cut = False
        for v in range(nm):
            for w in range(v+1, nm):
                for x in range(w+1, nm):
                    if not is_connected_without({v, w, x}):
                        min_cut = 3; found_3cut = True; break
                if found_3cut: break
            if found_3cut: break

    print("  Vertex connectivity >= %d" % min_cut)
    print("  (For d-polytope, need >= d by Balinski's theorem)")
    print()

# ============================================================
banner("PART 3: THE H-LATTICE STRUCTURE")
# ============================================================

# Is the H-ordering on G_n/Z_2 a LATTICE?
# The H-ordering: mi <= mj iff H(mi) <= H(mj)
# This is a PARTIAL order (multiple classes can have same H)

# For it to be a lattice, every pair must have a meet (greatest lower bound)
# and join (least upper bound) IN THE POSET.

# But the H-ordering is a TOTAL PREORDER (every pair comparable by H value).
# If H values are distinct for all classes, it's a total order = a chain = trivially a lattice.
# If H values repeat, it's not a total order.

# Actually, the H-ordering is NOT a partial order on G_n/Z_2 in the usual sense.
# The DAG structure uses H as a WEIGHT, not as the ordering relation.
# The actual ordering we should look at is the REACHABILITY ordering:
# mi <= mj iff there is a DIRECTED path from mi to mj in the H-increasing DAG.

# Let me look at the H-DAG more carefully.
# Actually, from S239, the H-DAG is: orient each edge of G_n/Z_2 from
# lower H to higher H. Level edges (same H) are undirected.

# For the merged meta-graph, let me check if the H-DAG has lattice structure.

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
cmap = {}; sizes = {}; reps = {}; tc = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap:
        cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        reps[cid] = [row[:] for row in A]
    sizes[cmap[c]] += 1; tc[bits] = cmap[c]
nc = len(cmap)

info5 = {}
for cid in range(nc):
    A = reps[cid]
    comp_c = canon(complement(A, n), n)
    info5[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c]}

merged5 = {}; mid5 = {}; mid = 0
for cid in range(nc):
    if cid in mid5: continue
    comp = info5[cid]['comp']
    merged5[mid] = (cid,) if comp == cid else (cid, comp)
    mid5[cid] = mid
    if comp != cid: mid5[comp] = mid
    mid += 1
nm5 = mid

F5 = np.zeros((nc, nc), dtype=int)
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m): F5[ci][tc[bits^(1<<k)]] += 1

adj5 = [[False]*nm5 for _ in range(nm5)]
for ci in range(nc):
    mi = mid5[ci]
    for cj in range(nc):
        if F5[ci][cj] > 0:
            mj = mid5[cj]
            if mi != mj:
                adj5[mi][mj] = True; adj5[mj][mi] = True

# H values for merged nodes
H5 = {mi: info5[merged5[mi][0]]['H'] for mi in range(nm5)}

# Build H-DAG: directed edges from lower H to higher H
# For level edges (same H), we need to handle specially
print("  H-DAG at n=5:")
print("  Merged nodes by H:")
for mi in sorted(range(nm5), key=lambda m: H5[m]):
    nbrs = [mj for mj in range(nm5) if adj5[mi][mj]]
    up = [mj for mj in nbrs if H5[mj] > H5[mi]]
    down = [mj for mj in nbrs if H5[mj] < H5[mi]]
    level = [mj for mj in nbrs if H5[mj] == H5[mi]]
    print("    mid=%d H=%d: up=%s down=%s level=%s" %
          (mi, H5[mi], up, down, level))

# Check if the H-directed graph is a lattice (Hasse diagram of a lattice)
# For this, every pair of nodes must have a unique meet and join.
# With the H-gradient: the reachability order is a partial order.

# Compute reachability
reach = [[False]*nm5 for _ in range(nm5)]
for mi in range(nm5): reach[mi][mi] = True

# Topological order by H
order = sorted(range(nm5), key=lambda m: H5[m])
for mi in order:
    for mj in range(nm5):
        if adj5[mi][mj] and H5[mj] > H5[mi]:
            # mi -> mj directed edge
            for mk in range(nm5):
                if reach[mk][mi]: reach[mk][mj] = True

print("\n  Reachability matrix (can node i reach node j via H-increasing path):")
for mi in order:
    reachable = [mj for mj in order if reach[mi][mj] and mi != mj]
    print("    mid=%d H=%d reaches: %s" % (mi, H5[mi], reachable))

# Check lattice property: for every pair (a,b), is there a unique join and meet?
is_lattice = True
for a in range(nm5):
    for b in range(a+1, nm5):
        # Find join = least upper bound
        upper_bounds = [c for c in range(nm5) if reach[a][c] and reach[b][c]]
        if not upper_bounds:
            # No upper bound at all — fails for pairs not in same component
            is_lattice = False
            print("  FAIL: no upper bound for (mid=%d, mid=%d)" % (a, b))
            break
        # Find minimal upper bounds
        min_ub = [c for c in upper_bounds if not any(reach[c][d] and d != c for d in upper_bounds)]
        if len(min_ub) != 1:
            is_lattice = False
            print("  FAIL: join not unique for (mid=%d H=%d, mid=%d H=%d): min_ub=%s" %
                  (a, H5[a], b, H5[b], [(m, H5[m]) for m in min_ub]))
            break

        # Find meet = greatest lower bound
        lower_bounds = [c for c in range(nm5) if reach[c][a] and reach[c][b]]
        if not lower_bounds:
            is_lattice = False
            print("  FAIL: no lower bound for (mid=%d, mid=%d)" % (a, b))
            break
        max_lb = [c for c in lower_bounds if not any(reach[d][c] and d != c for d in lower_bounds)]
        if len(max_lb) != 1:
            is_lattice = False
            print("  FAIL: meet not unique for (mid=%d H=%d, mid=%d H=%d): max_lb=%s" %
                  (a, H5[a], b, H5[b], [(m, H5[m]) for m in max_lb]))
            break
    if not is_lattice: break

print("\n  H-DAG at n=5 is a lattice: %s" % is_lattice)

# ============================================================
banner("PART 4: QUOTIENTOPE SPECTRUM POSITION")
# ============================================================

print("""
  THE QUOTIENTOPE SPECTRUM (Pilaud-Santos):

  PERMUTOHEDRON Π_n ← (weak order on S_n, n! vertices)
    |
    | lattice congruences contract edges
    |
  ASSOCIAHEDRON ← (sylvester congruence, Catalan number vertices)
    |
    | more contractions
    |
  m-CUBE Q_m ← (boolean congruence, 2^m vertices, m = n-1 for cubes)
    |
    | identity on cube
    |
  Q_m

  OUR CONSTRUCTION:
  Q_{C(n,2)} (the tournament cube)
    |
    | quotient by S_n action on coordinates
    |
  G_n (iso class graph, A000568 vertices)
    |
    | quotient by complement (Z_2 bit-flip)
    |
  G_n / Z_2 (merged meta-graph)

  G_n is NOT a quotientope in the Pilaud-Santos sense because:
  1. It's a quotient of Q_m, not of the permutohedron Π_n
  2. The S_n action on pairs is NOT a lattice congruence
  3. G_n has fewer vertices than ANY quotientope of S_n (which has ≥ n+1)

  HOWEVER, the CUBE Q_m is a polytope, and G_n is its orbit graph.
  The orbit polytope (convex hull of orbit representatives) IS a polytope,
  and G_n might be its 1-skeleton.

  THE ANALOGY:
  Quotientope         | Our construction
  --------------------|------------------
  Weak order on S_n   | Q_m = tournament cube
  Lattice congruence  | S_n × Z_2 group action
  Contraction         | Orbit identification
  Quotient lattice    | Orbit graph G_n/Z_2
  Polytope            | ??? (is G_n/Z_2 polytopal?)

  THE KEY QUESTION:
  Is G_n/Z_2 the 1-skeleton of a polytope?
  If so, what polytope? What are its higher faces?
""")

# ============================================================
banner("PART 5: SUBCUBE STRUCTURE")
# ============================================================

# In Q_m, a k-face is a k-subcube: fix m-k coordinates, vary k coordinates.
# Under the S_n quotient, a k-face of G_n corresponds to an orbit of k-subcubes.

# The number of k-faces of G_n should be:
# f_k(G_n) = #{orbits of k-subcubes of Q_m under S_n}

# A k-subcube is determined by:
# 1. A subset S of k coordinates (arcs to flip)
# 2. A partial orientation of the remaining m-k arcs

# Under S_n: orbits of (S, partial orientation) pairs.
# The orbit of S alone is determined by the "type" of the k-subset of arcs.

# For k=1: 1-subcubes = edges. f_1 = E(n) = our edge sequence.
# For k=0: 0-subcubes = vertices. f_0 = V(n) = A000568.

# For k=2: 2-subcubes = squares (4-cycles in Q_m).
# A 2-subcube is determined by choosing 2 arcs to flip.
# The 4 vertices of the square are: (T, T_e, T_f, T_{ef}).
# In G_n: the square maps to a 4-cycle or collapses.

# The number of squares in G_n is related to the "diamond" count:
# pairs of edges sharing both endpoints but via different arcs.

# Let me count squares at n=5

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
cmap = {}; tc = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k,(ii,jj) in enumerate(pairs):
        if (bits>>k)&1: A[ii][jj]=1
        else: A[jj][ii]=1
    c = canon(A, n)
    if c not in cmap: cmap[c] = len(cmap)
    tc[bits] = cmap[c]
nc = len(cmap)

# Count distinct 4-cycles (squares) in G_n at n=5
# A square in Q_m: choose 2 arc positions k1, k2.
# Vertices: bits, bits^k1, bits^k2, bits^k1^k2.
# In G_n: map each to its class. The square is (c0, c1, c2, c3).

# Count distinct square orbits
square_types = set()
for bits in range(2**m):
    c0 = tc[bits]
    for k1 in range(m):
        c1 = tc[bits ^ (1<<k1)]
        for k2 in range(k1+1, m):
            c2 = tc[bits ^ (1<<k2)]
            c3 = tc[bits ^ (1<<k1) ^ (1<<k2)]
            # The square is the 4-tuple (c0, c1, c2, c3)
            # Normalize: order doesn't matter for the FACE, but the
            # set of 4 classes does
            face = frozenset({c0, c1, c2, c3})
            if len(face) == 4:  # a proper square (all 4 vertices distinct)
                square_types.add(face)

print("  n=5: distinct 4-element faces (squares) in G_n: %d" % len(square_types))

# How many are proper 4-cycles (C_4)?
# A proper square has edges: c0-c1, c0-c2, c1-c3, c2-c3
# (from the cube structure: flipping k1 or k2)
# Check: are all 4 edges present in G_n?
F5b = np.zeros((nc, nc), dtype=int)
for bits in range(2**m):
    ci = tc[bits]
    for k in range(m):
        F5b[ci][tc[bits^(1<<k)]] += 1
adj5b = (F5b > 0).astype(int)
np.fill_diagonal(adj5b, 0)

proper_squares = 0
degenerate = 0
for face in square_types:
    verts = list(face)
    # Check all 6 possible edges among 4 vertices
    edges_present = sum(1 for i in range(4) for j in range(i+1, 4) if adj5b[verts[i]][verts[j]])
    if edges_present >= 4:
        proper_squares += 1
    else:
        degenerate += 1

print("  Proper squares (≥4 edges among 4 vertices): %d" % proper_squares)
print("  Degenerate (some edges missing): %d" % degenerate)

# Also: how many of the 4-vertex sets form K_4 (complete = tetrahedra)?
k4_count = sum(1 for face in square_types
               if all(adj5b[v1][v2] for v1 in face for v2 in face if v1 != v2))
print("  Of which are K_4 (complete): %d" % k4_count)

# ============================================================
banner("PART 6: THE TOURNAMENT POLYTOPE")
# ============================================================

print("""
  THE TOURNAMENT POLYTOPE P_n:
  Convex hull of all tournament adjacency vectors in R^m.

  Each tournament T gives a 0/1 vector x_T ∈ {0,1}^m:
    x_T[k] = A[i][j] where (i,j) = pairs[k]

  P_n = conv(x_T : T tournament on n vertices)

  Since x_T + x_{T^c} = (1,1,...,1) (complement flips all bits),
  P_n is centrally symmetric about (1/2, ..., 1/2).

  P_n is a SUBSET of the m-cube [0,1]^m.

  The ORBIT POLYTOPE: conv(one representative per iso class)
  This is the tournament polytope modulo S_n.

  KEY FACT: The orbit polytope has V(n) = A000568 vertices
  (one per iso class). Its 1-skeleton should relate to G_n.

  DIMENSION of P_n:
  For n≥3, P_n is full-dimensional (dim m = C(n,2)):
  The m edge orientations are "almost" independent (the only constraint
  is that they form a tournament, but since ALL 2^m orientations are
  tournaments, the polytope is the entire m-cube).

  Wait — EVERY 0/1 vector in {0,1}^m IS a tournament adjacency vector!
  (Given any orientation of all C(n,2) pairs, we get a tournament.)

  So P_n = the m-cube! Its orbit polytope is Q_m / S_n.

  THE ORBIT POLYTOPE has dimension:
  dim = m - dim(fixed subspace of S_n action)
  The S_n action on R^m permutes coordinates.
  Fixed subspace: vectors invariant under all permutations.
  The only invariant vector is the all-ones vector (and scalar multiples).
  Actually, the orbits of {pairs} under S_n have sizes:
  there are just a few types of pairs... no, all pairs are in the same S_n orbit.
  So the fixed subspace is 1-dimensional (the all-ones direction).
  Orbit polytope dimension = m - 1.

  Hmm, that's not quite right. The orbit polytope lives in R^m / S_n
  which is not a subspace of R^m. Let me think about this differently.

  The NUMBER OF ORBITS of S_n on the m coordinates is:
  #{orbits of pairs {i,j} under S_n} = 1 (since S_n acts transitively on pairs).
  So the orbit of any coordinate is the entire set of m coordinates.
  This means the "type" of a tournament is determined by:
  the multiset of values on each orbit.
  Since there's only 1 orbit, the type is just the number of 1s = sum of out-degrees / 2.
  Wait, that gives only the score sum, which is always m = C(n,2).

  No — the S_n action permutes pairs, but the ORBIT POLYTOPE uses the
  orbit-counting coordinates. The coordinates of the orbit polytope
  are the sums over each orbit.

  Since all pairs are in one orbit, the orbit polytope has just 1 coordinate:
  the total number of 1s = sum of out-degrees = m (constant for all tournaments!).
  This collapses everything to a single point — clearly wrong.

  The issue: I'm confusing the orbit polytope with something else.

  The correct object: G_n is NOT a polytope quotient in the standard sense.
  It's the orbit GRAPH (adjacency pattern among orbits), not the orbit polytope.
""")

print("Done. opus-2026-03-23-S240")
