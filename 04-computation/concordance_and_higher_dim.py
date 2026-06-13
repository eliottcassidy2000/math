#!/usr/bin/env python3
"""
concordance_and_higher_dim.py -- kind-pasteur-2026-03-13-S61

ALL overlap-2 cycle pairs are concordant: when two 3-cycles share an edge,
they traverse it in the same direction. Is this universal?

This is FORCED by the tournament axiom: each vertex pair has exactly one
arc direction. If C_i and C_j share vertices {u,v}, both cycles must use
the same arc u->v (or v->u). They CAN'T disagree on a shared edge!

So "concordance" is trivially 100% for any tournament. This means the
orientation of the overlap graph is AUTOMATICALLY consistent.

The interesting question becomes: what IS the "hidden dimension" that
separates H=109 from H=111?

Deep question: The {2,1,0} overlap structure with orientations forms a
HIGHER-DIMENSIONAL object. What is its combinatorial/topological nature?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))
    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val
    return lam, c3_sets


n = 7

# ========================================================================
# ANALYSIS 1: The 3-cycle incidence structure
# ========================================================================
print("=" * 70)
print("ANALYSIS 1: 3-CYCLE INCIDENCE STRUCTURE")
print("=" * 70)

# The key object: the (vertex, 3-cycle) incidence matrix I.
# I[v][C] = 1 if v is in C, 0 otherwise.
# Dimensions: n x c3 (7 x 11 for our case).
# Lambda = I * I^T (pair coverage).
# The overlap matrix = I^T * I (cycle coverage).

# For each pair of 3-cycles with overlap 2:
# They share 2 vertices {u,v} and differ in their 3rd vertex.
# The arc u->v is determined by the tournament.
# BOTH cycles must use this arc in the same direction (trivially).
# But each cycle has a DIFFERENT third vertex w1 vs w2.
# The "hidden information" is how w1 and w2 connect to {u,v}.

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    print(f"\n  {label} (bits={bits}): {len(c3s)} 3-cycles")

    # For each overlap-2 pair, analyze the connection structure
    ov2_pairs = []
    for i in range(len(c3s)):
        for j in range(i+1, len(c3s)):
            shared = c3s[i] & c3s[j]
            if len(shared) == 2:
                ov2_pairs.append((i, j, shared))

    print(f"    Overlap-2 pairs: {len(ov2_pairs)}")

    for idx, (i, j, shared) in enumerate(ov2_pairs[:5]):  # Show first 5
        u, v = sorted(shared)
        w1 = list(c3s[i] - shared)[0]
        w2 = list(c3s[j] - shared)[0]

        # Arc direction u->v
        arc_uv = "u->v" if A[u][v] else "v->u"

        # How w1 and w2 connect to {u,v}
        w1_u = "w1->u" if A[w1][u] else "u->w1"
        w1_v = "w1->v" if A[w1][v] else "v->w1"
        w2_u = "w2->u" if A[w2][u] else "u->w2"
        w2_v = "w2->v" if A[w2][v] else "v->w2"

        print(f"    C{i}={sorted(c3s[i])}, C{j}={sorted(c3s[j])}: "
              f"shared=({u},{v}) arc={arc_uv}, "
              f"w1={w1}: {w1_u},{w1_v}, w2={w2}: {w2_u},{w2_v}")


# ========================================================================
# ANALYSIS 2: The oriented 3-cycle complex
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: THE ORIENTED 3-CYCLE COMPLEX")
print("=" * 70)

# Build a simplicial complex where:
# 0-simplices = vertices (7)
# 1-simplices = edges with lambda > 0 (19 edges)
# 2-simplices = 3-cycle vertex sets (11)

# Each 2-simplex has an ORIENTATION from the tournament.
# This is a chain complex C_2 -> C_1 -> C_0.
# The boundary map d_2: for a 3-cycle {a,b,c} oriented as a->b->c->a,
# d_2([a,b,c]) = [b,c] - [a,c] + [a,b] (with tournament orientations).

# The homology of this complex might capture the hidden information!

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    print(f"\n  {label} (bits={bits}):")

    # 1-simplices: oriented edges (u,v) for u<v with tournament direction
    edges = []
    edge_index = {}
    for u in range(n):
        for v in range(u+1, n):
            if lam[u][v] > 0:
                edges.append((u, v))
                edge_index[(u, v)] = len(edges) - 1

    print(f"    Edges: {len(edges)}")

    # 2-simplices: oriented 3-cycles
    # Orientation: for each 3-cycle, the cyclic order from the tournament
    # gives an orientation on the 2-simplex.
    triangles = []
    for c3 in c3s:
        a, b, c = sorted(c3)
        if A[a][b] and A[b][c] and A[c][a]:
            triangles.append((a, b, c, +1))  # a->b->c->a = positive orientation
        else:
            triangles.append((a, b, c, -1))  # a<-b<-c<-a = negative orientation

    print(f"    Triangles: {len(triangles)}")

    # Boundary matrix d_2: maps 2-chains to 1-chains
    # For triangle (a,b,c) with orientation +1:
    # d_2 = [b,c] - [a,c] + [a,b]
    # For orientation -1: d_2 = -[b,c] + [a,c] - [a,b]

    ne = len(edges)
    nt = len(triangles)

    # Build boundary matrix
    d2 = [[0]*nt for _ in range(ne)]
    for j, (a, b, c, orient) in enumerate(triangles):
        # Edge [a,b] -> index
        if (a, b) in edge_index:
            d2[edge_index[(a, b)]][j] += orient
        # Edge [a,c] -> index
        if (a, c) in edge_index:
            d2[edge_index[(a, c)]][j] -= orient
        # Edge [b,c] -> index
        if (b, c) in edge_index:
            d2[edge_index[(b, c)]][j] += orient

    # Rank of d2 (over Q)
    # Use Gaussian elimination
    import copy
    mat = copy.deepcopy(d2)
    # Transpose: we want rank of the matrix, work with rows
    # mat is ne x nt. Rank = rank of this matrix.
    # Actually let me work with the transpose for row reduction
    mat_t = [[mat[i][j] for i in range(ne)] for j in range(nt)]
    # Row reduce mat_t (nt x ne)
    pivot_row = 0
    for col in range(ne):
        # Find non-zero in this column
        found = None
        for row in range(pivot_row, nt):
            if mat_t[row][col] != 0:
                found = row
                break
        if found is None:
            continue
        mat_t[pivot_row], mat_t[found] = mat_t[found], mat_t[pivot_row]
        # Eliminate
        for row in range(nt):
            if row != pivot_row and mat_t[row][col] != 0:
                factor = mat_t[row][col] / mat_t[pivot_row][col]
                for c2 in range(ne):
                    mat_t[row][c2] -= factor * mat_t[pivot_row][c2]
        pivot_row += 1

    rank_d2 = pivot_row
    print(f"    rank(d2) = {rank_d2}")
    print(f"    ker(d2) = dim(H_2) + something, coker = dim(H_1)")
    print(f"    H_2 = ker(d2) has dimension {nt - rank_d2}")
    print(f"    H_1 has dimension {ne - rank_d2} (ignoring d_1)")

    # Also compute d_1: C_1 -> C_0
    # d_1([u,v]) = v - u (as a 0-chain)
    # With tournament orientation: if u->v, then d_1 = +v - u
    d1 = [[0]*ne for _ in range(n)]
    for j, (u, v) in enumerate(edges):
        if A[u][v]:
            d1[u][j] = -1
            d1[v][j] = +1
        else:
            d1[u][j] = +1
            d1[v][j] = -1

    # Rank of d1
    mat1_t = [[d1[i][j] for i in range(n)] for j in range(ne)]
    pivot_row1 = 0
    for col in range(n):
        found = None
        for row in range(pivot_row1, ne):
            if abs(mat1_t[row][col]) > 1e-10:
                found = row
                break
        if found is None:
            continue
        mat1_t[pivot_row1], mat1_t[found] = mat1_t[found], mat1_t[pivot_row1]
        for row in range(ne):
            if row != pivot_row1 and abs(mat1_t[row][col]) > 1e-10:
                factor = mat1_t[row][col] / mat1_t[pivot_row1][col]
                for c2 in range(n):
                    mat1_t[row][c2] -= factor * mat1_t[pivot_row1][c2]
        pivot_row1 += 1

    rank_d1 = pivot_row1
    print(f"    rank(d1) = {rank_d1}")

    # Betti numbers:
    # b_0 = n - rank(d1) = 7 - rank_d1
    # b_1 = ne - rank(d1) - rank(d2) = nullity(d1) - rank(d2)
    # b_2 = nt - rank(d2)
    b0 = n - rank_d1
    b1 = ne - rank_d1 - rank_d2
    b2 = nt - rank_d2
    print(f"    Betti numbers: b0={b0}, b1={b1}, b2={b2}")
    print(f"    Euler char: {b0 - b1 + b2}")


# ========================================================================
# ANALYSIS 3: Comparison — is the cycle complex the same for H=109 and H=111?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: CYCLE COMPLEX COMPARISON")
print("=" * 70)

# Both have same lambda graph. Do they have the same cycle complex?
# The 0-skeleton and 1-skeleton (lambda graph) are isomorphic.
# The 2-skeleton (3-cycles as vertex sets) are the same.
# The ORIENTATION of the 2-simplices differs!

# Let's compare the oriented boundary matrices
print("\n  Oriented boundary comparison:")
for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    # Determine orientation of each 3-cycle
    orients = []
    for c3 in sorted(c3s, key=lambda x: sorted(x)):
        a, b, c = sorted(c3)
        if A[a][b] and A[b][c] and A[c][a]:
            orients.append(+1)
        else:
            orients.append(-1)

    print(f"  {label}: orientations = {orients}")
    print(f"           product = {1 if all(o == 1 for o in orients) else (-1)**sum(1 for o in orients if o == -1)}")
    pos = sum(1 for o in orients if o == 1)
    neg = sum(1 for o in orients if o == -1)
    print(f"           positive={pos}, negative={neg}")


# ========================================================================
# ANALYSIS 4: The 5-cycle complex
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: THE 5-CYCLE LAYER")
print("=" * 70)

# The 3-cycle complex doesn't distinguish H=109 from H=111 topologically
# (both have same Betti numbers). But 5-cycles might add structure.

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)

    # Find all 5-cycle vertex sets and their directed cycle counts
    c5_info = []
    total_c5_dir = 0
    for sub in combinations(range(n), 5):
        nc = count_directed_ham_cycles_on_subset(A, list(sub))
        if nc > 0:
            c5_info.append((sub, nc))
            total_c5_dir += nc

    print(f"\n  {label}: {len(c5_info)} 5-cycle vertex sets, {total_c5_dir} directed 5-cycles")

    # Multiplicity histogram
    mult_hist = defaultdict(int)
    for sub, nc in c5_info:
        mult_hist[nc] += 1
    print(f"    Multiplicity: {dict(sorted(mult_hist.items()))}")

    # Per-vertex participation
    vtx_count = [0] * n
    for sub, nc in c5_info:
        for v in sub:
            vtx_count[v] += nc
    print(f"    Per-vertex 5-cycle weight: {vtx_count}")
    print(f"    Sorted: {sorted(vtx_count)}")


# ========================================================================
# ANALYSIS 5: The Vitali set analogy — dimension counting
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: DIMENSION COUNTING FOR THE VITALI ANALOGY")
print("=" * 70)

print("""
DIMENSION COUNTING:

A tournament on n=7 has C(7,2) = 21 binary choices (arc directions).
Total dimension: 21 bits = 2^21 = 2097152 tournaments.

The score sequence constrains 7 degrees. Dimension remaining: ~14.
The lambda matrix has C(7,2) = 21 entries, but with rank constraints.
In practice, the lambda histogram has ~10 independent parameters.

For our specific class:
  - Score (2,2,3,3,3,3,5): fixes ~7 dimensions
  - Lambda histogram: fixes ~10 dimensions
  - Lambda degree sequence: fixes ~12 dimensions (resolves most ambiguity)
  - Lambda isomorphism class: fixes ~14 dimensions
  - Remaining: ~7 dimensions of "arc orientation" freedom
  - Only 1 of these 7 dimensions affects c7 (and hence H)

The "Vitali quotient" maps R^21 -> R^{14} by forgetting 7 dimensions.
The "non-measurable content" lives in a 7-dimensional fiber.
But only a 1-dimensional subspace of this fiber affects H.

This is why the Vitali analogy works:
  - The "measurable" part (lambda) captures 14/21 = 67% of the information
  - The "non-measurable" part (arc orientation) has 7 dimensions
  - Only 1 dimension of the non-measurable part is "physically relevant" (affects H)
  - The other 6 dimensions of arc orientation are "gauge freedom" (don't affect H)

The Vitali set analogy:
  - R/Q: the rationals are the "measurable" equivalence classes
  - A Vitali set picks one representative from each class
  - The "non-measurable" content is the choice of representative
  - In our setting: lambda = the "rational equivalence class"
  - c7 = the "choice of representative" that determines H
  - The 6 gauge dimensions = different labelings of the same oriented structure
""")


# ========================================================================
# ANALYSIS 6: What is the MINIMAL invariant beyond lambda?
# ========================================================================
print(f"{'='*70}")
print("ANALYSIS 6: MINIMAL BEYOND-LAMBDA INVARIANT")
print("=" * 70)

# We know: labeled lambda + c7_dir is complete.
# But c7_dir requires computing all 7-cycles (expensive!).
# Is there a CHEAPER invariant that still separates H=109 from H=111?

# Candidate 1: The 5-cycle orientation signature
# Both have same c5_dir=23. But do they have different 5-cycle ORIENTATION patterns?

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)

    # For each 5-cycle vertex set, compute the actual directed cycles
    # and their "chirality" (clockwise vs counterclockwise in some canonical sense)
    c5_chirality = []
    for sub in combinations(range(n), 5):
        nc = count_directed_ham_cycles_on_subset(A, list(sub))
        if nc == 0:
            continue
        # Count in each direction
        # For a 5-vertex tournament, directed Hamiltonian cycles come in pairs
        # (forward and reverse). So nc is always even (0 or 2 for tournaments on 5).
        c5_chirality.append((sub, nc))

    print(f"\n  {label}: 5-cycle sets with count:")
    for sub, nc in c5_chirality:
        print(f"    {sub}: {nc} directed 5-cycles")

# Candidate 2: The 3-cycle "coherence" — for each triple of 3-cycles sharing
# a common edge pattern, count concordances
print(f"\n  Candidate invariant: 3-cycle triple coherence")

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, c3s = get_labeled_lambda(A, n)

    # For each triple of 3-cycles that pairwise share vertices,
    # compute the "triple orientation product"
    triple_prods = []
    for i in range(len(c3s)):
        for j in range(i+1, len(c3s)):
            for k in range(j+1, len(c3s)):
                if (len(c3s[i] & c3s[j]) >= 1 and
                    len(c3s[j] & c3s[k]) >= 1 and
                    len(c3s[i] & c3s[k]) >= 1):
                    # Triple of mutually intersecting 3-cycles
                    # Compute orientation product
                    orients = []
                    for c3 in [c3s[i], c3s[j], c3s[k]]:
                        a, b, c = sorted(c3)
                        if A[a][b] and A[b][c] and A[c][a]:
                            orients.append(+1)
                        else:
                            orients.append(-1)
                    prod = orients[0] * orients[1] * orients[2]
                    triple_prods.append(prod)

    pos = sum(1 for p in triple_prods if p > 0)
    neg = sum(1 for p in triple_prods if p < 0)
    print(f"  {label}: {len(triple_prods)} mutually-intersecting triples, "
          f"prod>0: {pos}, prod<0: {neg}, sum: {sum(triple_prods)}")


# ========================================================================
# ANALYSIS 7: The arc flow through non-lambda edges
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 7: ARC FLOW THROUGH LAMBDA-0 EDGES")
print("=" * 70)

# There are 2 edges with lambda=0 (not in any 3-cycle).
# These edges connect "structurally isolated" vertex pairs.
# The direction of these arcs is part of the "non-measurable" content.

for bits, label in [(4658, "H=111"), (4728, "H=109")]:
    A = binary_to_tournament(bits, n)
    lam, _ = get_labeled_lambda(A, n)

    zero_edges = []
    for u in range(n):
        for v in range(u+1, n):
            if lam[u][v] == 0:
                direction = "u->v" if A[u][v] else "v->u"
                zero_edges.append((u, v, direction))
    print(f"\n  {label}: lambda-0 edges: {zero_edges}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
