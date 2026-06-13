#!/usr/bin/env python3
"""
ocf_homology_bridge.py — Can we connect Ω(T) structure to path homology?

TWO PARALLEL STORIES in this project:
  A. H(T) = I(Ω(T), 2) — the OCF (proved by Grinberg-Stanley)
  B. β_2(T) = 0 for all tournaments — the universal vanishing conjecture

Both concern CYCLES in tournaments. Can we connect them?

KEY OBSERVATION:
  Ω(T) is the conflict graph of odd directed cycles.
  β_2(T) measures whether there are "holes" in the path complex.
  β_2 = 0 means: every 2-cycle in the path complex bounds a 3-chain.

For circulant tournaments:
  Ω(T) has circulant symmetry (vertex-transitive)
  Path complex has circulant symmetry too
  β_k(T) is computable from the eigenspace decomposition (THM-125)

QUESTION: Is there a FORMULA relating I(Ω(T), 2) and β_*(T)?

From the data:
  T_7 (Paley): H = 189, β = [1,0,0,0,6,0,0]
  T_7 (Interval): H = 175, β = ???
  T_11 (Paley): H = 95095, β = [1,0,0,0,0,15,0,0,0,0,0]
  T_11 (Interval): H = 93027, β = ???

The Paley tournament T_p has β_{(p-1)/2 - 1} = p-2 (from THM-125).
What about the Interval tournament?

Author: opus-2026-03-12-S65
"""

import numpy as np
from itertools import combinations

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def tournament_adjacency(sigma, p):
    m = (p-1)//2
    n = p
    A = np.zeros((n, n), dtype=int)
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def compute_path_homology_small(A, n, max_deg=4):
    """Compute path homology of a tournament for small n.
    Returns dict of Betti numbers."""

    # Allowed p-paths: sequences (v_0,...,v_p) with A[v_i][v_{i+1}]=1
    # and all vertices distinct (elementary paths)

    betti = {}

    for deg in range(max_deg + 1):
        # Basis for Omega_deg: "allowed" paths that are in ker(face maps sum)
        # For deg=0: all vertices, dim = n, rank(d_0) = 0
        # For deg=1: all edges i→j with A[i][j]=1
        # For deg=2: paths i→j→k with all edges present and "allowed"

        if deg == 0:
            # Omega_0 = R^n
            # d_0 = 0 (no boundary of 0-paths)
            # beta_0 = dim(ker d_0 restricted to Omega_0) = n - rank(im d_1 in Omega_0)
            # Actually beta_0 = dim Omega_0 / im d_1
            pass

        # For simplicity, compute beta_1 only for small n
        if deg == 1 and n <= 12:
            # 1-paths = edges
            edges = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]]

            # d_1: maps edge (i,j) to j - i
            # Omega_1 = {paths p : d_0(p) = 0 in Omega_0}
            # Actually Omega_1 = regular 1-paths (all edges)
            # d_1(i→j) = j - i (as formal sum of vertices)

            # Boundary matrix d_1 (edges → vertices)
            d1 = np.zeros((n, len(edges)), dtype=int)
            for idx, (i, j) in enumerate(edges):
                d1[j][idx] += 1
                d1[i][idx] -= 1

            # rank(d_1)
            rank_d1 = np.linalg.matrix_rank(d1)

            # For beta_1, need rank(d_2) where d_2 maps 2-paths to edges

            # 2-paths = allowed elementary paths of length 2
            paths_2 = []
            for i in range(n):
                for j in range(n):
                    if i == j or not A[i][j]: continue
                    for k in range(n):
                        if k == i or k == j or not A[j][k]: continue
                        # Check "allowed": the path (i,j,k) is in Omega_2
                        # iff d_1(i→j→k) ∈ Omega_1
                        # d_1(i→j→k) = (j→k) - (i→k) + (i→j)
                        # This is in Omega_1 iff all three terms are edges
                        # i.e., A[i][k] = 1 or A[k][i] = 1
                        # Actually the "allowed" condition for GLMY:
                        # A 2-path v0→v1→v2 is allowed if EVERY face is an allowed 1-path
                        # faces: (v1,v2), (v0,v2), (v0,v1)
                        # "allowed" 1-path = directed edge. So need A[v0][v2] = 1 OR A[v2][v0]=1
                        # For tournaments: ALWAYS true (one of the two exists)
                        paths_2.append((i, j, k))

            if len(paths_2) > 0:
                # Boundary d_2: maps 2-path (i,j,k) to edges
                edge_idx = {e: idx for idx, e in enumerate(edges)}
                d2 = np.zeros((len(edges), len(paths_2)), dtype=int)
                for pidx, (i, j, k) in enumerate(paths_2):
                    # d_2(i→j→k) = (j→k) - (i→k) + (i→j)
                    # But (i→k) might not be an edge!
                    # In GLMY, the boundary includes ALL faces
                    if (j, k) in edge_idx:
                        d2[edge_idx[(j, k)]][pidx] += 1
                    if (i, k) in edge_idx:
                        d2[edge_idx[(i, k)]][pidx] -= 1
                    elif (k, i) in edge_idx:
                        # If i→k not edge, then k→i is. But this is a DIFFERENT edge.
                        # In GLMY, the face is still (i,k) — which might not be allowed.
                        # Actually for tournaments, exactly one of i→k or k→i exists.
                        # The boundary d_2(i→j→k) = (j,k) - (i,k) + (i,j)
                        # where (i,k) is the UNIQUE directed edge between i and k.
                        # If A[k][i]=1, then the edge is k→i, and the face is (k,i)
                        d2[edge_idx[(k, i)]][pidx] += 1  # sign depends on convention
                    if (i, j) in edge_idx:
                        d2[edge_idx[(i, j)]][pidx] += 1

                rank_d2 = np.linalg.matrix_rank(d2)
            else:
                rank_d2 = 0

            # beta_1 = dim(ker d_1 restricted to Omega_1) - rank(d_2)
            # dim(ker d_1) = len(edges) - rank_d1
            # But we need Omega_1 subspace... for GLMY, Omega_1 = ker(partial_0 restricted to A_1)
            # For tournaments, Omega_1 = all edges (since all faces of 1-paths are vertices)

            beta_1_approx = len(edges) - rank_d1 - rank_d2
            # This isn't quite right for GLMY but gives an approximation
            print(f"    n={n}: |edges|={len(edges)}, rank(d1)={rank_d1}, "
                  f"|2-paths|={len(paths_2)}, rank(d2)={rank_d2}, "
                  f"beta_1 ≈ {beta_1_approx}")

    return betti

print("=" * 70)
print("CONNECTING Ω(T) TO PATH HOMOLOGY")
print("=" * 70)

# Instead of computing path homology (complex), let's focus on a
# cleaner connection: the relationship between α_k(Ω) and cycle counts

print("""
STRUCTURAL RELATIONSHIP:

The Ω graph encodes which odd cycles can coexist in an independent set.
The independence polynomial I(Ω, λ) at λ=2 gives H.

For PATH HOMOLOGY: the boundary maps involve PATHS, not cycles.
β_k measures obstructions to extending (k-1)-dimensional path homology.

KEY QUESTION: Do the independent sets of Ω(T) relate to the
homology generators?

At n=7 (Paley, p=7):
  β_4 = 6 = p-1 (from eigenspace decomposition)
  H = 189 = 27 × 7 = 3³ × 7
  α_1(Ω) = 80 (all odd cycles), α_2 = 7
  I(Ω, 2) = 1 + 160 + 28 = 189 ✓

The 6 generators of H_4(T_7) correspond to 6 linearly independent
4-dimensional "holes" in the path complex.

These are NOT the same as the 7 independent pairs in α_2(Ω).
But they live in complementary spaces:
  - α_2 captures cycle DISJOINTNESS (vertex-based)
  - β_4 captures path INDEPENDENCE (chain-based)

Can the β_4 generators be expressed in terms of Ω structure?
""")

# At p=7, compute some relevant quantities
p = 7
m = 3

# Interval tournament
sigma_int = tuple(1 for _ in range(m))
A_int = tournament_adjacency(sigma_int, p)

# Paley tournament
sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))
A_pal = tournament_adjacency(sigma_pal, p)

for name, A in [("Interval", A_int), ("Paley", A_pal)]:
    print(f"\n--- {name} tournament at p={p} ---")

    # Score sequence
    scores = sorted([sum(A[i]) for i in range(p)])
    print(f"  Score sequence: {scores}")

    # Regular? (all scores = m)
    is_regular = all(s == m for s in scores)
    print(f"  Regular: {is_regular}")

    # Number of 3-cycles
    c3 = 0
    for i in range(p):
        for j in range(i+1, p):
            for k in range(j+1, p):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    print(f"  c3 = {c3}")

    # Approximate path homology
    compute_path_homology_small(A, p, max_deg=1)

# Connection: BOTH Paley and Interval at p=7 are regular tournaments
# (all scores = 3). But they have different path homology!

print("\n" + "=" * 70)
print("THE PALEY VS INTERVAL HOMOLOGY QUESTION")
print("=" * 70)
print("""
Both Paley and Interval are REGULAR tournaments at p=7 (scores all = 3).
Both have the same c3 = 14 (3-cycle count).

But they have DIFFERENT H values (189 vs 175) and presumably
different path homology.

The Paley tournament T_7 has:
  β = (1, 0, 0, 0, 6, 0, 0)
  Euler char χ = 1 - 0 + 0 - 0 + 6 - 0 + 0 = 7

What about the Interval tournament? Key differences:
  - Interval has chirality 1.0 (breaks reflection)
  - Paley has chirality 1.0 at p≡3 mod 4 too
  - Both have |Aut| = p·m = 21

The β_4 = 6 for Paley comes from the eigenspace decomposition.
For Interval, we would need to compute the eigenspace dims of Ω.

SPECULATION: Perhaps β_4(Interval) < β_4(Paley) = 6?
If so, the "fewer homology generators" could reflect the
"more structured / less tangled" cycle structure of Interval.

This would connect:
  Lower β_4 ↔ more cycle disjointness ↔ higher α_2 ↔ higher H
  (AT LARGE p, not at p=7 where Paley still wins)
""")

# Let's compute the Euler characteristic for both via the
# alternating sum of Omega dimensions
print("\n--- Omega dimensions comparison ---")
for name, A in [("Interval", A_int), ("Paley", A_pal)]:
    print(f"\n  {name}:")
    # Omega_k = dim of k-allowed paths modulo boundaries
    # For a circulant tournament, compute Omega_k dims by brute force

    # k=0: n vertices
    omega_0 = p
    print(f"    Ω_0 = {omega_0}")

    # k=1: edges = n(n-1)/2 directed edges in tournament
    omega_1 = sum(1 for i in range(p) for j in range(p) if i != j and A[i][j])
    print(f"    Ω_1 = {omega_1} (= n(n-1)/2 = {p*(p-1)//2})")

    # k=2: allowed 2-paths (need face condition)
    # In GLMY: p-path (v0,...,vp) is allowed if removing any vertex
    # gives an allowed (p-1)-path. For p=2: (v0,v1,v2) is allowed if
    # (v0,v1), (v0,v2), (v1,v2) are all allowed 1-paths.
    # A 1-path (u,v) is "regular" iff there exists a directed edge u→v.
    # For tournaments: between any u,v, exactly one direction exists.
    # So (u,v) is a valid 1-path iff A[u][v]=1.
    # A 2-path (a,b,c) needs: A[a][b]=1, A[b][c]=1, and either A[a][c]=1 or A[c][a]=1.
    # For tournaments: this is ALWAYS satisfied (one of A[a][c], A[c][a] = 1).
    # So Omega_2 = all 2-step directed paths with distinct vertices.
    omega_2 = sum(1 for i in range(p) for j in range(p) if i != j and A[i][j]
                  for k in range(p) if k != i and k != j and A[j][k])
    print(f"    Ω_2 = {omega_2}")

    # The full Omega dims for circulant tournament at p=7 should be:
    # [7, 21, ?, ?, ?, ?, 1] (the 7 is all vertices, 1 is Hamiltonian)
    # Actually Ω_k counts the NUMBER of allowed k-paths, which is the
    # dimension of the path chain group A_k.
    # For Omega_k (GLMY), we need the kernel of the boundary restricted to A_k.
    # That's more complex.

    # For now just count the chain group sizes
    for k in [3, 4, 5, 6]:
        count = 0
        # k-paths: sequences v_0→v_1→...→v_k with all edges and distinct vertices
        # For k=3 at n=7: feasible to enumerate
        if k <= 4:
            from itertools import permutations
            for perm in permutations(range(p), k+1):
                valid = True
                for step in range(k):
                    if not A[perm[step]][perm[step+1]]:
                        valid = False
                        break
                if valid:
                    count += 1
            print(f"    A_{k} = {count}")

print("\nDONE.")
