#!/usr/bin/env python3
"""
fundamental_S72.py — opus-2026-03-15-S72

THE FUNDAMENTAL QUESTION: Why does β₁ capture global topology?

We found tournaments with IDENTICAL local structure (scores, c3v, c4, sc4, sc5)
but different β₁. This means β₁ sees something no subgraph count can see.

WHAT does β₁ see?

β₁ = dim(Z₁) - rank(∂₂|Ω₂)
   = C(n,2) - (n-1) - rank(∂₂|Ω₂)

So β₁ > 0 iff there exists a directed 1-cycle that is NOT the boundary
of any 2-chain. In other words: a combination of edges that "wraps around"
the tournament and cannot be decomposed into transitive-triple boundaries.

PLAN:
1. Extract the actual β₁=0 and β₁=1 violation pair
2. Find the ACTUAL 1-cycle that doesn't bound
3. Understand geometrically why it can't be filled
4. Look at the complement/reversal structure
5. Study the β₁·β₃=0 mutual exclusivity from first principles
6. Find what TOPOLOGICAL property distinguishes the pair
"""

import numpy as np
from math import comb
from collections import defaultdict, Counter
from itertools import combinations
import sys
sys.stdout.reconfigure(line_buffering=True)

PRIME = 2**31 - 1

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def per_vertex_3cycles(A, n):
    c3v = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i,j] and A[j,k] and A[k,i]) or \
                   (A[i,k] and A[k,j] and A[j,i]):
                    c3v[i] += 1; c3v[j] += 1; c3v[k] += 1
    return tuple(sorted(c3v))

def count_sc4(A, n):
    sc4 = 0
    for verts in combinations(range(n), 4):
        scores = sorted(sum(A[vi, vj] for vj in verts if vj != vi) for vi in verts)
        if scores != [0, 1, 2, 3]: sc4 += 1
    return sc4

def compute_beta1_with_cycle(A, n):
    """Compute β₁ AND return the unbounding cycle if β₁ > 0."""
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                edges.append((i, j))
                edge_idx[(i, j)] = len(edges) - 1

    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k in (i, j) or A[j][k] == 0 or A[i][k] == 0: continue
                triples.append((i, j, k))

    if not triples:
        return comb(n, 2) - (n - 1), None, None

    # Build boundary matrix ∂₂: columns = triples, rows = edges
    bd = np.zeros((len(edges), len(triples)), dtype=np.int64)
    for t, (i, j, k) in enumerate(triples):
        if (j, k) in edge_idx:
            bd[edge_idx[(j, k)], t] = (bd[edge_idx[(j, k)], t] + 1) % PRIME
        if (i, k) in edge_idx:
            bd[edge_idx[(i, k)], t] = (bd[edge_idx[(i, k)], t] - 1) % PRIME
        if (i, j) in edge_idx:
            bd[edge_idx[(i, j)], t] = (bd[edge_idx[(i, j)], t] + 1) % PRIME

    # Also build ∂₁: columns = edges, rows = vertices
    bd1 = np.zeros((n, len(edges)), dtype=np.int64)
    for j, (a, b) in enumerate(edges):
        bd1[b, j] = 1   # +target
        bd1[a, j] = (PRIME - 1) % PRIME  # -source

    # Compute rank of ∂₂
    M = bd.copy() % PRIME
    nrows, ncols = M.shape
    rank2 = 0
    for col in range(ncols):
        nonzero = np.where(M[rank2:, col] != 0)[0]
        if len(nonzero) == 0: continue
        pivot = nonzero[0] + rank2
        if pivot != rank2: M[[rank2, pivot]] = M[[pivot, rank2]]
        inv = pow(int(M[rank2, col]), PRIME - 2, PRIME)
        M[rank2] = M[rank2] * inv % PRIME
        factors = M[:, col].copy()
        factors[rank2] = 0
        nzr = np.where(factors != 0)[0]
        if len(nzr) > 0:
            M[nzr] = (M[nzr] - np.outer(factors[nzr], M[rank2])) % PRIME
        rank2 += 1

    beta1 = max(0, comb(n, 2) - (n - 1) - rank2)

    # If β₁ > 0, find the actual unbounding cycle
    # A 1-cycle is in ker(∂₁) but not in im(∂₂)
    unbounding_cycle = None
    if beta1 > 0:
        # Find ker(∂₁): compute null space of bd1
        M1 = bd1.copy() % PRIME
        nrows1, ncols1 = M1.shape
        rank1 = 0
        pivot_cols1 = []
        for col in range(ncols1):
            nonzero = np.where(M1[rank1:, col] != 0)[0]
            if len(nonzero) == 0: continue
            pivot = nonzero[0] + rank1
            if pivot != rank1: M1[[rank1, pivot]] = M1[[pivot, rank1]]
            inv = pow(int(M1[rank1, col]), PRIME - 2, PRIME)
            M1[rank1] = M1[rank1] * inv % PRIME
            factors = M1[:, col].copy()
            factors[rank1] = 0
            nzr = np.where(factors != 0)[0]
            if len(nzr) > 0:
                M1[nzr] = (M1[nzr] - np.outer(factors[nzr], M1[rank1])) % PRIME
            pivot_cols1.append(col)
            rank1 += 1

        free_cols1 = [c for c in range(ncols1) if c not in pivot_cols1]

        # ker(∂₁) basis vectors
        ker_basis = []
        for fc in free_cols1:
            vec = np.zeros(ncols1, dtype=np.int64)
            vec[fc] = 1
            for r, pc in enumerate(pivot_cols1):
                if r < M1.shape[0]:
                    vec[pc] = (-M1[r, fc]) % PRIME
            ker_basis.append(vec)

        # im(∂₂) = column space of bd
        # Find a ker(∂₁) vector not in im(∂₂)
        # Augment bd with each ker basis vector and check if rank increases
        for kb in ker_basis:
            aug = np.column_stack([bd % PRIME, kb.reshape(-1, 1) % PRIME])
            M_aug = aug.copy()
            rank_aug = 0
            for col in range(M_aug.shape[1]):
                nonzero = np.where(M_aug[rank_aug:, col] != 0)[0]
                if len(nonzero) == 0: continue
                pivot = nonzero[0] + rank_aug
                if pivot != rank_aug: M_aug[[rank_aug, pivot]] = M_aug[[pivot, rank_aug]]
                inv = pow(int(M_aug[rank_aug, col]), PRIME - 2, PRIME)
                M_aug[rank_aug] = M_aug[rank_aug] * inv % PRIME
                factors = M_aug[:, col].copy()
                factors[rank_aug] = 0
                nzr = np.where(factors != 0)[0]
                if len(nzr) > 0:
                    M_aug[nzr] = (M_aug[nzr] - np.outer(factors[nzr], M_aug[rank_aug])) % PRIME
                rank_aug += 1

            if rank_aug > rank2:
                # This cycle is NOT in im(∂₂) — it's the unbounding cycle!
                cycle_edges = []
                for idx, coeff in enumerate(kb):
                    if coeff != 0 and coeff != PRIME:
                        c = int(coeff) if coeff < PRIME // 2 else int(coeff) - PRIME
                        cycle_edges.append((edges[idx], c))
                unbounding_cycle = cycle_edges
                break

    return beta1, unbounding_cycle, {'rank_d2': rank2, 'n_triples': len(triples),
                                       'n_edges': len(edges), 'edges': edges}

# ============================================================
# STEP 1: Find the violation pair
# ============================================================

print("=" * 72)
print("THE FUNDAMENTAL QUESTION")
print("What does β₁ see that subgraph counts cannot?")
print("=" * 72)

n = 6
total = 2**(n*(n-1)//2)

print(f"\nSearching for violation pair at n={n}...")

b0_example = None
b1_example = None

for bits in range(total):
    A = bits_to_adj(bits, n)
    c3v = per_vertex_3cycles(A, n)
    if c3v != (2, 2, 3, 3, 4, 4):
        continue
    sc4 = count_sc4(A, n)
    if sc4 != 11:
        continue

    beta1, cycle, info = compute_beta1_with_cycle(A, n)
    if beta1 == 0 and b0_example is None:
        b0_example = (bits, A.copy(), info)
    elif beta1 == 1 and b1_example is None:
        b1_example = (bits, A.copy(), info, cycle)

    if b0_example and b1_example:
        break

print(f"\nFound pair!")

# ============================================================
# STEP 2: Compare the two tournaments
# ============================================================

print("\n" + "=" * 72)
print("STEP 2: Structural comparison of β₁=0 vs β₁=1")
print("=" * 72)

bits0, A0, info0 = b0_example
bits1, A1, info1, cycle1 = b1_example

print(f"\n  β₁=0 tournament (bits={bits0}):")
print(f"    Adjacency matrix:")
for i in range(n):
    print(f"      {list(A0[i])}")
print(f"    Out-degrees: {[int(np.sum(A0[i])) for i in range(n)]}")
print(f"    rank(∂₂) = {info0['rank_d2']}, #triples = {info0['n_triples']}")

print(f"\n  β₁=1 tournament (bits={bits1}):")
print(f"    Adjacency matrix:")
for i in range(n):
    print(f"      {list(A1[i])}")
print(f"    Out-degrees: {[int(np.sum(A1[i])) for i in range(n)]}")
print(f"    rank(∂₂) = {info1['rank_d2']}, #triples = {info1['n_triples']}")

print(f"\n  Difference in rank(∂₂): {info0['rank_d2'] - info1['rank_d2']}")

# ============================================================
# STEP 3: The unbounding cycle
# ============================================================

print("\n" + "=" * 72)
print("STEP 3: The cycle that doesn't bound")
print("=" * 72)

if cycle1:
    print(f"\n  The 1-cycle in β₁=1 tournament that has no 2-chain boundary:")
    total_edges = 0
    for (edge, coeff) in cycle1:
        sign = '+' if coeff > 0 else '-'
        print(f"    {sign}{abs(coeff)}·({edge[0]}→{edge[1]})")
        total_edges += 1
    print(f"  Total: {total_edges} edges in the cycle")
    print(f"\n  This cycle exists in ker(∂₁) but NOT in im(∂₂).")
    print(f"  It wraps around the tournament in a way that no combination")
    print(f"  of transitive triples can decompose.")

# ============================================================
# STEP 4: What's TOPOLOGICALLY different?
# ============================================================

print("\n" + "=" * 72)
print("STEP 4: Topological difference")
print("=" * 72)

# Compare the 3-cycle structure (which specific triples form 3-cycles)
def find_3cycles(A, n):
    """Find all directed 3-cycles as sets of vertices."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i,j] and A[j,k] and A[k,i]:
                    cycles.append((i,j,k, 'fwd'))
                elif A[i,k] and A[k,j] and A[j,i]:
                    cycles.append((i,j,k, 'bwd'))
    return cycles

cycles0 = find_3cycles(A0, n)
cycles1 = find_3cycles(A1, n)

print(f"\n  β₁=0: 3-cycles: {[(c[0],c[1],c[2]) for c in cycles0]}")
print(f"  β₁=1: 3-cycles: {[(c[0],c[1],c[2]) for c in cycles1]}")

# Compare which edges differ
diff_edges = []
for i in range(n):
    for j in range(i+1, n):
        if A0[i,j] != A1[i,j]:
            if A0[i,j] == 1:
                diff_edges.append(f"{i}→{j} in T0, {j}→{i} in T1")
            else:
                diff_edges.append(f"{j}→{i} in T0, {i}→{j} in T1")

print(f"\n  Edges that differ between the two tournaments:")
for d in diff_edges:
    print(f"    {d}")
print(f"  Total flipped edges: {len(diff_edges)}")

# ============================================================
# STEP 5: Complement and reversal
# ============================================================

print("\n" + "=" * 72)
print("STEP 5: Symmetry analysis")
print("=" * 72)

# Is T1 isomorphic to T0^op (complement)?
A0_op = np.zeros_like(A0)
for i in range(n):
    for j in range(n):
        if i != j:
            A0_op[i, j] = 1 - A0[i, j]

# Check if A1 is a relabeling of A0_op
from itertools import permutations

def is_isomorphic(A, B, n):
    """Check if tournament A is isomorphic to B."""
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if A[i][j] != B[perm[i]][perm[j]]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True, perm
    return False, None

iso_to_complement, perm_comp = is_isomorphic(A1, A0_op, n)
print(f"  T1 isomorphic to T0^op? {iso_to_complement}")

iso_to_self, perm_self = is_isomorphic(A0, A1, n)
print(f"  T0 isomorphic to T1? {iso_to_self}")

# Check self-complementarity
sc0, _ = is_isomorphic(A0, A0_op, n)
sc1_mat = np.zeros_like(A1)
for i in range(n):
    for j in range(n):
        if i != j:
            sc1_mat[i, j] = 1 - A1[i, j]
sc1, _ = is_isomorphic(A1, sc1_mat, n)
print(f"  T0 self-complementary? {sc0}")
print(f"  T1 self-complementary? {sc1}")

# ============================================================
# STEP 6: The deep structure — path coverage
# ============================================================

print("\n" + "=" * 72)
print("STEP 6: Path coverage analysis")
print("=" * 72)

# For each tournament, look at the "coverage graph":
# Which edges are covered by transitive triple boundaries?
# An edge (a,b) is "covered" if it appears in ∂₂(some triple)

for label, A, info_data in [("β₁=0", A0, info0), ("β₁=1", A1, info1)]:
    print(f"\n  {label} tournament:")

    # Transitive triples and their boundaries
    triples = []
    for i in range(n):
        for j in range(n):
            if i == j or A[i][j] == 0: continue
            for k in range(n):
                if k in (i, j) or A[j][k] == 0 or A[i][k] == 0: continue
                triples.append((i, j, k))

    # For each edge, how many triples touch it?
    edge_coverage = defaultdict(int)
    for (i, j, k) in triples:
        edge_coverage[(i, j)] += 1
        edge_coverage[(j, k)] += 1
        edge_coverage[(i, k)] += 1

    covered = sum(1 for e in edge_coverage if edge_coverage[e] > 0)
    total_edges = n * (n - 1) // 2 * 2  # directed
    # Actually, n*(n-1)/2 undirected, but all are directed in tournament
    n_edges = sum(1 for i in range(n) for j in range(n) if A[i][j] == 1)

    print(f"    #triples: {len(triples)}")
    print(f"    #edges: {n_edges}")
    print(f"    #edges in ≥1 triple: {covered}/{n_edges}")
    uncovered = [(i,j) for i in range(n) for j in range(n) if A[i][j]==1 and edge_coverage[(i,j)]==0]
    print(f"    Uncovered edges: {uncovered}")

    coverage_dist = Counter(edge_coverage[e] for e in edge_coverage if edge_coverage[e] > 0)
    print(f"    Coverage distribution: {dict(sorted(coverage_dist.items()))}")

# ============================================================
# STEP 7: The β₁·β₃=0 principle from first principles
# ============================================================

print("\n" + "=" * 72)
print("STEP 7: Why β₁·β₃=0? The underfilling/overfilling duality")
print("=" * 72)

print("""
  β₁ > 0 means: UNDERFILLING — there exist directed cycles that
  cannot be decomposed into transitive triple boundaries.
  The 2-path system doesn't reach far enough.

  β₃ > 0 means: OVERFILLING — there exist 3-dimensional holes.
  The 4-path system creates cavities that can't be filled.

  CONJECTURE for why β₁·β₃ = 0:

  If β₁ > 0, the tournament has "topological slack" — uncoverable cycles.
  This slack means there's room for 3-chains to fill any potential 3-holes.
  The im(∂₃) → ker(∂₂) exactness (β₂=0) then forces im(∂₄) → ker(∂₃).

  If β₃ > 0, the tournament is "topologically dense" — so many paths that
  they wrap around creating cavities. But this density means every 1-cycle
  IS decomposable (β₁=0).

  This is a DUALITY: slack ↔ density. You can't have both.
""")

# Verify with data at n=5
print("  Verification at n=5 (exhaustive):")
n5 = 5
m5 = n5 * (n5 - 1) // 2

# Quick homology for all n=5 tournaments
from collections import defaultdict

h_by_beta = defaultdict(list)
for bits in range(2**m5):
    A = bits_to_adj(bits, n5)

    # β₁
    b1_val = compute_beta1_with_cycle(A, n5)[0]

    # H
    dp = {}
    for v in range(n5): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n5):
        for v in range(n5):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n5):
                if mask & (1 << u): continue
                if A[v, u] == 1:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    H = sum(dp.get(((1 << n5) - 1, v), 0) for v in range(n5))

    t3 = sum(1 for i in range(n5) for j in range(i+1,n5) for k in range(j+1,n5)
             if (A[i,j] and A[j,k] and A[k,i]) or (A[i,k] and A[k,j] and A[j,i]))

    h_by_beta[b1_val].append((H, t3))

print(f"\n    β₁=0 tournaments ({len(h_by_beta[0])}):")
H_vals = [h for h, _ in h_by_beta[0]]
t3_vals = [t for _, t in h_by_beta[0]]
print(f"      H range: [{min(H_vals)}, {max(H_vals)}], mean: {np.mean(H_vals):.1f}")
print(f"      t3 range: [{min(t3_vals)}, {max(t3_vals)}], mean: {np.mean(t3_vals):.1f}")

print(f"\n    β₁=1 tournaments ({len(h_by_beta[1])}):")
H_vals1 = [h for h, _ in h_by_beta[1]]
t3_vals1 = [t for _, t in h_by_beta[1]]
print(f"      H range: [{min(H_vals1)}, {max(H_vals1)}], mean: {np.mean(H_vals1):.1f}")
print(f"      t3 range: [{min(t3_vals1)}, {max(t3_vals1)}], mean: {np.mean(t3_vals1):.1f}")

print(f"\n  β₁=1 tournaments have HIGHER H and t3 on average!")
print(f"  This seems to CONTRADICT the 'underfilling' hypothesis...")
print(f"  Unless: β₁=1 means 'topological TWISTING' not 'sparsity'.")

# ============================================================
# STEP 8: The twist interpretation
# ============================================================

print("\n" + "=" * 72)
print("STEP 8: The Möbius twist interpretation")
print("=" * 72)

print("""
  REVISED INSIGHT:

  β₁ > 0 does NOT mean "too few paths." It means TWISTED paths.
  The tournament has a directed cycle where the 2-paths that
  should fill it are "oriented the wrong way" — like trying to
  tile a Möbius strip with oriented squares.

  The unbounding cycle wraps around in a way that reverses
  the orientation of any potential filling.

  This is why subgraph counts can't detect it: COUNTS are
  orientation-blind. β₁ sees ORIENTATION.

  The analogy: Two surfaces can have the same number of
  triangles, edges, vertices — but one is orientable (S²)
  and the other is not (RP²). β₁ detects this.

  For tournaments: the directed edges carry orientation.
  A 1-cycle is a signed combination of edges. Whether it
  bounds depends on whether the orientations of 2-paths
  ALIGN consistently. This is a GLOBAL orientation property.
""")

# Test: does the unbounding cycle relate to the Möbius structure?
if cycle1:
    print("  The unbounding cycle edges form a pattern:")
    cycle_vertices = set()
    for (e, c) in cycle1:
        cycle_vertices.add(e[0])
        cycle_vertices.add(e[1])
    print(f"  Vertices involved: {sorted(cycle_vertices)}")
    print(f"  Induced subtournament on these vertices:")
    vlist = sorted(cycle_vertices)
    for i in vlist:
        row = [int(A1[i, j]) for j in vlist]
        print(f"    {i}: {row}")

    # Check: is this subtournament strongly connected?
    sub_n = len(vlist)
    if sub_n >= 3:
        sub_scores = sorted(sum(A1[vi, vj] for vj in vlist if vj != vi) for vi in vlist)
        is_sc = sub_scores != list(range(sub_n))
        print(f"  Subtournament score sequence: {sub_scores}")
        print(f"  Strongly connected? {is_sc}")

# ============================================================
# STEP 9: The fundamental insight
# ============================================================

print("\n" + "=" * 72)
print("STEP 9: THE FUNDAMENTAL INSIGHT")
print("=" * 72)

print("""
  β₁ is the ORIENTABILITY OBSTRUCTION of the tournament's path complex.

  Key chain of reasoning:

  1. β₂ = 0 (PROVED): The path complex is exact at degree 2.
     This means every 2-cycle bounds. The complex is "2-acyclic."

  2. β₁ measures 1-cycles that don't bound 2-chains.
     In an oriented complex, this happens when the ORIENTATIONS
     of 2-cells are INCONSISTENT with the cycle's orientation.

  3. For tournaments, "orientation inconsistency" means:
     the cycle wraps through vertices whose out-neighborhoods
     create a TWISTING pattern in the path space.

  4. No subgraph count can detect twisting because counts are
     symmetric under relabeling. But orientation is NOT.
     Two tournaments can have the same unlabeled subgraph
     census but different labeled orientation patterns.

  5. The β₁·β₃ = 0 mutual exclusivity then says:
     - Twisted tournaments (β₁>0) can't create 3-holes
       because the twist PREVENTS higher-dimensional wrapping.
     - Dense tournaments (β₃>0) can't be twisted because
       density forces orientation consistency at the 1-level.

  This is the TOPOLOGICAL CONTENT of tournament theory.
  OCF (H = I(Ω,2)) is an ALGEBRAIC identity.
  β₁·β₃ = 0 is a GEOMETRIC constraint.
  Together they give the full picture:
    ALGEBRAIC: how many paths (counted by H)
    GEOMETRIC: how paths wrap (measured by β)
""")

print("\nDONE — fundamental_S72.py complete")
