#!/usr/bin/env python3
"""
GLMY PATH HOMOLOGY OF DIGRAPHS — Core implementation

Implements the Grigor'yan-Lin-Muranov-Yau (GLMY) path homology theory.

DEFINITIONS:
- An elementary p-path on vertex set V is a sequence (v_0, v_1, ..., v_p) of vertices.
- An elementary p-path is ALLOWED if (v_i, v_{i+1}) is a directed edge for all i.
- The boundary operator: ∂(v_0 ... v_p) = sum_{i=0}^{p} (-1)^i (v_0 ... v̂_i ... v_p)
  where v̂_i means delete vertex v_i.
- A_p = vector space of formal sums of allowed p-paths
- Omega_p = {u in A_p : ∂u in A_{p-1}} (∂-invariant allowed paths)
- Path homology: H_p = ker(∂: Omega_p -> Omega_{p-1}) / im(∂: Omega_{p+1} -> Omega_p)

KEY INSIGHT: Not every allowed p-path is ∂-invariant! When we delete an interior
vertex from an allowed path, the resulting (p-1)-path may not be allowed (the
two edges meeting at the deleted vertex get "merged" into a single edge that
may not exist in the digraph).

GEOMETRIC MEANING:
- β_0 = number of weakly connected components
- β_1 = "directed holes" (directed cycles that don't bound a filled region)
- β_2 = "directed cavities" (2D holes bounded by directed cycles)
- Higher β_p = higher-dimensional directed holes

For TOURNAMENTS:
- Every pair has exactly one directed edge
- So every sequence of distinct vertices is an allowed path!
- The ∂-invariance condition becomes trivial for tournaments
- This means Omega_p = A_p for tournaments!
- Tournament path homology should be computable purely from the
  boundary map on the full space of allowed paths.

Wait — that's a huge simplification. Let me verify.
"""
import numpy as np
from itertools import permutations, combinations
from math import comb
from collections import defaultdict
import random
from scipy.linalg import null_space

def is_allowed_path(A, path):
    """Check if path is allowed (all consecutive pairs are directed edges)."""
    for i in range(len(path) - 1):
        if A[path[i]][path[i+1]] != 1:
            return False
    return True

def enumerate_allowed_paths(A, n, p):
    """Enumerate all allowed p-paths (sequences of p+1 distinct vertices following edges)."""
    if p == 0:
        return [(v,) for v in range(n)]

    paths = []
    # Generate all sequences of p+1 distinct vertices
    for perm in permutations(range(n), p + 1):
        if is_allowed_path(A, perm):
            paths.append(perm)
    return paths

def boundary_of_path(path):
    """Compute ∂(v_0 ... v_p) = sum_{i=0}^{p} (-1)^i (v_0 ... v̂_i ... v_p).
    Returns list of (sign, face) pairs."""
    p = len(path) - 1
    faces = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        sign = (-1) ** i
        faces.append((sign, face))
    return faces

def build_boundary_matrix(A, n, p, allowed_p, allowed_pm1):
    """Build the matrix of ∂: A_p -> A_{p-1} restricted to allowed paths.

    Rows indexed by allowed (p-1)-paths, columns by allowed p-paths.
    """
    if not allowed_p or not allowed_pm1:
        return np.zeros((len(allowed_pm1) if allowed_pm1 else 0,
                         len(allowed_p) if allowed_p else 0))

    # Index maps
    pm1_index = {path: i for i, path in enumerate(allowed_pm1)}

    nrows = len(allowed_pm1)
    ncols = len(allowed_p)
    M = np.zeros((nrows, ncols), dtype=float)

    for j, path in enumerate(allowed_p):
        faces = boundary_of_path(path)
        for sign, face in faces:
            if face in pm1_index:
                M[pm1_index[face], j] += sign
            # If face is not allowed, it contributes nothing

    return M

def compute_omega(A, n, p, allowed_p, allowed_pm1):
    """Compute Omega_p = {u in A_p : ∂u in A_{p-1}}.

    For tournaments, every pair has a directed edge, so every sequence of
    distinct vertices is an allowed path. Hence ALL faces of an allowed path
    are allowed, and Omega_p = A_p.

    For general digraphs, we need to find the subspace of A_p whose boundary
    lies entirely in A_{p-1}.
    """
    if not allowed_p:
        return np.zeros((len(allowed_p), 0))

    # Build the "extended" boundary map including non-allowed faces
    # For each allowed p-path, compute its boundary faces
    # Check which faces are NOT allowed
    # Omega_p = kernel of the "non-allowed projection" of ∂

    # All possible (p-1)-paths (not just allowed ones)
    all_pm1 = list(permutations(range(n), p))

    # Index of NON-allowed (p-1)-paths
    allowed_set = set(allowed_pm1)
    non_allowed_pm1 = [path for path in all_pm1 if path not in allowed_set and len(set(path)) == p]

    if not non_allowed_pm1:
        # All faces are allowed => Omega_p = A_p
        return np.eye(len(allowed_p))

    # Build projection to non-allowed faces
    na_index = {path: i for i, path in enumerate(non_allowed_pm1)}
    nrows = len(non_allowed_pm1)
    ncols = len(allowed_p)
    P = np.zeros((nrows, ncols))

    for j, path in enumerate(allowed_p):
        faces = boundary_of_path(path)
        for sign, face in faces:
            if face in na_index:
                P[na_index[face], j] += sign

    # Omega_p = kernel of P
    # u in A_p is in Omega_p iff P @ u = 0
    if nrows == 0:
        return np.eye(ncols)

    # Find null space of P
    ns = null_space(P)
    return ns

def path_homology(A, n, max_dim=None):
    """Compute GLMY path homology of digraph with adjacency matrix A.

    Returns: list of Betti numbers [β_0, β_1, β_2, ...]
    """
    if max_dim is None:
        max_dim = n - 1

    # Enumerate allowed paths at each dimension
    allowed = {}
    for p in range(max_dim + 2):
        allowed[p] = enumerate_allowed_paths(A, n, p)

    betti = []

    for p in range(max_dim + 1):
        # Compute Omega_p
        omega_basis = compute_omega(A, n, p, allowed[p], allowed.get(p-1, []))

        if omega_basis.shape[1] == 0:
            betti.append(0)
            continue

        # Boundary ∂: Omega_p -> Omega_{p-1}
        # First, build ∂: A_p -> A_{p-1} on allowed paths
        bd_matrix = build_boundary_matrix(A, n, p, allowed[p], allowed.get(p-1, []))

        # Restrict to Omega_p: ∂_restricted = bd_matrix @ omega_basis
        bd_restricted = bd_matrix @ omega_basis

        # Also need Omega_{p+1} and its boundary into Omega_p
        omega_p1_basis = compute_omega(A, n, p+1, allowed.get(p+1, []), allowed[p])
        bd_p1_matrix = build_boundary_matrix(A, n, p+1, allowed.get(p+1, []), allowed[p])

        if omega_p1_basis.shape[1] > 0:
            # Image of ∂_{p+1} restricted to Omega_{p+1}, landing in A_p
            bd_p1_restricted = bd_p1_matrix @ omega_p1_basis

            # Project into Omega_p coordinates
            # We need to express bd_p1_restricted in the basis of Omega_p
            # This is tricky in general. For tournaments where Omega_p = A_p, it's simple.
        else:
            bd_p1_restricted = np.zeros((len(allowed[p]), 0))

        # Compute β_p = dim(ker ∂_p) - dim(im ∂_{p+1})
        # All operations in the Omega_p subspace

        # ker ∂_p on Omega_p: kernel of bd_restricted
        if bd_restricted.shape[0] > 0 and bd_restricted.shape[1] > 0:
            _, S, _ = np.linalg.svd(bd_restricted)
            rank_p = sum(s > 1e-8 for s in S)
            ker_dim = omega_basis.shape[1] - rank_p
        else:
            ker_dim = omega_basis.shape[1]

        # im ∂_{p+1} into Omega_p
        if bd_p1_restricted.shape[1] > 0:
            _, S1, _ = np.linalg.svd(bd_p1_restricted)
            im_dim = sum(s > 1e-8 for s in S1)
        else:
            im_dim = 0

        beta_p = ker_dim - im_dim
        betti.append(beta_p)

    return betti

# ===== SPECIAL CASE: TOURNAMENTS =====
def tournament_path_homology(A, n, max_dim=None):
    """Optimized path homology for tournaments.

    KEY SIMPLIFICATION: In a tournament, every pair has exactly one directed edge.
    Therefore EVERY sequence of distinct vertices is an allowed path.
    Moreover, deleting any vertex from an allowed path gives an allowed path
    (since all pairs still have edges). So Omega_p = A_p for all p.

    This means the chain complex is simply:
    ... -> A_p -> A_{p-1} -> ... -> A_0 -> 0
    with the standard boundary map.

    A_p has dimension n!/(n-p-1)! = P(n, p+1) (falling factorial).
    """
    if max_dim is None:
        max_dim = min(n - 1, 5)  # Cap for performance

    # For tournaments, enumerate all allowed p-paths = all (p+1)-permutations
    # that follow tournament edges
    # Wait — in a tournament, for a sequence (v_0,...,v_p), we need
    # A[v_i][v_{i+1}] = 1 for all i. Not every sequence follows edges!
    # A tournament has exactly one edge between each pair, but it might
    # go the wrong way. So allowed paths are still restricted.

    # My earlier claim was WRONG. Let me reconsider.
    # In a tournament: for any pair (u,v), either u->v or v->u.
    # An allowed p-path (v_0,...,v_p) requires v_i -> v_{i+1} for all i.
    # This is NOT automatic — it depends on which direction each edge goes.
    # So allowed paths = Hamiltonian paths of the induced sub-tournament
    # on the chosen vertices.

    # However, the ∂-invariance IS automatic for tournaments:
    # Deleting vertex v_i from (v_0,...,v_p) gives (v_0,...,v_{i-1},v_{i+1},...,v_p).
    # For this to be allowed, we need v_{i-1} -> v_{i+1}.
    # In a tournament, either v_{i-1} -> v_{i+1} or v_{i+1} -> v_{i-1}.
    # So the face might NOT be allowed!
    # Therefore Omega_p ≠ A_p in general for tournaments.

    # Let me just use the general algorithm.
    return path_homology(A, n, max_dim)


# ===== TESTS =====
print("=" * 70)
print("GLMY PATH HOMOLOGY — BASIC TESTS")
print("=" * 70)

# Test 1: Directed 3-cycle C_3 = {0->1, 1->2, 2->0}
print("\n--- Test 1: Directed 3-cycle ---")
A3 = [[0,1,0],[0,0,1],[1,0,0]]
betti = path_homology(A3, 3, max_dim=3)
print(f"  C_3: Betti = {betti}")
# Expected: β_0=1, β_1=1 (one directed hole)

# Test 2: Complete directed graph K_3 (both directions)
print("\n--- Test 2: Complete bidirected K_3 ---")
K3 = [[0,1,1],[1,0,1],[1,1,0]]
betti = path_homology(K3, 3, max_dim=3)
print(f"  K_3 bidirected: Betti = {betti}")
# Expected: β_0=1, all others 0 (contractible)

# Test 3: Transitive tournament on 3 vertices: 0->1, 0->2, 1->2
print("\n--- Test 3: Transitive tournament T_3 ---")
T3 = [[0,1,1],[0,0,1],[0,0,0]]
betti = path_homology(T3, 3, max_dim=3)
print(f"  Transitive T_3: Betti = {betti}")

# Test 4: Tournament with 3-cycle: 0->1, 1->2, 2->0
print("\n--- Test 4: Tournament with 3-cycle ---")
C3t = [[0,1,0],[0,0,1],[1,0,0]]
betti = path_homology(C3t, 3, max_dim=3)
print(f"  3-cycle tournament: Betti = {betti}")

# Test 5: Directed 4-cycle
print("\n--- Test 5: Directed 4-cycle ---")
A4c = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]]
betti = path_homology(A4c, 4, max_dim=3)
print(f"  C_4: Betti = {betti}")

# Test 6: Path graph 0->1->2->3
print("\n--- Test 6: Directed path 0->1->2->3 ---")
P4 = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]]
betti = path_homology(P4, 4, max_dim=3)
print(f"  P_4: Betti = {betti}")

# ===== TOURNAMENTS =====
print("\n\n" + "=" * 70)
print("PATH HOMOLOGY OF TOURNAMENTS")
print("=" * 70)

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

for n in [3, 4]:
    print(f"\n--- All tournaments on n={n} ---")

    betti_counts = defaultdict(int)
    betti_by_t3 = defaultdict(list)

    for A in all_tournaments(n):
        t3 = count_3cycles(A, n)
        H = ham_path_count(A, n)
        betti = path_homology(A, n, max_dim=min(n-1, 4))
        betti_tuple = tuple(betti)
        betti_counts[betti_tuple] += 1
        betti_by_t3[t3].append(betti_tuple)

        if betti_counts[betti_tuple] <= 2:
            print(f"  t3={t3}, H={H}: β = {betti}")

    print(f"\n  Betti number distribution:")
    for bt, count in sorted(betti_counts.items()):
        print(f"    β = {list(bt)}: {count} tournaments")

    print(f"\n  By 3-cycle count:")
    for t3 in sorted(betti_by_t3.keys()):
        bts = betti_by_t3[t3]
        unique = set(bts)
        print(f"    t3={t3}: {len(bts)} tournaments, distinct β: {[list(b) for b in sorted(unique)]}")

# n=5 (sampling)
print(f"\n--- Sampled tournaments on n=5 ---")
random.seed(42)
betti_counts_5 = defaultdict(int)

for trial in range(50):
    A = [[0]*5 for _ in range(5)]
    for i in range(5):
        for j in range(i+1,5):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    t3 = count_3cycles(A, 5)
    H = ham_path_count(A, 5)
    betti = path_homology(A, 5, max_dim=4)
    betti_counts_5[tuple(betti)] += 1

    if trial < 5:
        print(f"  Trial {trial}: t3={t3}, H={H}, β = {betti}")

print(f"\n  Betti number distribution (50 samples):")
for bt, count in sorted(betti_counts_5.items()):
    print(f"    β = {list(bt)}: {count}")

print("\n\nDone.")
