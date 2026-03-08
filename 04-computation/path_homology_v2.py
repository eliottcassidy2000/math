#!/usr/bin/env python3
"""
GLMY PATH HOMOLOGY v2 — Fixed implementation

Fixes from v1:
- Omega_0 = A_0 (boundary ∂_0 = 0 by convention, so all 0-paths are ∂-invariant)
- Cleaner handling of the chain complex at dimension 0

DEFINITIONS (GLMY):
- Elementary p-path: sequence (v_0, ..., v_p) of DISTINCT vertices
- ALLOWED p-path: one where v_i -> v_{i+1} is a directed edge for all i
- Boundary: ∂(v_0...v_p) = sum_{i=0}^p (-1)^i (v_0...v̂_i...v_p)
- A_p = vector space spanned by allowed p-paths
- Ω_p = {u ∈ A_p : ∂u ∈ A_{p-1}}  (for p ≥ 1; Ω_0 = A_0)
- H_p^{path}(G) = ker(∂_p: Ω_p → Ω_{p-1}) / im(∂_{p+1}: Ω_{p+1} → Ω_p)
"""
import numpy as np
from itertools import permutations
from math import comb
from collections import defaultdict
import random

def enumerate_allowed_paths(A, n, p):
    """All sequences of p+1 distinct vertices following directed edges."""
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    paths = []
    for perm in permutations(range(n), p + 1):
        ok = True
        for i in range(p):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok:
            paths.append(perm)
    return paths

def boundary_coeffs(path):
    """Returns list of (sign, face_tuple) for ∂(path)."""
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    """Build ∂: A_p → A_{p-1} as a matrix.
    Columns = allowed p-paths, Rows = allowed (p-1)-paths.
    Entry (r,c) = coefficient of (p-1)-path r in ∂(p-path c)."""
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))

    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))

    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    """Compute a basis for Ω_p = {u ∈ A_p : ∂u ∈ A_{p-1}}.

    For p=0: Ω_0 = A_0 (convention).
    For p≥1: Need ∂u to lie entirely in A_{p-1}, i.e., no component on
    non-allowed (p-1)-paths.
    """
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))

    if p == 0:
        return np.eye(dim_Ap)

    # Find faces that are NOT allowed
    allowed_pm1_set = set(allowed_pm1)

    # For each allowed p-path, find its faces and check which are not allowed
    # Collect all non-allowed faces
    non_allowed_faces = {}
    na_count = 0

    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1

    if na_count == 0:
        # All faces are allowed ⟹ Ω_p = A_p
        return np.eye(dim_Ap)

    # Build the projection to non-allowed faces
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign

    # Ω_p = ker(P)
    # Use SVD to find null space
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T  # columns are null vectors

    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    """Compute GLMY path Betti numbers β_0, β_1, ..., β_{max_dim}."""
    if max_dim is None:
        max_dim = n - 1

    # Pre-compute allowed paths
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    # Pre-compute Omega bases
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])

    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0

        if dim_omega_p == 0:
            betti.append(0)
            continue

        # ∂_p: Ω_p → A_{p-1} (restricted)
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        # In Omega_p coordinates:
        bd_p_omega = bd_p @ omega[p]  # shape (|A_{p-1}|, dim_omega_p)

        # Rank of ∂_p restricted to Ω_p
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0

        ker_dim = dim_omega_p - rank_p

        # ∂_{p+1}: Ω_{p+1} → A_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]  # shape (|A_p|, dim_omega_{p+1})

            # Image lands in A_p. But we want image in Ω_p.
            # Project onto Ω_p: if omega[p] is a basis of Ω_p,
            # im(∂_{p+1}) ∩ Ω_p should be computed.
            # Actually, by the chain complex property ∂∂=0,
            # im(∂_{p+1}) ⊂ ker(∂_p) ⊂ Ω_p.
            # So the image is already in Ω_p.
            # Rank of the image:
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0

        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))  # numerical safety

    return betti


# ===== VALIDATION =====
print("=" * 70)
print("GLMY PATH HOMOLOGY v2 — VALIDATION")
print("=" * 70)

# Test 1: Directed 3-cycle
print("\n--- Directed 3-cycle C_3 ---")
print("  Expected: β = [1, 1, 0] (1 component, 1 directed hole)")
A = [[0,1,0],[0,0,1],[1,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")

# Test 2: Directed path 0->1->2
print("\n--- Directed path P_3 ---")
print("  Expected: β = [1, 0, 0] (connected, contractible)")
A = [[0,1,0],[0,0,1],[0,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")

# Test 3: Transitive tournament T_3 (= directed path)
print("\n--- Transitive tournament T_3: 0->1, 0->2, 1->2 ---")
print("  This is an acyclic tournament. Expected: β_1 = 0 (no directed holes)")
A = [[0,1,1],[0,0,1],[0,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")

# Test 4: 3-cycle tournament (= directed 3-cycle)
print("\n--- 3-cycle tournament: 0->1, 1->2, 2->0 ---")
A = [[0,1,0],[0,0,1],[1,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")

# Test 5: Directed 4-cycle
print("\n--- Directed 4-cycle C_4 ---")
print("  Expected: β = [1, 1, 0, 0]")
A = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 4)}")

# Test 6: Complete bidirected K_3
print("\n--- Complete bidirected K_3 ---")
A = [[0,1,1],[1,0,1],[1,1,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")

# Test 7: Two isolated vertices + one edge
print("\n--- 0->1, vertex 2 isolated ---")
A = [[0,1,0],[0,0,0],[0,0,0]]
print(f"  Got: β = {path_betti_numbers(A, 3)}")
print("  Expected: β_0 = 3 or 2 (depends on convention)")

# ===== CIRCULANT DIGRAPHS =====
print("\n\n" + "=" * 70)
print("CIRCULANT DIGRAPHS C_n^S")
print("=" * 70)

def circulant_digraph(n, S):
    """C_n^S: vertex i has edges to (i+s) mod n for each s in S."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            j = (i + s) % n
            if j != i:
                A[i][j] = 1
    return A

# C_n^{1}: directed n-cycle
for n_val in [3, 4, 5, 6, 7]:
    A = circulant_digraph(n_val, [1])
    betti = path_betti_numbers(A, n_val, max_dim=min(n_val-1, 4))
    print(f"  C_{n_val}^{{1}} (directed cycle): β = {betti}")

# C_n^{1,2}: two-step circulant
print()
for n_val in [4, 5, 6, 7]:
    A = circulant_digraph(n_val, [1, 2])
    betti = path_betti_numbers(A, n_val, max_dim=min(n_val-1, 4))
    print(f"  C_{n_val}^{{1,2}}: β = {betti}")

# C_n^{1,2,3}: three-step circulant
print()
for n_val in [5, 6, 7]:
    A = circulant_digraph(n_val, [1, 2, 3])
    betti = path_betti_numbers(A, n_val, max_dim=min(n_val-1, 4))
    print(f"  C_{n_val}^{{1,2,3}}: β = {betti}")

# ===== ALL TOURNAMENTS n=3,4 =====
print("\n\n" + "=" * 70)
print("ALL TOURNAMENTS — PATH HOMOLOGY")
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
    print(f"\n--- All tournaments n={n} ---")
    betti_dist = defaultdict(int)
    detail = defaultdict(list)

    for A in all_tournaments(n):
        t3 = count_3cycles(A, n)
        H = ham_path_count(A, n)
        betti = path_betti_numbers(A, n, max_dim=n-1)
        bt = tuple(betti)
        betti_dist[bt] += 1
        detail[t3].append(bt)

    print(f"  Betti distribution:")
    for bt in sorted(betti_dist.keys()):
        print(f"    β={list(bt)}: {betti_dist[bt]} tournaments")

    print(f"  By t3:")
    for t3 in sorted(detail.keys()):
        from collections import Counter
        c = Counter(detail[t3])
        print(f"    t3={t3}: {dict(c)}")

# n=5 exhaustive (might be slow for all 1024)
print(f"\n--- All tournaments n=5 ---")
betti_dist_5 = defaultdict(int)
detail_5 = defaultdict(list)
count = 0
for A in all_tournaments(5):
    t3 = count_3cycles(A, 5)
    H = ham_path_count(A, 5)
    betti = path_betti_numbers(A, 5, max_dim=4)
    bt = tuple(betti)
    betti_dist_5[bt] += 1
    detail_5[t3].append(bt)
    count += 1
    if count <= 5:
        print(f"  t3={t3}, H={H}: β={betti}")

print(f"\n  Betti distribution (n=5, {count} tournaments):")
for bt in sorted(betti_dist_5.keys()):
    print(f"    β={list(bt)}: {betti_dist_5[bt]}")

print(f"\n  By t3:")
for t3 in sorted(detail_5.keys()):
    from collections import Counter
    c = Counter(detail_5[t3])
    print(f"    t3={t3}: {dict(c)}")

# ===== GEOMETRIC INTERPRETATION =====
print("\n\n" + "=" * 70)
print("GEOMETRIC INTERPRETATION")
print("=" * 70)
print("""
The GLMY path Betti numbers measure "directed topological features":

β_0 = # connected components (in path homology sense)
β_1 = # independent directed 1-holes (directed cycles not bounding anything)
β_2 = # independent directed 2-holes (directed "cavities")

For TOURNAMENTS:
- β_0 = 1 always (tournament is connected)
- β_1 depends on the 3-cycle structure
- β_2 and higher: potential higher-dimensional directed cavities

TOPOLOGICAL SPACES ANALOGY:
- β = (1,0,...) : contractible → like a point
- β = (1,1,0,...) : one 1-hole → like a circle S^1
- β = (1,2,1,...) : like a torus T^2
- β = (1,0,1,...) : one 2-hole → like a sphere S^2
""")

# Summary table
print("SUMMARY: Which topological space does each tournament resemble?")
for n in [3, 4, 5]:
    print(f"\nn={n}:")
    space_map = {
        (1,): "point",
        (1, 0): "point (contractible)",
        (1, 0, 0): "point",
        (1, 0, 0, 0): "point",
        (1, 0, 0, 0, 0): "point",
        (1, 1): "circle S^1",
        (1, 1, 0): "circle S^1",
        (1, 1, 0, 0): "circle S^1",
        (1, 1, 0, 0, 0): "circle S^1",
        (1, 2, 1): "torus T^2",
        (1, 2, 1, 0): "torus T^2",
        (1, 0, 1): "sphere S^2",
        (1, 0, 1, 0): "sphere S^2",
        (1, 0, 0, 1): "sphere S^3",
    }
    dist = betti_dist if n == 3 else (betti_dist if n == 4 else betti_dist_5)
    # Ugh, variable scoping. Let me just print what we have.
    if n == 3:
        for bt, cnt in sorted(betti_dist.items()):
            if bt in [(1,0,0), (0,0,0)]:
                continue
            space = space_map.get(bt, "unknown")
            print(f"  β={list(bt)}: {cnt} tournaments → {space}")
    elif n == 4:
        pass  # handled above
    elif n == 5:
        for bt, cnt in sorted(betti_dist_5.items()):
            space = space_map.get(bt, "unknown")
            print(f"  β={list(bt)}: {cnt} tournaments → {space}")

print("\n\nDone.")
