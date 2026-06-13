#!/usr/bin/env python3
"""
n9_eigenspace_betti.py — Eigenspace decomposition of β₅=10 at n=9

The n=7 Paley maximizer has β₄=6 = (p-1), from Z/7Z eigenspace decomposition:
  - k=0 (trivial): contributes β₄=0 (just β₀=1)
  - k=1,...,6 (non-trivial): each contributes β₄=1

At n=9, the maximizer is NOT Paley (9≡1 mod 4), but circulant maximizers
exist on Z₉ with connection set S={1,5,6,7}. β₅=10.

Question: How does Z/9Z decompose β₅=10 across eigenspaces?
  - Z/9Z has 9 irreducible reps over C (characters ω^k, ω=e^{2πi/9})
  - Over Q, these group into: trivial (dim 1), order-3 (dim 2), order-9 (dim 6)
  - If each of 9 eigenspaces contributes 1: total=9≠10. Something else is going on.

Also tests the Betti formula conjecture:
  - C((n+1)/2, 2) = C(5,2) = 10 (matches!)
  - Predicts n=11: β₆=15

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from scipy import linalg
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def build_circulant_tournament(n, S):
    """Build circulant tournament on Z_n with connection set S."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S:
                A[i][j] = 1
    return A

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def get_allowed_paths(A, n, length):
    """Get all allowed p-paths (sequences of distinct vertices with edges)."""
    if length == 0:
        return [(v,) for v in range(n)]

    paths = []
    # Generate all sequences of (length+1) distinct vertices
    def generate(path, used):
        if len(path) == length + 1:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in used:
                continue
            if A[last][v]:  # Must have edge from last to v
                generate(path + [v], used | {v})

    for start in range(n):
        generate([start], {start})

    return paths

def boundary_matrix(A, n, paths_p, paths_pm1, p):
    """Build the boundary matrix ∂_p : Ω_p → Ω_{p-1}.

    The GLMY boundary is: ∂_p(v_0,...,v_p) = sum_{i=1}^{p-1} (-1)^i (v_0,...,v̂_i,...,v_p)
    where v̂_i means delete v_i, BUT only if the resulting path is allowed.
    """
    path_to_idx = {p: i for i, p in enumerate(paths_pm1)}
    m = len(paths_p)
    n_cols = len(paths_pm1)

    if m == 0 or n_cols == 0:
        return np.zeros((n_cols, m), dtype=float)

    mat = np.zeros((n_cols, m), dtype=float)

    for j, path in enumerate(paths_p):
        for i in range(1, p):  # Delete interior vertices (index 1 to p-1)
            # Check if deleting vertex at position i gives an allowed path
            new_path = path[:i] + path[i+1:]
            # Need to check: is new_path allowed? i.e., edges exist between consecutive vertices
            allowed = True
            for k in range(len(new_path)-1):
                if not A[new_path[k]][new_path[k+1]]:
                    allowed = False
                    break
            if allowed and new_path in path_to_idx:
                mat[path_to_idx[new_path], j] += (-1)**i

    return mat

def shift_action_on_paths(paths, n):
    """Compute the Z/nZ shift action: σ(v_0,...,v_p) = (v_0+1,...,v_p+1) mod n.
    Returns a permutation matrix."""
    path_to_idx = {p: i for i, p in enumerate(paths)}
    m = len(paths)
    perm = np.zeros((m, m), dtype=float)

    for j, path in enumerate(paths):
        shifted = tuple((v + 1) % n for v in path)
        if shifted in path_to_idx:
            perm[path_to_idx[shifted], j] = 1.0

    return perm

# ===== Main computation =====
print("=" * 70)
print("EIGENSPACE DECOMPOSITION OF n=9 CIRCULANT MAXIMIZER")
print("=" * 70)

n = 9
S = {1, 5, 6, 7}  # Known maximizer connection set
A = build_circulant_tournament(n, S)
H = H_tournament(A, n)
print(f"\nCirculant Z_{n}, S={sorted(S)}, H={H}")

# Verify it's a tournament and regular
for i in range(n):
    for j in range(i+1, n):
        assert A[i][j] + A[j][i] == 1, f"Not a tournament at ({i},{j})"
scores = [sum(A[i]) for i in range(n)]
print(f"Scores: {scores}")
assert all(s == (n-1)//2 for s in scores), "Not regular"

# Full Betti numbers for verification
t0 = time.time()
beta_full = path_betti_numbers(A, n, max_dim=7)
beta_list = [int(beta_full[k]) if k < len(beta_full) else 0 for k in range(8)]
print(f"Full Betti: {beta_list} ({time.time()-t0:.1f}s)")

# ===== Eigenspace decomposition =====
print("\n" + "=" * 70)
print("CHAIN COMPLEX STRUCTURE")
print("=" * 70)

# Build chain complex up to dimension 7
max_dim = 7
all_paths = {}
for d in range(max_dim + 1):
    t1 = time.time()
    paths = get_allowed_paths(A, n, d)
    all_paths[d] = paths
    print(f"  dim {d}: {len(paths)} allowed paths ({time.time()-t1:.1f}s)")

# Build boundary matrices
boundaries = {}
for d in range(1, max_dim + 1):
    if len(all_paths[d]) > 0 and len(all_paths[d-1]) > 0:
        t1 = time.time()
        B = boundary_matrix(A, n, all_paths[d], all_paths[d-1], d)
        boundaries[d] = B
        print(f"  ∂_{d}: {B.shape[0]}×{B.shape[1]}, rank={np.linalg.matrix_rank(B)} ({time.time()-t1:.1f}s)")
    else:
        boundaries[d] = np.zeros((len(all_paths[d-1]) if d-1 in all_paths else 0,
                                   len(all_paths[d])), dtype=float)
        print(f"  ∂_{d}: trivial")

# ===== Shift action on each chain group =====
print("\n" + "=" * 70)
print("Z/9Z EIGENSPACE DECOMPOSITION")
print("=" * 70)

omega = np.exp(2j * np.pi / n)  # Primitive 9th root of unity

for d in range(max_dim + 1):
    if len(all_paths[d]) == 0:
        print(f"\n  dim {d}: empty")
        continue

    t1 = time.time()
    sigma = shift_action_on_paths(all_paths[d], n)

    # Check that sigma is a permutation and sigma^9 = I
    sigma_9 = np.linalg.matrix_power(sigma.astype(int), n)
    assert np.allclose(sigma_9, np.eye(len(all_paths[d]))), f"σ^{n} ≠ I at dim {d}"

    # Eigenvalues of sigma
    evals = np.linalg.eigvals(sigma)

    # Count eigenvalues by 9th root of unity
    eigenspace_dims = {}
    for k in range(n):
        target = omega**k
        count = sum(1 for e in evals if abs(e - target) < 0.01)
        if count > 0:
            eigenspace_dims[k] = count

    print(f"\n  dim {d}: total={len(all_paths[d])}, eigenspace dims: {eigenspace_dims} ({time.time()-t1:.1f}s)")

# ===== Per-eigenspace Betti numbers =====
print("\n" + "=" * 70)
print("PER-EIGENSPACE BETTI NUMBERS (key: dim 5)")
print("=" * 70)

# Focus on the chain complex around dimension 5: ∂₆: Ω₆→Ω₅ and ∂₅: Ω₅→Ω₄
# β₅ = dim(ker ∂₅) - dim(im ∂₆)

# For each eigenspace k, project the chain complex
for k in range(n):
    # Build projector onto eigenspace k for each chain group
    projectors = {}
    projected_dims = {}

    for d in range(max_dim + 1):
        if len(all_paths[d]) == 0:
            projectors[d] = None
            projected_dims[d] = 0
            continue

        sigma = shift_action_on_paths(all_paths[d], n)
        m = len(all_paths[d])

        # Projector: P_k = (1/n) * sum_{j=0}^{n-1} omega^{-kj} sigma^j
        P = np.zeros((m, m), dtype=complex)
        for j in range(n):
            P += (omega**(-k*j)) * np.linalg.matrix_power(sigma, j)
        P /= n

        projectors[d] = P
        projected_dims[d] = round(np.real(np.trace(P)))

    # For eigenspace k, compute β_d for each d
    betti_k = []
    for d in range(max_dim + 1):
        if projected_dims[d] == 0:
            betti_k.append(0)
            continue

        # ker(∂_d) in eigenspace k
        if d == 0 or d not in boundaries or boundaries[d].shape[1] == 0:
            ker_dim = projected_dims[d]
        else:
            # Project boundary matrix: P_{d-1} @ ∂_d @ P_d
            B = boundaries[d].astype(complex)
            P_d = projectors[d]
            P_dm1 = projectors[d-1] if projectors[d-1] is not None else np.eye(B.shape[0], dtype=complex)
            B_proj = P_dm1 @ B @ P_d
            rank_B = np.linalg.matrix_rank(B_proj, tol=1e-6)
            ker_dim = projected_dims[d] - rank_B

        # im(∂_{d+1}) in eigenspace k
        if d+1 > max_dim or d+1 not in boundaries or boundaries[d+1].shape[1] == 0:
            im_dim = 0
        else:
            B_next = boundaries[d+1].astype(complex)
            P_dp1 = projectors[d+1] if d+1 in projectors and projectors[d+1] is not None else np.eye(B_next.shape[1], dtype=complex)
            P_d_cur = projectors[d] if projectors[d] is not None else np.eye(B_next.shape[0], dtype=complex)
            B_next_proj = P_d_cur @ B_next @ P_dp1
            im_dim = np.linalg.matrix_rank(B_next_proj, tol=1e-6)

        betti_k.append(max(0, int(round(ker_dim - im_dim))))

    if any(b > 0 for b in betti_k):
        print(f"  k={k}: Ω_dims={list(projected_dims.values())}, β={betti_k}")

# ===== Summary =====
print("\n" + "=" * 70)
print("BETTI FORMULA CONJECTURE TEST")
print("=" * 70)

print("\nOdd-n maximizer Betti data:")
data = [(3, 1, 1), (5, 1, 1), (7, 4, 6), (9, 5, 10)]
for nn, d, b in data:
    predicted = d*(d-1)//2  # C(d,2)
    print(f"  n={nn}: β_{d}={b}, C({d},2)={predicted}, match={'✓' if b==predicted else '✗'}")

print(f"\nConjectured formula for n≥7: β_{{(n+1)/2}} = C((n+1)/2, 2) = (n+1)(n-1)/8")
for nn in [11, 13, 15]:
    d = (nn+1)//2
    b = d*(d-1)//2
    print(f"  n={nn}: predict β_{d} = C({d},2) = {b}")

print(f"\nOpus conjecture (Paley): β_{{p-3}} = p-1")
for p in [7, 11, 13]:
    d = p-3
    b = p-1
    print(f"  p={p}: predict β_{d} = {b}")
print(f"\n  At n=11 (Paley): this gives β_8=10, our formula gives β_6=15")
print(f"  These are INCOMPATIBLE — n=11 would discriminate")

print(f"\nTotal time: {time.time()-t0:.1f}s")
