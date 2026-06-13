#!/usr/bin/env python3
"""
eigenspace_boundary_structure.py — opus-2026-03-12-S69

Investigate WHY Q_k determines the topology.

For a circulant tournament C_p^S, eigenspace k has:
  d_m: Ω_m(k) → Ω_{m-1}(k)

The boundary map in eigenspace k should depend on ω^k and S_hat(k).

Key insight: For circulant graphs, the adjacency matrix A has eigenvalue
  λ_k = S_hat(k) = Σ_{s∈S} ω^{sk}
with eigenvector v_k = (1, ω^k, ω^{2k}, ..., ω^{(p-1)k}).

Q_k = |λ_k|² = λ_k * conj(λ_k).

Question: Does d_m in eigenspace k depend on λ_k or just Q_k?
If on λ_k, then Q_k alone shouldn't determine Ω. But THM-145 says Q determines Ω...

This means either:
(a) d_m depends only on |λ_k|, not arg(λ_k), OR
(b) The ESFs e_j(Q) somehow capture the relevant information about the phases too.

Let's investigate by looking at the actual eigenspace boundary matrices.
"""

import sys
import numpy as np
from numpy.linalg import matrix_rank
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def compute_lambda(p, S, k):
    """Compute eigenvalue λ_k = S_hat(k)."""
    omega = np.exp(2j * np.pi / p)
    return sum(omega**(s*k) for s in S)

# Start with p=7
p = 7
m = (p-1)//2

# Two orbits at p=7
orientations = [
    frozenset({1, 2, 3}),  # Interval
    frozenset({1, 2, 5}),  # Non-interval (Paley at p=7 ≡ 3 mod 4)
]

print(f"{'='*80}")
print(f"EIGENSPACE BOUNDARY STRUCTURE at p = {p}")
print(f"{'='*80}")

for S in orientations:
    print(f"\n--- S = {sorted(S)} ---")

    # Compute all eigenvalues
    print(f"\n  Eigenvalues λ_k = S_hat(k):")
    for k in range(1, m+1):
        lam = compute_lambda(p, S, k)
        Q = abs(lam)**2
        arg = np.angle(lam)
        print(f"    k={k}: λ = {lam:.4f} = {abs(lam):.4f} * e^(i*{arg:.4f}), Q = {Q:.4f}")

    # Key: λ_k and λ_{p-k} are conjugate (since S_hat is DFT of real-like data)
    print(f"\n  Conjugacy check:")
    for k in range(1, m+1):
        lam_k = compute_lambda(p, S, k)
        lam_pk = compute_lambda(p, S, p-k)
        print(f"    λ_{k} = {lam_k:.4f}, λ_{p-k} = {lam_pk:.4f}, product = {lam_k*lam_pk:.4f}")

# The key realization: for eigenspace decomposition of the BOUNDARY MAP,
# what matters is how the boundary operator acts on ω^k-eigenvectors.
# Let's look at this more carefully.

print(f"\n{'='*80}")
print(f"BOUNDARY MAP IN EIGENSPACE k")
print(f"{'='*80}")

# For an allowed m-path (v_0, v_1, ..., v_m) in T, the boundary is:
# d_m(path) = Σ (-1)^i (path with v_i removed)
#
# In eigenspace k, an m-path starting at vertex 0 contributes:
# ω^{k*0} * (contribution from remaining vertices)
# After projecting to eigenspace k, each path (a_0,...,a_m) becomes
# ω^{k*a_0} * path-weight.

# Actually, the key insight is that for a CIRCULANT tournament,
# the boundary map d_m: Ω_m → Ω_{m-1} commutes with the cyclic group action.
# So in each eigenspace k, d_m restricts to a linear map
# d_m(k): Ω_m(k) → Ω_{m-1}(k)
# whose rank determines the Betti numbers.

# Let's verify this by explicitly computing d_m(k) for small cases.

p = 7
m = 3
S = frozenset({1, 2, 3})  # Interval

print(f"\np={p}, S={sorted(S)}")
h = CirculantHomology(n=p, S=S)
h._ensure_enumerated(p)

# Get Ω_m for each m
for mm in range(p):
    omega_dim = h.omega_dims(max_degree=mm)
    if mm < len(omega_dim):
        print(f"  Ω_{mm} = {omega_dim[mm]}")

# For eigenspace analysis, we need the full boundary matrices
# Let's compute them and project to eigenspace k

# Actually, let me approach this differently.
# The CIRCULANT structure means we can work with DIFFERENCE SEQUENCES.
# An allowed m-path is (v_0, v_1, ..., v_m) where v_i → v_{i+1} in T.
# For circulant T, this means (v_{i+1} - v_i) mod p ∈ S.
# By translation invariance, the path starting at 0 is determined by
# the DIFFERENCE SEQUENCE (d_1, ..., d_m) where d_i = v_i - v_{i-1} ∈ S.

# BUT we also need all vertices DISTINCT. This is the non-trivial constraint.
# In fact, the ALLOWED paths starting at 0 are exactly the sequences
# (d_1, ..., d_m) with d_i ∈ S and all partial sums 0, d_1, d_1+d_2, ... distinct mod p.

# In eigenspace k, such a path (0, s_1, s_1+s_2, ...) contributes
# a weight ω^{k*(s_1 + (s_1+s_2) + ... + (s_1+...+s_m))} ??? No.
# Actually, the eigenvector projection maps vertex v to ω^{kv}.
# The path (v_0, ..., v_m) in eigenspace k has coefficient ω^{k*v_0}.
# Since we quotient by cyclic group, paths starting at different v_0
# give the same path in the quotient, so eigenspace k projects the path
# (v_0, ..., v_m) to ω^{k*v_0} times the "shape" of the path.

# For the BOUNDARY MAP in eigenspace k:
# d_m(path) = Σ (-1)^i (path with v_i removed)
# Each resulting (m-1)-path, projected to eigenspace k, gets weight ω^{k*v_0}.
# So d_m(k) maps m-path shapes to (m-1)-path shapes.

# The MATRIX of d_m(k) has entries involving ω^k, and its rank determines β_m(k).

# Let me compute this explicitly for p=7.

# Step 1: enumerate all path "shapes" (difference sequences)
def allowed_paths(p, S, length):
    """Enumerate allowed path shapes (difference sequences) of given length."""
    S_list = sorted(S)
    results = []

    def backtrack(diffs, partial_sums):
        if len(diffs) == length:
            results.append(tuple(diffs))
            return
        for s in S_list:
            new_sum = (partial_sums[-1] + s) % p
            if new_sum not in set(partial_sums):
                backtrack(diffs + [s], partial_sums + [new_sum])

    backtrack([], [0])
    return results

print(f"\n{'='*80}")
print(f"PATH SHAPES (difference sequences) for p={p}, S={sorted(S)}")
print(f"{'='*80}")

for length in range(1, p):
    paths = allowed_paths(p, S, length)
    if paths:
        print(f"  {length}-paths: {len(paths)} (should be Ω_{length}/p... wait, Ω_{length} total)")
        # Actually Ω_m is the TOTAL number of allowed (m+1)-tuples (including all starting vertices)
        # Per eigenspace: Ω_m / p paths per eigenspace
        # But the number of SHAPES is Ω_m / p (by translation invariance)
        # Actually for allowed paths, we need p * (number of shapes) since any starting vertex works
        # No — some shapes might have repeated partial sums for certain starting points.
        # For p prime and circulant, by shift: every shape starting at 0 is valid,
        # and there are p translates.
        # So #shapes = Ω_{length} / p if there's no overcounting.
        # But Ω_m counts (m+1)-tuples where ALL vertices are in {0,...,p-1} and distinct.
        # Hmm, actually Ω_m from the enumeration includes all starting vertices.

print(f"\n  Verification: path shapes * p = total paths?")
omega = h.omega_dims(max_degree=p-1)
for length in range(1, p):
    paths = allowed_paths(p, S, length)
    expected = omega[length] if length < len(omega) else 0
    actual = len(paths) * p
    print(f"    length={length}: shapes={len(paths)}, shapes*p={actual}, Ω_{length}={expected} {'✓' if actual == expected else '✗'}")

# Step 2: Boundary map in eigenspace k
# For eigenspace k, shape (d_1,...,d_m) maps to:
# d_m(shape) = Σ_{i=0}^{m} (-1)^i * (shape with vertex i removed)
# But removing vertex i from (0, s_1, s_1+s_2, ...) gives a different shape
# AND a phase factor.

# Actually, let me think about this more carefully.
# The boundary of (v_0, v_1, ..., v_m) is Σ_i (-1)^i (v_0,...,v_{i-1},v_{i+1},...,v_m)
# But (v_0,...,v_{i-1},v_{i+1},...,v_m) is a valid (m-1)-path only if
# v_{i-1} → v_{i+1} in T (i.e., (v_{i+1}-v_{i-1}) mod p ∈ S) for i=1,...,m-1.
# For i=0 or i=m, the resulting path is always valid (just shorter).

# In eigenspace k, the path (v_0,...,v_m) corresponds to shape (d_1,...,d_m)
# with phase ω^{k*v_0}. The boundary term removing v_i:
# - For i=0: path (v_1,...,v_m) with shape (d_2,...,d_m) and phase ω^{k*v_1} = ω^{k*(v_0+d_1)}
#   Relative to starting vertex v_0, this contributes ω^{k*d_1} * shape(d_2,...,d_m)
# - For i=m: path (v_0,...,v_{m-1}) with shape (d_1,...,d_{m-1}), same phase
# - For 0<i<m: shape (d_1,...,d_{i-1}, d_i+d_{i+1}, d_{i+2},...,d_m)
#   Valid only if d_i+d_{i+1} ∈ S ∪ {0,...,p-1}\{0} (actually, just that the resulting
#   vertex sequence is valid, which requires d_i+d_{i+1} mod p to keep vertices distinct)

# This is getting complex. Let me just compute the boundary MATRIX numerically.
omega_k = np.exp(2j * np.pi / p)

print(f"\n{'='*80}")
print(f"EIGENSPACE BOUNDARY MATRICES for p={p}, S={sorted(S)}")
print(f"{'='*80}")

for k in range(p):
    lam_k = compute_lambda(p, S, k)
    Q_k = abs(lam_k)**2

    # For each degree m, compute boundary matrix restricted to eigenspace k
    # Basis: path shapes (d_1,...,d_m)
    max_deg = min(6, p-1)

    ranks = []
    for mm in range(1, max_deg+1):
        paths_m = allowed_paths(p, S, mm)
        paths_m1 = allowed_paths(p, S, mm-1) if mm > 1 else [()]

        if not paths_m or not paths_m1:
            ranks.append(0)
            continue

        # Build boundary matrix: d_m(k)[i,j] = coefficient of j-th (m-1)-path
        # in boundary of i-th m-path, projected to eigenspace k

        path_m1_idx = {p: i for i, p in enumerate(paths_m1)}
        if mm == 1:
            # 0-paths are just single vertices, represented by ()
            path_m1_idx = {(): 0}

        n_source = len(paths_m)
        n_target = len(paths_m1)
        D = np.zeros((n_target, n_source), dtype=complex)

        for j, shape in enumerate(paths_m):
            # Boundary: remove each vertex
            for i in range(mm + 1):
                sign = (-1)**i

                if i == 0:
                    # Remove first vertex: shape becomes (d_2,...,d_m)
                    # Phase shift: ω^{k*d_1}
                    new_shape = shape[1:]
                    phase = omega_k**(k * shape[0])
                    if new_shape in path_m1_idx:
                        D[path_m1_idx[new_shape], j] += sign * phase
                elif i == mm:
                    # Remove last vertex: shape becomes (d_1,...,d_{m-1})
                    new_shape = shape[:mm-1]
                    phase = 1.0
                    if new_shape in path_m1_idx:
                        D[path_m1_idx[new_shape], j] += sign * phase
                else:
                    # Remove interior vertex i: merge d_i and d_{i+1}
                    merged = (shape[i-1] + shape[i]) % p
                    new_shape = shape[:i-1] + (merged,) + shape[i+1:]
                    phase = 1.0
                    # Check if new_shape is a valid path
                    if new_shape in path_m1_idx:
                        D[path_m1_idx[new_shape], j] += sign * phase

        r = matrix_rank(D, tol=1e-8)
        ranks.append(r)

    print(f"  k={k}: Q={Q_k:.4f}, boundary ranks = {ranks}")

print("\nDONE.")
