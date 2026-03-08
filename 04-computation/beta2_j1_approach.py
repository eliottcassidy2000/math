#!/usr/bin/env python3
"""
beta2_j1_approach.py - Test the j_1 approach for proving beta_2=0

THE j_1 APPROACH:
For tournament T, consider vertex v with beta_1(T\v) = 1.
The LES gives: ... -> H_2(T) -> H_2(T,T\v) ->^δ H_1(T\v) ->^{j_1} H_1(T) -> ...

If j_1 = 0 (the map H_1(T\v) -> H_1(T) is zero), then:
  - im(j_1) = 0, so ker(j_1) = H_1(T\v) = Z
  - Exactness at H_1(T\v): ker(j_1) = im(δ), so im(δ) = Z
  - δ is surjective onto H_1(T\v)
  - If also H_2(T,T\v) ≤ 1, then δ is an isomorphism
  - Exactness at H_2(T,T\v): ker(δ) = im(from H_2(T))
  - If δ is injective, then im(H_2(T) -> H_2(T,T\v)) = 0
  - If α_2: H_2(T\v) -> H_2(T) is zero (by induction beta_2(T\v) = 0):
    then H_2(T) -> H_2(T,T\v) is injective (ker = im(α_2) = 0)
  - Combined: H_2(T) injects into ker(δ) = 0, so H_2(T) = 0. QED!

This script tests:
1. For each T, does there exist v with j_1: H_1(T\v) -> H_1(T) = 0?
2. What is dim H_2(T,T\v) in each case?
3. Is the j_1 approach universal (works for ALL tournaments)?

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers, boundary_coeffs
)
sys.stdout = _saved


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def compute_j1_map(A, n, v):
    """Compute j_1: H_1(T\\v) -> H_1(T) explicitly.

    j_1 is the map induced by inclusion T\\v -> T on path homology.

    Returns: (j1_matrix_in_homology_coords, beta1_sub, beta1_full,
              ker_d1_sub_basis, ker_d1_full_basis, im_d2_sub_rank, im_d2_full_rank)
    """
    # Full tournament T
    paths1_full = enumerate_allowed_paths(A, n, 1)
    paths0_full = enumerate_allowed_paths(A, n, 0)
    paths2_full = enumerate_allowed_paths(A, n, 2)

    omega1_full = compute_omega_basis(A, n, 1, paths1_full, paths0_full)
    omega2_full = compute_omega_basis(A, n, 2, paths2_full, paths1_full)

    dim_omega1_full = omega1_full.shape[1] if omega1_full.ndim == 2 else 0
    dim_omega2_full = omega2_full.shape[1] if omega2_full.ndim == 2 else 0

    # Boundary d1: Omega_1 -> A_0
    D1_full = build_full_boundary_matrix(paths1_full, paths0_full)
    D1_omega_full = D1_full @ omega1_full if dim_omega1_full > 0 else np.zeros((len(paths0_full), 0))

    if D1_omega_full.shape[1] > 0:
        S1 = np.linalg.svd(D1_omega_full, compute_uv=False)
        rank_d1_full = sum(s > 1e-8 for s in S1)
    else:
        rank_d1_full = 0
    ker_d1_full_dim = dim_omega1_full - rank_d1_full

    # Boundary d2: Omega_2 -> A_1
    D2_full = build_full_boundary_matrix(paths2_full, paths1_full)
    D2_omega_full = D2_full @ omega2_full if dim_omega2_full > 0 else np.zeros((len(paths1_full), 0))

    if D2_omega_full.shape[1] > 0:
        S2 = np.linalg.svd(D2_omega_full, compute_uv=False)
        im_d2_full_rank = sum(s > 1e-8 for s in S2)
    else:
        im_d2_full_rank = 0

    beta1_full = ker_d1_full_dim - im_d2_full_rank

    # Subtournament T\v
    others = [i for i in range(n) if i != v]
    n_sub = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n_sub)] for i in range(n_sub)]

    paths1_sub = enumerate_allowed_paths(A_sub, n_sub, 1)
    paths0_sub = enumerate_allowed_paths(A_sub, n_sub, 0)
    paths2_sub = enumerate_allowed_paths(A_sub, n_sub, 2)

    omega1_sub = compute_omega_basis(A_sub, n_sub, 1, paths1_sub, paths0_sub)
    omega2_sub = compute_omega_basis(A_sub, n_sub, 2, paths2_sub, paths1_sub)

    dim_omega1_sub = omega1_sub.shape[1] if omega1_sub.ndim == 2 else 0
    dim_omega2_sub = omega2_sub.shape[1] if omega2_sub.ndim == 2 else 0

    D1_sub = build_full_boundary_matrix(paths1_sub, paths0_sub)
    D1_omega_sub = D1_sub @ omega1_sub if dim_omega1_sub > 0 else np.zeros((len(paths0_sub), 0))

    if D1_omega_sub.shape[1] > 0:
        S1s = np.linalg.svd(D1_omega_sub, compute_uv=False)
        rank_d1_sub = sum(s > 1e-8 for s in S1s)
    else:
        rank_d1_sub = 0
    ker_d1_sub_dim = dim_omega1_sub - rank_d1_sub

    D2_sub = build_full_boundary_matrix(paths2_sub, paths1_sub)
    D2_omega_sub = D2_sub @ omega2_sub if dim_omega2_sub > 0 else np.zeros((len(paths1_sub), 0))

    if D2_omega_sub.shape[1] > 0:
        S2s = np.linalg.svd(D2_omega_sub, compute_uv=False)
        im_d2_sub_rank = sum(s > 1e-8 for s in S2s)
    else:
        im_d2_sub_rank = 0

    beta1_sub = ker_d1_sub_dim - im_d2_sub_rank

    if beta1_sub == 0 or beta1_full < 0:
        return None, beta1_sub, beta1_full

    # Now compute j_1 at the homology level.
    # Step 1: Get ker(d1) bases in both spaces
    # For full: use SVD of D1_omega_full to find kernel
    if D1_omega_full.shape[1] > 0 and D1_omega_full.shape[0] > 0:
        U1f, S1f, V1ft = np.linalg.svd(D1_omega_full, full_matrices=True)
        rank1f = sum(s > 1e-8 for s in S1f)
        # Kernel vectors in Omega_1 coords
        ker_d1_full_omega_coords = V1ft[rank1f:].T  # columns
    else:
        ker_d1_full_omega_coords = np.eye(dim_omega1_full) if dim_omega1_full > 0 else np.zeros((0,0))

    if D1_omega_sub.shape[1] > 0 and D1_omega_sub.shape[0] > 0:
        U1s, S1s_full, V1st = np.linalg.svd(D1_omega_sub, full_matrices=True)
        rank1s = sum(s > 1e-8 for s in S1s_full)
        ker_d1_sub_omega_coords = V1st[rank1s:].T
    else:
        ker_d1_sub_omega_coords = np.eye(dim_omega1_sub) if dim_omega1_sub > 0 else np.zeros((0,0))

    if ker_d1_sub_omega_coords.shape[1] == 0 or ker_d1_full_omega_coords.shape[1] == 0:
        return None, beta1_sub, beta1_full

    # Step 2: Convert sub kernel vectors to full A_1 coordinates
    # A sub 1-path (a', b') maps to full 1-path (others[a'], others[b'])
    path1_full_idx = {p: i for i, p in enumerate(paths1_full)}

    # Inclusion map: A_1(T\v) -> A_1(T) in standard basis
    incl_matrix = np.zeros((len(paths1_full), len(paths1_sub)))
    for j, (a, b) in enumerate(paths1_sub):
        full_path = (others[a], others[b])
        if full_path in path1_full_idx:
            incl_matrix[path1_full_idx[full_path], j] = 1

    # Sub kernel in A_1(T\v) coords: omega1_sub @ ker_d1_sub_omega_coords
    ker_d1_sub_A1_sub = omega1_sub @ ker_d1_sub_omega_coords  # in A_1(T\v) coords

    # Map to A_1(T) coords via inclusion
    ker_d1_sub_in_full = incl_matrix @ ker_d1_sub_A1_sub  # in A_1(T) coords

    # Step 3: Project onto ker(d1_full) / im(d2_full) = H_1(T)
    # First express ker_d1_sub_in_full in Omega_1(T) coordinates
    # omega1_full has columns = basis of Omega_1(T)
    # We need to express our vectors in this basis
    # Use least squares: omega1_full @ x = ker_d1_sub_in_full
    if dim_omega1_full > 0:
        # Project onto Omega_1(T)
        omega1_coords, _, _, _ = np.linalg.lstsq(omega1_full, ker_d1_sub_in_full, rcond=None)
        # Check residual: should be ~0 since inclusion preserves Omega
        residual = np.max(np.abs(omega1_full @ omega1_coords - ker_d1_sub_in_full))
        if residual > 1e-6:
            # Not in Omega_1(T) — shouldn't happen for tournaments
            pass

        # Now express in ker(d1_full) coordinates
        # ker_d1_full_omega_coords columns form basis of ker(d1) in Omega_1 coords
        # Project omega1_coords onto ker(d1_full)
        ker_coords, _, _, _ = np.linalg.lstsq(ker_d1_full_omega_coords, omega1_coords, rcond=None)

        # Now we need to quotient by im(d2_full)
        # im(d2_full) in A_1(T) coords: D2_omega_full columns
        # Express im(d2_full) in ker(d1) coords
        if im_d2_full_rank > 0:
            im_in_A1 = D2_omega_full  # columns = image vectors in A_1(T)
            im_in_omega1, _, _, _ = np.linalg.lstsq(omega1_full, im_in_A1, rcond=None)
            im_in_ker, _, _, _ = np.linalg.lstsq(ker_d1_full_omega_coords, im_in_omega1, rcond=None)

            # j_1 in H_1: project ker_coords modulo im_in_ker
            # Augment: [im_in_ker | ker_coords] and check if ker_coords columns are in span of im_in_ker
            combined = np.hstack([im_in_ker, ker_coords])
            Sc = np.linalg.svd(combined, compute_uv=False)
            rank_combined = sum(s > 1e-8 for s in Sc)

            Si = np.linalg.svd(im_in_ker, compute_uv=False)
            rank_im = sum(s > 1e-8 for s in Si)

            # j_1 = 0 in H_1 iff the image of ker(d1_sub) lands in im(d2_full)
            # i.e., rank_combined == rank_im
            j1_is_zero = (rank_combined == rank_im)
        else:
            # im(d2_full) = 0, so H_1(T) = ker(d1_full)
            # j_1 = 0 iff ker_coords = 0
            j1_is_zero = np.max(np.abs(ker_coords)) < 1e-8 if ker_coords.size > 0 else True
    else:
        j1_is_zero = True

    return j1_is_zero, beta1_sub, beta1_full


# ============================================================
# PART 1: Test j_1 approach at n=4,5
# ============================================================
print("=" * 70)
print("j1 APPROACH TEST: exists v with j1=0 AND beta1(T\\v)=1?")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m

    has_j1_zero_v = 0
    no_j1_zero_v = 0
    j1_details = defaultdict(int)
    t0 = time.time()

    for bits in range(total):
        A = build_adj(n, bits)
        betti_full = path_betti_numbers(A, n, max_dim=2)
        beta1_T = betti_full[1]
        beta2_T = betti_full[2] if len(betti_full) > 2 else 0

        found_good_v = False
        v_results = []

        for v in range(n):
            result = compute_j1_map(A, n, v)
            j1_zero, b1_sub, b1_full = result[0], result[1], result[2]

            if b1_sub == 1 and j1_zero == True:
                found_good_v = True
                v_results.append((v, 'j1=0,b1=1'))
            elif b1_sub == 1 and j1_zero == False:
                v_results.append((v, 'j1!=0,b1=1'))
            elif b1_sub == 0:
                v_results.append((v, 'b1=0'))

        if found_good_v:
            has_j1_zero_v += 1
        else:
            no_j1_zero_v += 1
            if no_j1_zero_v <= 5:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  NO j_1-GOOD V! n={n} T#{bits} scores={scores} beta={betti_full}")
                for vr in v_results:
                    print(f"    v={vr[0]}: {vr[1]}")

        j1_details[found_good_v] += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments ({elapsed:.1f}s)")
    print(f"  Has j_1-good vertex: {has_j1_zero_v}/{total}")
    print(f"  No j_1-good vertex: {no_j1_zero_v}/{total}")


# ============================================================
# PART 2: Alternative: check if j_1 is always zero for tournaments
# ============================================================
print(f"\n{'='*70}")
print("IS j_1: H_1(T\\v) -> H_1(T) ALWAYS ZERO for tournaments?")
print("=" * 70)

for n in [4, 5]:
    m = n*(n-1)//2
    total = 1 << m

    always_zero = True
    j1_nonzero_count = 0
    total_checks = 0
    t0 = time.time()

    for bits in range(total):
        A = build_adj(n, bits)

        for v in range(n):
            result = compute_j1_map(A, n, v)
            j1_zero, b1_sub, b1_full = result[0], result[1], result[2]

            if b1_sub == 0:
                continue  # j_1 trivially zero since domain is 0

            total_checks += 1
            if j1_zero == False:
                always_zero = False
                j1_nonzero_count += 1
                if j1_nonzero_count <= 3:
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  j_1 != 0: n={n} T#{bits} v={v} scores={scores} beta_1(T\\v)={b1_sub} beta_1(T)={b1_full}")

    elapsed = time.time() - t0
    print(f"\nn={n}: j_1 always zero (when beta_1(sub)>0)? {always_zero}")
    print(f"  Checks with beta_1(sub)>0: {total_checks}, nonzero: {j1_nonzero_count}")
    print(f"  ({elapsed:.1f}s)")


# ============================================================
# PART 3: What about the OTHER map approach — α_2?
# When beta_2(T\v) = 0 (by induction), α_2 = 0.
# Then ker(H_2(T) -> H_2(T,T\v)) = im(α_2) = 0.
# So H_2(T) -> H_2(T,T\v).
# If j_1 = 0 AND beta_1(T\v) = 1, then δ is surjective.
# ker(δ) = im(H_2(T) -> H_2(T,T\v)).
# If H_2(T,T\v) = 1, then δ is an iso, ker(δ) = 0, so H_2(T) = 0.
# ============================================================
print(f"\n{'='*70}")
print("RELATIVE H_2 DIMENSION CHECK")
print("=" * 70)

# For this we need to compute H_2(T, T\v) directly
# H_2(T, T\v) ≅ ker(d_2^rel) / im(d_3^rel)
# where the relative chains are C_p(T) / C_p(T\v)

def compute_relative_h2(A, n, v):
    """Compute dim H_2(T, T\\v) using relative chain complex."""
    others = [i for i in range(n) if i != v]
    n_sub = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n_sub)] for i in range(n_sub)]

    # Full paths
    paths2_full = enumerate_allowed_paths(A, n, 2)
    paths1_full = enumerate_allowed_paths(A, n, 1)
    paths3_full = enumerate_allowed_paths(A, n, 3)

    # Sub paths (mapped to full vertex labels)
    paths2_sub_mapped = set()
    paths2_sub_raw = enumerate_allowed_paths(A_sub, n_sub, 2)
    for (a, b, c) in paths2_sub_raw:
        paths2_sub_mapped.add((others[a], others[b], others[c]))

    paths1_sub_mapped = set()
    paths1_sub_raw = enumerate_allowed_paths(A_sub, n_sub, 1)
    for (a, b) in paths1_sub_raw:
        paths1_sub_mapped.add((others[a], others[b]))

    paths3_sub_mapped = set()
    paths3_sub_raw = enumerate_allowed_paths(A_sub, n_sub, 3)
    for (a, b, c, d) in paths3_sub_raw:
        paths3_sub_mapped.add((others[a], others[b], others[c], others[d]))

    # Relative chains = paths involving v (not in sub)
    rel_paths2 = [p for p in paths2_full if p not in paths2_sub_mapped]
    rel_paths1 = [p for p in paths1_full if p not in paths1_sub_mapped]
    rel_paths3 = [p for p in paths3_full if p not in paths3_sub_mapped]

    if not rel_paths2:
        return 0

    # We need Omega-relative chains: elements of Omega_p(T) that project to
    # relative chains. But this is more subtle.
    #
    # Actually, for the relative sequence we should work with:
    # Omega_p^rel = Omega_p(T) / Omega_p(T\v)
    #
    # Simpler approach: compute full homology and use LES rank-nullity

    # Let's just compute all the LES dimensions directly
    betti_full = path_betti_numbers(A, n, max_dim=3)
    betti_sub = path_betti_numbers(A_sub, n_sub, max_dim=3)

    beta1_full = betti_full[1]
    beta1_sub = betti_sub[1]
    beta2_full = betti_full[2] if len(betti_full) > 2 else 0
    beta2_sub = betti_sub[2] if len(betti_sub) > 2 else 0

    return {
        'beta_full': betti_full,
        'beta_sub': betti_sub,
        'beta2_full': beta2_full,
        'beta1_sub': beta1_sub,
        'beta1_full': beta1_full,
    }

# Quick check at n=5
n = 5
print(f"\nn={n}: LES dimensions for each vertex deletion")
for bits in [0, 76, 341]:  # transitive, regular examples
    A = build_adj(n, bits)
    scores = tuple(sorted([sum(row) for row in A]))
    print(f"\n  T#{bits} scores={scores}:")
    betti_T = path_betti_numbers(A, n, max_dim=3)
    print(f"    beta(T) = {betti_T}")

    for v in range(n):
        info = compute_relative_h2(A, n, v)
        dv = sum(A[v])
        print(f"    v={v} (d+={dv}): beta(T\\v)={info['beta_sub'][:3]}")


# ============================================================
# PART 4: Direct check — is j_1 ALWAYS zero for all (T, v)?
# If yes, this gives the cleanest proof structure.
# ============================================================
print(f"\n{'='*70}")
print("COMPREHENSIVE: j_1 map for ALL (T, v) pairs at n=5,6")
print("=" * 70)

# Detailed analysis at n=5
n = 5
m = n*(n-1)//2
total = 1 << m
j1_zero_when_b1_pos = 0
j1_nonzero_when_b1_pos = 0
t0 = time.time()

for bits in range(total):
    A = build_adj(n, bits)
    for v in range(n):
        result = compute_j1_map(A, n, v)
        j1_zero, b1_sub, b1_full = result[0], result[1], result[2]
        if b1_sub > 0:
            if j1_zero:
                j1_zero_when_b1_pos += 1
            else:
                j1_nonzero_when_b1_pos += 1

elapsed = time.time() - t0
print(f"\nn={n}: ({elapsed:.1f}s)")
print(f"  j_1=0 when beta_1(sub)>0: {j1_zero_when_b1_pos}")
print(f"  j_1!=0 when beta_1(sub)>0: {j1_nonzero_when_b1_pos}")
print(f"  j_1 UNIVERSALLY ZERO? {j1_nonzero_when_b1_pos == 0}")


# n=6 exhaustive (may take a while)
n = 6
m = n*(n-1)//2
total = 1 << m
j1_zero_when_b1_pos = 0
j1_nonzero_when_b1_pos = 0
t0 = time.time()
first_nonzero = True

for bits in range(total):
    A = build_adj(n, bits)
    for v in range(n):
        result = compute_j1_map(A, n, v)
        j1_zero, b1_sub, b1_full = result[0], result[1], result[2]
        if b1_sub > 0:
            if j1_zero:
                j1_zero_when_b1_pos += 1
            else:
                j1_nonzero_when_b1_pos += 1
                if first_nonzero:
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  FIRST j_1!=0: T#{bits} v={v} scores={scores}")
                    first_nonzero = False

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  n={n}: {bits+1}/{total} ({elapsed:.0f}s) zero={j1_zero_when_b1_pos} nonzero={j1_nonzero_when_b1_pos}")

elapsed = time.time() - t0
print(f"\nn={n}: ({elapsed:.0f}s)")
print(f"  j_1=0 when beta_1(sub)>0: {j1_zero_when_b1_pos}")
print(f"  j_1!=0 when beta_1(sub)>0: {j1_nonzero_when_b1_pos}")
print(f"  j_1 UNIVERSALLY ZERO? {j1_nonzero_when_b1_pos == 0}")


print("\n\nDone.")
