#!/usr/bin/env python3
"""
beta2_homotopy_proof.py - Test a contracting homotopy for beta_2 = 0

IDEA: For a tournament T on n vertices, can we build a chain homotopy
s: Omega_2 -> Omega_3 such that d_3 * s = id on ker(d_2|Omega_2)?

If yes, then ker(d_2) = im(d_3), so beta_2 = 0.

APPROACH 1: Source cone
For transitive T with source v: s(a,b,c) = (v,a,b,c) works perfectly.
The cone from the source gives a contracting homotopy.

APPROACH 2: Averaged cone
For general T: s = (1/n) * sum_v s_v where s_v(a,b,c) = (v,a,b,c) if v->a.
Does d_3*s = id on ker(d_2)? We check this explicitly.

APPROACH 3: Signed cone
Maybe s_v needs to be adjusted: s_v(a,b,c) = epsilon_v * (v,a,b,c).
Can we find signs epsilon_v such that d_3*(sum epsilon_v s_v) = id?

This script tests these approaches at n=4,5.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, boundary_coeffs
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


def test_cone_homotopy(A, n, verbose=False):
    """Test whether the averaged cone homotopy gives d_3*s = id on ker(d_2).

    The cone from vertex v maps 2-path (a,b,c) to 3-path (v,a,b,c)
    if v->a and v not in {a,b,c}.

    d_3(v,a,b,c) = (a,b,c) - (v,b,c) + (v,a,c) - (v,a,b)

    For the averaged cone: s(sigma) = (1/K) * sum_{v: v->a, v not in sigma} (v,sigma)
    where K = |{v: v->a_0}| = number of in-neighbors of a_0 not in sigma.

    Actually, let's use the UNIFORM average: s = sum_v s_v / n.
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)
    paths0 = enumerate_allowed_paths(A, n, 0)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)

    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        return True, 0, "trivial"

    # Boundary d_2: Omega_2 -> A_1
    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2

    # Kernel of d_2 on Omega_2
    if D2_omega.shape[1] > 0 and D2_omega.shape[0] > 0:
        U, S, Vt = np.linalg.svd(D2_omega, full_matrices=True)
        rank_d2 = sum(s > 1e-8 for s in S)
        ker_d2_basis = Vt[rank_d2:].T  # columns = kernel vectors in Omega_2 coords
    else:
        ker_d2_basis = np.eye(dim_O2) if dim_O2 > 0 else np.zeros((0,0))

    ker_dim = ker_d2_basis.shape[1] if ker_d2_basis.ndim == 2 else 0

    if ker_dim == 0:
        return True, 0, "empty kernel"

    # Convert kernel vectors to A_2 coordinates
    ker_in_A2 = omega2 @ ker_d2_basis  # shape (|A_2|, ker_dim)

    # Build cone map: for each vertex v, s_v maps A_2 -> A_3
    # s_v(a,b,c) = (v,a,b,c) if v->a and v not in {a,b,c}, else 0
    path3_idx = {p: i for i, p in enumerate(paths3)}
    path2_idx = {p: i for i, p in enumerate(paths2)}

    # For each v, build cone matrix: |A_3| x |A_2|
    cone_matrices = {}
    for v in range(n):
        S_v = np.zeros((len(paths3), len(paths2)))
        for j, (a, b, c) in enumerate(paths2):
            if v in (a, b, c):
                continue
            if A[v][a] == 1:  # v -> a
                target = (v, a, b, c)
                if target in path3_idx:
                    S_v[path3_idx[target], j] = 1
        cone_matrices[v] = S_v

    # Averaged cone: S_avg = (1/n) * sum_v S_v
    S_avg = sum(cone_matrices[v] for v in range(n)) / n

    # Now compute d_3 * S_avg on ker(d_2)
    # d_3: A_3 -> A_2
    D3 = build_full_boundary_matrix(paths3, paths2)

    # d_3 * S_avg on ker vectors (in A_2 coords)
    # Result should be in A_2 coords
    result = D3 @ S_avg @ ker_in_A2  # shape (|A_2|, ker_dim)

    # Check: does result = ker_in_A2? (i.e., d_3 * s = id on ker(d_2))
    diff = result - ker_in_A2
    max_error = np.max(np.abs(diff)) if diff.size > 0 else 0

    if verbose and max_error > 1e-8:
        # What fraction of the identity does the cone recover?
        # Project result onto ker space
        proj = np.linalg.lstsq(ker_in_A2, result, rcond=None)[0]
        return False, max_error, proj

    return max_error < 1e-8, max_error, None


def test_weighted_cone(A, n):
    """Test whether WEIGHTED cone sum_v w_v * s_v gives d_3*s = id on ker(d_2).

    Solve for weights w_v such that d_3 * (sum w_v s_v) * z = z for all z in ker(d_2).
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        return True, None

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2

    U, S, Vt = np.linalg.svd(D2_omega, full_matrices=True)
    rank_d2 = sum(s > 1e-8 for s in S)
    ker_d2_basis = Vt[rank_d2:].T
    ker_dim = ker_d2_basis.shape[1] if ker_d2_basis.ndim == 2 else 0

    if ker_dim == 0:
        return True, None

    ker_in_A2 = omega2 @ ker_d2_basis

    path3_idx = {p: i for i, p in enumerate(paths3)}
    path2_idx = {p: i for i, p in enumerate(paths2)}
    D3 = build_full_boundary_matrix(paths3, paths2)

    # For each v, compute d_3 * s_v restricted to ker
    # d_3 * s_v * ker_in_A2 gives |A_2| x ker_dim matrix
    # We want sum_v w_v * (d_3 * s_v * ker_in_A2) = ker_in_A2

    # Stack: M_v = D3 @ S_v @ ker_in_A2 for each v
    # This is |A_2| x ker_dim for each v
    # We want sum_v w_v M_v = ker_in_A2
    # Vectorize: reshape M_v to (|A_2|*ker_dim) vector for each v
    # Then we have a linear system in n unknowns (the weights)

    M_vecs = []
    for v in range(n):
        S_v = np.zeros((len(paths3), len(paths2)))
        for j, (a, b, c) in enumerate(paths2):
            if v in (a, b, c):
                continue
            if A[v][a] == 1:
                target = (v, a, b, c)
                if target in path3_idx:
                    S_v[path3_idx[target], j] = 1
        M_v = D3 @ S_v @ ker_in_A2
        M_vecs.append(M_v.flatten())

    # target = ker_in_A2.flatten()
    target = ker_in_A2.flatten()

    # System: sum w_v M_v = target
    # M has shape (|A_2|*ker_dim, n)
    M = np.column_stack(M_vecs)

    # Solve via least squares
    w, residuals, rank, sv = np.linalg.lstsq(M, target, rcond=None)

    error = np.max(np.abs(M @ w - target))
    if error < 1e-8:
        return True, w
    else:
        return False, w


def test_position_cone(A, n):
    """Test cone from ALL positions, not just position 0.

    s_v^{(pos)}(a,b,c) inserts v at position pos:
    pos=0: (v,a,b,c)
    pos=1: (a,v,b,c) if a->v->b
    pos=2: (a,b,v,c) if b->v->c
    pos=3: (a,b,c,v) if c->v

    Does the combined map d_3 * sum_{v,pos} w_{v,pos} s_v^{(pos)} = id on ker(d_2)?
    """
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0

    if dim_O2 == 0:
        return True, None, 0

    D2 = build_full_boundary_matrix(paths2, paths1)
    D2_omega = D2 @ omega2

    U, S, Vt = np.linalg.svd(D2_omega, full_matrices=True)
    rank_d2 = sum(s > 1e-8 for s in S)
    ker_d2_basis = Vt[rank_d2:].T
    ker_dim = ker_d2_basis.shape[1] if ker_d2_basis.ndim == 2 else 0

    if ker_dim == 0:
        return True, None, 0

    ker_in_A2 = omega2 @ ker_d2_basis

    path3_idx = {p: i for i, p in enumerate(paths3)}
    D3 = build_full_boundary_matrix(paths3, paths2)

    M_vecs = []
    labels = []

    for v in range(n):
        for pos in range(4):  # Insert at position 0,1,2,3
            S_vp = np.zeros((len(paths3), len(paths2)))
            for j, (a, b, c) in enumerate(paths2):
                if v in (a, b, c):
                    continue
                if pos == 0:
                    target = (v, a, b, c)
                    ok = A[v][a] == 1
                elif pos == 1:
                    target = (a, v, b, c)
                    ok = A[a][v] == 1 and A[v][b] == 1
                elif pos == 2:
                    target = (a, b, v, c)
                    ok = A[b][v] == 1 and A[v][c] == 1
                elif pos == 3:
                    target = (a, b, c, v)
                    ok = A[c][v] == 1

                if ok and target in path3_idx:
                    S_vp[path3_idx[target], j] = 1

            M_vp = D3 @ S_vp @ ker_in_A2
            M_vecs.append(M_vp.flatten())
            labels.append((v, pos))

    target = ker_in_A2.flatten()
    M = np.column_stack(M_vecs)

    w, _, rank, _ = np.linalg.lstsq(M, target, rcond=None)
    error = np.max(np.abs(M @ w - target))

    return error < 1e-8, w if error < 1e-8 else None, ker_dim


# ============================================================
# PART 1: Test averaged cone at n=4,5
# ============================================================
print("=" * 70)
print("AVERAGED CONE HOMOTOPY TEST")
print("=" * 70)

for n_val in [4, 5]:
    m = n_val*(n_val-1)//2
    total = 1 << m
    works = 0
    fails = 0

    for bits in range(total):
        A = build_adj(n_val, bits)
        ok, err, info = test_cone_homotopy(A, n_val)
        if ok:
            works += 1
        else:
            fails += 1
            if fails <= 3:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  FAIL: n={n_val} T#{bits} scores={scores} error={err:.4f}")

    print(f"\nn={n_val}: averaged cone works: {works}/{total}, fails: {fails}/{total}")


# ============================================================
# PART 2: Weighted cone at n=4,5
# ============================================================
print(f"\n{'='*70}")
print("WEIGHTED CONE TEST (solve for vertex weights)")
print("=" * 70)

for n_val in [4, 5]:
    m = n_val*(n_val-1)//2
    total = 1 << m
    works = 0
    fails = 0
    weight_examples = []

    for bits in range(total):
        A = build_adj(n_val, bits)
        ok, w = test_weighted_cone(A, n_val)
        if ok:
            works += 1
            if len(weight_examples) < 5 and w is not None:
                scores = tuple(sorted([sum(row) for row in A]))
                weight_examples.append((bits, scores, w))
        else:
            fails += 1
            if fails <= 3:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  FAIL: n={n_val} T#{bits} scores={scores}")

    print(f"\nn={n_val}: weighted cone works: {works}/{total}, fails: {fails}/{total}")

    if weight_examples:
        print(f"  Example weights:")
        for bits, scores, w in weight_examples[:3]:
            w_str = ', '.join(f'{wi:.3f}' for wi in w)
            print(f"    T#{bits} scores={scores}: w=[{w_str}]")


# ============================================================
# PART 3: Position cone (all 4 positions) at n=4,5
# ============================================================
print(f"\n{'='*70}")
print("POSITION CONE TEST (insert v at all positions)")
print("=" * 70)

for n_val in [4, 5]:
    m = n_val*(n_val-1)//2
    total = 1 << m
    works = 0
    fails = 0

    for bits in range(total):
        A = build_adj(n_val, bits)
        ok, w, ker_dim = test_position_cone(A, n_val)
        if ok:
            works += 1
        else:
            fails += 1
            if fails <= 3:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  FAIL: n={n_val} T#{bits} scores={scores} ker_dim={ker_dim}")

    print(f"\nn={n_val}: position cone works: {works}/{total}, fails: {fails}/{total}")


print("\n\nDone.")
