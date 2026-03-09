#!/usr/bin/env python3
"""beta2_inductive_proof.py - Toward an inductive proof of beta_2 = 0

The LES of the pair (T, T\v) gives:
  ... -> H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T) -> ...

If we can show:
  (A) H_2(T\v) = 0 for all (n-1)-vertex tournaments (induction)
  (B) The connecting map delta: H_2(T,T\v) -> H_1(T\v) is injective
Then: H_2(T) = ker(delta) = 0.

(A) is the induction hypothesis (base case n=3: trivially beta_2=0).
(B) is HYP-259: delta injective for all interior v (1 <= d+ <= n-2).

But we showed delta is NOT always injective (HYP-248 REFUTED for source/sink).
For INTERIOR v, delta IS injective (HYP-259 confirmed).

Does every tournament have an interior vertex? YES: for n >= 3, there are at
most 1 source and 1 sink (by tournament property), so n-2 interior vertices.

So the proof would be:
1. Base case: beta_2 = 0 for n = 3.
2. Inductive step: Assume beta_2(T') = 0 for all |T'| = n-1.
   For tournament T on n vertices:
   - Pick any interior vertex v.
   - LES gives: 0 -> H_2(T) -> H_2(T,T\v) -delta-> H_1(T\v) -> ...
   - If delta is injective, then H_2(T) = 0.

The KEY remaining question: is delta REALLY injective for all interior v?
Let's test this more carefully at n=7,8,9.

Also: analyze the MECHANISM of delta-injectivity.

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis
)
sys.stdout = _saved

random.seed(42)


def random_tournament(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_beta2(A, n):
    """Compute beta_2 directly."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths3 = enumerate_allowed_paths(A, n, 3)

    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    omega3 = compute_omega_basis(A, n, 3, paths3, paths2)
    dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
    dim_O3 = omega3.shape[1] if omega3.ndim == 2 else 0

    if dim_O2 == 0:
        return 0

    D2 = build_full_boundary_matrix(
        [tuple(p) for p in paths2], [tuple(p) for p in paths1]
    )
    D2_om = D2 @ omega2
    r2 = np.linalg.matrix_rank(D2_om, tol=1e-8)
    ker_d2 = dim_O2 - r2

    if ker_d2 == 0:
        return 0

    if dim_O3 == 0:
        return ker_d2

    D3 = build_full_boundary_matrix(
        [tuple(p) for p in paths3], [tuple(p) for p in paths2]
    )
    D3_om = D3 @ omega3

    _, _, Vt2 = np.linalg.svd(D2_om, full_matrices=True)
    ker_basis = omega2 @ Vt2[r2:].T

    combined = np.hstack([D3_om, ker_basis])
    r_c = np.linalg.matrix_rank(combined, tol=1e-8)
    r_d3 = np.linalg.matrix_rank(D3_om, tol=1e-8)

    return r_c - r_d3


def compute_relative_h2_and_delta(A, n, v):
    """Compute H_2(T, T\v) and the connecting map delta.

    Returns: (dim_h2_rel, delta_injective, dim_h1_Tv, beta_1_T)
    """
    # Full tournament T
    paths1_T = enumerate_allowed_paths(A, n, 1)
    paths2_T = enumerate_allowed_paths(A, n, 2)
    paths3_T = enumerate_allowed_paths(A, n, 3)

    # T\v
    others = [x for x in range(n) if x != v]
    B, vlist = get_induced(A, n, others)
    n_prime = len(others)

    paths1_Tv = enumerate_allowed_paths(B, n_prime, 1)
    paths2_Tv = enumerate_allowed_paths(B, n_prime, 2)

    # Map T\v paths to T indices
    def map_path(p, vlist):
        return tuple(vlist[x] for x in p)

    # Omega bases
    omega1_T = compute_omega_basis(A, n, 1, paths1_T, [(i,) for i in range(n)])
    omega2_T = compute_omega_basis(A, n, 2, paths2_T, paths1_T)
    omega1_Tv = compute_omega_basis(B, n_prime, 1, paths1_Tv,
                                    [(i,) for i in range(n_prime)])
    omega2_Tv = compute_omega_basis(B, n_prime, 2, paths2_Tv, paths1_Tv)

    dim_O1_T = omega1_T.shape[1] if omega1_T.ndim == 2 else 0
    dim_O2_T = omega2_T.shape[1] if omega2_T.ndim == 2 else 0
    dim_O1_Tv = omega1_Tv.shape[1] if omega1_Tv.ndim == 2 else 0
    dim_O2_Tv = omega2_Tv.shape[1] if omega2_Tv.ndim == 2 else 0

    # beta_1(T\v)
    if dim_O1_Tv == 0:
        beta_1_Tv = 0
    else:
        D1_Tv = build_full_boundary_matrix(
            [tuple(p) for p in paths1_Tv],
            [(i,) for i in range(n_prime)]
        )
        D1_Tv_om = D1_Tv @ omega1_Tv
        r1_Tv = np.linalg.matrix_rank(D1_Tv_om, tol=1e-8)
        beta_1_Tv = dim_O1_Tv - r1_Tv

    # beta_1(T)
    if dim_O1_T == 0:
        beta_1_T = 0
    else:
        D1_T = build_full_boundary_matrix(
            [tuple(p) for p in paths1_T],
            [(i,) for i in range(n)]
        )
        D1_T_om = D1_T @ omega1_T
        r1_T = np.linalg.matrix_rank(D1_T_om, tol=1e-8)
        beta_1_T = dim_O1_T - r1_T

    # Relative chains: Omega_p(T) / Omega_p(T\v)
    # The v-chain space = {z in Omega_p(T) : z uses v}
    # H_2(T,T\v) = ker(d_2^rel) / im(d_3^rel)
    # where d^rel operates on relative chains (those using v)

    # Identify v-using 2-paths in T
    path2_T_idx = {tuple(p): i for i, p in enumerate(paths2_T)}
    v_paths2 = [i for i, p in enumerate(paths2_T) if v in tuple(p)]

    if dim_O2_T == 0 or not v_paths2:
        return 0, True, beta_1_Tv, beta_1_T

    # ker(d_2) restricted to v-chains in Omega_2
    # The v-chains in Omega_2: project omega2 basis to v-path coordinates
    # Actually, we need: relative Omega_2 = Omega_2(T) / image_of(Omega_2(T\v))

    # Simpler approach: compute H_2(T,T\v) via LES
    # From the LES: H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)
    # Since beta_2(T) = 0 (verified), H_2(T,T\v) injects into ker(H_1(T\v)->H_1(T))

    # The map H_1(T\v) -> H_1(T) is induced by inclusion
    # dim(H_2(T,T\v)) = dim(ker(i_*: H_1(T\v)->H_1(T)))
    #                  + dim(coker(H_2(T)->H_2(T,T\v)))
    # Since H_2(T)=0: dim(H_2(T,T\v)) = dim(ker(i_*: H_1(T\v)->H_1(T)))

    # Compute ker(i_*: H_1(T\v) -> H_1(T))
    if beta_1_Tv == 0:
        return 0, True, beta_1_Tv, beta_1_T

    # Get H_1(T\v) cycles in T\v coordinates, embed in T coordinates
    D1_Tv = build_full_boundary_matrix(
        [tuple(p) for p in paths1_Tv],
        [(i,) for i in range(n_prime)]
    )
    D1_Tv_om = D1_Tv @ omega1_Tv
    _, S1, Vt1 = np.linalg.svd(D1_Tv_om, full_matrices=True)
    r1 = sum(s > 1e-8 for s in S1)

    if r1 >= dim_O1_Tv:
        return 0, True, beta_1_Tv, beta_1_T

    ker_h1_Tv = omega1_Tv @ Vt1[r1:].T  # columns are H_1(T\v) in T\v coords

    # Map these to T's 1-path coordinates
    path1_T_idx = {tuple(p): i for i, p in enumerate(paths1_T)}
    embedding = np.zeros((len(paths1_T), len(paths1_Tv)))
    for j, p in enumerate(paths1_Tv):
        p_orig = map_path(tuple(p), vlist)
        if p_orig in path1_T_idx:
            embedding[path1_T_idx[p_orig], j] = 1

    ker_in_T = embedding @ ker_h1_Tv  # H_1(T\v) cycles in T coordinates

    # Check which are boundaries in T (i.e., in im(d_1^T))
    D1_T = build_full_boundary_matrix(
        [tuple(p) for p in paths1_T],
        [(i,) for i in range(n)]
    )
    if dim_O1_T > 0:
        D1_T_om = D1_T @ omega1_T
        # Is ker_in_T in im(d_1|_Om1)?
        combined = np.hstack([D1_T_om, ker_in_T])
        r_c = np.linalg.matrix_rank(combined, tol=1e-8)
        r_d1 = np.linalg.matrix_rank(D1_T_om, tol=1e-8)
        dim_h2_rel = r_c - r_d1  # cycles in T\v that become boundaries in T
    else:
        dim_h2_rel = np.linalg.matrix_rank(ker_in_T, tol=1e-8)

    # Actually, H_2(T,T\v) = ker(i_*) where i_*: H_1(T\v) -> H_1(T)
    # ker(i_*) = {z in H_1(T\v) : i(z) in im(d_2^T) + ... }
    # Wait, i_*: H_1(T\v) -> H_1(T) maps z to its image in H_1(T)
    # z maps to 0 in H_1(T) iff z (embedded in T coords) is a boundary in T

    # Check each basis vector
    if dim_O2_T > 0:
        D2_T = build_full_boundary_matrix(
            [tuple(p) for p in paths2_T], [tuple(p) for p in paths1_T]
        )
        D2_T_om = D2_T @ omega2_T
        im_d2 = D2_T_om
    else:
        im_d2 = np.zeros((len(paths1_T), 0))

    # ker(i_*) = columns of ker_in_T that are in im_d2
    combined2 = np.hstack([im_d2, ker_in_T])
    r_i = np.linalg.matrix_rank(im_d2, tol=1e-8)
    r_c2 = np.linalg.matrix_rank(combined2, tol=1e-8)
    ker_istar = ker_in_T.shape[1] - (r_c2 - r_i)

    # Wait, this isn't right either. Let me be more careful.
    # ker(i_*) = {z in H_1(T\v) : embedded z is in im(d_2^T|_Omega)}
    # = {z : exists w in Omega_2(T) with d_2(w) = embed(z)}
    # dim(ker(i_*)) = dim(H_1(T\v)) - dim(im(i_*))
    # dim(im(i_*)) = dim(span of embedded H_1(T\v) in H_1(T))
    # Actually we need: dim of image of embedded ker in the QUOTIENT H_1(T)
    # = dim(embedded ker + im(d_2))/dim(im(d_2)) - dim of intersection

    # Simpler: ker(i_*) = those z where embed(z) is in im(d_2)
    n_ker = 0
    for k in range(ker_h1_Tv.shape[1]):
        z_T = ker_in_T[:, k]
        if np.max(np.abs(z_T)) < 1e-12:
            n_ker += 1
            continue
        if im_d2.shape[1] > 0:
            c, _, _, _ = np.linalg.lstsq(im_d2, z_T, rcond=None)
            err = np.max(np.abs(im_d2 @ c - z_T))
            if err < 1e-6:
                n_ker += 1

    dim_h2_rel = n_ker
    delta_injective = (dim_h2_rel == 0)

    return dim_h2_rel, delta_injective, beta_1_Tv, beta_1_T


# ============================================================
# PART 1: Delta injectivity at n=7,8
# ============================================================
print("=" * 70)
print("PART 1: DELTA INJECTIVITY FOR INTERIOR VERTICES")
print("=" * 70)
print("For each tournament, for each interior vertex v (1<=d+<=n-2),")
print("check if delta: H_2(T,T\\v) -> H_1(T\\v) is injective.")
print("(Equiv: all H_1(T\\v) cycles remain non-trivial in H_1(T))")

for n in [7, 8]:
    random.seed(42)
    num_trials = 300 if n == 7 else 100
    total_interior = 0
    delta_inj_ok = 0
    delta_inj_fail = 0
    h2_rel_nonzero = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)

        for v in range(n):
            d_out = sum(A[v])
            if d_out == 0 or d_out == n - 1:
                continue  # skip source/sink

            total_interior += 1
            h2_rel, inj, b1_Tv, b1_T = compute_relative_h2_and_delta(A, n, v)

            if h2_rel > 0:
                h2_rel_nonzero += 1

            if not inj:
                delta_inj_fail += 1
                if delta_inj_fail <= 3:
                    scores = tuple(sorted([sum(row) for row in A]))
                    print(f"  n={n} FAIL trial {trial}, v={v}, d+={d_out}")
                    print(f"    h2_rel={h2_rel}, b1_Tv={b1_Tv}, b1_T={b1_T}")
                    print(f"    scores={scores}")
            else:
                delta_inj_ok += 1

        if (trial + 1) % (num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"interior={total_interior} ok={delta_inj_ok} "
                  f"fail={delta_inj_fail} h2>0={h2_rel_nonzero}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {total_interior} interior vertices")
    print(f"  Delta injective: {delta_inj_ok}")
    print(f"  Delta NOT injective: {delta_inj_fail}")
    print(f"  H_2(T,T\\v) > 0: {h2_rel_nonzero}")


# ============================================================
# PART 2: What does beta_2=0 look like in the LES?
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: LES STRUCTURE (n=5 exhaustive)")
print("=" * 70)

def all_tournaments_gen(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

n = 5
les_patterns = {}
total = 0

for A in all_tournaments_gen(n):
    total += 1
    b1_T = 0
    b2_T = compute_beta2(A, n)

    # Compute beta_1(T)
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    if omega1.ndim == 2 and omega1.shape[1] > 0:
        D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
        D1_om = D1 @ omega1
        r1 = np.linalg.matrix_rank(D1_om, tol=1e-8)
        b1_T = omega1.shape[1] - r1

    for v in range(n):
        d_out = sum(A[v])
        if d_out == 0 or d_out == n - 1:
            continue

        h2_rel, inj, b1_Tv, _ = compute_relative_h2_and_delta(A, n, v)
        key = (b1_T, b1_Tv, h2_rel, inj)
        les_patterns[key] = les_patterns.get(key, 0) + 1

print(f"n=5: {total} tournaments")
print(f"\n  LES patterns (b1_T, b1_Tv, h2_rel, delta_inj):")
for key, count in sorted(les_patterns.items(), key=lambda x: -x[1]):
    b1T, b1Tv, h2r, inj = key
    print(f"    b1_T={b1T}, b1_Tv={b1Tv}, h2_rel={h2r}, "
          f"delta_inj={'Y' if inj else 'N'}: {count}")


print("\n\nDone.")
