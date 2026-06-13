#!/usr/bin/env python3
"""beta2_les_test.py - Test LES injectivity for beta_2=0 proof

KEY FACT from LES of (T, T\v):
  0 = H_2(T\v) -> H_2(T) -> H_2(T,T\v) -> H_1(T\v) -> H_1(T)

where H_2(T\v) = 0 by induction.

So: H_2(T) injects into H_2(T,T\v), and H_2(T,T\v) = ker(delta).
delta: H_2(T,T\v) -> H_1(T\v) is the connecting map.
H_2(T) = 0 iff delta is injective.

But actually H_2(T) = ker(j_*: H_2(T) -> H_2(T,T\v)) where j_* comes
from the LES. Since H_2(T\v) -> H_2(T) is injective and H_2(T\v)=0,
H_2(T) injects into H_2(T,T\v).

So beta_2(T) = 0 iff ker(delta) = 0 iff delta: H_2(T,T\v) -> H_1(T\v)
is injective. And by LES exactness, im(delta) = ker(i_*: H_1(T\v)->H_1(T)).

So H_2(T,T\v) = ker(delta) + im(j_*) where H_2(T) -j_*-> H_2(T,T\v).
Since H_2(T\v) = 0, j_* is injective, so H_2(T) embeds in H_2(T,T\v).
And H_2(T) = ker(delta composed with j_*^{-1}) on im(j_*).

Bottom line: beta_2(T) <= dim(ker(i_*: H_1(T\v) -> H_1(T)))

TEST: For each interior v, is H_1(T\v) -> H_1(T) injective?
If yes for ALL interior v, then H_2(T,T\v) = 0 for all, hence H_2(T) = 0.

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


def check_h1_injectivity(A, n, v):
    """Check if i_*: H_1(T\v) -> H_1(T) is injective.

    Returns (b1_Tv, b1_T, ker_istar_dim, is_injective)
    """
    # T\v
    others = [x for x in range(n) if x != v]
    B, vlist = get_induced(A, n, others)
    n_prime = len(others)

    # 1-paths and boundaries for T\v
    paths1_Tv = enumerate_allowed_paths(B, n_prime, 1)
    paths2_Tv = enumerate_allowed_paths(B, n_prime, 2)

    if not paths1_Tv:
        return 0, 0, 0, True

    # Omega_1(T\v)
    paths0_Tv = [(i,) for i in range(n_prime)]
    omega1_Tv = compute_omega_basis(B, n_prime, 1, paths1_Tv, paths0_Tv)
    dim_O1_Tv = omega1_Tv.shape[1] if omega1_Tv.ndim == 2 else 0

    if dim_O1_Tv == 0:
        return 0, 0, 0, True

    # d_1 for T\v
    D1_Tv = build_full_boundary_matrix(
        [tuple(p) for p in paths1_Tv], paths0_Tv
    )
    D1_Tv_om = D1_Tv @ omega1_Tv
    r1_Tv = np.linalg.matrix_rank(D1_Tv_om, tol=1e-8)
    b1_Tv = dim_O1_Tv - r1_Tv

    if b1_Tv == 0:
        return 0, 0, 0, True

    # Get H_1(T\v) cycle basis (in T\v coordinates)
    _, _, Vt1 = np.linalg.svd(D1_Tv_om, full_matrices=True)
    ker_cycles_Tv = omega1_Tv @ Vt1[r1_Tv:].T  # cols are cycles in T\v A_1

    # 1-paths for T (full)
    paths1_T = enumerate_allowed_paths(A, n, 1)
    paths2_T = enumerate_allowed_paths(A, n, 2)
    paths0_T = [(i,) for i in range(n)]

    omega1_T = compute_omega_basis(A, n, 1, paths1_T, paths0_T)
    dim_O1_T = omega1_T.shape[1] if omega1_T.ndim == 2 else 0

    D1_T = build_full_boundary_matrix(
        [tuple(p) for p in paths1_T], paths0_T
    )

    # Embed T\v 1-paths into T's A_1 space
    path1_T_idx = {tuple(p): i for i, p in enumerate(paths1_T)}
    path1_Tv_list = [tuple(p) for p in paths1_Tv]

    embedding = np.zeros((len(paths1_T), len(paths1_Tv)))
    for j, p_Tv in enumerate(path1_Tv_list):
        p_T = tuple(vlist[x] for x in p_Tv)
        if p_T in path1_T_idx:
            embedding[path1_T_idx[p_T], j] = 1

    # Map cycles to T coordinates
    cycles_in_T = embedding @ ker_cycles_Tv  # cols are H_1(T\v) cycles in T's A_1

    # Check if these are cycles in T (they should be: d_1 commutes with inclusion)
    d1_check = D1_T @ cycles_in_T
    cycle_err = np.max(np.abs(d1_check))

    # Compute d_2 for T
    omega2_T = compute_omega_basis(A, n, 2, paths2_T, paths1_T)
    dim_O2_T = omega2_T.shape[1] if omega2_T.ndim == 2 else 0

    if dim_O2_T > 0:
        D2_T = build_full_boundary_matrix(
            [tuple(p) for p in paths2_T], [tuple(p) for p in paths1_T]
        )
        D2_T_om = D2_T @ omega2_T  # im(d_2) in A_1(T)
    else:
        D2_T_om = np.zeros((len(paths1_T), 0))

    # For each cycle z in H_1(T\v), check if its embedding is in im(d_2(T))
    # If z_T in im(d_2) => z maps to 0 in H_1(T) => in ker(i_*)
    ker_istar = 0
    for k in range(b1_Tv):
        z_T = cycles_in_T[:, k]
        if np.max(np.abs(z_T)) < 1e-12:
            ker_istar += 1
            continue

        if D2_T_om.shape[1] > 0:
            c, _, _, _ = np.linalg.lstsq(D2_T_om, z_T, rcond=None)
            err = np.max(np.abs(D2_T_om @ c - z_T))
            if err < 1e-6:
                ker_istar += 1

    # beta_1(T)
    if dim_O1_T > 0:
        D1_T_om = D1_T @ omega1_T
        r1_T = np.linalg.matrix_rank(D1_T_om, tol=1e-8)
        b1_T = dim_O1_T - r1_T
    else:
        b1_T = 0

    return b1_Tv, b1_T, ker_istar, (ker_istar == 0)


# ============================================================
# PART 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("PART 1: H_1 INJECTIVITY at n=5 (EXHAUSTIVE)")
print("=" * 70)

n = 5
total_interior = 0
inj_ok = 0
inj_fail = 0
has_interior_inj = 0
total_tours = 0

for A in all_tournaments_gen(n):
    total_tours += 1
    tour_ok = False

    for v in range(n):
        d_out = sum(A[v])
        if d_out == 0 or d_out == n - 1:
            continue

        total_interior += 1
        b1_Tv, b1_T, ker_dim, inj = check_h1_injectivity(A, n, v)

        if inj:
            inj_ok += 1
            tour_ok = True
        else:
            inj_fail += 1

    if tour_ok:
        has_interior_inj += 1

print(f"n=5: {total_tours} tournaments, {total_interior} interior (T,v) pairs")
print(f"  i_* injective: {inj_ok}/{total_interior}")
print(f"  i_* NOT injective: {inj_fail}/{total_interior}")
print(f"  Tournaments with SOME interior v having inj i_*: "
      f"{has_interior_inj}/{total_tours}")


# ============================================================
# PART 2: n=6 exhaustive
# ============================================================
print(f"\n{'=' * 70}")
print("PART 2: H_1 INJECTIVITY at n=6 (EXHAUSTIVE)")
print("=" * 70)

n = 6
total_interior = 0
inj_ok = 0
inj_fail = 0
has_interior_inj = 0
total_tours = 0
t0 = time.time()

for tidx, A in enumerate(all_tournaments_gen(n)):
    total_tours += 1
    tour_ok = False

    for v in range(n):
        d_out = sum(A[v])
        if d_out == 0 or d_out == n - 1:
            continue

        total_interior += 1
        b1_Tv, b1_T, ker_dim, inj = check_h1_injectivity(A, n, v)

        if inj:
            inj_ok += 1
            tour_ok = True
        else:
            inj_fail += 1
            if inj_fail <= 5:
                scores = tuple(sorted([sum(row) for row in A]))
                print(f"  FAIL T#{tidx} v={v}: b1_Tv={b1_Tv}, b1_T={b1_T}, "
                      f"ker={ker_dim}, scores={scores}")

    if tour_ok:
        has_interior_inj += 1

    if (tidx + 1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {tidx + 1}/32768 ({elapsed:.0f}s) "
              f"inj_ok={inj_ok} inj_fail={inj_fail}")

elapsed = time.time() - t0
print(f"\nn=6 ({elapsed:.0f}s): {total_tours} tournaments, "
      f"{total_interior} interior pairs")
print(f"  i_* injective: {inj_ok}/{total_interior}")
print(f"  i_* NOT injective: {inj_fail}/{total_interior}")
print(f"  Tournaments with SOME interior v having inj i_*: "
      f"{has_interior_inj}/{total_tours}")


# ============================================================
# PART 3: n=7,8 sampled
# ============================================================
print(f"\n{'=' * 70}")
print("PART 3: H_1 INJECTIVITY at n=7,8 (sampled)")
print("=" * 70)

for n in [7, 8]:
    random.seed(42)
    num_trials = 500 if n == 7 else 200
    total_interior = 0
    inj_ok = 0
    inj_fail = 0
    has_interior_inj = 0
    no_interior_inj = 0

    t0 = time.time()
    for trial in range(num_trials):
        A = random_tournament(n)
        tour_ok = False

        for v in range(n):
            d_out = sum(A[v])
            if d_out == 0 or d_out == n - 1:
                continue

            total_interior += 1
            b1_Tv, b1_T, ker_dim, inj = check_h1_injectivity(A, n, v)

            if inj:
                inj_ok += 1
                tour_ok = True
            else:
                inj_fail += 1

        if tour_ok:
            has_interior_inj += 1
        else:
            no_interior_inj += 1

        if (trial + 1) % (num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial + 1}/{num_trials} ({elapsed:.0f}s) "
                  f"ok={inj_ok} fail={inj_fail}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {total_interior} interior pairs")
    print(f"  i_* injective: {inj_ok}/{total_interior}")
    print(f"  i_* NOT injective: {inj_fail}/{total_interior}")
    print(f"  ALL interior v inj: {has_interior_inj}/{num_trials}")
    print(f"  NO interior v inj: {no_interior_inj}/{num_trials}")


print("\n\nDone.")
