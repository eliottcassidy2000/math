#!/usr/bin/env python3
"""
beta2_dt_cancel_n7.py — Test DT+cancel filling of Z₂ at n=7 (sampled)

Known results:
- n=4,5: DT alone fills Z₂ (100%)
- n=6: DT+cancel fills Z₂ (100%), DT alone fails 3%
- n=7: β₂=0 verified exhaustively, but does DT+cancel suffice?

Author: opus-2026-03-08-S49
"""
import sys, time, random
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved


def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


print("=" * 70)
print("DT+CANCEL FILLING AT n=7 (SAMPLED)")
print("=" * 70)

random.seed(42)

for n in [7, 8]:
    m = n*(n-1)//2
    total = 1 << m
    n_samples = 1000

    samples = random.sample(range(total), min(n_samples, total))

    dt_ok = 0
    dt_fail = 0
    cancel_ok = 0
    cancel_fail = 0
    trivial = 0

    deficit_dist = Counter()

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)

        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0 or not ap3:
            trivial += 1
            continue

        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk2
        if z2_dim == 0:
            trivial += 1
            continue

        # Build DT vectors
        bd3 = build_full_boundary_matrix(ap3, ap2)
        dt_indices = [i for i, p in enumerate(ap3) if A[p[0]][p[2]] and A[p[1]][p[3]]]

        # DT filling rank
        if dt_indices:
            bd3_dt = np.column_stack([bd3[:, i] for i in dt_indices])
            coords_dt = np.linalg.lstsq(om2, bd3_dt, rcond=None)[0]
            rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)
        else:
            rk_dt = 0

        deficit = z2_dim - rk_dt
        deficit_dist[deficit] += 1

        if deficit == 0:
            dt_ok += 1
        else:
            dt_fail += 1

            # Build cancel vectors
            bad_groups = defaultdict(list)
            for i, p in enumerate(ap3):
                a, b, c, d = p
                if not A[a][c]: bad_groups[('02', a, c)].append(i)
                if not A[b][d]: bad_groups[('13', b, d)].append(i)

            cancel_vecs = []
            for key, indices in bad_groups.items():
                for j in range(1, len(indices)):
                    v = np.zeros(len(ap3))
                    v[indices[0]] = 1; v[indices[j]] = -1
                    cancel_vecs.append(v)

            # DT + cancel
            all_vecs = []
            for i in dt_indices:
                v = np.zeros(len(ap3)); v[i] = 1
                all_vecs.append(v)
            all_vecs.extend(cancel_vecs)

            if all_vecs:
                V = np.column_stack(all_vecs)
                bd3_all = bd3 @ V
                coords_all = np.linalg.lstsq(om2, bd3_all, rcond=None)[0]
                rk_all = np.linalg.matrix_rank(coords_all, tol=1e-8)
            else:
                rk_all = 0

            if rk_all >= z2_dim:
                cancel_ok += 1
            else:
                cancel_fail += 1
                scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
                print(f"  CANCEL FAIL: bits={bits}, scores={scores}, Z₂={z2_dim}, "
                      f"rk_DT={rk_dt}, rk_DT+cancel={rk_all}, deficit={deficit}")

        if (idx + 1) % 200 == 0:
            elapsed = time.time() - t0
            print(f"  ... {idx+1}/{len(samples)} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(samples)} samples in {elapsed:.0f}s")
    print(f"  Trivial (no Z₂): {trivial}")
    print(f"  DT alone fills Z₂: {dt_ok}")
    print(f"  DT fails, cancel fixes: {cancel_ok}")
    print(f"  DT+cancel FAILS: {cancel_fail}")
    print(f"  DT deficit distribution: {dict(sorted(deficit_dist.items()))}")

    if cancel_fail == 0:
        print(f"\n  ✓ DT+cancel ALWAYS fills Z₂ at n={n} ({dt_ok+cancel_ok} nontrivial tests)")
    else:
        print(f"\n  ✗ DT+cancel FAILS {cancel_fail} times at n={n}")

    # Summary rate
    total_nontrivial = dt_ok + dt_fail
    if total_nontrivial > 0:
        dt_rate = 100 * dt_ok / total_nontrivial
        print(f"\n  DT alone success rate: {dt_rate:.1f}%")

print("\nDone.")
