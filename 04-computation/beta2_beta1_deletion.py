#!/usr/bin/env python3
"""
beta2_beta1_deletion.py — Test: does every tournament have interior v with β₁(T\v) = 0?

If true, this immediately gives h₂_rel(T,T\v) = 0 → β₂(T) = 0 by LES.

From earlier data: at n=5, there exist 24 tournaments (and at n=6, 960)
where ALL interior vertices have β₁(T\v) > 0. So β₁=0 is NOT sufficient!

But: we also know h₂_rel can be 0 even when β₁(T\v) > 0.
So we need a different approach.

NEW APPROACH: Count how many vertices have h₂_rel > 0 directly.
If h₂_rel > 0 requires β₁(T\v) > 0, and β₁(T\v) > 0 requires T\v
to have 3-cycles, we can bound h₂_rel through the cycle structure.

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


def compute_beta1(A, n):
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1); d2 = dim_om(om2)
    bd1 = build_full_boundary_matrix(ap1, ap0)
    rk = np.linalg.matrix_rank(bd1 @ om1, tol=1e-8)
    z1 = d1 - rk
    if d2 > 0:
        bd2om = np.linalg.lstsq(om1, build_full_boundary_matrix(ap2, ap1) @ om2, rcond=None)[0]
        b1 = np.linalg.matrix_rank(bd2om, tol=1e-8)
    else:
        b1 = 0
    return z1 - b1


print("=" * 70)
print("β₁(T\\v) = 0 EXISTENCE TEST")
print("=" * 70)

random.seed(42)

for n in [5, 6, 7]:
    m = n*(n-1)//2
    total = 1 << m

    if n <= 6:
        samples = list(range(total))
    else:
        samples = random.sample(range(total), 500)

    # Count: for how many tournaments does EVERY interior v have β₁(T\v) > 0?
    all_beta1_pos = 0
    # Count: for how many tournaments does SOME interior v have β₁(T\v) = 0?
    some_beta1_zero = 0

    beta1_zero_count_dist = Counter()  # #(interior v with β₁=0) distribution
    beta1_values_dist = Counter()  # (β₁(T), #(int v with β₁(T\v)=0)) distribution

    t0 = time.time()
    for idx, bits in enumerate(samples):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

        beta1_T = compute_beta1(A, n)

        count_beta1_zero = 0
        for v in range(n):
            if scores[v] == 0 or scores[v] == n-1:
                continue
            others = [i for i in range(n) if i != v]
            n1 = n-1
            A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
            b1_sub = compute_beta1(A_sub, n1)
            if b1_sub == 0:
                count_beta1_zero += 1

        beta1_zero_count_dist[count_beta1_zero] += 1
        beta1_values_dist[(beta1_T, count_beta1_zero)] += 1

        if count_beta1_zero == 0:
            all_beta1_pos += 1
        else:
            some_beta1_zero += 1

        if (idx + 1) % 5000 == 0:
            print(f"  n={n}: {idx+1}/{len(samples)} ({time.time()-t0:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n}: {len(samples)} tournaments in {elapsed:.0f}s")
    print(f"  ALL interior v have β₁(T\\v) > 0: {all_beta1_pos}")
    print(f"  SOME interior v has β₁(T\\v) = 0: {some_beta1_zero}")
    print(f"\n  #(interior v with β₁(T\\v)=0): {dict(sorted(beta1_zero_count_dist.items()))}")
    print(f"\n  (β₁(T), #β₁-zero-deletions):")
    for key in sorted(beta1_values_dist):
        print(f"    {key}: {beta1_values_dist[key]}")

print("\nDone.")
