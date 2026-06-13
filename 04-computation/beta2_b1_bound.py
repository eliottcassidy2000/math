#!/usr/bin/env python3
"""beta2_b1_bound.py - Is b1(T) always <= 1 for tournaments?

DISCOVERY: At n=4,5,6 (exhaustive), b1(T) in {0, 1} always.
This would be a strong structural result if true for all n.

For the inductive proof of beta_2 = 0:
  If b1(T) <= 1 for all tournaments, and beta_2(T\v) = 0 by induction,
  then the induction step only needs: exists v with b1(T\v) <= b1(T).
  Since b1 in {0,1}, this means: exists v with b1(T\v) = 0 when b1(T) = 0,
  or b1(T\v) <= 1 when b1(T) = 1.

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
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


def compute_b1(A, n):
    """Compute beta_1 of tournament."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)

    if not paths1:
        return 0, 0, 0

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return 0, 0, 0

    D1 = build_full_boundary_matrix([tuple(p) for p in paths1], paths0)
    D1_om = D1 @ omega1
    sv = np.linalg.svd(D1_om, compute_uv=False)
    rk_d1 = int(sum(s > 1e-8 for s in sv))

    if paths2:
        omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
        dim_O2 = omega2.shape[1] if omega2.ndim == 2 else 0
        if dim_O2 > 0:
            D2 = build_full_boundary_matrix([tuple(p) for p in paths2],
                                            [tuple(p) for p in paths1])
            D2_om = D2 @ omega2
            sv2 = np.linalg.svd(D2_om, compute_uv=False)
            rk_d2 = int(sum(s > 1e-8 for s in sv2))
        else:
            rk_d2 = 0
    else:
        rk_d2 = 0

    z1_dim = dim_O1 - rk_d1
    b1 = z1_dim - rk_d2
    return b1, dim_O1, rk_d2


# ============================================================
print("=" * 70)
print("b1(T) BOUND FOR TOURNAMENTS")
print("=" * 70)

for n in [7, 8, 9, 10, 11, 12, 15, 20]:
    random.seed(42)
    num_trials = {7: 1000, 8: 500, 9: 300, 10: 200, 11: 100,
                  12: 50, 15: 20, 20: 10}[n]
    b1_dist = Counter()
    max_b1 = 0
    t0 = time.time()

    for trial in range(num_trials):
        A = random_tournament(n)
        b1, dim_O1, rk_d2 = compute_b1(A, n)
        b1_dist[b1] += 1
        if b1 > max_b1:
            max_b1 = b1

    elapsed = time.time() - t0
    z1_dim = (n - 1) * (n - 2) // 2
    print(f"\nn={n} ({elapsed:.1f}s, {num_trials} samples):")
    print(f"  dim(Z_1) = {z1_dim}")
    print(f"  b1 distribution: {dict(sorted(b1_dist.items()))}")
    print(f"  max b1 = {max_b1}")
    print(f"  b1 <= 1: {'YES' if max_b1 <= 1 else 'NO -- COUNTEREXAMPLE'}")


# ============================================================
# Part 2: b1 = 0 rate (what fraction of tournaments are acyclic in H_1?)
# ============================================================
print(f"\n{'=' * 70}")
print("b1 = 0 RATE")
print("=" * 70)

# At n=4: 40/64 = 62.5% have b1=0
# At n=5: 720/1024 = 70.3% have b1=0
# At n=6: 27968/32768 = 85.4% have b1=0

for n in [7, 8, 9, 10]:
    random.seed(42)
    num_trials = {7: 2000, 8: 1000, 9: 500, 10: 200}[n]
    b1_zero = 0
    b1_one = 0
    b1_more = 0

    for trial in range(num_trials):
        A = random_tournament(n)
        b1, _, _ = compute_b1(A, n)
        if b1 == 0:
            b1_zero += 1
        elif b1 == 1:
            b1_one += 1
        else:
            b1_more += 1

    print(f"n={n}: b1=0: {b1_zero}/{num_trials} ({100*b1_zero/num_trials:.1f}%), "
          f"b1=1: {b1_one}/{num_trials} ({100*b1_one/num_trials:.1f}%), "
          f"b1>1: {b1_more}/{num_trials}")


# ============================================================
# Part 3: For b1=1 tournaments, which vertices give b1(T\v)=0?
# ============================================================
print(f"\n{'=' * 70}")
print("b1=1 TOURNAMENTS: WHICH VERTICES REDUCE b1?")
print("=" * 70)

n = 7
random.seed(42)

def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist

b1_one_cases = 0
good_v_stats = Counter()

for trial in range(2000):
    A = random_tournament(n)
    b1, _, _ = compute_b1(A, n)
    if b1 != 1:
        continue

    b1_one_cases += 1
    good_vs = 0
    out_degs = [sum(A[v]) for v in range(n)]

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        b1v, _, _ = compute_b1(B_sub, len(others))
        if b1v <= 1:  # monotone (b1 doesn't increase)
            good_vs += 1

    good_v_stats[good_vs] += 1

    if b1_one_cases <= 3:
        print(f"  b1=1 case: scores={sorted(out_degs)}, "
              f"#good_v (b1 non-increasing)={good_vs}")

print(f"\nn={n}: {b1_one_cases} tournaments with b1=1")
print(f"  Distribution of #good_vertices: {dict(sorted(good_v_stats.items()))}")


print("\n\nDone.")
