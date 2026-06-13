#!/usr/bin/env python3
"""beta2_b1_monotonicity.py - Test b1 monotonicity under vertex deletion

KEY HYPOTHESIS (HYP-278): For every tournament T on n >= 5 vertices,
there exists a vertex v such that b1(T\v) <= b1(T).

EQUIVALENCE: b1(T\v) <= b1(T) iff i_*: H_1(T\v) -> H_1(T) is injective
             iff H_2(T, T\v) = 0.

If HYP-278 + H_2(T\v) = 0 (by induction), then H_2(T) = 0.

This would PROVE beta_2(T) = 0 for all tournaments.

ANALYSIS PLAN:
1. Exhaustive verification at n=5,6 (done: both 100%)
2. Sampled verification at n=7,8,9,10
3. Which vertices are monotone-safe? (Not just existence but structure)
4. Algebraic understanding of b1 formula

ALGEBRAIC FACTS:
- dim(Z_1(T)) = (n-1)(n-2)/2 for ALL tournaments (constant!)
- b1(T) = dim(Z_1) - rk(d_2) = (n-1)(n-2)/2 - rk(d_2|Omega_2)
- dim(Z_1(T\v)) = (n-2)(n-3)/2 (constant for (n-1)-vertex tournaments)
- b1(T\v) = (n-2)(n-3)/2 - rk(d_2|Omega_2(T\v))

So: b1(T\v) <= b1(T) iff
  (n-2)(n-3)/2 - rk(d_2(T\v)) <= (n-1)(n-2)/2 - rk(d_2(T))
  iff rk(d_2(T)) - rk(d_2(T\v)) <= (n-1)(n-2)/2 - (n-2)(n-3)/2 = n-2
  iff rk(d_2(T)) - rk(d_2(T\v)) <= n-2

So: b1(T\v) <= b1(T) iff the rank of d_2 drops by at most n-2 upon
deleting vertex v.

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


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def compute_b1(A, n):
    """Compute beta_1 of tournament."""
    paths1 = enumerate_allowed_paths(A, n, 1)
    paths0 = [(i,) for i in range(n)]
    paths2 = enumerate_allowed_paths(A, n, 2)

    if not paths1:
        return 0

    omega1 = compute_omega_basis(A, n, 1, paths1, paths0)
    dim_O1 = omega1.shape[1] if omega1.ndim == 2 else 0
    if dim_O1 == 0:
        return 0

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

    return dim_O1 - rk_d1 - rk_d2


def check_b1_monotonicity(A, n):
    """For each vertex v, check if b1(T\v) <= b1(T)."""
    b1_T = compute_b1(A, n)

    results = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        b1_Tv = compute_b1(B_sub, len(others))

        d_out = sum(A[v])
        results.append({
            'v': v, 'd_out': d_out,
            'b1_T': b1_T, 'b1_Tv': b1_Tv,
            'delta': b1_Tv - b1_T,
            'good': b1_Tv <= b1_T,
        })

    has_good = any(r['good'] for r in results)
    return b1_T, results, has_good


# ============================================================
# Part 1: Large n verification
# ============================================================
print("=" * 70)
print("b1 MONOTONICITY VERIFICATION")
print("=" * 70)

for n in [7, 8, 9, 10, 11, 12]:
    random.seed(42)
    num_trials = {7: 500, 8: 300, 9: 200, 10: 100, 11: 50, 12: 20}[n]
    has_good = 0
    no_good = 0
    max_delta = -999
    min_good_count = 999
    t0 = time.time()

    for trial in range(num_trials):
        A = random_tournament(n)
        b1_T, results, hg = check_b1_monotonicity(A, n)

        if hg:
            has_good += 1
        else:
            no_good += 1
            print(f"  NO GOOD v: n={n}, trial={trial}, b1={b1_T}")

        deltas = [r['delta'] for r in results]
        max_delta = max(max_delta, max(deltas))
        good_count = sum(1 for r in results if r['good'])
        min_good_count = min(min_good_count, good_count)

        if (trial + 1) % max(1, num_trials // 4) == 0:
            elapsed = time.time() - t0
            print(f"  n={n}: {trial+1}/{num_trials} ({elapsed:.0f}s) "
                  f"good={has_good} bad={no_good}")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {has_good}/{num_trials} have monotone v")
    print(f"  Max delta (worst increase): {max_delta}")
    print(f"  Min #good vertices in any tournament: {min_good_count}")


# ============================================================
# Part 2: Detailed delta distribution at n=7
# ============================================================
print(f"\n{'=' * 70}")
print("DELTA DISTRIBUTION at n=7")
print("=" * 70)

n = 7
random.seed(42)
delta_dist = Counter()
delta_by_degree = {}

for trial in range(500):
    A = random_tournament(n)
    b1_T, results, _ = check_b1_monotonicity(A, n)

    for r in results:
        delta_dist[r['delta']] += 1
        d = r['d_out']
        if d not in delta_by_degree:
            delta_by_degree[d] = Counter()
        delta_by_degree[d][r['delta']] += 1

print(f"Overall delta distribution:")
for d in sorted(delta_dist.keys()):
    print(f"  delta={d}: {delta_dist[d]}")

print(f"\nDelta by out-degree:")
for d in sorted(delta_by_degree.keys()):
    dist = delta_by_degree[d]
    total = sum(dist.values())
    neg = sum(v for k, v in dist.items() if k < 0)
    zero = dist.get(0, 0)
    pos = sum(v for k, v in dist.items() if k > 0)
    print(f"  d_out={d}: neg={neg} zero={zero} pos={pos} "
          f"(good_rate={(neg+zero)/total:.3f})")


# ============================================================
# Part 3: Rank drop formula
# ============================================================
print(f"\n{'=' * 70}")
print("RANK DROP ANALYSIS: rk(d_2(T)) - rk(d_2(T\\v))")
print("=" * 70)

print("""
b1(T\\v) <= b1(T) iff rk(d_2(T)) - rk(d_2(T\\v)) <= n-2.

What is rk(d_2(T)) - rk(d_2(T\\v)) typically?
""")

n = 7
random.seed(42)
rank_drops = []
good_drops = []
bad_drops = []

for trial in range(200):
    A = random_tournament(n)

    # rk(d_2(T))
    paths2 = enumerate_allowed_paths(A, n, 2)
    paths1 = enumerate_allowed_paths(A, n, 1)
    omega2 = compute_omega_basis(A, n, 2, paths2, paths1)
    if omega2.ndim < 2 or omega2.shape[1] == 0:
        continue
    D2 = build_full_boundary_matrix([tuple(p) for p in paths2], [tuple(p) for p in paths1])
    D2_om = D2 @ omega2
    sv = np.linalg.svd(D2_om, compute_uv=False)
    rk_d2_T = int(sum(s > 1e-8 for s in sv))

    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, vlist = get_induced(A, n, others)
        np2 = len(others)

        p2tv = enumerate_allowed_paths(B_sub, np2, 2)
        p1tv = enumerate_allowed_paths(B_sub, np2, 1)

        if p2tv and p1tv:
            o2tv = compute_omega_basis(B_sub, np2, 2, p2tv, p1tv)
            do2tv = o2tv.shape[1] if o2tv.ndim == 2 else 0
            if do2tv > 0:
                D2tv = build_full_boundary_matrix([tuple(p) for p in p2tv],
                                                  [tuple(p) for p in p1tv])
                D2tv_o = D2tv @ o2tv
                sv2 = np.linalg.svd(D2tv_o, compute_uv=False)
                rk_d2_Tv = int(sum(s > 1e-8 for s in sv2))
            else:
                rk_d2_Tv = 0
        else:
            rk_d2_Tv = 0

        drop = rk_d2_T - rk_d2_Tv
        rank_drops.append(drop)
        d_out = sum(A[v])
        # b1(T\v) <= b1(T) iff drop <= n-2 = 5
        if drop <= n - 2:
            good_drops.append((d_out, drop))
        else:
            bad_drops.append((d_out, drop))

print(f"n=7: {len(rank_drops)} (T,v) pairs")
print(f"  Rank drop distribution:")
for d in sorted(set(rank_drops)):
    cnt = rank_drops.count(d)
    good_bad = "GOOD" if d <= n - 2 else "BAD"
    print(f"    drop={d}: {cnt} ({good_bad})")

print(f"\n  Threshold: drop <= {n-2} is good")
print(f"  Good: {len(good_drops)}, Bad: {len(bad_drops)}")

# By degree
print(f"\n  Good drops by d_out:")
g_by_d = Counter(d for d, _ in good_drops)
b_by_d = Counter(d for d, _ in bad_drops)
for d in range(n):
    g = g_by_d.get(d, 0)
    b = b_by_d.get(d, 0)
    t = g + b
    if t > 0:
        print(f"    d_out={d}: good={g}/{t} ({100*g/t:.0f}%)")


print("\n\nDone.")
