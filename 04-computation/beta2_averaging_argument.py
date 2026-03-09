#!/usr/bin/env python3
"""beta2_averaging_argument.py - Test averaging argument for HYP-278

KEY IDEA: If sum_v (drop(v) - (n-2)) has a useful bound, then at least one
vertex must have drop(v) <= n-2, proving HYP-278.

RECALL:
- drop(v) = rk(d_2(T)) - rk(d_2(T\v))
- b1(T\v) <= b1(T) iff drop(v) <= n-2
- drop(v) = n-2 + dim(ker(i_*)) where i_*: H_1(T\v) -> H_1(T)
- When b1(T)=1: ALL vertices are good (drop <= n-2), so only b1=0 matters

APPROACH 1: Sum of drops
- Sum_v drop(v) = Sum_v [rk(d_2(T)) - rk(d_2(T\v))] = n * rk(d_2(T)) - Sum_v rk(d_2(T\v))
- If Sum_v drop(v) <= n*(n-2), then at least one v has drop <= n-2

APPROACH 2: Average b1 change
- Sum_v [b1(T\v) - b1(T)] = Sum_v b1(T\v) - n*b1(T)
- For b1=0: Sum_v b1(T\v). Need this < n (i.e., not all T\v have b1=1)

APPROACH 3: Exact formulas
- rk(d_2(T)) = dim(Z_1) - b1 = (n-1)(n-2)/2 - b1
- rk(d_2(T\v)) = (n-2)(n-3)/2 - b1(T\v)
- drop(v) = [(n-1)(n-2)/2 - b1] - [(n-2)(n-3)/2 - b1(T\v)]
          = (n-2) + b1(T\v) - b1
- So drop(v) <= n-2 iff b1(T\v) <= b1(T). (We knew this.)
- Sum_v drop(v) = n(n-2) + Sum_v b1(T\v) - n*b1(T)
- For b1=0: Sum_v drop(v) = n(n-2) + Sum_v b1(T\v)
- Average drop = (n-2) + (1/n)*Sum_v b1(T\v)
- So: at least one drop <= (n-2) iff at least one b1(T\v) = 0 (since b1 in {0,1})
- EQUIVALENTLY: HYP-278 (for b1=0 case) iff NOT ALL T\v have b1=1

So the question reduces to: Can a b1=0 tournament T have b1(T\v)=1 for ALL v?

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


# ============================================================
# Part 1: Test if ALL T\v can have b1=1 when b1(T)=0
# ============================================================
print("=" * 70)
print("CAN b1(T)=0 WITH ALL b1(T\\v)=1?")
print("=" * 70)

print("""
If b1(T)=0, HYP-278 reduces to: NOT ALL of b1(T\\v) are 1.
i.e., there exists v with b1(T\\v) = 0.

Testing this exhaustively at n=5,6 and sampled at larger n.
""")

for n in [5, 6]:
    t0 = time.time()
    all_good = 0
    all_b1_1 = 0  # b1(T)=0 but all T\v have b1=1
    b1_0_count = 0
    b1_1_count = 0
    sum_b1_tv_dist = Counter()

    total = 2 ** (n * (n - 1) // 2)
    for bits in range(total):
        # Build tournament
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        b1_T = compute_b1(A, n)
        if b1_T == 0:
            b1_0_count += 1
        else:
            b1_1_count += 1

        if b1_T > 0:
            continue  # only check b1=0 case

        # Check all T\v
        b1_tvs = []
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            b1_tvs.append(compute_b1(B_sub, n - 1))

        sum_b1 = sum(b1_tvs)
        sum_b1_tv_dist[sum_b1] += 1
        num_b1_1 = sum(1 for x in b1_tvs if x == 1)
        if num_b1_1 == n:
            all_b1_1 += 1
            print(f"  COUNTEREXAMPLE at n={n}, bits={bits}: all T\\v have b1=1!")
        if all(x == 0 for x in b1_tvs):
            all_good += 1

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s):")
    print(f"  b1=0 tournaments: {b1_0_count}/{total}")
    print(f"  b1=1 tournaments: {b1_1_count}/{total}")
    print(f"  Among b1=0: ALL T\\v have b1=0: {all_good}/{b1_0_count}")
    print(f"  Among b1=0: ALL T\\v have b1=1: {all_b1_1}/{b1_0_count} {'COUNTEREXAMPLE!' if all_b1_1 > 0 else '(none)'}")
    print(f"  Sum(b1(T\\v)) distribution for b1=0 tournaments:")
    for s in sorted(sum_b1_tv_dist.keys()):
        print(f"    Sum={s}: {sum_b1_tv_dist[s]}")


# ============================================================
# Part 2: Sampled at larger n
# ============================================================
print(f"\n{'=' * 70}")
print("SAMPLED VERIFICATION at n=7-12")
print("=" * 70)

for n in [7, 8, 9, 10, 11, 12]:
    random.seed(42)
    num_trials = {7: 500, 8: 200, 9: 100, 10: 50, 11: 30, 12: 15}[n]
    b1_0_tested = 0
    b1_0_all_b1tv_1 = 0
    b1_0_some_b1tv_0 = 0
    min_frac_good = 1.0
    t0 = time.time()

    for trial in range(num_trials):
        A = random_tournament(n)
        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue
        b1_0_tested += 1

        b1_tvs = []
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            b1_tvs.append(compute_b1(B_sub, n - 1))

        num_good = sum(1 for x in b1_tvs if x == 0)
        frac = num_good / n
        min_frac_good = min(min_frac_good, frac)

        if all(x == 1 for x in b1_tvs):
            b1_0_all_b1tv_1 += 1
            print(f"  COUNTEREXAMPLE n={n} trial={trial}")
        else:
            b1_0_some_b1tv_0 += 1

        if (trial + 1) % max(1, num_trials // 4) == 0:
            print(f"  n={n}: trial {trial+1}/{num_trials} ({time.time()-t0:.0f}s)")

    elapsed = time.time() - t0
    print(f"n={n} ({elapsed:.0f}s): b1=0 tested: {b1_0_tested}")
    print(f"  All b1(T\\v)=1: {b1_0_all_b1tv_1}/{b1_0_tested}")
    print(f"  Some b1(T\\v)=0: {b1_0_some_b1tv_0}/{b1_0_tested}")
    print(f"  Min fraction of good vertices: {min_frac_good:.3f}")


# ============================================================
# Part 3: What fraction of vertices are good when b1=0?
# ============================================================
print(f"\n{'=' * 70}")
print("FRACTION OF GOOD VERTICES WHEN b1(T)=0")
print("=" * 70)
print("At what rate do b1(T\\v)=0 vertices appear in b1(T)=0 tournaments?")

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    frac_dist = Counter()
    count_b1_0 = 0

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue
        count_b1_0 += 1

        good = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            if compute_b1(B_sub, n - 1) == 0:
                good += 1
        frac_dist[good] += 1

    print(f"\nn={n}: {count_b1_0} tournaments with b1=0")
    print(f"  #good vertices | count | fraction")
    for g in sorted(frac_dist.keys()):
        print(f"  {g}/{n}            | {frac_dist[g]:5d} | {frac_dist[g]/count_b1_0:.4f}")
    min_g = min(frac_dist.keys())
    print(f"  MINIMUM good vertices: {min_g}/{n}")


# ============================================================
# Part 4: Algebraic formula for Sum_v b1(T\v)
# ============================================================
print(f"\n{'=' * 70}")
print("SUM_v b1(T\\v) FORMULA")
print("=" * 70)
print("""
Sum_v b1(T\\v) = Sum_v [(n-2)(n-3)/2 - rk(d_2(T\\v))]
              = n*(n-2)(n-3)/2 - Sum_v rk(d_2(T\\v))

rk(d_2(T)) = (n-1)(n-2)/2 - b1(T)

Sum_v drop(v) = n*rk(d_2(T)) - Sum_v rk(d_2(T\\v))
              = n*[(n-1)(n-2)/2 - b1] - [n*(n-2)(n-3)/2 - Sum_v b1(T\\v)]
              = n(n-2)[(n-1)-(n-3)]/2 - n*b1 + Sum_v b1(T\\v)
              = n(n-2) - n*b1 + Sum_v b1(T\\v)

So Sum_v b1(T\\v) = Sum_v drop(v) - n(n-2) + n*b1

For b1=0: Sum_v b1(T\\v) = Sum_v drop(v) - n(n-2) = #bad_vertices

Let's verify and see if there's a nice formula for Sum_v b1(T\\v).
""")

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    sum_b1_tv_values = []
    score_to_sum = {}

    for bits in range(total):
        A = [[0] * n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1

        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue

        b1_sum = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            b1_sum += compute_b1(B_sub, n - 1)

        sum_b1_tv_values.append(b1_sum)
        score = tuple(sorted(sum(A[i]) for i in range(n)))
        if score not in score_to_sum:
            score_to_sum[score] = []
        score_to_sum[score].append(b1_sum)

    print(f"\nn={n}: Sum_v b1(T\\v) distribution for b1(T)=0:")
    dist = Counter(sum_b1_tv_values)
    for s in sorted(dist.keys()):
        print(f"  Sum={s}: {dist[s]}")
    print(f"  Mean: {np.mean(sum_b1_tv_values):.3f}")
    print(f"  Max: {max(sum_b1_tv_values)}")

    print(f"\n  By score sequence:")
    for score in sorted(score_to_sum.keys()):
        vals = score_to_sum[score]
        print(f"    {score}: sum in {set(vals)} (count={len(vals)})")


# ============================================================
# Part 5: 3-cycle count and b1 relationship
# ============================================================
print(f"\n{'=' * 70}")
print("3-CYCLE COUNT AND b1")
print("=" * 70)

for n in [5, 6, 7]:
    if n <= 6:
        total = 2 ** (n * (n - 1) // 2)
        num_iter = total
    else:
        num_iter = 500

    c3_b1 = {}
    random.seed(42)

    for trial in range(num_iter):
        if n <= 6:
            bits = trial
            A = [[0] * n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> idx) & 1:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    idx += 1
        else:
            A = random_tournament(n)

        # Count 3-cycles
        c3 = 0
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    cyc = 0
                    if A[i][j] and A[j][k] and A[k][i]:
                        cyc += 1
                    if A[i][k] and A[k][j] and A[j][i]:
                        cyc += 1
                    c3 += cyc

        b1 = compute_b1(A, n)
        if c3 not in c3_b1:
            c3_b1[c3] = Counter()
        c3_b1[c3][b1] += 1

    print(f"\nn={n}:")
    print(f"  c3 | b1=0 | b1=1 | frac_b1_0")
    for c3 in sorted(c3_b1.keys()):
        b0 = c3_b1[c3].get(0, 0)
        b1 = c3_b1[c3].get(1, 0)
        total = b0 + b1
        print(f"  {c3:3d} | {b0:5d} | {b1:5d} | {b0/total:.4f}")


print("\n\nDone.")
