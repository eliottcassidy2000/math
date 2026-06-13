#!/usr/bin/env python3
"""beta2_regular_b1.py - b1 of regular and near-regular tournaments

KEY QUESTION: For the HARD case (kappa>=2 + SC + b1=0),
all n=6 examples have near-regular score (2,2,2,3,3,3).

What about regular tournaments at odd n?
- n=5 regular: ALL b1=1 (from data)
- n=7 regular: b1=0 or b1=1?
- If regular at n=7 has b1=1: then kappa>=2 + b1=0 would need near-regular.
- If regular at n=7 has b1=0: then we need to verify good vertex existence.

Also check: what fraction of kappa>=2 + b1=0 at n=7 have near-regular scores?

Author: kind-pasteur-2026-03-08-S43
"""
import sys, os, random, time
import numpy as np
from collections import Counter
from itertools import combinations
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


def random_regular_tournament(n):
    """Generate random regular tournament on odd n."""
    assert n % 2 == 1
    while True:
        A = random_tournament(n)
        scores = [sum(A[i]) for i in range(n)]
        if all(s == (n - 1) // 2 for s in scores):
            return A


def get_induced(A, n, vertices):
    vlist = sorted(vertices)
    m = len(vlist)
    B = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            B[i][j] = A[vlist[i]][vlist[j]]
    return B, vlist


def is_sc(A, n):
    if n <= 1:
        return True
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[u][v] and v not in visited:
                visited.add(v)
                stack.append(v)
    if len(visited) < n:
        return False
    visited = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if A[v][u] and v not in visited:
                visited.add(v)
                stack.append(v)
    return len(visited) == n


def compute_b1(A, n):
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
# Part 1: b1 of regular tournaments
# ============================================================
print("=" * 70)
print("b1 OF REGULAR TOURNAMENTS")
print("=" * 70)

# n=5: all regular have b1=1 (known)
# n=7: check
for n in [5, 7, 9]:
    random.seed(42)
    if n == 5:
        # Exhaustive
        total = 2 ** (n * (n - 1) // 2)
        b1_dist = Counter()
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
            scores = [sum(A[i]) for i in range(n)]
            if not all(s == (n - 1) // 2 for s in scores):
                continue
            b1 = compute_b1(A, n)
            b1_dist[b1] += 1
        print(f"\nn={n} (exhaustive):")
        for b, cnt in sorted(b1_dist.items()):
            print(f"  b1={b}: {cnt}")
    else:
        # Random regular
        trials = {7: 200, 9: 50}[n]
        b1_dist = Counter()
        attempts = 0
        found = 0
        t0 = time.time()
        while found < trials and attempts < trials * 1000:
            attempts += 1
            A = random_tournament(n)
            scores = [sum(A[i]) for i in range(n)]
            if not all(s == (n - 1) // 2 for s in scores):
                continue
            found += 1
            b1 = compute_b1(A, n)
            b1_dist[b1] += 1
            if found % 50 == 0:
                print(f"  n={n}: {found}/{trials} ({time.time()-t0:.0f}s)")

        print(f"\nn={n} ({found} regular tournaments sampled):")
        for b, cnt in sorted(b1_dist.items()):
            print(f"  b1={b}: {cnt} ({100*cnt/found:.1f}%)")


# ============================================================
# Part 2: b1 monotonicity for regular tournaments
# ============================================================
print(f"\n{'=' * 70}")
print("b1 MONOTONICITY FOR REGULAR TOURNAMENTS")
print("=" * 70)

for n in [7]:
    random.seed(42)
    checked = 0
    all_good = 0
    all_bad = 0
    sum_b1_dist = Counter()
    t0 = time.time()
    attempts = 0

    while checked < 200 and attempts < 100000:
        attempts += 1
        A = random_tournament(n)
        scores = [sum(A[i]) for i in range(n)]
        if not all(s == (n - 1) // 2 for s in scores):
            continue
        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue  # Only check b1=0 regulars
        checked += 1

        s = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B, _ = get_induced(A, n, others)
            s += compute_b1(B, n - 1)
        sum_b1_dist[s] += 1

        if s == n:
            all_bad += 1
            print(f"  ALL BAD at trial {attempts}!")
        else:
            all_good += 1

        if checked % 50 == 0:
            print(f"  checked {checked} ({time.time()-t0:.0f}s)")

    print(f"\nn={n}: {checked} regular b1=0 tournaments checked ({time.time()-t0:.0f}s)")
    print(f"  All have good vertex: {all_good}/{checked}")
    print(f"  Sum distribution:")
    for s in sorted(sum_b1_dist.keys()):
        print(f"    Sum={s}: {sum_b1_dist[s]}")


# ============================================================
# Part 3: Score distribution of kappa>=2 + b1=0 at n=7
# ============================================================
print(f"\n{'=' * 70}")
print("KAPPA>=2 + b1=0 SCORE DISTRIBUTION at n=7")
print("=" * 70)

n = 7
random.seed(42)
score_dist = Counter()
kappa_dist = Counter()
num_found = 0
t0 = time.time()

for trial in range(3000):
    A = random_tournament(n)
    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue

    # Quick kappa check: try deleting each vertex
    all_sc = True
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        if not is_sc(B, n - 1):
            all_sc = False
            break

    if not all_sc:
        continue  # kappa = 1

    # kappa >= 2 (all deletions keep SC)
    num_found += 1
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    score_dist[scores] += 1

    if trial % 500 == 0:
        print(f"  trial {trial}, found {num_found} ({time.time()-t0:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=7 ({elapsed:.0f}s): {num_found} kappa>=2 + SC + b1=0")
print(f"  Score distribution:")
for s, cnt in sorted(score_dist.items(), key=lambda x: -x[1]):
    print(f"    {s}: {cnt}")


# ============================================================
# Part 4: NEW APPROACH - vertex with minimum c3(v)
# ============================================================
print(f"\n{'=' * 70}")
print("VERTEX WITH MIN c3(v) - IS IT ALWAYS GOOD?")
print("=" * 70)

# At n=6 kappa>=2 case: c3(T)=8, c3(T\v) = 8 - c3(v).
# If c3(T\v) <= 3 at n=5, then b1(T\v) = 0 (from our data).
# So we need c3(v) >= 5 for some v. Average = 4. Might not always hold.
# Let's check: does the vertex with MAX c3(v) always give b1(T\v)=0?

n = 6
total = 2 ** (n * (n - 1) // 2)
max_c3v_good = 0
max_c3v_bad = 0
checked = 0

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

    b1 = compute_b1(A, n)
    if b1 != 0:
        continue

    # Compute c3(v) for each v
    c3_vs = []
    for v in range(n):
        c3v = 0
        for i in range(n):
            if i == v:
                continue
            for j in range(n):
                if j == v or j == i:
                    continue
                if A[v][i] and A[i][j] and A[j][v]:
                    c3v += 1
        c3_vs.append(c3v // 2)  # each cycle counted twice

    # Vertex with maximum c3(v)
    max_v = max(range(n), key=lambda v: c3_vs[v])
    others = [x for x in range(n) if x != max_v]
    B, _ = get_induced(A, n, others)
    b1v = compute_b1(B, n - 1)

    checked += 1
    if b1v == 0:
        max_c3v_good += 1
    else:
        max_c3v_bad += 1

print(f"\nn=6: vertex with MAX c3(v) always good?")
print(f"  Good: {max_c3v_good}/{checked} ({100*max_c3v_good/checked:.1f}%)")
print(f"  Bad: {max_c3v_bad}/{checked} ({100*max_c3v_bad/checked:.1f}%)")


# Try: vertex whose deletion gives SMALLEST c3(T\v) = c3(T) - c3(v)
# I.e., vertex with LARGEST c3(v) (same as above)
# Also try: vertex with smallest degree variance in T\v

print(f"\nAlternative: vertex minimizing c3(T\\v) = vertex maximizing c3(v)")
print(f"  (Same as above)")


# Try: vertex with median degree
median_good = 0
median_bad = 0
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

    b1 = compute_b1(A, n)
    if b1 != 0:
        continue

    degs = [sum(A[v]) for v in range(n)]
    target = (n - 1) / 2
    # Sort by distance to median, pick closest
    sorted_vs = sorted(range(n), key=lambda v: abs(degs[v] - target))
    v = sorted_vs[0]
    others = [x for x in range(n) if x != v]
    B, _ = get_induced(A, n, others)
    if compute_b1(B, n - 1) == 0:
        median_good += 1
    else:
        median_bad += 1

print(f"\nVertex closest to median degree:")
print(f"  Good: {median_good}/{checked} ({100*median_good/checked:.1f}%)")
print(f"  Bad: {median_bad}/{checked} ({100*median_bad/checked:.1f}%)")


print("\n\nDone.")
