#!/usr/bin/env python3
"""beta2_kappa2_analysis.py - Analyze the hard case: kappa>=2 + SC + b1=0

At n=5: this case is EMPTY (no such tournaments exist)
At n=6: 1680 tournaments, all have good vertices

KEY QUESTION: What structural property guarantees a good vertex?

PLAN:
1. Characterize the 1680 kappa>=2 + SC + b1=0 tournaments at n=6
2. Which vertices are good vs bad? Pattern?
3. Check: does some vertex always have T\v "close to transitive"?
4. Check: does min-degree or max-degree vertex always work?
5. At n=7: sample and check

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


def count_c3(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


def vertex_connectivity(A, n):
    if not is_sc(A, n):
        return 0
    if n <= 2:
        return n - 1
    for k in range(1, n):
        for subset in combinations(range(n), k):
            remaining = [v for v in range(n) if v not in subset]
            B, _ = get_induced(A, n, remaining)
            if not is_sc(B, len(remaining)):
                return k
    return n - 1


# ============================================================
# Part 1: Characterize kappa>=2 + SC + b1=0 at n=6
# ============================================================
print("=" * 70)
print("KAPPA>=2 + SC + b1=0 TOURNAMENTS at n=6")
print("=" * 70)

n = 6
total = 2 ** (n * (n - 1) // 2)
data = []
t0 = time.time()

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

    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    c3 = count_c3(A, n)

    # Check each vertex
    good_vs = []
    bad_vs = []
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        b1v = compute_b1(B, n - 1)
        d_out_v = sum(A[v])
        if b1v == 0:
            good_vs.append((v, d_out_v))
        else:
            bad_vs.append((v, d_out_v))

    data.append({
        'bits': bits, 'scores': scores, 'c3': c3,
        'kappa': kappa, 'good_vs': good_vs, 'bad_vs': bad_vs,
        'num_good': len(good_vs), 'num_bad': len(bad_vs)
    })

elapsed = time.time() - t0
print(f"\n{len(data)} tournaments found ({elapsed:.0f}s)")

# Score distribution
score_dist = Counter(d['scores'] for d in data)
print(f"\nScore distribution:")
for s, count in sorted(score_dist.items(), key=lambda x: -x[1]):
    print(f"  {s}: {count}")

# Kappa distribution
kappa_dist = Counter(d['kappa'] for d in data)
print(f"\nKappa distribution:")
for k, count in sorted(kappa_dist.items()):
    print(f"  kappa={k}: {count}")

# c3 distribution
c3_dist = Counter(d['c3'] for d in data)
print(f"\nc3 distribution:")
for c, count in sorted(c3_dist.items()):
    print(f"  c3={c}: {count}")

# Number of good/bad vertices
num_good_dist = Counter(d['num_good'] for d in data)
print(f"\n#good vertices distribution:")
for g, count in sorted(num_good_dist.items()):
    print(f"  {g}/{n}: {count}")

# Degree of good vs bad vertices
good_deg_dist = Counter()
bad_deg_dist = Counter()
for d in data:
    for v, deg in d['good_vs']:
        good_deg_dist[deg] += 1
    for v, deg in d['bad_vs']:
        bad_deg_dist[deg] += 1

print(f"\nDegree of good vertices:")
for deg in sorted(set(list(good_deg_dist.keys()) + list(bad_deg_dist.keys()))):
    g = good_deg_dist.get(deg, 0)
    b = bad_deg_dist.get(deg, 0)
    t = g + b
    print(f"  d_out={deg}: good={g}, bad={b}, rate={g/t:.3f}" if t > 0 else f"  d_out={deg}: none")


# ============================================================
# Part 2: Is min-degree or max-degree vertex always good?
# ============================================================
print(f"\n{'=' * 70}")
print("MIN/MAX DEGREE VERTEX ANALYSIS")
print("=" * 70)

min_good = 0
max_good = 0
median_good = 0

for d in data:
    scores_list = [sum(d_row) for d_row in
                   [[0]*n]]  # dummy, need actual A
    # Re-extract from bits
    bits = d['bits']
    A = [[0] * n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    degs = [sum(A[v]) for v in range(n)]
    min_v = degs.index(min(degs))
    max_v = degs.index(max(degs))
    # Median: vertex closest to (n-1)/2
    target = (n - 1) / 2
    median_v = min(range(n), key=lambda v: abs(degs[v] - target))

    good_set = {v for v, _ in d['good_vs']}
    if min_v in good_set:
        min_good += 1
    if max_v in good_set:
        max_good += 1
    if median_v in good_set:
        median_good += 1

print(f"  Min-degree vertex is good: {min_good}/{len(data)} ({100*min_good/len(data):.1f}%)")
print(f"  Max-degree vertex is good: {max_good}/{len(data)} ({100*max_good/len(data):.1f}%)")
print(f"  Median-degree vertex is good: {median_good}/{len(data)} ({100*median_good/len(data):.1f}%)")


# ============================================================
# Part 3: Check at n=7 (sampled, kappa >= 2)
# ============================================================
print(f"\n{'=' * 70}")
print("KAPPA>=2 + SC + b1=0 at n=7 (sampled)")
print("=" * 70)

n = 7
random.seed(42)
num_found = 0
num_good = 0
num_bad = 0
t0 = time.time()

for trial in range(2000):
    A = random_tournament(n)
    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    num_found += 1
    has_good = False
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        if compute_b1(B, n - 1) == 0:
            has_good = True
            break

    if has_good:
        num_good += 1
    else:
        num_bad += 1
        print(f"  COUNTEREXAMPLE at trial {trial}!")

    if num_found % 50 == 0:
        print(f"  Found {num_found} ({time.time()-t0:.0f}s)")

elapsed = time.time() - t0
print(f"\nn=7 ({elapsed:.0f}s): {num_found} kappa>=2 + SC + b1=0 found")
print(f"  Has good vertex: {num_good}/{num_found}")


# ============================================================
# Part 4: Key question - does kappa>=2 + b1=0 IMPLY b1=0
# for ALL but a bounded number of vertices?
# ============================================================
print(f"\n{'=' * 70}")
print("SUM b1(T\\v) FOR KAPPA>=2 + SC + b1=0")
print("=" * 70)

n = 6
total = 2 ** (n * (n - 1) // 2)
sum_dist = Counter()

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

    if not is_sc(A, n):
        continue
    b1 = compute_b1(A, n)
    if b1 != 0:
        continue
    kappa = vertex_connectivity(A, n)
    if kappa < 2:
        continue

    s = 0
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B, _ = get_induced(A, n, others)
        s += compute_b1(B, n - 1)
    sum_dist[s] += 1

print(f"\nn=6: Sum_v b1(T\\v) for kappa>=2:")
for s in sorted(sum_dist.keys()):
    print(f"  Sum={s}: {sum_dist[s]}")


print("\n\nDone.")
