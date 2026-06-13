#!/usr/bin/env python3
"""beta2_sum_b1_bound.py - Investigate bound on Sum_v b1(T\v)

KEY DISCOVERY: For b1(T)=0 tournaments:
  Sum_v b1(T\v) <= 3 at n=5,6 (exhaustive)

Conjecture: Sum_v b1(T\v) <= floor(n/2) or something similar?

If Sum_v b1(T\v) < n, then HYP-278 follows (exists good vertex).

We also investigate: WHICH vertices are bad? Is there a vertex-local
condition that guarantees goodness?

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


def count_3cycles(A, n):
    """Count total directed 3-cycles."""
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    c3 += 1
    return c3


def c3_through_v(A, n, v):
    """Count 3-cycles through vertex v."""
    c3v = 0
    for i in range(n):
        if i == v:
            continue
        for j in range(n):
            if j == v or j == i:
                continue
            # Check cycle v -> i -> j -> v
            if A[v][i] and A[i][j] and A[j][v]:
                c3v += 1
    # Each cycle counted 2 times (once per ordering of the other 2 vertices)
    return c3v // 2


# ============================================================
# Part 1: Sum_v b1(T\v) vs tournament invariants
# ============================================================
print("=" * 70)
print("SUM_v b1(T\\v) vs TOURNAMENT INVARIANTS")
print("=" * 70)

for n in [5, 6]:
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

        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue

        # Compute b1(T\v) for each v
        bad_vertices = []
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            if compute_b1(B_sub, n - 1) == 1:
                bad_vertices.append(v)

        sum_b1tv = len(bad_vertices)

        # Tournament invariants
        scores = sorted(sum(A[i]) for i in range(n))
        c3 = count_3cycles(A, n)

        # c3 through each bad vertex
        c3_bad = [c3_through_v(A, n, v) for v in bad_vertices] if bad_vertices else []
        # c3 through each good vertex
        good_vs = [v for v in range(n) if v not in bad_vertices]
        c3_good = [c3_through_v(A, n, v) for v in good_vs]

        # Degree of bad/good vertices
        deg_bad = [sum(A[v]) for v in bad_vertices]
        deg_good = [sum(A[v]) for v in good_vs]

        data.append({
            'bits': bits, 'scores': tuple(scores), 'c3': c3,
            'sum_b1tv': sum_b1tv, 'bad_vertices': bad_vertices,
            'c3_bad': c3_bad, 'c3_good': c3_good,
            'deg_bad': deg_bad, 'deg_good': deg_good
        })

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {len(data)} b1=0 tournaments")

    # Max sum
    max_sum = max(d['sum_b1tv'] for d in data)
    print(f"  Max Sum_v b1(T\\v): {max_sum}")
    print(f"  n-2 = {n-2}")
    print(f"  floor(n/2) = {n//2}")

    # Distribution
    sum_dist = Counter(d['sum_b1tv'] for d in data)
    for s in sorted(sum_dist.keys()):
        print(f"  Sum={s}: {sum_dist[s]}")

    # Bad vertex characterization
    print(f"\n  Bad vertex degree distribution:")
    all_bad_degs = []
    for d in data:
        all_bad_degs.extend(d['deg_bad'])
    if all_bad_degs:
        for deg in sorted(set(all_bad_degs)):
            print(f"    d_out={deg}: {all_bad_degs.count(deg)}")

    print(f"\n  Good vertex degree distribution:")
    all_good_degs = []
    for d in data:
        all_good_degs.extend(d['deg_good'])
    if all_good_degs:
        for deg in sorted(set(all_good_degs)):
            print(f"    d_out={deg}: {all_good_degs.count(deg)}")

    # c3 through bad vs good
    print(f"\n  c3 through bad vertices: {Counter(x for d in data for x in d['c3_bad'])}")
    print(f"  c3 through good vertices: {Counter(x for d in data for x in d['c3_good'])}")

    # Check: is there a LOCAL condition for goodness?
    # Hypothesis: vertex v is good if c3(v) >= threshold?
    print(f"\n  c3(v) vs goodness:")
    c3v_good = Counter()
    c3v_bad = Counter()
    for d in data:
        for v in d['bad_vertices']:
            c3v_bad[c3_through_v_fast(A, n, v) if False else d['c3_bad'][d['bad_vertices'].index(v)]] += 1
        for i, v in enumerate([x for x in range(n) if x not in d['bad_vertices']]):
            c3v_good[d['c3_good'][i]] += 1

    for c3v in sorted(set(list(c3v_good.keys()) + list(c3v_bad.keys()))):
        g = c3v_good.get(c3v, 0)
        b = c3v_bad.get(c3v, 0)
        tot = g + b
        print(f"    c3(v)={c3v}: good={g}, bad={b}, good_rate={g/tot:.3f}")


# ============================================================
# Part 2: Check bound at n=7 (sampled)
# ============================================================
print(f"\n{'=' * 70}")
print("SUM_v b1(T\\v) BOUND at n=7,8 (sampled)")
print("=" * 70)

for n in [7, 8, 9, 10]:
    random.seed(42)
    trials = {7: 1000, 8: 500, 9: 200, 10: 100}[n]
    max_sum = 0
    sum_dist = Counter()
    tested = 0
    t0 = time.time()

    for trial in range(trials):
        A = random_tournament(n)
        b1_T = compute_b1(A, n)
        if b1_T > 0:
            continue
        tested += 1

        s = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            b1tv = compute_b1(B_sub, n - 1)
            s += b1tv

        sum_dist[s] += 1
        max_sum = max(max_sum, s)

        if (trial + 1) % max(1, trials // 4) == 0:
            print(f"  n={n}: trial {trial+1}/{trials} ({time.time()-t0:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nn={n} ({elapsed:.0f}s): {tested} b1=0 tournaments tested")
    print(f"  Max Sum_v b1(T\\v): {max_sum}")
    print(f"  n-2 = {n-2}, floor(n/2) = {n//2}")
    print(f"  Distribution:")
    for s in sorted(sum_dist.keys()):
        print(f"    Sum={s}: {sum_dist[s]} ({sum_dist[s]/tested*100:.1f}%)")


# ============================================================
# Part 3: Relationship between Sum and c3
# ============================================================
print(f"\n{'=' * 70}")
print("SUM_v b1(T\\v) vs c3 (3-CYCLE COUNT)")
print("=" * 70)

for n in [5, 6]:
    total = 2 ** (n * (n - 1) // 2)
    c3_sum_pairs = []

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

        c3 = count_3cycles(A, n)
        s = 0
        for v in range(n):
            others = [x for x in range(n) if x != v]
            B_sub, _ = get_induced(A, n, others)
            s += compute_b1(B_sub, n - 1)

        c3_sum_pairs.append((c3, s))

    print(f"\nn={n}:")
    c3_to_sums = {}
    for c3, s in c3_sum_pairs:
        if c3 not in c3_to_sums:
            c3_to_sums[c3] = Counter()
        c3_to_sums[c3][s] += 1

    print(f"  c3 | Sum distribution")
    for c3 in sorted(c3_to_sums.keys()):
        vals = c3_to_sums[c3]
        parts = [f"{s}:{vals[s]}" for s in sorted(vals.keys())]
        print(f"  {c3:3d} | {', '.join(parts)}")


# ============================================================
# Part 4: Key structural question - what makes v bad?
# ============================================================
print(f"\n{'=' * 70}")
print("STRUCTURAL ANALYSIS OF BAD VERTICES at n=6")
print("=" * 70)

n = 6
total = 2 ** (n * (n - 1) // 2)
bad_v_patterns = Counter()

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

    # Find bad vertices
    for v in range(n):
        others = [x for x in range(n) if x != v]
        B_sub, _ = get_induced(A, n, others)
        b1tv = compute_b1(B_sub, n - 1)
        if b1tv == 1:
            # Characterize: degree of v, score of T\v, c3 of T\v
            d_out_v = sum(A[v])
            # Score of T\v
            scores_tv = tuple(sorted(sum(B_sub[i]) for i in range(n - 1)))
            # c3 of T\v
            c3_tv = count_3cycles(B_sub, n - 1)
            bad_v_patterns[(d_out_v, scores_tv, c3_tv)] += 1

print(f"Bad vertex patterns (d_out(v), score(T\\v), c3(T\\v)):")
for pattern, count in sorted(bad_v_patterns.items(), key=lambda x: -x[1])[:20]:
    print(f"  {pattern}: {count}")

# Check: is b1(T\v) = 1 iff T\v has "many" 3-cycles?
print(f"\n  b1(T\\v)=1 always correlates with c3(T\\v) >= ?")
c3_values = [p[2] for p in bad_v_patterns.keys()]
print(f"  c3(T\\v) range for bad: {min(c3_values)} to {max(c3_values)}")


print("\n\nDone.")
