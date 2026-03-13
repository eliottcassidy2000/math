#!/usr/bin/env python3
"""
p13_k7_mechanism.py -- The k=7 mechanism at p=13

At p=13, 78% of H comes from k>=7 cycles. Since max independent set
size = floor(13/3) = 4 (from four disjoint 3-cycles), and
3+3+7 = 13 = p (no room), (3,7) pairs can't be disjoint.

But wait: 5+7=12 < 13, so (5,7) disjoint pairs ARE possible!
And 7 < 13 but 7+7=14 > 13, so no (7,7) disjoint pairs.

So the FULL alpha structure at p=13 involves:
  alpha_1 from k=3,5,7,9,11,13
  alpha_2 from (3,3), (3,5), (5,5), (5,7) -- only these fit
  alpha_3 from (3,3,3), (3,3,5) -- 3+3+5=11 fits, 3+5+5=13 is tight!
  alpha_4 from (3,3,3,3) -- 12 < 13, fits!

Tight constraint: 3+5+5=13 means a (3,5,5) triple uses ALL vertices.
These are "perfect covers" of Z_13 by disjoint odd cycles.

This script computes ALL these contributions for each distinct H-value class.

Author: kind-pasteur-2026-03-12-S59c
"""

import time
from itertools import combinations
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def count_directed_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
    start = 0
    dp = {}
    dp[(1 << start, start)] = 1
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0 and A[verts[v]][verts[start]]:
            total += dp[key]
    return total


def full_alpha_decomposition(A, p, max_k=7):
    """Compute alpha_1 through alpha_4 with full cross-length terms."""
    m = (p - 1) // 2

    # Enumerate all directed cycles up to max_k
    all_dc = []  # list of (frozenset, k, dc_index_within_vertex_set)
    dc_by_k = {}

    for k in range(3, max_k + 1, 2):
        cycles_k = []
        for subset in combinations(range(p), k):
            d = count_directed_ham_cycles(A, list(subset))
            for idx in range(d):
                cycles_k.append(frozenset(subset))
        dc_by_k[k] = cycles_k
        all_dc.extend([(fs, k) for fs in cycles_k])

    n = len(all_dc)

    # alpha_1
    alpha_1 = n
    alpha_1_by_k = {k: len(dc_by_k[k]) for k in dc_by_k}

    # alpha_2: disjoint pairs
    alpha_2 = 0
    alpha_2_by_type = defaultdict(int)
    for i in range(n):
        for j in range(i + 1, n):
            if not (all_dc[i][0] & all_dc[j][0]):
                alpha_2 += 1
                kk = tuple(sorted([all_dc[i][1], all_dc[j][1]]))
                alpha_2_by_type[kk] += 1

    # alpha_3: mutually disjoint triples
    alpha_3 = 0
    alpha_3_by_type = defaultdict(int)

    # Pre-compute non-adjacent list
    non_adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if not (all_dc[i][0] & all_dc[j][0]):
                non_adj[i].append(j)

    for i in range(n):
        for j in non_adj[i]:
            for l in non_adj[j]:
                if l <= j:
                    continue
                if not (all_dc[i][0] & all_dc[l][0]):
                    alpha_3 += 1
                    ks = tuple(sorted([all_dc[i][1], all_dc[j][1], all_dc[l][1]]))
                    alpha_3_by_type[ks] += 1

    # alpha_4: mutually disjoint quadruples (only from 3-cycles at p=13)
    alpha_4 = 0
    alpha_4_by_type = defaultdict(int)

    for i in range(n):
        for j in non_adj[i]:
            for l in non_adj[j]:
                if l <= j:
                    continue
                if all_dc[i][0] & all_dc[l][0]:
                    continue
                for m_val in non_adj[l]:
                    if m_val <= l:
                        continue
                    if all_dc[i][0] & all_dc[m_val][0]:
                        continue
                    if all_dc[j][0] & all_dc[m_val][0]:
                        continue
                    alpha_4 += 1
                    ks = tuple(sorted([all_dc[i][1], all_dc[j][1],
                                      all_dc[l][1], all_dc[m_val][1]]))
                    alpha_4_by_type[ks] += 1

    return {
        'alpha_1': alpha_1,
        'alpha_1_by_k': alpha_1_by_k,
        'alpha_2': alpha_2,
        'alpha_2_by_type': dict(sorted(alpha_2_by_type.items())),
        'alpha_3': alpha_3,
        'alpha_3_by_type': dict(sorted(alpha_3_by_type.items())),
        'alpha_4': alpha_4,
        'alpha_4_by_type': dict(sorted(alpha_4_by_type.items())),
    }


def main():
    print("=" * 70)
    print("p=13 FULL ALPHA DECOMPOSITION (k<=7)")
    print("=" * 70)

    p = 13
    m = (p - 1) // 2

    # Get one representative from each H-value class
    S_int = list(range(1, m + 1))
    pairs = [(s, p - s) for s in range(1, m + 1)]
    H_classes = defaultdict(list)

    for bits in range(1 << m):
        S = []
        for i, (a, b) in enumerate(pairs):
            S.append(a if (bits & (1 << i)) else b)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        H_classes[H].append(S)

    # Process each class
    results = {}
    for H_val in sorted(H_classes, reverse=True):
        S_rep = H_classes[H_val][0]
        A = build_adj(p, S_rep)

        is_int = S_rep == S_int

        print(f"\n{'='*60}")
        print(f"H = {H_val}, S = {S_rep}{' [INTERVAL]' if is_int else ''}")
        print(f"{'='*60}")

        t0 = time.time()
        data = full_alpha_decomposition(A, p, max_k=7)
        t1 = time.time()

        print(f"\n  Computed in {t1-t0:.1f}s")

        print(f"\n  alpha_1 = {data['alpha_1']}")
        for k, cnt in sorted(data['alpha_1_by_k'].items()):
            print(f"    k={k}: {cnt}")

        print(f"\n  alpha_2 = {data['alpha_2']}")
        for kk, cnt in sorted(data['alpha_2_by_type'].items()):
            print(f"    {kk}: {cnt}")

        print(f"\n  alpha_3 = {data['alpha_3']}")
        for ks, cnt in sorted(data['alpha_3_by_type'].items()):
            print(f"    {ks}: {cnt}")

        print(f"\n  alpha_4 = {data['alpha_4']}")
        for ks, cnt in sorted(data['alpha_4_by_type'].items()):
            print(f"    {ks}: {cnt}")

        # OCF check
        H_check = 1 + 2*data['alpha_1'] + 4*data['alpha_2'] + 8*data['alpha_3'] + 16*data['alpha_4']
        remainder = H_val - H_check
        print(f"\n  H_check (k<=7, through alpha_4) = {H_check}")
        print(f"  Actual H = {H_val}")
        print(f"  Remainder from k>=9 = {remainder}")

        results[H_val] = data

    # ====== COMPARISON TABLE ======
    print(f"\n{'='*70}")
    print("COMPARISON: INTERVAL vs OTHERS")
    print("=" * 70)

    H_sorted = sorted(results.keys(), reverse=True)

    # Header
    print(f"\n  {'H':>12} {'a1':>8} {'a2':>8} {'a3':>8} {'a4':>6} "
          f"{'2a1':>10} {'4a2':>10} {'8a3':>10} {'16a4':>8} {'sum':>10} {'rem':>10}")
    print(f"  {'-'*108}")

    for H_val in H_sorted:
        d = results[H_val]
        a1, a2, a3, a4 = d['alpha_1'], d['alpha_2'], d['alpha_3'], d['alpha_4']
        s = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4
        rem = H_val - s

        print(f"  {H_val:>12} {a1:>8} {a2:>8} {a3:>8} {a4:>6} "
              f"{2*a1:>10} {4*a2:>10} {8*a3:>10} {16*a4:>8} {s:>10} {rem:>10}")

    # ====== DELTA TABLE (relative to Interval) ======
    H_int = max(H_sorted)  # Interval is maximum
    d_int = results[H_int]

    print(f"\n  DELTA relative to Interval (H={H_int}):")
    print(f"  {'H':>12} {'da1':>8} {'da2':>8} {'da3':>8} {'da4':>6} "
          f"{'eff1':>10} {'eff2':>10} {'eff3':>10} {'eff4':>8} {'dH':>10}")
    print(f"  {'-'*98}")

    for H_val in H_sorted:
        d = results[H_val]
        da1 = d['alpha_1'] - d_int['alpha_1']
        da2 = d['alpha_2'] - d_int['alpha_2']
        da3 = d['alpha_3'] - d_int['alpha_3']
        da4 = d['alpha_4'] - d_int['alpha_4']
        dH = H_val - H_int

        print(f"  {H_val:>12} {da1:>8} {da2:>8} {da3:>8} {da4:>6} "
              f"{2*da1:>10} {4*da2:>10} {8*da3:>10} {16*da4:>8} {dH:>10}")

    # ====== BY-TYPE BREAKDOWN ======
    print(f"\n  alpha_2 by type for max and min H:")
    H_max = max(H_sorted)
    H_min = min(H_sorted)
    print(f"    H={H_max} (Interval): {results[H_max]['alpha_2_by_type']}")
    print(f"    H={H_min} (minimum):  {results[H_min]['alpha_2_by_type']}")

    print(f"\n  alpha_3 by type for max and min H:")
    print(f"    H={H_max} (Interval): {results[H_max]['alpha_3_by_type']}")
    print(f"    H={H_min} (minimum):  {results[H_min]['alpha_3_by_type']}")


if __name__ == '__main__':
    main()
