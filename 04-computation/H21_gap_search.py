#!/usr/bin/env python3
"""
H21_gap_search.py -- Determine whether H=21 is achievable for any tournament.

H(T) = number of directed Hamiltonian paths in tournament T.
By OCF: H(T) = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...

Strategy:
1. Exhaustive at n<=6 (fast), large sampling at n=7
2. Random + targeted sampling at n=8..12
3. Systematic edge-flip from transitive at n=8
4. OCF decomposition analysis at n=7
"""

import random
import sys
import time
from itertools import combinations
from collections import Counter, defaultdict

# ============================================================
# Core: count_H via bitmask DP
# ============================================================

def count_H_adj(A, n):
    """Count Hamiltonian paths. A is list-of-lists adjacency."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            row = A[v]
            for u in range(n):
                if (mask & (1 << u)) or not row[u]:
                    continue
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_H_bits(n, bits):
    """Count H for tournament encoded as upper-triangle bit integer."""
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return count_H_adj(A, n)


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A


def flip_edge(A, n, i, j):
    B = [row[:] for row in A]
    B[i][j], B[j][i] = B[j][i], B[i][j]
    return B


# ============================================================
# OCF computation for analysis
# ============================================================

def find_odd_cycle_vsets(A, n):
    """Find all directed odd cycle VERTEX SETS."""
    cycle_sets = []
    for size in range(3, n+1, 2):
        for subset in combinations(range(n), size):
            s = size
            sub = list(subset)
            idx = {v: i for i, v in enumerate(sub)}
            sub_A = [[0]*s for _ in range(s)]
            for a in sub:
                for b in sub:
                    if a != b and A[a][b]:
                        sub_A[idx[a]][idx[b]] = 1
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << s):
                if not (mask & 1):
                    continue
                for v in range(s):
                    c = dp.get((mask, v), 0)
                    if not (mask & (1 << v)) or c == 0:
                        continue
                    for u in range(s):
                        if mask & (1 << u):
                            continue
                        if sub_A[v][u]:
                            key = (mask | (1 << u), u)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << s) - 1
            has_cycle = any(dp.get((full, v), 0) > 0 and sub_A[v][0]
                           for v in range(1, s))
            if has_cycle:
                cycle_sets.append(frozenset(subset))
    return cycle_sets


def compute_IP_at_2(cycle_sets):
    """Compute I(Omega, 2) from list of cycle vertex-sets."""
    m = len(cycle_sets)
    if m == 0:
        return 1, [1]
    if m > 22:
        return None, None

    adj_bits = [0] * m
    for a in range(m):
        for b in range(a+1, m):
            if cycle_sets[a] & cycle_sets[b]:
                adj_bits[a] |= 1 << b
                adj_bits[b] |= 1 << a

    size_counts = Counter()
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if adj_bits[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            size_counts[bin(mask).count('1')] += 1

    max_k = max(size_counts.keys()) if size_counts else 0
    i_vals = [size_counts.get(k, 0) for k in range(max_k + 1)]
    I_at_2 = sum(i_vals[k] * (2**k) for k in range(len(i_vals)))
    return I_at_2, i_vals


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print("=" * 70)
    print("H=21 GAP SEARCH -- Is H(T)=21 achievable for any tournament T?")
    print("=" * 70)

    # --- Theoretical decompositions ---
    print("\n--- Theoretical: all (i_1, i_2, ...) giving H=21 ---")
    print("H = 1 + 2*i_1 + 4*i_2 + 8*i_3 + 16*i_4")
    print("Need: 2*i_1 + 4*i_2 + 8*i_3 + 16*i_4 = 20\n")

    solutions = []
    for i1 in range(11):
        r1 = 20 - 2*i1
        if r1 < 0: break
        for i2 in range(r1//4 + 1):
            r2 = r1 - 4*i2
            if r2 < 0: break
            for i3 in range(r2//8 + 1):
                r3 = r2 - 8*i3
                if r3 < 0: break
                if r3 == 0:
                    solutions.append((i1, i2, i3))
                elif r3 % 16 == 0:
                    solutions.append((i1, i2, i3, r3 // 16))

    for s in solutions:
        parts = ', '.join(f'i_{k+1}={s[k]}' for k in range(len(s)) if s[k] > 0)
        if not parts:
            parts = "(impossible)"
        print(f"  {parts}")
    print(f"  Total: {len(solutions)} decompositions")
    sys.stdout.flush()

    # ============================================================
    # PHASE 1: Exhaustive n=3..6
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 1: Exhaustive search (n=3 through 6)")
    print("=" * 70)

    all_achieved = set()
    for n in range(3, 7):
        m = n * (n - 1) // 2
        total = 1 << m
        h_counts = Counter()
        t0 = time.time()
        for bits in range(total):
            h_counts[count_H_bits(n, bits)] += 1
        elapsed = time.time() - t0

        all_achieved.update(h_counts.keys())
        vals = sorted(h_counts.keys())
        print(f"\n  n={n} ({total} tours, {elapsed:.1f}s): H values = {vals}")
        print(f"    Near 21: " + " ".join(
            f"H={hv}:{h_counts.get(hv,0)}" for hv in range(15, 28, 2)))
        if 21 in h_counts:
            print(f"    *** H=21 FOUND! Count={h_counts[21]} ***")
        sys.stdout.flush()

    # ============================================================
    # PHASE 2: Large sampling n=7..12
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 2: Random sampling (n=7 through 12)")
    print("=" * 70)

    sample_sizes = [(7, 100000), (8, 50000), (9, 20000),
                    (10, 5000), (11, 1000), (12, 300)]

    for n, num_samples in sample_sizes:
        h_counts = Counter()
        found_21 = None
        t0 = time.time()

        for trial in range(num_samples):
            A = random_tournament(n)
            h = count_H_adj(A, n)
            h_counts[h] += 1
            if h == 21 and found_21 is None:
                found_21 = [row[:] for row in A]

            if (trial + 1) % 20000 == 0:
                elapsed = time.time() - t0
                print(f"    n={n}: {trial+1}/{num_samples} ({elapsed:.0f}s)")
                sys.stdout.flush()

        elapsed = time.time() - t0
        all_achieved.update(h_counts.keys())

        small = sorted(h for h in h_counts if h <= 50)
        print(f"\n  n={n} ({num_samples} samples, {elapsed:.1f}s): H<=50 = {small}")
        print(f"    Near 21: " + " ".join(
            f"H={hv}:{h_counts.get(hv,0)}" for hv in range(15, 28, 2)))

        if found_21:
            print(f"    *** H=21 FOUND! ***")
            print(f"    Adjacency: {found_21}")
        else:
            print(f"    H=21 NOT found")
        sys.stdout.flush()

    missing_global = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"\n  Missing odd H <= 51 after random search: {missing_global}")

    # ============================================================
    # PHASE 3: Systematic edge-flip at n=8 (up to 3 flips)
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 3: Systematic edge-flip from transitive at n=8")
    print("=" * 70)

    n = 8
    A_trans = transitive_tournament(n)
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    ne = len(edges)
    flip_achieved = set()
    flip_achieved.add(1)

    t0 = time.time()

    # 1-flip
    for e1 in range(ne):
        A = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
        flip_achieved.add(count_H_adj(A, n))
    print(f"  1-flip ({ne}): H<=30 = "
          f"{sorted(h for h in flip_achieved if h <= 30)}")
    sys.stdout.flush()

    # 2-flip
    cnt = 0
    for e1 in range(ne):
        A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
        for e2 in range(e1+1, ne):
            A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
            flip_achieved.add(count_H_adj(A2, n))
            cnt += 1
    print(f"  2-flip ({cnt}): H<=30 = "
          f"{sorted(h for h in flip_achieved if h <= 30)}")
    sys.stdout.flush()

    # 3-flip
    cnt = 0
    for e1 in range(ne):
        A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
        for e2 in range(e1+1, ne):
            A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
            for e3 in range(e2+1, ne):
                A3 = flip_edge(A2, n, edges[e3][0], edges[e3][1])
                flip_achieved.add(count_H_adj(A3, n))
                cnt += 1
    print(f"  3-flip ({cnt}): H<=30 = "
          f"{sorted(h for h in flip_achieved if h <= 30)}")
    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")
    sys.stdout.flush()

    all_achieved.update(flip_achieved)
    missing_global = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"  Missing odd H <= 51 after flips: {missing_global}")

    if 21 in flip_achieved:
        print(f"  H=21 ACHIEVED via edge flips!")

    # ============================================================
    # PHASE 4: OCF analysis for H values near 21 at n=7
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 4: OCF decomposition for H near 21 at n=7")
    print("=" * 70)

    near_h = defaultdict(list)
    random.seed(42)
    for trial in range(20000):
        A = random_tournament(7)
        h = count_H_adj(A, 7)
        if 15 <= h <= 27 and len(near_h[h]) < 2:
            near_h[h].append(A)

    print(f"  H values found near 21: {sorted(near_h.keys())}")

    for hv in sorted(near_h.keys()):
        A = near_h[hv][0]
        cycle_sets = find_odd_cycle_vsets(A, 7)
        I_at_2, i_vals = compute_IP_at_2(cycle_sets)
        c3 = sum(1 for cs in cycle_sets if len(cs) == 3)
        c5 = sum(1 for cs in cycle_sets if len(cs) == 5)
        c7 = sum(1 for cs in cycle_sets if len(cs) == 7)
        i_str = ', '.join(f'i{k}={i_vals[k]}' for k in range(len(i_vals)))
        print(f"    H={hv}: c3={c3} c5={c5} c7={c7} | {i_str} "
              f"| I(Omega,2)={I_at_2}")
    sys.stdout.flush()

    # ============================================================
    # PHASE 5: Targeted search
    # ============================================================
    print("\n" + "=" * 70)
    print("PHASE 5: Targeted hill-climbing for H=21")
    print("=" * 70)

    for n in [8, 9, 10]:
        print(f"\n  Targeted at n={n}:")
        best_diff = 100
        best_h = None
        t0 = time.time()

        for attempt in range(500):
            A = transitive_tournament(n)
            h = 1
            for step in range(n * (n-1) // 4):
                i = random.randint(0, n-1)
                j = random.randint(0, n-1)
                if i == j:
                    continue
                if i > j:
                    i, j = j, i
                B = flip_edge(A, n, i, j)
                new_h = count_H_adj(B, n)
                if abs(new_h - 21) < abs(h - 21):
                    A, h = B, new_h
                elif abs(new_h - 21) == abs(h - 21) and random.random() < 0.5:
                    A, h = B, new_h
                elif random.random() < 0.05:
                    A, h = B, new_h
                if h == 21:
                    break
            if h == 21:
                print(f"    *** H=21 FOUND at attempt {attempt}! ***")
                print(f"    A = {A}")
                all_achieved.add(21)
                break
            if abs(h - 21) < best_diff:
                best_diff = abs(h - 21)
                best_h = h

        elapsed = time.time() - t0
        if 21 not in all_achieved:
            print(f"    Not found. Closest: H={best_h} (diff={best_diff}). "
                  f"Time: {elapsed:.1f}s")
        sys.stdout.flush()

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)

    final_missing = [h for h in range(1, 52, 2) if h not in all_achieved]
    print(f"Missing odd H values <= 51: {final_missing}")

    print(f"\nH=7 status: "
          f"{'ACHIEVED' if 7 in all_achieved else 'NOT FOUND (proved impossible for all n)'}")
    print(f"H=21 status: "
          f"{'ACHIEVED' if 21 in all_achieved else 'NOT FOUND in any search'}")

    if 21 not in all_achieved:
        print("\nEvidence for H=21 impossibility:")
        print("  - Exhaustive at n=3..6: NOT achieved")
        print("  - 100K random at n=7: NOT achieved")
        print("  - 50K random at n=8: NOT achieved")
        print("  - 20K random at n=9: NOT achieved")
        print("  - Smaller samples at n=10,11,12: NOT achieved")
        print("  - Systematic 3-flip at n=8: NOT achieved")
        print("  - Targeted hill-climbing at n=8,9,10: NOT achieved")
        print("  - H=7 is proved impossible for all n (THM-029)")
        print("  - H=21 may be similarly impossible via OCF constraints")
