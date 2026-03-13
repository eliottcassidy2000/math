#!/usr/bin/env python3
"""
trh_d2_t3_test.py — Test whether d²≠0 ⟺ t3=2 at n=5, and explore n=6,7.

Key discovery from trh_d2_structure.py:
At n=5, d²=0 FAILS for EXACTLY the tournaments with t3=2 (240/1024).
This is a PERFECT characterization with 0 exceptions.

Questions:
1. Why does t3=2 specifically break d²=0?
2. What is the analogous condition at n=6?
3. What specific t3 values cause failure at n=6?
"""

import numpy as np
from itertools import combinations
from collections import Counter

def get_regular_paths(A, m):
    n = A.shape[0]
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

def interior_boundary_matrix(paths_m, paths_m1):
    if not paths_m or not paths_m1:
        return np.zeros((len(paths_m1) if paths_m1 else 0,
                         len(paths_m) if paths_m else 0), dtype=int)
    idx = {p: i for i, p in enumerate(paths_m1)}
    m = len(paths_m[0]) - 1
    B = np.zeros((len(paths_m1), len(paths_m)), dtype=int)
    for j, path in enumerate(paths_m):
        for i in range(1, m):
            face = path[:i] + path[i+1:]
            if face in idx:
                B[idx[face], j] += (-1)**i
    return B

def check_d2(A):
    n = A.shape[0]
    all_paths = {m: get_regular_paths(A, m) for m in range(n)}
    for m in range(2, n):
        p_m = all_paths[m]
        p_m1 = all_paths[m-1]
        p_m2 = all_paths[m-2]
        if p_m and p_m1 and p_m2:
            B_m = interior_boundary_matrix(p_m, p_m1)
            B_m1 = interior_boundary_matrix(p_m1, p_m2)
            d2 = B_m1 @ B_m
            if np.max(np.abs(d2)) > 0:
                return False
    return True

def check_d2_by_degree(A):
    """Return dict: degree -> whether d²=0 at that degree."""
    n = A.shape[0]
    all_paths = {m: get_regular_paths(A, m) for m in range(n)}
    results = {}
    for m in range(2, n):
        p_m = all_paths[m]
        p_m1 = all_paths[m-1]
        p_m2 = all_paths[m-2]
        if p_m and p_m1 and p_m2:
            B_m = interior_boundary_matrix(p_m, p_m1)
            B_m1 = interior_boundary_matrix(p_m1, p_m2)
            d2 = B_m1 @ B_m
            results[m] = np.max(np.abs(d2)) == 0
        else:
            results[m] = True
    return results

def score_sequence(A):
    n = A.shape[0]
    return tuple(sorted([int(sum(A[v])) for v in range(n)]))

def count_3cycles(A):
    n = A.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[i][k] and A[k][j]):
                    count += 1
    return count

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    for bits in range(2**len(pairs)):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ============================================================
print("=" * 70)
print("VERIFICATION: n=5, d²≠0 ⟺ t3=2")
print("=" * 70)

n = 5
t3_vs_d2 = {}  # t3 -> (pass_count, fail_count)
for A in all_tournaments(n):
    t3 = count_3cycles(A)
    d2_ok = check_d2(A)
    if t3 not in t3_vs_d2:
        t3_vs_d2[t3] = [0, 0]
    t3_vs_d2[t3][0 if d2_ok else 1] += 1

print("\nt3 | d²=0 | d²≠0")
print("-" * 30)
for t3 in sorted(t3_vs_d2):
    p, f = t3_vs_d2[t3]
    print(f"  {t3} | {p:4d} | {f:4d}  {'*** FAIL' if f > 0 else ''}")

print(f"\nPERFECT: d²≠0 ⟺ t3=2 at n=5: {t3_vs_d2[2][0] == 0 and t3_vs_d2[2][1] == 240}")

# ============================================================
print(f"\n{'=' * 70}")
print("n=6: COMPLETE ENUMERATION — t3 vs d² (may take a few minutes)")
print("=" * 70)

n = 6
t3_vs_d2_6 = {}
fail_by_degree_6 = Counter()
total_6 = 0
for A in all_tournaments(n):
    t3 = count_3cycles(A)
    d2_ok = check_d2(A)
    if t3 not in t3_vs_d2_6:
        t3_vs_d2_6[t3] = [0, 0]
    t3_vs_d2_6[t3][0 if d2_ok else 1] += 1

    if not d2_ok:
        deg_results = check_d2_by_degree(A)
        for m in sorted(deg_results):
            if not deg_results[m]:
                fail_by_degree_6[m] += 1
                break

    total_6 += 1
    if total_6 % 5000 == 0:
        print(f"  processed {total_6}/32768...")

print(f"\nTotal n=6 tournaments: {total_6}")
print("\nt3 | d²=0 | d²≠0 | fail%")
print("-" * 45)
for t3 in sorted(t3_vs_d2_6):
    p, f = t3_vs_d2_6[t3]
    pct = f / (p + f) * 100
    marker = ""
    if f > 0 and p > 0:
        marker = " (MIXED)"
    elif f > 0 and p == 0:
        marker = " (ALL FAIL)"
    print(f"  {t3:2d} | {p:5d} | {f:5d} | {pct:5.1f}%{marker}")

total_pass = sum(v[0] for v in t3_vs_d2_6.values())
total_fail = sum(v[1] for v in t3_vs_d2_6.values())
print(f"\nTotal: {total_pass} pass, {total_fail} fail ({total_fail/total_6*100:.1f}%)")
print(f"First failing degree distribution: {dict(fail_by_degree_6)}")

# Check if there's a clean t3-based characterization
pure_fail_t3 = [t3 for t3, (p, f) in t3_vs_d2_6.items() if p == 0 and f > 0]
pure_pass_t3 = [t3 for t3, (p, f) in t3_vs_d2_6.items() if f == 0 and p > 0]
mixed_t3 = [t3 for t3, (p, f) in t3_vs_d2_6.items() if p > 0 and f > 0]

print(f"\nPure FAIL t3 values: {pure_fail_t3}")
print(f"Pure PASS t3 values: {pure_pass_t3}")
print(f"Mixed t3 values: {mixed_t3}")

if mixed_t3:
    print("\n*** t3 alone does NOT characterize d²=0 at n=6 ***")
    print("Need additional invariant. Checking score sequence + t3...")

    # Check if (score_seq, t3) is sufficient
    ss_t3_vs_d2 = {}
    for A in all_tournaments(n):
        ss = score_sequence(A)
        t3 = count_3cycles(A)
        key = (ss, t3)
        d2_ok = check_d2(A)
        if key not in ss_t3_vs_d2:
            ss_t3_vs_d2[key] = [0, 0]
        ss_t3_vs_d2[key][0 if d2_ok else 1] += 1

    mixed_ss_t3 = [(k, v) for k, v in ss_t3_vs_d2.items() if v[0] > 0 and v[1] > 0]
    if mixed_ss_t3:
        print(f"  (score_seq, t3) has {len(mixed_ss_t3)} mixed classes:")
        for (ss, t3), (p, f) in sorted(mixed_ss_t3):
            print(f"    ss={ss}, t3={t3}: {p} pass, {f} fail")
    else:
        print(f"  (score_seq, t3) PERFECTLY characterizes d²=0 at n=6!")

# ============================================================
print(f"\n{'=' * 70}")
print("WHY t3=2 AT n=5?")
print("=" * 70)

print("""
At n=5 with 5 vertices, C(5,3)=10 triples.
Each triple is either a 3-cycle or a transitive triple.
t3=2 means exactly 2 of the 10 triples form 3-cycles.

The two iso classes with t3=2 are:
  ss=(0,2,2,3,3): "near-transitive with two cycles"
  ss=(1,1,2,2,4): "dominant vertex + two cycles"

Both have a specific structural feature: the two 3-cycles share
exactly one edge (they share 2 vertices), creating a "butterfly" pattern.
""")

# Verify: do the two 3-cycles always share 2 vertices?
print("Checking: do the two 3-cycles share vertices at t3=2?")
for A in all_tournaments(5):
    if count_3cycles(A) != 2:
        continue
    cycles = []
    for i in range(5):
        for j in range(i+1, 5):
            for k in range(j+1, 5):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[i][k] and A[k][j]):
                    cycles.append({i, j, k})
    overlap = len(cycles[0] & cycles[1])
    print(f"  ss={score_sequence(A)}: cycles {cycles[0]} & {cycles[1]}, overlap={overlap}")
    break  # just one example per class, they're isomorphic

# Check all
overlaps = Counter()
for A in all_tournaments(5):
    if count_3cycles(A) != 2:
        continue
    cycles = []
    for i in range(5):
        for j in range(i+1, 5):
            for k in range(j+1, 5):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[i][k] and A[k][j]):
                    cycles.append({i, j, k})
    overlap = len(cycles[0] & cycles[1])
    overlaps[overlap] += 1

print(f"\nOverlap distribution for t3=2 tournaments: {dict(overlaps)}")

# ============================================================
print(f"\n{'=' * 70}")
print("ADDITIONAL n=6 ANALYSIS: WHAT DISTINGUISHES PASS/FAIL AT MIXED t3?")
print("=" * 70)

if mixed_t3:
    for t3_val in mixed_t3[:3]:
        print(f"\n--- t3={t3_val} (mixed) ---")
        pass_count = 0
        fail_count = 0
        pass_4cycles = Counter()
        fail_4cycles = Counter()

        for A in all_tournaments(6):
            if count_3cycles(A) != t3_val:
                continue
            d2_ok = check_d2(A)
            # Count 4-cycles
            t4 = 0
            for quad in combinations(range(6), 4):
                B = np.zeros((4, 4), dtype=int)
                for ii, vi in enumerate(quad):
                    for jj, vj in enumerate(quad):
                        B[ii][jj] = A[vi][vj]
                # 4-cycle: 4 edges forming a cycle
                # Just count directed 4-cycles
                for perm in [(0,1,2,3), (0,1,3,2), (0,2,1,3), (0,2,3,1), (0,3,1,2), (0,3,2,1)]:
                    if all(B[perm[i]][perm[(i+1)%4]] for i in range(4)):
                        t4 += 1
            t4 //= 4  # each 4-cycle counted 4 times (cyclic rotations... actually n times starting points)

            if d2_ok:
                pass_count += 1
                pass_4cycles[t4] += 1
            else:
                fail_count += 1
                fail_4cycles[t4] += 1

        print(f"  pass={pass_count}, fail={fail_count}")
        if pass_count <= 20 and fail_count <= 20:
            print(f"  4-cycle counts in pass: {dict(pass_4cycles)}")
            print(f"  4-cycle counts in fail: {dict(fail_4cycles)}")

print("\nDONE.")
