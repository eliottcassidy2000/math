#!/usr/bin/env python3
"""
root_profile_determinacy_s18.py — kind-pasteur-2026-03-21-S18

KEY DISCOVERY FROM petersen_lie_bridge: at n=5, the sorted root cycle
profile uniquely determines H(T). Test at n=6 (exhaustive) and n=7 (sampled).

Root cycle profile: for each arc (positive root of A_{n-1}), count the
number of directed 3-cycles through that arc. Sort the resulting vector.

Also test: does the FULL root cycle profile (including 5-cycle counts,
7-cycle counts, etc.) do better?

And: joint (weight norm, Casimir, root profile) classification.
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
import random

sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def directed_3cycle_count_per_arc(A, n):
    """For each arc (i,j) with i<j, count directed 3-cycles through it."""
    profile = []
    for i in range(n):
        for j in range(i+1, n):
            c = 0
            for k in range(n):
                if k == i or k == j:
                    continue
                # 3-cycle using arc i->j (if it exists): need j->k->i
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                # 3-cycle using arc j->i (if it exists): need i->k->j
                if A[j][i] and A[i][k] and A[k][j]:
                    c += 1
            profile.append(c)
    return tuple(sorted(profile))

def directed_all_cycle_profile(A, n):
    """For each arc, count directed odd cycles of ALL lengths through it."""
    # For efficiency, only do 3-cycles and 5-cycles
    profile_3 = []
    profile_5 = []

    for i in range(n):
        for j in range(i+1, n):
            c3 = 0
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[i][j] and A[j][k] and A[k][i]:
                    c3 += 1
                if A[j][i] and A[i][k] and A[k][j]:
                    c3 += 1
            profile_3.append(c3)

            # 5-cycles through arc (i,j)
            if n >= 5:
                c5 = 0
                for subset in combinations(range(n), 5):
                    if i not in subset or j not in subset:
                        continue
                    # Check if there's a directed 5-cycle on this subset
                    sub = list(subset)
                    for perm in permutations(sub):
                        is_cycle = True
                        for idx in range(5):
                            if not A[perm[idx]][perm[(idx+1) % 5]]:
                                is_cycle = False
                                break
                        if is_cycle:
                            c5 += 1
                # Each 5-cycle is counted 5 times (once per starting vertex)
                # and each arc appears in the cycle, so divide by 5
                c5 = c5 // 5
                profile_5.append(c5)

    return tuple(sorted(profile_3)), tuple(sorted(profile_5)) if n >= 5 else ()

print("=" * 72)
print("  ROOT CYCLE PROFILE DETERMINACY TEST")
print("  kind-pasteur-2026-03-21-S18")
print("=" * 72)

# ========================================================================
# n=5: CONFIRM (should be True)
# ========================================================================
print("\n--- n=5 (exhaustive, confirmation) ---")
n = 5
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

profile_to_H_5 = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    H = count_hp(A, n)
    prof = directed_3cycle_count_per_arc(A, n)
    profile_to_H_5[prof].add(H)

determines_5 = all(len(v) == 1 for v in profile_to_H_5.values())
print(f"  Profiles: {len(profile_to_H_5)}, determines H: {determines_5}")
for prof, Hs in sorted(profile_to_H_5.items()):
    print(f"    {list(prof)} -> H = {sorted(Hs)}")

# ========================================================================
# n=6: EXHAUSTIVE
# ========================================================================
print("\n--- n=6 (exhaustive) ---")
n = 6
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

profile_to_H_6 = defaultdict(set)
profile_count_6 = defaultdict(int)
total_6 = 0
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    H = count_hp(A, n)
    prof = directed_3cycle_count_per_arc(A, n)
    profile_to_H_6[prof].add(H)
    profile_count_6[prof] += 1
    total_6 += 1

    if total_6 % 100000 == 0:
        print(f"  ... processed {total_6}/{2**m} ({100*total_6/2**m:.1f}%)")

determines_6 = all(len(v) == 1 for v in profile_to_H_6.values())
print(f"  Profiles: {len(profile_to_H_6)}, determines H: {determines_6}")

if not determines_6:
    print(f"\n  COUNTEREXAMPLES (profile maps to multiple H values):")
    failures = [(prof, Hs) for prof, Hs in profile_to_H_6.items() if len(Hs) > 1]
    failures.sort(key=lambda x: -len(x[1]))
    for prof, Hs in failures[:20]:
        count = profile_count_6[prof]
        print(f"    {list(prof)} -> H in {sorted(Hs)} (count={count})")
    print(f"  Total failures: {len(failures)} / {len(profile_to_H_6)} profiles")
else:
    # Show all profiles
    for prof, Hs in sorted(profile_to_H_6.items()):
        count = profile_count_6[prof]
        print(f"    {list(prof)} -> H = {sorted(Hs)} (count={count})")

# ========================================================================
# n=6: JOINT (3-cycle profile, weight norm) CLASSIFICATION
# ========================================================================
print("\n--- n=6: Joint (profile, ||w||^2) classification ---")
n = 6

C_A5 = np.array([
    [ 2, -1,  0,  0,  0],
    [-1,  2, -1,  0,  0],
    [ 0, -1,  2, -1,  0],
    [ 0,  0, -1,  2, -1],
    [ 0,  0,  0, -1,  2]
])

joint_to_H_6 = defaultdict(set)
for bits in range(2**m):
    A = np.zeros((n, n), dtype=int)
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    H = count_hp(A, n)
    prof = directed_3cycle_count_per_arc(A, n)

    # Weight norm
    w = np.zeros(n-1, dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            eps = 1 if A[i][j] else -1
            for k_idx in range(i, j):
                w[k_idx] += eps
    wt_norm = int(w @ C_A5 @ w)

    joint_to_H_6[(prof, wt_norm)].add(H)

determines_joint_6 = all(len(v) == 1 for v in joint_to_H_6.values())
print(f"  Joint classes: {len(joint_to_H_6)}, determines H: {determines_joint_6}")

if not determines_joint_6:
    failures_j = [(k, v) for k, v in joint_to_H_6.items() if len(v) > 1]
    print(f"  Remaining failures: {len(failures_j)} / {len(joint_to_H_6)}")
    for (prof, wn), Hs in sorted(failures_j, key=lambda x: -len(x[1]))[:10]:
        print(f"    profile sum={sum(prof)}, ||w||^2={wn}: H in {sorted(Hs)}")

# ========================================================================
# n=7: SAMPLED
# ========================================================================
print("\n--- n=7 (sampled, 50K random tournaments) ---")
n = 7
random.seed(42)

profile_to_H_7 = defaultdict(set)
for trial in range(50000):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    H = count_hp(A, n)
    prof = directed_3cycle_count_per_arc(A, n)
    profile_to_H_7[prof].add(H)

    if (trial + 1) % 10000 == 0:
        fails = sum(1 for v in profile_to_H_7.values() if len(v) > 1)
        print(f"  ... {trial+1}: {len(profile_to_H_7)} profiles, {fails} multi-H")

determines_7 = all(len(v) == 1 for v in profile_to_H_7.values())
print(f"  Profiles: {len(profile_to_H_7)}, determines H: {determines_7}")

if not determines_7:
    failures_7 = [(prof, Hs) for prof, Hs in profile_to_H_7.items() if len(Hs) > 1]
    print(f"  Failures: {len(failures_7)} / {len(profile_to_H_7)}")
    for prof, Hs in sorted(failures_7, key=lambda x: -len(x[1]))[:15]:
        print(f"    sum={sum(prof)} {list(prof)} -> H in {sorted(Hs)}")

print("\n" + "=" * 72)
print("  SUMMARY")
print("=" * 72)
