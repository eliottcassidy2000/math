#!/usr/bin/env python3
"""
Test NONHAM=0 for Paley tournament T_13.

QR mod 13: 1^2=1, 2^2=4, 3^2=9, 4^2=3, 5^2=12, 6^2=10
So QR = {1,3,4,9,10,12}, QNR = {2,5,6,7,8,11}.

Check pair (0,2) where T[0,2]=0.
U has 11 vertices, 2^11 = 2048 subsets.
Larger subsets up to 12! paths — may be slow. Use early termination.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import time
from functools import lru_cache

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def count_paths_ending_at(T, verts, a):
    """Count Ham paths in T[verts] ending at a. Uses DP for efficiency."""
    verts = tuple(sorted(verts))
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == a else 0

    # DP: dp[mask][last] = # paths visiting exactly the vertices in mask, ending at last
    idx = {v: i for i, v in enumerate(verts)}
    if a not in idx:
        return 0
    a_idx = idx[a]

    dp = {}
    for i, v in enumerate(verts):
        dp[(1 << i, i)] = 1

    for step in range(1, n):
        new_dp = {}
        for (mask, last), cnt in dp.items():
            if cnt == 0: continue
            for i, v in enumerate(verts):
                if mask & (1 << i): continue
                if T.get((verts[last], v), 0) == 1:
                    key = (mask | (1 << i), i)
                    new_dp[key] = new_dp.get(key, 0) + cnt
        dp.update(new_dp)

    full_mask = (1 << n) - 1
    return dp.get((full_mask, a_idx), 0)

def count_paths_starting_at(T, verts, b):
    """Count Ham paths in T[verts] starting at b. Uses DP."""
    verts = tuple(sorted(verts))
    n = len(verts)
    if n == 1:
        return 1 if verts[0] == b else 0

    idx = {v: i for i, v in enumerate(verts)}
    if b not in idx:
        return 0
    b_idx = idx[b]

    # DP: dp[mask][last] = # paths visiting mask, starting at b, ending at last
    dp = {((1 << b_idx), b_idx): 1}

    for step in range(1, n):
        new_dp = {}
        for (mask, last), cnt in dp.items():
            if cnt == 0: continue
            for i, v in enumerate(verts):
                if mask & (1 << i): continue
                if T.get((verts[last], v), 0) == 1:
                    key = (mask | (1 << i), i)
                    new_dp[key] = new_dp.get(key, 0) + cnt
        dp.update(new_dp)

    full_mask = (1 << n) - 1
    total = sum(dp.get((full_mask, i), 0) for i in range(n))
    return total


def check_nonham(T, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    t_start = time.time()

    for mask in range(1 << len(U)):
        S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
        R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
        sign = (-1)**len(S_list)
        S_set = sorted(set(S_list) | {a})
        R_set = sorted(set(R) | {b})

        ea = count_paths_ending_at(T, S_set, a)
        bb = count_paths_starting_at(T, R_set, b)

        total += sign * ea * bb

        if (mask + 1) % 256 == 0:
            elapsed = time.time() - t_start
            print(f"    mask {mask+1}/{1<<len(U)}: total = {total} ({elapsed:.1f}s)", flush=True)

    return total


print("=" * 70)
print("n=13: NONHAM=0 for Paley T_13?")
print("=" * 70)

n = 13
QR = {1, 3, 4, 9, 10, 12}
QNR = {2, 5, 6, 7, 8, 11}
print(f"  QR = {sorted(QR)}")
print(f"  QNR = {sorted(QNR)}")

T = circulant_tournament(n, QR)
print(f"  T[0,1]={T[(0,1)]}, T[0,2]={T[(0,2)]}")

a, b = 0, 2
print(f"\n  Checking NONHAM({a},{b}) where T[{a},{b}]={T[(a,b)]}")

t0 = time.time()
nh = check_nonham(T, n, a, b)
elapsed = time.time() - t0

print(f"\n  NONHAM({a},{b}) = {nh}")
print(f"  Time: {elapsed:.1f}s")

if nh == 0:
    print("\n  CONFIRMED: NONHAM=0 for Paley T_13")
else:
    print(f"\n  COUNTEREXAMPLE!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
