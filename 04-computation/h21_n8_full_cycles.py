#!/usr/bin/env python3
"""
At n=8: check if alpha_1=10 with i_2=0 (using 3+5 cycles only) changes
when we add 7-cycles to the count.

The previous script found alpha_1(3+5)=10, i_2=0 at n=8.
But the full alpha_1 includes 7-cycles, which would increase alpha_1
and potentially add new disjoint pairs.

Key question: when alpha_1(3+5)=10 and i_2(3+5)=0, does adding 7-cycles
push alpha_1 above 10 (making H > 21) or add disjoint pairs?

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter
import random
import time

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def find_5cycles_dp(adj, n):
    cycles = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 32):
            for i in range(5):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(5):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 5):
            if (31, j) in dp and adj[v[j]][v[0]]:
                count += dp[(31, j)]
        num = count // 5
        for _ in range(num):
            cycles.append(frozenset(verts))
    return cycles

def find_7cycles_dp(adj, n):
    cycles = []
    for verts in combinations(range(n), 7):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 128):
            for i in range(7):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(7):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 7):
            if (127, j) in dp and adj[v[j]][v[0]]:
                count += dp[(127, j)]
        num = count // 7
        for _ in range(num):
            cycles.append(frozenset(verts))
    return cycles

def count_i2(all_cycles):
    alpha1 = len(all_cycles)
    i2 = 0
    for a in range(alpha1):
        for b in range(a+1, alpha1):
            if not (all_cycles[a] & all_cycles[b]):
                i2 += 1
    return i2

def held_karp(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[S | (1 << u)][u] += dp[S][v]
    return sum(dp[full])


def main():
    n = 8
    print(f"n={n}: Checking if alpha_1(3+5)=10, i_2=0 survives adding 7-cycles")
    print(f"Also checking alpha_1(3+5)=8, i_2=1")

    random.seed(42)
    found_10_0 = []
    found_8_1 = []

    t0 = time.time()
    checked = 0
    for trial in range(100000):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        c3 = find_3cycles(adj, n)
        if len(c3) > 12:
            continue

        c5 = find_5cycles_dp(adj, n)
        partial = c3 + c5
        alpha1_35 = len(partial)

        if alpha1_35 == 10:
            i2_35 = count_i2(partial)
            if i2_35 == 0:
                # NOW add 7-cycles
                c7 = find_7cycles_dp(adj, n)
                full = c3 + c5 + c7
                alpha1_full = len(full)
                i2_full = count_i2(full)
                H = held_karp(adj, n)
                found_10_0.append({
                    't3': len(c3), 't5': len(c5), 't7': len(c7),
                    'alpha1_35': alpha1_35, 'i2_35': i2_35,
                    'alpha1_full': alpha1_full, 'i2_full': i2_full,
                    'H': H,
                })
                if len(found_10_0) <= 5:
                    print(f"  10,0 candidate: t3={len(c3)},t5={len(c5)},t7={len(c7)}, "
                          f"alpha_full={alpha1_full}, i2_full={i2_full}, H={H}")

        if alpha1_35 == 8:
            i2_35 = count_i2(partial)
            if i2_35 == 1:
                c7 = find_7cycles_dp(adj, n)
                full = c3 + c5 + c7
                alpha1_full = len(full)
                i2_full = count_i2(full)
                H = held_karp(adj, n)
                found_8_1.append({
                    't3': len(c3), 't5': len(c5), 't7': len(c7),
                    'alpha1_35': alpha1_35, 'i2_35': i2_35,
                    'alpha1_full': alpha1_full, 'i2_full': i2_full,
                    'H': H,
                })
                if len(found_8_1) <= 5:
                    print(f"  8,1 candidate: t3={len(c3)},t5={len(c5)},t7={len(c7)}, "
                          f"alpha_full={alpha1_full}, i2_full={i2_full}, H={H}")

        checked += 1

    elapsed = time.time() - t0
    print(f"\nChecked {checked} in {elapsed:.1f}s")

    print(f"\nalpha_1(3+5)=10, i2(3+5)=0: {len(found_10_0)} found")
    if found_10_0:
        # Do the 7-cycles push alpha_1 > 10?
        alpha_full_dist = Counter(d['alpha1_full'] for d in found_10_0)
        i2_full_dist = Counter(d['i2_full'] for d in found_10_0)
        H_dist = Counter(d['H'] for d in found_10_0)
        t7_dist = Counter(d['t7'] for d in found_10_0)
        print(f"  alpha1_full: {dict(sorted(alpha_full_dist.items()))}")
        print(f"  i2_full: {dict(sorted(i2_full_dist.items()))}")
        print(f"  t7: {dict(sorted(t7_dist.items()))}")
        print(f"  H: {dict(sorted(H_dist.items()))}")

        h21_count = sum(1 for d in found_10_0 if d['H'] == 21)
        print(f"  H=21 count: {h21_count}")

    print(f"\nalpha_1(3+5)=8, i2(3+5)=1: {len(found_8_1)} found")
    if found_8_1:
        alpha_full_dist = Counter(d['alpha1_full'] for d in found_8_1)
        i2_full_dist = Counter(d['i2_full'] for d in found_8_1)
        H_dist = Counter(d['H'] for d in found_8_1)
        t7_dist = Counter(d['t7'] for d in found_8_1)
        print(f"  alpha1_full: {dict(sorted(alpha_full_dist.items()))}")
        print(f"  i2_full: {dict(sorted(i2_full_dist.items()))}")
        print(f"  t7: {dict(sorted(t7_dist.items()))}")
        print(f"  H: {dict(sorted(H_dist.items()))}")

        h21_count = sum(1 for d in found_8_1 if d['H'] == 21)
        print(f"  H=21 count: {h21_count}")


if __name__ == "__main__":
    main()
