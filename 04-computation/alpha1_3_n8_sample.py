#!/usr/bin/env python3
"""
Check if alpha_1=3 occurs at n=8 (via sampling, since 2^28 is too large).

kind-pasteur-2026-03-06
"""
import sys, os, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from itertools import combinations

def c3_from_scores(scores):
    n = len(scores)
    return (n*(n-1)*(n-2)//6) - sum(s*(s-1)//2 for s in scores)

def generate_random_tournament(n):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1
    return T

def count_directed_cycles_on_subset(T, verts):
    k = len(verts)
    if k < 3:
        return 0
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for vi in range(k):
            if not (mask & (1 << vi)):
                continue
            c = dp.get((mask, vi), 0)
            if c == 0:
                continue
            for ui in range(k):
                if mask & (1 << ui):
                    continue
                if T[verts[vi]][verts[ui]]:
                    key = (mask | (1 << ui), ui)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << k) - 1
    total = 0
    for vi in range(1, k):
        c = dp.get((full, vi), 0)
        if c > 0 and T[verts[vi]][verts[0]]:
            total += c
    return total

def count_ck(T, n, k):
    total = 0
    for combo in combinations(range(n), k):
        total += count_directed_cycles_on_subset(T, list(combo))
    return total

random.seed(42)
n = 8
SAMPLES = 200000

print(f"n={n}: Sampling {SAMPLES} random tournaments for alpha_1=3...")
found = 0
c3_3_count = 0

for trial in range(SAMPLES):
    if trial % 50000 == 0:
        print(f"  {trial}/{SAMPLES}...", flush=True)
    T = generate_random_tournament(n)
    scores = [sum(T[i]) for i in range(n)]
    c3 = c3_from_scores(scores)

    if c3 > 3:
        continue
    if c3 == 3:
        c3_3_count += 1

    # Need c3 + c5 + c7 = 3
    target = 3 - c3
    c5 = count_ck(T, n, 5)
    if c5 > target:
        continue
    c7 = count_ck(T, n, 7)
    if c5 + c7 > target:
        continue
    # No need to check higher odd cycles since c3+c5+c7 <= 3 already
    alpha1 = c3 + c5 + c7
    if alpha1 == 3:
        found += 1
        if found <= 5:
            print(f"  FOUND: c3={c3}, c5={c5}, c7={c7}")

print(f"\nResults: {found} tournaments with alpha_1=3 out of {SAMPLES}")
print(f"  c3=3 found: {c3_3_count}")
print(f"  (If found > 0, alpha_1=3 exists at n={n})")
print("Done.")
