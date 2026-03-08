#!/usr/bin/env python3
"""
betti_landscape_n9.py — How does β₅ vary near the n=9 maximum?

From beta5_n9_check.py:
- H=3255 (regular): β₃=1 (S-phase!)
- H=3075: β₄=4, β₅=2 (!)
- H=3069: β₄=1
- H=3067: contractible

Questions:
1. Does the H=3357 maximizer ALWAYS have β₅=10? (check multiple)
2. What Betti vectors appear at H=3333 (second-highest regular)?
3. Is the β₁/β₃/β₅ odd-Betti-only pattern preserved?

Strategy: quick screening via just the high-H regular tournaments.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, random, os
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def random_regular_tournament(n):
    assert n % 2 == 1
    target = (n - 1) // 2
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    for _ in range(1000):
        scores = [sum(A[i]) for i in range(n)]
        if max(abs(s - target) for s in scores) == 0:
            break
        high_v = max(range(n), key=lambda v: scores[v])
        low_v = min(range(n), key=lambda v: scores[v])
        if high_v == low_v:
            break
        if A[high_v][low_v]:
            A[high_v][low_v] = 0
            A[low_v][high_v] = 1
        else:
            candidates = [w for w in range(n) if w != high_v and w != low_v
                         and A[high_v][w] and A[w][low_v]]
            if candidates:
                w = random.choice(candidates)
                A[high_v][w] = 0
                A[w][high_v] = 1
                A[w][low_v] = 0
                A[low_v][w] = 1
    return A

# ===== Check 3 distinct H=3357 maximizers =====
print("=" * 70)
print("BETTI CONSISTENCY CHECK: Multiple H=3357 maximizers")
print("=" * 70)

n = 9
t0 = time.time()

checked = 0
for trial in range(3):
    while True:
        A = random_regular_tournament(n)
        scores = [sum(A[i]) for i in range(n)]
        if max(scores) - min(scores) > 0:
            continue
        H = H_tournament(A, n)
        if H == 3357:
            break

    beta = path_betti_numbers(A, n, max_dim=5)
    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    print(f"  Maximizer #{trial+1}: β = {beta_list} ({time.time()-t0:.1f}s)")
    checked += 1

# ===== Check H=3333 (second-best regular H observed) =====
print("\n" + "=" * 70)
print("SECOND-TIER REGULAR TOURNAMENTS")
print("=" * 70)

h_values = {}
print("Searching for non-max regular tournaments...")
for trial in range(30000):
    A = random_regular_tournament(n)
    scores = [sum(A[i]) for i in range(n)]
    if max(scores) - min(scores) > 0:
        continue
    H = H_tournament(A, n)
    if H not in h_values and H != 3357 and H >= 3200:
        h_values[H] = [row[:] for row in A]
        print(f"  Found H={H} at trial {trial}")

print(f"\nDistinct regular H values >= 3200 (excluding 3357): {sorted(h_values.keys(), reverse=True)}")

# Compute Betti for top 3 non-max regular H values
for H in sorted(h_values.keys(), reverse=True)[:3]:
    A = h_values[H]
    t1 = time.time()
    beta = path_betti_numbers(A, n, max_dim=5)
    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    print(f"  H={H}: β = {beta_list} ({time.time()-t1:.1f}s)")

# ===== Check: does β₃ ever coexist with β₅ at n=9? =====
print("\n" + "=" * 70)
print("β₃ AND β₅ COEXISTENCE CHECK")
print("=" * 70)
print("From prior data: H=3255 has β₃=1, H=3075 has β₅=2. Are these β₃ AND β₅ at once?")
print("β₃ would need S-phase (3-cycle structure), β₅ is higher-dimensional.")
print("From beta5_n9_check: H=3075 has β=[1,0,0,0,4,2] — β₃=0, β₅=2. So NO coexistence.")
print("H=3255 has β=[1,0,0,1,0,0] — β₃=1, β₅=0. Also NO coexistence.")
print("This preserves mutual exclusivity of odd Betti numbers at n=9.")

print(f"\nTotal time: {time.time()-t0:.1f}s")
