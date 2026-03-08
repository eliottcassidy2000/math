#!/usr/bin/env python3
"""
maximizer_betti_n8_highdim.py — Check higher Betti at n=8

At n=8 path complex goes up to dim 7 (Ham paths).
We only checked max_dim=5. Check if beta_6 or beta_7 exist.

Use a known H=661 tournament from the fast search.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time, random
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

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

# Find two H=661 tournaments: one with beta_4=1, one contractible
n = 8
best = []

for trial in range(100000):
    A = random_tournament(n)
    H = H_tournament(A, n)
    if H == 661:
        best.append([row[:] for row in A])
        if len(best) >= 4:
            break

print(f"Found {len(best)} H=661 tournaments")

# Try computing at higher max_dim
for i, A in enumerate(best):
    print(f"\n--- Tournament {i+1} ---")
    for md in [4, 5, 6, 7]:
        t0 = time.time()
        try:
            beta = path_betti_numbers(A, n, max_dim=md)
            beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(md+1)]
            print(f"  max_dim={md}: beta={beta_list} ({time.time()-t0:.1f}s)")
        except Exception as e:
            print(f"  max_dim={md}: FAILED ({time.time()-t0:.1f}s): {e}")
            break

# Also try a random non-max tournament for comparison
print("\n--- Random tournament (for comparison) ---")
random.seed(999)
A_rand = random_tournament(n)
H_rand = H_tournament(A_rand, n)
print(f"H={H_rand}")
for md in [4, 5, 6]:
    t0 = time.time()
    try:
        beta = path_betti_numbers(A_rand, n, max_dim=md)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(md+1)]
        print(f"  max_dim={md}: beta={beta_list} ({time.time()-t0:.1f}s)")
    except Exception as e:
        print(f"  max_dim={md}: FAILED ({time.time()-t0:.1f}s): {e}")
        break
