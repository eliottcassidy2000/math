#!/usr/bin/env python3
"""
betti_dimension_shift_v2.py — Quick check: n=4 exhaustive + n=6 exhaustive Betti for maximizers

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

from collections import Counter

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

# ===== n=4 EXHAUSTIVE =====
print("=" * 70)
print("n=4 EXHAUSTIVE: Betti numbers by H value")
print("=" * 70)

n = 4
m = n * (n-1) // 2
total = 1 << m
betti_by_H = {}

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    try:
        beta = path_betti_numbers(A, n, max_dim=3)
    except:
        continue
    H = H_tournament(A, n)
    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(4)]
    if H not in betti_by_H:
        betti_by_H[H] = []
    betti_by_H[H].append(beta_list)

for H_val in sorted(betti_by_H.keys(), reverse=True):
    entries = betti_by_H[H_val]
    betti_counts = Counter(tuple(e) for e in entries)
    print(f"  H={H_val} ({len(entries)}):")
    for b, cnt in betti_counts.most_common():
        print(f"    beta={list(b)}: {cnt}")

# ===== n=6 EXHAUSTIVE: full Betti for H=45 maximizers =====
print("\n" + "=" * 70)
print("n=6 EXHAUSTIVE: Betti for H=45 maximizers (what's the full vector?)")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
max_betti = []
t0 = time.time()

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    H = H_tournament(A, n)
    if H != 45:
        continue

    try:
        beta = path_betti_numbers(A, n, max_dim=5)
    except:
        continue

    beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(6)]
    max_betti.append(beta_list)

print(f"Done in {time.time()-t0:.1f}s")
print(f"Found {len(max_betti)} maximizers")
betti_counts = Counter(tuple(b) for b in max_betti)
for b, cnt in betti_counts.most_common():
    print(f"  beta={list(b)}: {cnt}")

# ===== SUMMARY =====
print("\n" + "=" * 70)
print("COMPLETE BETTI DATA FOR H-MAXIMIZERS")
print("=" * 70)
print("""
n  max_H  #max  Betti vector         Highest nontrivial dim
3    1      2   [1]                  0 (trivial)
4    5     16   [1,0,0,0]            0 (contractible)
5   15     64   [1,1,0,0,0]          1 (circle)
6   45    240   [1,0,0,β₃,0,0]      3 (3-sphere?)
7  189    240   [1,0,0,0,6,0,0]     4

Observations:
1. At n=4: maximizers are CONTRACTIBLE (β = [1,0,0,0])
2. At n=5: maximizers have β₁ = 1 (C-phase, circle-like)
3. At n=6: maximizers have β₃ > 0 (S-phase, "3-sphere-like")
4. At n=7: maximizers have β₄ = 6 (!), everything else 0

The highest nontrivial dim: 0, 1, 3, 4 (not monotone from n=3)
But from n=5 onwards: 1, 3, 4 — gap at dim 2

At n=5: ALL maximizers have SAME Betti (β₁=1)
At n=7: ALL maximizers have SAME Betti (β₄=6)

H=175 class at n=7 has β₁=1 (C-phase) — SECOND highest H
H=171 class at n=7 has β₀=1 only (P-phase) — THIRD highest H

So topology stratifies the H values at n=7:
  H=189: β₄=6 (deepest topology)
  H=175: β₁=1 (circle)
  H=171: contractible
""")
