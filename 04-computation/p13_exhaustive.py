#!/usr/bin/env python3
"""
Exhaustive enumeration of ALL circulant tournaments at p=13.
There are 2^6 = 64 connection sets. Find the global H-maximizer.

Also: landscape analysis at p=13 — how does H vary across all orientations?
This is the PHASE TRANSITION point.

opus-2026-03-12-S62c
"""

import numpy as np
from itertools import combinations

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def count_H_dp(A):
    n = len(A)
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i][i] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for u in range(n):
            if not (mask & (1 << u)):
                continue
            if dp[mask][u] == 0:
                continue
            for v in range(n):
                if mask & (1 << v):
                    continue
                if A[u][v]:
                    dp[mask | (1 << v)][v] += dp[mask][u]
    return int(sum(dp[full]))

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

p = 13
m = (p - 1) // 2  # = 6
QR = get_QR(p)

print(f"p = {p}, m = {m}")
print(f"QR = {QR}")
print(f"Interval = {list(range(1, m+1))}")
print(f"Total circulant tournaments: 2^{m} = {2**m}")
print()

# Each tournament is determined by choosing one from each pair {j, p-j}
pairs = [(j, p-j) for j in range(1, m+1)]
print(f"Pairs: {pairs}")
print()

# Enumerate all 2^m tournaments
results = []
for bits in range(2**m):
    S = []
    for i in range(m):
        if bits & (1 << i):
            S.append(pairs[i][1])
        else:
            S.append(pairs[i][0])
    S = sorted(S)
    A = make_tournament(p, S)
    H = count_H_dp(A)

    # Identify known tournaments
    name = ""
    if S == sorted(list(range(1, m+1))):
        name = " *** INTERVAL ***"
    elif S == sorted(QR):
        name = " *** PALEY ***"

    results.append((H, S, name))

# Sort by H
results.sort(key=lambda x: -x[0])

print("=" * 72)
print("ALL CIRCULANT TOURNAMENTS AT p=13, RANKED BY H")
print("=" * 72)
print()
print(f"{'Rank':>4s}  {'H':>12s}  {'S':>40s}  {'Notes'}")
print("-" * 72)
for rank, (H, S, name) in enumerate(results, 1):
    print(f"{rank:>4d}  {H:>12d}  {str(S):>40s}  {name}")

print()
H_int = [H for H, S, _ in results if 'INTERVAL' in _][0]
paley_results = [H for H, S, _ in results if 'PALEY' in _]

print(f"NOTE: p=13 ≡ 1 mod 4, so NO Paley tournament exists!")
print(f"  QR = {QR} is symmetric (j and p-j both in QR), not a valid tournament set.")
print()
print(f"Interval rank: {[i+1 for i, (H,S,n) in enumerate(results) if 'INTERVAL' in n][0]} / {len(results)}")
print(f"H(Interval) = {H_int}")

# Global max
H_max, S_max, _ = results[0]
print(f"\nGlobal maximum: H = {H_max}")
print(f"  Achieved by S = {S_max}")
print(f"  Is interval: {S_max == sorted(list(range(1, m+1)))}")

# Check affine equivalence of top entries to interval
print()
print("Affine equivalences of top 10:")
for rank, (H, S, _) in enumerate(results[:10], 1):
    equiv_to = None
    for a in range(1, p):
        aS = sorted([(a * s) % p for s in S])
        if aS == sorted(list(range(1, m+1))):
            equiv_to = f"interval * {a}"
            break
        if aS == sorted(QR):
            equiv_to = f"Paley * {a}"
            break
    print(f"  Rank {rank}: S={S}, H={H}, equiv: {equiv_to or 'unique'}")

# Spectral analysis of top 5
print()
print("=" * 72)
print("SPECTRAL ANALYSIS OF TOP 5")
print("=" * 72)
print()

omega = np.exp(2j * np.pi / p)
for rank, (H, S, name) in enumerate(results[:5], 1):
    eigs = [abs(sum(omega**(j*s) for s in S)) for j in range(1, p)]
    top3 = sorted(eigs, reverse=True)[:3]
    concentration = top3[0]**2 / sum(e**2 for e in eigs) * 100
    print(f"  Rank {rank}: S={S}, H={H}")
    print(f"    Top eigenvalue magnitudes: {', '.join(f'{e:.3f}' for e in top3)}")
    print(f"    Concentration (top-1 |μ|²): {concentration:.1f}%")
    print(f"    Stab size: {sum(1 for a in range(1,p) if sorted([(a*s)%p for s in S]) == S)}")
    print()

# Distribution of H values
print("=" * 72)
print("H DISTRIBUTION AT p=13")
print("=" * 72)
print()
H_vals = [H for H, _, _ in results]
print(f"  Min H:    {min(H_vals)}")
print(f"  Max H:    {max(H_vals)}")
print(f"  Mean H:   {np.mean(H_vals):.0f}")
print(f"  Median H: {np.median(H_vals):.0f}")
print(f"  Std H:    {np.std(H_vals):.0f}")
print(f"  CV:       {np.std(H_vals)/np.mean(H_vals)*100:.1f}%")

# How many are above interval? Above Paley?
above_int = sum(1 for H, _, _ in results if H > H_int)
above_pal = 0  # no Paley at p=13
print(f"  Above Interval: {above_int}")
print(f"  Above Paley: {above_pal}")

# Hamming distance to interval for each tournament
print()
print("=" * 72)
print("H vs HAMMING DISTANCE FROM INTERVAL")
print("=" * 72)
print()
S_int_set = set(range(1, m+1))
print(f"  {'Hamming':>7s}  {'Count':>5s}  {'Mean H':>12s}  {'Max H':>12s}  {'Best S'}")
for d in range(m+1):
    subset = [(H, S) for H, S, _ in results if len(set(S) - S_int_set) == d]
    if subset:
        mean_H = np.mean([H for H, _ in subset])
        max_H = max(H for H, _ in subset)
        best_S = max(subset, key=lambda x: x[0])[1]
        print(f"  {d:>7d}  {len(subset):>5d}  {mean_H:>12.0f}  {max_H:>12d}  {best_S}")

print()
print("DONE.")
