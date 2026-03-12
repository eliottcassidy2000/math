#!/usr/bin/env python3
"""
Landscape analysis at the TRUE crossover primes p = 19, 23.

For p ≡ 3 mod 4: p=7, 11 (Paley wins), p=19, 23 (Interval wins).
The crossover is between p=11 and p=19.

We can't do exhaustive enumeration at p=19 (2^9 = 512 tournaments, each
requiring DP on n=19 which is borderline). But we CAN:

1. Check ALL single swaps from both Interval and Paley at p=19
2. Verify the landscape has a SINGLE basin of attraction at large p
3. Quantify the "distance" between Paley and Interval in orientation space
4. Check: is the crossover at p=19, or is there p=17 (≡ 1 mod 4, skip)?

The sequence of p ≡ 3 mod 4 primes: 3, 7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83

opus-2026-03-12-S62c
"""

import numpy as np
from itertools import combinations
import time

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=int)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def count_H_dp(A):
    """Count Hamiltonian paths via DP bitmask. Works up to n~20."""
    n = len(A)
    if n > 20:
        return -1  # too slow
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

def eigenvalues(p, S):
    omega = np.exp(2j * np.pi / p)
    return [sum(omega**(j*s) for s in S) for j in range(1, p)]

print("=" * 72)
print("CROSSOVER LANDSCAPE: p ≡ 3 mod 4 primes")
print("=" * 72)
print()

# First: verify which primes have Paley tournament
for p in [3, 7, 11, 19, 23, 31]:
    if p % 4 != 3:
        print(f"  p={p}: ≡ {p%4} mod 4, SKIP (no Paley)")
        continue
    m = (p - 1) // 2
    QR = get_QR(p)
    # Verify QR is a valid tournament connection set
    # For each pair {j, p-j}, exactly one should be in QR
    valid = all((j in QR) != ((p-j) in QR) for j in range(1, m+1))
    print(f"  p={p}: ≡ 3 mod 4, m={m}, QR valid tournament: {valid}")

print()

# p=19 analysis
print("=" * 72)
print("p = 19 DETAILED ANALYSIS")
print("=" * 72)
print()

p = 19
m = (p - 1) // 2  # = 9
QR = get_QR(p)
S_int = list(range(1, m + 1))
pairs = [(j, p-j) for j in range(1, m+1)]

print(f"p = {p}, m = {m}")
print(f"QR = {QR}")
print(f"Interval = {S_int}")
print(f"Pairs: {pairs}")
print()

# Hamming distance between Paley and Interval orientations
paley_bits = tuple(0 if j in QR else 1 for j, _ in pairs)
int_bits = tuple(0 for _ in pairs)
hamming = sum(a != b for a, b in zip(paley_bits, int_bits))
print(f"Paley orientation bits: {paley_bits}")
print(f"Interval orientation bits: {int_bits}")
print(f"Hamming distance: {hamming} out of {m}")
print()

# Compute H for Paley and Interval
t0 = time.time()
A_pal = make_tournament(p, QR)
H_pal = count_H_dp(A_pal)
t1 = time.time()
print(f"H(Paley) = {H_pal} (computed in {t1-t0:.1f}s)")

t0 = time.time()
A_int = make_tournament(p, S_int)
H_int = count_H_dp(A_int)
t1 = time.time()
print(f"H(Interval) = {H_int} (computed in {t1-t0:.1f}s)")
print(f"Gap: H(Int) - H(Pal) = {H_int - H_pal:+d} ({(H_int-H_pal)/H_pal*100:+.3f}%)")
print()

# Single-swap analysis for Interval
print("Single-swap analysis from Interval:")
int_swaps = []
for a in S_int:
    for b in range(1, p):
        if b not in S_int:
            S_new = sorted([x for x in S_int if x != a] + [b])
            if len(S_new) == m and len(set(S_new)) == m:
                # Check it's a valid tournament set
                valid = all((j in S_new) != ((p-j) in S_new) for j in range(1, m+1))
                if valid:
                    A_new = make_tournament(p, S_new)
                    H_new = count_H_dp(A_new)
                    delta = H_new - H_int
                    int_swaps.append((a, b, H_new, delta, delta/H_int*100))

int_swaps.sort(key=lambda x: -x[3])
print(f"  Total valid single swaps: {len(int_swaps)}")
print(f"  {'a→b':>8s}  {'H_new':>15s}  {'ΔH':>15s}  {'ΔH/H':>8s}")
for a, b, h, d, pct in int_swaps[:10]:
    print(f"  {a}→{b:>2d}  {h:>15d}  {d:>+15d}  {pct:>+7.4f}%")
all_neg = all(s[3] <= 0 for s in int_swaps)
print(f"  All swaps negative: {all_neg}")
print(f"  → Interval is {'LOCAL MAX' if all_neg else 'NOT local max'} at p=19")
print()

# Single-swap analysis for Paley
print("Single-swap analysis from Paley:")
pal_swaps = []
for a in QR:
    for b in range(1, p):
        if b not in QR:
            S_new = sorted([x for x in QR if x != a] + [b])
            if len(S_new) == m and len(set(S_new)) == m:
                valid = all((j in S_new) != ((p-j) in S_new) for j in range(1, m+1))
                if valid:
                    A_new = make_tournament(p, S_new)
                    H_new = count_H_dp(A_new)
                    delta = H_new - H_pal
                    pal_swaps.append((a, b, H_new, delta, delta/H_pal*100))

pal_swaps.sort(key=lambda x: -x[3])
print(f"  Total valid single swaps: {len(pal_swaps)}")
print(f"  {'a→b':>8s}  {'H_new':>15s}  {'ΔH':>15s}  {'ΔH/H':>8s}")
for a, b, h, d, pct in pal_swaps[:10]:
    print(f"  {a}→{b:>2d}  {h:>15d}  {d:>+15d}  {pct:>+7.4f}%")
pal_local_max = all(s[3] <= 0 for s in pal_swaps)
print(f"  All swaps negative: {pal_local_max}")
print(f"  → Paley is {'LOCAL MAX' if pal_local_max else 'NOT local max'} at p=19")
print()

# Is Paley STILL a local max at p=19?
if not pal_local_max:
    best_swap = pal_swaps[0]
    print(f"  Best swap from Paley: {best_swap[0]}→{best_swap[1]}")
    print(f"    Brings Paley TOWARD interval? ", end="")
    new_S = sorted([x for x in QR if x != best_swap[0]] + [best_swap[1]])
    ham_new = sum(1 for j in range(1, m+1) if (j in new_S) != (j in S_int))
    print(f"Hamming to Int: {ham_new} (was {hamming})")
    if ham_new < hamming:
        print(f"    YES — the gradient from Paley points TOWARD Interval!")
    else:
        print(f"    NO — gradient points away from Interval")

# Spectral comparison
print()
print("=" * 72)
print("SPECTRAL LANDSCAPE")
print("=" * 72)
print()

eigs_int = eigenvalues(p, S_int)
eigs_pal = eigenvalues(p, QR)

int_mags = sorted([abs(e) for e in eigs_int], reverse=True)
pal_mags = sorted([abs(e) for e in eigs_pal], reverse=True)

print(f"  Interval |μ| (top 5): {', '.join(f'{v:.4f}' for v in int_mags[:5])}")
print(f"  Paley |μ| (top 5):    {', '.join(f'{v:.4f}' for v in pal_mags[:5])}")
print(f"  |μ_1| ratio: {int_mags[0]/pal_mags[0]:.4f}")
print(f"  Spectral entropy: Int={-sum(v**2/sum(u**2 for u in int_mags)*np.log(v**2/sum(u**2 for u in int_mags)) for v in int_mags if v > 1e-10):.4f}, Pal={-sum(v**2/sum(u**2 for u in pal_mags)*np.log(v**2/sum(u**2 for u in pal_mags)) for v in pal_mags if v > 1e-10):.4f}")

# The gradient flow: starting from random tournaments, does H increase toward interval?
print()
print("=" * 72)
print("GRADIENT ASCENT FROM RANDOM STARTING POINTS")
print("=" * 72)
print()

import random
random.seed(42)

n_starts = 5
for trial in range(n_starts):
    # Random connection set
    S = []
    for j, pj in pairs:
        if random.random() < 0.5:
            S.append(j)
        else:
            S.append(pj)
    S = sorted(S)

    # Greedy ascent
    path = [S[:]]
    H_curr = count_H_dp(make_tournament(p, S))
    steps = 0
    while True:
        best_next = None
        best_H = H_curr
        for a in S:
            for b in range(1, p):
                if b not in S:
                    S_new = sorted([x for x in S if x != a] + [b])
                    if len(set(S_new)) == m:
                        valid = all((j in S_new) != ((p-j) in S_new) for j in range(1, m+1))
                        if valid:
                            H_new = count_H_dp(make_tournament(p, S_new))
                            if H_new > best_H:
                                best_H = H_new
                                best_next = S_new
        if best_next is None:
            break
        S = best_next
        H_curr = best_H
        path.append(S[:])
        steps += 1
        if steps > 20:
            break

    ham_to_int = sum(1 for j in range(1, m+1) if (j in S) != (j in S_int))
    print(f"  Trial {trial+1}: start S={path[0]}")
    print(f"    Steps: {steps}, Final S={S}, H={H_curr}")
    print(f"    Hamming to Interval: {ham_to_int}")
    is_int_equiv = any(sorted([(a*s)%p for s in S]) == S_int for a in range(1, p))
    print(f"    Equivalent to Interval: {is_int_equiv}")
    print()

print("=" * 72)
print("CRITICAL CONCLUSION")
print("=" * 72)
print()
print(f"""
  AT p=19:
    H(Interval) = {H_int}
    H(Paley) = {H_pal}
    Gap: {(H_int-H_pal)/H_pal*100:+.3f}%

  The orientation landscape at p=19 shows:
  - Interval is local max: {all_neg}
  - Paley is local max: {pal_local_max}

  PHASE TRANSITION SUMMARY (p ≡ 3 mod 4):
  p=7:  Paley wins, Paley is global max, Interval is NOT local max
  p=11: Paley wins, Paley is global max, Interval is NOT local max
  p=19: Interval wins, both may be local maxima
  p=23: Interval wins (kind-pasteur: 0/36 beat Interval)

  The crossover from Paley→Interval happens between p=11 and p=19.
  (p=13 ≡ 1 mod 4 is outside the Paley domain.)
""")

print("DONE.")
