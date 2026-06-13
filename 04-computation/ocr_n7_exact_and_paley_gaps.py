#!/usr/bin/env python3
"""
ocr_n7_exact_and_paley_gaps.py — kind-pasteur-2026-03-21-S16

THREE INVESTIGATIONS IN ONE:

1. OCR(7) exact: exhaustive computation over all 2^21 = 2,097,152 tournaments
   to find the exact rational R^2(S2, H).

2. PALEY MAXIMIZER: Analyze why Paley tournaments maximize H.
   Key idea: Paley T_p has the flattest eigenvalue distribution AND
   the most symmetric cycle structure. Both contribute to max H.
   NEW ANGLE: Use OCF. H = 1 + 2*alpha_1 + 4*alpha_2 + ...
   Paley maximizes the weighted sum. Is alpha_1 or alpha_2 dominant?

3. PERMANENT GAPS: H=7 and H=21.
   From OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
   H=7 requires 1+2a1+4a2+...=7, so 2a1+4a2+...=6, so a1+2a2+...=3.
   H=21 requires a1+2a2+...=10.
   WHY are these impossible?
   NEW ANGLE: The OCF decomposition constrains which (a1,a2,...) are achievable.
   Is there a PACKING constraint?

Author: kind-pasteur-2026-03-21-S16
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
from math import comb, gcd, factorial
from fractions import Fraction
import time

def build_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def ham_paths_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    key = (mask | (1 << w), w)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# ============================================================
# INVESTIGATION 1: OCR(7) EXACT
# ============================================================

print("=" * 72)
print("  INVESTIGATION 1: OCR(7) EXACT")
print("=" * 72)

n = 7
m = n*(n-1)//2  # 21
total = 1 << m   # 2,097,152

print(f"  Enumerating all {total:,} tournaments on n={n} vertices...")
t0 = time.time()

# We need: sum(H), sum(S2), sum(H^2), sum(S2^2), sum(H*S2)
# Use Python ints to avoid overflow
sum_H = 0
sum_S2 = 0
sum_H2 = 0
sum_S22 = 0
sum_HS2 = 0
count = 0

# Also track H distribution for gap analysis
H_counter = Counter()

# Track Paley tournament specifically
# Paley T_7: QR mod 7 = {1, 2, 4}, connection set S = {1, 2, 4}
# T_7[i][j] = 1 iff (j-i) mod 7 in {1, 2, 4}
paley_bits = 0
paley_A = [[0]*7 for _ in range(7)]
QR7 = {1, 2, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j - i) % 7 in QR7:
            paley_A[i][j] = 1

for bits in range(total):
    A = build_tournament(n, bits)
    scores = [sum(A[i]) for i in range(n)]
    S2 = sum((s - 3)**2 for s in scores)  # (n-1)/2 = 3
    H = ham_paths_dp(A, n)

    sum_H += H
    sum_S2 += S2
    sum_H2 += H * H
    sum_S22 += S2 * S2
    sum_HS2 += H * S2
    count += 1

    H_counter[H] += 1

    if count % 500000 == 0:
        elapsed = time.time() - t0
        rate = count / elapsed
        remaining = (total - count) / rate
        print(f"    {count:,}/{total:,} ({100*count/total:.1f}%) "
              f"elapsed {elapsed:.0f}s, remaining ~{remaining:.0f}s")

elapsed = time.time() - t0
print(f"  Done in {elapsed:.1f}s ({count/elapsed:.0f} tournaments/sec)")
print()

# Exact R^2
N = total
numer_cov = N * sum_HS2 - sum_H * sum_S2
var_S2 = N * sum_S22 - sum_S2**2
var_H = N * sum_H2 - sum_H**2

r2_numer = numer_cov**2
r2_denom = var_S2 * var_H

g = gcd(r2_numer, r2_denom)
p, q = r2_numer // g, r2_denom // g

print(f"  EXACT R^2(S2, H) = {p} / {q}")
print(f"  1 - R^2 = {q - p} / {q}")
print(f"  OCR(7) = {p/q:.12f}")
print()

# Mean H
mean_H = Fraction(sum_H, N)
exp_random = factorial(n) // (1 << (n-1))
print(f"  Mean H = {mean_H} = {float(mean_H):.6f}")
print(f"  Expected random = {exp_random}")
print()

# ============================================================
# THE OCR SEQUENCE
# ============================================================

print("  THE EXACT OCR SEQUENCE:")
print(f"  n=3: OCR = 1/1")
print(f"  n=4: OCR = 1/1")
print(f"  n=5: OCR = 18/19")
print(f"  n=6: OCR = 12/13")
print(f"  n=7: OCR = {p}/{q}")
print()

# Analyze the denominators
denoms = [1, 1, 19, 13, q]
print(f"  Denominators: {denoms}")
print(f"  1-OCR denominators: {denoms}")

# Check various patterns
print(f"\n  Pattern checks for denominator sequence:")
if q - p == 1:
    print(f"  1-OCR(7) has numerator 1 => 1/{q}")
else:
    print(f"  1-OCR(7) = {q-p}/{q}")

# ============================================================
# INVESTIGATION 2: H=7 and H=21 GAP ANALYSIS VIA OCF
# ============================================================

print()
print("=" * 72)
print("  INVESTIGATION 2: PERMANENT GAPS VIA OCF DECOMPOSITION")
print("=" * 72)

# From OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
# H=7: 2a1+4a2+8a3+...=6 => a1+2a2+4a3+...=3
# Possibilities: (a1,a2,a3,...) = (3,0,0,...), (1,1,0,...), (0,0,0,...,[higher])
# But alpha_1 >= alpha_2 always (can't have disjoint pairs without individual cycles)
# Actually NO: alpha_2 counts pairs, alpha_1 counts total.
# Need alpha_1 >= 2*alpha_2 (each pair uses 2 cycles)

# H=21: 2a1+4a2+8a3+...=20 => a1+2a2+4a3+...=10
# Many possible decompositions.

print("\n  H=7 requires a1+2a2+4a3+...=3")
print("  Possible (a1,a2): (3,0), (1,1)")
print("  H=21 requires a1+2a2+4a3+...=10")
print("  Possible: (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)")
print()

# Let's compute (alpha_1, alpha_2) for ALL n=7 tournaments and see
# which (a1,a2) pairs actually occur.

print("  Computing (alpha_1, alpha_2) distribution at n=7...")
print("  (Using score-class sampling for speed)")

# For each score class at n=7, sample a few tournaments and compute alpha_1
# alpha_1 at n=7 = c3 + c5_directed + c7_directed
# Actually, let's compute c3 (= from scores) and track H values

# From the exhaustive H count, let's analyze what (c3, H) pairs exist
# c3 is determined by scores, so we group by c3

score_H = defaultdict(lambda: defaultdict(int))
c3_H = defaultdict(lambda: defaultdict(int))

print("  Re-scanning for c3-H distribution...")
for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    c3 = comb(n, 3) - sum(s*(s-1)//2 for s in scores)
    H = ham_paths_dp(A, n)
    c3_H[c3][H] += 1

print()
print("  c3 -> H distribution (showing H values near 7 and 21):")
for c3_val in sorted(c3_H.keys()):
    H_vals = c3_H[c3_val]
    H_min = min(H_vals.keys())
    H_max = max(H_vals.keys())
    n_distinct = len(H_vals)
    # Check if H=7 or H=21 appears
    has_7 = 7 in H_vals
    has_21 = 21 in H_vals
    marker = ""
    if has_7: marker += " H=7!"
    if has_21: marker += " H=21!"
    print(f"    c3={c3_val:3d}: H in [{H_min}, {H_max}], {n_distinct} distinct values{marker}")

print()
print("  Full H spectrum at n=7 (confirming gaps):")
all_H_values = sorted(H_counter.keys())
print(f"  {len(all_H_values)} distinct H values")
print(f"  Gaps in [1..50]:", end=" ")
for h in range(1, 51, 2):
    if h not in H_counter:
        print(h, end=" ")
print()

# ============================================================
# INVESTIGATION 3: PALEY MAXIMIZER
# ============================================================

print()
print("=" * 72)
print("  INVESTIGATION 3: PALEY MAXIMIZER ANALYSIS")
print("=" * 72)

# Find all tournaments achieving max H
max_H = max(H_counter.keys())
print(f"\n  Max H at n=7: {max_H}")
print(f"  Number of max-H tournaments: {H_counter[max_H]}")

# Compute H of Paley T_7
H_paley = ham_paths_dp(paley_A, 7)
print(f"  H(Paley T_7) = {H_paley}")
print(f"  Paley achieves max: {H_paley == max_H}")

# What score does Paley have?
paley_scores = tuple(sorted([sum(paley_A[i]) for i in range(7)]))
paley_c3 = comb(7, 3) - sum(s*(s-1)//2 for s in paley_scores)
print(f"  Paley scores: {paley_scores}")
print(f"  Paley c3: {paley_c3}")

# How many regular tournaments at n=7?
reg_count = 0
reg_max_H = 0
reg_H_vals = Counter()
for bits in range(total):
    A = build_tournament(n, bits)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))
    if scores == (3,3,3,3,3,3,3):
        H = ham_paths_dp(A, n)
        reg_count += 1
        reg_H_vals[H] += 1
        if H > reg_max_H:
            reg_max_H = H

print(f"\n  Regular tournaments (all scores 3): {reg_count}")
print(f"  Regular H distribution: {dict(reg_H_vals)}")
print(f"  Paley ratio: {H_counter[max_H]} / {reg_count} = "
      f"{H_counter[max_H]/reg_count:.4f} of regular tournaments achieve max H")

# ============================================================
# THE DEEPEST QUESTION: Why does the OCR formula give these rationals?
# ============================================================

print()
print("=" * 72)
print("  SYNTHESIS: THE OCR FORMULA")
print("=" * 72)
print()

# Key facts:
# E[H] = n!/2^{n-1} (by linearity of expectation over permutations)
# E[S2] = (n-1)(n-2)/4 (by score variance for random tournaments)
# Cov(H, S2) and Var(H) determine R^2.

# Can we compute E[H^2] theoretically?
# H = sum_{permutations sigma} indicator(sigma is a HP of T)
# H^2 = sum_{sigma, tau} indicator(both are HPs)
# E[H^2] = sum_{sigma,tau} P(both sigma and tau are HPs of random T)

# For random tournament: each arc independent, P(arc in T) = 1/2.
# P(sigma is HP) = (1/2)^{n-1} (need n-1 specific arcs)
# P(sigma AND tau both HPs) = (1/2)^{|arcs(sigma) union arcs(tau)|}

# This is a COMBINATORIAL quantity! Let me compute it.

print("  THEORETICAL E[H^2] COMPUTATION:")
print("  E[H^2] = sum_{sigma,tau in S_n} 2^{-|E(sigma) union E(tau)|}")
print("  where E(sigma) = arcs used by permutation sigma as a HP")
print()

for n_th in [3, 4, 5]:
    total_th = 0
    # Each permutation sigma uses arcs (sigma(0)->sigma(1), ..., sigma(n-2)->sigma(n-1))
    perms = list(permutations(range(n_th)))
    for sigma in perms:
        arcs_sigma = set()
        for k in range(n_th - 1):
            arcs_sigma.add((sigma[k], sigma[k+1]))

        for tau in perms:
            arcs_tau = set()
            for k in range(n_th - 1):
                arcs_tau.add((tau[k], tau[k+1]))

            union_size = len(arcs_sigma | arcs_tau)
            total_th += 2**(-union_size)

    # E[H^2] = total_th / (n!)^2 * (n!)^2 = total_th
    # Wait: E[H^2] = sum_{sigma,tau} P(both HPs) = sum_{sigma,tau} (1/2)^{|union|}
    E_H2 = total_th

    # E[H] = n! * (1/2)^{n-1} = n!/2^{n-1}
    E_H = factorial(n_th) / 2**(n_th - 1)

    # Var(H) = E[H^2] - E[H]^2
    Var_H_th = E_H2 - E_H**2

    print(f"  n={n_th}: E[H^2] = {E_H2:.6f}, E[H] = {E_H:.6f}, Var(H) = {Var_H_th:.6f}")

    # Check against exhaustive
    # 2^m * Var(H) = sum(H^2) - sum(H)^2 / 2^m
    # Total Var = (sum_H2/N - (sum_H/N)^2) where N=2^m
    # In our exact computations earlier:
    # sum_H / N = E[H], sum_H2 / N = E[H^2]
    # So this should match

print()
print("  The key: E[H^2] depends on the OVERLAP structure of permutation pairs.")
print("  This is a purely combinatorial quantity independent of tournament structure.")
print("  The OCR denominator (Var(H)) is thus a fixed combinatorial number.")
print("  The OCR numerator (Cov^2) involves how S2 correlates with HP overlap.")
