#!/usr/bin/env python3
"""
alpha_crossover_smalln.py -- Compute exact alpha decompositions for small n
to map the crossover landscape between α₁, α₂, α₃ dominance.

Computes both interval (cyclic) and Paley tournaments for all odd n ≤ 15.
Uses direct enumeration for n ≤ 13, DP for larger.

Author: opus-2026-04-16-S1
"""

import numpy as np
from itertools import combinations
import time

# ─── Tournament construction ──────────────────────────────────────────────────

def make_interval(n):
    """Interval/cyclic tournament on n vertices (S={1,...,(n-1)/2})."""
    m = (n - 1) // 2
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for k in range(1, m + 1):
            j = (i + k) % n
            A[i][j] = 1
    return A

def make_paley(p):
    """Paley tournament on p vertices (p ≡ 3 mod 4 prime)."""
    assert p % 4 == 3, f"{p} must be ≡ 3 mod 4"
    QR = set(pow(a, 2, p) for a in range(1, p))
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in QR:
                A[i][j] = 1
    return A

# ─── Odd cycle enumeration ────────────────────────────────────────────────────

def find_all_odd_directed_cycles(A):
    """Find all directed odd cycles in tournament A. Returns list of frozensets."""
    n = len(A)
    cycles = set()
    # DFS-based cycle finding
    for start in range(n):
        stack = [(start, [start], 1 << start)]
        while stack:
            v, path, visited = stack.pop()
            for u in range(n):
                if A[v][u]:
                    if u == start and len(path) >= 3 and len(path) % 2 == 1:
                        cycles.add(frozenset(path))
                    elif not (visited & (1 << u)):
                        stack.append((u, path + [u], visited | (1 << u)))
    return list(cycles)

def alpha_decomposition_small(A):
    """
    Compute alpha_k = number of independent sets of size k in Omega(T).
    Omega(T) = conflict graph of odd directed cycles (adjacent iff share vertex).

    Returns dict {k: alpha_k} and H.
    """
    n = len(A)
    cycles = find_all_odd_directed_cycles(A)

    # Build vertex bitmasks for each cycle
    cycle_masks = [sum(1 << v for v in cyc) for cyc in cycles]
    m = len(cycle_masks)

    # Build adjacency: two cycles conflict iff they share a vertex
    # (i.e., their masks have a nonzero AND)

    # Count independent sets by DP
    # dp[S] = 1 if S is an independent set in Omega
    # We want count by size

    alphas = {}
    kmax = n // 3

    # For small m, enumerate directly
    for k in range(1, kmax + 1):
        count = 0
        for combo in combinations(range(m), k):
            # Check all pairs are disjoint (no vertex sharing)
            union = 0
            ok = True
            for idx in combo:
                if union & cycle_masks[idx]:
                    ok = False
                    break
                union |= cycle_masks[idx]
            if ok:
                count += 1
        if count > 0:
            alphas[k] = count
        else:
            break

    H = 1 + sum(2**k * alphas.get(k, 0) for k in range(1, kmax + 1))
    return alphas, H

# ─── Main analysis ────────────────────────────────────────────────────────────

print("=" * 76)
print("ALPHA CROSSOVER LANDSCAPE — Small n exact decompositions")
print("=" * 76)
print()

# Gather data
all_data = {}

for n in [3, 5, 7, 9, 11, 13, 15]:
    print(f"n={n}:")

    # Interval tournament
    A_int = make_interval(n)
    t0 = time.time()
    alphas_int, H_int = alpha_decomposition_small(A_int)
    t_int = time.time() - t0
    all_data[('int', n)] = (alphas_int, H_int)

    a1 = alphas_int.get(1, 0)
    a2 = alphas_int.get(2, 0)
    a3 = alphas_int.get(3, 0)
    a4 = alphas_int.get(4, 0)
    a5 = alphas_int.get(5, 0)

    terms = [2**k * alphas_int.get(k, 0) for k in range(1, 10) if k in alphas_int]
    kmax_nonzero = max(alphas_int.keys()) if alphas_int else 0

    print(f"  Interval: H={H_int:,} ({t_int:.1f}s)")
    print(f"    α₁={a1}, α₂={a2}, α₃={a3}, α₄={a4}, α₅={a5}")
    for k, v in sorted(alphas_int.items()):
        pct = 100 * (2**k * v) / H_int if H_int > 1 else 0
        print(f"    2^{k}·α_{k} = {2**k * v:>15,}  ({pct:.1f}%)")

    if a2 > 0:
        print(f"    α₁/(2α₂)  = {a1/(2*a2):.4f}")
    if a2 > 0 and a3 > 0:
        print(f"    α₃/α₂     = {a3/a2:.4f}")

    # Paley (if n ≡ 3 mod 4 and prime)
    from sympy import isprime
    if n % 4 == 3 and isprime(n):
        A_pal = make_paley(n)
        t0 = time.time()
        alphas_pal, H_pal = alpha_decomposition_small(A_pal)
        t_pal = time.time() - t0
        all_data[('paley', n)] = (alphas_pal, H_pal)

        b1 = alphas_pal.get(1, 0)
        b2 = alphas_pal.get(2, 0)
        b3 = alphas_pal.get(3, 0)

        print(f"  Paley:    H={H_pal:,} ({t_pal:.1f}s)")
        print(f"    α₁={b1}, α₂={b2}, α₃={b3}")
        for k, v in sorted(alphas_pal.items()):
            pct = 100 * (2**k * v) / H_pal if H_pal > 1 else 0
            print(f"    2^{k}·α_{k} = {2**k * v:>15,}  ({pct:.1f}%)")

        winner = "Paley" if H_pal > H_int else ("Interval" if H_int > H_pal else "TIE")
        print(f"  → Winner: {winner} (H_pal={H_pal:,}, H_int={H_int:,})")

        # Contribution differences
        delta_H = H_pal - H_int
        delta_2a1 = 2*(b1 - a1)
        delta_4a2 = 4*(b2 - a2)
        delta_8a3 = 8*(b3 - a3)
        print(f"  → ΔH = {delta_H:+,}")
        print(f"     2·Δα₁ = {delta_2a1:+,}")
        print(f"     4·Δα₂ = {delta_4a2:+,}")
        print(f"     8·Δα₃ = {delta_8a3:+,}")

    print()

# ─── Crossover summary table ──────────────────────────────────────────────────
print("=" * 76)
print("CROSSOVER SUMMARY TABLE (interval tournament)")
print("=" * 76)
print(f"{'n':>4} {'H':>20} {'α₁':>15} {'α₂':>15} {'α₃':>12} {'α₁/2α₂':>9} {'α₃/α₂':>8} {'dom_k':>6}")
print("-" * 100)

# Known large-n data
large_n_data = {
    17: {'alpha1': 1651334601, 'alpha2': 1482234998, 'alpha3': 458011858,
         'alpha4': 45997104, 'alpha5': 1800368, 'H': 13689269499},
    19: {'alpha1': 126443605257, 'alpha2': 122111579294, 'alpha3': 42960731622,
         'alpha4': 5521030944, 'alpha5': 331078344, 'alpha6': 4100656, 'H': 1184212824763},
    21: {'alpha1': 12030499746751, 'alpha2': 12330182836208, 'alpha3': 4796354751404,
         'alpha4': 738531326288, 'alpha5': 58868297768, 'alpha6': 1454221328,
         'alpha7': 12571712, 'H': 125547534942879},
    23: {'alpha1': 1391602826199187, 'alpha2': 1499656616321278, 'alpha3': 632921002322216,
         'alpha4': 111796734828336, 'alpha5': 10945293151712, 'alpha6': 412282843184,
         'alpha7': 7454017376, 'H': 16011537490557279},
}

for n in [3, 5, 7, 9, 11, 13, 15]:
    if ('int', n) in all_data:
        alphas, H = all_data[('int', n)]
        a1 = alphas.get(1, 0)
        a2 = alphas.get(2, 0)
        a3 = alphas.get(3, 0)
        ratio_12 = a1/(2*a2) if a2 > 0 else float('inf')
        ratio_32 = a3/a2 if a2 > 0 else 0
        # Dominant k: find k that maximizes 2^k * alpha_k
        dom_k = max(alphas, key=lambda k: 2**k * alphas[k]) if alphas else 0
        print(f"{n:>4} {H:>20,} {a1:>15,} {a2:>15,} {a3:>12,} {ratio_12:>9.4f} {ratio_32:>8.4f} {dom_k:>6}")

for n in sorted(large_n_data.keys()):
    d = large_n_data[n]
    a1, a2, a3, H = d['alpha1'], d['alpha2'], d['alpha3'], d['H']
    ratio_12 = a1/(2*a2)
    ratio_32 = a3/a2
    # Find dominant k among known terms
    terms = {k: 2**k * d[f'alpha{k}'] for k in range(1, 8) if f'alpha{k}' in d}
    dom_k = max(terms, key=terms.get)
    print(f"{n:>4} {H:>20,} {a1:>15,} {a2:>15,} {a3:>12,} {ratio_12:>9.4f} {ratio_32:>8.4f} {dom_k:>6}")

print()
print("NOTE: dom_k = index k with largest 2^k·αₖ contribution to H")
print("      α₁/(2α₂) < 0.5 means 4α₂ > 2α₁ (α₂ term dominates α₁ term)")
print("      α₃/α₂   > 0.5 means 8α₃ > 4α₂ (α₃ term dominates α₂ term)")

# ─── Extrapolation of crossover n ────────────────────────────────────────────
print()
print("=" * 76)
print("CROSSOVER PREDICTION — Analytic extrapolation")
print("=" * 76)

# Fit α₁/(2α₂) vs n using known data
# n=17: 0.5570, n=19: 0.5177, n=21: 0.4878, n=23: 0.4640
# Crossover at ratio=0.5: between n=19 and n=21

r12_data = [(17, 0.557042), (19, 0.51774), (21, 0.487848), (23, 0.463974)]
print("\nFirst crossover (4α₂ becomes dominant over 2α₁):")
print("  α₁/(2α₂) = 0.5 means the α₂ term surpasses the α₁ term")
for n, r in r12_data:
    flag = "← crossover region" if r < 0.5 else ""
    print(f"  n={n}: α₁/(2α₂) = {r:.4f} {flag}")

# Fit linear model to estimate crossover
# Use last two bracketing points
n1, r1 = 19, 0.51774
n2, r2 = 21, 0.487848
n_cross_12 = n1 + (n2 - n1) * (r1 - 0.5) / (r1 - r2)
print(f"  → Linear interpolation: crossover at n ≈ {n_cross_12:.1f}")

# Fit α₃/α₂ vs n
r32_data = [(17, 0.309001), (19, 0.35182), (21, 0.388993), (23, 0.422044)]
print("\nSecond crossover (8α₃ becomes dominant over 4α₂):")
print("  α₃/α₂ = 0.5 means the α₃ term surpasses the α₂ term")
for n, r in r32_data:
    print(f"  n={n}: α₃/α₂ = {r:.4f}")

# Fit linear and check increment deceleration
increments = [r32_data[i+1][1] - r32_data[i][1] for i in range(len(r32_data)-1)]
print(f"  Increments: {[f'{x:.4f}' for x in increments]}")

# Simple linear extrapolation from n=21 to n=23
slope = (r32_data[-1][1] - r32_data[-2][1]) / 2  # per unit n
n_last, r_last = r32_data[-1]

# If increments decelerating by ~0.003 per step (2 units):
# r(n+2) ≈ r(n) + Δ, with Δ decreasing
# Use geometric decay: Δ_k = Δ_0 * q^k
delta_ns = [r32_data[i+1][0]-r32_data[i][0] for i in range(len(r32_data)-1)]  # all 2
delta_rs = increments

# Fit: r(n) = r_inf - (r_inf - r0) * exp(-lambda * n)
import numpy as np
ns = np.array([x[0] for x in r32_data])
rs = np.array([x[1] for x in r32_data])
# Try r(n) = A - B * exp(-c*n)
# From data, looks like approaching ~0.7-0.8 asymptotically
# Use least squares with A, B, c

from scipy.optimize import curve_fit

def model(n, A, B, c):
    return A - B * np.exp(-c * n)

try:
    popt, _ = curve_fit(model, ns, rs, p0=[0.7, 5.0, 0.08], maxfev=10000)
    A, B, c = popt
    n_cross_32 = -np.log((A - 0.5) / B) / c
    print(f"  Fitted: r(n) = {A:.4f} - {B:.4f}·exp(-{c:.4f}·n)")
    print(f"  Asymptote: r(∞) = {A:.4f}")
    print(f"  → Predicted crossover at n ≈ {n_cross_32:.1f}")

    # Predictions
    for n_pred in [25, 27, 29, 31, 33, 35, 41]:
        r_pred = model(n_pred, A, B, c)
        print(f"  n={n_pred}: predicted α₃/α₂ ≈ {r_pred:.4f}")
except Exception as e:
    print(f"  Fit failed: {e}")
    # Fallback linear
    slope_32 = (r32_data[-1][1] - r32_data[0][1]) / (r32_data[-1][0] - r32_data[0][0])
    n_cross_32 = r32_data[-1][0] + (0.5 - r32_data[-1][1]) / slope_32
    print(f"  → Linear extrapolation: crossover at n ≈ {n_cross_32:.1f}")

# α₄/α₃ crossover
r43_data = [
    (17, 45997104/458011858),
    (19, 5521030944/42960731622),
    (21, 738531326288/4796354751404),
    (23, 111796734828336/632921002322216),
]
print("\nThird crossover (16α₄ becomes dominant over 8α₃):")
print("  α₄/α₃ = 0.5 means the α₄ term surpasses the α₃ term")
for n, r in r43_data:
    print(f"  n={n}: α₄/α₃ = {r:.4f}")

try:
    ns43 = np.array([x[0] for x in r43_data])
    rs43 = np.array([x[1] for x in r43_data])
    popt43, _ = curve_fit(model, ns43, rs43, p0=[0.5, 5.0, 0.15], maxfev=10000)
    A43, B43, c43 = popt43
    print(f"  Fitted: r(n) = {A43:.4f} - {B43:.4f}·exp(-{c43:.4f}·n)")
    print(f"  Asymptote: r(∞) = {A43:.4f}")
    if A43 > 0.5:
        n_cross_43 = -np.log((A43 - 0.5) / B43) / c43
        print(f"  → Predicted crossover at n ≈ {n_cross_43:.1f}")
    else:
        print(f"  → No crossover predicted (asymptote {A43:.4f} < 0.5)")
except Exception as e:
    print(f"  Fit failed: {e}")

print()
print("=" * 76)
print("PALEY vs INTERVAL: Why Paley loses at large n")
print("=" * 76)
# At n=7: Paley wins because 2·Δα₁ > deficit from α₂
# At n=11: Paley wins but margin shrinking
# At n=19: Interval wins
# At n=23: Need Paley data

paley_data = {
    7: {'alpha1': 80, 'alpha2': 7, 'H': 189},
    11: {'alpha1': 21169, 'alpha2': 10879, 'alpha3': 1155, 'H': 95095},
}
interval_data = {
    7: {'alpha1': 59, 'alpha2': 14, 'H': 175},
    11: {'alpha1': 18397, 'alpha2': 11110, 'alpha3': 1474, 'H': 93027},
    17: {'alpha1': 1651334601, 'H': 13689269499},
    19: {'alpha1': 126443605257, 'alpha2': 122111579294, 'H': 1184212824763},
}

print("\nPaley primes p ≡ 3 mod 4:")
print(f"{'p':>5} {'H(Paley)':>20} {'H(Int)':>20} {'Winner':>10} {'ΔH':>15}")
print("-" * 75)

p_series = [(7, 189, 175), (11, 95095, 93027), (19, 1172695746915, 1184212824763)]
for p, HP, HI in p_series:
    winner = "Paley" if HP > HI else "Interval"
    print(f"{p:>5} {HP:>20,} {HI:>20,} {winner:>10} {HP-HI:>+15,}")

print()
print("At p=23: H(Interval)=16,011,537,490,557,279; H(Paley) = COMPUTING...")
print("  [Run: python3 alpha_from_cc_bin.py 23 paley]")

print()
print("KEY INSIGHT: Paley has MORE α₁ but interval has MORE α₂,α₃,...")
print("Paley wins when 2·Δα₁ > 4·|Δα₂| + 8·|Δα₃| + ...")
print("Interval wins when the accumulated higher-order surplus exceeds the α₁ lead.")
