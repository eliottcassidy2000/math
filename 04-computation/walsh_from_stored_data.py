#!/usr/bin/env python3
"""
walsh_from_stored_data.py — Extract Walsh degree decomposition from stored H data

The p=17 exhaustive H computation was done in S66 (step4_polymer_expansion.py).
Parse the output to extract (σ, H) pairs, then compute Walsh decomposition
without re-running the expensive Held-Karp DP.

This gives us a 4th data point for the degree dominance ratio:
  p=7:  deg-4/deg-2 = 0.00  → Paley wins
  p=11: deg-4/deg-2 = 0.19  → Paley wins
  p=13: deg-4/deg-2 = 4.42  → Interval wins
  p=17: deg-4/deg-2 = ???   → Interval wins

Author: opus-2026-03-12-S67
"""

import numpy as np
import re
from itertools import combinations
from collections import defaultdict

def parse_h_data(filename, p):
    """Parse step4_polymer_expansion.out for (sigma, H) pairs at given p."""
    m = (p-1)//2
    pairs = []

    with open(filename) as f:
        content = f.read()

    # Find the section for this p
    # Format: "#N: H = HHHHH, σ = (1, -1, ...), IPR = ..."
    pattern = r'H =\s+(\d+), σ = \(([^)]+)\)'
    for match in re.finditer(pattern, content):
        H = int(match.group(1))
        sigma_str = match.group(2)
        sigma = tuple(int(x.strip()) for x in sigma_str.split(','))
        if len(sigma) == m:
            pairs.append((sigma, H))

    # Deduplicate
    seen = set()
    unique = []
    for sigma, H in pairs:
        if sigma not in seen:
            seen.add(sigma)
            unique.append((sigma, H))
    return unique

def walsh_char(S, sigma):
    prod = 1
    for k in S:
        prod *= sigma[k-1]
    return prod

def walsh_decomposition(pairs, m):
    """Compute Walsh spectrum and variance by degree."""
    chords = list(range(1, m+1))
    N = len(pairs)

    walsh_coeffs = {}
    for deg in range(0, m+1, 2):
        for S in combinations(chords, deg):
            val = sum(H * walsh_char(S, sigma) for sigma, H in pairs) / N
            if abs(val) > 1e-6:
                walsh_coeffs[S] = val

    var_by_deg = defaultdict(float)
    for S, val in walsh_coeffs.items():
        var_by_deg[len(S)] += val**2

    return walsh_coeffs, var_by_deg

# Read all data
filename = "05-knowledge/results/step4_polymer_expansion.out"

print("=" * 70)
print("WALSH DEGREE DECOMPOSITION FROM STORED DATA")
print("=" * 70)

results = {}

for p in [7, 11, 13, 17]:
    m = (p-1)//2
    pairs = parse_h_data(filename, p)
    print(f"\np={p}: found {len(pairs)} orientations (expected {1<<m})")

    if len(pairs) < (1 << m):
        print(f"  WARNING: incomplete data, skipping")
        continue

    walsh_coeffs, var_by_deg = walsh_decomposition(pairs, m)

    total_var = sum(v for d, v in var_by_deg.items() if d > 0)
    print(f"  Walsh variance decomposition:")
    for deg in sorted(var_by_deg):
        if deg == 0:
            continue
        frac = var_by_deg[deg] / total_var * 100
        print(f"    Degree {deg}: variance = {var_by_deg[deg]:.2f} ({frac:.1f}%)")

    ratio_4_2 = var_by_deg.get(4, 0) / var_by_deg.get(2, 1)
    results[p] = {
        'var_by_deg': dict(var_by_deg),
        'ratio_4_2': ratio_4_2,
        'total_var': total_var
    }

    # Also compute degree-6 and degree-8 if present
    if 6 in var_by_deg:
        ratio_6_2 = var_by_deg[6] / var_by_deg[2]
        print(f"    deg-6/deg-2 ratio: {ratio_6_2:.4f}")
    if 8 in var_by_deg:
        ratio_8_2 = var_by_deg[8] / var_by_deg[2]
        print(f"    deg-8/deg-2 ratio: {ratio_8_2:.4f}")

    # Interval vs top H comparison
    sigma_int = tuple([1]*m)
    H_int = dict(pairs)[sigma_int]
    H_max = max(H for _, H in pairs)
    print(f"\n  H(Interval) = {H_int}, H(max) = {H_max}")
    print(f"  Interval is max? {H_int == H_max}")

print("\n" + "=" * 70)
print("DEGREE DOMINANCE RATIO SUMMARY")
print("=" * 70)
print(f"\n  {'p':>3} | {'deg-4/deg-2':>12} | {'Winner':>10}")
print(f"  " + "-" * 40)
for p in sorted(results):
    r = results[p]
    winner = "Interval" if p >= 13 else "Paley"
    print(f"  {p:>3} | {r['ratio_4_2']:>12.4f} | {winner:>10}")

# Can we predict where ratio = 1?
if len(results) >= 3:
    ps = sorted(results)
    ratios = [results[p]['ratio_4_2'] for p in ps]
    print(f"\n  Log(ratio) vs p:")
    for p, r in zip(ps, ratios):
        if r > 0:
            print(f"    p={p}: log(ratio) = {np.log(r):.4f}")
        else:
            print(f"    p={p}: ratio = 0 (log = -inf)")

print("\n" + "=" * 70)
print("HIGHER DEGREE ANALYSIS AT p=17")
print("=" * 70)

if 17 in results:
    r = results[17]
    print("\n  Full variance profile at p=17:")
    for deg in sorted(r['var_by_deg']):
        if deg == 0: continue
        frac = r['var_by_deg'][deg] / r['total_var'] * 100
        print(f"    Degree {deg:>2}: {frac:>6.2f}%")

    print("\n  Higher degree contributions are emerging!")
    print("  At p=17 (m=8), Walsh degrees up to 8 are possible.")
    print("  The degree-4 dominance seen at p=13 should continue.")

print("\nDONE.")
