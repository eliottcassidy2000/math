#!/usr/bin/env python3
"""
p11_topological_dichotomy.py — Extend topological dichotomy to p=11

At p=7: Paley β=[1,0,0,0,6,0,0] χ=7 vs ALL others β=[1,1,0,0,0,0,0] χ=0.
Question: Does this persist at p=11? Does the Interval (H-maximizer at p≥13)
have different topology from Paley?

Uses CirculantHomology for efficient computation via eigenspace decomposition.

Author: opus-2026-03-12-S68
"""
import sys
sys.path.insert(0, '04-computation')

from circulant_homology import CirculantHomology, PaleyHomology
import numpy as np
import time

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

print("=" * 70)
print("TOPOLOGICAL DICHOTOMY AT p=11: ALL CIRCULANT TOURNAMENTS")
print("=" * 70)

p = 11
m = (p - 1) // 2  # = 5

# Enumerate all 2^m = 32 circulant orientations
sigma_int = tuple([1]*m)
sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

print(f"\np={p}, m={m}")
print(f"Interval σ = {sigma_int}")
print(f"Paley    σ = {sigma_pal}")
print()

# Compute H values using path_homology_v2 approach (or just from stored data)
# For now, focus on Betti numbers

results = {}
for bits in range(1 << m):
    sigma = tuple(1 if (bits >> j) & 1 else -1 for j in range(m))
    S = set()
    for k in range(1, m+1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    S = frozenset(S)

    tag = ""
    if sigma == sigma_int: tag = " [INTERVAL]"
    elif sigma == sigma_pal: tag = " [PALEY]"
    elif sigma == tuple(-x for x in sigma_pal): tag = " [anti-PALEY]"
    elif sigma == tuple(-x for x in sigma_int): tag = " [anti-INTERVAL]"

    print(f"  σ={sigma}, S={sorted(S)}{tag}")

    t0 = time.time()
    h = CirculantHomology(n=p, S=S)

    # Compute Betti numbers - p=11 should be feasible
    try:
        betti = h.betti_numbers(max_degree=p-1, verbose=False)
        chi = sum((-1)**i * b for i, b in enumerate(betti))
        dt = time.time() - t0
        print(f"    β = {betti}, χ = {chi}  ({dt:.1f}s)")
        results[sigma] = (betti, chi, sorted(S))
    except Exception as e:
        dt = time.time() - t0
        print(f"    Error after {dt:.1f}s: {e}")
        # Try with lower max_degree
        try:
            betti = h.betti_numbers(max_degree=6, verbose=False)
            chi_partial = sum((-1)**i * b for i, b in enumerate(betti))
            print(f"    β (partial, deg≤6) = {betti}, partial χ = {chi_partial}")
            results[sigma] = (betti, chi_partial, sorted(S))
        except Exception as e2:
            print(f"    Also failed at deg≤6: {e2}")

    sys.stdout.flush()

# Summary
print("\n" + "=" * 70)
print("SUMMARY: TOPOLOGICAL TYPES AT p=11")
print("=" * 70)

# Group by Betti type
from collections import Counter
betti_groups = Counter()
for sigma, (betti, chi, S) in results.items():
    betti_groups[tuple(betti)] += 1

for betti_type, count in betti_groups.most_common():
    print(f"\n  β = {list(betti_type)} (χ = {sum((-1)**i * b for i, b in enumerate(betti_type))})")
    print(f"    Count: {count}/{len(results)} orientations")
    for sigma, (betti, chi, S) in results.items():
        if tuple(betti) == betti_type:
            tag = ""
            if sigma == sigma_int: tag = " ← INTERVAL"
            elif sigma == sigma_pal: tag = " ← PALEY"
            print(f"      σ={sigma}, S={S}{tag}")

print("\n" + "=" * 70)
print("KEY QUESTION: Is Paley STILL unique at p=11?")
print("=" * 70)
print("""
At p=7: Paley is UNIQUE — the only orientation with β_4=6, χ=7.
All 7 others have β_1=1, χ=0.

If p=11 shows the same pattern (Paley uniquely has high-dim homology,
while Interval has β_1=1), this is a TOPOLOGICAL characterization of
the Paley tournament that goes beyond spectral properties.

If Interval acquires new topology at p=11 (near the phase transition
p≈13), this suggests topology ALSO undergoes a phase transition.
""")

print("DONE.")
