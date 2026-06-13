#!/usr/bin/env python3
"""
circulant_topology_survey.py — Survey eigenspace anomaly depth across circulant orientations at p=11

Tests whether the k=0 anomaly pattern at m=2 (Δ₂ = r₂^(0) - r₂^(k≥1)) is nonzero,
which determines whether the tournament has "deep anomaly" (Paley-like) or not.

Since checking m=2 is fast, we can survey all 32 orientations.

Key insight: If Δ₂ = 0, the tournament likely has β=[1,1,0,...,0] (Interval-like).
If Δ₂ ≠ 0, it may have high-dimensional Betti.

Author: opus-2026-03-12-S68
"""
import sys, time
from itertools import combinations
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology, find_nth_root_of_unity

p = 11
m = (p - 1) // 2  # = 5

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

# Generate all m-element subsets of {1,...,p-1} where S ∩ (p-S) = ∅
# i.e., for each pair {d, p-d}, exactly one is in S
def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = set()
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.add(a)
            else:
                S.add(b)
        yield frozenset(S)

# Classify each orientation
S_paley = frozenset(d for d in range(1, p) if legendre(d, p) == 1)
S_interval = frozenset(range(1, m + 1))

print(f"p = {p}, m = {m}")
print(f"Paley S = {sorted(S_paley)}")
print(f"Interval S = {sorted(S_interval)}")
print(f"Total orientations: {2**m}")
print()

results = []
for S in gen_orientations(p):
    h = CirculantHomology(n=p, S=S)
    h._ensure_enumerated(3)  # Only need up to degree 3
    prime_f = h.prime
    omega_p = find_nth_root_of_unity(p, prime_f)

    # Compute ranks at m=1 and m=2 for k=0 and k=1
    omega_k0 = 1
    omega_k1 = omega_p

    r1_k0 = h._boundary_rank_k(1, omega_k0)
    r1_k1 = h._boundary_rank_k(1, omega_k1)
    r2_k0 = h._boundary_rank_k(2, omega_k0)
    r2_k1 = h._boundary_rank_k(2, omega_k1)
    r3_k0 = h._boundary_rank_k(3, omega_k0)
    r3_k1 = h._boundary_rank_k(3, omega_k1)

    delta_1 = r1_k0 - r1_k1
    delta_2 = r2_k0 - r2_k1
    delta_3 = r3_k0 - r3_k1

    omega_dims = h.omega_dims(max_degree=3)

    label = ""
    if S == S_paley:
        label = " [PALEY]"
    elif S == S_interval:
        label = " [INTERVAL]"
    elif S == frozenset(p - d for d in S_paley):
        label = " [ANTI-PALEY]"
    elif S == frozenset(p - d for d in S_interval):
        label = " [ANTI-INTERVAL]"

    anomaly_depth = 0
    if delta_1 != 0: anomaly_depth = 1
    if delta_2 != 0: anomaly_depth = 2
    if delta_3 != 0: anomaly_depth = 3

    results.append({
        'S': sorted(S),
        'label': label,
        'delta': (delta_1, delta_2, delta_3),
        'depth': anomaly_depth,
        'omega': omega_dims[:4],
        'r0': (r1_k0, r2_k0, r3_k0),
        'r1': (r1_k1, r2_k1, r3_k1),
    })

# Group by anomaly depth
from collections import Counter
depth_counts = Counter(r['depth'] for r in results)
print(f"Anomaly depth distribution (from m=1-3 check):")
for d in sorted(depth_counts):
    print(f"  D >= {d}: {depth_counts[d]} orientations")

print(f"\n{'='*80}")
print(f"DETAIL: Orientations with deep anomaly (D >= 2)")
print(f"{'='*80}")
for r in sorted(results, key=lambda x: -x['depth']):
    if r['depth'] >= 2:
        print(f"  S={r['S']}{r['label']}")
        print(f"    Δ = {r['delta']}, Ω = {r['omega']}")

print(f"\n{'='*80}")
print(f"DETAIL: Orientations with D = 1 (Interval-like)")
print(f"{'='*80}")
for r in sorted(results, key=lambda x: str(x['S'])):
    if r['depth'] <= 1:
        print(f"  S={r['S']}{r['label']}")
        print(f"    Δ = {r['delta']}")
