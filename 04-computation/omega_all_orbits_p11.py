#!/usr/bin/env python3
"""
omega_all_orbits_p11.py — opus-2026-03-12-S69

Compute complete Ω_m profiles for all 4 orbit types at p=11.
chi_per = Σ(-1)^m Ω_m / p.

This directly determines the Euler characteristic from path enumeration alone.
"""

import sys
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

p = 11

# Orbit representatives
orbit_reps = {
    'Interval': {1,2,3,4,5},
    'Paley': {1,3,4,5,9},
    'Orbit_A': {1,2,3,4,6},
    'Orbit_B': {1,6,7,8,9},
}

topo_data = {
    'Interval': {'chi': 0, 'chi_per': 0},
    'Paley': {'chi': 11, 'chi_per': 1},
    'Orbit_A': {'chi': 22, 'chi_per': 2},
    'Orbit_B': {'chi': 0, 'chi_per': 0},
}

print(f"{'='*80}")
print(f"COMPLETE Ω PROFILES FOR ALL ORBIT TYPES AT p={p}")
print(f"{'='*80}")

for name, S in orbit_reps.items():
    h = CirculantHomology(n=p, S=frozenset(S))
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)

    chi_total = sum((-1)**m * w for m, w in enumerate(omega))
    chi_per_computed = chi_total / p

    td = topo_data[name]

    print(f"\n{'─'*60}")
    print(f"  {name}: S = {sorted(S)}")
    print(f"  Ω = {omega}")
    print(f"  |Ω| = {sum(omega)}")
    print(f"  chi = Σ(-1)^m Ω_m = {chi_total}")
    print(f"  chi_per = chi/p = {chi_per_computed}")
    print(f"  Expected chi_per = {td['chi_per']}")
    print(f"  {'✓' if abs(chi_per_computed - td['chi_per']) < 0.01 else '✗'}")

    # Per-eigenspace Ω
    omega_per = [w // p for w in omega]
    print(f"  Ω/p = {omega_per}")

    # Check which degrees contribute to chi_per
    print(f"\n  Degree contributions to chi_per:")
    for m, w in enumerate(omega):
        contrib = (-1)**m * w / p
        if abs(contrib) > 0.01:
            print(f"    m={m}: (-1)^{m} * {w}/{p} = {contrib:+.1f}")

# Cross-comparison: what makes Ω different?
print(f"\n{'='*80}")
print(f"CROSS-COMPARISON OF Ω PROFILES")
print(f"{'='*80}")

all_omega = {}
for name, S in orbit_reps.items():
    h = CirculantHomology(n=p, S=frozenset(S))
    h._ensure_enumerated(p)
    omega = h.omega_dims(max_degree=p-1)
    all_omega[name] = omega

print(f"\n  {'m':>3} {'Interval':>10} {'Paley':>10} {'Orbit_A':>10} {'Orbit_B':>10}  {'Int-Pal':>10} {'A-B':>10}")
for m in range(p):
    vals = [all_omega[name][m] if m < len(all_omega[name]) else 0 for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']]
    diff1 = vals[0] - vals[1]
    diff2 = vals[2] - vals[3]
    print(f"  {m:3d} {vals[0]:10d} {vals[1]:10d} {vals[2]:10d} {vals[3]:10d}  {diff1:+10d} {diff2:+10d}")

# Total
for name in ['Interval', 'Paley', 'Orbit_A', 'Orbit_B']:
    print(f"  Total Ω({name}) = {sum(all_omega[name])}")

# Key question: do the differences Ω_m(A) - Ω_m(B) have structure?
print(f"\n{'='*80}")
print(f"DIFFERENCES IN Ω — WHY ORBIT_A ≠ ORBIT_B?")
print(f"{'='*80}")

print(f"\n  Orbit_A - Orbit_B (both have S-S = Z_p*, prod Q = 23):")
for m in range(p):
    diff = all_omega['Orbit_A'][m] - all_omega['Orbit_B'][m]
    if diff != 0:
        print(f"    m={m}: Δ = {diff:+d} (A={all_omega['Orbit_A'][m]}, B={all_omega['Orbit_B'][m]})")

print(f"\n  Paley - Interval:")
for m in range(p):
    diff = all_omega['Paley'][m] - all_omega['Interval'][m]
    if diff != 0:
        print(f"    m={m}: Δ = {diff:+d} (P={all_omega['Paley'][m]}, I={all_omega['Interval'][m]})")

# Now try to connect Ω_m differences to the Q-spectrum
print(f"\n{'='*80}")
print(f"CONNECTION: Ω_m FROM DIFFERENCE SEQUENCES")
print(f"{'='*80}")
print(f"  Ω_m = # of m-tuples (d_1,...,d_m) from S with all prefix sums distinct nonzero")
print(f"  This is a PURELY COMBINATORIAL property of the set S mod p")
print()

# Compute A_m counts (number of allowed diff sequences) directly
for name, S in orbit_reps.items():
    S_sorted = sorted(S)
    # Direct enumeration of m-paths for small m
    print(f"  {name}: S = {S_sorted}")

    # Additive structure
    from itertools import product as cart_product

    # How many 2-paths? (d1, d2) ∈ S² with d1, d1+d2 distinct nonzero and ≠ each other
    count2 = 0
    for d1 in S:
        for d2 in S:
            s1 = d1 % p
            s2 = (d1 + d2) % p
            if s1 != 0 and s2 != 0 and s1 != s2:
                count2 += 1
    print(f"    Ω_2 = {count2} (via direct enum), expected = {all_omega[name][2]}")

    # Sumset structure: S+S mod p
    sumset = set()
    for a in S:
        for b in S:
            sumset.add((a + b) % p)
    sumset.discard(0)  # Remove 0 if present
    print(f"    |S+S| = {len(sumset)} (out of {p-1} nonzero residues)")

    # 2S (doubling set)
    double_set = set((2*s) % p for s in S)
    print(f"    |2S| = {len(double_set)}")

    # Higher sumsets
    S3 = set()
    for a in S:
        for b in S:
            for c in S:
                S3.add((a+b+c) % p)
    S3.discard(0)
    print(f"    |S+S+S| = {len(S3)}")

print("\nDONE.")
