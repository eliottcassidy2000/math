#!/usr/bin/env python3
"""
H_chi_bridge.py — opus-2026-03-12-S69

Bridge between Hamilton path count H and GLMY Euler characteristic chi_per
for circulant tournaments.

Both H and chi_per are determined by the Q-spectrum (THM-145 for chi_per).
Question: Is there a direct relationship H = f(chi_per, p, other invariants)?

For the Interval:
  H/p is related to Morgan-Voyce polynomials
  chi_per = 0 at p≡3 mod 4, 1 at p≡1 mod 4

For Paley (p≡3 mod 4):
  Q_k = (p+1)/4 for all k (flat spectrum)
  chi_per = 1
"""

import sys
import numpy as np
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        yield S

def zpstar_orbit_rep(p, S):
    best = sorted(S)
    for a in range(2, p):
        S_a = frozenset((a * s) % p for s in S)
        if sorted(S_a) < best:
            best = sorted(S_a)
    return tuple(best)

def count_hamiltonian_paths(p, S):
    """Count Hamilton paths in circulant tournament C_p^S."""
    # Use the adjacency matrix approach
    # For small p, enumerate directly
    adj = {}
    for v in range(p):
        adj[v] = set()
        for s in S:
            adj[v].add((v + s) % p)

    count = 0
    # Try all starting vertices (by symmetry, multiply by p/p = 1... no, need to count all)
    # Actually for circulant: all starting vertices give same count, so count from 0 and multiply by p
    # But we want TOTAL H, so count from 0 and multiply by p

    def backtrack(path, visited):
        nonlocal count
        if len(path) == p:
            count += 1
            return
        v = path[-1]
        for u in adj[v]:
            if u not in visited:
                visited.add(u)
                path.append(u)
                backtrack(path, visited)
                path.pop()
                visited.remove(u)

    # Count from vertex 0 only (multiply by p at end)
    backtrack([0], {0})
    return count  # This is H/p per starting vertex... actually this IS H_from_0

def compute_Q(p, S):
    m = (p-1)//2
    return np.array([abs(sum(np.exp(2j*np.pi*s*k/p) for s in S))**2 for k in range(1, m+1)])

def compute_esf(Q):
    m = len(Q)
    coeffs = np.polynomial.polynomial.polyfromroots(Q)
    return [round(float(np.real(coeffs[m-j] * (-1)**j))) for j in range(m+1)]

print("="*80)
print("H vs chi_per BRIDGE FOR CIRCULANT TOURNAMENTS")
print("="*80)

for p in [7, 11]:
    m = (p-1)//2
    print(f"\np = {p}, m = {m}")
    print(f"  {'S':20s} | {'e_j (non-univ)':25s} | {'chi_per':>8s} | {'H_from_0':>10s} | {'H_total':>10s} | {'H/p':>8s}")

    orbits = {}
    for S in gen_orientations(p):
        rep = zpstar_orbit_rep(p, S)
        if rep not in orbits:
            orbits[rep] = S

    for rep, S in sorted(orbits.items()):
        # Compute chi_per
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))

        # Compute H
        H_from_0 = count_hamiltonian_paths(p, S)
        H_total = H_from_0 * p  # By circulant symmetry

        # Compute Q/ESF
        Q = compute_Q(p, S)
        esf = compute_esf(Q)
        esf_key = esf[2:]

        print(f"  {str(sorted(S)):20s} | {str(esf_key):25s} | {chi_per:8d} | {H_from_0:10d} | {H_total:10d} | {H_total//p:8d}")

    # Compute some relationships
    print(f"\n  Relationships:")
    results = []
    for rep, S in sorted(orbits.items()):
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))
        H_from_0 = count_hamiltonian_paths(p, S)
        H_total = H_from_0 * p
        total_omega = sum(omega)
        Q = compute_Q(p, S)
        prod_1_Q = np.prod(1 + Q)
        sum_Q2 = sum(Q**2)

        results.append({
            'S': sorted(S), 'chi_per': chi_per, 'H': H_total,
            'total_omega': total_omega, 'prod_1_Q': prod_1_Q,
            'sum_Q2': sum_Q2
        })

        print(f"  S={sorted(S)}: H={H_total}, Σ Ω = {total_omega}, χ_per={chi_per}, "
              f"Π(1+Q)={prod_1_Q:.1f}, H/(p·Π(1+Q))={H_total/(p*prod_1_Q):.4f}")

    # Check: is H/p related to Σ Ω?
    print(f"\n  H/p vs Σ Ω:")
    for r in results:
        print(f"    S={r['S']}: H/p = {r['H']//p}, Σ Ω = {r['total_omega']}, ratio = {r['H']/(p*r['total_omega']):.4f}")

# For p=13, skip H computation (too slow), just show chi_per data
print(f"\n{'='*80}")
print(f"p = 13: chi_per ALREADY COMPUTED (H too expensive)")
print(f"{'='*80}")

# From thm145_p13_verify results
p13_data = [
    ([7,8,9,10,11,12], (70,84,45,11,1), 1),     # Interval
    ([1,2,7,8,9,10], (122,305,357,180,27), -3),
    ([1,2,4,7,8,10], (161,552,812,375,53), 4),
    ([1,3,7,8,9,11], (174,721,1566,1701,729), -17),  # quasi-Paley
    ([1,7,8,9,10,11], (148,448,604,362,79), -8),
    ([2,7,8,9,10,12], (161,552,851,570,131), 1),
]

print(f"\n  {'S':25s} | {'e_m (prod Q)':>12s} | {'chi_per':>8s} | {'chi_per * e_m':>14s}")
for S, esf, chi in p13_data:
    prod_Q = esf[-1]
    print(f"  {str(S):25s} | {prod_Q:12d} | {chi:8d} | {chi * prod_Q:14d}")

# Interesting: chi_per * prod(Q) for quasi-Paley = -17 * 729 = -12393
# For Interval: 1 * 1 = 1

# Check: sum of chi_per * orbit_size
print(f"\n  Sum of chi_per weighted by orbit size:")
# At p=13, Z_p*-orbits have sizes:
# Need to compute orbit sizes
from math import gcd

p = 13
m = (p-1)//2
orbit_sizes = {}
all_orbits = {}
for S in gen_orientations(p):
    rep = zpstar_orbit_rep(p, S)
    if rep not in all_orbits:
        all_orbits[rep] = 0
    all_orbits[rep] += 1

# Match to our data
chi_pers = {
    tuple([1,2,3,4,5,6]): 1,   # = [7,8,...,12] complement
    tuple([1,2,3,4,5,7]): -3,
    tuple([1,2,4,7,8,10]): 4,
    tuple([1,2,3,5,6,9]): -17,
    tuple([1,2,3,5,7,9]): -8,
    tuple([1,2,6,8,9,10]): 1,
}

total = 0
for rep, size in sorted(all_orbits.items()):
    chi = chi_pers.get(rep, '?')
    print(f"    rep={list(rep)}: orbit size = {size}, chi_per = {chi}")
    if chi != '?':
        total += chi * size

print(f"\n  Σ chi_per * |orbit| = {total}")
print(f"  Total orientations = {2**m} = {2**m}")
print(f"  Average chi_per = {total}/{2**m} = {total/(2**m):.4f}")

# KEY: Σ chi_per * |orbit| over all Z_p*-orbits = what?
# This equals Σ_T chi_per(T) where sum is over ALL orientations.
# By Burnside or character theory, this might simplify.

print("\nDONE.")
