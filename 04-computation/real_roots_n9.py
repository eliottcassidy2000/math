#!/usr/bin/env python3
"""
OPEN-Q-015: Real-rootedness of I(Omega(T), x) at n=9.

At n=9, Omega(T) can have claws (~86% of random tournaments).
Chudnovsky-Seymour no longer applies. Yet real-rootedness still seems
to hold empirically. This script tests exhaustively enough to be convincing
and looks for structural patterns.

Key questions:
1. Does I(Omega(T), x) have all real roots at n=9?
2. What is the typical degree of I(Omega(T), x) at n=9? (max = floor(9/3) = 3)
3. Can we bound the discriminant using tournament structure?

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament, find_odd_cycles, conflict_graph
import numpy as np
from collections import defaultdict

def indep_poly_coefficients(adj):
    """Compute independence polynomial coefficients [a_0, a_1, ..., a_k]."""
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    # Count independent sets by size
    max_indep = 0
    counts = defaultdict(int)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            size = bin(mask).count('1')
            counts[size] += 1
            max_indep = max(max_indep, size)
    return [counts.get(k, 0) for k in range(max_indep + 1)]

def check_real_roots(coeffs):
    """Check if polynomial has all real roots. Returns (all_real, roots)."""
    if len(coeffs) <= 1:
        return True, []
    # numpy roots expects highest degree first
    p = list(reversed(coeffs))
    if len(p) <= 1:
        return True, []
    roots = np.roots(p)
    all_real = all(abs(r.imag) < 1e-8 for r in roots)
    return all_real, sorted([r.real for r in roots])

def has_claw(adj):
    """Check if graph has K_{1,3} (claw) as induced subgraph."""
    m = len(adj)
    for v in range(m):
        # Find non-neighbors of v
        non_nbrs = [u for u in range(m) if u != v and not adj[v][u]]
        # Need 3 pairwise non-adjacent non-neighbors of v
        # But wait, claw = v adjacent to 3 vertices that are pairwise NON-adjacent
        nbrs = [u for u in range(m) if u != v and adj[v][u]]
        if len(nbrs) < 3:
            continue
        for i in range(len(nbrs)):
            for j in range(i+1, len(nbrs)):
                for k in range(j+1, len(nbrs)):
                    a, b, c = nbrs[i], nbrs[j], nbrs[k]
                    if not adj[a][b] and not adj[a][c] and not adj[b][c]:
                        return True
    return False

print("=" * 70)
print("REAL-ROOTEDNESS OF I(Omega(T), x) AT n=9")
print("=" * 70)

# Test at n=7,8,9 with increasing sample sizes
for n in [7, 8, 9]:
    num_samples = 500 if n <= 8 else 300
    failures = 0
    claw_count = 0
    degree_dist = defaultdict(int)
    all_coeffs = []
    discriminant_stats = []

    for trial in range(num_samples):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        if not cycles:
            degree_dist[0] += 1
            continue

        cg = conflict_graph(cycles)
        coeffs = indep_poly_coefficients(cg)
        deg = len(coeffs) - 1
        degree_dist[deg] += 1
        all_coeffs.append(coeffs)

        is_claw = has_claw(cg)
        if is_claw:
            claw_count += 1

        all_real, roots = check_real_roots(coeffs)
        if not all_real:
            failures += 1
            print(f"  FAILURE at n={n}, trial {trial}: coeffs={coeffs}")
            print(f"    roots={roots}")
            print(f"    claw={is_claw}, #cycles={len(cycles)}")

        # Discriminant for degree 2
        if deg == 2:
            a0, a1, a2 = coeffs
            disc = a1**2 - 4*a0*a2
            discriminant_stats.append(disc)

        # For degree 3: discriminant
        if deg == 3:
            a0, a1, a2, a3 = coeffs
            # Discriminant of a0 + a1*x + a2*x^2 + a3*x^3
            disc = 18*a0*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a0*a2**3 - 27*a0**2*a3**2
            discriminant_stats.append(disc)

    print(f"\n  n={n}: {num_samples} random tournaments")
    print(f"    Real-root failures: {failures}")
    print(f"    Claw count: {claw_count}/{num_samples} ({100*claw_count/num_samples:.1f}%)")
    print(f"    Degree distribution: {dict(sorted(degree_dist.items()))}")

    if discriminant_stats:
        min_disc = min(discriminant_stats)
        max_disc = max(discriminant_stats)
        neg_count = sum(1 for d in discriminant_stats if d < 0)
        print(f"    Discriminant: min={min_disc}, max={max_disc}, negative={neg_count}")

# Detailed n=9 analysis: look at the structure of degree-3 polynomials
print(f"\n{'='*70}")
print("DETAILED n=9 ANALYSIS: Degree-3 independence polynomials")
print("=" * 70)

n = 9
deg3_polys = []
for trial in range(500):
    T = random_tournament(n)
    cycles = find_odd_cycles(T)
    if not cycles:
        continue
    cg = conflict_graph(cycles)
    coeffs = indep_poly_coefficients(cg)
    if len(coeffs) - 1 == 3:
        deg3_polys.append(coeffs)

print(f"  Found {len(deg3_polys)} degree-3 polynomials")
if deg3_polys:
    a1_vals = [c[1] for c in deg3_polys]
    a2_vals = [c[2] for c in deg3_polys]
    a3_vals = [c[3] for c in deg3_polys]
    print(f"  a1 range: [{min(a1_vals)}, {max(a1_vals)}]")
    print(f"  a2 range: [{min(a2_vals)}, {max(a2_vals)}]")
    print(f"  a3 range: [{min(a3_vals)}, {max(a3_vals)}]")

    # Check log-concavity: a1^2 >= a0*a2 and a2^2 >= a1*a3
    lc1_fails = sum(1 for c in deg3_polys if c[1]**2 < c[0]*c[2])
    lc2_fails = sum(1 for c in deg3_polys if c[2]**2 < c[1]*c[3])
    print(f"  Log-concavity a1^2 >= a0*a2 failures: {lc1_fails}")
    print(f"  Log-concavity a2^2 >= a1*a3 failures: {lc2_fails}")

    # Ratio analysis: for real roots of 1 + a1*x + a2*x^2 + a3*x^3,
    # Newton's inequalities give a1^2 >= 2*a2 and a2^2 >= 3*a1*a3 (for normalized)
    # More precisely: for all-real-roots polynomial with positive coefficients,
    # the ratios a_k^2 / (a_{k-1}*a_{k+1}) >= (k+1)/(d-k+1) * d/(d+1) ... complicated
    # Simpler: for all-negative-real-roots, (a_k/C(d,k))^2 >= (a_{k-1}/C(d,k-1))*(a_{k+1}/C(d,k+1))
    # (ultra-log-concavity)
    ulc_fails = 0
    for c in deg3_polys:
        # Check (a1/3)^2 >= (a0/1)*(a2/3) = a2/3
        # i.e., a1^2/9 >= a2/3, i.e., a1^2 >= 3*a2
        if c[1]**2 < 3*c[2]:
            ulc_fails += 1
    print(f"  Ultra-log-concavity (a1^2 >= 3*a2) failures: {ulc_fails}")

    # Check strongest condition: discriminant
    all_disc = []
    for c in deg3_polys:
        a0, a1, a2, a3 = c
        disc = 18*a0*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a0*a2**3 - 27*a0**2*a3**2
        all_disc.append(disc)
    neg_disc = sum(1 for d in all_disc if d < 0)
    print(f"  Negative discriminant (=> complex roots): {neg_disc}/{len(all_disc)}")
    if all_disc:
        print(f"  Discriminant range: [{min(all_disc)}, {max(all_disc)}]")
        # Show some examples
        for c, d in sorted(zip(deg3_polys, all_disc), key=lambda x: x[1])[:5]:
            print(f"    coeffs={c}, disc={d}")

# Test at even larger n
print(f"\n{'='*70}")
print("LARGER n TESTS")
print("=" * 70)

for n in [10, 12, 15]:
    num_samples = 100 if n <= 12 else 50
    failures = 0
    degree_dist = defaultdict(int)

    for trial in range(num_samples):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        if not cycles:
            degree_dist[0] += 1
            continue

        cg = conflict_graph(cycles)
        coeffs = indep_poly_coefficients(cg)
        deg = len(coeffs) - 1
        degree_dist[deg] += 1

        all_real, roots = check_real_roots(coeffs)
        if not all_real:
            failures += 1
            print(f"  FAILURE at n={n}, trial {trial}: coeffs={coeffs}")

    print(f"\n  n={n}: {num_samples} random tournaments")
    print(f"    Real-root failures: {failures}")
    print(f"    Degree distribution: {dict(sorted(degree_dist.items()))}")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
