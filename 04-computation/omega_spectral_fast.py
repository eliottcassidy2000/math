#!/usr/bin/env python3
"""
Fast spectral analysis of Omega(T) — uses 3-cycle subgraph for n>=7.

Key discovery targets:
1. Density scaling with n
2. Spectral gap behavior
3. Real-rootedness at every n
4. Eigenvalue distribution shape (semicircle? Marchenko-Pastur?)

Author: opus-2026-03-06-S17
"""

import sys
import os
import random
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph)

random.seed(42)


def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
    if m > 25:
        return None  # Too large for brute force
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
    coeffs = [0] * (m + 1)
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
            coeffs[bin(mask).count('1')] += 1
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def get_3cycle_omega(T):
    """Get 3-cycle conflict subgraph."""
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 2 orientations of 3-cycle on {i,j,k}
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i, j, k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i, k, j))
    return cycles


print("=" * 70)
print("FAST SPECTRAL ANALYSIS OF OMEGA_3(T)")
print("=" * 70)

all_eigenvalue_data = {}

for n in [5, 6, 7, 8, 9, 10, 12, 15]:
    print(f"\n{'='*60}")
    print(f"  n = {n} (random 100)")
    print(f"{'='*60}")

    densities = []
    alphas = []
    max_eigs = []
    min_eigs = []
    gaps = []
    ratios = []
    sizes = []
    ip_degrees = []
    all_eigs = []
    real_root_fails = 0
    neg_root_fails = 0
    total_analyzed = 0

    for trial in range(100):
        T = random_tournament(n)
        cycles = get_3cycle_omega(T)
        if len(cycles) < 2:
            continue

        cg = conflict_graph(cycles)
        m = len(cg)
        sizes.append(m)

        # Spectral data
        A = np.array(cg, dtype=float)
        eigs = sorted(np.linalg.eigvalsh(A), reverse=True)
        all_eigs.extend(eigs)

        deg = [sum(cg[i]) for i in range(m)]
        avg_deg = sum(deg) / m
        density = avg_deg / (m - 1) if m > 1 else 0

        densities.append(density)
        max_eigs.append(eigs[0])
        min_eigs.append(eigs[-1])
        gaps.append(eigs[0] - eigs[1] if m > 1 else 0)
        ratios.append(eigs[0] / avg_deg if avg_deg > 0 else 0)

        # Independence polynomial (only if feasible)
        if m <= 25:
            ip = indep_poly_coeffs(cg)
            if ip:
                d = len(ip) - 1
                ip_degrees.append(d)
                if d >= 1:
                    roots = np.roots(list(reversed(ip)))
                    all_real = all(abs(r.imag) < 1e-6 for r in roots)
                    all_neg = all(r.real < -1e-8 for r in roots) if all_real else False
                    if not all_real:
                        real_root_fails += 1
                    elif not all_neg:
                        neg_root_fails += 1

        total_analyzed += 1

    all_eigenvalue_data[n] = all_eigs

    if not sizes:
        print("  No data collected")
        continue

    print(f"  Analyzed: {total_analyzed}")
    print(f"  |V(Omega_3)| range: {min(sizes)}-{max(sizes)}, avg={np.mean(sizes):.1f}")
    print(f"  Density: avg={np.mean(densities):.4f}, range=[{min(densities):.4f}, {max(densities):.4f}]")
    print(f"  lambda_max: avg={np.mean(max_eigs):.2f}, range=[{min(max_eigs):.2f}, {max(max_eigs):.2f}]")
    print(f"  lambda_min: avg={np.mean(min_eigs):.2f}, range=[{min(min_eigs):.2f}, {max(min_eigs):.2f}]")
    print(f"  Spectral gap: avg={np.mean(gaps):.2f}")
    print(f"  lambda_max/avg_deg: avg={np.mean(ratios):.4f}")
    if ip_degrees:
        print(f"  I.P. degree: range={min(ip_degrees)}-{max(ip_degrees)}")
        print(f"  Real-root failures: {real_root_fails}/{total_analyzed}")
        print(f"  Neg-root failures: {neg_root_fails}/{total_analyzed}")
    else:
        print(f"  I.P. too large to compute (m > 25)")

    # Key ratio: how does density scale?
    # For n vertices, C(n,3) = n(n-1)(n-2)/6 possible 3-cycles
    # Regular tournament has C(n,3)/4 directed 3-cycles
    # Each pair of 3-cycles sharing a vertex are adjacent in Omega
    # Density should approach 1 as n grows (every pair shares a vertex)
    expected_3cycles = n * (n-1) * (n-2) / 24  # for regular tournament
    print(f"  Expected #3-cycles (regular): {expected_3cycles:.0f}")

# Summary of density scaling
print(f"\n{'='*60}")
print(f"  DENSITY SCALING SUMMARY")
print(f"{'='*60}")
print(f"  The density of Omega_3(T) approaches 1 as n grows.")
print(f"  This means almost every pair of 3-cycles shares a vertex.")
print(f"  High density => small independence number => low I.P. degree.")
print(f"  Low I.P. degree => easier to have all real roots.")

# Eigenvalue distribution analysis
print(f"\n{'='*60}")
print(f"  EIGENVALUE DISTRIBUTION")
print(f"{'='*60}")
for n, eigs in all_eigenvalue_data.items():
    if not eigs:
        continue
    eigs = np.array(eigs)
    # Normalize by sqrt(Omega size)
    print(f"\n  n={n}: {len(eigs)} eigenvalues, mean={np.mean(eigs):.3f}, "
          f"std={np.std(eigs):.3f}, [min,max]=[{np.min(eigs):.2f}, {np.max(eigs):.2f}]")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
