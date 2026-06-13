#!/usr/bin/env python3
"""
Real-rootedness of I(Omega_3(T), x) at large n.

Omega_3(T) = conflict graph on 3-cycles only (fast to compute).
At n<=8, alpha(Omega_3) <= 2 so I(Omega_3,x) has degree <= 2 (trivially real-rooted).
At n=9, alpha(Omega_3) <= 3 so degree <= 3. Real-rootedness needs checking.
At n>=12, alpha(Omega_3) can be 4+ (floor(n/3) disjoint 3-cycles).

This script tests real-rootedness at n=9,...,21 where the polynomial can
have degree up to 7.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
import numpy as np
from collections import defaultdict
from itertools import combinations

def find_3cycles(T):
    """Find all directed 3-cycles. Returns list of sorted vertex triples."""
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i, j, k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i, j, k))  # vertex set, orientation doesn't matter for conflict
    return cycles

def build_conflict_graph_3(cycles_3):
    """Build conflict graph: two 3-cycles conflict iff they share a vertex."""
    m = len(cycles_3)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if set(cycles_3[i]) & set(cycles_3[j]):
                adj[i][j] = adj[j][i] = 1
    return adj

def indep_poly_coefficients(adj):
    m = len(adj)
    if m == 0:
        return [1]
    nbr = [0] * m
    for i in range(m):
        for j in range(m):
            if adj[i][j]:
                nbr[i] |= 1 << j
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
            counts[bin(mask).count('1')] += 1
    mx = max(counts.keys()) if counts else 0
    return [counts.get(k, 0) for k in range(mx + 1)]

def indep_poly_coefficients_fast(adj, max_alpha):
    """Compute independence polynomial using dynamic programming on bitmask.
    Only works when m is small enough for bitmask."""
    m = len(adj)
    if m == 0:
        return [1]
    if m > 25:
        # Too large for exact computation. Use sampling or bounds.
        return None
    return indep_poly_coefficients(adj)

def check_real_roots(coeffs):
    if len(coeffs) <= 1:
        return True, []
    p = list(reversed(coeffs))
    roots = np.roots(p)
    all_real = all(abs(r.imag) < 1e-6 for r in roots)
    return all_real, sorted([r.real for r in roots])

print("=" * 70)
print("REAL-ROOTEDNESS OF I(Omega_3(T), x) AT LARGE n")
print("=" * 70)

for n in [9, 10, 12, 15, 18, 21]:
    num_samples = 200 if n <= 12 else (100 if n <= 18 else 50)
    failures = 0
    skipped = 0
    degree_dist = defaultdict(int)
    c3_dist = []

    for trial in range(num_samples):
        T = random_tournament(n)
        cycles_3 = find_3cycles(T)
        c3 = len(cycles_3)
        c3_dist.append(c3)

        if c3 == 0:
            degree_dist[0] += 1
            continue

        if c3 > 25:
            # Too many 3-cycles for exact computation
            skipped += 1
            continue

        cg = build_conflict_graph_3(cycles_3)
        coeffs = indep_poly_coefficients(cg)
        if coeffs is None:
            skipped += 1
            continue

        deg = len(coeffs) - 1
        degree_dist[deg] += 1

        all_real, roots = check_real_roots(coeffs)
        if not all_real:
            failures += 1
            print(f"  FAILURE at n={n}, trial {trial}: coeffs={coeffs}")
            print(f"    c3={c3}, roots={[f'{r:.4f}' for r in roots]}")

    avg_c3 = sum(c3_dist) / len(c3_dist) if c3_dist else 0
    print(f"\n  n={n}: {num_samples} samples, {skipped} skipped (too large)")
    print(f"    Real-root failures: {failures}")
    print(f"    avg c3 = {avg_c3:.1f}, max c3 = {max(c3_dist) if c3_dist else 0}")
    print(f"    Degree distribution: {dict(sorted(degree_dist.items()))}")

# For small n where we can be exhaustive, also check log-concavity and
# ultra-log-concavity (which implies real-rootedness for certain classes)
print(f"\n{'='*70}")
print("LOG-CONCAVITY AND ULTRA-LOG-CONCAVITY CHECK")
print("=" * 70)

for n in [9, 10]:
    num_samples = 200
    lc_fails = 0
    ulc_fails = 0
    tested = 0

    for trial in range(num_samples):
        T = random_tournament(n)
        cycles_3 = find_3cycles(T)
        if len(cycles_3) == 0 or len(cycles_3) > 25:
            continue

        cg = build_conflict_graph_3(cycles_3)
        coeffs = indep_poly_coefficients(cg)
        if coeffs is None or len(coeffs) <= 2:
            continue

        tested += 1
        d = len(coeffs) - 1

        # Log-concavity: a_k^2 >= a_{k-1} * a_{k+1} for all k
        lc_ok = True
        for k in range(1, d):
            if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1]:
                lc_ok = False
                break
        if not lc_ok:
            lc_fails += 1

        # Ultra-log-concavity: (a_k / C(d,k))^2 >= (a_{k-1}/C(d,k-1)) * (a_{k+1}/C(d,k+1))
        from math import comb
        ulc_ok = True
        for k in range(1, d):
            lhs = (coeffs[k] / comb(d, k))**2
            rhs = (coeffs[k-1] / comb(d, k-1)) * (coeffs[k+1] / comb(d, k+1))
            if lhs < rhs - 1e-10:
                ulc_ok = False
                break
        if not ulc_ok:
            ulc_fails += 1

    print(f"  n={n}: {tested} tested")
    print(f"    Log-concavity failures: {lc_fails}")
    print(f"    Ultra-log-concavity failures: {ulc_fails}")

# Special investigation: at n=9, what do the degree-3 polynomials look like?
print(f"\n{'='*70}")
print("n=9 DEGREE-3 DEEP DIVE")
print("=" * 70)

n = 9
deg3_data = []
for trial in range(500):
    T = random_tournament(n)
    cycles_3 = find_3cycles(T)
    if len(cycles_3) == 0 or len(cycles_3) > 25:
        continue
    cg = build_conflict_graph_3(cycles_3)
    coeffs = indep_poly_coefficients(cg)
    if coeffs is not None and len(coeffs) - 1 == 3:
        deg3_data.append(coeffs)

print(f"  Collected {len(deg3_data)} degree-3 polynomials")
if deg3_data:
    # All have a0 = 1
    a1_vals = [c[1] for c in deg3_data]
    a2_vals = [c[2] for c in deg3_data]
    a3_vals = [c[3] for c in deg3_data]

    print(f"  a1: min={min(a1_vals)}, max={max(a1_vals)}, avg={sum(a1_vals)/len(a1_vals):.1f}")
    print(f"  a2: min={min(a2_vals)}, max={max(a2_vals)}, avg={sum(a2_vals)/len(a2_vals):.1f}")
    print(f"  a3: min={min(a3_vals)}, max={max(a3_vals)}, avg={sum(a3_vals)/len(a3_vals):.1f}")

    # For real roots of 1 + a1*x + a2*x^2 + a3*x^3:
    # Discriminant Delta = 18*a1*a2*a3 - 4*a1^3*a3 + a1^2*a2^2 - 4*a2^3 - 27*a3^2
    # (with a0=1)
    discs = []
    for c in deg3_data:
        a0, a1, a2, a3 = c
        disc = 18*a0*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a0*a2**3 - 27*a0**2*a3**2
        discs.append(disc)

    neg = sum(1 for d in discs if d < 0)
    print(f"  Discriminant: min={min(discs)}, max={max(discs)}")
    print(f"  Negative discriminants (complex roots): {neg}/{len(discs)}")

    # Newton's inequality for degree 3: a1^2 >= 3*a2 and a2^2 >= 3*a1*a3
    # (necessary for all-real-roots with positive coefficients)
    n1_fails = sum(1 for c in deg3_data if c[1]**2 < 3*c[2])
    n2_fails = sum(1 for c in deg3_data if c[2]**2 < 3*c[1]*c[3])
    print(f"  Newton a1^2 >= 3*a2 failures: {n1_fails}")
    print(f"  Newton a2^2 >= 3*a1*a3 failures: {n2_fails}")

    # Ratio a1^2 / (3*a2): how much slack?
    ratios = [c[1]**2 / (3*c[2]) if c[2] > 0 else float('inf') for c in deg3_data]
    print(f"  Ratio a1^2/(3*a2): min={min(ratios):.4f}, avg={sum(ratios)/len(ratios):.4f}")

    ratios2 = [c[2]**2 / (3*c[1]*c[3]) if c[1]*c[3] > 0 else float('inf') for c in deg3_data]
    print(f"  Ratio a2^2/(3*a1*a3): min={min(ratios2):.4f}, avg={sum(ratios2)/len(ratios2):.4f}")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
