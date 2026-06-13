#!/usr/bin/env python3
"""
Spectral analysis of Omega(T) — the odd-cycle conflict graph.

Investigates whether the adjacency spectrum of Omega(T) has special structure
that could explain the persistent real-rootedness of I(Omega(T), x).

Key questions:
1. Is Omega(T) always a "high-energy" graph (large spectral radius)?
2. Does the ratio lambda_max / n grow in a specific way?
3. Are there spectral gap properties?
4. How does the Lovász theta function relate to independence number?
5. Does the Gallai graph G(Omega) have special structure?

The connection: for claw-free graphs, Chudnovsky-Seymour works via
properties of the neighborhood of each vertex (quasi-line structure).
Even when Omega has claws (n>=9), the NEIGHBORHOOD structure may retain
enough regularity for real roots.

Author: opus-2026-03-06-S17
"""

import sys
import os
import random
import itertools
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, random_tournament,
                             find_odd_cycles, conflict_graph)

random.seed(42)

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("WARNING: numpy not available, skipping spectral analysis")


def indep_poly_coeffs(adj):
    m = len(adj)
    if m == 0:
        return [1]
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


def graph_properties(adj):
    """Compute key properties of a graph given as adjacency matrix."""
    m = len(adj)
    if m == 0:
        return {}

    props = {}
    props['vertices'] = m

    # Edge count
    edges = sum(adj[i][j] for i in range(m) for j in range(i+1, m))
    props['edges'] = edges
    props['density'] = 2 * edges / (m * (m - 1)) if m > 1 else 0

    # Degree sequence
    degrees = [sum(adj[i]) for i in range(m)]
    props['min_degree'] = min(degrees)
    props['max_degree'] = max(degrees)
    props['avg_degree'] = sum(degrees) / m

    # Independence number (brute force for small m)
    if m <= 20:
        alpha = 0
        nbr = [0] * m
        for i in range(m):
            for j in range(m):
                if adj[i][j]:
                    nbr[i] |= 1 << j
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
                alpha = max(alpha, bin(mask).count('1'))
        props['alpha'] = alpha

    # Spectral properties
    if HAS_NUMPY:
        A = np.array(adj, dtype=float)
        eigenvalues = sorted(np.linalg.eigvalsh(A), reverse=True)
        props['lambda_max'] = eigenvalues[0]
        props['lambda_min'] = eigenvalues[-1]
        props['spectral_gap'] = eigenvalues[0] - eigenvalues[1] if m > 1 else 0
        props['eigenvalues'] = eigenvalues

        # Lovász theta (approximate via SDP relaxation using eigenvalues)
        # theta(G) = -n * lambda_min / (lambda_max - lambda_min)  [for regular graphs]
        # General bound: theta(G) >= -lambda_min * n / (lambda_max - lambda_min)
        if eigenvalues[0] > eigenvalues[-1]:
            theta_bound = 1 - eigenvalues[0] / eigenvalues[-1]
            props['lovasz_theta_bound'] = theta_bound

    return props


def check_neighborhood_structure(adj, v):
    """Analyze the structure of N(v) in the graph."""
    m = len(adj)
    neighbors = [j for j in range(m) if adj[v][j]]
    k = len(neighbors)
    if k <= 1:
        return {'type': 'trivial', 'size': k}

    # Count edges within N(v)
    internal_edges = sum(1 for i in range(k) for j in range(i+1, k)
                        if adj[neighbors[i]][neighbors[j]])

    # Check if N(v) is a union of at most 2 cliques
    # (quasi-line condition: each neighborhood is a union of at most 2 cliques)
    max_edges = k * (k - 1) // 2

    return {
        'size': k,
        'internal_edges': internal_edges,
        'max_edges': max_edges,
        'density': 2 * internal_edges / (k * (k - 1)) if k > 1 else 1,
        'is_clique': internal_edges == max_edges,
    }


# ============================================================
# Main analysis
# ============================================================
print("=" * 70)
print("SPECTRAL ANALYSIS OF OMEGA(T)")
print("=" * 70)

for n in [5, 6, 7, 8, 9]:
    print(f"\n{'='*60}")
    print(f"  n = {n}")
    print(f"{'='*60}")

    if n <= 6:
        mode = 'exhaustive'
        m = n * (n - 1) // 2
        count = 1 << m
    else:
        mode = 'random'
        count = 200 if n <= 8 else 100

    # Collect statistics
    stats = defaultdict(list)

    for idx in range(count):
        if mode == 'exhaustive':
            T = tournament_from_bits(n, idx)
        else:
            T = random_tournament(n)

        cycles = find_odd_cycles(T)
        if not cycles:
            continue

        cg = conflict_graph(cycles)
        props = graph_properties(cg)

        for key, val in props.items():
            if isinstance(val, (int, float)):
                stats[key].append(val)

        # Independence polynomial
        ip = indep_poly_coeffs(cg)
        stats['ip_degree'].append(len(ip) - 1)

        # Neighborhood analysis (sample a few vertices)
        m_omega = len(cg)
        sample_verts = range(min(5, m_omega))
        clique_nbhd_count = 0
        high_density_count = 0
        for v in sample_verts:
            ns = check_neighborhood_structure(cg, v)
            if ns.get('is_clique'):
                clique_nbhd_count += 1
            if ns.get('density', 0) > 0.8:
                high_density_count += 1

    # Report
    if not stats.get('vertices'):
        print("  No tournaments with odd cycles at this n")
        continue

    print(f"  Tournaments analyzed: {count}")
    print(f"  |V(Omega)| range: {min(stats['vertices'])}-{max(stats['vertices'])}, "
          f"avg={sum(stats['vertices'])/len(stats['vertices']):.1f}")
    print(f"  Edge density range: {min(stats['density']):.3f}-{max(stats['density']):.3f}, "
          f"avg={sum(stats['density'])/len(stats['density']):.3f}")
    print(f"  alpha range: {min(stats.get('alpha', [0]))}-{max(stats.get('alpha', [0]))}")
    print(f"  I.P. degree range: {min(stats['ip_degree'])}-{max(stats['ip_degree'])}")

    if HAS_NUMPY and 'lambda_max' in stats:
        print(f"  lambda_max range: {min(stats['lambda_max']):.3f}-{max(stats['lambda_max']):.3f}")
        print(f"  lambda_min range: {min(stats['lambda_min']):.3f}-{max(stats['lambda_min']):.3f}")
        print(f"  Spectral gap range: {min(stats['spectral_gap']):.3f}-{max(stats['spectral_gap']):.3f}")

        if 'lovasz_theta_bound' in stats:
            theta_vals = [v for v in stats['lovasz_theta_bound'] if v is not None]
            if theta_vals:
                print(f"  Lovász theta bound: {min(theta_vals):.3f}-{max(theta_vals):.3f}")

        # Key ratio: lambda_max / avg_degree (regularity indicator)
        if stats['avg_degree']:
            ratios = [l / d if d > 0 else 0
                     for l, d in zip(stats['lambda_max'], stats['avg_degree'])]
            print(f"  lambda_max / avg_degree: {min(ratios):.3f}-{max(ratios):.3f}, "
                  f"avg={sum(ratios)/len(ratios):.3f}")
            # For regular graphs this ratio = 1. Closer to 1 = more regular.

    print(f"  Min degree range: {min(stats['min_degree'])}-{max(stats['min_degree'])}")
    print(f"  Max degree range: {min(stats['max_degree'])}-{max(stats['max_degree'])}")

# ============================================================
# Deep dive: eigenvalue distribution at n=7,8
# ============================================================
if HAS_NUMPY:
    print(f"\n{'='*60}")
    print(f"  EIGENVALUE DISTRIBUTION DEEP DIVE")
    print(f"{'='*60}")

    for n in [7, 8]:
        print(f"\n  --- n = {n} ---")
        all_eigs = []
        for _ in range(500):
            T = random_tournament(n)
            cycles = find_odd_cycles(T)
            if not cycles:
                continue
            cg = conflict_graph(cycles)
            A = np.array(cg, dtype=float)
            eigs = np.linalg.eigvalsh(A)
            all_eigs.extend(eigs.tolist())

        all_eigs = np.array(all_eigs)
        print(f"  Total eigenvalues collected: {len(all_eigs)}")
        print(f"  Mean: {np.mean(all_eigs):.4f}")
        print(f"  Std:  {np.std(all_eigs):.4f}")
        print(f"  Min:  {np.min(all_eigs):.4f}")
        print(f"  Max:  {np.max(all_eigs):.4f}")

        # Histogram bins
        bins = np.linspace(np.min(all_eigs), np.max(all_eigs), 20)
        hist, _ = np.histogram(all_eigs, bins=bins)
        print(f"  Distribution shape (20 bins):")
        max_bar = max(hist)
        for i in range(len(hist)):
            bar = '#' * int(40 * hist[i] / max_bar) if max_bar > 0 else ''
            print(f"    [{bins[i]:6.2f}, {bins[i+1]:6.2f}): {bar} ({hist[i]})")

        # Check: fraction of eigenvalues that are negative
        neg_frac = np.sum(all_eigs < -0.01) / len(all_eigs)
        print(f"  Fraction negative: {neg_frac:.3f}")

        # Semicircle law comparison
        print(f"  (Compare: Wigner semicircle has mean=0, extends to ±2*sqrt(n))")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
