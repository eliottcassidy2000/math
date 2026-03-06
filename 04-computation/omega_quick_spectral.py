#!/usr/bin/env python3
"""
Quick spectral analysis of Omega(T) via random sampling.
Focused on getting key structural data at n=5,6,7,8,9.

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


for n in [5, 6, 7, 8, 9]:
    print(f"\n{'='*60}")
    print(f"  n = {n} (random 100)")
    print(f"{'='*60}")

    results = []
    for _ in range(100):
        T = random_tournament(n)
        cycles = find_odd_cycles(T)
        if not cycles:
            continue

        # Use only 3-cycles for n>=8 (full is too slow)
        if n >= 8:
            cycles = [c for c in cycles if len(c) == 3]
            if not cycles:
                continue

        cg = conflict_graph(cycles)
        m = len(cg)
        ip = indep_poly_coeffs(cg)
        d = len(ip) - 1  # degree of I.P.

        A = np.array(cg, dtype=float)
        eigs = sorted(np.linalg.eigvalsh(A), reverse=True)

        deg = [sum(cg[i]) for i in range(m)]
        avg_deg = sum(deg) / m if m > 0 else 0
        density = sum(sum(row) for row in cg) / (m * (m - 1)) if m > 1 else 0

        # Check real roots of I.P.
        if d >= 1:
            roots = np.roots(list(reversed(ip)))
            all_real = all(abs(r.imag) < 1e-6 for r in roots)
            all_neg = all(r.real < -1e-8 for r in roots) if all_real else False
        else:
            all_real = True
            all_neg = True

        # Check if eigenvalues of complement-adjacency give PSD matrix
        # M = max_eig*I - A has eigenvalues max_eig - lambda_i >= 0
        # This is always PSD. Its eigenvalues are max_eig - eigs[i].
        max_eig = eigs[0]
        complement_eigs = [max_eig - e for e in eigs]

        results.append({
            'm': m, 'd': d, 'ip': ip, 'density': density,
            'avg_deg': avg_deg, 'max_eig': max_eig, 'min_eig': eigs[-1],
            'spectral_gap': eigs[0] - eigs[1] if m > 1 else 0,
            'all_real': all_real, 'all_neg': all_neg,
            'ratio': max_eig / avg_deg if avg_deg > 0 else 0,
        })

    real_root_fails = sum(1 for r in results if not r['all_real'])
    neg_root_fails = sum(1 for r in results if r['all_real'] and not r['all_neg'])

    print(f"  Analyzed: {len(results)}")
    print(f"  |V(Omega)| range: {min(r['m'] for r in results)}-{max(r['m'] for r in results)}")
    print(f"  I.P. degree range: {min(r['d'] for r in results)}-{max(r['d'] for r in results)}")
    print(f"  Density: avg={np.mean([r['density'] for r in results]):.3f}")
    print(f"  lambda_max: avg={np.mean([r['max_eig'] for r in results]):.2f}, "
          f"range=[{min(r['max_eig'] for r in results):.2f}, {max(r['max_eig'] for r in results):.2f}]")
    print(f"  lambda_min: avg={np.mean([r['min_eig'] for r in results]):.2f}, "
          f"range=[{min(r['min_eig'] for r in results):.2f}, {max(r['min_eig'] for r in results):.2f}]")
    print(f"  lambda_max / avg_degree: avg={np.mean([r['ratio'] for r in results]):.3f}")
    print(f"  Real-root failures: {real_root_fails}/{len(results)}")
    print(f"  Negative-root failures: {neg_root_fails}/{len(results)}")

    # Show unique I.P. types
    ip_types = {}
    for r in results:
        key = tuple(r['ip'])
        ip_types[key] = ip_types.get(key, 0) + 1
    print(f"  Unique I.P. types: {len(ip_types)}")
    for ip, cnt in sorted(ip_types.items(), key=lambda x: -x[1])[:8]:
        print(f"    {list(ip)}: {cnt} times")

print("\nDONE")
