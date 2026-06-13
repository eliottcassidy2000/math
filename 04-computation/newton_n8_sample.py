#!/usr/bin/env python3
"""
Sample n=8 tournaments and check Newton's inequality for full Omega(T).
The question: does the Newton failure at n=9 also appear at n=8?

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament
from itertools import combinations, permutations
import numpy as np

def find_all_directed_odd_cycles(T, n, max_len=None):
    if max_len is None:
        max_len = n
    cycles = set()
    for length in range(3, max_len + 1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                if all(T[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
    return list(cycles)

def independence_poly_dc(adj_dict):
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return [1]
        v = max(verts, key=lambda u: len(adj_dict[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj_dict[v] & verts) - {v})
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]
        memo[verts] = result
        return result
    return solve(frozenset(adj_dict.keys()))

n = 8
print(f"Sampling n={n} tournaments for Newton's inequality on full Omega...")

newton_fails = 0
real_root_fails = 0
total_checked = 0
samples = 200

for trial in range(samples):
    T = random_tournament(n)
    cycles = find_all_directed_odd_cycles(T, n)
    m = len(cycles)
    if m < 2:
        continue

    cycle_vsets = [frozenset(c) for c in cycles]

    # Build adjacency
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and cycle_vsets[i] & cycle_vsets[j])

    # Quick Newton check: a1^2 >= 2*a2
    a1 = m
    a2 = sum(1 for i in range(m) for j in range(i+1, m) if not (cycle_vsets[i] & cycle_vsets[j]))

    total_checked += 1
    if a1**2 < 2 * a2:
        newton_fails += 1
        print(f"  Trial {trial}: NEWTON FAIL a1={a1}, a2={a2}, a1^2={a1**2}, 2*a2={2*a2}")

    # For a subset, check full real-rootedness
    if trial < 50 and m <= 40:
        coeffs = independence_poly_dc(adj)
        deg = len(coeffs) - 1
        while deg > 0 and coeffs[deg] == 0:
            deg -= 1
        if deg >= 2:
            poly_c = [coeffs[deg-i] for i in range(deg+1)]
            roots = np.roots(poly_c)
            if not all(abs(r.imag) < 1e-8 for r in roots):
                real_root_fails += 1
                print(f"  Trial {trial}: REAL-ROOT FAIL, coeffs={coeffs[:deg+1]}")

    if (trial + 1) % 50 == 0:
        print(f"  ... {trial+1}/{samples} done")

print(f"\nResults for n={n}:")
print(f"  Checked: {total_checked}")
print(f"  Newton failures: {newton_fails}")
print(f"  Real-root failures: {real_root_fails} (checked {min(50, total_checked)} tournaments)")
