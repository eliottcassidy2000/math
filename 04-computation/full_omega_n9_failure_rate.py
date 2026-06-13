#!/usr/bin/env python3
"""
Measure failure rate of real-rootedness for FULL Omega(T) at n=9.
Focus on tournaments with skewed score sequences (where Omega_3 fails).

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

n = 9

# Strategy: generate tournaments with HIGH score variance (skewed)
# These are more likely to fail.
# Score variance = sum (s_i - 4)^2 / 9 for n=9 (mean out-degree = 4)

print("Testing full Omega real-rootedness at n=9...")
print("Targeting high-variance score sequences.\n")

total = 0
fail_full = 0
fail_omega3 = 0
tested_full = 0

# Also test the known counterexample variants
# The counterexample has scores (1,1,3,4,4,4,6,6,7), var = 4.22

for trial in range(10000):
    T = random_tournament(n)
    scores = sorted(sum(row) for row in T)
    var = sum((s - 4)**2 for s in scores) / 9

    # Focus on high-variance tournaments
    if var < 3.0:
        continue

    total += 1

    # Check Omega_3 first (fast)
    cycles3 = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles3.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles3.append((i,j,k))

    a1 = len(cycles3)
    if a1 < 3:
        continue

    cycle_sets3 = [frozenset(c) for c in cycles3]
    a2 = sum(1 for i in range(a1) for j in range(i+1, a1)
             if not (cycle_sets3[i] & cycle_sets3[j]))
    a3 = 0
    for i in range(a1):
        for j in range(i+1, a1):
            if cycle_sets3[i] & cycle_sets3[j]:
                continue
            uij = cycle_sets3[i] | cycle_sets3[j]
            for k in range(j+1, a1):
                if not (cycle_sets3[k] & uij):
                    a3 += 1

    if a3 > 0:
        disc = 18*a1*a2*a3 - 4*a1**3*a3 + a1**2*a2**2 - 4*a2**3 - 27*a3**2
        if disc < 0:
            fail_omega3 += 1
            print(f"  Trial {trial}: Omega_3 FAIL! scores={tuple(scores)}, var={var:.2f}, (a1,a2,a3)=({a1},{a2},{a3}), disc={disc}")

    # Check full Omega (slower)
    all_cycles = find_all_directed_odd_cycles(T, n)
    cycle_list = list(all_cycles)
    m = len(cycle_list)
    cycle_vsets = [frozenset(c) for c in cycle_list]

    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and cycle_vsets[i] & cycle_vsets[j])

    coeffs = independence_poly_dc(adj)
    deg = len(coeffs) - 1
    while deg > 0 and coeffs[deg] == 0:
        deg -= 1

    tested_full += 1

    if deg >= 2:
        # Check Newton's inequalities
        newton_fail = False
        for k in range(1, deg):
            if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1] * (k+1) / k:
                newton_fail = True
                break

        if newton_fail:
            fail_full += 1
            val2 = sum(coeffs[i] * 2**i for i in range(deg+1))
            print(f"  Trial {trial}: FULL Omega FAIL! scores={tuple(scores)}, var={var:.2f}, {m} cycles, coeffs={coeffs[:deg+1]}, I(2)={val2}")

    if tested_full % 100 == 0:
        print(f"  ... {tested_full} full-Omega tests done ({fail_full} failures)")

print(f"\n{'='*70}")
print(f"Results:")
print(f"  High-variance tournaments tested: {total}")
print(f"  Full Omega tests: {tested_full}")
print(f"  Omega_3 failures: {fail_omega3}")
print(f"  Full Omega failures: {fail_full}")
if tested_full > 0:
    print(f"  Full Omega failure rate: {fail_full/tested_full*100:.2f}%")
