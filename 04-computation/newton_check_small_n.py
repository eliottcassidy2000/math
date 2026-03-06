#!/usr/bin/env python3
"""
Check Newton's first inequality for I(Omega(T), x) at small n.
Newton requires a1^2 >= 2*a0*a2 (i.e., a1^2 >= 2*a2 since a0=1).

Uses directed odd cycles as vertices of Omega.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits
from itertools import combinations, permutations

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

for n in range(3, 9):
    total = 2**(n*(n-1)//2)
    fail_count = 0
    checked = 0
    worst = None
    for bits in range(total):
        T = tournament_from_bits(n, bits)
        cycles = find_all_directed_odd_cycles(T, n)
        m = len(cycles)
        if m < 2:
            continue
        cycle_vsets = [frozenset(c) for c in cycles]
        a1 = m
        a2 = 0
        for i in range(m):
            for j in range(i+1, m):
                if not (cycle_vsets[i] & cycle_vsets[j]):
                    a2 += 1
        # Newton: a1^2 >= 2*a2
        if a1**2 < 2 * a2:
            fail_count += 1
            if worst is None or (2*a2 - a1**2) > worst[0]:
                worst = (2*a2 - a1**2, a1, a2)
        checked += 1
    msg = f"n={n}: {checked}/{total} checked, Newton failures: {fail_count}"
    if worst:
        msg += f" (worst: a1={worst[1]}, a2={worst[2]}, gap={worst[0]})"
    print(msg)
