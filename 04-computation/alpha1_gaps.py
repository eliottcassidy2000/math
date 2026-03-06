#!/usr/bin/env python3
"""
Which alpha_1 values (total directed odd cycles) are achievable?
The gap at alpha_1=3 explains why H=7 is impossible.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations

def count_cycles_by_length(T, n):
    """Count directed odd cycles by length."""
    by_length = {}
    for k in range(3, n+1, 2):
        total = 0
        for combo in combinations(range(n), k):
            verts = list(combo)
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    total += c
            by_length[k] = total
    return by_length

# Collect achievable alpha_1 = total directed odd cycles
for n in range(3, 7):
    m = n*(n-1)//2
    achievable = set()
    # Also collect achievable (c3, c5, c7) tuples
    achievable_tuples = set()
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cycles = count_cycles_by_length(T, n)
        alpha_1 = sum(cycles.values())
        achievable.add(alpha_1)
        tup = tuple(cycles.get(k, 0) for k in range(3, n+1, 2))
        achievable_tuples.add(tup)

    all_vals = set(range(max(achievable) + 1))
    gaps = sorted(all_vals - achievable)
    print(f"\nn={n}: achievable alpha_1 = {sorted(achievable)}")
    print(f"  Gaps in alpha_1: {gaps}")
    print(f"  Number of distinct (c3,c5,...) tuples: {len(achievable_tuples)}")

    # Show cycle length breakdown
    print(f"  Achievable c3 values: {sorted(set(t[0] for t in achievable_tuples))}")
    if n >= 5:
        print(f"  Achievable c5 values: {sorted(set(t[1] for t in achievable_tuples))}")
    if n >= 7:
        print(f"  Achievable c7 values: {sorted(set(t[2] for t in achievable_tuples))}")

    # Can c3=3 occur?
    c3_vals = set(t[0] for t in achievable_tuples)
    print(f"  c3=3 achievable? {3 in c3_vals}")

# KEY QUESTION: Is alpha_1=3 impossible because c3=3 is impossible?
# At n=4: can have 0, 1, 2, or 4 cyclic triples. c3=3 seems impossible.
# Why? Because the total number of cyclic triples has a parity constraint!
# Moon's formula: c3 = C(n,3) - (1/2)*sum(s_i choose 2) where s_i = score of vertex i
# Since C(n,3) is fixed and sum(s_i choose 2) has a specific parity...

print("\n\n=== MOON's FORMULA ANALYSIS ===")
for n in range(3, 8):
    m = n*(n-1)//2
    c3_vals = set()
    for bits in range(min(1 << m, 1 << 15)):  # Cap at n=6 exhaustive
        T = tournament_from_bits(n, bits)
        scores = sorted(sum(T[i]) for i in range(n))
        # Moon: c3 = C(n,3) - sum(s_i * (s_i - 1) / 2) / 1
        # Actually: 3*c3 = C(n,3) - sum C(s_i, 2) + ... no
        # Moon's formula: c3 = C(n,3) - sum_{i} C(s_i, 2)
        # No: the correct formula is:
        # In any tournament, number of 3-cycles = C(n,3) - sum_{i} C(s_i^+, 2)
        # where s_i^+ = out-degree of vertex i
        s_plus = [sum(T[i]) for i in range(n)]
        c3_moon = n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in s_plus)
        c3_vals.add(c3_moon)
    c3_sorted = sorted(c3_vals)
    gaps = [x for x in range(max(c3_sorted)+1) if x not in c3_vals]
    print(f"n={n}: achievable c3 = {c3_sorted}")
    if gaps:
        print(f"  c3 gaps: {gaps}")

print("\nDone.")
