#!/usr/bin/env python3
"""
Analyze how c5 and c7 individually depend on disj_33 at n=7 regular.

The cascade: 3-cycle overlap → determines c5 and c7 individually.
c5 + c7 = (disj^2 - 27*disj + 272)/2 (quadratic in disj).

Do c5 and c7 also have quadratic formulas?

opus-2026-03-13-S71b
"""

import itertools
from collections import Counter

def all_tournaments_n7():
    n = 7
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj, bits

def is_regular(adj, n):
    return all(sum(adj[i]) == (n-1)//2 for i in range(n))

def count_directed_cycles(adj, n, length):
    """Count directed cycles of given length."""
    count = 0
    for combo in itertools.combinations(range(n), length):
        for perm in itertools.permutations(combo):
            valid = True
            for i in range(length):
                if not adj[perm[i]][perm[(i+1) % length]]:
                    valid = False
                    break
            if valid:
                count += 1
    return count // length  # each cycle counted 'length' times

def count_disjoint_3pairs(adj, n):
    """Count disjoint 3-cycle pairs."""
    # Find all 3-cycles
    three_cycles = []
    for combo in itertools.combinations(range(n), 3):
        for perm in itertools.permutations(combo):
            valid = True
            for i in range(3):
                if not adj[perm[i]][perm[(i+1) % 3]]:
                    valid = False
                    break
            if valid:
                # Canonical form: start from min vertex
                min_idx = perm.index(min(perm))
                canonical = perm[min_idx:] + perm[:min_idx]
                three_cycles.append(set(canonical))
                break  # only one canonical direction per triple

    # Actually, count unique 3-cycles properly
    seen = set()
    unique_3 = []
    for combo in itertools.combinations(range(n), 3):
        a, b, c = combo
        # Check both directions
        if adj[a][b] and adj[b][c] and adj[c][a]:
            unique_3.append(set(combo))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            unique_3.append(set(combo))

    # Count disjoint pairs
    disj = 0
    for i in range(len(unique_3)):
        for j in range(i+1, len(unique_3)):
            if not (unique_3[i] & unique_3[j]):
                disj += 1
    return disj, len(unique_3)

n = 7
print(f"Cycle cascade analysis at n={n}")
print(f"{'disj':>6s} {'c3':>5s} {'c5':>5s} {'c7':>5s} {'c5+c7':>7s} {'count':>7s}")

data = {}
count = 0
for adj, bits in all_tournaments_n7():
    if not is_regular(adj, n):
        continue
    count += 1

    c3 = count_directed_cycles(adj, n, 3)
    c5 = count_directed_cycles(adj, n, 5)
    c7 = count_directed_cycles(adj, n, 7)
    disj, n3 = count_disjoint_3pairs(adj, n)

    key = (disj, c3, c5, c7)
    if key not in data:
        data[key] = 0
        print(f"{disj:6d} {c3:5d} {c5:5d} {c7:5d} {c5+c7:7d}", flush=True)
    data[key] += 1

print(f"\nTotal regular tournaments: {count}")
print(f"\nDistinct (disj, c3, c5, c7) classes: {len(data)}")

# Check if c5 and c7 individually depend only on disj
by_disj = {}
for (disj, c3, c5, c7), cnt in data.items():
    if disj not in by_disj:
        by_disj[disj] = []
    by_disj[disj].append((c5, c7, cnt))

print(f"\nPer-disj values:")
for disj in sorted(by_disj.keys()):
    vals = by_disj[disj]
    c5_vals = set(v[0] for v in vals)
    c7_vals = set(v[1] for v in vals)
    total = sum(v[2] for v in vals)
    print(f"  disj={disj}: c5 ∈ {sorted(c5_vals)}, c7 ∈ {sorted(c7_vals)}, count={total}")
    for c5, c7, cnt in sorted(vals):
        c5_c7_pred = (disj**2 - 27*disj + 272) // 2
        print(f"    c5={c5}, c7={c7}, c5+c7={c5+c7} (pred={c5_c7_pred}), n={cnt}")

# Check: is c5 a linear function of disj? quadratic?
print("\nQuadratic fit for c5 and c7:")
import numpy as np
disj_vals = sorted(by_disj.keys())
c5_vals = [by_disj[d][0][0] for d in disj_vals]
c7_vals = [by_disj[d][0][1] for d in disj_vals]

# c5 = a*disj^2 + b*disj + c
A = np.array([[d**2, d, 1] for d in disj_vals])
b5 = np.array(c5_vals)
b7 = np.array(c7_vals)

coeffs_5 = np.linalg.solve(A, b5)
coeffs_7 = np.linalg.solve(A, b7)
print(f"  c5 = {coeffs_5[0]:.1f}*disj^2 + {coeffs_5[1]:.1f}*disj + {coeffs_5[2]:.1f}")
print(f"  c7 = {coeffs_7[0]:.1f}*disj^2 + {coeffs_7[1]:.1f}*disj + {coeffs_7[2]:.1f}")

# Verify
for d in disj_vals:
    c5_pred = coeffs_5[0]*d**2 + coeffs_5[1]*d + coeffs_5[2]
    c7_pred = coeffs_7[0]*d**2 + coeffs_7[1]*d + coeffs_7[2]
    actual = by_disj[d][0]
    print(f"  disj={d}: c5_pred={c5_pred:.1f} (actual={actual[0]}), c7_pred={c7_pred:.1f} (actual={actual[1]})")

# Connection to H
print("\n" + "="*60)
print("FULL DECOMPOSITION")
print("="*60)
for d in disj_vals:
    c5, c7, cnt = by_disj[d][0]
    H = d**2 - 23*d + 301
    alpha_1 = 14 + c5 + c7
    alpha_2 = d
    print(f"disj={d}: H={H}, α₁={alpha_1}, α₂={alpha_2}")
    print(f"  Check: 1 + 2*{alpha_1} + 4*{alpha_2} = {1+2*alpha_1+4*alpha_2}")
    print(f"  c5={c5}, c7={c7}")
    print(f"  3-cycle overlap edges = C(14,2) - disj = {14*13//2 - d}")
