#!/usr/bin/env python3
"""
efwd4_invariants.py — What determines E[fwd^4]?

E[fwd^4] is NOT determined by t3 alone (proved at n=5,6).
What additional invariant is needed?

Candidates:
- t5: 5-cycle count
- alpha_2: disjoint 3-cycle pairs
- #fwd5path: directed 5-path count
- sum of squares of out-degrees (s2)
- specific 4-vertex subgraph counts

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict
import numpy as np

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def count_3cycles(adj, n):
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

def count_kcycles(adj, n, k):
    count = 0
    for combo in combinations(range(n), k):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[(i+1)%k]] for i in range(k)):
                count += 1
    return count // k

def count_fwd_paths(adj, n, length):
    """Count forward (length)-step paths = directed paths on (length+1) vertices."""
    count = 0
    for combo in combinations(range(n), length + 1):
        for perm in permutations(combo):
            if all(adj[perm[i]][perm[i+1]] for i in range(length)):
                count += 1
    return count

print("=" * 60)
print("WHAT DETERMINES E[fwd^4]?")
print("=" * 60)

for n in [5, 6]:
    m_vals = n*(n-1)//2
    seen = set()
    data = []

    for bits in range(1 << m_vals):
        adj = tournament_from_bits(n, bits)
        F = compute_F(adj, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        t3 = count_3cycles(adj, n)
        t5 = count_kcycles(adj, n, 5) if n >= 5 else 0

        # Disjoint 3-cycle pairs
        cycles3 = []
        for i, j, k in combinations(range(n), 3):
            if (adj[i][j] and adj[j][k] and adj[k][i]) or \
               (adj[i][k] and adj[k][j] and adj[j][i]):
                cycles3.append((i, j, k))
        alpha_2 = sum(1 for a in range(len(cycles3))
                      for b in range(a+1, len(cycles3))
                      if set(cycles3[a]).isdisjoint(set(cycles3[b])))

        # Count forward 4-paths and 5-paths
        fwd4 = count_fwd_paths(adj, n, 3)
        fwd5 = count_fwd_paths(adj, n, 4)

        # Count specific 4-vertex subgraph types
        # For 4 vertices, classify tournaments by # of 3-cycles
        # Each 4-tournament has 0, 1, or 4 three-cycles
        # W4 = # of 4-vertex subsets with exactly 1 three-cycle
        # R4 = # of 4-vertex subsets with 4 three-cycles (regular tournament = C_4)
        W4 = 0
        R4 = 0
        for combo in combinations(range(n), 4):
            sub_t3 = 0
            for triple in combinations(combo, 3):
                i, j, k = triple
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    sub_t3 += 1
            if sub_t3 == 1:
                W4 += 1
            elif sub_t3 == 4:
                R4 += 1

        total = sum(F)
        m4 = Fraction(sum(k**4 * F[k] for k in range(n)), total)

        data.append({
            't3': t3, 't5': t5, 'alpha_2': alpha_2,
            'fwd4': fwd4, 'fwd5': fwd5,
            'W4': W4, 'R4': R4,
            'E_fwd4': float(m4)
        })

    print(f"\nn={n}: {len(data)} F-classes")

    # Try various combinations
    y = np.array([d['E_fwd4'] for d in data])

    combos = [
        ('t3', lambda d: [1, d['t3']]),
        ('t3,t5', lambda d: [1, d['t3'], d['t5']]),
        ('t3,alpha_2', lambda d: [1, d['t3'], d['alpha_2']]),
        ('t3,R4', lambda d: [1, d['t3'], d['R4']]),
        ('t3,W4', lambda d: [1, d['t3'], d['W4']]),
        ('t3,fwd5', lambda d: [1, d['t3'], d['fwd5']]),
        ('t3,t3^2', lambda d: [1, d['t3'], d['t3']**2]),
        ('t3,t5,alpha_2', lambda d: [1, d['t3'], d['t5'], d['alpha_2']]),
        ('t3,t5,R4', lambda d: [1, d['t3'], d['t5'], d['R4']]),
        ('t3,R4,alpha_2', lambda d: [1, d['t3'], d['R4'], d['alpha_2']]),
    ]

    for name, feat in combos:
        X = np.array([feat(d) for d in data], dtype=float)
        c, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        err = max(abs(X @ c - y))
        exact = "EXACT" if err < 1e-8 else f"err={err:.6f}"
        print(f"  {name:>20}: {exact}  coeffs=[{', '.join(f'{x:.4f}' for x in c)}]")

    # Print table for n=5 to see what varies
    if n == 5:
        print(f"\n  Detailed table:")
        print(f"  {'t3':>3} {'t5':>3} {'a2':>3} {'R4':>3} {'W4':>3} {'fwd5':>6} {'E[fwd^4]':>12}")
        for d in sorted(data, key=lambda x: (x['t3'], x['t5'])):
            print(f"  {d['t3']:>3} {d['t5']:>3} {d['alpha_2']:>3} {d['R4']:>3} {d['W4']:>3} {d['fwd5']:>6} {d['E_fwd4']:>12.4f}")
