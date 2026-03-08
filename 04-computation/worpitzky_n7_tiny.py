#!/usr/bin/env python3
"""
worpitzky_n7_tiny.py — Minimal n=7 Worpitzky test.

Use numpy instead of Fraction for speed, small sample.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from collections import defaultdict
import random
import numpy as np
import sys

n = 7
# Pre-compute permutations list
ALL_PERMS = list(permutations(range(n)))

def tournament_from_bits(bits):
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

def compute_F(adj):
    F = [0]*n
    for P in ALL_PERMS:
        fwd = 0
        for i in range(n-1):
            fwd += adj[P[i]][P[i+1]]
        F[fwd] += 1
    return F

def worpitzky_a(F, m):
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

def count_3cycles(adj):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

m_edges = n*(n-1)//2
random.seed(42)
seen_F = {}
all_data = []

print(f"n={n}: Sampling tournaments...", flush=True)

# Small sample: 5000
for trial in range(5000):
    bits = random.getrandbits(m_edges)
    adj = tournament_from_bits(bits)
    F = compute_F(adj)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F[key] = True

    t3 = count_3cycles(adj)

    # Use numpy polyfit for speed
    m_pts = list(range(n + 2))
    a_vals = [worpitzky_a(F, m) for m in m_pts]
    coeffs = np.polyfit(m_pts, a_vals, n - 1)  # high to low

    # Convert to delta from binomial
    # coeffs[0] = coeff of m^6, coeffs[1] = m^5, ...
    # c_j (coeff of m^j) = coeffs[n-1-j]
    delta = [round(coeffs[n-1-j] - comb(n, j)) for j in range(n)]

    all_data.append({
        't3': t3, 'delta': delta, 'F': F, 'H': F[n-1]
    })

    if len(all_data) % 10 == 0:
        print(f"  {len(all_data)} distinct F-vectors (trial {trial})", file=sys.stderr, flush=True)

# Also add transitive
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
F_trans = compute_F(adj_trans)
key_trans = tuple(F_trans)
if key_trans not in seen_F:
    seen_F[key_trans] = True
    t3 = 0
    m_pts = list(range(n + 2))
    a_vals = [worpitzky_a(F_trans, m) for m in m_pts]
    coeffs = np.polyfit(m_pts, a_vals, n - 1)
    delta = [round(coeffs[n-1-j] - comb(n, j)) for j in range(n)]
    all_data.append({'t3': t3, 'delta': delta, 'F': F_trans, 'H': F_trans[n-1]})

print(f"\nTotal distinct F-vectors: {len(all_data)}", flush=True)

# Check levels
print("\nUNIVERSAL LEVELS:", flush=True)
for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) == 1 and 0 in vals:
        print(f"  delta_{j}: UNIVERSAL = 0")
    else:
        print(f"  delta_{j}: {len(vals)} distinct values, range [{min(vals)}, {max(vals)}]")

# Check t3 dependence
print("\nt3 DEPENDENCE:", flush=True)
for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) <= 1:
        continue
    t3_to_d = defaultdict(set)
    for d in all_data:
        t3_to_d[d['t3']].add(d['delta'][j])
    determined = all(len(v) == 1 for v in t3_to_d.values())
    if determined:
        items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_d.items())]
        slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0]) if len(items) >= 2 else 0
        print(f"  delta_{j}: determined by t3, slope ~ {slope:.4f}")
    else:
        n_amb = sum(1 for v in t3_to_d.values() if len(v) > 1)
        print(f"  delta_{j}: NOT determined by t3 ({n_amb} ambiguous)")

# Verify c0 = H
all_match = all(round(d['delta'][0]) == d['H'] - 1 for d in all_data)
print(f"\nc0 = H(T): {'ALL MATCH' if all_match else 'MISMATCH'}")

# E[fwd^3] = A(n) + 6/n * t3
print("\nE[fwd^3] check:", flush=True)
total = factorial(n)
t3_to_m3 = defaultdict(set)
for d in all_data:
    m3 = sum(k**3 * d['F'][k] for k in range(n)) / total
    t3_to_m3[d['t3']].add(round(m3, 10))
det_m3 = all(len(v) == 1 for v in t3_to_m3.values())
print(f"  E[fwd^3] determined by t3: {det_m3}")
if det_m3:
    items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_m3.items())]
    if len(items) >= 2:
        slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0])
        print(f"  slope = {slope:.8f}, 6/n = {6/n:.8f}, match = {abs(slope - 6/n) < 0.0001}")

# Compact table
print(f"\nCOMPACT TABLE ({len(all_data)} classes):", flush=True)
print(f"{'t3':>3} {'H':>5} {'d5':>4} {'d4':>6} {'d3':>6} {'d2':>8} {'d1':>10} {'d0':>6}")
for d in sorted(all_data, key=lambda x: (x['t3'], x['H'])):
    dd = d['delta']
    print(f"{d['t3']:>3} {d['H']:>5} {dd[5]:>4} {dd[4]:>6} {dd[3]:>6} {dd[2]:>8} {dd[1]:>10} {dd[0]:>6}")
