#!/usr/bin/env python3
"""
worpitzky_n7_fast.py — Fast n=7 Worpitzky test.

Reduced sample + skip expensive invariants until needed.

Author: opus-2026-03-07-S46c
"""
from itertools import permutations, combinations
from math import comb, factorial
from fractions import Fraction
from collections import defaultdict
import random
import sys

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

def compute_F_fast(adj, n):
    """Compute F-vector using cached permutations."""
    F = [0]*n
    for P in ALL_PERMS:
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def worpitzky_a(F, n, m):
    return sum(F[k] * comb(m + n - 1 - k, n - 1) for k in range(n))

def exact_worpitzky_coeffs(F, n):
    N = n
    pts = list(range(N))
    vals = [Fraction(worpitzky_a(F, n, m)) for m in pts]
    mat = [[Fraction(m**j) for j in range(N)] + [vals[i]] for i, m in enumerate(pts)]
    for col in range(N):
        for row in range(col, N):
            if mat[row][col] != 0:
                mat[col], mat[row] = mat[row], mat[col]
                break
        pivot = mat[col][col]
        for row in range(col+1, N):
            factor = mat[row][col] / pivot
            for k in range(N+1):
                mat[row][k] -= factor * mat[col][k]
    coeffs = [Fraction(0)] * N
    for row in range(N-1, -1, -1):
        val = mat[row][N]
        for k in range(row+1, N):
            val -= mat[row][k] * coeffs[k]
        coeffs[row] = val / mat[row][row]
    return coeffs

def count_3cycles(adj, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or \
           (adj[i][k] and adj[k][j] and adj[j][i]):
            t3 += 1
    return t3

n = 7
ALL_PERMS = list(permutations(range(n)))
m_edges = n*(n-1)//2

random.seed(42)
seen_F = {}
all_data = []

print(f"n={n}: Sampling tournaments...")

# Sample
for trial in range(30000):
    bits = random.getrandbits(m_edges)
    adj = tournament_from_bits(n, bits)
    F = compute_F_fast(adj, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F[key] = True

    t3 = count_3cycles(adj, n)

    coeffs = exact_worpitzky_coeffs(F, n)
    delta = [coeffs[j] - Fraction(comb(n, j)) for j in range(n)]

    all_data.append({
        't3': t3, 'coeffs': coeffs, 'delta': delta,
        'F': F, 'H': F[n-1], 'adj': adj
    })

    if len(all_data) % 10 == 0:
        print(f"  {len(all_data)} distinct F-vectors (trial {trial})", file=sys.stderr)

# Also add transitive
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
F_trans = compute_F_fast(adj_trans, n)
key_trans = tuple(F_trans)
if key_trans not in seen_F:
    seen_F[key_trans] = True
    t3 = count_3cycles(adj_trans, n)
    coeffs = exact_worpitzky_coeffs(F_trans, n)
    delta = [coeffs[j] - Fraction(comb(n, j)) for j in range(n)]
    all_data.append({'t3': t3, 'coeffs': coeffs, 'delta': delta, 'F': F_trans, 'H': F_trans[n-1], 'adj': adj_trans})

print(f"\nTotal distinct F-vectors: {len(all_data)}")

# Check universal levels
print("\n" + "=" * 60)
print("CHECKING WHICH LEVELS ARE UNIVERSAL")
print("=" * 60)

for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) == 1 and Fraction(0) in vals:
        print(f"  delta_{j} (coeff of m^{j}): UNIVERSAL = 0")
    else:
        print(f"  delta_{j} (coeff of m^{j}): {len(vals)} distinct values")

# Check t3 dependence for each level
print("\n" + "=" * 60)
print("CHECKING t3 DEPENDENCE")
print("=" * 60)

for j in range(n-1, -1, -1):
    vals = set(d['delta'][j] for d in all_data)
    if len(vals) <= 1:
        continue

    t3_to_delta = defaultdict(set)
    for d in all_data:
        t3_to_delta[d['t3']].add(d['delta'][j])
    determined = all(len(v) == 1 for v in t3_to_delta.values())

    if determined:
        items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_delta.items())]
        if len(items) >= 2:
            slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0])
            intercept = items[0][1] - slope * items[0][0]
            ok = all(v == intercept + slope * t for t, v in items)
            if ok:
                print(f"  delta_{j}: {intercept} + {slope}*t3  [EXACT LINEAR]")
            else:
                print(f"  delta_{j}: determined by t3 but NOT linear")
    else:
        n_amb = sum(1 for v in t3_to_delta.values() if len(v) > 1)
        print(f"  delta_{j}: NOT determined by t3 alone ({n_amb} ambiguous)")

# Verify c0 = H(T)
print("\n" + "=" * 60)
print("VERIFYING c0 = H(T)")
print("=" * 60)
all_match = all(d['coeffs'][0] == Fraction(d['H']) for d in all_data)
print(f"  c0 = H(T): {'ALL MATCH' if all_match else 'MISMATCH'} ({len(all_data)} classes)")

# Check E[fwd^3] = A + 6/n * t3
print("\n" + "=" * 60)
print("E[fwd^3] = A(n) + 6/n * t3?")
print("=" * 60)

t3_to_m3 = defaultdict(set)
for d in all_data:
    total = factorial(n)
    m3 = Fraction(sum(k**3 * d['F'][k] for k in range(n)), total)
    t3_to_m3[d['t3']].add(m3)

determined = all(len(v) == 1 for v in t3_to_m3.values())
if determined:
    items = [(t3, list(vals)[0]) for t3, vals in sorted(t3_to_m3.items())]
    slope = (items[1][1] - items[0][1]) / (items[1][0] - items[0][0])
    A_val = items[0][1] - slope * items[0][0]
    pred_slope = Fraction(6, n)
    print(f"  E[fwd^3] = {A_val} + {slope}*t3")
    print(f"  Predicted slope: 6/{n} = {pred_slope}")
    print(f"  Match: {slope == pred_slope}")
else:
    print(f"  E[fwd^3] NOT determined by t3 alone!")

# Print compact table
print("\n" + "=" * 60)
print(f"COMPACT TABLE ({len(all_data)} classes)")
print("=" * 60)
print(f"{'t3':>3} {'H':>5} {'d5':>4} {'d4':>6} {'d3':>6} {'d2':>8} {'d1':>10} {'d0':>6}")
for d in sorted(all_data, key=lambda x: (x['t3'], x['H'])):
    dd = d['delta']
    print(f"{d['t3']:>3} {d['H']:>5} {int(dd[5]):>4} {int(dd[4]):>6} {int(dd[3]):>6} {int(dd[2]):>8} {int(dd[1]):>10} {int(dd[0]):>6}")
