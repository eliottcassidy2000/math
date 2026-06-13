#!/usr/bin/env python3
"""
Verify deformed Eulerian number formula at n=9.

At n=9, the invariants are: t3 (f=6), t5 (f=4), t7 (f=2), t9 (f=0),
bc (f=4, parts=2), bc35 (f=2, parts=2), bc37 (f=0, parts=2), a3 (f=2, parts=3).

a_k(T) = A(9,k) + sum_I 2^{parts(I)} * c_k^{(f_I, 8)} * I(T)

opus-2026-03-07-S32
"""
from itertools import combinations
from collections import defaultdict
from math import comb
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def inflated_eulerian(f, d, k):
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    """Pairs of vertex-disjoint 3-cycles."""
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def count_bc35(A, n):
    """Pairs: one 3-cycle and one disjoint 5-cycle."""
    cyc3 = []
    for t in combinations(range(n), 3):
        sub = [[A[t[i]][t[j]] for j in range(3)] for i in range(3)]
        if sub[0][1]*sub[1][2]*sub[2][0] or sub[0][2]*sub[2][1]*sub[1][0]:
            cyc3.append(set(t))

    cyc5 = []
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for m in range(1, 1 << 5):
            for v in range(5):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(5):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 5) - 1
        ct = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        if ct > 0:
            cyc5.append((set(verts), ct))

    total = 0
    for c3 in cyc3:
        for c5_verts, c5_ct in cyc5:
            if c3.isdisjoint(c5_verts):
                total += c5_ct
    return total

def count_bc37(A, n):
    """Pairs: one 3-cycle and one disjoint 7-cycle."""
    cyc3 = []
    for t in combinations(range(n), 3):
        sub = [[A[t[i]][t[j]] for j in range(3)] for i in range(3)]
        if sub[0][1]*sub[1][2]*sub[2][0] or sub[0][2]*sub[2][1]*sub[1][0]:
            cyc3.append(set(t))

    cyc7 = []
    for verts in combinations(range(n), 7):
        sub = [[A[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for m in range(1, 1 << 7):
            for v in range(7):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(7):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 7) - 1
        ct = sum(dp[full][v] for v in range(1, 7) if sub[v][0])
        if ct > 0:
            cyc7.append((set(verts), ct))

    total = 0
    for c3 in cyc3:
        for c7_verts, c7_ct in cyc7:
            if c3.isdisjoint(c7_verts):
                total += c7_ct
    return total

def count_alpha3(A, n):
    """Triples of vertex-disjoint 3-cycles."""
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    total = 0
    for i in range(len(cyc3)):
        for j in range(i+1, len(cyc3)):
            if not cyc3[i].isdisjoint(cyc3[j]): continue
            for k in range(j+1, len(cyc3)):
                if cyc3[k].isdisjoint(cyc3[i]) and cyc3[k].isdisjoint(cyc3[j]):
                    total += 1
    return total

def forward_edge_dist_dp(A, n):
    """Count permutations by forward edges using Hamiltonian path DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

# ====================================================================
# Setup
# ====================================================================
n = 9
d = n - 1  # = 8

# Invariants at n=9
invariants = [
    ('t3', 6, 1),   # f = n-1-2 = 6
    ('t5', 4, 1),   # f = n-1-4 = 4
    ('t7', 2, 1),   # f = n-1-6 = 2
    ('t9', 0, 1),   # f = n-1-8 = 0
    ('bc', 4, 2),   # f = n-1-2*2 = 4
    ('bc35', 2, 2), # f = n-1-2-4 = 2  (weighted by cycle count product)
    ('bc37', 0, 2), # f = n-1-2-6 = 0
    ('a3', 2, 3),   # f = n-1-3*2 = 2
]

print(f"DEFORMED EULERIAN VERIFICATION at n={n}")
print("=" * 70)

print(f"\nInflated coefficients c_k^(f, {d}):")
for name, f, parts in invariants:
    coeffs = [2**parts * inflated_eulerian(f, d, k) for k in range(d + 1)]
    print(f"  {name:5s} (f={f}, 2^{parts}): {coeffs}")

print(f"\nBaseline A({n}, k): {[eulerian_number(n, k) for k in range(n)]}")

# ====================================================================
# Verify
# ====================================================================
print(f"\n{'=' * 70}")
print("VERIFICATION (5 random n=9 tournaments)")
print("=" * 70)

all_match = True
for trial in range(5):
    print(f"\n  Trial {trial}:")
    A = random_tournament(n, n * 777 + trial)

    # Compute actual distribution
    actual = forward_edge_dist_dp(A, n)
    print(f"    Actual dist: {dict(sorted(actual.items()))}")

    # Compute invariants
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    t9 = count_directed_cycles(A, n, 9)
    bc = count_bc(A, n)
    bc35 = count_bc35(A, n)
    bc37 = count_bc37(A, n)
    a3 = count_alpha3(A, n)

    inv_vals = {'t3': t3, 't5': t5, 't7': t7, 't9': t9,
                'bc': bc, 'bc35': bc35, 'bc37': bc37, 'a3': a3}

    print(f"    Invariants: t3={t3}, t5={t5}, t7={t7}, t9={t9}, bc={bc}, bc35={bc35}, bc37={bc37}, a3={a3}")

    # Predict
    errors = []
    for k in range(n):
        predicted = eulerian_number(n, k)
        for name, f, parts in invariants:
            predicted += 2**parts * inflated_eulerian(f, d, k) * inv_vals[name]
        actual_k = actual.get(k, 0)
        if predicted != actual_k:
            errors.append((k, predicted, actual_k))

    if errors:
        print(f"    ERRORS: {errors}")
        all_match = False
    else:
        # Also verify H via OCF
        H = actual.get(n-1, 0)
        H_ocf = 1 + 2*(t3+t5+t7+t9) + 4*(bc+bc35+bc37) + 8*a3
        print(f"    ALL {n} coefficients MATCH. H={H}, OCF={H_ocf}, OCF_ok={H==H_ocf}")

print(f"\n{'=' * 70}")
print(f"RESULT: {'ALL PASS' if all_match else 'SOME FAIL'}")
print("=" * 70)
