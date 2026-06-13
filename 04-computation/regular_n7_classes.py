#!/usr/bin/env python3
"""
Find the exact invariants of all regular n=7 tournament classes.
Verify p_0(2) and W(i/2) for each class.

kind-pasteur-2026-03-07-S28
"""
from itertools import combinations
import random

def make_tournament(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
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

def count_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            c = dp.get((mask, v), 0)
            if c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[(mask|(1<<u), u)] = dp.get((mask|(1<<u), u), 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_bc(A, n):
    cyc3 = []
    for t in combinations(range(n), 3):
        if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]:
            cyc3.append(set(t))
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3)) if cyc3[i].isdisjoint(cyc3[j]))

def compute_W(A, n, r):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = complex(1, 0)
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + val * (r + A[v][u] - 0.5)
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


n = 7
classes = {}
examples = {}
rng = random.Random(42)
for trial in range(10000):
    bits = rng.randint(0, (1 << 21) - 1)
    A = make_tournament(n, bits)
    scores = sorted([sum(A[i]) for i in range(n)])
    if scores != [3] * 7:
        continue
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)
    H = count_H(A, n)
    key = (H, t3, t5, t7, bc)
    if key not in classes:
        classes[key] = 0
        examples[key] = A
    classes[key] += 1

print("Regular n=7 tournament classes (from 10K random samples):")
print("=" * 80)
for key in sorted(classes.keys()):
    H, t3, t5, t7, bc = key
    alpha1 = t3 + t5 + t7
    # p_0(2) formula from opus's reduced polynomial
    p0 = 2176 + 2 * (-128 * t3 + 16 * t5 - 8 * t7) + 64 * bc
    W_formula = -p0 / 8  # W(i/2) = (-1)^3/2^3 * p_0(2)

    # Verify W(i/2) directly
    A = examples[key]
    W_direct = compute_W(A, n, 0.5j)

    print(f"H={H}: t3={t3}, t5={t5}, t7={t7}, bc={bc}")
    print(f"  alpha1={alpha1}, alpha2={bc}")
    print(f"  p_0(2)={p0}, W(i/2) formula={W_formula:.1f}, W(i/2) direct={W_direct.real:.1f}")
    print(f"  match: {abs(W_formula - W_direct.real) < 0.1}")
    print(f"  count in sample: {classes[key]}")
    print()
