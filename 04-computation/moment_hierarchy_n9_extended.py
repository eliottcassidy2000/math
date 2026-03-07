#!/usr/bin/env python3
"""
MOMENT HIERARCHY at n=9 — extended with t7, t9, and bc variants.

At n=9, m6 is NOT determined by (t3, t5, bc, H).
Need to find what additional invariants enter.

Candidates:
  - t7 = #directed 7-cycles
  - t9 = #directed 9-cycles (= Ham cycles)
  - bc4 = sum over 8-vertex subsets of #(disjoint 4-vertex pair with both H>0)
    (analogous to bc for disjoint 3-cycle pairs on 6-vertex subsets)

opus-2026-03-06-S28
"""
from itertools import permutations, combinations
from math import factorial, comb
import random
import numpy as np

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

def count_k_cycles_dp(A, n, k):
    """Count directed k-cycles using DP over k-vertex subsets."""
    count = 0
    for verts in combinations(range(n), k):
        sub = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
        dp = [[0]*k for _ in range(1 << k)]
        dp[1][0] = 1
        for mask in range(1, 1 << k):
            for v in range(k):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(k):
                    if mask & (1 << u): continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << k) - 1
        hc = sum(dp[full][v] for v in range(1, k) if sub[v][0])
        count += hc
    return count

def ham_count_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if (mask & (1 << u)) or not A[v][u]: continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

def ham_count_sub(A, verts):
    """Ham path count for subtournament on given vertices."""
    k = len(verts)
    if k <= 1: return 1
    sub = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    return ham_count_dp(sub, k)

def is_3cycle(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total(A, n):
    """bc = sum over 6-vertex subsets of #complementary-cyclic-triple-partitions."""
    total = 0
    for S in combinations(range(n), 6):
        S_list = list(S)
        for T in combinations(S_list, 3):
            T_comp = tuple(v for v in S_list if v not in T)
            if T < T_comp:
                if is_3cycle(A, T) and is_3cycle(A, T_comp):
                    total += 1
    return total

def bc4_total(A, n):
    """bc4 = sum over 8-vertex subsets of #(unordered pairs of disjoint 4-vertex groups
    where both have at least one directed Ham path)."""
    total = 0
    for S in combinations(range(n), 8):
        S_list = list(S)
        for G1 in combinations(S_list, 4):
            G2 = tuple(v for v in S_list if v not in G1)
            if G1 < G2:
                h1 = ham_count_sub(A, list(G1))
                h2 = ham_count_sub(A, list(G2))
                if h1 > 0 and h2 > 0:
                    total += 1
    return total

def hp_product_sum(A, n):
    """sum over 8-vertex subsets S, over unordered (G1,G2) partition of S into
    two 4-vertex groups, of H(T[G1]) * H(T[G2])."""
    total = 0
    for S in combinations(range(n), 8):
        S_list = list(S)
        for G1 in combinations(S_list, 4):
            G2 = tuple(v for v in S_list if v not in G1)
            if G1 < G2:
                h1 = ham_count_sub(A, list(G1))
                h2 = ham_count_sub(A, list(G2))
                total += h1 * h2
    return total

def compute_moments_fast(A, n, max_j=8):
    moments = [0] * (max_j + 1)
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]])
        fpow = 1
        for j in range(max_j + 1):
            moments[j] += fpow
            fpow *= f
    return moments

# =====================================================================
n = 9
print("=" * 70)
print(f"EXTENDED MOMENT HIERARCHY at n={n}")
print("=" * 70)

data = []
for trial in range(12):
    random.seed(trial * 73 + 9)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_k_cycles_dp(A, n, 5)
    t7 = count_k_cycles_dp(A, n, 7)
    t9 = count_k_cycles_dp(A, n, 9)
    H = ham_count_dp(A, n)
    bc = both_cyclic_total(A, n)
    hpp = hp_product_sum(A, n)

    moments = compute_moments_fast(A, n, max_j=8)

    d = {
        't3': t3, 't5': t5, 't7': t7, 't9': t9,
        'H': H, 'bc': bc, 'hpp': hpp,
        'moments': moments
    }
    data.append(d)
    print(f"  T{trial:2d}: t3={t3:3d} t5={t5:4d} t7={t7:5d} t9={t9:5d} "
          f"bc={bc:4d} hpp={hpp:6d} H={H:5d}")

# =====================================================================
print("\n" + "=" * 70)
print("REGRESSION WITH EXTENDED INVARIANTS")
print("=" * 70)

for j in [6, 7, 8]:
    y = np.array([d['moments'][j] for d in data], dtype=float)

    # (t3, t5, bc, H, 1) — same as n=7
    X_old = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    c_old, _, _, _ = np.linalg.lstsq(X_old, y, rcond=None)
    err_old = max(abs(X_old @ c_old - y))

    # Add t7
    X_t7 = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['t7'] for d in data],
        [d['bc'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    c_t7, _, _, _ = np.linalg.lstsq(X_t7, y, rcond=None)
    err_t7 = max(abs(X_t7 @ c_t7 - y))

    # Add t7, hpp
    X_full = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['t7'] for d in data],
        [d['bc'] for d in data],
        [d['hpp'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    c_full, _, _, _ = np.linalg.lstsq(X_full, y, rcond=None)
    err_full = max(abs(X_full @ c_full - y))

    # Add t7, t9, hpp
    X_all = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['t7'] for d in data],
        [d['t9'] for d in data],
        [d['bc'] for d in data],
        [d['hpp'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    c_all, _, _, _ = np.linalg.lstsq(X_all, y, rcond=None)
    err_all = max(abs(X_all @ c_all - y))

    print(f"\n  j={j}: sum_P f^{j}")
    print(f"    (t3,t5,bc,H): err={err_old:.2f}")
    print(f"    + t7:          err={err_t7:.2f}")
    print(f"    + t7,hpp:      err={err_full:.2f}")
    print(f"    + t7,t9,hpp:   err={err_all:.2f}")

    if err_full < 1:
        names = ['t3', 't5', 't7', 'bc', 'hpp', 'H', '1']
        cstr = []
        for i, name in enumerate(names):
            if abs(c_full[i]) > 0.5:
                cstr.append(f"{c_full[i]:.0f}*{name}")
        print(f"    EXACT: {' + '.join(cstr)}")
    elif err_all < 1:
        names = ['t3', 't5', 't7', 't9', 'bc', 'hpp', 'H', '1']
        cstr = []
        for i, name in enumerate(names):
            if abs(c_all[i]) > 0.5:
                cstr.append(f"{c_all[i]:.0f}*{name}")
        print(f"    EXACT: {' + '.join(cstr)}")

# =====================================================================
# Check m8 = 8!*H + ...
# =====================================================================
print("\n" + "=" * 70)
print(f"m8 - 8!*H regression")
print("=" * 70)
y8_adj = np.array([d['moments'][8] - factorial(8) * d['H'] for d in data], dtype=float)
X_ext = np.column_stack([
    [d['t3'] for d in data],
    [d['t5'] for d in data],
    [d['t7'] for d in data],
    [d['t9'] for d in data],
    [d['bc'] for d in data],
    [d['hpp'] for d in data],
    [1 for _ in data],
])
c_ext, _, _, _ = np.linalg.lstsq(X_ext, y8_adj, rcond=None)
err_ext = max(abs(X_ext @ c_ext - y8_adj))
print(f"  m8 - 8!*H = f(t3,t5,t7,t9,bc,hpp): err={err_ext:.2f}")
if err_ext < 1:
    names = ['t3', 't5', 't7', 't9', 'bc', 'hpp', '1']
    for i, name in enumerate(names):
        if abs(c_ext[i]) > 0.5:
            print(f"    {c_ext[i]:.0f}*{name}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
