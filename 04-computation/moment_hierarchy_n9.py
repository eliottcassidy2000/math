#!/usr/bin/env python3
"""
MOMENT HIERARCHY at n=9: what new invariants enter?

At n=9, n-1=8 edges per permutation. Expected hierarchy:
  - m0, m1: universal
  - m2, m3: f(t3)
  - m4, m5: f(t3, t5, bc)
  - m6, m7: f(t3, t5, t7, bc, ?, H)?
  - m8: f(t3, ..., H) with coefficient 8! for H

9! = 362880 permutations — feasible.

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

def count_5_cycles_dp(A, n):
    count = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for mask in range(1, 1 << 5):
            for v in range(5):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(5):
                    if mask & (1 << u): continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << 5) - 1
        hc = sum(dp[full][v] for v in range(1, 5) if sub[v][0])
        count += hc
    return count

def count_7_cycles_dp(A, n):
    count = 0
    for verts in combinations(range(n), 7):
        sub = [[A[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for mask in range(1, 1 << 7):
            for v in range(7):
                if not (mask & (1 << v)) or dp[mask][v] == 0: continue
                for u in range(7):
                    if mask & (1 << u): continue
                    if sub[v][u]:
                        dp[mask | (1 << u)][u] += dp[mask][v]
        full = (1 << 7) - 1
        hc = sum(dp[full][v] for v in range(1, 7) if sub[v][0])
        count += hc
    return count

def count_9_cycles_dp(A, n):
    """Count 9-cycles = Hamiltonian cycles."""
    dp = [[0]*n for _ in range(1 << n)]
    dp[1][0] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(1, n) if A[v][0]) // n

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

def is_3cycle(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total(A, n):
    """bc = sum over (n-1)-vertex subsets of #complementary-cyclic-triple partitions."""
    total = 0
    for S in combinations(range(n), 6):
        S_list = list(S)
        for T in combinations(S_list, 3):
            T_comp = tuple(v for v in S_list if v not in T)
            if T < T_comp:
                if is_3cycle(A, T) and is_3cycle(A, T_comp):
                    total += 1
    return total

def compute_moments_fast(A, n, max_j=8):
    """Compute moments using DP over permutations."""
    moments = [0] * (max_j + 1)
    # We need to iterate over all permutations — use a generator
    for p in permutations(range(n)):
        f = 0
        for i in range(n-1):
            f += A[p[i]][p[i+1]]
        fpow = 1
        for j in range(max_j + 1):
            moments[j] += fpow
            fpow *= f
    return moments

# =====================================================================
n = 9
print("=" * 70)
print(f"MOMENT HIERARCHY at n={n}")
print("=" * 70)

data = []
for trial in range(15):
    random.seed(trial * 73 + 9)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    H = ham_count_dp(A, n)
    bc = both_cyclic_total(A, n)

    # Moments (fast for n=9 — 362880 perms)
    moments = compute_moments_fast(A, n, max_j=8)

    d = {
        't3': t3, 't5': t5, 'H': H, 'bc': bc,
        'moments': moments
    }
    data.append(d)
    print(f"  T{trial:2d}: t3={t3:3d} t5={t5:4d} bc={bc:4d} H={H:5d}")

# =====================================================================
# Regression: what determines each moment?
# =====================================================================
print("\n" + "=" * 70)
print("MOMENT REGRESSION")
print("=" * 70)

for j in range(2, 9):
    y = np.array([d['moments'][j] for d in data], dtype=float)

    # Linear in (t3, t5, bc, H, 1)
    X = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coeffs
    max_err = max(abs(pred - y))

    print(f"\n  j={j}: sum_P f^{j}")
    print(f"    f(t3,t5,bc,H): max_err={max_err:.2f}")
    if max_err < 1:
        cstr = []
        names = ['t3', 't5', 'bc', 'H', '1']
        for i, name in enumerate(names):
            if abs(coeffs[i]) > 0.5:
                cstr.append(f"{coeffs[i]:.0f}*{name}")
        print(f"    EXACT: {' + '.join(cstr)}")

    # Also try t3-only
    X_t3 = np.column_stack([[d['t3'] for d in data], [1 for _ in data]])
    c_t3, _, _, _ = np.linalg.lstsq(X_t3, y, rcond=None)
    err_t3 = max(abs(X_t3 @ c_t3 - y))
    if err_t3 < 1:
        print(f"    [t3 only: {c_t3[0]:.0f}*t3 + {c_t3[1]:.0f}]")

    # Try (t3, t5, bc) without H
    X_no_H = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [1 for _ in data],
    ])
    c_no_H, _, _, _ = np.linalg.lstsq(X_no_H, y, rcond=None)
    err_no_H = max(abs(X_no_H @ c_no_H - y))
    if err_no_H < 1 and max_err > 1:
        print(f"    WITHOUT H: max_err={err_no_H:.2f} — EXACT")

    if max_err > 1:
        # Not determined by (t3, t5, bc, H)
        # Try adding t3^2, t3*t5, t3*bc
        X_quad = np.column_stack([
            [d['t3']**2 for d in data],
            [d['t3']*d['t5'] for d in data],
            [d['t3'] for d in data],
            [d['t5'] for d in data],
            [d['bc'] for d in data],
            [d['H'] for d in data],
            [1 for _ in data],
        ])
        c_quad, _, _, _ = np.linalg.lstsq(X_quad, y, rcond=None)
        err_quad = max(abs(X_quad @ c_quad - y))
        print(f"    Quadratic (t3^2, t3*t5, ...): max_err={err_quad:.2f}")

# =====================================================================
# Check: at n=9, does the (n-1)! coefficient for H hold in m8?
# =====================================================================
print("\n" + "=" * 70)
print(f"COEFFICIENT OF H IN m8: should be 8! = {factorial(8)}")
print("=" * 70)

y8 = np.array([d['moments'][8] for d in data], dtype=float)
y8_adj = y8 - factorial(8) * np.array([d['H'] for d in data], dtype=float)

# Does m8 - 8!*H depend on (t3, t5, bc)?
X_no_H_2 = np.column_stack([
    [d['t3'] for d in data],
    [d['t5'] for d in data],
    [d['bc'] for d in data],
    [1 for _ in data],
])
c_adj, _, _, _ = np.linalg.lstsq(X_no_H_2, y8_adj, rcond=None)
err_adj = max(abs(X_no_H_2 @ c_adj - y8_adj))
print(f"  m8 - 8!*H = a*t3 + b*t5 + c*bc + d")
print(f"  Coefficients: {c_adj}")
print(f"  Max error: {err_adj:.2f}")

if err_adj > 1:
    print(f"  NOT exact with (t3, t5, bc) — needs more invariants at n=9!")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
