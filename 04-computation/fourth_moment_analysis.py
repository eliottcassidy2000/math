#!/usr/bin/env python3
"""
FOURTH MOMENT OF FORWARD-ARC COUNT: sum_P f_P^4

THM-055 tells us tr(c_{n-5}) depends on sum_P f^4 (the 4th moment).
At n=7, tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231 (kind-pasteur).

Question: What is sum_P f^4 in terms of t_3, t_5, bc, and other invariants?

The approach: compute sum_P f^j for j=0..6 across many n=7 tournaments
and regress against (t_3, t_5, bc, H).

opus-2026-03-06-S28
"""

from itertools import permutations, combinations
from math import factorial, comb
import random

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

def ham_count_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1: continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

def is_3cycle(A, triple):
    a, b, c = triple
    return (A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]) > 0

def both_cyclic_total(A, n):
    total = 0
    for v_del in range(n):
        S = [v for v in range(n) if v != v_del]
        for T in combinations(S, 3):
            T_comp = tuple(v for v in S if v not in T)
            if T < T_comp:
                if is_3cycle(A, T) and is_3cycle(A, T_comp):
                    total += 1
    return total

def compute_moments(A, n, max_j=6):
    """Compute sum_P f^j for j=0..max_j."""
    moments = [0] * (max_j + 1)
    for p in permutations(range(n)):
        f = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        fpow = 1
        for j in range(max_j + 1):
            moments[j] += fpow
            fpow *= f
    return moments

# =====================================================================
n = 7
print("=" * 70)
print(f"FOURTH MOMENT ANALYSIS at n={n}")
print("=" * 70)

data = []
for trial in range(25):
    random.seed(trial * 41 + 7)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_3_cycles(A, n)
    t5 = count_5_cycles_dp(A, n)
    bc = both_cyclic_total(A, n)
    H = ham_count_dp(A, n)
    moments = compute_moments(A, n, max_j=6)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    data.append({
        't3': t3, 't5': t5, 'bc': bc, 'H': H, 'scores': scores,
        'moments': moments
    })
    if trial < 10:
        print(f"  T{trial:2d}: t3={t3:2d} t5={t5:2d} bc={bc:2d} H={H:3d} "
              f"m4={moments[4]:10d} m3={moments[3]:9d}")

# =====================================================================
# Regression: sum_P f^4 = a*t_3^2 + b*t_3 + c*t_5 + d*bc + e*H + f_const?
# =====================================================================
print("\n" + "=" * 70)
print("REGRESSION: sum_P f^j against tournament invariants")
print("=" * 70)

import numpy as np

# Build feature matrix
for j in range(2, 7):
    y = np.array([d['moments'][j] for d in data], dtype=float)

    # Try linear in (t3, t5, bc, H)
    X_linear = np.column_stack([
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    try:
        coeffs_lin, res_lin, _, _ = np.linalg.lstsq(X_linear, y, rcond=None)
        max_err_lin = max(abs(X_linear @ coeffs_lin - y))
    except:
        max_err_lin = float('inf')
        coeffs_lin = [0]*5

    # Try with t3^2 added
    X_quad = np.column_stack([
        [d['t3']**2 for d in data],
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [d['H'] for d in data],
        [1 for _ in data],
    ])
    try:
        coeffs_quad, res_quad, _, _ = np.linalg.lstsq(X_quad, y, rcond=None)
        max_err_quad = max(abs(X_quad @ coeffs_quad - y))
    except:
        max_err_quad = float('inf')
        coeffs_quad = [0]*6

    print(f"\n  j={j}: sum_P f^{j}")
    print(f"    Linear (t3,t5,bc,H,1): max_err={max_err_lin:.4f}")
    if max_err_lin < 0.01:
        print(f"    EXACT: {coeffs_lin[0]:.1f}*t3 + {coeffs_lin[1]:.1f}*t5 + "
              f"{coeffs_lin[2]:.1f}*bc + {coeffs_lin[3]:.1f}*H + {coeffs_lin[4]:.1f}")
    print(f"    Quadratic (t3^2,...): max_err={max_err_quad:.4f}")
    if max_err_quad < 0.01:
        print(f"    EXACT: {coeffs_quad[0]:.1f}*t3^2 + {coeffs_quad[1]:.1f}*t3 + "
              f"{coeffs_quad[2]:.1f}*t5 + {coeffs_quad[3]:.1f}*bc + "
              f"{coeffs_quad[4]:.1f}*H + {coeffs_quad[5]:.1f}")

# =====================================================================
# Focus: what does sum_P f^4 depend on that sum_P f^2 and f^3 don't?
# =====================================================================
print("\n" + "=" * 70)
print("WHAT DETERMINES sum_P f^4 BEYOND t_3?")
print("=" * 70)

# Group by t3 value and see how m4 varies
from collections import defaultdict
by_t3 = defaultdict(list)
for d in data:
    by_t3[d['t3']].append(d)

for t3_val in sorted(by_t3.keys()):
    group = by_t3[t3_val]
    if len(group) >= 2:
        m4_vals = [d['moments'][4] for d in group]
        if len(set(m4_vals)) > 1:
            print(f"\n  t3={t3_val}: {len(group)} tournaments, m4 varies!")
            for d in group:
                print(f"    t5={d['t5']:2d} bc={d['bc']:2d} H={d['H']:3d} m4={d['moments'][4]}")

# =====================================================================
# Direct: is m4 = a*t3^2 + b*t3 + c*t5 + d*bc + e?
# (Excluding H since we want to understand what m4 captures intrinsically)
# =====================================================================
print("\n" + "=" * 70)
print("m4 WITHOUT H")
print("=" * 70)

y = np.array([d['moments'][4] for d in data], dtype=float)
X = np.column_stack([
    [d['t3']**2 for d in data],
    [d['t3'] for d in data],
    [d['t5'] for d in data],
    [d['bc'] for d in data],
    [1 for _ in data],
])
coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
pred = X @ coeffs
max_err = max(abs(pred - y))
print(f"  m4 = {coeffs[0]:.2f}*t3^2 + {coeffs[1]:.2f}*t3 + {coeffs[2]:.2f}*t5 + "
      f"{coeffs[3]:.2f}*bc + {coeffs[4]:.2f}")
print(f"  Max error: {max_err:.4f}")

# Try t3^2, t3*t5, t3, t5, bc
X2 = np.column_stack([
    [d['t3']**2 for d in data],
    [d['t3']*d['t5'] for d in data],
    [d['t3'] for d in data],
    [d['t5'] for d in data],
    [d['bc'] for d in data],
    [1 for _ in data],
])
coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
pred2 = X2 @ coeffs2
max_err2 = max(abs(pred2 - y))
print(f"\n  With t3*t5: max_err={max_err2:.4f}")
if max_err2 < 0.01:
    print(f"  EXACT: {coeffs2[0]:.1f}*t3^2 + {coeffs2[1]:.1f}*t3*t5 + "
          f"{coeffs2[2]:.1f}*t3 + {coeffs2[3]:.1f}*t5 + "
          f"{coeffs2[4]:.1f}*bc + {coeffs2[5]:.1f}")

# =====================================================================
# Verify tr(c_2) formula from kind-pasteur
# =====================================================================
print("\n" + "=" * 70)
print("VERIFY tr(c_2) = 24*bc - 60*t_3 + 12*t_5 + 231")
print("=" * 70)

# At n=7, tr(c_2) = sum_P e_4(s_P)
# e_4(f) = f^4/24 - f^3/2 + 47f^2/24 - 11f/4 + 15/16
# So tr(c_2) = m4/24 - m3/2 + 47*m2/24 - 11*m1/4 + 15*m0/16

for d in data[:10]:
    m = d['moments']
    tr_c2_moments = m[4]/24 - m[3]/2 + 47*m[2]/24 - 11*m[1]/4 + 15*m[0]/16
    tr_c2_formula = 24*d['bc'] - 60*d['t3'] + 12*d['t5'] + 231
    diff = abs(tr_c2_moments - tr_c2_formula)
    print(f"  t3={d['t3']:2d} t5={d['t5']:2d} bc={d['bc']:2d}: "
          f"from_moments={tr_c2_moments:.2f} formula={tr_c2_formula} diff={diff:.4f}")

# =====================================================================
# THE KEY QUESTION: Can we express m4 purely in terms of m0..m3 + bc?
# That is, does m4 = polynomial(t_3) + linear(bc)?
# =====================================================================
print("\n" + "=" * 70)
print("KEY: m4 = polynomial(t_3) + alpha*bc?")
print("=" * 70)

# We know m2 = f(t_3), m3 = g(t_3). So try m4 = h(t_3) + alpha*bc
# i.e. m4 = a*t3^2 + b*t3 + c*bc + d
X_key = np.column_stack([
    [d['t3']**2 for d in data],
    [d['t3'] for d in data],
    [d['bc'] for d in data],
    [1 for _ in data],
])
coeffs_key, _, _, _ = np.linalg.lstsq(X_key, y, rcond=None)
pred_key = X_key @ coeffs_key
max_err_key = max(abs(pred_key - y))
print(f"  m4 = {coeffs_key[0]:.2f}*t3^2 + {coeffs_key[1]:.2f}*t3 + "
      f"{coeffs_key[2]:.2f}*bc + {coeffs_key[3]:.2f}")
print(f"  Max error: {max_err_key:.4f}")

if max_err_key > 0.1:
    print(f"  NOT exact — m4 needs more than (t3, bc)")
    # Try adding t5
    X_key2 = np.column_stack([
        [d['t3']**2 for d in data],
        [d['t3'] for d in data],
        [d['t5'] for d in data],
        [d['bc'] for d in data],
        [1 for _ in data],
    ])
    coeffs_key2, _, _, _ = np.linalg.lstsq(X_key2, y, rcond=None)
    pred_key2 = X_key2 @ coeffs_key2
    max_err_key2 = max(abs(pred_key2 - y))
    print(f"  With t5: max_err={max_err_key2:.4f}")
    if max_err_key2 < 0.01:
        print(f"  EXACT: m4 = {coeffs_key2[0]:.1f}*t3^2 + {coeffs_key2[1]:.1f}*t3 + "
              f"{coeffs_key2[2]:.1f}*t5 + {coeffs_key2[3]:.1f}*bc + {coeffs_key2[4]:.1f}")
    elif max_err_key2 > 0.1:
        print(f"  STILL not exact — m4 needs even more")
        # Try t5^2
        X_key3 = np.column_stack([
            [d['t3']**2 for d in data],
            [d['t3'] for d in data],
            [d['t5'] for d in data],
            [d['t5']**2 for d in data],
            [d['bc'] for d in data],
            [1 for _ in data],
        ])
        coeffs_key3, _, _, _ = np.linalg.lstsq(X_key3, y, rcond=None)
        pred_key3 = X_key3 @ coeffs_key3
        max_err_key3 = max(abs(pred_key3 - y))
        print(f"  With t5^2: max_err={max_err_key3:.4f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
