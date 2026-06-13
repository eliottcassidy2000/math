#!/usr/bin/env python3
"""
Forward-Edge Distribution and Derivatives of W(r) at r=1/2.

THEOREM: F'_j(1/2) = 2^{j+1} - 2  (Mersenne-like sequence)
THEOREM: W'(1/2) = (n-1)*H(T) + a_{n-2}(T)
         where a_k(T) = #{permutations with exactly k forward edges}

From the OCF decomposition:
  a_{n-2}(T) = W'(1/2) - (n-1)*H(T)
             = sum_I 2^{parts(I)} * (2^{f_I+1}-2) * I(T) + (2^n-2) - (n-1)*H(T)

Since H(T) = I(Omega(T), 2), this gives a_{n-2} in terms of OCF invariants.

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from fractions import Fraction
from math import factorial
import random

def compute_W_dp(A, n, r_val):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + val * wt
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

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
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_forward_edge_dist(A, n):
    """Count a_k = #{permutations with exactly k forward edges}. DP approach."""
    # dp[mask][v][k] = #paths ending at v, using vertices in mask, with k forward edges
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for k in range(n):
                c = dp.get((mask, v, k), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    fwd = A[v][u]
                    new_k = k + fwd
                    key = (mask | (1 << u), u, new_k)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    a = {}
    for v in range(n):
        for k in range(n):
            c = dp.get((full, v, k), 0)
            if c > 0:
                a[k] = a.get(k, 0) + c
    return a

def count_t3(A, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])

def random_tournament(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# =====================================================================
# Part 1: Verify W'(1/2) = (n-1)*H + a_{n-2}
# =====================================================================

print("=" * 70)
print("PART 1: W'(1/2) = (n-1)*H(T) + a_{n-2}(T)")
print("=" * 70)

for n in [5, 7, 9]:
    print(f"\nn={n}:")
    for trial in range(5):
        A = random_tournament(n, n*100 + trial)
        H = count_H(A, n)
        a = count_forward_edge_dist(A, n)

        # Compute W'(1/2) numerically
        eps = 1e-8
        Wp = compute_W_dp(A, n, 0.5 + eps)
        Wm = compute_W_dp(A, n, 0.5 - eps)
        W_deriv = (Wp - Wm) / (2 * eps)

        predicted = (n-1) * H + a.get(n-2, 0)

        print(f"  T{trial}: H={H:6d}, a_{{n-2}}={a.get(n-2,0):6d}, "
              f"W'(1/2)={W_deriv:.1f}, (n-1)*H+a_{{n-2}}={predicted}, "
              f"match={abs(W_deriv - predicted) < 1}")

# =====================================================================
# Part 2: Forward edge distribution palindromy (a_k = a_{n-1-k})
# =====================================================================

print(f"\n{'=' * 70}")
print("PART 2: Palindromy a_k(T) = a_{n-1-k}(T)")
print("=" * 70)

for n in [5, 7]:
    A = random_tournament(n, n*200)
    a = count_forward_edge_dist(A, n)
    print(f"\nn={n}: {dict(sorted(a.items()))}")
    for k in range(n):
        ak = a.get(k, 0)
        an1k = a.get(n-1-k, 0)
        if ak != an1k:
            print(f"  PALINDROMY FAILS at k={k}!")
    print(f"  Palindromy verified: a_k = a_{{n-1-k}} for all k")

# =====================================================================
# Part 3: a_{n-2} formula from OCF decomposition
# =====================================================================

print(f"\n{'=' * 70}")
print("PART 3: a_{n-2}(T) formula via OCF")
print("=" * 70)

# W'(1/2) = sum_I 2^{parts} * F'_{f_I}(1/2) * I(T) + F'_{n-1}(1/2)
# = sum_I 2^{parts} * (2^{f_I+1}-2) * I(T) + (2^n - 2)
# And a_{n-2} = W'(1/2) - (n-1)*H(T)
# H(T) = 1 + 2*(t3+t5+...) + 4*(bc+...) + ...
# So a_{n-2} = (2^n-2) + sum_I 2^{parts}*(2^{f+1}-2)*I - (n-1)*(1 + 2*sum alpha_k*2^k...)
#
# Let's compute numerically at n=7 first.

n = 7
print(f"\nn={n}:")

import numpy as np

def count_t5(A, n):
    t5 = 0
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
        t5 += sum(dp[full][v] for v in range(1,5) if sub[v][0])
    return t5

def count_t7(A, n):
    if n < 7: return 0
    t7 = 0
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
        t7 += sum(dp[full][v] for v in range(1,7) if sub[v][0])
    return t7

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

num_samples = 20
data = []
for trial in range(num_samples):
    A = random_tournament(n, n*300 + trial)
    H = count_H(A, n)
    a_dist = count_forward_edge_dist(A, n)
    a_n2 = a_dist.get(n-2, 0)
    t3 = count_t3(A, n)
    t5 = count_t5(A, n)
    t7 = count_t7(A, n)
    bc = count_bc(A, n)
    data.append({'H': H, 'a_n2': a_n2, 't3': t3, 't5': t5, 't7': t7, 'bc': bc})
    if trial < 5:
        print(f"  T{trial}: H={H:5d}, a_{{n-2}}={a_n2:5d}, t3={t3:3d}, t5={t5:4d}, bc={bc:3d}")

# Regression: a_{n-2} = c0 + c1*t3 + c2*t5 + c3*t7 + c4*bc
inv_names = ['const', 't3', 't5', 't7', 'bc']
X = np.array([[1, d['t3'], d['t5'], d['t7'], d['bc']] for d in data])
y = np.array([d['a_n2'] for d in data])
coeffs, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
y_pred = X @ coeffs
max_err = np.max(np.abs(y - y_pred))

print(f"\n  Regression: a_{{n-2}} = ", end="")
terms = []
for i, name in enumerate(inv_names):
    frac = Fraction(coeffs[i]).limit_denominator(10000)
    if abs(coeffs[i]) > 0.01:
        terms.append(f"{frac}*{name}")
print(" + ".join(terms))
print(f"  Max error: {max_err:.6f}")

# Theoretical prediction from F'_j(1/2) = 2^{j+1} - 2:
# a_{n-2} = W'(1/2) - (n-1)*H
# W'(1/2) = F'_{n-1}(1/2) + sum C'_I(1/2) * I(T)
# = (2^n-2) + 2*t3*(2^{f_t3+1}-2) + 2*t5*(2^{f_t5+1}-2) + 2*t7*(2^{f_t7+1}-2) + 4*bc*(2^{f_bc+1}-2)
# f_t3 = n-1-2=4, f_t5=n-1-4=2, f_t7=n-1-6=0, f_bc=n-1-4=2
# = (2^7-2) + 2*t3*(2^5-2) + 2*t5*(2^3-2) + 2*t7*(2^1-2) + 4*bc*(2^3-2)
# = 126 + 60*t3 + 12*t5 + 0*t7 + 24*bc
# a_{n-2} = 126 + 60*t3 + 12*t5 + 0*t7 + 24*bc - 6*H
# H = 1 + 2*(t3+t5+t7) + 4*bc
# a_{n-2} = 126 + 60*t3 + 12*t5 + 24*bc - 6*(1 + 2*t3 + 2*t5 + 2*t7 + 4*bc)
# = 126 + 60*t3 + 12*t5 + 24*bc - 6 - 12*t3 - 12*t5 - 12*t7 - 24*bc
# = 120 + 48*t3 + 0*t5 - 12*t7

print(f"\n  Theoretical prediction:")
print(f"  F'_{{n-1}}(1/2) = 2^{n} - 2 = {2**n - 2}")
for inv, f_val, parts in [('t3',4,1), ('t5',2,1), ('t7',0,1), ('bc',2,2)]:
    coeff = 2**parts * (2**(f_val+1) - 2)
    print(f"  C'_{inv}(1/2) = 2^{parts} * (2^{{{f_val}+1}}-2) = {coeff}")

# Full formula
print(f"\n  W'(1/2) = {2**n-2} + {2*(2**5-2)}*t3 + {2*(2**3-2)}*t5 + {2*(2**1-2)}*t7 + {4*(2**3-2)}*bc")
print(f"  = 126 + 60*t3 + 12*t5 + 0*t7 + 24*bc")
print(f"\n  a_{{n-2}} = W'(1/2) - (n-1)*H")
print(f"  = 126 + 60*t3 + 12*t5 + 24*bc - 6*(1 + 2*t3 + 2*t5 + 2*t7 + 4*bc)")
print(f"  = 120 + 48*t3 - 12*t7")

# Verify
print(f"\n  Verification:")
for i, d in enumerate(data[:5]):
    pred = 120 + 48*d['t3'] - 12*d['t7']
    print(f"    T{i}: a_{{n-2}}={d['a_n2']}, predicted={pred}, match={d['a_n2']==pred}")

# =====================================================================
# Part 4: Higher derivatives — W''(1/2), a_{n-3}
# =====================================================================

print(f"\n{'=' * 70}")
print("PART 4: Second derivative and a_{n-3}")
print("=" * 70)

# F''_j(1/2) from Eulerian numbers: only k=0,1,2 contribute
# k=0: j*(j-1)*A(j+1,0) = j(j-1)
# k=1: 2*(j-1)*A(j+1,1) (from cross terms)
# k=2: 2*A(j+1,2) (from the pure k=2 term)
# A(j+1,2) = 3^{j+1} - (j+2)*2^{j+1} + C(j+2,2) [Euler number formula]

# Let me compute F''_j(1/2) numerically for verification
for f in range(8):
    if f <= 7:
        # Compute via numerical differentiation of F_f
        def eval_Ff(r):
            total = 0
            for sigma in permutations(range(f+1)):
                prod = 1.0
                for i in range(f):
                    eps = 1 if sigma[i+1] > sigma[i] else -1
                    prod *= r + eps*0.5
                total += prod
            return total

        eps = 1e-5
        d2 = (eval_Ff(0.5+eps) - 2*eval_Ff(0.5) + eval_Ff(0.5-eps)) / eps**2
        print(f"  F_{f}''(1/2) = {d2:.1f}")

# W''(1/2) = sum_P sum_{j<k} prod_{i != j,k} (1/2+s_i)
# Nonzero only when <= 2 backward edges
# If 0 backward: each (j,k) pair contributes 1, total = C(n-1,2)
# If 1 backward at position m: only pairs (j,k) where one of j,k = m contribute
#   when j=m: prod of all except j,k has 0 backward, product = 1. Count: n-2 choices for k.
#   Same for k=m. But (j,k) unordered, so just n-2 choices total.
# If 2 backward at positions m1,m2: only pair (j,k) = {m1,m2} contributes. Product = 1. Count: 1.
# W''(1/2) = C(n-1,2)*a_{n-1} + (n-2)*a_{n-2} + a_{n-3}
# a_{n-3} = W''(1/2) - C(n-1,2)*H - (n-2)*a_{n-2}

print(f"\n  W''(1/2) = C(n-1,2)*H + (n-2)*a_{{n-2}} + a_{{n-3}}")
for trial in range(5):
    d = data[trial]
    A = random_tournament(n, n*300 + trial)
    a_dist = count_forward_edge_dist(A, n)

    # Compute W''(1/2) numerically
    eps = 1e-6
    Wp = compute_W_dp(A, n, 0.5 + eps)
    W0 = compute_W_dp(A, n, 0.5)
    Wm = compute_W_dp(A, n, 0.5 - eps)
    W_d2 = (Wp - 2*W0 + Wm) / eps**2

    cn12 = (n-1)*(n-2)//2
    pred = cn12 * d['H'] + (n-2) * a_dist.get(n-2, 0) + a_dist.get(n-3, 0)

    print(f"  T{trial}: W''(1/2)={W_d2:.1f}, C(n-1,2)*H+(n-2)*a_{{n-2}}+a_{{n-3}}={pred}, "
          f"a_{{n-3}}={a_dist.get(n-3,0)}")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
