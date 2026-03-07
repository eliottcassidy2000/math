#!/usr/bin/env python3
"""
W-polynomial hierarchy at EVEN n (n=4,6,8).

At even n, W(r) has only ODD powers: W(r) = w_1*r + w_3*r^3 + ... + w_{n-1}*r^{n-1}.
The hierarchy should parallel odd n: w_{n-1} = n!, w_{n-3} = f(t_3), etc.

Key question: does the shift principle carry over to even n?
  C_{t_3}(r) at even n=4: 4r
  C_{t_3}(r) at even n=6: -8r + 48r^3
  C_{t_3}(r) at even n=8: ?

And does C_{t_5}(r) at n=6 = C_{t_3}(r) at n=4?
  C_{t_5}(r) at n=6 = 4r.  C_{t_3}(r) at n=4 = 4r.  YES!

opus-2026-03-06-S30
"""
from itertools import combinations
from math import factorial, comb
from fractions import Fraction
import random
import numpy as np

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

def extract_W_coeffs_even(A, n):
    """Extract odd-power coefficients at even n."""
    num_coeffs = n // 2
    r_sample = [0.1 * (k+1) for k in range(num_coeffs)]
    W_vals = [compute_W_dp(A, n, r) for r in r_sample]
    V = np.array([[r**(2*k+1) for k in range(num_coeffs)] for r in r_sample])
    return np.linalg.solve(V, W_vals)

def count_t3(A, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])

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

def random_tournament(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# =====================================================================
print("=" * 70)
print("W-POLYNOMIAL AT EVEN n")
print("=" * 70)

for n in [4, 6, 8]:
    print(f"\n{'='*50}")
    print(f"n = {n}")
    print(f"{'='*50}")

    num_samples = 30 if n <= 6 else 20
    data = []

    for trial in range(num_samples):
        A = random_tournament(n, n*1000 + trial)
        wc = extract_W_coeffs_even(A, n)
        t3 = count_t3(A, n)
        H = count_H(A, n)
        d = {'wc': wc, 't3': t3, 'H': H}

        if n >= 6:
            d['t5'] = count_t5(A, n)
            d['bc'] = count_bc(A, n)
        if n >= 8:
            d['t7'] = count_t7(A, n)

        data.append(d)

    # Print first few
    for trial in range(min(5, num_samples)):
        d = data[trial]
        coeff_str = " ".join(f"w_{2*k+1}={d['wc'][k]:.1f}" for k in range(n//2))
        print(f"  T{trial}: t3={d['t3']:3d} H={d['H']:5d} {coeff_str}")

    # Verify W(1/2) = H
    print(f"\n  W(1/2) = H check:")
    for d in data[:3]:
        W_half = sum(d['wc'][k] * 0.5**(2*k+1) for k in range(n//2))
        print(f"    W(1/2)={W_half:.2f}  H={d['H']}  diff={abs(W_half - d['H']):.6f}")

    # Regression for each coefficient
    if n == 4:
        inv_names = ['const', 't3']
        X = np.array([[1, d['t3']] for d in data])
    elif n == 6:
        inv_names = ['const', 't3', 't5', 'bc']
        X = np.array([[1, d['t3'], d['t5'], d['bc']] for d in data])
    elif n == 8:
        inv_names = ['const', 't3', 't5', 'bc', 't7']
        X = np.array([[1, d['t3'], d['t5'], d['bc'], d['t7']] for d in data])

    print(f"\n  Coefficient formulas:")
    for k in range(n//2):
        y = np.array([d['wc'][k] for d in data])
        coeffs, res, _, _ = np.linalg.lstsq(X, y, rcond=None)
        y_pred = X @ coeffs
        err = np.max(np.abs(y - y_pred))
        terms = []
        for i, name in enumerate(inv_names):
            frac = Fraction(coeffs[i]).limit_denominator(10000)
            if abs(coeffs[i]) > 0.01:
                terms.append(f"{frac}*{name}")
        power = 2*k + 1
        print(f"    w_{power} = {' + '.join(terms)}  (err={err:.6f})")

    # Per-invariant r-polynomial
    print(f"\n  Per-invariant r-polynomials (ODD powers):")
    for k in range(n//2):
        y = np.array([d['wc'][k] for d in data])
        coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        for i, name in enumerate(inv_names):
            if abs(coeffs[i]) > 0.01:
                frac = Fraction(coeffs[i]).limit_denominator(10000)

    # Check shift principle for even n
    if n == 6:
        # C_{t3}(r) at n=4: 4r
        # C_{t5}(r) at n=6: should also be 4r if shift principle holds
        # C_{t3}(r) at n=6: -8r + 48r^3
        y_w1 = np.array([d['wc'][0] for d in data])
        y_w3 = np.array([d['wc'][1] for d in data])
        y_w5 = np.array([d['wc'][2] for d in data])

        coeffs_w1, _, _, _ = np.linalg.lstsq(X, y_w1, rcond=None)
        coeffs_w3, _, _, _ = np.linalg.lstsq(X, y_w3, rcond=None)

        t5_in_w1 = coeffs_w1[2]
        t5_in_w3 = coeffs_w3[2]
        print(f"\n  SHIFT PRINCIPLE CHECK (even n):")
        print(f"    C_t5(r) at n=6 = {Fraction(t5_in_w1).limit_denominator(100)}*r + {Fraction(t5_in_w3).limit_denominator(100)}*r^3")
        print(f"    C_t3(r) at n=4 = 4*r")
        print(f"    Match: {abs(t5_in_w1 - 4) < 0.01 and abs(t5_in_w3) < 0.01}")

    if n == 8:
        y_coeffs = {}
        for k in range(4):
            y = np.array([d['wc'][k] for d in data])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            y_coeffs[k] = coeffs

        print(f"\n  SHIFT PRINCIPLE CHECK (n=8):")
        # C_{t5}(r) at n=8 should equal C_{t3}(r) at n=6 = -8r + 48r^3
        t5_idx = 2  # index in inv_names for t5
        t5_coeffs = [y_coeffs[k][t5_idx] for k in range(4)]
        print(f"    C_t5(r) at n=8 = {Fraction(t5_coeffs[0]).limit_denominator(100)}*r "
              f"+ {Fraction(t5_coeffs[1]).limit_denominator(100)}*r^3 "
              f"+ {Fraction(t5_coeffs[2]).limit_denominator(100)}*r^5 "
              f"+ {Fraction(t5_coeffs[3]).limit_denominator(100)}*r^7")
        print(f"    C_t3(r) at n=6 = -8*r + 48*r^3")
        match = abs(t5_coeffs[0] - (-8)) < 0.01 and abs(t5_coeffs[1] - 48) < 0.01
        print(f"    Match: {match}")

        # C_{t7}(r) at n=8 should equal C_{t5}(r) at n=6 = 4r
        t7_idx = 4  # index in inv_names for t7
        t7_coeffs = [y_coeffs[k][t7_idx] for k in range(4)]
        print(f"    C_t7(r) at n=8 = {Fraction(t7_coeffs[0]).limit_denominator(100)}*r "
              f"+ {Fraction(t7_coeffs[1]).limit_denominator(100)}*r^3")
        print(f"    C_t5(r) at n=6 = 4*r")
        match2 = abs(t7_coeffs[0] - 4) < 0.01
        print(f"    Match: {match2}")

        # bc shift check
        bc_idx = 3
        bc_coeffs = [y_coeffs[k][bc_idx] for k in range(4)]
        print(f"\n    C_bc(r) at n=8 = {Fraction(bc_coeffs[0]).limit_denominator(100)}*r "
              f"+ {Fraction(bc_coeffs[1]).limit_denominator(100)}*r^3 "
              f"+ {Fraction(bc_coeffs[2]).limit_denominator(100)}*r^5")
        print(f"    C_bc(r) at n=6 = 8*r")

# =====================================================================
print(f"\n{'='*70}")
print("UNIVERSAL PATTERN: even vs odd n")
print("="*70)

# Collect all per-invariant polynomials
print("""
MASTER TABLE OF PER-INVARIANT R-POLYNOMIALS:

ODD n:
  C_t3(r):
    n=5:  -1 + 12r^2
    n=7:  2 - 60r^2 + 240r^4
    n=9:  -17/2 + 462r^2 - 4200r^4 + 10080r^6
  C_t5(r):
    n=5:  2
    n=7:  -1 + 12r^2     [= C_t3(r) at n=5]
    n=9:  2 - 60r^2 + 240r^4  [= C_t3(r) at n=7]
  C_t7(r):
    n=7:  2
    n=9:  -1 + 12r^2     [= C_t3(r) at n=5]
  C_bc(r):
    n=7:  -2 + 24r^2
    n=9:  4 - 120r^2 + 480r^4

EVEN n:
  C_t3(r):
    n=4:  4r
    n=6:  -8r + 48r^3
    n=8:  (computed above)
  C_t5(r):
    n=6:  4r             [= C_t3(r) at n=4]
    n=8:  (should = C_t3(r) at n=6)
  C_bc(r):
    n=6:  8r
    n=8:  (computed above)

SHIFT PRINCIPLE: C_{t_{2j+1}}(r) at n = C_{t_{2j-1}}(r) at n-2
  (same parity of n preserved)
""")

# =====================================================================
print(f"\n{'='*70}")
print("THE MASTER SEQUENCE M_k(r)")
print("="*70)

print("""
Define the master sequence for ODD n:
  M_0(r) = 2                              [cycle using all n vertices]
  M_1(r) = -1 + 12r^2                     [cycle using n-2 vertices]
  M_2(r) = 2 - 60r^2 + 240r^4            [cycle using n-4 vertices]
  M_3(r) = -17/2 + 462r^2 - 4200r^4 + 10080r^6

Properties:
  - M_k(1/2) = 2 for all k  (from OCF: each cycle contributes 2 to H)
  - Leading coefficient of M_k = 2*(2k+1)!
  - M_k has degree 2k in r (only even powers)
  - C_{t_{2j+1}}(r) at n = M_{(n-2j-1)/2}(r)

For EVEN n, the master sequence uses odd powers:
  M_0^e(r) = 4r                           [cycle using all n vertices]
  M_1^e(r) = -8r + 48r^3                  [cycle using n-2 vertices]
  M_2^e(r) = ...                          [to be computed at n=8]

Properties:
  - M_k^e(1/2) = 2 for all k
  - Leading coefficient of M_k^e = 2*(2k+1)! (conjectured)
""")

# Verify M_0^e(1/2) = 2
print(f"M_0^e(1/2) = 4*0.5 = {4*0.5}")
print(f"M_1^e(1/2) = -8*0.5 + 48*0.5^3 = {-8*0.5 + 48*0.125}")

print(f"\n{'='*70}")
print("DONE")
print("="*70)
