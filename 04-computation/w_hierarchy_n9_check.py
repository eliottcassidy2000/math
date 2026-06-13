#!/usr/bin/env python3
"""
Verify the W-coefficient hierarchy at n=9:
  w_8 = 9! = 362880 (universal)
  w_6 = 2*7!*t3 - C(9,3)*7!/2 = 10080*t3 - C(9,3)*5040/2 (predicted)
      = 10080*t3 - 84*2520 = 10080*t3 - 211680

Also verify: at n=9, does w_6 depend ONLY on t3?

kind-pasteur-2026-03-06-S25g
"""
import random
import numpy as np
from itertools import permutations, combinations
from math import factorial, comb

def count_t3(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j]*A[j][k]*A[k][i]: count += 1
                if A[i][k]*A[k][j]*A[j][i]: count += 1
    return count

n = 9
random.seed(42)

print(f"W-COEFFICIENT HIERARCHY at n={n}")
print(f"Predicted w_{{n-3}} = w_6 = 2*(n-2)!*t3 - C(n,3)*(n-2)!/2")
pred_a = 2 * factorial(n-2)
pred_b = -comb(n, 3) * factorial(n-2) // 2
print(f"  = {pred_a}*t3 + {pred_b}")
print("=" * 60)

# At n=9, we can't enumerate all permutations (9! = 362880 per tournament).
# But we can compute W(r) at several r values and fit.

data = []
for trial in range(15):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_t3(A, n)

    # W(r) = sum_P prod(r + s_e) at several r values
    # Use DP to compute W(r) = sum of Hamiltonian path weights
    # DP on (visited_set, last_vertex)
    r_vals = [0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]

    W_vals = []
    for r in r_vals:
        # DP: dp[(S, v)] = sum of prod(r + s_e) over paths ending at v visiting S
        dp = {}
        for v in range(n):
            dp[(1 << v, v)] = 1.0

        for size in range(2, n+1):
            for S in range(1 << n):
                if bin(S).count('1') != size:
                    continue
                for v in range(n):
                    if not (S & (1 << v)):
                        continue
                    S_prev = S ^ (1 << v)
                    total = 0.0
                    for u in range(n):
                        if (S_prev & (1 << u)):
                            s_e = A[u][v] - 0.5
                            val = dp.get((S_prev, u), 0.0)
                            if val != 0:
                                total += val * (r + s_e)
                    if total != 0:
                        dp[(S, v)] = total

            # Clear entries for smaller sizes to save memory
            if size >= 3:
                for S in range(1 << n):
                    if bin(S).count('1') == size - 2:
                        for v in range(n):
                            key = (S, v)
                            if key in dp:
                                del dp[key]

        full = (1 << n) - 1
        W = sum(dp.get((full, v), 0.0) for v in range(n))
        W_vals.append(W)

    # Fit W(r) = w0 + w2*r^2 + w4*r^4 + w6*r^6 + w8*r^8
    # (only even powers for odd n)
    V = np.array([[r**(2*j) for j in range(5)] for r in r_vals])
    coeffs = np.linalg.lstsq(V, np.array(W_vals), rcond=None)[0]
    w0, w2, w4, w6, w8 = coeffs

    data.append({'t3': t3, 'w0': w0, 'w2': w2, 'w4': w4, 'w6': w6, 'w8': w8})

    pred_w6 = pred_a * t3 + pred_b
    err = abs(w6 - pred_w6)

    if trial < 8:
        print(f"  trial={trial:2d}: t3={t3:3d}, w8={w8:.0f}, w6={w6:.1f}, pred_w6={pred_w6}, err={err:.1f}")

# Check w8 = 9!
w8_vals = [d['w8'] for d in data]
print(f"\nw8: min={min(w8_vals):.0f}, max={max(w8_vals):.0f}, expected={factorial(n)}")

# Check w6 = pred_a * t3 + pred_b
max_err_w6 = max(abs(d['w6'] - (pred_a * d['t3'] + pred_b)) for d in data)
print(f"Max |w6 - ({pred_a}*t3 + {pred_b})| = {max_err_w6:.2f}")

if max_err_w6 < 1:
    print(f"\nw_{{n-3}} = (n-2)! * [2*t3 - C(n,3)/2] VERIFIED at n={n}!")
else:
    print(f"\nw_{{n-3}} formula does NOT match at n={n}. Fitting...")
    y = np.array([d['w6'] for d in data])
    X = np.array([[d['t3'], 1] for d in data])
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    print(f"  w6 = {c[0]:.1f}*t3 + {c[1]:.1f}")

# Also check w4 pattern
print(f"\nw4 analysis:")
y4 = np.array([d['w4'] for d in data])
X4 = np.array([[d['t3'], 1] for d in data])
c4 = np.linalg.lstsq(X4, y4, rcond=None)[0]
err4 = max(abs(y4 - X4 @ c4))
print(f"  w4 = {c4[0]:.1f}*t3 + {c4[1]:.1f}  (max err = {err4:.2f})")
if err4 > 1:
    print(f"  w4 at n={n} depends on MORE than t3!")

print("\nDONE")
