#!/usr/bin/env python3
"""
CREATIVE EXPLORATION: The w_0 penalty shift pattern.

Key discovery: H - w_0 penalizes DIFFERENT cycle types at different n:
  n=5: H - w_0 = 3*t_3                     (penalizes smallest cycles)
  n=7: H - w_0 = 3*t_5 + 6*alpha_2 + 21/4  (penalizes MIDDLE cycles, not smallest)

This is reminiscent of:
1. Euler characteristic: chi = V - E + F (alternating signs by dimension)
2. Jones polynomial at t=0: determinant (alternating sum of cycle data)
3. Mobius function: alternating sum over poset lattice

Question: Is there a UNIVERSAL pattern?
  H - w_0 = sum_k c_k(n) * alpha_k  for some coefficients c_k(n)?

Let me check: does W(r) at r=0 admit a clean formula in terms of OCF alpha's?

W(1/2) = H = sum_k alpha_k * 2^k
W(0) = w_0

Define: Delta = H - w_0 = W(1/2) - W(0)

This is the INTEGRAL of W'(r) from 0 to 1/2.
Since W(r) is a polynomial of degree n-1, this integral has a closed form.

But we also know W(r) = sum_{j=0}^{n-1} w_j * r^j (even j only for odd n).

So Delta = sum_{j>0} w_j * (1/2)^j = W(1/2) - w_0.

The question is: can we express Delta in terms of cycle data?

At n=5: Delta = 3*t_3
  w_2 = 12*t_3 - 30, w_4 = 120
  Delta = w_2/4 + w_4/16 = (12*t_3 - 30)/4 + 120/16 = 3*t_3 - 7.5 + 7.5 = 3*t_3. CHECK!

At n=7: Delta = 3*t_5 + 6*alpha_2 + 21/4
  W(r) has coefficients w_0, w_2, w_4, w_6
  Delta = w_2/4 + w_4/16 + w_6/64

Let me compute w_2, w_4, w_6 at n=7 and verify.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
from fractions import Fraction
import random

def count_directed_k_cycles(A, n, k):
    count = 0
    for verts in combinations(range(n), k):
        seen = set()
        for p in permutations(verts):
            if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                start = min(range(k), key=lambda i: p[i])
                canon = tuple(p[start:] + p[:start])
                if canon not in seen:
                    seen.add(canon)
                    count += 1
    return count

def count_alpha2(A, n):
    all_cycles = []
    for k in range(3, n+1, 2):
        for verts in combinations(range(n), k):
            seen = set()
            for p in permutations(verts):
                if all(A[p[i]][p[(i+1)%k]] == 1 for i in range(k)):
                    start = min(range(k), key=lambda i: p[i])
                    canon = tuple(p[start:] + p[:start])
                    if canon not in seen:
                        seen.add(canon)
                        all_cycles.append(frozenset(verts))
    count = 0
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            if all_cycles[i].isdisjoint(all_cycles[j]):
                count += 1
    return count

n = 7
random.seed(42)

print("=" * 70)
print("W(r) COEFFICIENT STRUCTURE AT n=7")
print("=" * 70)

import numpy as np

for trial in range(10):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_directed_k_cycles(A, n, 3)
    t5 = count_directed_k_cycles(A, n, 5)
    t7 = count_directed_k_cycles(A, n, 7)
    a2 = count_alpha2(A, n)

    # Compute W(r) at several points to extract coefficients
    r_vals = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    W_vals = []
    for r in r_vals:
        W = 0.0
        for p in permutations(range(n)):
            prod = 1.0
            for i in range(n-1):
                s_e = A[p[i]][p[i+1]] - 0.5
                prod *= (r + s_e)
            W += prod
        W_vals.append(W)

    # Fit polynomial: W(r) = w0 + w2*r^2 + w4*r^4 + w6*r^6
    # (odd coefficients = 0 for odd n)
    V = np.array([[r**(2*j) for j in range(4)] for r in r_vals])
    coeffs = np.linalg.lstsq(V, np.array(W_vals), rcond=None)[0]
    w0, w2, w4, w6 = coeffs

    # Compute Delta = W(1/2) - W(0)
    H_computed = W_vals[r_vals.index(0.5)]
    delta = H_computed - w0
    delta_pred = w2/4 + w4/16 + w6/64

    if trial < 5:
        print(f"\n  trial={trial}: t3={t3}, t5={t5}, t7={t7}, a2={a2}")
        print(f"    W coefficients: w0={w0:.4f}, w2={w2:.4f}, w4={w4:.4f}, w6={w6:.4f}")
        print(f"    H={H_computed:.0f}, Delta={delta:.4f}, Delta_pred={delta_pred:.4f}")
        print(f"    3*t5+6*a2+21/4 = {3*t5+6*a2+21/4:.4f}")

# Now fit w2 as function of cycle data
random.seed(42)
data = []
for trial in range(20):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3 = count_directed_k_cycles(A, n, 3)
    t5 = count_directed_k_cycles(A, n, 5)
    t7 = count_directed_k_cycles(A, n, 7)
    a2 = count_alpha2(A, n)

    r_vals = [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
    W_vals = []
    for r in r_vals:
        W = 0.0
        for p in permutations(range(n)):
            prod = 1.0
            for i in range(n-1):
                s_e = A[p[i]][p[i+1]] - 0.5
                prod *= (r + s_e)
            W += prod
        W_vals.append(W)

    V = np.array([[r**(2*j) for j in range(4)] for r in r_vals])
    coeffs = np.linalg.lstsq(V, np.array(W_vals), rcond=None)[0]
    data.append({
        't3': t3, 't5': t5, 't7': t7, 'a2': a2,
        'w0': coeffs[0], 'w2': coeffs[1], 'w4': coeffs[2], 'w6': coeffs[3]
    })

# Fit w2 = a*t3 + b*t5 + c*t7 + d*a2 + e
y = np.array([d['w2'] for d in data])
X = np.array([[d['t3'], d['t5'], d['t7'], d['a2'], 1] for d in data])
c = np.linalg.lstsq(X, y, rcond=None)[0]
err = max(abs(y - X @ c))
print(f"\n{'='*70}")
print(f"w2 = {c[0]:.4f}*t3 + {c[1]:.4f}*t5 + {c[2]:.4f}*t7 + {c[3]:.4f}*a2 + {c[4]:.4f}")
print(f"Max error: {err:.6f}")

# Fit w4
y4 = np.array([d['w4'] for d in data])
c4 = np.linalg.lstsq(X, y4, rcond=None)[0]
err4 = max(abs(y4 - X @ c4))
print(f"\nw4 = {c4[0]:.4f}*t3 + {c4[1]:.4f}*t5 + {c4[2]:.4f}*t7 + {c4[3]:.4f}*a2 + {c4[4]:.4f}")
print(f"Max error: {err4:.6f}")

# w6 should be universal (= 7! = 5040)
w6_vals = [d['w6'] for d in data]
print(f"\nw6 values: min={min(w6_vals):.1f}, max={max(w6_vals):.1f} (should be 5040)")

print(f"\n{'='*70}")
print("HIERARCHY OF W(r) COEFFICIENTS AT n=7")
print(f"{'='*70}")
print(f"  w6 = 7! = 5040 (universal)")
print(f"  w4 depends on: {['t3' if abs(c4[0])>0.1 else '', 't5' if abs(c4[1])>0.1 else '', 't7' if abs(c4[2])>0.1 else '', 'a2' if abs(c4[3])>0.1 else '']}")
print(f"  w2 depends on: {['t3' if abs(c[0])>0.1 else '', 't5' if abs(c[1])>0.1 else '', 't7' if abs(c[2])>0.1 else '', 'a2' if abs(c[3])>0.1 else '']}")
print(f"  w0 depends on: t3, t5, t7, a2")

print("\nDONE")
