#!/usr/bin/env python3
"""
Test whether W(r) coefficients stratify by odd-cycle complexity.

W(r) = sum_P prod_{e in P} (r + s_e) where s_e = A[u,v] - 1/2

Expanding: W(r) = sum_k w_k * r^k

Known:
  w_{n-1} = H (path count)
  w_{n-2} = 0 (by symmetry)
  w_{n-3} = 2*(n-2)!*t_3 - const (opus S27)

Conjecture: w_{n-5} depends on t_3 and t_5 (5-cycle count)
            w_{n-7} depends on t_3, t_5, t_7

kind-pasteur-2026-03-06-S25f
"""

from itertools import permutations, combinations
from math import factorial, comb
import numpy as np
import random

def tournament_random(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_kcycles(A, k):
    """Count directed k-cycles in tournament A."""
    n = len(A)
    count = 0
    for verts in combinations(range(n), k):
        # Check all cyclic orderings
        from itertools import permutations as perms
        for p in perms(verts):
            # Check if p[0]->p[1]->...->p[k-1]->p[0] is a directed cycle
            is_cycle = True
            for i in range(k):
                if A[p[i]][p[(i+1)%k]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    # Each k-cycle is counted k times (cyclic rotations)
    return count // k

def W_coefficients(A):
    """Compute all coefficients of W(r) = sum_P prod(r + s_e)."""
    n = len(A)
    # W(r) is a polynomial of degree n-1
    coeffs = [0.0] * n  # coeffs[k] = coefficient of r^k

    for p in permutations(range(n)):
        # s_i = A[p[i]][p[i+1]] - 0.5
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]

        # Product (r+s_0)(r+s_1)...(r+s_{n-2})
        # Use polynomial multiplication
        poly = [1.0]  # Start with 1
        for si in s:
            new_poly = [0.0] * (len(poly) + 1)
            for j, c in enumerate(poly):
                new_poly[j+1] += c      # r * current
                new_poly[j] += c * si   # s_i * current
            poly = new_poly

        for k in range(len(poly)):
            if k < n:
                coeffs[k] += poly[k]

    return coeffs

print("=" * 70)
print("W(r) COEFFICIENT STRATIFICATION BY ODD-CYCLE COMPLEXITY")
print("=" * 70)

# Test at n=5
n = 5
print(f"\n--- n={n} ---")
print(f"W(r) has degree {n-1} = {n-1}")
print(f"Coefficients: w_0, w_1, w_2, w_3, w_4")
print(f"Expected: w_4 = H, w_3 = 0, w_2 ~ t_3, w_1 = 0, w_0 ~ ???")

data = []
for seed in range(50):
    A = tournament_random(n, seed * 13)
    coeffs = W_coefficients(A)
    t3 = count_kcycles(A, 3)
    t5 = count_kcycles(A, 5)
    H = int(round(coeffs[n-1]))
    data.append({
        'seed': seed, 'H': H, 't3': t3, 't5': t5,
        'w': [round(c, 6) for c in coeffs]
    })

# Check w_{n-2} = w_3 = 0
print(f"\nw_3 (should be 0): {[d['w'][3] for d in data[:5]]}")

# Fit w_2 = a*t3 + b
X = np.array([[d['t3'], 1] for d in data])
y2 = np.array([d['w'][2] for d in data])
c2 = np.linalg.lstsq(X, y2, rcond=None)[0]
err2 = max(abs(y2 - X @ c2))
print(f"\nw_2 = {c2[0]:.4f}*t3 + {c2[1]:.4f}  (max_err={err2:.6f})")
print(f"  Expected: 2*(n-2)!*t3 - (n-2)!*C(n,3)/2 = {2*factorial(n-2)}*t3 + {-factorial(n-2)*comb(n,3)//2}")

# Check w_1 = 0
print(f"\nw_1 (should be 0): {[d['w'][1] for d in data[:5]]}")

# Fit w_0 = a*t3 + b*t5 + c  or  w_0 = a*t3^2 + b*t3 + c*t5 + d
print(f"\nw_0 analysis:")
y0 = np.array([d['w'][0] for d in data])

# Try w_0 = a*t3 + b
X_t3 = np.array([[d['t3'], 1] for d in data])
c_t3 = np.linalg.lstsq(X_t3, y0, rcond=None)[0]
err_t3 = max(abs(y0 - X_t3 @ c_t3))
print(f"  w_0 = {c_t3[0]:.4f}*t3 + {c_t3[1]:.4f}  (max_err={err_t3:.4f})")

# Try w_0 = a*t3 + b*t5 + c
X_t3t5 = np.array([[d['t3'], d['t5'], 1] for d in data])
c_t3t5 = np.linalg.lstsq(X_t3t5, y0, rcond=None)[0]
err_t3t5 = max(abs(y0 - X_t3t5 @ c_t3t5))
print(f"  w_0 = {c_t3t5[0]:.4f}*t3 + {c_t3t5[1]:.4f}*t5 + {c_t3t5[2]:.4f}  (max_err={err_t3t5:.4f})")

# Try w_0 = a*t3^2 + b*t3 + c
X_t3sq = np.array([[d['t3']**2, d['t3'], 1] for d in data])
c_t3sq = np.linalg.lstsq(X_t3sq, y0, rcond=None)[0]
err_t3sq = max(abs(y0 - X_t3sq @ c_t3sq))
print(f"  w_0 = {c_t3sq[0]:.6f}*t3^2 + {c_t3sq[1]:.4f}*t3 + {c_t3sq[2]:.4f}  (max_err={err_t3sq:.4f})")

# Try w_0 = a*t3^2 + b*t3 + c*t5 + d
X_full = np.array([[d['t3']**2, d['t3'], d['t5'], 1] for d in data])
c_full = np.linalg.lstsq(X_full, y0, rcond=None)[0]
err_full = max(abs(y0 - X_full @ c_full))
print(f"  w_0 = {c_full[0]:.6f}*t3^2 + {c_full[1]:.4f}*t3 + {c_full[2]:.4f}*t5 + {c_full[3]:.4f}  (max_err={err_full:.4f})")

if err_full < 0.01:
    print(f"  EXACT FIT with t3^2, t3, t5!")
    if abs(c_full[2]) > 0.01:
        print(f"  t5 CONTRIBUTES to w_0!")
    else:
        print(f"  t5 does NOT contribute (t3^2 + t3 suffice)")

# Now do n=7
n = 7
print(f"\n--- n={n} ---")
print(f"W(r) has degree {n-1} = 6")
print(f"Coefficients: w_0 through w_6")

data7 = []
for seed in range(20):
    A = tournament_random(n, seed * 17)
    coeffs = W_coefficients(A)
    t3 = count_kcycles(A, 3)
    t5 = count_kcycles(A, 5)
    H = int(round(coeffs[n-1]))
    data7.append({
        'seed': seed, 'H': H, 't3': t3, 't5': t5,
        'w': [round(c, 4) for c in coeffs]
    })
    print(f"  seed={seed}: H={H}, t3={t3}, t5={t5}, w=[{', '.join(f'{c:.1f}' for c in coeffs)}]")

# Check odd-indexed coeffs are 0
print(f"\nOdd coefficients (should be 0):")
print(f"  w_5: {[d['w'][5] for d in data7[:5]]}")
print(f"  w_3: {[d['w'][3] for d in data7[:5]]}")
print(f"  w_1: {[d['w'][1] for d in data7[:5]]}")

# Fit w_4 = a*t3 + b
y4 = np.array([d['w'][4] for d in data7])
X7 = np.array([[d['t3'], 1] for d in data7])
c4 = np.linalg.lstsq(X7, y4, rcond=None)[0]
err4 = max(abs(y4 - X7 @ c4))
print(f"\nw_4 = {c4[0]:.4f}*t3 + {c4[1]:.4f}  (max_err={err4:.4f})")
print(f"  Expected: 2*(n-2)!*t3 - ... = {2*factorial(n-2)}*t3 + ...")

# Fit w_2 with t3, t5, t3^2
y2_7 = np.array([d['w'][2] for d in data7])
X7_full = np.array([[d['t3']**2, d['t3'], d['t5'], 1] for d in data7])
c2_7 = np.linalg.lstsq(X7_full, y2_7, rcond=None)[0]
err2_7 = max(abs(y2_7 - X7_full @ c2_7))
print(f"\nw_2 = {c2_7[0]:.4f}*t3^2 + {c2_7[1]:.4f}*t3 + {c2_7[2]:.4f}*t5 + {c2_7[3]:.4f}  (max_err={err2_7:.4f})")

if abs(c2_7[2]) > 0.01:
    print(f"  t5 CONTRIBUTES to w_2! Coefficient = {c2_7[2]:.4f}")
    print(f"  This confirms w_{n-5} depends on 5-cycle count!")
else:
    print(f"  t5 does NOT contribute to w_2")

# Fit w_0 similarly
y0_7 = np.array([d['w'][0] for d in data7])

# Try various models
# w_0 = a*t3^3 + b*t3^2 + c*t3*t5 + d*t3 + e*t5 + f*t7 + g
# First check t3-only models
X7_t3 = np.array([[d['t3'], 1] for d in data7])
c0_t3 = np.linalg.lstsq(X7_t3, y0_7, rcond=None)[0]
err0_t3 = max(abs(y0_7 - X7_t3 @ c0_t3))
print(f"\nw_0 = {c0_t3[0]:.4f}*t3 + {c0_t3[1]:.4f}  (max_err={err0_t3:.4f})")

X7_t3t5 = np.array([[d['t3'], d['t5'], 1] for d in data7])
c0_t3t5 = np.linalg.lstsq(X7_t3t5, y0_7, rcond=None)[0]
err0_t3t5 = max(abs(y0_7 - X7_t3t5 @ c0_t3t5))
print(f"w_0 = {c0_t3t5[0]:.4f}*t3 + {c0_t3t5[1]:.4f}*t5 + {c0_t3t5[2]:.4f}  (max_err={err0_t3t5:.4f})")

X7_quad = np.array([[d['t3']**2, d['t3'], d['t5'], 1] for d in data7])
c0_quad = np.linalg.lstsq(X7_quad, y0_7, rcond=None)[0]
err0_quad = max(abs(y0_7 - X7_quad @ c0_quad))
print(f"w_0 = {c0_quad[0]:.6f}*t3^2 + {c0_quad[1]:.4f}*t3 + {c0_quad[2]:.4f}*t5 + {c0_quad[3]:.4f}  (max_err={err0_quad:.4f})")

X7_cubic = np.array([[d['t3']**3, d['t3']**2, d['t3']*d['t5'], d['t3'], d['t5'], 1] for d in data7])
c0_cubic = np.linalg.lstsq(X7_cubic, y0_7, rcond=None)[0]
err0_cubic = max(abs(y0_7 - X7_cubic @ c0_cubic))
print(f"w_0 = {c0_cubic[0]:.8f}*t3^3 + {c0_cubic[1]:.6f}*t3^2 + {c0_cubic[2]:.6f}*t3*t5 + ...")
print(f"      {c0_cubic[3]:.4f}*t3 + {c0_cubic[4]:.4f}*t5 + {c0_cubic[5]:.4f}  (max_err={err0_cubic:.4f})")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key question: Does w_{n-2k-1} depend on cycle counts t_3, t_5, ..., t_{2k+1}?

If YES: W(r) is a "cycle-stratified generating function" where each
        coefficient level introduces the next odd-cycle count.
        This would mean the W(r) polynomial ENCODES the full odd-cycle
        structure of the tournament, stratified by "complexity."

Connection to OCF: H = W(1/2) = I(Omega(T), 2)
  The OCF says H decomposes into alpha_k * 2^k.
  The W(r) expansion says H decomposes into w_k * (1/2)^k.
  These are DUAL decompositions of the same number!
""")
print("DONE")
