#!/usr/bin/env python3
"""
ocr_hidden_geometry.py — kind-pasteur-2026-03-21-S16d

IS THERE A HIDDEN PARABOLA?

The OCR data:
  n:     5       6       7       8       9       10      11
  OCR:   0.9474  0.9231  0.9160  0.9160  0.9189  0.9228  0.9267
  1-OCR: 0.0526  0.0769  0.0840  0.0840  0.0811  0.0772  0.0733

The 1-OCR values look like they could trace a parabola with minimum at n=7.5.

Key exact values:
  1-OCR(5) = 1/19
  1-OCR(6) = 1/13
  1-OCR(7) = 11/131
  1-OCR(8) = 11/131

And we proved: Cov(H, S2) / E[H] = -(n-1)/2 (exact for n=4..7).

R^2 = Cov^2 / (Var_S2 * Var_H)
    = [(n-1)/2 * E[H]]^2 / (Var_S2 * Var_H)
    = (n-1)^2/4 * E[H]^2 / (Var_S2 * Var_H)

So: 1 - R^2 = 1 - (n-1)^2 * E[H]^2 / (4 * Var_S2 * Var_H)

Now: E[H] = n!/2^{n-1}
     Var_S2 = ? (need formula)

If Cov(H,S2)/E[H] = -(n-1)/2 holds for ALL n, then:
  R^2 = (n-1)^2 * E[H]^2 / (4 * Var_S2 * Var_H)

  1 - R^2 = 1 - (n-1)^2 / (4 * CV_S2^2 * CV_H^2)

where CV_X^2 = Var(X)/E[X]^2.

This means the V-shape comes from the ratio of CV^2 values.

Let me compute EVERYTHING as exact rationals and look for the parabola.
"""
from fractions import Fraction
from math import factorial, comb, gcd
import sys

sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  IS THERE A HIDDEN PARABOLA IN 1-OCR?")
print("  kind-pasteur-2026-03-21-S16d")
print("=" * 72)

# Exact data
exact = {
    5:  {'E_H': Fraction(15,2), 'Var_H': Fraction(285,16),
         'Cov': Fraction(-45,4), 'OCR': Fraction(18,19)},
    6:  {'E_H': Fraction(45,2), 'Var_H': Fraction(585,4),
         'Cov': Fraction(-45,1), 'OCR': Fraction(12,13)},
    7:  {'E_H': Fraction(315,4), 'Var_H': Fraction(206325,128),
         'Cov': Fraction(-1575,8), 'OCR': Fraction(120,131)},
    8:  {'E_H': Fraction(315,1), 'Var_H': Fraction(371385,16),
         'Cov': None, 'OCR': Fraction(120,131)},
}

# Compute Var_S2 from the exact R^2 and other moments
# R^2 = Cov^2 / (Var_S2 * Var_H)
# => Var_S2 = Cov^2 / (R^2 * Var_H)

print("\n  PART 1: Extract exact Var(S2) from R^2 identity")
print()

for n in [5, 6, 7]:
    d = exact[n]
    Var_S2 = d['Cov']**2 / (d['OCR'] * d['Var_H'])
    print(f"  n={n}: Var(S2) = {Var_S2} = {float(Var_S2):.6f}")
    # Also compute from Var(S2) = n(n-1)(n-3)/8 (my earlier wrong formula)?
    # Let me just record the exact values.
    exact[n]['Var_S2'] = Var_S2

# At n=8, we need Var_S2 and Cov separately. From the C output:
# sum_H = 84557168640, sum_4S2 = 15032385536
# sum_H2 = 32866314485760, sum_16S22 = 1022202216448
# sum_4HS2 = 3720515420160
# N = 2^28 = 268435456

N8 = 268435456
sum_H_8 = 84557168640
sum_4S2_8 = 15032385536
sum_H2_8 = 32866314485760
sum_16S22_8 = 1022202216448
sum_4HS2_8 = 3720515420160

E_H_8 = Fraction(sum_H_8, N8)
E_4S2_8 = Fraction(sum_4S2_8, N8)
Var_H_8 = Fraction(sum_H2_8, N8) - E_H_8**2
Var_4S2_8 = Fraction(sum_16S22_8, N8) - E_4S2_8**2
Cov_H_4S2_8 = Fraction(sum_4HS2_8, N8) - E_H_8 * E_4S2_8

# True Var(S2) = Var(4*S2)/16, Cov(H,S2) = Cov(H,4*S2)/4
Var_S2_8 = Var_4S2_8 / 16
Cov_HS2_8 = Cov_H_4S2_8 / 4

exact[8]['Cov'] = Cov_HS2_8
exact[8]['Var_S2'] = Var_S2_8

print(f"\n  n=8: Var(S2) = {Var_S2_8} = {float(Var_S2_8):.6f}")
print(f"  n=8: Cov(H,S2) = {Cov_HS2_8} = {float(Cov_HS2_8):.6f}")
print(f"  n=8: E[H] = {E_H_8} = {float(E_H_8):.6f}")

# Check Cov/E[H] = -(n-1)/2 at n=8
cov_ratio_8 = Cov_HS2_8 / E_H_8
print(f"  n=8: Cov/E[H] = {cov_ratio_8} = {float(cov_ratio_8):.6f}")
print(f"  Expected -(n-1)/2 = {-(8-1)/2}")

# ============================================================
# PART 2: The key ratio — CV^2(H) as a function of n
# ============================================================

print("\n" + "=" * 72)
print("  PART 2: CV^2(H) = Var(H) / E[H]^2")
print("=" * 72)
print()

cv2_H = {}
for n in [5, 6, 7, 8]:
    d = exact[n]
    cv2 = d['Var_H'] / d['E_H']**2
    cv2_H[n] = cv2
    print(f"  n={n}: CV^2(H) = {cv2} = {float(cv2):.10f}")

# Also: CV^2(S2) = Var(S2) / E[S2]^2
print()
cv2_S2 = {}
for n in [5, 6, 7, 8]:
    d = exact[n]
    E_S2 = Fraction(n*(n-1), 4)
    cv2 = d['Var_S2'] / E_S2**2
    cv2_S2[n] = cv2
    print(f"  n={n}: CV^2(S2) = {cv2} = {float(cv2):.10f}")

# ============================================================
# PART 3: R^2 in terms of CV ratios
# ============================================================

print("\n" + "=" * 72)
print("  PART 3: R^2 = (n-1)^2 / [4 * n^2(n-1)^2/16 * CV^2(S2) * ... ]")
print("=" * 72)

# R^2 = Cov^2 / (Var_S2 * Var_H)
# If Cov = -(n-1)/2 * E[H]:
#   Cov^2 = (n-1)^2/4 * E[H]^2
#   R^2 = (n-1)^2/4 * E[H]^2 / (Var_S2 * Var_H)
#       = 1 / (4/(n-1)^2 * Var_S2/E[H]^2 * Var_H)
# Hmm, let me be more careful.

# R^2 = [(n-1)/2]^2 * E[H]^2 / (Var_S2 * Var_H)
# 1-R^2 = 1 - [(n-1)/2]^2 / (Var_S2/E[H]^2 * Var_H)
# Wait, this doesn't simplify cleanly. Let me just compute 1-R^2 directly.

print()
for n in [5, 6, 7, 8]:
    d = exact[n]
    one_minus = 1 - d['OCR']
    print(f"  n={n}: 1-OCR = {one_minus} = {float(one_minus):.10f}")

# ============================================================
# PART 4: Fit a parabola to 1-OCR
# ============================================================

print("\n" + "=" * 72)
print("  PART 4: Fit parabola to 1-OCR")
print("=" * 72)

# Points: (5, 1/19), (6, 1/13), (7, 11/131), (8, 11/131)
# Plus sampled: (9, 0.0811), (10, 0.0772), (11, 0.0733)

import numpy as np

ns = np.array([5, 6, 7, 8, 9, 10, 11], dtype=float)
residuals = np.array([1/19, 1/13, 11/131, 11/131, 0.0811, 0.0772, 0.0733])

# Fit quadratic: 1-OCR = a*(n-c)^2 + d
# = a*n^2 - 2ac*n + ac^2 + d
# Equivalently: y = A*n^2 + B*n + C

X = np.column_stack([ns**2, ns, np.ones(7)])
coefs = np.linalg.lstsq(X, residuals, rcond=None)[0]
A, B, C = coefs

pred = X @ coefs
ss_res = np.sum((residuals - pred)**2)
ss_tot = np.sum((residuals - np.mean(residuals))**2)
r2_fit = 1 - ss_res / ss_tot

vertex_n = -B / (2*A)
vertex_y = C - B**2 / (4*A)

print(f"\n  Quadratic fit: 1-OCR = {A:.6f}*n^2 + {B:.6f}*n + {C:.6f}")
print(f"  R^2 of fit = {r2_fit:.8f}")
print(f"  Vertex at n = {vertex_n:.4f}, 1-OCR = {vertex_y:.6f}")
print(f"  Residuals: {residuals - pred}")
print()

# The vertex should be near n=7.5 (midpoint of 7 and 8)
print(f"  Parabola vertex at n = {vertex_n:.4f}")
print(f"  This is {'close to' if abs(vertex_n - 7.5) < 0.5 else 'NOT near'} n = 7.5")

# Extrapolate
print(f"\n  Extrapolated 1-OCR:")
for n_ext in [12, 15, 20, 50, 100]:
    pred_val = A * n_ext**2 + B * n_ext + C
    print(f"    n={n_ext:4d}: 1-OCR = {pred_val:.6f}, OCR = {1-pred_val:.6f}")

# ============================================================
# PART 5: Is the parabola in a TRANSFORMED variable?
# ============================================================

print("\n" + "=" * 72)
print("  PART 5: Parabola in transformed coordinates?")
print("=" * 72)

# Maybe 1-OCR is quadratic in 1/n, or in log(n), or in m = C(n,2)

# Try 1/n
inv_ns = 1.0 / ns
X2 = np.column_stack([inv_ns**2, inv_ns, np.ones(7)])
coefs2 = np.linalg.lstsq(X2, residuals, rcond=None)[0]
pred2 = X2 @ coefs2
r2_2 = 1 - np.sum((residuals - pred2)**2) / ss_tot
print(f"\n  Quadratic in 1/n: R^2 = {r2_2:.8f}")

# Try m = C(n,2)
ms = ns * (ns - 1) / 2
X3 = np.column_stack([ms**2, ms, np.ones(7)])
coefs3 = np.linalg.lstsq(X3, residuals, rcond=None)[0]
pred3 = X3 @ coefs3
r2_3 = 1 - np.sum((residuals - pred3)**2) / ss_tot
print(f"  Quadratic in m=C(n,2): R^2 = {r2_3:.8f}")

# Try (n-7.5)^2 directly
centered = (ns - 7.5)**2
X4 = np.column_stack([centered, np.ones(7)])
coefs4 = np.linalg.lstsq(X4, residuals, rcond=None)[0]
pred4 = X4 @ coefs4
r2_4 = 1 - np.sum((residuals - pred4)**2) / ss_tot
print(f"  Linear in (n-7.5)^2: R^2 = {r2_4:.8f}")
print(f"    1-OCR = {coefs4[0]:.6f} * (n-7.5)^2 + {coefs4[1]:.6f}")

# Try (n-7)*(n-8)
prod78 = (ns - 7) * (ns - 8)
X5 = np.column_stack([prod78, np.ones(7)])
coefs5 = np.linalg.lstsq(X5, residuals, rcond=None)[0]
pred5 = X5 @ coefs5
r2_5 = 1 - np.sum((residuals - pred5)**2) / ss_tot
print(f"  Linear in (n-7)(n-8): R^2 = {r2_5:.8f}")
print(f"    1-OCR = {coefs5[0]:.6f} * (n-7)(n-8) + {coefs5[1]:.6f}")

# Try CV^2(H) directly — we know these exactly
# CV^2(H) denominators: 285/16 / (15/2)^2 = 285/16 * 4/225 = 285/(16*225/4) = 19/60
# So CV^2(H)(5) = 19/60, CV^2(H)(6) = 13/45, CV^2(H)(7) = 131/504, CV^2(H)(8) = 131/560

print("\n" + "=" * 72)
print("  PART 6: The denominators of CV^2(H)")
print("=" * 72)

for n in [3, 4, 5, 6, 7, 8]:
    E_H = Fraction(factorial(n), 2**(n-1))
    if n in exact:
        Var_H = exact[n]['Var_H']
    elif n == 3:
        Var_H = Fraction(3, 4)
    elif n == 4:
        Var_H = Fraction(3)
    else:
        continue
    cv2 = Var_H / E_H**2
    print(f"  n={n}: CV^2(H) = {cv2}")

    # Factor numerator and denominator
    def factor(x):
        if x == 0: return {}
        f = {}
        d = 2
        x = abs(x)
        while d*d <= x:
            while x % d == 0:
                f[d] = f.get(d,0)+1
                x //= d
            d += 1
        if x > 1: f[x] = f.get(x,0)+1
        return f

    print(f"    Numerator {cv2.numerator}: {factor(cv2.numerator)}")
    print(f"    Denominator {cv2.denominator}: {factor(cv2.denominator)}")

# ============================================================
# PART 7: The EXACT formula for R^2
# ============================================================

print("\n" + "=" * 72)
print("  PART 7: If Cov = -(n-1)/2 * E[H], what is R^2?")
print("=" * 72)

# R^2 = Cov^2 / (Var_S2 * Var_H)
#      = (n-1)^2/4 * E[H]^2 / (Var_S2 * Var_H)
#      = (n-1)^2/4 / (Var_S2/E[H]^2 * Var_H/1)   -- wrong
# Let me redo:
# R^2 = [(n-1)^2/4 * (n!/2^{n-1})^2] / [Var_S2 * Var_H]

# Now: Var_S2 depends only on score distribution (computable from arc independence)
# And Var_H depends on permutation pair structure.

# The OCR denominator (prime p) satisfies:
# p | Var_H but p ∤ [(n-1)^2 * E[H]^2 / Var_S2]
# So p must divide Var_H and NOT divide the rest.

# Let's verify: Var_S2 at each n
for n in [5, 6, 7, 8]:
    d = exact[n]
    print(f"\n  n={n}:")
    print(f"    Var_S2 = {d['Var_S2']}")
    numer = (n-1)**2 * d['E_H']**2
    print(f"    (n-1)^2 * E[H]^2 = {numer}")
    print(f"    (n-1)^2 * E[H]^2 / 4 = {numer / 4}")
    print(f"    Var_S2 * Var_H = {d['Var_S2'] * d['Var_H']}")
    print(f"    R^2 = {numer/4} / {d['Var_S2'] * d['Var_H']} = {d['OCR']}")

# ============================================================
# PART 8: Var(S2) formula
# ============================================================

print("\n" + "=" * 72)
print("  PART 8: What is Var(S2) exactly?")
print("=" * 72)

for n in [5, 6, 7, 8]:
    d = exact[n]
    VS = d['Var_S2']
    print(f"  n={n}: Var(S2) = {VS}")

# Check: Var(S2) = n(n-1)(n-3)/8?
# n=5: 5*4*2/8 = 5, but actual Var_S2(5) = ?
# From R^2 = 18/19: Var_S2 = Cov^2 / (R^2 * Var_H) = (45/4)^2 / (18/19 * 285/16)
# = 2025/16 / (5130/304) = 2025/16 * 304/5130 = 2025*19/(16*270) = 2025*19/4320
# = 38475/4320 = 8.90625... hmm

VS5 = Fraction(45,4)**2 / (Fraction(18,19) * Fraction(285,16))
print(f"\n  Computed Var(S2) at n=5: {VS5} = {float(VS5):.6f}")
print(f"  n(n-1)(n-3)/8 at n=5: {5*4*2/8}")

# These don't match. The issue is that S2 = sum (s_i - (n-1)/2)^2
# and the scores s_i are NOT independent (they share arcs).
# So Var(S2) is not given by the naive formula.

# Let me just record the exact values.
print("\n  Exact Var(S2) values:")
for n in [5, 6, 7]:
    d = exact[n]
    print(f"  n={n}: {d['Var_S2']}")

# At n=8:
print(f"  n=8: {exact[8]['Var_S2']}")

# Factor these
for n in [5, 6, 7, 8]:
    d = exact[n]
    VS = d['Var_S2']
    def factor(x):
        if x <= 0: return {}
        f = {}
        d_val = 2
        x = abs(int(x))
        while d_val*d_val <= x:
            while x % d_val == 0:
                f[d_val] = f.get(d_val,0)+1
                x //= d_val
            d_val += 1
        if x > 1: f[x] = f.get(x,0)+1
        return f

    print(f"  n={n}: Var(S2) = {VS}, num={VS.numerator}, den={VS.denominator}")
    print(f"    num factors: {factor(VS.numerator)}")
    print(f"    den factors: {factor(VS.denominator)}")
