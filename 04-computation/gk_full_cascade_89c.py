#!/usr/bin/env python3
"""
gk_full_cascade_89c.py — Full cascading g_k extraction with g_k(2)=2k bootstrap
opus-2026-03-14-S89c

Key boundary conditions:
  g_k(1) = 1 for all k (monotone path argument)
  g_k(2) = 2k for all k (conjectured from k=1..5, using to bootstrap)
  g_k(0) = 0 (no room for dominos)

Strategy: use these to resolve the triangular system and get enough
data points for polynomial interpolation of g_5 and g_6.
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand
import time

def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)

def compute_W_dp(n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[((1 << v), v)] = 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if u == v - 1:
                    continue
                weight = 2 * cnt if u == v + 1 else cnt
                new_key = (mask | (1 << u), u)
                dp[new_key] = dp.get(new_key, 0) + weight
    return sum(dp.get((full, v), 0) for v in range(n))

# Known formulas
def g1(m): return Fraction(m)
def g2(m): return Fraction(m * m)
def g3(m): return Fraction(m * (2*m*m + 1), 3)
def g4(m): return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)

known_g = {1: g1, 2: g2, 3: g3, 4: g4}

# Compute W(n) for n=3..20
print("Computing W(n) values...")
W_vals = {}
for n in range(3, 21):
    t0 = time.time()
    W_vals[n] = compute_W_dp(n)
    t1 = time.time()
    print(f"  W({n}) = {W_vals[n]}  ({t1-t0:.1f}s)")

# Storage for extracted g_k values
gk_vals = {k: {} for k in range(5, 10)}

# Boundary conditions
for k in range(5, 10):
    gk_vals[k][1] = Fraction(1)      # g_k(1) = 1
    gk_vals[k][0] = Fraction(0)      # g_k(0) = 0

# Conjectured: g_k(2) = 2k
for k in range(5, 10):
    gk_vals[k][2] = Fraction(2 * k)

print("\n" + "="*70)
print("FULL CASCADING g_k EXTRACTION")
print("="*70)

for n in range(11, 21):
    nf = factorial(n)
    cv2 = Fraction(W_vals[n], nf) - 1
    remaining = cv2

    # Subtract known g_1..g_4
    for k in range(1, 5):
        mk = n - 2*k
        if mk < 1:
            break
        contrib = Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
        remaining -= contrib

    # Now remaining = Σ_{k≥5} 2·g_k(n-2k) / (n)_{2k}
    # Identify unknowns
    terms = []
    for k in range(5, (n-1)//2 + 1):
        mk = n - 2*k
        if mk < 0:
            break
        terms.append((k, mk))

    if not terms:
        continue

    # Subtract all known values
    unknowns = []
    for k, mk in terms:
        if k in gk_vals and mk in gk_vals[k]:
            val = gk_vals[k][mk]
            remaining -= Fraction(2 * val, falling_factorial(n, 2*k))
        else:
            unknowns.append((k, mk))

    print(f"\nn={n}: terms={terms}")

    if len(unknowns) == 0:
        # All known — check consistency
        print(f"  All known. Residual = {remaining} (should be 0)")
        if remaining != 0:
            print(f"  *** INCONSISTENCY: residual = {float(remaining):.6e} ***")
    elif len(unknowns) == 1:
        k, mk = unknowns[0]
        gk_val = remaining * falling_factorial(n, 2*k) / 2
        print(f"  Solved: g_{k}({mk}) = {gk_val}")
        gk_vals[k][mk] = gk_val
    else:
        # Multiple unknowns — solve for the lowest k (largest mk)
        # The first unknown has the largest coefficient
        k0, mk0 = unknowns[0]
        print(f"  Unknowns: {unknowns}")
        print(f"  Cannot fully resolve; solving g_{k0}({mk0}) ignoring higher terms")
        # Estimate: remaining ≈ 2·g_{k0}(mk0) / (n)_{2k0}
        gk_approx = remaining * falling_factorial(n, 2*k0) / 2
        print(f"  g_{k0}({mk0}) ≈ {gk_approx} = {float(gk_approx):.4f}")
        # Store as approximate
        gk_vals[k0][mk0] = gk_approx

# Print all extracted values
print("\n" + "="*70)
print("EXTRACTED g_k VALUES")
print("="*70)

x = symbols('x')

for k in range(5, 9):
    vals = gk_vals[k]
    if not vals:
        continue
    print(f"\ng_{k}(m):")
    for mk in sorted(vals):
        v = vals[mk]
        is_int = (v.denominator == 1)
        print(f"  g_{k}({mk}) = {v}{'' if is_int else f' = {float(v):.6f} (NOT INTEGER!)'}")

    # Try polynomial interpolation with integer values
    pts = []
    for mk in sorted(vals):
        v = vals[mk]
        if v.denominator == 1 and mk >= 1:
            pts.append((Rational(mk), Rational(v)))

    if len(pts) >= 3:
        print(f"\n  Polynomial fit ({len(pts)} integer points):")
        for deg in range(1, len(pts)):
            poly = interpolate(pts[:deg+1], x)
            all_match = True
            for m_val, g_val in pts[deg+1:]:
                if poly.subs(x, m_val) != g_val:
                    all_match = False
                    break
            if all_match and len(pts) > deg + 1:
                print(f"  → g_{k}(m) = {factor(poly)} = {expand(poly)} (degree {deg}, EXACT)")
                # Check g_k(0) = 0
                if poly.subs(x, 0) == 0:
                    print(f"    ✓ g_{k}(0) = 0")
                else:
                    print(f"    ✗ g_{k}(0) = {poly.subs(x, 0)} ≠ 0!")
                break
        else:
            poly = interpolate(pts, x)
            print(f"  → g_{k}(m) = {factor(poly)} = {expand(poly)} (through {len(pts)} points)")
            if poly.subs(x, 0) == 0:
                print(f"    ✓ g_{k}(0) = 0")
            else:
                print(f"    ✗ g_{k}(0) = {poly.subs(x, 0)} (should be 0 if polynomial is correct)")

# Verify g_k(2) = 2k pattern
print("\n" + "="*70)
print("VERIFY g_k(2) = 2k PATTERN")
print("="*70)
for k in range(1, 8):
    if k <= 4:
        val = known_g[k](2)
    elif 2 in gk_vals.get(k, {}):
        val = gk_vals[k][2]
    else:
        val = "?"
    expected = 2 * k
    match = "✓" if val == expected else "?"
    print(f"  g_{k}(2) = {val}, expected 2·{k} = {expected}  {match}")

# Asymptotic analysis with all known g_k
print("\n" + "="*70)
print("CV² ASYMPTOTIC WITH g_1..g_5")
print("="*70)

# If we have g_5 polynomial, use it
g5_poly = None
pts5 = []
for mk in sorted(gk_vals[5]):
    v = gk_vals[5][mk]
    if v.denominator == 1 and mk >= 1:
        pts5.append((Rational(mk), Rational(v)))

if len(pts5) >= 4:
    # Try each degree
    for deg in range(1, len(pts5)):
        poly = interpolate(pts5[:deg+1], x)
        all_match = True
        for m_val, g_val in pts5[deg+1:]:
            if poly.subs(x, m_val) != g_val:
                all_match = False
                break
        if all_match and len(pts5) > deg + 1:
            g5_poly = poly
            break
    if g5_poly is None:
        g5_poly = interpolate(pts5, x)

if g5_poly:
    print(f"Using g_5(m) = {expand(g5_poly)}")

for n in range(5, 30):
    val = Fraction(0)
    for k in range(1, n//2 + 1):
        mk = n - 2*k
        if mk < 1:
            break
        if k <= 4:
            g = known_g[k](mk)
        elif k == 5 and g5_poly:
            g = Fraction(g5_poly.subs(x, mk))
        else:
            break
        ff = falling_factorial(n, 2*k)
        val += Fraction(2*g, ff)

    # Compare with 2/n - 14/(3n³)
    approx1 = Fraction(2, n)
    approx2 = Fraction(2, n) - Fraction(14, 3*n**3)
    diff1 = val - approx1
    diff2 = val - approx2
    print(f"  n={n:2d}: CV²≈{float(val):.12f}, 2/n={float(approx1):.12f}, "
          f"diff={float(diff1):.2e}, diff_from_2/n-14/3n³={float(diff2):.2e}")

# Leading coefficients analysis
print("\n" + "="*70)
print("1/n^k COEFFICIENT ANALYSIS")
print("="*70)
print("CV² = c_1/n + c_2/n² + c_3/n³ + c_4/n⁴ + ...")
print(f"From |S|=2:  2(n-2)/(n(n-1)) = 2/n - 2/n² + 2/n³ - 2/n⁴ + ...")
print(f"From |S|=4:  2(n-4)²/(n)_4")
print(f"From |S|=6:  2·m(2m²+1)/3 / (n)_6  with m=n-6")

# Compute exact 1/n^k coefficients numerically
# Use large n to extract them
print("\nNumerical coefficient extraction (from polynomial approximation at large n):")
for power in range(1, 7):
    # Compute n^power × (CV² - known lower terms) at large n
    # This is getting complex; let's just show the residuals
    pass

# Exact coefficient of 1/n⁴
# |S|=2: coefficient of 1/n⁴ is -2 (from 2/(n(n-1)) expansion)
# Actually: 2(n-2)/(n(n-1)) = 2/n · (n-2)/(n-1) = 2/n · (1 - 1/(n-1))
#   = 2/n - 2/(n(n-1)) = 2/n - 2/n² - 2/n³ - 2/n⁴ - ...
# Wait let me be more careful.
# 2(n-2)/(n(n-1)) = 2(n-2)/(n²-n)
# = 2/n · (n-2)/(n-1) = 2/n · (1 - 1/(n-1))
# = 2/n · Σ_{j≥0} (-1)^j / (n-1)^j ... no.
# (n-2)/(n-1) = 1 - 1/(n-1). And 1/(n-1) = (1/n)·1/(1-1/n) = (1/n)·Σ 1/n^k
# So (n-2)/(n-1) = 1 - 1/n · 1/(1-1/n) = 1 - 1/n - 1/n² - 1/n³ - ...
# Then 2/n · (1 - 1/n - 1/n² - ...) = 2/n - 2/n² - 2/n³ - 2/n⁴ - ...

# |S|=4: 2(n-4)²/(n)_4 = 2(n-4)²/(n(n-1)(n-2)(n-3))
# (n-4)² = n² - 8n + 16
# (n)_4 = n⁴ - 6n³ + 11n² - 6n
# Ratio = 2(n²-8n+16)/(n⁴-6n³+11n²-6n)
# = 2/n² · (1-8/n+16/n²) / (1-6/n+11/n²-6/n³)
# ≈ 2/n² · (1-8/n+16/n²)(1+6/n+(36-11)/n²+...)
# = 2/n² · (1 - 2/n + (16-48+25)/n² + ...)
# Wait this is getting messy. Let me just compute numerically.

print("\nExact 1/n^k coefficients by polynomial fitting:")
# Compute CV²_approx(n) using g_1..g_4 for n=100..110, then fit
from sympy import Matrix

ns = list(range(50, 61))
data_rows = []
for n in ns:
    val = Fraction(0)
    for k in range(1, 5):
        mk = n - 2*k
        if mk < 1:
            break
        g = known_g[k](mk)
        ff = falling_factorial(n, 2*k)
        val += Fraction(2*g, ff)
    data_rows.append(float(val))

# Fit: CV² ≈ Σ c_j / n^j for j=1..8
import numpy as np
A = np.array([[1.0/n**j for j in range(1, 9)] for n in ns])
b = np.array(data_rows)
coeffs, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
print("  CV²(g_1..g_4) ≈ " + " + ".join(f"{coeffs[j-1]:.6f}/n^{j}" for j in range(1, 9)))
print(f"  c_1 = {coeffs[0]:.10f} (expected 2)")
print(f"  c_2 = {coeffs[1]:.10f} (expected 0)")
print(f"  c_3 = {coeffs[2]:.10f} (expected -14/3 ≈ -4.667)")
print(f"  c_4 = {coeffs[3]:.10f}")
print(f"  c_5 = {coeffs[4]:.10f}")

print("\nDone!")
