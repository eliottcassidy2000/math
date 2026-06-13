#!/usr/bin/env python3
"""
gk_extend_89c.py — Extend g_k computation to larger n via DP
opus-2026-03-14-S89c

Use bitmask DP to compute W(n), then extract g_k values by
subtracting known lower-order contributions.
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand
import time

def falling_factorial(n, k):
    return reduce(lambda a,b: a*b, range(n, n-k, -1), 1)

def compute_W_dp(n):
    """Bitmask DP for W(n) = Σ_{σ∈NUD(n)} 2^{adj1(σ)}."""
    full = (1 << n) - 1
    dp = {}

    for v in range(n):
        dp[((1 << v), v)] = 1

    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            cnt = dp.get(key, 0)
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

    W = sum(dp.get((full, v), 0) for v in range(n))
    return W

# Known g_k formulas
def g1(m):
    return m

def g2(m):
    return m * m

def g3(m):
    return m * (2*m*m + 1) // 3

# g4: conjecture from 4 data points
# g_4(m) = (10m³-33m²+50m-24)/3
def g4(m):
    val = 10*m**3 - 33*m**2 + 50*m - 24
    assert val % 3 == 0
    return val // 3

# Compute W(n) for n=11..16 using DP
print("="*70)
print("W(n) AND g_k EXTRACTION")
print("="*70)

gk_data = {1: [], 2: [], 3: [], 4: [], 5: [], 6: []}

for n in range(3, 17):
    t0 = time.time()
    W = compute_W_dp(n)
    t1 = time.time()
    nf = factorial(n)
    cv2 = Fraction(W, nf) - 1
    print(f"\nn={n}: W={W}, time={t1-t0:.1f}s")
    print(f"  CV² = {float(cv2):.14f}")

    # Subtract known contributions
    remaining = cv2
    m = n - 1  # max positions

    for k in range(1, n // 2 + 1):
        size = 2 * k
        if size > m:
            break

        mk = n - size  # m value for this k

        if k == 1:
            contrib = Fraction(2 * g1(mk), falling_factorial(n, size))
        elif k == 2 and mk >= 1:
            contrib = Fraction(2 * g2(mk), falling_factorial(n, size))
        elif k == 3 and mk >= 1:
            contrib = Fraction(2 * g3(mk), falling_factorial(n, size))
        elif k == 4 and mk >= 1:
            contrib = Fraction(2 * g4(mk), falling_factorial(n, size))
        else:
            # Unknown g_k, this is what we're extracting
            # Record remaining as upper bound
            ff = falling_factorial(n, size)
            gk_approx = remaining * ff / 2
            gk_data[k].append((n, mk, gk_approx, remaining))
            print(f"  g_{k}({mk}) ≈ {gk_approx} = {float(gk_approx):.6f} (remaining={float(remaining):.2e})")
            break

        remaining -= contrib

    if remaining != 0 and n - 2*(n//2) == 0:
        # n is even, remaining should be 0 if we computed everything
        pass

# Now analyze g_5 and g_6
print("\n" + "="*70)
print("g_5 AND g_6 ANALYSIS")
print("="*70)

x = symbols('x')

for k in [5, 6]:
    data = gk_data.get(k, [])
    if not data:
        continue
    print(f"\ng_{k}:")
    for n, mk, g_approx, rem in data:
        # Check if g_approx is close to integer
        g_rounded = round(float(g_approx))
        diff = abs(float(g_approx) - g_rounded)
        print(f"  m={mk}: g ≈ {float(g_approx):.6f} (nearest int: {g_rounded}, diff={diff:.2e})")

    # If values look like integers, extract them
    integer_pts = []
    for n, mk, g_approx, rem in data:
        g_rounded = round(float(g_approx))
        if abs(float(g_approx) - g_rounded) < 0.01:
            integer_pts.append((Rational(mk), Rational(g_rounded)))

    if len(integer_pts) >= 2:
        poly = interpolate(integer_pts[:min(len(integer_pts), 5)], x)
        print(f"  Interpolation: g_{k}(m) = {factor(poly)} = {expand(poly)}")

# Also verify g_4 with extended data
print("\n" + "="*70)
print("g_4 VERIFICATION")
print("="*70)

g4_verified = [(1, 1), (2, 8), (3, 33), (4, 96)]
for n, mk, g_approx, rem in gk_data.get(4, []):
    g_rounded = round(float(g_approx))
    pred = g4(mk)
    print(f"  m={mk}: predicted g_4={pred}, computed≈{float(g_approx):.6f} (diff={abs(pred-float(g_approx)):.2e})")
    if abs(pred - float(g_approx)) < 0.5:
        print(f"    ✓ Matches!")
    else:
        print(f"    ✗ MISMATCH!")

# Summary of all g_k polynomials
print("\n" + "="*70)
print("SUMMARY OF g_k POLYNOMIALS")
print("="*70)
print(f"g_1(m) = m")
print(f"g_2(m) = m²")
print(f"g_3(m) = m(2m²+1)/3")
print(f"g_4(m) = (5m-4)(2m²-5m+6)/3 = (10m³-33m²+50m-24)/3")

# Evaluate the full CV² asymptotic:
print("\n" + "="*70)
print("FULL ASYMPTOTIC ANALYSIS")
print("="*70)

for n in range(5, 30):
    val = Fraction(0)
    for k in range(1, n//2 + 1):
        mk = n - 2*k
        if mk < 1:
            break
        if k == 1:
            g = g1(mk)
        elif k == 2:
            g = g2(mk)
        elif k == 3:
            g = g3(mk)
        elif k == 4:
            g = g4(mk)
        else:
            break  # Unknown g_k
        ff = falling_factorial(n, 2*k)
        val += Fraction(2*g, ff)

    exact_val = val
    print(f"  n={n}: CV²(approx) = {float(exact_val):.12f}, 2/n = {2.0/n:.12f}, diff = {float(exact_val - Fraction(2,n)):.2e}")

print("\nDone!")
