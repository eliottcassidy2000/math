#!/usr/bin/env python3
"""
gk_coefficient_analysis_89c.py — Analyze coefficients of degree-3 g_k polynomials
opus-2026-03-14-S89c

All g_k(m) for k≥3 are degree 3: g_k(m) = (a_k m³ + b_k m² + c_k m + d_k)/3

Extract a_k, b_k, c_k, d_k and look for patterns/recurrences.
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand, Poly
import time

x = symbols('x')

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

# Known polynomial formulas
def g_formula(k, m):
    """Evaluate g_k(m) for k=1..4."""
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k == 3: return Fraction(m * (2*m*m + 1), 3)
    if k == 4: return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)
    return None

# Compute W(n)
print("Computing W(n)...")
W = {}
for n in range(3, 24):
    t0 = time.time()
    W[n] = compute_W_dp(n)
    t1 = time.time()
    print(f"  W({n}) = {W[n]}  ({t1-t0:.1f}s)")

# Cascading extraction
print("\n" + "="*70)
print("CASCADING EXTRACTION")
print("="*70)

# Store g_k polynomials as they're found
g_polys = {}  # k -> sympy polynomial
g_data = {}   # k -> {m: value}

# Initialize with known polynomials
for k in range(1, 5):
    g_data[k] = {}

# Boundary conditions for k ≥ 5
for k in range(5, 12):
    g_data[k] = {1: Fraction(1), 2: Fraction(2*k)}

def g_eval(k, m):
    """Evaluate g_k(m) using best available info."""
    if k <= 4:
        return g_formula(k, m)
    if k in g_polys:
        return Fraction(int(g_polys[k].subs(x, m)))
    if k in g_data and m in g_data[k]:
        return g_data[k][m]
    return None

def try_fit_polynomial(k):
    """Try to fit degree-3 polynomial to g_k data."""
    data = g_data[k]
    int_pts = [(m, data[m]) for m in sorted(data) if data[m].denominator == 1 and m >= 1]
    if len(int_pts) < 4:
        return None
    pts = [(Rational(m), Rational(v)) for m, v in int_pts[:4]]
    poly = interpolate(pts, x)
    # Verify remaining points
    for m, v in int_pts[4:]:
        if poly.subs(x, m) != Rational(v):
            return None  # Doesn't fit!
    return poly

# Iterative cascade
for iteration in range(8):
    max_known = max(k for k in range(1, 12) if k <= 4 or k in g_polys)
    print(f"\nIteration {iteration+1}: polynomials known for k=1..{max_known}")

    new_data = False
    for n in sorted(W):
        nf = factorial(n)
        cv2 = Fraction(W[n], nf) - 1
        remaining = cv2

        # Subtract all known contributions
        for k in range(1, (n-1)//2 + 1):
            mk = n - 2*k
            if mk < 1:
                break
            val = g_eval(k, mk)
            if val is not None:
                remaining -= Fraction(2 * val, falling_factorial(n, 2*k))
            else:
                # This k is the first unknown — can we solve?
                # Check how many unknowns remain
                unknowns = []
                for k2 in range(k, (n-1)//2 + 1):
                    mk2 = n - 2*k2
                    if mk2 < 1:
                        break
                    if g_eval(k2, mk2) is None:
                        unknowns.append((k2, mk2))

                if len(unknowns) == 1:
                    ku, mu = unknowns[0]
                    val = remaining * falling_factorial(n, 2*ku) / 2
                    if val.denominator == 1:
                        g_data[ku][mu] = val
                        new_data = True
                break

    # Try to fit polynomials
    for k in range(5, 12):
        if k not in g_polys:
            poly = try_fit_polynomial(k)
            if poly is not None:
                g_polys[k] = poly
                print(f"  Found polynomial for g_{k}: {expand(poly)}")
                new_data = True

    if not new_data:
        print("  No new data extracted, stopping.")
        break

# Print all polynomials
print("\n" + "="*70)
print("ALL g_k POLYNOMIALS (degree 3 for k≥3)")
print("="*70)

# Extract coefficients
coefficients = {}
for k in range(3, 12):
    if k <= 4:
        if k == 3:
            poly = Rational(2,3)*x**3 + Rational(1,3)*x
        else:
            poly = Rational(10,3)*x**3 - 11*x**2 + Rational(50,3)*x - 8
    elif k in g_polys:
        poly = g_polys[k]
    else:
        continue

    p = Poly(3*expand(poly), x)
    coeffs = p.all_coeffs()
    while len(coeffs) < 4:
        coeffs = [0] + coeffs
    a, b, c, d = [int(c) for c in coeffs]
    coefficients[k] = (a, b, c, d)
    print(f"  g_{k}(m) = ({a}m³ + ({b})m² + ({c})m + ({d})) / 3")
    print(f"        = {factor(poly)}")

# Analyze coefficient sequences
print("\n" + "="*70)
print("COEFFICIENT SEQUENCES")
print("="*70)

ks = sorted(coefficients)
print(f"\n{'k':>3} {'a_k':>15} {'b_k':>15} {'c_k':>15} {'d_k':>15}")
print("-"*65)
for k in ks:
    a, b, c, d = coefficients[k]
    print(f"{k:>3} {a:>15} {b:>15} {c:>15} {d:>15}")

# Check: a_k + b_k + c_k + d_k = 3 (since g_k(1) = 1 → 3g_k(1) = 3)
print("\nVerification: a_k + b_k + c_k + d_k = 3 (from g_k(1) = 1):")
for k in ks:
    a, b, c, d = coefficients[k]
    s = a + b + c + d
    print(f"  k={k}: {a} + {b} + {c} + {d} = {s} {'✓' if s == 3 else '✗'}")

# Check: 8a + 4b + 2c + d = 6k (from g_k(2) = 2k → 3g_k(2) = 6k)
print("\nVerification: 8a + 4b + 2c + d = 6k (from g_k(2) = 2k):")
for k in ks:
    a, b, c, d = coefficients[k]
    s = 8*a + 4*b + 2*c + d
    print(f"  k={k}: {s}, expected {6*k} {'✓' if s == 6*k else '✗'}")

# Ratios of leading coefficients
print("\nLeading coefficient a_k ratios:")
prev = None
for k in ks:
    a = coefficients[k][0]
    if prev:
        ratio = a / prev
        print(f"  a_{k}/a_{k-1} = {a}/{prev} = {ratio:.6f}")
    prev = a

# Try to find recurrence for a_k
print("\n" + "="*70)
print("RECURRENCE SEARCH FOR a_k")
print("="*70)
a_vals = [coefficients[k][0] for k in ks]
print(f"a_k = {a_vals}")

# Check a_k = c1 * k * a_{k-1} + c2 * a_{k-1}
# a_4 = c1*4*a_3 + c2*a_3 = a_3(4c1+c2) → 10 = 2(4c1+c2)
# a_5 = a_4(5c1+c2) → 388 = 10(5c1+c2)
# From first: 4c1+c2 = 5
# From second: 5c1+c2 = 38.8
# → c1 = 33.8, c2 = 5 - 4*33.8 = -130.2 → not integers
print("Linear recurrence a_k = (αk+β)·a_{k-1}:")
if len(a_vals) >= 3:
    # r_k = a_k/a_{k-1} should be linear in k
    for i in range(1, len(a_vals)):
        r = a_vals[i] / a_vals[i-1]
        print(f"  r_{ks[i]} = a_{ks[i]}/a_{ks[i]-1} = {r:.6f}")

# Try a_k = α·(2k)!·f(k) or a_k = α·k!·f(k)
print("\nNormalized by factorials:")
for k in ks:
    a = coefficients[k][0]
    fk = factorial(k)
    f2k = factorial(2*k)
    print(f"  k={k}: a_k/k! = {a/fk:.6f}, a_k/(2k)! = {a/f2k:.10f}, a_k/(2k-2)! = {a/factorial(2*k-2):.6f}")

# Check if d_k = -3·g_k(0) relates to something
print("\n" + "="*70)
print("g_k(0) VALUES (= d_k/3)")
print("="*70)
for k in ks:
    d = coefficients[k][3]
    print(f"  g_{k}(0) = {d}/3 = {Fraction(d,3)}")

# The -d_k sequence: 0, 24, 1776, 342960, 100014720, ...
print("\n-d_k sequence:")
for k in ks:
    print(f"  k={k}: -d_k = {-coefficients[k][3]}")

# Search for relation: d_k = -a_k + b_k - c_k + d_k... no.
# d_k from a_k + b_k + c_k + d_k = 3: d_k = 3 - a_k - b_k - c_k

# Sum relations
print("\n" + "="*70)
print("g_k(3) SEQUENCE")
print("="*70)
for k in ks:
    a, b, c, d = coefficients[k]
    g3 = (27*a + 9*b + 3*c + d) / 3
    print(f"  g_{k}(3) = {int(g3)}")

print("\nDone!")
