#!/usr/bin/env python3
"""
gk_coefficients_v2_89c.py — Extract all g_k polynomials via manual cascade
opus-2026-03-15-S89c

Strategy: for each k, use known polynomials for k'<k to extract 4 data points
g_k(1)=1, g_k(2)=2k, g_k(3), g_k(4), then fit degree-3 polynomial.
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

# Compute W
print("Computing W(n)...")
W = {}
for n in range(3, 24):
    t0 = time.time()
    W[n] = compute_W_dp(n)
    t1 = time.time()
    print(f"  W({n}) = {W[n]}  ({t1-t0:.1f}s)")

# Known g_k formulas
g_formulas = {}  # k -> lambda m: Fraction

g_formulas[1] = lambda m: Fraction(m)
g_formulas[2] = lambda m: Fraction(m * m)
g_formulas[3] = lambda m: Fraction(m * (2*m*m + 1), 3)
g_formulas[4] = lambda m: Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)

def g_eval(k, m):
    if k in g_formulas:
        return g_formulas[k](m)
    # Boundary conditions for all k
    if m == 1:
        return Fraction(1)
    if m == 2:
        return Fraction(2 * k)
    return None

def extract_gk_at_m(target_k, target_m):
    """Extract g_{target_k}(target_m) from the appropriate n = target_m + 2*target_k."""
    n = target_m + 2 * target_k
    if n not in W:
        return None

    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2

    # Subtract all known contributions
    for k in range(1, (n-1)//2 + 1):
        mk = n - 2*k
        if mk < 1:
            break
        if k == target_k and mk == target_m:
            continue  # This is what we're solving for
        val = g_eval(k, mk)
        if val is not None:
            remaining -= Fraction(2 * val, falling_factorial(n, 2*k))
        else:
            # Unknown — can't fully extract
            return None

    # remaining = 2·g_{target_k}(target_m) / (n)_{2·target_k}
    result = remaining * falling_factorial(n, 2*target_k) / 2
    return result

# Sequential extraction: for each k, extract g_k(3) then g_k(4), fit polynomial
print("\n" + "="*70)
print("SEQUENTIAL EXTRACTION")
print("="*70)

g_polys = {}
coefficients = {}

# Known polynomials
for k in [3, 4]:
    if k == 3:
        poly = Rational(2,3)*x**3 + Rational(1,3)*x
    else:
        poly = (Rational(10)*x**3 - 33*x**2 + 50*x - 24) / 3
    g_polys[k] = expand(poly)
    p = Poly(3*expand(poly), x)
    c = p.all_coeffs()
    while len(c) < 4:
        c = [0] + c
    coefficients[k] = tuple(int(ci) for ci in c)

for target_k in range(5, 12):
    print(f"\n--- Extracting g_{target_k} ---")

    # We know g_{target_k}(1) = 1 and g_{target_k}(2) = 2·target_k
    data = {1: Fraction(1), 2: Fraction(2 * target_k)}

    # Try to extract g_{target_k}(3) from n = 3 + 2·target_k
    g3_val = extract_gk_at_m(target_k, 3)
    if g3_val is not None:
        if g3_val.denominator == 1:
            data[3] = g3_val
            print(f"  g_{target_k}(3) = {g3_val} ✓")
        else:
            print(f"  g_{target_k}(3) = {g3_val} = {float(g3_val):.4f} NOT INTEGER!")
            break

    # Try to extract g_{target_k}(4) from n = 4 + 2·target_k
    g4_val = extract_gk_at_m(target_k, 4)
    if g4_val is not None:
        if g4_val.denominator == 1:
            data[4] = g4_val
            print(f"  g_{target_k}(4) = {g4_val} ✓")
        else:
            print(f"  g_{target_k}(4) = {g4_val} = {float(g4_val):.4f} NOT INTEGER!")
            break

    if len(data) >= 4:
        # Fit degree-3 polynomial through (1,1), (2,2k), (3,g3), (4,g4)
        pts = [(Rational(m), Rational(data[m])) for m in [1, 2, 3, 4]]
        poly = interpolate(pts, x)
        poly = expand(poly)
        g_polys[target_k] = poly

        # Create formula
        def make_formula(p):
            def f(m):
                return Fraction(p.subs(x, m))
            return f
        g_formulas[target_k] = make_formula(poly)

        # Extract coefficients (multiply by 3 to clear denominators)
        p3 = Poly(3*poly, x)
        c = p3.all_coeffs()
        while len(c) < 4:
            c = [0] + c
        coefficients[target_k] = tuple(int(ci) for ci in c)

        print(f"  g_{target_k}(m) = {poly}")
        print(f"  = {factor(poly)}")
        print(f"  3·g_{target_k}(m) = {coefficients[target_k][0]}m³ + ({coefficients[target_k][1]})m² + ({coefficients[target_k][2]})m + ({coefficients[target_k][3]})")

        # Verify by extracting g_{target_k}(5) if possible
        g5_val = extract_gk_at_m(target_k, 5)
        if g5_val is not None:
            pred = Fraction(poly.subs(x, 5))
            match = "✓" if g5_val == pred else "✗ MISMATCH!"
            print(f"  Verification: g_{target_k}(5) = {g5_val}, predicted {pred} {match}")
            if g5_val != pred:
                print(f"  *** DEGREE 3 IS WRONG FOR k={target_k}! ***")
                break
    else:
        print(f"  Insufficient data for g_{target_k}")
        break

# Print coefficient table
print("\n" + "="*70)
print("COEFFICIENT TABLE: 3·g_k(m) = a_k·m³ + b_k·m² + c_k·m + d_k")
print("="*70)

print(f"\n{'k':>3} {'a_k':>20} {'b_k':>20} {'c_k':>20} {'d_k':>20}")
print("-"*85)
for k in sorted(coefficients):
    a, b, c, d = coefficients[k]
    print(f"{k:>3} {a:>20} {b:>20} {c:>20} {d:>20}")

# Leading coefficient sequence
print("\n" + "="*70)
print("LEADING COEFFICIENT a_k ANALYSIS")
print("="*70)
a_seq = [coefficients[k][0] for k in sorted(coefficients)]
ks = sorted(coefficients)
print(f"a_k = {a_seq}")

# Ratios
print("\nRatios a_k / a_{k-1}:")
for i in range(1, len(a_seq)):
    r = a_seq[i] / a_seq[i-1]
    print(f"  a_{ks[i]}/a_{ks[i]-1} = {r:.6f}")

# Ratios of ratios
print("\nSecond ratios:")
ratios = [a_seq[i]/a_seq[i-1] for i in range(1, len(a_seq))]
for i in range(1, len(ratios)):
    print(f"  r_{ks[i+1]}/r_{ks[i]} = {ratios[i]/ratios[i-1]:.6f}")

# a_k / (2k)!
print("\na_k / (2k)!:")
for k in ks:
    a = coefficients[k][0]
    print(f"  k={k}: {a} / {factorial(2*k)} = {a/factorial(2*k):.15f}")

# a_k / (2k-2)!
print("\na_k / (2k-2)!:")
for k in ks:
    a = coefficients[k][0]
    print(f"  k={k}: {a} / {factorial(2*k-2)} = {a/factorial(2*k-2):.10f}")

# Factor the leading coefficients
print("\nFactorizations of a_k:")
for k in ks:
    a = coefficients[k][0]
    from sympy import factorint
    print(f"  a_{k} = {a} = {dict(factorint(abs(a)))}")

# Constant term d_k analysis
print("\n" + "="*70)
print("CONSTANT TERM d_k ANALYSIS")
print("="*70)
d_seq = [coefficients[k][3] for k in ks]
print(f"d_k = {d_seq}")
print("\n-d_k factorizations:")
for k in ks:
    d = -coefficients[k][3]
    if d > 0:
        from sympy import factorint
        print(f"  -d_{k} = {d} = {dict(factorint(d))}")
    else:
        print(f"  -d_{k} = {d}")

# g_k(3) sequence
print("\n" + "="*70)
print("g_k(3) SEQUENCE")
print("="*70)
g3_seq = []
for k in ks:
    a, b, c, d = coefficients[k]
    g3 = (27*a + 9*b + 3*c + d) // 3
    g3_seq.append(g3)
    print(f"  g_{k}(3) = {g3}")

# g_k(m) for small m
print("\n" + "="*70)
print("g_k(m) TABLE")
print("="*70)
print(f"{'k':>3}", end="")
for m in range(1, 8):
    print(f" {'m='+str(m):>15}", end="")
print()
for k in ks:
    a, bc, c, d = coefficients[k]
    print(f"{k:>3}", end="")
    for m in range(1, 8):
        val = (a*m**3 + bc*m**2 + c*m + d) // 3
        print(f" {val:>15}", end="")
    print()

print("\nDone!")
