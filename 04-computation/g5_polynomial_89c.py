#!/usr/bin/env python3
"""
g5_polynomial_89c.py — Find g_5 polynomial using boundary condition g_5(0)=0
opus-2026-03-14-S89c

Data points:
  g_5(0) = 0   (boundary: no room for dominos)
  g_5(1) = 1   (monotone path)
  g_5(2) = 10  (g_k(2) = 2k pattern)
  g_5(3) = 211 (exact from cascading extraction)
  g_5(4) = 1380 (exact, using g_6(2) = 12)

5 points → degree 4 polynomial.
"""

from sympy import symbols, interpolate, Rational, factor, expand, simplify

x = symbols('x')

# Data points
pts = [
    (Rational(0), Rational(0)),
    (Rational(1), Rational(1)),
    (Rational(2), Rational(10)),
    (Rational(3), Rational(211)),
    (Rational(4), Rational(1380)),
]

poly = interpolate(pts, x)
print("g_5 polynomial through 5 points (degree 4):")
print(f"  g_5(m) = {expand(poly)}")
print(f"  g_5(m) = {factor(poly)}")
print()

# Check values
for m, expected in pts:
    val = poly.subs(x, m)
    print(f"  g_5({m}) = {val} (expected {expected}) {'✓' if val == expected else '✗'}")

# Predict g_5(5), g_5(6), g_5(7)
print()
for m in range(5, 11):
    val = poly.subs(x, m)
    print(f"  g_5({m}) = {val}")

# Now check: does g_5 have degree 5?
# For degree 5, we need one more data point.
# We can use the n=15 extraction to get g_5(5) + g_6(3) constraint,
# then check consistency.
print("\n" + "="*70)
print("CHECKING IF g_5 IS DEGREE 5 (needs g_6 data)")
print("="*70)

# From gk_full_cascade output at n=15:
# remaining after g_1..g_4 and g_7(1)=1 subtraction:
# 2·g_5(5)/(15)_10 + 2·g_6(3)/(15)_12 = remaining
# If g_5 is degree 4, g_5(5) is predicted. Then we can extract g_6(3).

from fractions import Fraction
from math import factorial
from functools import reduce

def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)

def g1(m): return Fraction(m)
def g2(m): return Fraction(m * m)
def g3(m): return Fraction(m * (2*m*m + 1), 3)
def g4(m): return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)
known_g = {1: g1, 2: g2, 3: g3, 4: g4}

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

# Precompute W
W = {}
for n in range(3, 21):
    W[n] = compute_W_dp(n)

# g_5 as fraction polynomial
g5_coeffs = expand(poly)  # degree 4 polynomial in x
print(f"\nDegree-4 g_5(m) = {g5_coeffs}")

def g5_frac(m):
    """Evaluate the degree-4 g_5 polynomial."""
    return Fraction(int(poly.subs(x, m)))

# Try extracting g_6 values using degree-4 g_5
g6_vals = {0: Fraction(0), 1: Fraction(1), 2: Fraction(12)}

print("\nExtracting g_6 values using degree-4 g_5:")
for n in range(15, 21):
    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2

    # Subtract g_1..g_4
    for k in range(1, 5):
        mk = n - 2*k
        if mk < 1:
            break
        remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))

    # Subtract g_5 (degree-4 assumption)
    mk5 = n - 10
    if mk5 >= 1:
        remaining -= Fraction(2 * g5_frac(mk5), falling_factorial(n, 10))

    # Subtract known g_k for k≥7
    for k in range(7, (n-1)//2 + 1):
        mk = n - 2*k
        if mk == 1:
            remaining -= Fraction(2, falling_factorial(n, 2*k))
        elif mk == 2:
            remaining -= Fraction(2 * 2*k, falling_factorial(n, 2*k))
        elif mk == 0:
            pass  # g_k(0) = 0
        else:
            break  # unknown

    # Now remaining should be 2·g_6(n-12)/(n)_12 + higher unknown g_6 terms
    mk6 = n - 12
    if mk6 >= 1:
        # Check if there are other unknown terms
        other_unknowns = []
        for k in range(7, (n-1)//2 + 1):
            mk = n - 2*k
            if mk > 2:
                other_unknowns.append((k, mk))

        if not other_unknowns:
            g6_val = remaining * falling_factorial(n, 12) / 2
            is_int = g6_val.denominator == 1
            print(f"  n={n}: g_6({mk6}) = {g6_val}{'' if is_int else f' = {float(g6_val):.6f} NOT INT'}")
            if is_int:
                g6_vals[mk6] = g6_val
        else:
            print(f"  n={n}: g_6({mk6}) has unresolved higher terms: {other_unknowns}")
            # Approximate
            g6_approx = remaining * falling_factorial(n, 12) / 2
            print(f"    g_6({mk6}) ≈ {float(g6_approx):.6f}")

print("\n" + "="*70)
print("g_6 VALUES")
print("="*70)
for mk in sorted(g6_vals):
    print(f"  g_6({mk}) = {g6_vals[mk]}")

# Fit g_6 polynomial
pts6 = [(Rational(mk), Rational(g6_vals[mk])) for mk in sorted(g6_vals) if mk >= 0]
if len(pts6) >= 4:
    poly6 = interpolate(pts6, x)
    print(f"\n  g_6(m) = {expand(poly6)}")
    print(f"  g_6(m) = {factor(poly6)}")

# Now check: is degree-4 g_5 correct? If so, g_6 values should be integers.
# If NOT integers, g_5 is higher degree and our values are wrong.

print("\n" + "="*70)
print("CONSISTENCY CHECK: Does degree-4 g_5 produce integer g_6?")
print("="*70)

all_int = True
for mk in sorted(g6_vals):
    v = g6_vals[mk]
    if v.denominator != 1:
        all_int = False
        print(f"  g_6({mk}) = {v} = {float(v):.6f} — NOT INTEGER!")
    else:
        print(f"  g_6({mk}) = {v} — integer ✓")

if all_int:
    print("\n✓ Degree-4 g_5 is CONSISTENT (all g_6 values are integers)")
else:
    print("\n✗ Degree-4 g_5 produces non-integer g_6 — g_5 may be higher degree")

# Try degree 5 for g_5 (if we have enough data)
# degree 5 needs 6 points. We have (0,0),(1,1),(2,10),(3,211),(4,1380)
# and potentially g_5(5) from n=15 if we can get g_6(3)
print("\n" + "="*70)
print("ATTEMPTING DEGREE 5 g_5")
print("="*70)

# With degree 4 assumption, we got g_6(3). Now use THAT to get g_5(5) directly from n=15.
# At n=15: remaining_after_g1..g4 = 2·g_5(5)/(15)_10 + 2·g_6(3)/(15)_12 + 2·g_7(1)/(15)_14
# If g_6(3) from above is integer, great.
if 3 in g6_vals and g6_vals[3].denominator == 1:
    print(f"Using g_6(3) = {g6_vals[3]} to extract g_5(5) independently")

    n = 15
    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2

    for k in range(1, 5):
        mk = n - 2*k
        if mk < 1:
            break
        remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))

    # Subtract g_6(3)
    remaining -= Fraction(2 * g6_vals[3], falling_factorial(15, 12))
    # Subtract g_7(1) = 1
    remaining -= Fraction(2, falling_factorial(15, 14))

    g5_5 = remaining * falling_factorial(15, 10) / 2
    print(f"g_5(5) = {g5_5}")

    if g5_5.denominator == 1:
        # Now fit degree 5 through 6 points
        pts5_ext = pts + [(Rational(5), Rational(g5_5))]
        poly5 = interpolate(pts5_ext, x)
        print(f"\nDegree-5 g_5(m) = {expand(poly5)}")
        print(f"Degree-5 g_5(m) = {factor(poly5)}")

        # Compare predictions
        print("\nDegree-4 vs Degree-5 predictions:")
        for m in range(6, 11):
            v4 = poly.subs(x, m)
            v5 = poly5.subs(x, m)
            print(f"  g_5({m}): deg4={v4}, deg5={v5}, diff={v5-v4}")

        # Check g_5(0) = 0 for degree 5
        print(f"\n  g_5(0) = {poly5.subs(x, 0)} (should be 0)")

print("\nDone!")
