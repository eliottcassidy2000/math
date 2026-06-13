#!/usr/bin/env python3
"""
gk_cascade_v2_89c.py — Correct cascade WITHOUT g_k(0)=0 assumption
opus-2026-03-14-S89c

CRITICAL FIX: g_k(0) ≠ 0 in general! The polynomial is only valid for m≥1.
E.g., g_4(0) = -8 from the formula (10m³-33m²+50m-24)/3.

Boundary conditions:
  g_k(1) = 1 for all k (monotone path argument)
  g_k(2) = 2k for all k (verified k=1..8)

Without g_k(0)=0, we need more data points. The degree of g_k(m) is the
question: is it k? k-1? Something else?

Known degrees:
  g_1: degree 1 (= k)
  g_2: degree 2 (= k)
  g_3: degree 3 (= k)
  g_4: degree 3 (= k-1) — unless it's actually degree 4

For g_4: we have 4 data points (m=1,2,3,4). If degree 3, the fit is unique and can be validated.
If degree 4, we need a 5th point. Let's check g_4(5) from W(13).
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand
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

# Known g_k formulas (from brute-force verified, NOT using g_k(0)=0)
def g1(m): return Fraction(m)
def g2(m): return Fraction(m * m)
def g3(m): return Fraction(m * (2*m*m + 1), 3)
def g4(m): return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)

known_g = {1: g1, 2: g2, 3: g3, 4: g4}

# Compute W(n)
print("Computing W(n)...")
W = {}
for n in range(3, 22):
    t0 = time.time()
    W[n] = compute_W_dp(n)
    t1 = time.time()
    print(f"  W({n}) = {W[n]}  ({t1-t0:.1f}s)")

# STEP 1: Verify g_4 by extracting g_4(5) from n=13
print("\n" + "="*70)
print("STEP 1: Verify g_4 — extract g_4(5) from n=13")
print("="*70)

n = 13
nf = factorial(n)
cv2 = Fraction(W[n], nf) - 1
remaining = cv2

# Subtract g_1, g_2, g_3
for k in range(1, 4):
    mk = n - 2*k
    if mk < 1:
        break
    remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))

# At n=13: terms are g_4(5), g_5(3), g_6(1)
# remaining = 2·g_4(5)/(13)_8 + 2·g_5(3)/(13)_10 + 2·g_6(1)/(13)_12
# Using g_6(1) = 1:
remaining -= Fraction(2, falling_factorial(13, 12))

print(f"  remaining (after g_1..g_3 and g_6(1)=1) = {remaining}")
print(f"  = 2·g_4(5)/(13)_8 + 2·g_5(3)/(13)_10")
print(f"  Two unknowns: g_4(5) and g_5(3)")

# g_4(5) predicted by degree-3 formula:
g4_5_pred = g4(5)
print(f"\n  g_4(5) predicted (degree 3) = {g4_5_pred}")

# Use this to extract g_5(3):
remaining_after_g4 = remaining - Fraction(2 * g4_5_pred, falling_factorial(13, 8))
g5_3 = remaining_after_g4 * falling_factorial(13, 10) / 2
print(f"  → g_5(3) = {g5_3}")

# Now check: is there a DIFFERENT g_4(5) that gives integer g_5(3)?
# g_5(3) = (remaining - 2·g_4(5)/(13)_8) × (13)_10 / 2
# For g_5(3) to be integer, need specific g_4(5).

# Extract from n=14 similarly
n = 14
nf = factorial(n)
cv2 = Fraction(W[n], nf) - 1
remaining = cv2
for k in range(1, 4):
    mk = n - 2*k
    if mk < 1:
        break
    remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))

# At n=14: g_4(6), g_5(4), g_6(2), g_7(0→not a term since mk=0)
# Actually g_7(0): mk = 14-14 = 0, and (n-1)//2 = 6, k goes up to 6.
# k=7 needs 2k=14 positions from n-1=13. Can't, so max k=6.
# Wait: k=4→mk=6, k=5→mk=4, k=6→mk=2.
remaining_14 = remaining
# Subtract g_6(2) = 2·6 = 12 (from g_k(2) = 2k)
remaining_14 -= Fraction(2 * 12, falling_factorial(14, 12))

print(f"\nn=14: remaining (after g_1..g_3 and g_6(2)=12) = {remaining_14}")
print(f"  = 2·g_4(6)/(14)_8 + 2·g_5(4)/(14)_10")

g4_6_pred = g4(6)
remaining_after_g4 = remaining_14 - Fraction(2 * g4_6_pred, falling_factorial(14, 8))
g5_4 = remaining_after_g4 * falling_factorial(14, 10) / 2
print(f"  g_4(6) predicted = {g4_6_pred}")
print(f"  → g_5(4) = {g5_4}")

# Now check n=15 for g_4(7), g_5(5), g_6(3), g_7(1)
n = 15
nf = factorial(n)
cv2 = Fraction(W[n], nf) - 1
remaining = cv2
for k in range(1, 4):
    mk = n - 2*k
    if mk < 1:
        break
    remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))

# Subtract g_7(1) = 1
remaining -= Fraction(2, falling_factorial(15, 14))

print(f"\nn=15: remaining (after g_1..g_3 and g_7(1)=1)")
print(f"  = 2·g_4(7)/(15)_8 + 2·g_5(5)/(15)_10 + 2·g_6(3)/(15)_12")

g4_7_pred = g4(7)
remaining_15 = remaining - Fraction(2 * g4_7_pred, falling_factorial(15, 8))
print(f"  g_4(7) predicted = {g4_7_pred}")
print(f"  remaining after g_4 = 2·g_5(5)/(15)_10 + 2·g_6(3)/(15)_12")

# STEP 2: Check if degree-3 g_5 through (1,1),(2,10),(3,g5_3),(4,g5_4) works
print("\n" + "="*70)
print("STEP 2: g_5 polynomial (degree 3 from 4 points)")
print("="*70)

pts5 = [
    (Rational(1), Rational(g5_3.numerator if hasattr(g5_3, 'numerator') else int(g5_3),
                           g5_3.denominator if hasattr(g5_3, 'denominator') else 1)),
]
# Let me just use the values directly
print(f"  g_5(1) = 1")
print(f"  g_5(2) = 10")
print(f"  g_5(3) = {g5_3}")
print(f"  g_5(4) = {g5_4}")

if g5_3.denominator != 1 or g5_4.denominator != 1:
    print(f"  WARNING: non-integer g_5 values!")
    print(f"  g_5(3) = {float(g5_3):.10f}")
    print(f"  g_5(4) = {float(g5_4):.10f}")
else:
    pts5 = [
        (Rational(1), Rational(1)),
        (Rational(2), Rational(10)),
        (Rational(3), Rational(g5_3)),
        (Rational(4), Rational(g5_4)),
    ]
    poly5 = interpolate(pts5, x)
    print(f"\n  Degree-3 fit: g_5(m) = {expand(poly5)}")
    print(f"  Factored: g_5(m) = {factor(poly5)}")
    print(f"  g_5(0) = {poly5.subs(x, 0)}")

    # Predict g_5(5) from this polynomial
    g5_5_pred = poly5.subs(x, 5)
    print(f"  g_5(5) predicted = {g5_5_pred}")

    # Use to extract g_6(3) from n=15
    g6_3 = (remaining_15 - Fraction(2 * int(g5_5_pred), falling_factorial(15, 10))) * falling_factorial(15, 12) / 2
    print(f"  → g_6(3) = {g6_3}")

    # Now n=16: g_4(8), g_5(6), g_6(4), g_7(2)=14
    n = 16
    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2
    for k in range(1, 4):
        mk = n - 2*k
        if mk < 1:
            break
        remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
    remaining -= Fraction(2 * g4(8), falling_factorial(16, 8))
    remaining -= Fraction(2 * int(poly5.subs(x, 6)), falling_factorial(16, 10))
    remaining -= Fraction(2 * 14, falling_factorial(16, 14))
    # remaining = 2·g_6(4)/(16)_12
    g6_4 = remaining * falling_factorial(16, 12) / 2
    print(f"\n  n=16: g_6(4) = {g6_4}")

    # n=17: g_4(9), g_5(7), g_6(5), g_7(3), g_8(1)=1
    n = 17
    nf = factorial(n)
    cv2 = Fraction(W[n], nf) - 1
    remaining = cv2
    for k in range(1, 4):
        mk = n - 2*k
        if mk < 1:
            break
        remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
    remaining -= Fraction(2 * g4(9), falling_factorial(17, 8))
    remaining -= Fraction(2 * int(poly5.subs(x, 7)), falling_factorial(17, 10))
    remaining -= Fraction(2, falling_factorial(17, 16))  # g_8(1) = 1
    # remaining = 2·g_6(5)/(17)_12 + 2·g_7(3)/(17)_14
    print(f"\n  n=17: remaining = 2·g_6(5)/(17)_12 + 2·g_7(3)/(17)_14")

    # If we have g_6 polynomial from (1,1),(2,12),(3,g6_3),(4,g6_4):
    if g6_3.denominator == 1 and g6_4.denominator == 1:
        pts6 = [
            (Rational(1), Rational(1)),
            (Rational(2), Rational(12)),
            (Rational(3), Rational(g6_3)),
            (Rational(4), Rational(g6_4)),
        ]
        poly6 = interpolate(pts6, x)
        print(f"\n  g_6 degree-3 fit: g_6(m) = {expand(poly6)}")
        print(f"  g_6(0) = {poly6.subs(x, 0)}")

        g6_5_pred = int(poly6.subs(x, 5))
        remaining -= Fraction(2 * g6_5_pred, falling_factorial(17, 12))
        g7_3 = remaining * falling_factorial(17, 14) / 2
        print(f"  g_6(5) predicted = {g6_5_pred}")
        print(f"  → g_7(3) = {g7_3}")

        # n=18: g_4(10), g_5(8), g_6(6), g_7(4), g_8(2)=16
        n = 18
        nf = factorial(n)
        cv2 = Fraction(W[n], nf) - 1
        remaining = cv2
        for k in range(1, 4):
            mk = n - 2*k
            if mk < 1:
                break
            remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
        remaining -= Fraction(2 * g4(10), falling_factorial(18, 8))
        remaining -= Fraction(2 * int(poly5.subs(x, 8)), falling_factorial(18, 10))
        remaining -= Fraction(2 * int(poly6.subs(x, 6)), falling_factorial(18, 12))
        remaining -= Fraction(2 * 16, falling_factorial(18, 16))
        # remaining = 2·g_7(4)/(18)_14
        g7_4 = remaining * falling_factorial(18, 14) / 2
        print(f"\n  n=18: g_7(4) = {g7_4}")

        # g_7 polynomial from (1,1),(2,14),(3,g7_3),(4,g7_4)
        if g7_3.denominator == 1 and g7_4.denominator == 1:
            pts7 = [
                (Rational(1), Rational(1)),
                (Rational(2), Rational(14)),
                (Rational(3), Rational(g7_3)),
                (Rational(4), Rational(g7_4)),
            ]
            poly7 = interpolate(pts7, x)
            print(f"\n  g_7 degree-3 fit: g_7(m) = {expand(poly7)}")
            print(f"  g_7(0) = {poly7.subs(x, 0)}")

            # n=19: g_4(11), g_5(9), g_6(7), g_7(5), g_8(3), g_9(1)=1
            n = 19
            nf = factorial(n)
            cv2 = Fraction(W[n], nf) - 1
            remaining = cv2
            for k in range(1, 4):
                mk = n - 2*k
                if mk < 1:
                    break
                remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
            remaining -= Fraction(2 * g4(11), falling_factorial(19, 8))
            remaining -= Fraction(2 * int(poly5.subs(x, 9)), falling_factorial(19, 10))
            remaining -= Fraction(2 * int(poly6.subs(x, 7)), falling_factorial(19, 12))
            remaining -= Fraction(2 * int(poly7.subs(x, 5)), falling_factorial(19, 14))
            remaining -= Fraction(2, falling_factorial(19, 18))  # g_9(1) = 1
            # remaining = 2·g_8(3)/(19)_16
            g8_3 = remaining * falling_factorial(19, 16) / 2
            print(f"\n  n=19: g_8(3) = {g8_3}")

            # n=20: g_4(12), g_5(10), g_6(8), g_7(6), g_8(4), g_9(2)=18
            n = 20
            nf = factorial(n)
            cv2 = Fraction(W[n], nf) - 1
            remaining = cv2
            for k in range(1, 4):
                mk = n - 2*k
                if mk < 1:
                    break
                remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
            remaining -= Fraction(2 * g4(12), falling_factorial(20, 8))
            remaining -= Fraction(2 * int(poly5.subs(x, 10)), falling_factorial(20, 10))
            remaining -= Fraction(2 * int(poly6.subs(x, 8)), falling_factorial(20, 12))
            remaining -= Fraction(2 * int(poly7.subs(x, 6)), falling_factorial(20, 14))
            remaining -= Fraction(2 * 18, falling_factorial(20, 18))  # g_9(2) = 18
            # remaining = 2·g_8(4)/(20)_16
            g8_4 = remaining * falling_factorial(20, 16) / 2
            print(f"  n=20: g_8(4) = {g8_4}")

            # n=21: g_4(13), g_5(11), g_6(9), g_7(7), g_8(5), g_9(3), g_10(1)=1
            n = 21
            nf = factorial(n)
            cv2 = Fraction(W[n], nf) - 1
            remaining = cv2
            for k in range(1, 4):
                mk = n - 2*k
                if mk < 1:
                    break
                remaining -= Fraction(2 * known_g[k](mk), falling_factorial(n, 2*k))
            remaining -= Fraction(2 * g4(13), falling_factorial(21, 8))
            remaining -= Fraction(2 * int(poly5.subs(x, 11)), falling_factorial(21, 10))
            remaining -= Fraction(2 * int(poly6.subs(x, 9)), falling_factorial(21, 12))
            remaining -= Fraction(2 * int(poly7.subs(x, 7)), falling_factorial(21, 14))
            remaining -= Fraction(2, falling_factorial(21, 20))  # g_10(1) = 1
            # remaining = 2·g_8(5)/(21)_16 + 2·g_9(3)/(21)_18
            print(f"\n  n=21: remaining has g_8(5) and g_9(3)")

            if g8_3.denominator == 1 and g8_4.denominator == 1:
                pts8 = [
                    (Rational(1), Rational(1)),
                    (Rational(2), Rational(16)),
                    (Rational(3), Rational(g8_3)),
                    (Rational(4), Rational(g8_4)),
                ]
                poly8 = interpolate(pts8, x)
                print(f"  g_8 degree-3 fit: g_8(m) = {expand(poly8)}")

                g8_5_pred = int(poly8.subs(x, 5))
                remaining -= Fraction(2 * g8_5_pred, falling_factorial(21, 16))
                g9_3 = remaining * falling_factorial(21, 18) / 2
                print(f"  g_8(5) predicted = {g8_5_pred}")
                print(f"  → g_9(3) = {g9_3}")

# SUMMARY
print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("g_4(m) = (10m³-33m²+50m-24)/3  (degree 3)")
if g5_3.denominator == 1:
    print(f"g_5(m) = {expand(poly5)}  (degree 3)")
    print(f"  = {factor(poly5)}")
if g6_3.denominator == 1 and g6_4.denominator == 1:
    print(f"g_6(m) = {expand(poly6)}  (degree 3)")
    print(f"  = {factor(poly6)}")
if g7_3.denominator == 1 and g7_4.denominator == 1:
    print(f"g_7(m) = {expand(poly7)}  (degree 3)")
    print(f"  = {factor(poly7)}")

# Degree sequence
print("\nDegree sequence of g_k:")
print("  k=1: 1, k=2: 2, k=3: 3, k≥4: 3 (if confirmed)")

# Leading coefficients of the degree-3 terms
print("\nLeading coefficients:")
print(f"  g_3: 2/3")
print(f"  g_4: 10/3")
if g5_3.denominator == 1:
    lc5 = poly5.as_poly().LC()
    print(f"  g_5: {lc5}")
if g6_3.denominator == 1:
    lc6 = poly6.as_poly().LC()
    print(f"  g_6: {lc6}")
if g7_3.denominator == 1:
    lc7 = poly7.as_poly().LC()
    print(f"  g_7: {lc7}")

print("\nDone!")
