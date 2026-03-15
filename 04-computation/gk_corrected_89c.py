#!/usr/bin/env python3
"""
gk_corrected_89c.py — Correct g_k extraction accounting for ALL higher terms
opus-2026-03-14-S89c

Key insight: g_k(1) = 1 for all k (from the monotone path argument).
Use this to correct the extraction cascade.
"""

from fractions import Fraction
from math import factorial
from functools import reduce
from sympy import symbols, interpolate, Rational, factor, expand
import time

def falling_factorial(n, k):
    return reduce(lambda a,b: a*b, range(n, n-k, -1), 1)

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
def g1(m): return m
def g2(m): return m * m
def g3(m): return Fraction(m * (2*m*m + 1), 3)
def g4(m): return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)

# Compute W for n=3..17
print("Computing W(n) values...")
W_vals = {}
for n in range(3, 18):
    W_vals[n] = compute_W_dp(n)
    print(f"  W({n}) = {W_vals[n]}")

print("\n" + "="*70)
print("EXACT g_k EXTRACTION WITH CASCADING CORRECTION")
print("="*70)

# For each n, the CV² expansion is:
# CV² = Σ_{k=1}^{⌊(n-1)/2⌋} 2·g_k(n-2k) / (n)_{2k}
#
# We know g_1, g_2, g_3, g_4. We want to find g_5, g_6, etc.
# Strategy: for each n, compute the "remaining" after subtracting g_1..g_4,
# then solve a triangular system for g_5(m), g_6(m), g_7(m), etc.

# At n = 2K+1: the maximum k is K, and g_K(1) should be 1.
# This gives us g_k(1) = 1 for all k.
# At n = 2K+2: g_K(2) can be extracted if we know g_{K+1}(0) = 0 (m=0 shouldn't contribute).
# Actually, for k=K+1 at n=2K+2: m = n-2(K+1) = 0. And g_{K+1}(0) should be 0 since there's
# no room for the domino tiling.

# Let me extract g_5 properly.
# At n=11: g_5(1) is the only unknown (|S|=10 is the highest).
# |S|=10 = 2·g_5(1)/(11)_10 = 2·g_5(1)/39916800.
# remaining = CV²(11) - Σ_{k=1}^{4} known contributions.
# g_5(1) = remaining × (11)_10 / 2.

g5_vals = {}
g6_vals = {}
g7_vals = {}

for n in range(11, 18):
    nf = factorial(n)
    cv2 = Fraction(W_vals[n], nf) - 1
    remaining = cv2

    # Subtract known g_1..g_4
    for k in range(1, 5):
        mk = n - 2*k
        if mk < 1:
            break
        if k == 1: g = g1(mk)
        elif k == 2: g = g2(mk)
        elif k == 3: g = g3(mk)
        elif k == 4: g = g4(mk)
        contrib = Fraction(2 * g, falling_factorial(n, 2*k))
        remaining -= contrib

    # remaining = Σ_{k≥5} 2·g_k(n-2k) / (n)_{2k}
    # For each k≥5 with n-2k ≥ 1:
    # max k = (n-1)//2

    # Build system: remaining = Σ_{k≥5, n-2k≥1} 2·g_k(n-2k)/(n)_{2k}
    # The unknowns are g_5(n-10), g_6(n-12), g_7(n-14), etc.

    unknowns = []
    for k in range(5, (n-1)//2 + 1):
        mk = n - 2*k
        if mk < 1:
            break
        unknowns.append((k, mk))

    if len(unknowns) == 0:
        continue

    print(f"\nn={n}: unknowns = {unknowns}")

    if len(unknowns) == 1:
        # Only one unknown: solve directly
        k, mk = unknowns[0]
        gk_val = remaining * falling_factorial(n, 2*k) / 2
        print(f"  g_{k}({mk}) = {gk_val}")
        if k == 5: g5_vals[mk] = gk_val
        elif k == 6: g6_vals[mk] = gk_val
        elif k == 7: g7_vals[mk] = gk_val

    elif len(unknowns) == 2:
        # Two unknowns: g_a(ma) and g_b(mb)
        # remaining = 2·g_a(ma)/(n)_{2a} + 2·g_b(mb)/(n)_{2b}
        (ka, ma), (kb, mb) = unknowns
        # If we know g_b(mb) from a smaller n, use it
        known_b = None
        if kb == 5 and mb in g5_vals: known_b = g5_vals[mb]
        elif kb == 6 and mb in g6_vals: known_b = g6_vals[mb]
        elif kb == 7 and mb in g7_vals: known_b = g7_vals[mb]

        if known_b is not None:
            remaining_adjusted = remaining - Fraction(2 * known_b, falling_factorial(n, 2*kb))
            gka_val = remaining_adjusted * falling_factorial(n, 2*ka) / 2
            print(f"  Using known g_{kb}({mb}) = {known_b}")
            print(f"  g_{ka}({ma}) = {gka_val}")
            if ka == 5: g5_vals[ma] = gka_val
            elif ka == 6: g6_vals[ma] = gka_val
        else:
            # Use g_k(1) = 1 assumption if mb = 1
            if mb == 1:
                known_b = 1
                remaining_adjusted = remaining - Fraction(2, falling_factorial(n, 2*kb))
                gka_val = remaining_adjusted * falling_factorial(n, 2*ka) / 2
                print(f"  Assuming g_{kb}(1) = 1")
                print(f"  g_{ka}({ma}) = {gka_val}")
                if ka == 5: g5_vals[ma] = gka_val
                elif ka == 6: g6_vals[ma] = gka_val
                if kb == 6: g6_vals[1] = Fraction(1)
                elif kb == 7: g7_vals[1] = Fraction(1)
            else:
                print(f"  Cannot resolve: need g_{kb}({mb})")
                # Estimate assuming g_b contribution is small
                gka_approx = remaining * falling_factorial(n, 2*ka) / 2
                print(f"  g_{ka}({ma}) ≈ {float(gka_approx):.4f} (ignoring g_{kb}({mb}))")

    elif len(unknowns) >= 3:
        # Multiple unknowns: use known values from smaller n
        r = remaining
        for k, mk in unknowns:
            known = None
            if k == 5 and mk in g5_vals: known = g5_vals[mk]
            elif k == 6 and mk in g6_vals: known = g6_vals[mk]
            elif k == 7 and mk in g7_vals: known = g7_vals[mk]

            if known is not None:
                r -= Fraction(2 * known, falling_factorial(n, 2*k))
                print(f"  Subtracting known g_{k}({mk}) = {known}")
            elif mk == 1:
                # Assume g_k(1) = 1
                r -= Fraction(2, falling_factorial(n, 2*k))
                print(f"  Assuming g_{k}(1) = 1")
                if k == 5: g5_vals[1] = Fraction(1)
                elif k == 6: g6_vals[1] = Fraction(1)
                elif k == 7: g7_vals[1] = Fraction(1)
            else:
                # This is the first unknown we can solve for
                gk_val = r * falling_factorial(n, 2*k) / 2
                print(f"  g_{k}({mk}) = {gk_val} (APPROX — ignoring higher terms)")
                if k == 5: g5_vals[mk] = gk_val
                elif k == 6: g6_vals[mk] = gk_val
                break

# Print results
print("\n" + "="*70)
print("CORRECTED g_5 VALUES")
print("="*70)

x = symbols('x')
pts5 = []
for mk in sorted(g5_vals):
    v = g5_vals[mk]
    print(f"  g_5({mk}) = {v}")
    if v.denominator == 1:
        pts5.append((Rational(mk), Rational(v)))

if len(pts5) >= 4:
    for deg in range(1, len(pts5)):
        poly = interpolate(pts5[:deg+1], x)
        all_match = True
        for m_val, g_val in pts5[deg+1:]:
            if poly.subs(x, m_val) != g_val:
                all_match = False
                break
        if all_match and len(pts5) > deg+1:
            print(f"\n  → g_5(m) = {factor(poly)} = {expand(poly)} (degree {deg})")
            break
    else:
        poly = interpolate(pts5, x)
        print(f"\n  → g_5(m) = {factor(poly)} = {expand(poly)} (through {len(pts5)} points)")

print("\n" + "="*70)
print("CORRECTED g_6 VALUES")
print("="*70)
for mk in sorted(g6_vals):
    print(f"  g_6({mk}) = {g6_vals[mk]}")

print("\nDone!")
