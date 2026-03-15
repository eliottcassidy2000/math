#!/usr/bin/env python3
"""
gk_verify_89c.py — Cross-verify g_k polynomials against exact W(n) computation
opus-2026-03-15-S89c

For each n, compute CV² from W(n) and compare with the sum using all g_k polynomials.
"""

from fractions import Fraction
from math import factorial
from functools import reduce

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

# All g_k polynomials: 3·g_k(m) = a_k m³ + b_k m² + c_k m + d_k
gk_coeffs = {
    1: (0, 0, 3, 0),      # g_1(m) = m
    2: (0, 3, 0, 0),      # g_2(m) = m²
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
    10: (1707771898925208, -9582908872031805, 16794323323619016, -8919186350512416),
}

def g_eval(k, m):
    if k in gk_coeffs:
        a, b, c, d = gk_coeffs[k]
        return Fraction(a*m**3 + b*m**2 + c*m + d, 3)
    if m == 1:
        return Fraction(1)
    if m == 2:
        return Fraction(2*k)
    return None

# Compute and verify
print("="*70)
print("CROSS-VERIFICATION: W(n)/n! - 1 vs Σ 2·g_k(n-2k)/(n)_{2k}")
print("="*70)

for n in range(3, 22):
    W = compute_W_dp(n)
    nf = factorial(n)
    cv2_exact = Fraction(W, nf) - 1

    cv2_sum = Fraction(0)
    max_k = (n - 1) // 2
    all_known = True
    for k in range(1, max_k + 1):
        mk = n - 2*k
        if mk < 1:
            break
        g = g_eval(k, mk)
        if g is None:
            all_known = False
            break
        cv2_sum += Fraction(2 * g, falling_factorial(n, 2*k))

    residual = cv2_exact - cv2_sum
    if residual == 0 and all_known:
        print(f"n={n:2d}: EXACT MATCH ✓ (max_k={max_k}, all terms accounted)")
    elif all_known:
        print(f"n={n:2d}: RESIDUAL = {residual} = {float(residual):.6e}")
    else:
        print(f"n={n:2d}: residual = {float(residual):.6e} (missing k>{max(gk_coeffs)} terms)")

print("\nDone!")
