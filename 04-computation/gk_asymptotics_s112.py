#!/usr/bin/env python3
"""
gk_asymptotics_s112.py — Full asymptotic expansion of CV^2 via Delannoy
kind-pasteur-2026-03-15-S112

CV^2 = sum_{k>=1} 2*g_k(n-2k)/(n)_{2k}
     = sum_{k>=1} sum_j C(k-1,j-1)*C(n-2k,j)*2^j / (n)_{2k}

Expand in x = 1/n to get CV^2 = a_1*x + a_2*x^2 + a_3*x^3 + ...
"""

from fractions import Fraction
from math import comb, factorial

# ============================================================
# Formal power series in x = 1/n, to order x^P
# ============================================================

P = 15  # number of terms in the expansion

def fps_mult(a, b):
    """Multiply two formal power series (lists of Fraction coefficients)."""
    result = [Fraction(0)] * P
    for i in range(min(P, len(a))):
        if a[i] == 0:
            continue
        for j in range(min(P - i, len(b))):
            result[i + j] += a[i] * b[j]
    return result

def fps_inv(a):
    """Inverse of a formal power series with a[0] != 0."""
    result = [Fraction(0)] * P
    result[0] = Fraction(1) / a[0]
    for i in range(1, P):
        s = Fraction(0)
        for j in range(1, i + 1):
            if j < len(a) and a[j] != 0:
                s += a[j] * result[i - j]
        result[i] = -s / a[0]
    return result

def fps_from_poly(coeffs):
    """Convert polynomial [c0, c1, ...] to FPS."""
    result = [Fraction(0)] * P
    for i, c in enumerate(coeffs):
        if i < P:
            result[i] = Fraction(c)
    return result

# ============================================================
# Compute the k-th contribution as FPS in x = 1/n
# ============================================================

# g_k(n-2k) = sum_j C(k-1,j-1)*C(n-2k,j)*2^{j-1}
# (n)_{2k} = n(n-1)...(n-2k+1)
#
# Let m = n-2k. In terms of x = 1/n: m = 1/x - 2k.
# C(m, j) = product_{i=0}^{j-1} (m-i) / j! = product (1/x - 2k - i) / j!
# (n)_{2k} = product_{i=0}^{2k-1} (1/x - i) = (1/x)^{2k} * product (1 - ix)
#
# 2*g_k(m) / (n)_{2k} = sum_j C(k-1,j-1) * 2^j * C(m,j) / (n)_{2k}
#
# C(m,j)/(n)_{2k} = [prod_{i=0}^{j-1}(m-i)/j!] / [prod_{i=0}^{2k-1}(n-i)]
#                  = prod_{i=0}^{j-1}(n-2k-i) / (j! * prod_{i=0}^{2k-1}(n-i))
#                  = 1/j! * prod_{i=0}^{j-1}(n-2k-i) / (n)_{2k}
#                  = 1/j! * 1 / [prod_{i=0}^{2k-1}(n-i) / prod_{i=0}^{j-1}(n-2k-i)]
#                  = 1/j! * 1 / (n)_{2k} * (n-2k)_j
#                  = 1/j! * (n)_{2k+j} / ((n)_{2k})^2 ... no

# Simpler: C(m,j) / (n)_{2k} where m=n-2k
# = (m)_j / (j! * (n)_{2k})
# = (n-2k)_j / (j! * (n)_{2k})

# As a FPS in x = 1/n:
# n = 1/x, (n-i) = (1-ix)/x
# (n)_{2k} = prod_{i=0}^{2k-1} (1-ix)/x = x^{-2k} * prod_{i=0}^{2k-1}(1-ix)
# (n-2k)_j = prod_{i=0}^{j-1} (1-(2k+i)x)/x = x^{-j} * prod_{i=0}^{j-1}(1-(2k+i)x)

# C(m,j)/(n)_{2k} = x^{-j} * prod(1-(2k+i)x) / (j! * x^{-2k} * prod(1-ix))
#                  = x^{2k-j} / j! * prod_{i=0}^{j-1}(1-(2k+i)x) / prod_{i=0}^{2k-1}(1-ix)

# The FPS contribution of the (k,j) term to CV^2:
# 2^j * C(k-1,j-1) * x^{2k-j} / j! * [prod(1-(2k+i)x, i=0..j-1)] / [prod(1-ix, i=0..2k-1)]

# Note: 2k-j >= 2k-k = k >= 1, so all terms start at x^k or higher.

def contribution_k(k):
    """FPS for 2*g_k(n-2k)/(n)_{2k} in x = 1/n."""
    total = [Fraction(0)] * P

    for j in range(1, k + 1):
        coeff = Fraction(2**j * comb(k - 1, j - 1), factorial(j))
        power = 2 * k - j  # starting power of x

        if power >= P:
            continue

        # Numerator: prod_{i=0}^{j-1} (1 - (2k+i)*x)
        num = fps_from_poly([Fraction(1)])
        for i in range(j):
            factor = [Fraction(1), Fraction(-(2*k + i))]
            num = fps_mult(num, fps_from_poly(factor))

        # Denominator: prod_{i=1}^{2k-1} (1 - i*x)  [skip i=0 which gives 1]
        # Actually prod_{i=0}^{2k-1} (1-ix) = 1 * (1-x) * (1-2x) * ... * (1-(2k-1)x)
        den = fps_from_poly([Fraction(1)])
        for i in range(1, 2*k):
            factor = [Fraction(1), Fraction(-i)]
            den = fps_mult(den, fps_from_poly(factor))

        # Ratio = num / den
        ratio = fps_mult(num, fps_inv(den))

        # Multiply by coeff * x^power
        for p in range(P):
            if p - power >= 0 and p - power < P:
                total[p] += coeff * ratio[p - power]

    return total

# ============================================================
# Sum over k
# ============================================================

print("="*70)
print("ASYMPTOTIC EXPANSION: CV^2 = sum a_p / n^p")
print("="*70)

cv2_fps = [Fraction(0)] * P

for k in range(1, P):
    ck = contribution_k(k)
    for p in range(P):
        cv2_fps[p] += ck[p]

print("\nCV^2 = ", end="")
terms = []
for p in range(1, P):
    if cv2_fps[p] != 0:
        terms.append(f"({cv2_fps[p]})/n^{p}")
print(" + ".join(terms[:10]))

print("\nCoefficients a_p = [x^p] CV^2:")
for p in range(1, min(P, 12)):
    print(f"  a_{p} = {cv2_fps[p]} = {float(cv2_fps[p]):.10f}")

# ============================================================
# Check against known values
# ============================================================

print("\n" + "="*70)
print("VERIFICATION against exact CV^2")
print("="*70)

# W(n) via bitmask DP
def compute_W(n):
    if n == 1: return 1
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for length in range(1, n):
        new_dp = {}
        for (mask, last), weight in dp.items():
            if bin(mask).count('1') != length:
                continue
            for v in range(n):
                if mask & (1 << v):
                    continue
                if v == last - 1:
                    continue
                new_weight = weight * (2 if v == last + 1 else 1)
                key = (mask | (1 << v), v)
                new_dp[key] = new_dp.get(key, 0) + new_weight
        for k2, w in new_dp.items():
            dp[k2] = dp.get(k2, 0) + w
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

for n in range(3, 16):
    W = compute_W(n)
    cv2_exact = Fraction(W, factorial(n)) - 1

    # Evaluate FPS at x = 1/n
    x = Fraction(1, n)
    cv2_approx = sum(cv2_fps[p] * x**p for p in range(1, P))

    err = cv2_exact - cv2_approx
    print(f"  n={n:2d}: exact={float(cv2_exact):.10f}, "
          f"FPS={float(cv2_approx):.10f}, err={float(err):.2e}")

# ============================================================
# Pattern analysis of coefficients
# ============================================================

print("\n" + "="*70)
print("COEFFICIENT ANALYSIS")
print("="*70)

print("\na_p sequence:")
for p in range(1, min(P, 12)):
    a = cv2_fps[p]
    print(f"  a_{p} = {a}")

# Check: is a_p related to Bernoulli numbers?
# B_0=1, B_1=-1/2, B_2=1/6, B_4=-1/30, B_6=1/42
print("\nRatios a_{p+1}/a_p:")
for p in range(1, min(P-1, 10)):
    if cv2_fps[p] != 0:
        ratio = cv2_fps[p+1] / cv2_fps[p]
        print(f"  a_{p+1}/a_{p} = {ratio} = {float(ratio):.6f}")

# Check: a_1 = 2, a_2 = 0, a_3 = -14/3
# Is a_3 = -14/3? From THM-216.
print(f"\na_3 = {cv2_fps[3]} (THM-216 predicts -14/3 = {Fraction(-14,3)})")

# Factor the denominators
print("\nDenominators of a_p:")
for p in range(1, min(P, 12)):
    a = cv2_fps[p]
    if a != 0:
        print(f"  a_{p}: denominator = {a.denominator}")

# Product p! * a_p
print("\np! * a_p:")
for p in range(1, min(P, 12)):
    val = factorial(p) * cv2_fps[p]
    print(f"  {p}! * a_{p} = {val}")

# Check if a_p = 2*(-1)^{p+1} * something / p!
print("\na_p * p! / 2:")
for p in range(1, min(P, 12)):
    val = cv2_fps[p] * factorial(p) / 2
    print(f"  a_{p} * {p}!/2 = {val}")

print("\nDone!")
