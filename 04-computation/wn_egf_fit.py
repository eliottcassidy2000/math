#!/usr/bin/env python3
"""
Try to determine the EGF of W(n) by fitting.

Known: NUD EGF = exp(-x)/(1-x)^2
W(n)/NUD(n) -> e, so W(n) ~ e * NUD(n) ~ e * n!/e = n!

The EGF f(x) = Σ W(n) x^n/n! has f(x) -> 1/(1-x) as the leading singularity.

Since CV²(n) = W(n)/n! - 1 ~ 2/n, and Σ (c/n) x^n = -c ln(1-x),
we expect f(x) = 1/(1-x) - 2 ln(1-x) + g(x) where g is analytic at x=1.

Let me test this: define R(n) = W(n)/n! - 1 - 2/n (removing the known terms).
Then Σ R(n) x^n should be analytic at x=1 or have milder singularity.
"""
from fractions import Fraction as F
import math

wn = {1: 1, 2: 2, 3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752, 9: 439670,
      10: 4327904, 11: 46963358, 12: 556953448, 13: 7166360054,
      14: 99428495088, 15: 1479600188798, 16: 23506712352248,
      17: 397095175477430, 18: 7107209383674112, 19: 134345623603516190,
      20: 2674426516381764744, 21: 55925062706620208438,
      22: 1225582324497129993488, 23: 28088043261491650347134,
      24: 671901551280362054066968, 25: 16746599265666151628174198,
      26: 434187457363955400414175008, 27: 11692423738081050318010736030}

# CV²(n) = W(n)/n! - 1
cv2 = {}
for n in range(1, 28):
    cv2[n] = F(wn[n], math.factorial(n)) - 1

# R(n) = CV²(n) - 2/n (remove the known leading term)
R = {}
for n in range(1, 28):
    R[n] = cv2[n] - F(2, n)

print("=== R(n) = CV²(n) - 2/n ===")
for n in range(1, 28):
    print(f"  n={n:2d}: {float(R[n]):.15f}")

# R(n) -> 0 but slowly. Check R(n) * n:
print("\n=== R(n) * n ===")
for n in range(3, 28):
    print(f"  n={n:2d}: {float(R[n] * n):.10f}")

# R(n) * n -> 0 also. R(n) * n^2:
print("\n=== R(n) * n^2 ===")
for n in range(3, 28):
    print(f"  n={n:2d}: {float(R[n] * n * n):.10f}")

# From earlier: (CV² - 2/n) * n² → 0 slowly (~-0.20 at n=27)
# Maybe R(n) ~ -c * ln(n)/n² ?
print("\n=== R(n) * n^2 / ln(n) ===")
for n in range(3, 28):
    val = float(R[n] * n * n) / math.log(n)
    print(f"  n={n:2d}: {val:.10f}")

# Still going to 0. R(n) * n^2 / ln(n)^2?
print("\n=== R(n) * n^2 / ln(n)^2 ===")
for n in range(3, 28):
    val = float(R[n] * n * n) / math.log(n)**2
    print(f"  n={n:2d}: {val:.10f}")

# Hmm, that's also going to 0. Let me try R(n) * n^1.5:
print("\n=== R(n) * n^1.5 ===")
for n in range(5, 28):
    val = float(R[n]) * n**1.5
    print(f"  n={n:2d}: {val:.10f}")

# Let me try a different decomposition. Since the g_k contribution gives
# CV²(n) = Σ_{k=1}^{...} g_k(n-2k) / (n)_{2k}
# and g_1(m) = m gives the leading 2/n term,
# the next term is g_2(m) = m(m+1)/2... wait, let me check.
# 3*g_k = a*m³+b*m²+c*m+d
# g_1: (0,0,3,0) → g_1(m) = m
# g_2: (0,3,0,0) → g_2(m) = m²/3... wait
# 3*g_2(m) = 3m², so g_2(m) = m²
# g_2(m)/(n)_4 where m = n-4, (n)_4 = n(n-1)(n-2)(n-3)
# g_2(n-4) = (n-4)²
# Contribution: (n-4)²/[n(n-1)(n-2)(n-3)]
# For large n: ~ n²/n⁴ = 1/n² → 0.
# So the g_2 term decays as 1/n². Not enough to explain the -2/n² + ... pattern.
# Actually there are O(n/2) terms in the sum, so the total is more complex.

# Let me compute the first few g_k contributions explicitly
print("\n=== g_k contributions to CV²(n) ===")
# g_k coefficients: 3*g_k = a*m³+b*m²+c*m+d, m = n-2k
gk_coeffs = {
    1: (0, 0, 3, 0),
    2: (0, 3, 0, 0),
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
}

def gk(k, m):
    if k not in gk_coeffs:
        return None
    a, b, c, d = gk_coeffs[k]
    return F(a*m**3 + b*m**2 + c*m + d, 3)

def falling_fact(n, k):
    """(n)_k = n*(n-1)*...*(n-k+1)"""
    result = F(1)
    for i in range(k):
        result *= F(n - i)
    return result

# CV²(n) = Σ_k g_k(n-2k) / (n)_{2k}
for n in [10, 15, 20, 25, 27]:
    print(f"\n  n={n}:")
    total = F(0)
    for k in range(1, 6):
        m = n - 2*k
        if m < 0:
            break
        gval = gk(k, m)
        ff = falling_fact(n, 2*k)
        if ff != 0:
            contrib = gval / ff
            total += contrib
            print(f"    k={k}: g_{k}({m})={gval}, (n)_{2*k}={ff}, contrib={float(contrib):.15f}")
    # Remaining terms (k=6+) that we don't have g_k for
    cv2_exact = F(wn[n], math.factorial(n)) - 1
    remainder = cv2_exact - total
    print(f"    Sum(k=1..5)={float(total):.15f}")
    print(f"    CV²(exact) ={float(cv2_exact):.15f}")
    print(f"    Remainder  ={float(remainder):.15f}")
    print(f"    Remainder*n²={float(remainder*n*n):.10f}")

# The NUD formula gives NUD(n) via its EGF.
# W(n) = Σ_σ∈NUD 2^adj1(σ). We know the EGF of NUD is exp(-x)/(1-x)².
# The bivariate EGF F(x,u) = Σ_n Σ_j N(n,j) u^j x^n/n!
# with F(x,1) = exp(-x)/(1-x)².
# We want F(x,2).
#
# Key idea: F(x,u) might be related to exp((u-2)x)/(1-x)^2 * correction.
# At u=1: exp(-x)/(1-x)² ✓
# At u=0: should give Hertzsprung EGF.
#
# exp(-2x)/(1-x)² at n=1: (-2) + 2 = 0. But N(1,0) = 1, Hertz(1) = 1. So this fails.
# The simple exp((u-2)x)/(1-x)² doesn't work.

# Let me try to see if F(x,2) has a nice form.
# f(x) = F(x,2) = Σ W(n) x^n / n!
# Evaluate f at a few points to see structure:
print("\n=== Numerical evaluation of f(x) = EGF of W ===")
for x_val in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    f_val = sum(wn[n] / math.factorial(n) * x_val**n for n in range(1, 28))
    # Compare with some candidates:
    pole = 1/(1-x_val)
    log = -2 * math.log(1-x_val)
    nud_egf = math.exp(-x_val) / (1-x_val)**2
    print(f"  x={x_val:.1f}: f(x)={f_val:.10f}, 1/(1-x)={pole:.10f}, "
          f"-2ln(1-x)={log:.10f}, nud_egf={nud_egf:.10f}, "
          f"f-pole={f_val-pole:.10f}")
