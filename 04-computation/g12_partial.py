#!/usr/bin/env python3
"""
Compute g_12(3) from W(27) data.
CV²(n) = W(n)/n! - 1 = Σ_{k=1}^{⌊n/2⌋} g_k(n-2k) / (n)_{2k}
For n=27: k ranges 1..13. g_13(1) = 1.
g_12(3) = g_12(27-24) = g_12(3).

Known g_k coefficients for k=1..11.
Also: g_13(1) = 1 (universal), g_k(m) for m≤0 we don't need.
For k=14, m = 27-28 = -1 < 0, doesn't contribute (k ≤ ⌊n/2⌋ = 13).
"""
from fractions import Fraction as F
import math

# g_k coefficients: 3*g_k(m) = a*m³ + b*m² + c*m + d
gk_coeffs = {
    1: (0, 0, 3, 0),
    2: (0, 3, 0, 0),
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
    10: (1707771898925208, -9582908872031805, 16794323323619016, -8919186350512416),
    11: (1030345293948358170, -5802477398736520560, 10195015138571054553, -5422883033782892160),
}

def gk(k, m):
    """Compute g_k(m) exactly."""
    a, b, c, d = gk_coeffs[k]
    return F(a*m**3 + b*m**2 + c*m + d, 3)

def falling_fact(n, k):
    result = F(1)
    for i in range(k):
        result *= F(n - i)
    return result

W27 = 11692423738081050318010736030
n = 27
nfact = math.factorial(n)
CV2_27 = F(W27, nfact) - 1

print(f"W(27) = {W27}")
print(f"27! = {nfact}")
print(f"CV²(27) = W(27)/27! - 1 = {float(CV2_27):.15f}")

# Compute known g_k contributions for n=27
total_known = F(0)
for k in range(1, 12):  # k=1..11
    m = n - 2*k
    if m < 1:
        print(f"  k={k}: m={m} < 1, skipping")
        continue
    g_val = gk(k, m)
    ff = falling_fact(n, 2*k)
    contrib = g_val / ff
    total_known += contrib
    print(f"  k={k}: g_{k}({m}) = {g_val}, (27)_{{{2*k}}} = {ff}")
    print(f"         contrib = {float(contrib):.15e}")

# k=13: m = 27-26 = 1, g_13(1) = 1 (universal boundary)
g13_contrib = F(1) / falling_fact(27, 26)
print(f"  k=13: g_13(1) = 1, (27)_26 = {falling_fact(27, 26)}")
print(f"         contrib = {float(g13_contrib):.15e}")

# g_12(3) is unknown. k=12, m=27-24=3, (27)_24 = 27!/3!
g12_ff = falling_fact(27, 24)
print(f"\n  k=12: m=3, (27)_24 = {g12_ff}")

# CV²(27) = total_known + g_12(3)/g12_ff + g13_contrib
# g_12(3) = (CV²(27) - total_known - g13_contrib) * g12_ff

g12_at_3 = (CV2_27 - total_known - g13_contrib) * g12_ff
print(f"\n=== RESULT ===")
print(f"g_12(3) = {g12_at_3}")
print(f"g_12(3) = {float(g12_at_3):.6f}")

# Verify this is an integer (or at least rational with small denominator)
print(f"g_12(3) numerator = {g12_at_3.numerator}")
print(f"g_12(3) denominator = {g12_at_3.denominator}")

# From THM-217 binomial truncation:
# 3*g_12(m) = a_12*m³ + b_12*m² + c_12*m + d_12
# g_12(1) = 1 => a_12 + b_12 + c_12 + d_12 = 3
# g_12(2) = 24 => 8*a_12 + 4*b_12 + 2*c_12 + d_12 = 72
# g_12(3) = ? => 27*a_12 + 9*b_12 + 3*c_12 + d_12 = 3*g_12(3)
#
# From the binomial truncation structure:
# d_12 = 3*g_12(0)
# b_12 = 3*12 - 3 - 3*a_12 + 3*g_12(0)/2  [from constraint derivation]
# Actually: Δ¹ = 1 - Δ⁰, Δ² = Δ⁰ + 2(k-1), Δ³ = 2*a_k
# where Δ⁰ = g_k(0), Δ¹ = g_k(1)-g_k(0), Δ² = g_k(2)-2g_k(1)+g_k(0), Δ³ = g_k(3)-3g_k(2)+3g_k(1)-g_k(0)
#
# Δ³ = g_12(3) - 3*24 + 3*1 - g_12(0) = g_12(3) - 69 - g_12(0)
# And Δ³ = 2*a_12 / 3... wait, Δ³ is the third forward difference.
# In the binomial basis: g_k(m) = Δ⁰ C(m,0) + Δ¹ C(m,1) + Δ² C(m,2) + Δ³ C(m,3)
# With the constraints:
# g_k(0) = Δ⁰  (trivially)
# g_k(1) = Δ⁰ + Δ¹ = 1  =>  Δ¹ = 1 - g_k(0)
# g_k(2) = Δ⁰ + 2Δ¹ + Δ² = 2k  =>  Δ² = 2k - Δ⁰ - 2(1-Δ⁰) = 2k - 2 + Δ⁰
# So Δ² = g_k(0) + 2(k-1)
# g_k(3) = Δ⁰ + 3Δ¹ + 3Δ² + Δ³ = g_12(0) + 3(1-g_12(0)) + 3(g_12(0)+22) + Δ³
#         = g_12(0) + 3 - 3*g_12(0) + 3*g_12(0) + 66 + Δ³
#         = g_12(0) + 69 + Δ³

# So Δ³ = g_12(3) - g_12(0) - 69
# And Δ³ = 2*a_12 (the third binomial coefficient is related to the leading monomial coeff)

# Also, 3*g_k(m) = a_k*m³ + ..., and the third forward difference of 3*g_k is:
# Δ³(3*g_k) = 3*Δ³(g_k) = 6*a_k (since Δ³ of m³ = 6)
# Wait: Δ³(m³) = (m+3)³ - 3(m+2)³ + 3(m+1)³ - m³ at m=0: 27-24+3-0 = 6
# So Δ³(3g_k) = 6*a_k (from leading term), meaning Δ³(g_k) = 2*a_k

# Therefore: 2*a_12 = g_12(3) - g_12(0) - 69
# And we know g_12(3) from the computation above.
# This gives ONE equation with TWO unknowns: a_12, g_12(0).

print(f"\n=== Constraint from g_12(3) ===")
g12_3 = g12_at_3  # should be exact integer
print(f"2*a_12 = g_12(3) - g_12(0) - 69")
print(f"2*a_12 + g_12(0) = g_12(3) - 69 = {g12_3 - 69}")
D0_plus_D3 = g12_3 - 69
print(f"Δ⁰ + Δ³ = g_12(0) + 2*a_12 = {D0_plus_D3}")

# From pattern: g_k(0)/a_k → -2 slowly
# At k=11: ratio = -1.754
# Extrapolating: at k=12, ratio ≈ -1.765
# So g_12(0) ≈ -1.765 * a_12
# And 2*a_12 + (-1.765)*a_12 = D0_plus_D3
# 0.235 * a_12 = D0_plus_D3
# a_12 ≈ D0_plus_D3 / 0.235

print(f"\n=== Estimate a_12 from ratio pattern ===")
# Known ratios: g_k(0)/a_k for k=3..11
ratios = [0, -0.8, -1.526, -1.641, -1.682, -1.707, -1.725, -1.741, -1.754]
# k =        3,   4,     5,     6,     7,     8,     9,    10,    11
# Extrapolate: the increments are decreasing:
# -0.8, -0.726, -0.115, -0.041, -0.025, -0.018, -0.016, -0.013
# Next increment ≈ -0.011, so ratio_12 ≈ -1.765
for r_est in [-1.76, -1.765, -1.77]:
    # g_12(0) = r * a_12
    # g_12(0) + 2*a_12 = D0_plus_D3
    # r*a_12 + 2*a_12 = D0_plus_D3
    # a_12 = D0_plus_D3 / (r + 2)
    a12_est = float(D0_plus_D3) / (r_est + 2)
    g12_0_est = r_est * a12_est
    print(f"  ratio={r_est}: a_12 ≈ {a12_est:.0f}, g_12(0) ≈ {g12_0_est:.0f}")

# Also check: a_k/a_{k-1} pattern
# k=11: a_11/a_10 = 603.3
# Next: ~683 (from second difference extrapolation)
a12_from_ratio = 683 * 1030345293948358170
print(f"\n  From a_12/a_11 ≈ 683: a_12 ≈ {a12_from_ratio:.6e}")

# Cross-check
g12_0_cross = float(D0_plus_D3) - 2 * a12_from_ratio
ratio_cross = g12_0_cross / a12_from_ratio
print(f"  Implied g_12(0)/a_12 = {ratio_cross:.6f}")
