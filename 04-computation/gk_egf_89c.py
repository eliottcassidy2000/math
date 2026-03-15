#!/usr/bin/env python3
"""
gk_egf_89c.py — EGF analysis of g_k coefficients
opus-2026-03-15-S89c

KEY STRUCTURAL THEOREM:
  g_k(m) = g_k(0)·C(m-1,2) + m + 2(k-1)·C(m,2) + 2·a_k·C(m,3)

where:
  - g_k(0) = polynomial extrapolation to m=0 (free parameter)
  - a_k = leading cubic coefficient in 3·g_k = a_k·m³ + ...

The ENTIRE polynomial is determined by just (g_k(0), a_k).
Both sequences are new to OEIS.

g_k(0): 0, 0, 0, -8, -592, -114320, -33338240, -12475185560, -5629549881808
a_k:    0, 0, 2, 10, 388, 69660, 19826270, 7309726742, 3262687720240

This script investigates EGF/recurrence structure of these sequences.
"""

from fractions import Fraction
from math import factorial, comb
from sympy import symbols, Rational, factor, expand, Poly, solve

gk_coeffs = {
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}

def g_true(k, m):
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m*m)
    a, b, c, d = gk_coeffs[k]
    return Fraction(a*m**3 + b*m**2 + c*m + d, 3)

# Extract the two key sequences
ak_seq = {1: 0, 2: 0}
gk0_seq = {1: 0, 2: 0}
for k in range(1, 10):
    gk0_seq[k] = g_true(k, 0)
    if k >= 3:
        ak_seq[k] = gk_coeffs[k][0]
    elif k == 1:
        ak_seq[k] = 0
    elif k == 2:
        ak_seq[k] = 0

print("="*70)
print("KEY SEQUENCES")
print("="*70)

print("\na_k (leading cubic coeff in 3*g_k):")
for k in range(1, 10):
    print(f"  k={k}: a_k = {ak_seq[k]}")

print("\ng_k(0):")
for k in range(1, 10):
    print(f"  k={k}: g_k(0) = {gk0_seq[k]}")

# EGF: A(t) = Σ a_k · t^k / k!  or  Σ a_k · t^k / (2k)!  ?
print("\n" + "="*70)
print("EGF CANDIDATES FOR a_k")
print("="*70)

print("\na_k / k!:")
for k in range(3, 10):
    r = Fraction(ak_seq[k], factorial(k))
    print(f"  k={k}: {float(r):.15f}")

print("\na_k / (2k)!:")
for k in range(3, 10):
    r = Fraction(ak_seq[k], factorial(2*k))
    print(f"  k={k}: {float(r):.15e}")

print("\na_k / k!²:")
for k in range(3, 10):
    r = Fraction(ak_seq[k], factorial(k)**2)
    print(f"  k={k}: {float(r):.15f}")

print("\na_k / (k-1)!²:")
for k in range(3, 10):
    r = Fraction(ak_seq[k], factorial(k-1)**2)
    print(f"  k={k}: {float(r):.15f}")

# Look for recurrence: a_k = f(k) * a_{k-1} + g(k) * a_{k-2}?
print("\n" + "="*70)
print("RECURRENCE SEARCH FOR a_k")
print("="*70)

ak = [ak_seq[k] for k in range(3, 10)]
# Try a_k = α(k)·a_{k-1} + β(k)·a_{k-2}
# k=5: 388 = α·10 + β·2 → 5α + β = 194
# k=6: 69660 = α·388 + β·10 → 388α + 10β = 69660
# From first: β = 194 - 5α
# Sub: 388α + 10(194-5α) = 69660 → 388α + 1940 - 50α = 69660 → 338α = 67720 → α = 200.35...
# Not integer. Try with k-dependent coefficients.

# Try a_k = (αk² + βk + γ)·a_{k-1}
# k=4: 10 = (16α+4β+γ)·2
# k=5: 388 = (25α+5β+γ)·10
# k=6: 69660 = (36α+6β+γ)·388
# From k=4: 16α+4β+γ = 5
# From k=5: 25α+5β+γ = 38.8 → not integer either

# Try a_k = (αk+β)·a_{k-1} + (γk+δ)·a_{k-2}
# 4 unknowns, need 4 equations
# k=5: 388 = (5α+β)·10 + (5γ+δ)·2
# k=6: 69660 = (6α+β)·388 + (6γ+δ)·10
# k=7: 19826270 = (7α+β)·69660 + (7γ+δ)·388
# k=8: 7309726742 = (8α+β)·19826270 + (8γ+δ)·69660

# Let A = 5α+β, B = 5γ+δ:
# 388 = 10A + 2B → 5A + B = 194
# Let A' = 6α+β, B' = 6γ+δ (= A+α, B+γ):
# 69660 = 388A' + 10B'

# This is getting messy. Let me try numerically.
from sympy import Matrix, Rational as R

# System: a_k = p(k)·a_{k-1} + q(k)·a_{k-2}
# where p(k) = p1·k + p0, q(k) = q1·k + q0
# k=5: 388 = (5p1+p0)·10 + (5q1+q0)·2
# k=6: 69660 = (6p1+p0)·388 + (6q1+q0)·10
# k=7: 19826270 = (7p1+p0)·69660 + (7q1+q0)·388
# k=8: 7309726742 = (8p1+p0)·19826270 + (8q1+q0)·69660

M = Matrix([
    [50, 10, 10, 2],
    [6*388, 388, 60, 10],
    [7*69660, 69660, 7*388, 388],
    [8*19826270, 19826270, 8*69660, 69660],
])
b = Matrix([388, 69660, 19826270, 7309726742])
sol = M.solve(b)
p1, p0, q1, q0 = sol
print(f"\nLinear recurrence a_k = (p1·k+p0)·a_(k-1) + (q1·k+q0)·a_(k-2):")
print(f"  p1 = {p1}")
print(f"  p0 = {p0}")
print(f"  q1 = {q1}")
print(f"  q0 = {q0}")

# Verify
print("\nVerification:")
for k in range(5, 10):
    pred = (p1*k + p0) * ak_seq[k-1] + (q1*k + q0) * ak_seq[k-2]
    actual = ak_seq[k]
    print(f"  k={k}: pred={pred}, actual={actual}, match={pred==actual}")

# Try quadratic coefficients: a_k = (p2k²+p1k+p0)·a_{k-1} + (q2k²+q1k+q0)·a_{k-2}
print("\n" + "="*70)
print("QUADRATIC RECURRENCE SEARCH")
print("="*70)

M2 = Matrix([
    [25*10, 50, 10, 25*2, 10, 2],
    [36*388, 6*388, 388, 36*10, 60, 10],
    [49*69660, 7*69660, 69660, 49*388, 7*388, 388],
    [64*19826270, 8*19826270, 19826270, 64*69660, 8*69660, 69660],
])
# Overdetermined with 4 equations and 6 unknowns if we only have up to k=8
# Add k=9
M2_full = Matrix([
    [25*10, 50, 10, 25*2, 10, 2],
    [36*388, 6*388, 388, 36*10, 60, 10],
    [49*69660, 7*69660, 69660, 49*388, 7*388, 388],
    [64*19826270, 8*19826270, 19826270, 64*69660, 8*69660, 69660],
    [81*7309726742, 9*7309726742, 7309726742, 81*19826270, 9*19826270, 19826270],
])
# Wait, this needs k=10 which we don't have. Let me try only 3-term recurrence with quadratic coefficients.
# Actually with 6 unknowns we need 6 equations (k=5..10), but only have up to k=9.
# Let me try: does the simple linear recurrence already work? Check k=9.

if p1 is not None:
    pred9 = (p1*9 + p0) * ak_seq[8] + (q1*9 + q0) * ak_seq[7]
    print(f"k=9 prediction from linear recurrence: {pred9}")
    print(f"k=9 actual: {ak_seq[9]}")
    print(f"Match: {pred9 == ak_seq[9]}")

# Also try: a_k = f(k) * a_{k-1} alone (no a_{k-2} term)
print("\n" + "="*70)
print("RATIO a_k/a_{k-1} ANALYSIS")
print("="*70)

for k in range(4, 10):
    r = Fraction(ak_seq[k], ak_seq[k-1])
    print(f"  a_{k}/a_{k-1} = {r} = {float(r):.6f}")

# Is a_k/a_{k-1} a polynomial in k?
ratios = [Fraction(ak_seq[k], ak_seq[k-1]) for k in range(4, 10)]
print(f"\nRatios: {[float(r) for r in ratios]}")

# Differences of ratios
diffs = [ratios[i+1] - ratios[i] for i in range(len(ratios)-1)]
print(f"First diffs: {[float(d) for d in diffs]}")
diffs2 = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
print(f"Second diffs: {[float(d) for d in diffs2]}")

# Now look at g_k(0) similarly
print("\n" + "="*70)
print("g_k(0) ANALYSIS")
print("="*70)

print("\ng_k(0) / a_k ratios:")
for k in range(4, 10):
    r = Fraction(int(gk0_seq[k]), ak_seq[k])
    print(f"  k={k}: g_k(0)/a_k = {r} = {float(r):.10f}")

print("\ng_k(0) / g_{k-1}(0) ratios:")
for k in range(5, 10):
    r = Fraction(int(gk0_seq[k]), int(gk0_seq[k-1]))
    print(f"  k={k}: {float(r):.6f}")

# Relationship between g_k(0) and a_k?
# Recall: 3·g_k(m) = a_k·m³ + b_k·m² + c_k·m + d_k
# g_k(0) = d_k/3
# g_k(1) = 1 → a_k + b_k + c_k + d_k = 3
# g_k(2) = 2k → 8a_k + 4b_k + 2c_k + d_k = 6k
# So: d_k = 3·g_k(0), and a_k + b_k + c_k = 3 - 3·g_k(0)

print("\nb_k and c_k in terms of a_k and g_k(0):")
for k in range(3, 10):
    a, b, c, d = gk_coeffs[k]
    g0 = g_true(k, 0)
    print(f"  k={k}: a={a}, b={b}, c={c}, d={d}")
    print(f"    d/3 = {Fraction(d,3)} = g_k(0) = {g0}")
    # From g_k(1)=1: a+b+c+d=3 → b+c = 3-d-a
    # From g_k(2)=2k: 8a+4b+2c+d=6k → 4b+2c = 6k-d-8a → 2b+c = 3k-(d+8a)/2
    # So b = (2b+c) - (b+c) = [3k-(d+8a)/2] - [3-d-a]
    #      = 3k - d/2 - 4a - 3 + d + a = 3k + d/2 - 3a - 3
    b_formula = 3*k + d/2 - 3*a - 3
    print(f"    b formula check: {b_formula} = {b}? {b_formula == b}")

print("\n" + "="*70)
print("COMPACT FORMULA: b_k = 3k + 3g_k(0)/2 - 3a_k - 3")
print("="*70)

# Hmm wait. d = 3·g_k(0). So d/2 = 3·g_k(0)/2. Hmm, is 3·g_k(0) always even?
for k in range(3, 10):
    d = gk_coeffs[k][3]
    print(f"  k={k}: d={d}, d mod 2 = {d % 2}")

print("\nSo d is always even (d/2 is integer), confirming b_k = 3(k-1) + d/2 - 3a_k")
print("Or equivalently: b_k = 3(k-1) + 3g_k(0)/1·(3/2) - 3a_k")

# Let me just print the clean relationship
print("\n" + "="*70)
print("RELATIONSHIP BETWEEN b_k, a_k, g_k(0)")
print("="*70)
for k in range(3, 10):
    a, b, c, d = gk_coeffs[k]
    # b = 3k - 3 + d//2 - 3a
    expected_b = 3*k - 3 + d//2 - 3*a
    print(f"  k={k}: b={b}, 3(k-1)+d/2-3a = {expected_b}, match = {b == expected_b}")

print("\nDone!")
