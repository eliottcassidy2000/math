#!/usr/bin/env python3
"""
Derive tr(c_{n-7}) at general n using moments m0-m5.

tr(c_{n-7}) = sum_P e_6(s_P) where e_6 is degree 6 in f.
This requires m0,...,m6. But we only have closed forms for m0-m5.

However, m6 depends on invariants beyond (t3, t5, bc):
- At n=7: m6 = 720*H + 850080*t3 + 40320*t5 + 80640*bc + 8636880
  H is the Hamiltonian path count, a new invariant.
- At general n: m6 involves sub-tournament H-counts at size 7+.

So tr(c_{n-7}) will depend on MORE than (t3, t5, bc). It will bring in
t7 and higher bc variants.

But we can still compute the (t3, t5, bc) part of tr(c_{n-7}) from our
existing formulas, and determine what NEW invariants enter.

Actually: tr(c_{n-7}) only exists for n >= 9 (since c_{n-7} = c_2 for n=9,
c_4 for n=11, etc.).

Strategy: compute tr(c_{n-7}) numerically at n=9 and n=11 to see the
structure, then connect to the known n=9 formula from THM-055.

opus-2026-03-06-S28
"""
from math import factorial, comb
from fractions import Fraction
from itertools import combinations, permutations
import random

def compute_e2k_coefficients(n, max_k):
    """Compute e_{2k} as polynomial in g = f - (n-1)/2."""
    N = n - 1
    def p(k):
        if k % 2 == 0: return {0: Fraction(N, 4**(k//2))}
        else: return {1: Fraction(1, 4**((k-1)//2))}
    def poly_mult(a, b):
        result = {}
        for da, ca in a.items():
            for db, cb in b.items():
                result[da+db] = result.get(da+db, Fraction(0)) + ca*cb
        return {d: c for d, c in result.items() if c}
    def poly_add(a, b, sign=1):
        result = dict(a)
        for d, c in b.items():
            result[d] = result.get(d, Fraction(0)) + sign*c
        return {d: c for d, c in result.items() if c}
    def poly_scale(a, s):
        return {d: c*s for d, c in a.items() if c*s}
    e = [{0: Fraction(1)}]
    for k_val in range(1, 2*max_k+1):
        total = {}
        for i in range(1, k_val+1):
            term = poly_mult(e[k_val-i], p(i))
            total = poly_add(total, term, (-1)**(i+1))
        e.append(poly_scale(total, Fraction(1, k_val)))
    return [e[2*k] for k in range(max_k+1)]

# =====================================================================
# What e_6 looks like symbolically
# =====================================================================
print("=" * 70)
print("e_6 polynomial in g = f - (n-1)/2")
print("=" * 70)

for n in [7, 9, 11]:
    e = compute_e2k_coefficients(n, 3)
    e6 = e[3]
    print(f"\n  n={n}: e_6(g) =")
    for deg in sorted(e6.keys()):
        print(f"    g^{deg}: {e6[deg]}")

# =====================================================================
# Compute tr(c_{n-7}) = sum_P e_6(s_P) symbolically
# =====================================================================
print(f"\n{'='*70}")
print("tr(c_{n-7}) from moments m0-m5 + correction from m6")
print(f"{'='*70}")

def moment_formula(j, n):
    """Return (const, t3, t5, bc) for m_j, or None if unknown."""
    if j == 0: return (factorial(n), 0, 0, 0)
    if j == 1: return (factorial(n)*(n-1)//2, 0, 0, 0)
    if j == 2: return (factorial(n)*(3*n*n-5*n+4)//12, 4*factorial(n-2), 0, 0)
    if j == 3: return (factorial(n)*(n-1)*(n*n-n+2)//8, 6*factorial(n-1), 0, 0)
    if j == 4: return (factorial(n)*(15*n**4-30*n**3+65*n**2-82*n+48)//240,
                       2*factorial(n-2)*(3*n*n-5*n+4),
                       48*factorial(n-4), 96*factorial(n-4))
    if j == 5: return (factorial(n)*(n-1)*(3*n**4-2*n**3+13*n**2-14*n+16)//96,
                       5*(n-1)*(n*n-n+2)*factorial(n-2),
                       120*(n-1)*factorial(n-4), 240*(n-1)*factorial(n-4))
    return None  # m6+ not available in closed form

# For tr(c_{n-7}), we need m0-m6. Compute what we can (m0-m5 contribution)
# and leave the m6 term as symbolic.

for n in [9, 11, 13]:
    print(f"\n  n={n}:")
    e = compute_e2k_coefficients(n, 3)
    e6 = e[3]
    half_n1 = Fraction(n-1, 2)

    # tr(c_{n-7}) = sum_P e_6(g_P) = sum_{deg in e6} e6[deg] * sum_P g^deg
    # sum_P g^deg = sum_{l=0}^{deg} C(deg,l) * (-half_n1)^{deg-l} * m_l

    # Separate into known (m0-m5) and unknown (m6) contributions
    result_known = [Fraction(0)] * 4  # const, t3, t5, bc
    m6_coefficient = Fraction(0)

    for deg, coeff in e6.items():
        for l in range(min(deg, 5) + 1):
            m = moment_formula(l, n)
            binom = comb(deg, l)
            factor = coeff * binom * (-half_n1)**(deg - l)
            for i in range(4):
                result_known[i] += factor * m[i]

        # m6 term (l=6, only if deg >= 6)
        if deg >= 6:
            binom = comb(deg, 6)
            factor = coeff * binom * (-half_n1)**(deg - 6)
            m6_coefficient += factor

    print(f"    Known part (from m0-m5):")
    print(f"      const = {float(result_known[0]):.2f}")
    print(f"      t3 = {float(result_known[1]):.2f}")
    print(f"      t5 = {float(result_known[2]):.2f}")
    print(f"      bc = {float(result_known[3]):.2f}")
    print(f"    m6 coefficient = {float(m6_coefficient):.6f}")

    # At n=9: tr(c_2) is known from THM-055:
    # c_2 = 462*t3 - 60*t5 + 12*t7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640
    if n == 9:
        print(f"\n    Expected (THM-055 n=9):")
        print(f"      c_2 = -2640 + 462*t3 - 60*t5 + 12*t7 - 120*bc + 24*bc35_w + 48*a3")
        print(f"\n    m6 must provide the difference between known part and this.")
        print(f"    If m6 = A + B*t3 + C*t5 + D*bc + E*t7 + F*bc35_w + G*a3 + ...")
        print(f"    then the m6_coeff * m6 fills the gap.")

# =====================================================================
# At n=7: tr(c_0) from e_6
# =====================================================================
print(f"\n{'='*70}")
print("Special case n=7: tr(c_0)")
print(f"{'='*70}")

n = 7
e = compute_e2k_coefficients(n, 3)
e6 = e[3]
half_n1 = Fraction(n-1, 2)

# For n=7, m6 = 720*H + stuff
# tr(c_0) should give us the formula for c_0 in terms of invariants + H.

result_known = [Fraction(0)] * 4
m6_coefficient = Fraction(0)

for deg, coeff in e6.items():
    for l in range(min(deg, 5) + 1):
        m = moment_formula(l, n)
        binom = comb(deg, l)
        factor = coeff * binom * (-half_n1)**(deg - l)
        for i in range(4):
            result_known[i] += factor * m[i]
    if deg >= 6:
        binom = comb(deg, 6)
        factor = coeff * binom * (-half_n1)**(deg - 6)
        m6_coefficient += factor

print(f"  Known part (from m0-m5):")
for i, name in enumerate(["const", "t3", "t5", "bc"]):
    print(f"    {name} = {result_known[i]} = {float(result_known[i]):.4f}")
print(f"  m6 coefficient = {m6_coefficient} = {float(m6_coefficient):.6f}")

# m6 at n=7 = 8636880 + 850080*t3 + 40320*t5 + 80640*bc + 720*H
# So the m6 contribution to tr(c_0) is m6_coefficient * (above)
m6_n7 = [8636880, 850080, 40320, 80640]  # const, t3, t5, bc + 720*H
H_coeff = m6_coefficient * 720

total = [result_known[i] + m6_coefficient * m6_n7[i] for i in range(4)]

print(f"\n  Complete tr(c_0) at n=7:")
print(f"    const = {float(total[0]):.4f}")
print(f"    t3 = {float(total[1]):.4f}")
print(f"    t5 = {float(total[2]):.4f}")
print(f"    bc = {float(total[3]):.4f}")
print(f"    H_coeff = {float(H_coeff):.4f}")

# Expected: tr(c_0) = 2*t3 - t5 + 2*t7 - 2*bc + 253/4
# But H = 1 + 2*(t3+t5+t7) + 4*bc at n=7
# So without OCF: c_0 = H_coeff*H + rest
# THM-055: c_0 = H - 6*bc - 3*t5 + 249/4 ... or
# c_0 = 2*t3 - t5 + 2*t7 - 2*bc + 253/4 (expanded via OCF)
print(f"\n  Expected: c_0 = 253/4 + 2*t3 - t5 + 2*t7 - 2*bc")
print(f"  Or equivalently: c_0 = 249/4 + H - 6*bc - 3*t5")
print(f"  (where H = 1 + 2*t3 + 2*t5 + 2*t7 + 4*bc at n=7)")

# =====================================================================
# KEY INSIGHT: universality of t5 and bc in tr(c_{n-5})
# =====================================================================
print(f"\n{'='*70}")
print("UNIVERSALITY of t5 and bc coefficients")
print(f"{'='*70}")
print("""
THEOREM: For ALL odd n >= 7:

  tr(c_{n-5}) = (n-4)!*C(n,5)*(5n-13)/12
                - (n-4)!*C(n-2,3) * t3
                + 2*(n-4)! * t5
                + 4*(n-4)! * bc

The t5 and bc coefficients are INDEPENDENT of n (up to the (n-4)! scaling):
  t5_coeff / (n-4)! = 2 (universal)
  bc_coeff / (n-4)! = 4 (universal)

This universality arises because:
1. t5 enters only through 5-vertex H-counts (pattern (4,))
2. bc enters only through 6-vertex structure (patterns (2,2) and (5,))
3. The position-pattern counting and Stirling aggregation exactly cancel
   all n-dependent factors except (n-4)!
4. The bc/t5 ratio = 2 traces to SIGMA_k ratios being exactly 2 for k <= 5

CONJECTURE: Similar universality holds for tr(c_{n-7}):
  The coefficient of t7 in tr(c_{n-7}) is a simple multiple of (n-6)!
  (and similarly for bc35_w, alpha_3, and other invariants at this level).
""")

print("DONE")
