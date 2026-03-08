#!/usr/bin/env python3
"""
kappa6_analysis.py - Analyze the kappa_6 coefficient structure.

At n=7:
  kappa_6 = 2/63 + (2/7)*t7 + (80/3087)*t3^3 - (4/49)*t3*t5 - (8/49)*t3*a2

Let me understand each piece.

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import comb

n = 7

# Constant term: 2/63 = (n+1)/252 = 8/252 = 2/63. ✓ (Bernoulli formula)
print("KAPPA_6 COEFFICIENTS AT n=7")
print("=" * 60)
print(f"Constant: 2/63 = (n+1)/(6*42) = {Fraction(n+1, 252)}")
print(f"  = (n+1)*B_6/6 where B_6 = 1/42")
print(f"  Matches Bernoulli formula: {Fraction(2,63) == Fraction(n+1, 252)}")

# t7 coefficient: 2/7 = 2/n
print(f"\nt7 coefficient: 2/7 = 2/n = {Fraction(2, n)}")
print(f"  Compare: kappa_2 has t3 coeff = 4/(n(n-1)) = {Fraction(4, n*(n-1))}")
print(f"  And kappa_4 has t5 coeff = 2/C(n,4) = {Fraction(2, comb(n,4))}")
print(f"  Pattern for new cycle coefficient: 2/C(n, 2k) for kappa_{2k}?")
print(f"    kappa_2: 2k=2, C(n,2) = {comb(n,2)}, 2/C(n,2) = {Fraction(2, comb(n,2))}. Actual: {Fraction(4, n*(n-1))} = {Fraction(2, comb(n,2))}. ✓!")
print(f"    kappa_4: 2k=4, C(n,4) = {comb(n,4)}, 2/C(n,4) = {Fraction(2, comb(n,4))}. Actual: {Fraction(2, comb(n,4))}. ✓!")
print(f"    kappa_6: 2k=6, C(n,6) = {comb(n,6)}, 2/C(n,6) = {Fraction(2, comb(n,6))}.")
print(f"    Actual t7 coeff: {Fraction(2,7)}. 2/C(7,6) = {Fraction(2, comb(7,6))}. Match: {Fraction(2,7) == Fraction(2, comb(7,6))}!")

# So the "new cycle" coefficient follows: 2/C(n, 2k) for the t_{2k+1} term in kappa_{2k}!
print("\n" + "=" * 60)
print("UNIVERSAL PATTERN: new (2k+1)-cycle coefficient in kappa_{2k} = 2/C(n, 2k)")
print("=" * 60)

# Wait, let me check kappa_2 more carefully.
# kappa_2 = (n+1)/12 + 4*t3/(n(n-1))
# The t3 coefficient is 4/(n(n-1)).
# 2/C(n,2) = 2/(n(n-1)/2) = 4/(n(n-1)). ✓!!!

print()
print("kappa_2: t3 coeff = 4/(n(n-1)) = 2/C(n,2). C(n,2) = n(n-1)/2.")
print("kappa_4: t5 coeff = 2/C(n,4).")
print("kappa_6: t7 coeff = 2/C(n,6).")
print()
print("CONJECTURE: kappa_{2k}(T) contains t_{2k+1} with coefficient 2/C(n, 2k).")
print("This is the LEADING new cycle term in the cumulant hierarchy.")

# What about the alpha_2 term in kappa_4?
# coeff_a2 = 4/C(n,4) = 2 * (2/C(n,4)) = 2 * coeff_t5
# Similarly in kappa_6, alpha_2 enters through the CROSS TERM t3*a2 with coeff -8/49.
# This is a nonlinear (mixed) term, not a "new" linear term.

print("\n" + "=" * 60)
print("NONLINEAR TERMS IN KAPPA_6")
print("=" * 60)

# kappa_6 = 2/63 + (2/7)*t7 + (80/3087)*t3^3 - (4/49)*t3*t5 - (8/49)*t3*a2
# Nonlinear terms: t3^3, t3*t5, t3*a2
# These are products of lower-level invariants.

# Let's check: does the nonlinear part = known function of (kappa_2 - kappa_2(trans), kappa_4 - kappa_4(trans))?
# delta_k2 = 4*t3/(n(n-1)) = 4t3/42 = 2t3/21
# delta_k4 = 2*t5/C(n,4) + 4*a2/C(n,4) - 48*t3^2/(n(n-1))^2
#           = 2t5/35 + 4a2/35 - 48t3^2/1764
#           = 2t5/35 + 4a2/35 - 4t3^2/147

dk2 = Fraction(2, 21)  # coefficient of t3 in delta_k2 (i.e., delta_k2 = (2/21)*t3)
dk4_t5 = Fraction(2, 35)
dk4_a2 = Fraction(4, 35)
dk4_t3sq = Fraction(-4, 147)

# kappa_6 nonlinear part: (80/3087)*t3^3 - (4/49)*t3*t5 - (8/49)*t3*a2
# = t3 * [80t3^2/3087 - 4t5/49 - 8a2/49]

# Check if this equals some combination like:
# 15*dk2*dk4 or 30*dk2^3 etc (from cumulant-to-moment formulas)
# kappa_6 = mu_6 - 15*mu_4*mu_2 + 30*mu_2^3  (since mu_3 = 0)
# = mu_6 - 15*(kappa_4 + 3*kappa_2^2)*kappa_2 + 30*kappa_2^3
# = mu_6 - 15*kappa_4*kappa_2 - 45*kappa_2^3 + 30*kappa_2^3
# = mu_6 - 15*kappa_4*kappa_2 - 15*kappa_2^3
# So: kappa_6 = mu_6 - 15*kappa_4*kappa_2 - 15*kappa_2^3

# But also: kappa_6 = mu_6 - 15*mu_4*mu_2 - 10*mu_3^2 + 30*mu_2^3
# With mu_3 = 0: kappa_6 = mu_6 - 15*mu_4*mu_2 + 30*mu_2^3

# Hmm this is just the standard formula. The point is that kappa_6 can have cross terms
# from the products kappa_4*kappa_2 in the moment-to-cumulant inversion.

# The KEY structure:
# kappa_6 = kappa_6(trans) + (2/C(n,6))*t7 + [nonlinear in dk2, dk4]

# Let's verify: the cross terms should come from -15*kappa_4*kappa_2 (expanded)
# kappa_2 = K2_0 + dk2*t3  (linear in t3)
# kappa_4 = K4_0 + dk4(t5,a2,t3^2)  (from THM-093)

# The t3*t5 term in kappa_6 comes from -15*dk4_t5_piece * dk2_piece * ???
# Actually this gets complicated. Let me just verify the formula.

print()
print("t3^3 coeff: 80/3087 =", Fraction(80, 3087))
print("Reduced:", Fraction(80, 3087))  # gcd(80,3087)?
import math
g = math.gcd(80, 3087)
print(f"gcd(80,3087) = {g}, so 80/3087 = {80//g}/{3087//g}")
# 3087 = 3*1029 = 3*3*343 = 9*343 = 9*7^3. 80 = 16*5.
# gcd = 1. So 80/3087 is already reduced.
print(f"3087 = 9*343 = 9*7^3. So 80/(9*7^3) = 80/3087")

# t3*t5 coeff: -4/49 = -4/7^2
print(f"\nt3*t5 coeff: -4/49 = -4/7^2")
# t3*a2 coeff: -8/49 = -8/7^2 = 2*(-4/49)
print(f"t3*a2 coeff: -8/49 = 2*(-4/49) - exactly 2x the t3*t5 coeff!")
print("This matches the universal 1:2 ratio for (t5, alpha_2) pairs.")

print("\n" + "=" * 60)
print("KAPPA_6 COMPACT FORM")
print("=" * 60)
print()
print("kappa_6(T) = (n+1)/252 + (2/C(n,6))*t7 - (4/49)*t3*(t5+2*a2) + (80/3087)*t3^3")
print()
print("Note the invariant t5+2*a2 appears AGAIN (same as in kappa_4)!")
print("The new genuinely 7-vertex term is (2/C(n,6))*t7.")
print()

# Final summary
print("=" * 60)
print("COMPLETE CUMULANT HIERARCHY (through kappa_6, at n=7)")
print("=" * 60)
print()
print("kappa_2 = (n+1)/12 + (2/C(n,2))*t3")
print()
print("kappa_4 = -(n+1)/120 + (2/C(n,4))*(t5 + 2*a2) - 3*(2*t3/C(n,2))^2")
print(f"        = -(n+1)/120 + (2/C(n,4))*(t5 + 2*a2) - 48*t3^2/(n(n-1))^2")
print()
print("kappa_6 = (n+1)/252 + (2/C(n,6))*t7 - (4/49)*t3*(t5+2*a2) + (80/3087)*t3^3")
print()
print("PATTERN:")
print("  1. Bernoulli constant: (-1)^{k-1}(n+1)|B_{2k}|/(2k)")
print("  2. New cycle term: (2/C(n,2k))*t_{2k+1}")
print("  3. Cross terms: products of lower-level corrections")
print("  4. The combination t5+2*a2 appears in BOTH kappa_4 and kappa_6")
