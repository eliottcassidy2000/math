#!/usr/bin/env python3
"""
kappa4_coeff_pattern.py — Find closed-form pattern for kappa_4 coefficients.

Data:
  n=5: coeff_t5 = 2/5,  coeff_a2 = 0,    t3^2 = -3/25
  n=6: coeff_t5 = 2/15, coeff_a2 = 4/15, t3^2 = -4/75
  n=7: coeff_t5 = 2/35, coeff_a2 = 4/35, t3^2 = -4/147

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import comb

print("COEFFICIENT PATTERNS IN kappa_4")
print("=" * 60)

# t5 coefficient: 2/5, 2/15, 2/35
print("\nt5 coefficient:")
for n in [5, 6, 7]:
    if n == 5: c = Fraction(2, 5)
    elif n == 6: c = Fraction(2, 15)
    elif n == 7: c = Fraction(2, 35)
    print(f"  n={n}: {c} = 2/{1/float(c/2):.0f}")

# Denominators: 5, 15, 35
# = 5*1, 5*3, 5*7?
# Or: 1*5, 3*5, 5*7 — products of consecutive odds
# 1*5, 3*5, 5*7 = C(5,4), C(6,4)*?, ...
# Actually: C(n,4)*4!/something?
# C(5,4)=5, C(6,4)=15, C(7,4)=35. YES!

print("\n  C(n,4): C(5,4)=5, C(6,4)=15, C(7,4)=35")
print("  Pattern: coeff_t5 = 2/C(n,4)")

for n in [5, 6, 7]:
    predicted = Fraction(2, comb(n, 4))
    if n == 5: actual = Fraction(2, 5)
    elif n == 6: actual = Fraction(2, 15)
    elif n == 7: actual = Fraction(2, 35)
    print(f"  n={n}: predicted 2/C(n,4) = {predicted}, actual = {actual}, match = {predicted == actual}")

# alpha_2 coefficient: 0, 4/15, 4/35
print("\nalpha_2 coefficient:")
for n in [5, 6, 7]:
    if n == 5: c = Fraction(0)
    elif n == 6: c = Fraction(4, 15)
    elif n == 7: c = Fraction(4, 35)
    print(f"  n={n}: {c}")

# At n=5, alpha_2=0 always (can't have 2 disjoint 3-cycles on 5 vertices)
# So the coefficient is "undefined" or "0 trivially"
# For n=6,7: 4/15, 4/35 = 4/C(6,4), 4/C(7,4). Pattern: 4/C(n,4)?

print("\n  Pattern: coeff_a2 = 4/C(n,4)?")
for n in [6, 7]:
    predicted = Fraction(4, comb(n, 4))
    if n == 6: actual = Fraction(4, 15)
    elif n == 7: actual = Fraction(4, 35)
    print(f"  n={n}: predicted 4/C(n,4) = {predicted}, actual = {actual}, match = {predicted == actual}")

# So: coeff_t5 = 2/C(n,4), coeff_a2 = 4/C(n,4)
# Note: alpha_2 counts PAIRS of disjoint 3-cycles. Its coefficient is exactly 2x the t5 coefficient!

# t3^2 coefficient: -3/25, -4/75, -4/147
print("\nt3^2 coefficient:")
for n in [5, 6, 7]:
    coeff = Fraction(-48, (n*(n-1))**2)
    print(f"  n={n}: -48/(n(n-1))^2 = {coeff}")

# ALSO check M4 formula
# M4 = M4(trans) + dM4/dt3 * t3 + coeff_t5 * t5 + coeff_a2 * a2
# Since kappa_4 inherits t5 and a2 coefficients from M4 (M3, M2 don't depend on t5, a2),
# the M4 coefficients are:
# coeff_t5(M4) = coeff_t5(kappa_4) = 2/C(n,4)
# coeff_a2(M4) = coeff_a2(kappa_4) = 4/C(n,4)

print("\n" + "=" * 60)
print("COMPLETE kappa_4 FORMULA FOR GENERAL n")
print("=" * 60)
print()
print("kappa_4(T) = -(n+1)/120 + (2/C(n,4))*t5 + (4/C(n,4))*alpha_2 - 48/(n(n-1))^2 * t3^2")
print()
print("Equivalently:")
print("kappa_4(T) = -(n+1)/120 + (2/C(n,4))*(t5 + 2*alpha_2) - 48/(n(n-1))^2 * t3^2")
print()

# Note: t5 + 2*alpha_2 has a nice interpretation!
# alpha_2 = number of pairs of vertex-disjoint 3-cycles
# t5 = number of directed 5-cycles
# t5 + 2*alpha_2 = ?

print("INVARIANT: t5 + 2*alpha_2")
print("  This combines 5-cycles with disjoint 3-cycle pairs.")
print("  Both involve exactly 5 vertices (at n>=6, alpha_2 uses 6 vertices).")
print("  Wait... actually alpha_2 uses 6 vertices (two disjoint triples).")
print("  So at n=5, alpha_2=0 trivially.")
print("  But 2*alpha_2 has the same coefficient structure as t5.")
print()

# Actually let's think about this differently.
# In the OCF, alpha_2 contributes 4*alpha_2 to H.
# In kappa_4, 4*alpha_2/C(n,4).
# t5 contributes 2*t5 to H. In kappa_4, 2*t5/C(n,4).
# The ratio is preserved: coeff_a2/coeff_t5 = 2 in both OCF and kappa_4.
# This is significant: the kappa_4 dependence on (t5, alpha_2) mirrors the OCF dependence!

print("=" * 60)
print("CONNECTION TO OCF")
print("=" * 60)
print()
print("OCF: H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ... = 1 + 2*t3 + 2*t5 + 4*alpha_2 + ...")
print("       (at level of t5, alpha_2 contributions)")
print()
print("kappa_4 dependence on (t5, alpha_2):")
print("  = (1/C(n,4)) * (2*t5 + 4*alpha_2)")
print("  = (1/C(n,4)) * [OCF contribution from 5-vertex independent sets in Omega(T)]")
print()
print("This is remarkable: kappa_4 captures the 5-vertex OCF structure,")
print("normalized by C(n,4) = number of 4-element subsets.")

# E[fwd^4] formula at general n
print("\n" + "=" * 60)
print("E[fwd^4] FORMULA FOR GENERAL n")
print("=" * 60)
print()

# M4 = M4(trans) + (2(3n^2-5n+4)/(n(n-1)))*t3 + (2/C(n,4))*t5 + (4/C(n,4))*a2
# M4(trans) = fourth moment of Eulerian distribution
# Let me compute M4(trans) for several n

from math import factorial

for n in [5, 6, 7, 8]:
    # Eulerian numbers A(n,k) for k=0..n-1
    A = [0]*n
    A[0] = 1
    for step in range(2, n+1):
        new_A = [0]*n
        for k in range(n):
            new_A[k] = (k+1) * (A[k] if k < len(A) else 0) + (step - k) * (A[k-1] if k > 0 else 0)
        A = new_A

    total = factorial(n)
    M4_trans = sum(Fraction(k**4 * A[k], total) for k in range(n))

    dM4_dt3 = Fraction(2*(3*n**2 - 5*n + 4), n*(n-1))
    coeff_t5 = Fraction(2, comb(n, 4)) if n >= 5 else Fraction(0)
    coeff_a2 = Fraction(4, comb(n, 4)) if n >= 6 else Fraction(0)

    print(f"n={n}:")
    print(f"  M4(trans) = {M4_trans}")
    print(f"  dM4/dt3 = {dM4_dt3}")
    print(f"  coeff_t5 = {coeff_t5}")
    print(f"  coeff_a2 = {coeff_a2}")
    print(f"  M4 = {M4_trans} + ({dM4_dt3})*t3 + ({coeff_t5})*t5 + ({coeff_a2})*a2")
    print()
