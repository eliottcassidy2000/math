#!/usr/bin/env python3
"""
kappa4_general_n.py — Derive general-n kappa_4 formula.

KEY FINDINGS from n=5,6:
  kappa_4 = C(n) + 0*t3 + B(n)*t5 + D(n)*alpha_2 - 48/(n^2(n-1)^2)*t3^2

where:
  - t3^2 coefficient = -48/(n(n-1))^2 = -3*(4/(n(n-1)))^2  [from -3*Var^2]
  - t5 coefficient: 2/5 (n=5), 2/15 (n=6)
  - alpha_2 coefficient: 0 (n=5), 4/15 (n=6)
  - constant: -1/20 (n=5), -7/120 (n=6)

Let me understand each coefficient.

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction

# mu_4 = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4
# mu_4 is LINEAR in (t3, t5, alpha_2):
#   n=5: mu_4 = 7/10 + (3/5)*t3 + (2/5)*t5
#   n=6: mu_4 = 77/80 + (7/15)*t3 + (2/15)*t5 + (4/15)*alpha_2

# The t3 coefficient in mu_4:
# d(mu_4)/dt3 = d(M4)/dt3 - 4*mu*d(M3)/dt3 + 6*mu^2*d(M2)/dt3
# = d(M4)/dt3 - 4*(n-1)/2*(6/n) + 6*(n-1)^2/4*(4/(n(n-1)))
# = d(M4)/dt3 - 12(n-1)/n + 6(n-1)/n
# = d(M4)/dt3 - 6(n-1)/n

print("COEFFICIENTS OF mu_4")
print("=" * 60)

for n in [5, 6]:
    mu = Fraction(n-1, 2)

    # Known M4 t3 coefficient
    if n == 5:
        dM4_dt3 = Fraction(27, 5)
    else:
        dM4_dt3 = Fraction(82, 15)

    dM3_dt3 = Fraction(6, n)
    dM2_dt3 = Fraction(4, n*(n-1))

    dmu4_dt3 = dM4_dt3 - 4*mu*dM3_dt3 + 6*mu**2*dM2_dt3
    print(f"n={n}: dmu4/dt3 = {dM4_dt3} - {4*mu*dM3_dt3} + {6*mu**2*dM2_dt3} = {dmu4_dt3}")

    # Also: -6(n-1)/n shortcut
    shortcut = dM4_dt3 - Fraction(6*(n-1), n)
    print(f"  shortcut: dM4/dt3 - 6(n-1)/n = {dM4_dt3} - {Fraction(6*(n-1),n)} = {shortcut}")
    print(f"  match: {shortcut == dmu4_dt3}")

# Hmm, let me recompute more carefully
print("\n\nCARFUL RECOMPUTATION")
for n in [5, 6]:
    mu = Fraction(n-1, 2)

    if n == 5:
        dM4_dt3 = Fraction(27, 5)
    else:
        dM4_dt3 = Fraction(82, 15)

    term1 = dM4_dt3
    term2 = -4*mu*Fraction(6, n)
    term3 = 6*mu**2 * Fraction(4, n*(n-1))

    print(f"\nn={n}:")
    print(f"  dM4/dt3 = {term1}")
    print(f"  -4*mu*dM3/dt3 = {term2}")
    print(f"  6*mu^2*dM2/dt3 = {term3}")
    print(f"  sum = {term1 + term2 + term3}")

# Now for kappa_4 = mu_4 - 3*Var^2:
# d(kappa_4)/dt3 = d(mu_4)/dt3 - 6*Var*d(Var)/dt3
# At t3=0: d(kappa_4)/dt3 = d(mu_4)/dt3 - 6*(n+1)/12 * 4/(n(n-1))
#                          = d(mu_4)/dt3 - 2(n+1)/(n(n-1))

print("\n\nKAPPA_4 LINEAR T3 CHECK")
print("=" * 60)
for n in [5, 6]:
    mu = Fraction(n-1, 2)

    if n == 5:
        dmu4_dt3 = Fraction(3, 5)
    else:
        dmu4_dt3 = Fraction(7, 15)

    dVar_dt3 = Fraction(4, n*(n-1))
    Var_0 = Fraction(n+1, 12)

    dk4_dt3 = dmu4_dt3 - 6*Var_0*dVar_dt3
    print(f"n={n}: dmu4/dt3 = {dmu4_dt3}, 6*Var(0)*dVar/dt3 = {6*Var_0*dVar_dt3}")
    print(f"  dk4/dt3 at t3=0 = {dk4_dt3}")

# aha! This IS zero!
# dmu4/dt3 = 6*Var(0)*dVar/dt3 exactly!
# n=5: 3/5 = 6*(1/2)*(1/5) = 6/10 = 3/5. ✓
# n=6: 7/15 = 6*(7/12)*(2/15) = 6*14/180 = 84/180 = 7/15. ✓

# So: dmu4/dt3 = 6*((n+1)/12)*(4/(n(n-1))) = 2(n+1)/(n(n-1))
# And: dM4/dt3 = dmu4/dt3 + 4*mu*6/n - 6*mu^2*4/(n(n-1))
#              = 2(n+1)/(n(n-1)) + 12(n-1)/n - 6(n-1)/n
#              = 2(n+1)/(n(n-1)) + 6(n-1)/n

print("\n\nPREDICTED dM4/dt3 FOR GENERAL n")
print("=" * 60)
for n in [5, 6, 7, 8, 9]:
    mu = Fraction(n-1, 2)
    predicted = Fraction(2*(n+1), n*(n-1)) + Fraction(6*(n-1), n)
    print(f"n={n}: dM4/dt3 = 2(n+1)/(n(n-1)) + 6(n-1)/n = {Fraction(2*(n+1), n*(n-1))} + {Fraction(6*(n-1), n)} = {predicted}")

    if n == 5:
        known = Fraction(27, 5)
        print(f"  known = {known}, match = {predicted == known}")
    elif n == 6:
        known = Fraction(82, 15)
        print(f"  known = {known}, match = {predicted == known}")

# Simplify: 2(n+1)/(n(n-1)) + 6(n-1)/n = [2(n+1) + 6(n-1)^2] / (n(n-1))
# = [2n+2 + 6n^2-12n+6] / (n(n-1))
# = [6n^2 - 10n + 8] / (n(n-1))
# = 2[3n^2 - 5n + 4] / (n(n-1))

print("\nSimplified: dM4/dt3 = 2(3n^2 - 5n + 4) / (n(n-1))")
for n in [5, 6, 7, 8, 9]:
    val = Fraction(2*(3*n**2 - 5*n + 4), n*(n-1))
    print(f"  n={n}: {val} = {float(val):.6f}")

# Now: the kappa_4 constant term and the t5, alpha_2 coefficients.
# Constant term: kappa_4(transitive) = mu_4(trans) - 3*Var(trans)^2
# mu_4(trans) = M4(trans) - 4*mu*M3(trans) + 6*mu^2*M2(trans) - 3*mu^4
# For transitive: M_r = (1/n!) sum_sigma des(sigma)^r = r-th moment of Eulerian distribution
# These are universal polynomials in n.

# Actually, kappa_4 of the Eulerian distribution (= transitive tournament) is known.
# It's: kappa_4 = -(n+1)/120 (this is a classical result)
# Wait let me check: n=5: should be -6/120 = -1/20. ✓
# n=6: should be -7/120. ✓!

print("\n\nKAPPA_4 CONSTANT TERM")
print("=" * 60)
for n in [5, 6, 7, 8, 9, 10]:
    predicted_const = Fraction(-(n+1), 120)
    print(f"  n={n}: kappa_4(transitive) = -(n+1)/120 = {predicted_const}")

# ✓ This is the fourth cumulant of the Eulerian distribution!
# Classical result: kappa_r of des(sigma) = (Bernoulli numbers related)
# kappa_2 = (n+1)/12, kappa_4 = -(n+1)/120
# These are proportional to Bernoulli numbers: B_2 = 1/6, B_4 = -1/30

print("\nNote: kappa_2 = (n+1)/12, kappa_4 = -(n+1)/120")
print("Ratio: kappa_4/kappa_2 = (-1/120)/(1/12) = -1/10")
print("This ratio is INDEPENDENT of n!")

# t5 coefficient: 2/5 (n=5), 2/15 (n=6).
# These equal the t5 coefficient in M4 exactly, because M3, M2 don't depend on t5.
# Pattern: 2/5 = 2/C(5,2)? No, C(5,2)=10.
# 2/5 and 2/15. Ratio = 3. C(6,5)/C(5,5) = 6.
# Actually: 2/(n-1)! * something?
# 2/5 = 2/n, 2/15 = 2/(n(n-1)/2) = 4/(n(n-1))?
# No: 4/(5*4)=1/5 ≠ 2/5. 4/(6*5)=2/15. ✓ for n=6 but not n=5.
# Maybe: at n=5, t5 uses ALL 5 vertices, so the normalization is different.
# At n=5: there's C(5,5)=1 5-vertex subset, each contributing 2 directed 5-cycles (CW, CCW).
# At n=6: there are C(6,5)=6 5-vertex subsets.

# The coefficient of t5 in E[fwd^4] was found to be 2/5 and 2/15.
# These were computed by Cramer's rule fitting, not derived from first principles.
# To predict the general formula, we'd need to understand what produces the t5 dependence.

# alpha_2 coefficient: 0 (n=5, impossible to have disjoint 3-cycles on 5 vertices),
# 4/15 (n=6).

print("\n\nSUMMARY: KAPPA_4 AT GENERAL n")
print("=" * 60)
print("kappa_4 = -(n+1)/120 + coeff_t5(n)*t5 + coeff_a2(n)*alpha_2 - 48/(n(n-1))^2 * t3^2")
print("\nKnown coefficients:")
print(f"  n=5: coeff_t5 = 2/5, coeff_a2 = 0 (no disjoint pair possible)")
print(f"  n=6: coeff_t5 = 2/15, coeff_a2 = 4/15")
print(f"  t3^2: -48/(n(n-1))^2 for all n")
print(f"  constant: -(n+1)/120 for all n")
print(f"\nPrediction: coeff_t5(n) and coeff_a2(n) likely have clean formulas in n.")

# Let me check the mu_4 formula more carefully.
# mu_4 = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4
# The t3 part:
# dmu4/dt3 = dM4/dt3 - 4*mu*6/n + 6*mu^2*4/(n(n-1))
# We showed dmu4/dt3 = 2(n+1)/(n(n-1)) for all n.
# This gives: dM4/dt3 = 2(n+1)/(n(n-1)) + 12(n-1)/n - 6(n-1)/n
#                      = 2(n+1)/(n(n-1)) + 6(n-1)/n
# = 2(3n^2-5n+4)/(n(n-1))

# At n=7: dM4/dt3 = 2(147-35+4)/(7*6) = 2*116/42 = 232/42 = 116/21
print(f"\nPredicted E[fwd^4] t3 coefficient at n=7: {Fraction(2*(3*49-35+4), 42)} = {float(Fraction(232,42)):.6f}")

# Can verify at n=7 using the worpitzky_n7_tiny.py data
