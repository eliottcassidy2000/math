#!/usr/bin/env python3
"""
kappa4_pattern.py — Analyze patterns in kappa_4 coefficients.

From computation:
  n=5: kappa_4 = -1/20 + (2/5)*t5 - (3/25)*t3^2
  n=6: kappa_4 = -7/120 + (2/15)*t5 + (4/15)*a2 - (4/75)*t3^2

Pattern analysis:
  Constant term: -1/20, -7/120. Let's check:
    -1/20 = -1/(4*5) = -1/(4n)
    -7/120 = -7/(20*6) = -7/(20n)
    Hmm. Or: -1/20, -7/120 = -6/120. Ratio = 6.
    Actually: for n=5, (n-1)/2 = 2. For n=6, (n-1)/2 = 5/2.
    -1/(n(n-1)) = -1/20, -1/30. No.
    Maybe: -1/20 and -7/120 don't have an obvious n-pattern with just 2 points.

  t5 coefficient: 2/5 = 2/n, 2/15 = 2/(n*(n-1)/2). Hmm.
    2/5 = 2/C(5,1)? No. 2/5 and 2/15 = 2/C(6,2)? C(6,2)=15 yes!
    Wait, at n=5: 2/5. At n=6: 2/15.
    2/5 = 2/5 and 2/15 = 2/15.
    ratio = 3. And 15/5 = 3.
    C(n,2)/n = (n-1)/2. 15/5 = 3 = (6-1)/2? No, 5/2.
    Actually 2/C(n,2): C(5,2)=10, 2/10=1/5 ≠ 2/5.
    Let me try: 2/(n!/(n-4)!/something...)
    Or simpler: the coefficient of t5 in M4 is 2/5 (n=5) and 2/15 (n=6).
    From the moment formula, M4 = ... + coeff_t5 * t5.
    kappa_4 = M4 - ... so the t5 coefficient in kappa_4 equals coeff in M4.
    (Because M3 doesn't depend on t5, and Var doesn't depend on t5.)
    So coeff_t5(kappa_4) = coeff_t5(M4).

    From efwd4_exact.py: 2/5 (n=5), 2/15 (n=6).
    2/5 = 2/(n) and 2/15 = 2/(n(n-1)/2)? 2/(6*5/2) = 2/15. Yes!
    But 2/5 ≠ 2/(5*4/2) = 2/10 = 1/5. So no.

    Try: 2/C(n,2)? C(5,2)=10 → 2/10=1/5. No.
    2*C(n,5-specific)?
    At n=5, 5-cycles use ALL vertices.
    At n=6, 5-cycles use 5 out of 6 vertices.
    Number of 5-vertex subsets: C(5,5)=1, C(6,5)=6.
    Normalization: 2/(5*C(5,5)) = 2/5. 2/(5*C(6,5)) = 2/30 ≠ 2/15.
    Hmm. 2/(C(5,5)*5) = 2/5. 2/(C(6,5)*5/2) = 2/15.
    Pattern: 2 / (C(n,5) * something)?
    Let me just note the values and move on to n=7 computation.

  t3^2 coefficient: -3/25 = -3/5^2, -4/75 = -4/(3*25) = -4/75.
    At n=5: -3/n^2 = -3/25. ✓
    At n=6: -4/75. -4/(n^2*?) = -4/36*something.
    Actually: -3/25 and -4/75.
    -3/(5^2) and -4/(75) = -4/(3*5^2).
    Hmm. Or: from Var formula, Var = (n+1)/12 - 4t3/(n(n-1)).
    Var^2 has t3^2 coefficient: (4/(n(n-1)))^2 = 16/(n^2(n-1)^2).
    kappa_4 has -3*Var^2, so t3^2 coeff from that = -3*16/(n^2(n-1)^2) = -48/(n^2(n-1)^2).
    n=5: -48/(25*16) = -48/400 = -3/25. ✓
    n=6: -48/(36*25) = -48/900 = -4/75. ✓

    So the t3^2 coefficient is EXACTLY -48/(n^2(n-1)^2) = -3*(4/(n(n-1)))^2.
    This comes purely from -3*Var^2, meaning mu_4 itself has NO t3^2 term!

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction

print("PATTERN IN KAPPA_4 COEFFICIENTS")
print("=" * 60)

# t3^2 coefficient check
for n in [5, 6]:
    coeff = Fraction(-48, n**2 * (n-1)**2)
    print(f"n={n}: -48/(n^2*(n-1)^2) = {coeff}")

print("\nn=5: -3/25 matches -48/400 =", Fraction(-48, 400))
print("n=6: -4/75 matches -48/900 =", Fraction(-48, 900))

# The t3^2 term comes ENTIRELY from -3*Var^2
# So mu_4 (central 4th moment) has NO quadratic terms in t3
# mu_4 is LINEAR in the invariants!

print("\n" + "=" * 60)
print("KEY INSIGHT: mu_4 is LINEAR in (t3, t5, alpha_2)")
print("=" * 60)

# mu_4 = kappa_4 + 3*Var^2
# For n=5: mu_4 = (-1/20 + 2t5/5 - 3t3^2/25) + 3*(1/2 - t3/5)^2
#         = -1/20 + 2t5/5 - 3t3^2/25 + 3/4 - 6t3/5 + 3t3^2/25
#         = 7/10 - 6t3/5 + 2t5/5

for n in [5, 6]:
    print(f"\nn={n}:")
    if n == 5:
        # kappa_4 = -1/20 + 2t5/5 - 3t3^2/25
        # Var = 1/2 - t3/5
        # Var^2 = 1/4 - t3/5 + t3^2/25
        # mu_4 = kappa_4 + 3*Var^2 = -1/20 + 2t5/5 - 3t3^2/25 + 3/4 - 3t3/5 + 3t3^2/25
        #       = (-1/20 + 3/4) + (-3/5)*t3 + (2/5)*t5 + 0*t3^2
        const = Fraction(-1,20) + Fraction(3,4)
        print(f"  mu_4 = {const} + ({Fraction(-3,5)})*t3 + ({Fraction(2,5)})*t5")
        print(f"  mu_4 = {const} - (3/5)*t3 + (2/5)*t5")
    elif n == 6:
        # kappa_4 = -7/120 + 2t5/15 + 4a2/15 - 4t3^2/75
        # Var = 7/12 - 2t3/15
        # Var^2 = 49/144 - 28t3/180 + 4t3^2/225 = 49/144 - 7t3/45 + 4t3^2/225
        # 3*Var^2 = 49/48 - 7t3/15 + 4t3^2/75
        # mu_4 = -7/120 + 2t5/15 + 4a2/15 - 4t3^2/75 + 49/48 - 7t3/15 + 4t3^2/75
        const = Fraction(-7,120) + Fraction(49,48)
        t3_coeff = Fraction(-7, 15)
        print(f"  mu_4 = {const} + ({t3_coeff})*t3 + ({Fraction(2,15)})*t5 + ({Fraction(4,15)})*a2")

print("\n" + "=" * 60)
print("PATTERN IN mu_4 COEFFICIENTS")
print("=" * 60)

# n=5: mu_4 = 7/10 - 3t3/5 + 2t5/5
# n=6: mu_4 = 77/80 - 7t3/15 + 2t5/15 + 4a2/15
# Wait let me recompute n=6 constant
c6 = Fraction(-7, 120) + Fraction(49, 48)
print(f"n=6 constant: -7/120 + 49/48 = {c6}")

# n=5: constant = 7/10 = 14/20
# n=6: constant = 77/80

# t3 coefficient: -3/5 = -12/20, -7/15 = -28/60
# ratio: (7/15)/(3/5) = 7/9
# Or: at n=5, coeff_t3 = -6/(n(n-1)) * n/... hmm
# Actually from mu_4 = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4
# The t3 dependence comes from M3 (slope 6/n) and M2 (slope -4/(n(n-1)))
# and M4 (slope = what for t3?)

# M4 t3 coefficient: 27/5 (n=5), 82/15 (n=6)
# M3 t3 coefficient: 6/5 (n=5), 1 (=6/6, n=6)  -- this is 6/n
# M2 t3 coefficient: -1/5 (n=5), -2/15 (n=6)  -- this is -4/(n(n-1))

# mu_4 = M4 - 4*mu*M3 + 6*mu^2*M2 - 3*mu^4
# d(mu_4)/d(t3) = d(M4)/d(t3) - 4*mu*d(M3)/d(t3) + 6*mu^2*d(M2)/d(t3)

for n_val in [5, 6]:
    mu = Fraction(n_val-1, 2)
    dM4_dt3 = Fraction(27, 5) if n_val == 5 else Fraction(82, 15)
    dM3_dt3 = Fraction(6, n_val)
    dM2_dt3 = Fraction(-4, n_val*(n_val-1))

    dmu4_dt3 = dM4_dt3 - 4*mu*dM3_dt3 + 6*mu**2*dM2_dt3
    print(f"\nn={n_val}:")
    print(f"  dM4/dt3 = {dM4_dt3}")
    print(f"  -4*mu*dM3/dt3 = {-4*mu*dM3_dt3}")
    print(f"  6*mu^2*dM2/dt3 = {6*mu**2*dM2_dt3}")
    print(f"  dmu_4/dt3 = {dmu4_dt3}")

    # But kappa_4 has zero t3 coefficient!
    # d(kappa_4)/dt3 = dmu4/dt3 - 6*Var*dVar/dt3
    dVar_dt3 = Fraction(-4, n_val*(n_val-1))
    Var_trans = Fraction(n_val+1, 12)  # Var at t3=0
    dk4_dt3_at_t3_0 = dmu4_dt3 - 6*Var_trans*dVar_dt3
    print(f"  dk4/dt3 at t3=0 = {dk4_dt3_at_t3_0}")
    print(f"  (should be 0 since kappa_4 has no linear t3 term)")

# So the LINEAR t3 coefficient vanishes in kappa_4.
# This means: dM4/dt3 - 4*mu*(6/n) + 6*mu^2*(-4/(n(n-1))) + 6*(n+1)/12 * 4/(n(n-1)) = 0
# i.e., dM4/dt3 = 4*mu*6/n - 6*mu^2*4/(n(n-1)) - 6*(n+1)/12 * 4/(n(n-1))
#                = 24*mu/n - 24*mu^2/(n(n-1)) - 2(n+1)/(n(n-1))

print("\n" + "=" * 60)
print("PREDICTING dM4/dt3 from zero linear kappa_4 constraint")
print("=" * 60)

for n_val in [5, 6, 7, 8, 9, 10]:
    mu = Fraction(n_val-1, 2)
    predicted = 24*mu/n_val - 24*mu**2/(n_val*(n_val-1)) - Fraction(2*(n_val+1), n_val*(n_val-1))
    print(f"n={n_val}: predicted dM4/dt3 = {predicted} = {float(predicted):.6f}")

# Check against known:
print("\nKnown: n=5: 27/5 =", float(Fraction(27,5)))
print("Known: n=6: 82/15 =", float(Fraction(82,15)))

for n_val in [5, 6]:
    mu = Fraction(n_val-1, 2)
    predicted = 24*mu/n_val - 24*mu**2/(n_val*(n_val-1)) - Fraction(2*(n_val+1), n_val*(n_val-1))
    known = Fraction(27, 5) if n_val == 5 else Fraction(82, 15)
    print(f"n={n_val}: predicted={predicted}, known={known}, match={predicted==known}")
