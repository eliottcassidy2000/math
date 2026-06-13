#!/usr/bin/env python3
"""
cumulant_bernoulli.py — Investigate the cumulant-Bernoulli structure.

For the transitive tournament, fwd = des(sigma) (descent count).
The cumulant generating function of des(sigma) is related to Bernoulli numbers.

Carlitz 1959: The CGF of des(sigma)/sqrt(Var) has cumulants:
  kappa_{2k} = (-1)^{k+1} (n+1) B_{2k} / (2k)  for k >= 1

where B_{2k} are Bernoulli numbers: B_2=1/6, B_4=-1/30, B_6=1/42, B_8=-1/30

So:
  kappa_2 = (n+1)/12  [from B_2 = 1/6, k=1: (n+1)*1/6 / 2 = (n+1)/12]
  kappa_4 = -(n+1)/120  [from B_4 = -1/30, k=2: -(n+1)*(-1/30)/4 = (n+1)/120... hmm sign]

Wait, let me be more careful. The standard formula is:
  kappa_r(des) = (-1)^{r-1} * (n+1) * B_r / r  for r >= 2 even (odd cumulants = 0 by symmetry)

Check: kappa_2 = (-1)^1 * (n+1) * B_2 / 2 = -(n+1)*(1/6)/2 = -(n+1)/12? No, should be positive.

Actually the sign convention depends on the reference. Let me just compute directly.

Author: opus-2026-03-07-S46d
"""
from fractions import Fraction
from math import factorial, comb
from itertools import permutations

def compute_cumulants_eulerian(n, max_k=8):
    """Compute cumulants of the Eulerian (descent) distribution directly."""
    # Get Eulerian numbers
    A = [0]*n
    A[0] = 1
    for step in range(2, n+1):
        new_A = [0]*n
        for k in range(n):
            new_A[k] = (k+1) * A[k] + (step - k) * (A[k-1] if k > 0 else 0)
        A = new_A

    total = factorial(n)
    mu = Fraction(n-1, 2)

    # Compute central moments
    moments = {}
    for r in range(max_k+1):
        moments[r] = sum(Fraction((k - mu)**r * A[k], total) for k in range(n))

    # Compute cumulants from moments
    # kappa_1 = mu_1 = 0 (centered)
    # kappa_2 = mu_2
    # kappa_3 = mu_3
    # kappa_4 = mu_4 - 3*mu_2^2
    # kappa_6 = mu_6 - 15*mu_4*mu_2 - 10*mu_3^2 + 30*mu_2^3
    # (all odd ones are 0 by symmetry)
    cumulants = {}
    cumulants[2] = moments[2]
    cumulants[4] = moments[4] - 3*moments[2]**2
    cumulants[6] = moments[6] - 15*moments[4]*moments[2] - 10*moments[3]**2 + 30*moments[2]**3
    cumulants[8] = (moments[8] - 28*moments[6]*moments[2] - 56*moments[5]*moments[3]
                    - 35*moments[4]**2 + 420*moments[4]*moments[2]**2
                    + 560*moments[3]**2*moments[2] - 630*moments[2]**4)

    return cumulants

print("EULERIAN CUMULANTS (transitive tournament)")
print("=" * 60)

# Bernoulli numbers
B = {0: Fraction(1), 1: Fraction(-1,2), 2: Fraction(1,6), 4: Fraction(-1,30),
     6: Fraction(1,42), 8: Fraction(-1,30), 10: Fraction(5,66)}

for n in range(3, 12):
    cums = compute_cumulants_eulerian(n)
    print(f"\nn={n}:")
    for r in [2, 4, 6, 8]:
        kr = cums[r]
        print(f"  kappa_{r} = {kr} = {float(kr):.8f}")

        # Test formula: kappa_{2k} = (-1)^{k+1} * (n+1) * B_{2k} / (2k)
        k = r // 2
        predicted = (-1)**(k+1) * (n+1) * B[r] / r
        print(f"    formula: (-1)^{k+1}*(n+1)*B_{r}/{r} = {predicted}, match = {kr == predicted}")

print("\n" + "=" * 60)
print("CUMULANT RATIOS (universal for transitive, n-independent)")
print("=" * 60)

for n in [5, 7, 9, 11]:
    cums = compute_cumulants_eulerian(n)
    k2 = cums[2]
    k4 = cums[4]
    k6 = cums[6]
    k8 = cums[8]

    r42 = k4/k2 if k2 != 0 else None
    r62 = k6/k2 if k2 != 0 else None
    r82 = k8/k2 if k2 != 0 else None

    print(f"\nn={n}:")
    print(f"  kappa_4/kappa_2 = {r42}")
    print(f"  kappa_6/kappa_2 = {r62}")
    print(f"  kappa_8/kappa_2 = {r82}")

# The ratios should be: B_{2k}/(2k) / (B_2/2) = B_{2k}/(2k * B_2/2) = 2*B_{2k}/(2k*B_2) = B_{2k}/(k*B_2)
print("\n\nPredicted ratios from Bernoulli:")
for r in [4, 6, 8]:
    k = r // 2
    ratio = (-1)**(k+1) * B[r] / r / ((-1)**(2) * B[2] / 2)
    # = (-1)^{k+1} * B_{2k}/(2k) / (-B_2/2)
    # = (-1)^{k} * 2*B_{2k} / (2k*B_2)
    ratio = (-1)**(k) * 2 * B[r] / (r * B[2])
    print(f"  kappa_{r}/kappa_2 = {ratio} = {float(ratio):.8f}")

print("\n" + "=" * 60)
print("GENERATING FUNCTION INTERPRETATION")
print("=" * 60)
print()
print("The CGF of des(sigma) for n-permutations, at parameter t, is:")
print("  K(t) = sum_{k>=1} kappa_k * t^k / k!")
print()
print("For the Eulerian distribution, the probability generating function is:")
print("  P(x) = A_n(x)/n! where A_n(x) = sum A(n,k) x^k")
print()
print("The CGF satisfies: K(t) = log(P(e^t))")
print()
print("For large n, the Eulerian distribution approaches N(mu, sigma^2),")
print("and the normalized cumulants kappa_r/sigma^r approach the Gaussian values")
print("(0 for r >= 3). The exact cumulants kappa_r = O(n) for all r,")
print("while sigma^r = O(n^{r/2}), so kappa_r/sigma^r = O(n^{1-r/2}) -> 0.")
print()

# KEY INSIGHT: For tournament T, the cumulant structure is:
# kappa_{2k}(T) = kappa_{2k}(transitive) + correction(T)
# The correction depends on cycle invariants at the (2k+1)-vertex level.
# kappa_2: correction = 4*t3/(n(n-1))
# kappa_4: correction = (2/C(n,4))*(t5 + 2*alpha_2) - 48*t3^2/(n(n-1))^2

print("=" * 60)
print("TOURNAMENT CUMULANT CORRECTIONS")
print("=" * 60)
print()
print("kappa_2(T) = kappa_2(trans) + 4*t3/(n(n-1))")
print("           = (n+1)*B_2 + 4*t3/(n(n-1))")
print()
print("kappa_4(T) = kappa_4(trans) + (2/C(n,4))*(t5 + 2*alpha_2) - 48*t3^2/(n(n-1))^2")
print("           = -(n+1)/120 + (2/C(n,4))*(t5 + 2*alpha_2) - 3*(4*t3/(n(n-1)))^2")
print()
print("Note: the t3^2 term = -3*[kappa_2 correction]^2 = -3*(delta kappa_2)^2")
print("This means: kappa_4(T) = kappa_4(trans) + (2/C(n,4))*(t5+2a2) - 3*(kappa_2(T)-kappa_2(trans))^2")
print("i.e., kappa_4 correction = new 5-vertex part - 3*(3-vertex correction)^2")
print()
print("This is the CUMULANT HIERARCHY:")
print("  Each kappa_{2k} adds genuinely new cycle structure (on 2k+1 vertices)")
print("  plus nonlinear corrections from lower levels.")
