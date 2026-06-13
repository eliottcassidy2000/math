#!/usr/bin/env python3
"""
why_seven_ocr.py — kind-pasteur-2026-03-21-S16c

WHY DOES OCR STABILIZE AT n=7?

The exact OCR sequence:
  n: 3   4   5    6     7       8
  :  1   1   18/19 12/13 120/131 120/131

What changes at n=7 that makes it the transition?

STRUCTURAL THRESHOLDS AT EACH n:
- n=3: Only 3-cycles. H = 1 + 2c3. Trivially score-determined.
- n=4: Only 3-cycles. H = 1 + 2c3. Score-determined.
- n=5: 5-cycles appear! First non-score-determined contribution.
       alpha_2 = 0 (need 6 verts for disjoint pair). OCR drops to 18/19.
- n=6: alpha_2 > 0 first time (disjoint 3-3 pairs). OCR drops to 12/13.
- n=7: 7-cycles (Hamiltonian cycles) appear. OCR = 120/131.
- n=8: 5-3 disjoint pairs possible (3+5=8=n). OCR SAME = 120/131.

HYPOTHESIS: The stabilization at n=7 corresponds to a topological invariant.

n=7 is where:
1. ALL odd cycle lengths 3,5,7 up to n first coexist
2. The "full cycle spectrum" stabilizes
3. The dimension of the cycle space reaches its final growth regime

APPROACH: Compute the THEORETICAL R^2 by computing the three moments
E[S2^2], E[H*S2], E[H^2] from first principles as polynomials in n.
If these are polynomials, the ratio should simplify.

For E[S2]: we already know E[S2] = n(n-1)/4.
For Var(S2): computable from score covariances.
For Cov(H, S2): this is the key new quantity.
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb, factorial, gcd
from fractions import Fraction
import sys

sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  WHY DOES OCR STABILIZE AT n=7?")
print("  kind-pasteur-2026-03-21-S16c")
print("=" * 72)

# ============================================================
# PART 1: Theoretical Var(S2) as a function of n
# ============================================================

print("\n" + "=" * 72)
print("  PART 1: Var(S2) for random tournaments")
print("=" * 72)

# S2 = sum_i (s_i - (n-1)/2)^2
# Each s_i = sum_{j!=i} A[i][j], where A[i][j] ~ Bernoulli(1/2), independent.
# E[s_i] = (n-1)/2. Var(s_i) = (n-1)/4.
# X_i = (s_i - (n-1)/2)^2. E[X_i] = Var(s_i) = (n-1)/4.
# S2 = sum X_i. E[S2] = n(n-1)/4.

# For Var(S2) = sum Var(X_i) + 2 * sum_{i<j} Cov(X_i, X_j):

# Var(X_i) = E[X_i^2] - E[X_i]^2 = E[(s_i - mu)^4] - (Var(s_i))^2
# For Binomial(n-1, 1/2):
#   mu_4 = (n-1)(3n-7)/16  (fourth central moment)
#   Var(X_i) = mu_4 - sigma^4 = (n-1)(3n-7)/16 - (n-1)^2/16 = (n-1)(2n-6)/16 = (n-1)(n-3)/8

# For Cov(X_i, X_j): s_i and s_j share the arc (i,j).
# s_i = A[i][j] + sum_{k!=i,j} A[i][k]
# s_j = A[j][i] + sum_{k!=i,j} A[j][k] = (1-A[i][j]) + sum_{k!=i,j} A[j][k]
# The shared part: A[i][j] appears in s_i, (1-A[i][j]) appears in s_j.
# Cov(s_i, s_j) = Cov(A[i][j], 1-A[i][j]) = -Var(A[i][j]) = -1/4.
# So Cov(X_i, X_j) needs more careful computation.

# Let me compute Var(S2) directly for small n and look for a formula.

for n in range(3, 9):
    m = n*(n-1)//2
    total = 1 << m

    if total > 300000:
        # Sample
        np.random.seed(42)
        sample = 100000
        S2_vals = []
        for _ in range(sample):
            bits = np.random.randint(0, total)
            bits = int(bits)
            A = [[0]*n for _ in range(n)]
            pos = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1
            scores = [sum(A[i]) for i in range(n)]
            S2 = sum((s - (n-1)/2)**2 for s in scores)
            S2_vals.append(S2)
        E_S2 = np.mean(S2_vals)
        Var_S2 = np.var(S2_vals)
        method = "sampled"
    else:
        S2_vals = []
        for bits in range(total):
            A = [[0]*n for _ in range(n)]
            pos = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1
            scores = [sum(A[i]) for i in range(n)]
            S2 = sum((2*s - (n-1))**2 for s in scores)  # 4 * actual S2 to keep integer
            S2_vals.append(S2)
        E_S2 = Fraction(sum(S2_vals), total)
        Var_S2 = Fraction(sum(s**2 for s in S2_vals), total) - E_S2**2
        # Convert back: actual Var(S2) = Var(4*S2_actual)/16 = Var_S2/16
        # But for Var(S2) formula check, use the 4x-scaled version directly
        method = "exact"

    # Theoretical E[S2] = n(n-1)/4
    theo_E = Fraction(n*(n-1), 4)

    print(f"\n  n={n} ({method}):")
    print(f"    E[S2] = {float(E_S2 if isinstance(E_S2, Fraction) else E_S2):.6f} (theory: {float(theo_E):.6f})")
    print(f"    Var(S2) = {Var_S2 if isinstance(Var_S2, Fraction) else f'{Var_S2:.4f}'}")
    # Check formula: Var(S2) = n(n-1)(n-3)/8
    theo_var = Fraction(n*(n-1)*(n-3), 8) if n >= 3 else Fraction(0)
    print(f"    Theory n(n-1)(n-3)/8 = {float(theo_var):.6f}")
    if method == "exact":
        print(f"    Match: {Var_S2 == theo_var}")

# ============================================================
# PART 2: Exact moments E[H], Var(H), Cov(H, S2) as functions of n
# ============================================================

print("\n" + "=" * 72)
print("  PART 2: Exact moments as functions of n")
print("=" * 72)

# Known exact values:
exact_data = {
    3: {'E_H': Fraction(3,2), 'Var_H': Fraction(3,4),
        'Var_S2': Fraction(0), 'Cov': Fraction(0)},  # Var_S2=0 at n=3
    4: {'E_H': Fraction(3), 'Var_H': Fraction(3),
        'Var_S2': Fraction(1), 'Cov': Fraction(-3)},
    5: {'E_H': Fraction(15,2), 'Var_H': Fraction(285,16),
        'Var_S2': Fraction(15,2), 'Cov': Fraction(-45,4)},
    6: {'E_H': Fraction(45,2), 'Var_H': Fraction(585,4),
        'Var_S2': Fraction(15), 'Cov': Fraction(-45)},
    7: {'E_H': Fraction(315,4), 'Var_H': Fraction(206325,128),
        'Var_S2': Fraction(105,4), 'Cov': Fraction(-1575,8)},
    8: {'E_H': Fraction(315), 'Var_H': Fraction(371385,16),
        'Var_S2': None, 'Cov': None},  # Need to compute
}

# E[H] = n! / 2^{n-1}
print("\n  E[H] pattern:")
for n in range(3, 9):
    E_H = Fraction(factorial(n), 2**(n-1))
    print(f"    n={n}: E[H] = {E_H} = {float(E_H):.6f}")

# Var(H) pattern: let's look at Var(H) * 2^{something} / n!
print("\n  Var(H) scaled pattern:")
for n in [3, 4, 5, 6, 7, 8]:
    d = exact_data[n]
    VH = d['Var_H']
    EH = Fraction(factorial(n), 2**(n-1))
    ratio = VH / EH**2
    print(f"    n={n}: Var(H)/E[H]^2 = {ratio} = {float(ratio):.8f}")

# CV(H) = sqrt(Var(H)) / E[H]
print("\n  Coefficient of variation CV(H):")
for n in [3, 4, 5, 6, 7, 8]:
    d = exact_data[n]
    VH = d['Var_H']
    EH = Fraction(factorial(n), 2**(n-1))
    cv2 = VH / EH**2
    print(f"    n={n}: CV^2(H) = {cv2} = {float(cv2):.8f}, CV = {float(cv2)**0.5:.6f}")

# ============================================================
# PART 3: Cov(H, S2) pattern
# ============================================================

print("\n" + "=" * 72)
print("  PART 3: Cov(H, S2) pattern")
print("=" * 72)

print("\n  Cov(H, S2):")
for n in [4, 5, 6, 7]:
    d = exact_data[n]
    cov = d['Cov']
    print(f"    n={n}: Cov = {cov}")

# Ratio Cov / E[H]
print("\n  Cov(H,S2) / E[H]:")
for n in [4, 5, 6, 7]:
    d = exact_data[n]
    ratio = d['Cov'] / exact_data[n]['E_H']
    print(f"    n={n}: {ratio} = {float(ratio):.6f}")

# Ratio -Cov / (E[H] * E[S2])
print("\n  -Cov / (E[H] * n(n-1)/4):")
for n in [4, 5, 6, 7]:
    d = exact_data[n]
    ratio = -d['Cov'] / (d['E_H'] * Fraction(n*(n-1), 4))
    print(f"    n={n}: {ratio} = {float(ratio):.6f}")

# ============================================================
# PART 4: R^2 decomposition — what makes it 120/131?
# ============================================================

print("\n" + "=" * 72)
print("  PART 4: Decomposing R^2 = Cov^2 / (Var_S2 * Var_H)")
print("=" * 72)

for n in [5, 6, 7]:
    d = exact_data[n]
    r2 = d['Cov']**2 / (d['Var_S2'] * d['Var_H'])
    print(f"\n  n={n}: R^2 = {r2}")
    print(f"    Cov^2 = {d['Cov']**2}")
    print(f"    Var_S2 * Var_H = {d['Var_S2'] * d['Var_H']}")

    # The denominator minus numerator
    diff = d['Var_S2'] * d['Var_H'] - d['Cov']**2
    print(f"    Denominator - Numerator = {diff}")
    print(f"    This is Var_S2 * Var_eps_conditional")

# ============================================================
# PART 5: The structural meaning of 131
# ============================================================

print("\n" + "=" * 72)
print("  PART 5: What IS 131?")
print("=" * 72)

# Var(H) numerators: 3, 3, 285, 585, 206325, 371385
# Factored:
# n=3: 3
# n=4: 3
# n=5: 3 * 5 * 19
# n=6: 3^2 * 5 * 13
# n=7: 3^2 * 5^2 * 7 * 131
# n=8: 3^4 * 5 * 7 * 131

# E[H^2] numerators: 3, 12, 1185, 1305, 1000125, ...
# At n=8: E[H^2] = sum_H2 / N = 32866314485760 / 268435456
EH2_8 = Fraction(32866314485760, 268435456)
print(f"  E[H^2](n=8) = {EH2_8}")
VH_8 = EH2_8 - Fraction(315)**2
print(f"  Var(H)(n=8) = {VH_8} = {VH_8.numerator}/{VH_8.denominator}")

# Factor 371385
n_val = 371385
factors = []
temp = n_val
for f in range(2, 200):
    while temp % f == 0:
        factors.append(f)
        temp //= f
if temp > 1:
    factors.append(temp)
print(f"\n  371385 = {' * '.join(map(str, factors))}")

# What is 131 in terms of tournament parameters?
print("\n  131 in various expressions:")
print(f"    2^7 + 3 = {2**7 + 3}")
print(f"    2^7 + 2 + 1 = {2**7 + 3}")
print(f"    C(7,3)*4 - 5 = {comb(7,3)*4 - 5}")  # 140-5=135, no
print(f"    C(7,2)*6 + 5 = {comb(7,2)*6 + 5}")   # 126+5=131!
print(f"    *** C(7,2)*6 + 5 = 21*6 + 5 = 131 ***")
print()

# Check: does C(n,2)*6 + 5 give the denominator for other n?
print("  Testing C(n,2)*6 + 5:")
for n_test in [5, 6, 7, 8]:
    val = comb(n_test, 2) * 6 + 5
    print(f"    n={n_test}: C(n,2)*6+5 = {val}", end="")
    if n_test == 5:
        print(f"  (actual denom: 19, {val==19})")
    elif n_test == 6:
        print(f"  (actual denom: 13, {val==13})")
    elif n_test >= 7:
        print(f"  (actual denom: 131, {val==131})")

# Try other formulas
print("\n  Other attempts:")
# n=5: 19 = ?. n=6: 13 = ?. n=7,8: 131 = ?
# These are all different, so it's not a single formula of n.
# But maybe it's a formula of the Var(H) components.

# Let's look at Var(H) = A/B where B is a power of 2.
# Var(H) * B gives an integer. That integer contains the prime.
# The prime comes from E[H^2] * B - (E[H])^2 * B = integer.

# E[H^2] involves the permutation pair overlap structure.
# The prime comes from the IRREDUCIBLE part of this structure.

print("\n  KEY OBSERVATION:")
print("  131 = 6*C(7,2) + 5 = 6*21 + 5")
print("  But also: 19 ≠ 6*C(5,2) + 5 = 65")
print("  So this is NOT a universal formula.")
print()
print("  The denominators 19, 13, 131 may not come from a simple formula of n.")
print("  They come from the ARITHMETIC of the permutation pair overlap structure,")
print("  which depends on n in a complex way.")
print()
print("  The key structural question: why does the SAME prime (131)")
print("  appear in Var(H) at both n=7 and n=8?")
print()
print("  HYPOTHESIS: At n >= 7, the permutation pair overlap distribution")
print("  f(n, a) has a STABLE shape, and the resulting Var(H) always")
print("  contains 131 as its largest prime factor.")
print()
print("  This would mean OCR = 120/131 for ALL n >= 7.")
print("  The n=9, 10, 11 samplers will test this.")
