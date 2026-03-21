#!/usr/bin/env python3
"""
ocr_denominator_analysis.py — kind-pasteur-2026-03-21-S16

The OCR denominators are: 1, 1, 19, 13, 131 for n=3,4,5,6,7.
All nondegenerate ones are PRIME: 19, 13, 131.

Let me investigate these as elements of Var(H)/Var(S2).

From exact computations:
  n=5: Var(H) = 285/16,  Var(S2) = 15/8,  R^2 = 18/19
       => Var(H)*Var(S2) = 285*15/(16*8) = 4275/128
       => Cov^2 = R^2 * Var(H)*Var(S2) = 18/19 * 4275/128 = 77400 - 4050/128
       Wait, let me just work with the formula directly.

  R^2 = Cov(S2,H)^2 / (Var(S2)*Var(H))
  1-R^2 = 1 - Cov^2/(Var_S2*Var_H) = (Var_S2*Var_H - Cov^2) / (Var_S2*Var_H)

The denominator of 1-R^2 (= the denominator of OCR) is related to Var(H).

KEY: Var(H) is a rational number. Let's find its denominator structure.
"""
from fractions import Fraction

# Exact values from computation
data = {
    3: {'sum_H': 12, 'sum_S2': 8, 'sum_H2': 28, 'sum_S22': 16, 'sum_HS2': 8, 'N': 8},
    4: {'sum_H': 192, 'sum_S2': 192, 'sum_H2': 768, 'sum_S22': 640, 'sum_HS2': 384, 'N': 64},
    5: {'sum_H': 7680, 'sum_S2': 5120, 'sum_H2': 75840, 'sum_S22': 33280, 'sum_HS2': 26880, 'N': 1024},
    6: {'sum_H': 737280, 'sum_S2': 229376, 'sum_H2': 21381120, 'sum_S22': 2097152, 'sum_HS2': 3686400, 'N': 32768},
    7: {'sum_H': 165150720, 'sum_S2': 22020096, 'sum_H2': 16386048000, 'sum_S22': 286261248, 'sum_HS2': 1321205760, 'N': 2097152},
}

print("=" * 72)
print("  OCR DENOMINATOR ANALYSIS")
print("=" * 72)

for n, d in data.items():
    N = d['N']
    E_H = Fraction(d['sum_H'], N)
    E_S2 = Fraction(d['sum_S2'], N)
    Var_H = Fraction(d['sum_H2'], N) - E_H**2
    Var_S2 = Fraction(d['sum_S22'], N) - E_S2**2
    Cov = Fraction(d['sum_HS2'], N) - E_H * E_S2

    if Var_S2 == 0 or Var_H == 0:
        R2 = Fraction(1)
    else:
        R2 = Cov**2 / (Var_S2 * Var_H)

    one_minus = 1 - R2

    print(f"\n  n={n}:")
    print(f"    E[H] = {E_H}")
    print(f"    E[S2] = {E_S2}")
    print(f"    Var(H) = {Var_H}")
    print(f"    Var(S2) = {Var_S2}")
    print(f"    Cov(H,S2) = {Cov}")
    print(f"    R^2 = {R2}")
    print(f"    1-R^2 = {one_minus}")

    # Key ratio: Var(H) / Cov(H,S2)^2 * Var(S2)
    # = 1/R^2

    # Let me also compute Var(H) * N^2 to remove fractions
    VH_int = d['sum_H2'] * N - d['sum_H']**2
    VS_int = d['sum_S22'] * N - d['sum_S2']**2
    Cov_int = d['sum_HS2'] * N - d['sum_H'] * d['sum_S2']

    print(f"    N^2*Var(H) = {VH_int}")
    print(f"    N^2*Var(S2) = {VS_int}")
    print(f"    N*Cov = {Cov_int}")

# Now let's look at patterns
print("\n" + "=" * 72)
print("  PATTERN SEARCH")
print("=" * 72)

# The denominators of 1-R^2 when expressed in lowest terms:
# n=5: 1/19. 19 is prime.
# n=6: 1/13. 13 is prime.
# n=7: 11/131. 131 is prime.

# The numerators of 1-R^2:
# n=5: 1
# n=6: 1
# n=7: 11

# What are 19, 13, 131?
# 19 = 20-1 = 4*5-1 = C(5,2)*2-1
# 13 = 15-2 = C(6,2)-2
# 131 = 132-1 = 12*11-1
# Or: 19 = 2^5-13? No, 32-13=19. Hmm.
# 13 = 2^4-3
# 131 = 2^7+3

# Another angle: denominators * (1-R^2 numerator) = 19*1, 13*1, 131*11
# = 19, 13, 1441
# 1441 = 11*131

# The R^2 numerators: 18, 12, 120
# 18 = 2*3^2, 12 = 2^2*3, 120 = 2^3*3*5
# Powers of 2 in numerator: 1, 2, 3 — increasing by 1!
# Factor of 3: 2, 1, 1 — decreasing
# 120 = 5! = 5*4*3*2*1. Interesting.
# 18 = 3!/3 * 6? No.
# 18 = C(9,2)? No, C(9,2)=36.
# 12 = C(4,2)*2
# 120 = C(10,3) = C(16,2)? No, C(16,2)=120! Yes!

# C(16,2) = 120. And 16 = C(7,2) + 1 = 22. No, 16 ≠ 22.
# Hmm. 120 = 5! = P(5,5). Also 120 = C(10,3).

# Actually: 120/131 at n=7.
# Var(H) at n=7 = 206325/128.
# 206325 = 3*5^2*1751? Let me factor.
print("\nKey quantities:")
for n in [5, 6, 7]:
    d = data[n]
    N = d['N']
    E_H = Fraction(d['sum_H'], N)
    Var_H = Fraction(d['sum_H2'], N) - E_H**2
    Var_S2 = Fraction(d['sum_S22'], N) - Fraction(d['sum_S2'], N)**2
    Cov = Fraction(d['sum_HS2'], N) - E_H * Fraction(d['sum_S2'], N)

    print(f"  n={n}: Var(H) = {Var_H}, Var(S2) = {Var_S2}, Cov = {Cov}")
    print(f"    Var(H) * 2^{2*n-2} = {Var_H * (1 << (2*(d['N'].bit_length()-1)))}")
    # Actually Var(H) = (sum_H2*N - sum_H^2) / N^2
    # = VH_int / N^2
    VH = d['sum_H2'] * N - d['sum_H']**2
    VS = d['sum_S22'] * N - d['sum_S2']**2
    CV = d['sum_HS2'] * N - d['sum_H'] * d['sum_S2']
    print(f"    N^2*Var(H) = {VH}")
    print(f"    N^2*Var(S2) = {VS}")
    print(f"    N*Cov = {CV}")

    # R^2 = CV^2 / (VS * VH)
    from math import gcd
    r2_num = CV**2
    r2_den = VS * VH
    g = gcd(r2_num, r2_den)
    p, q = r2_num // g, r2_den // g
    print(f"    R^2 = {p}/{q}")

    # Check: does Var(H) have n! in numerator?
    # E[H] = n!/2^{n-1}
    # E[H^2] involves permutation pair overlaps
    # Var(H) = E[H^2] - (n!/2^{n-1})^2

print("\n\n  POSSIBLE FORMULA FOR DENOMINATORS:")
print("  19  = ?")
print("  13  = ?")
print("  131 = ?")
print()
print("  Testing: denominator(n) = Var(H)*N^2 / gcd_factor")
print("  n=5: denom from VH=4354560, VS=2621440, CV=-3932160")
print("    R^2 = CV^2/(VS*VH) = 15461882265600/(11408506880000)")
from math import gcd
for n in [5, 6, 7]:
    d = data[n]
    N = d['N']
    VH = d['sum_H2'] * N - d['sum_H']**2
    VS = d['sum_S22'] * N - d['sum_S2']**2
    CV = d['sum_HS2'] * N - d['sum_H'] * d['sum_S2']
    r2_num = CV**2
    r2_den = VS * VH
    g = gcd(r2_num, r2_den)
    p, q = r2_num // g, r2_den // g
    print(f"  n={n}: R^2 = {p}/{q}, denom = {q}, num = {p}")
    # Factor the denominator
    q_temp = q
    factors = []
    for f in range(2, 200):
        while q_temp % f == 0:
            factors.append(f)
            q_temp //= f
    if q_temp > 1:
        factors.append(q_temp)
    print(f"    denom factors: {factors}")
