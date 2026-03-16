#!/usr/bin/env python3
"""
AMPLITUDE GLITCH ANALYSIS
===========================
The n=9 row of the THM-080 amplitude table is wrong.
Let's figure out WHERE the wrong values came from.

opus-2026-03-16-S73
"""
from math import factorial
from fractions import Fraction

print("=" * 70)
print("GLITCH ANALYSIS: Where did the wrong n=9 table values come from?")
print("=" * 70)

# The wrong table values for n=9:
wrong = {
    (1, 0): Fraction(3, 2),
    (3, 0): Fraction(3, 8),
    (3, 1): Fraction(3, 4),
    (5, 0): Fraction(1, 16),
    (7, 0): Fraction(1, 128),
}

# The correct formula values:
def correct_M(n, d, s):
    return Fraction(2**s * factorial(n - 2 - d), 2**(n-2))

# The H formula values:
def correct_H(n, d, r):
    return Fraction(2**r * factorial(n - d), 2**(n-1))

print("\nHypothesis 1: Table values are H amplitudes at degree d+1?")
for (d, s), tv in sorted(wrong.items()):
    for r in range(5):
        hv = correct_H(9, d+1, r)
        if hv == tv:
            print(f"  MATCH: H(n=9, d={d+1}, r={r}) = {hv} == table(d={d},s={s})")
            break
    else:
        print(f"  No H match for table(d={d},s={s})={tv}")

print("\nHypothesis 2: Someone shifted n by 2 (used n=7 formula for n=9)?")
for (d, s), tv in sorted(wrong.items()):
    if 7 - 2 - d >= 0:
        mv = correct_M(7, d, s)
        ratio = tv / mv if mv != 0 else None
        print(f"  d={d}, s={s}: M(n=7)={mv}, table(n=9)={tv}, ratio={ratio}")

print("\nHypothesis 3: Table uses (n-d)!/2^{n-1} instead of (n-2-d)!/2^{n-2}?")
for (d, s), tv in sorted(wrong.items()):
    alt = Fraction(2**s * factorial(9 - d), 2**8)
    print(f"  d={d}, s={s}: alt_formula={alt}, table={tv}, match={alt==tv}")

print("\nHypothesis 4: Table entries are M values with wrong s count?")
for (d, s), tv in sorted(wrong.items()):
    for s_try in range(5):
        mv = correct_M(9, d, s_try)
        if mv == tv:
            print(f"  MATCH: M(n=9, d={d}, s={s_try}) = {mv} == table(d={d},s={s})")
            break
    else:
        print(f"  No s-match for table(d={d},s={s})={tv}")

print("\nHypothesis 5: DIRECT comparison — what s value produces each table entry?")
print("  I.e., solve 2^s * (n-2-d)!/2^{n-2} = table_value for s")
for (d, s_claimed), tv in sorted(wrong.items()):
    base = Fraction(factorial(9 - 2 - d), 2**7)
    if base > 0:
        ratio = tv / base
        # Is ratio a power of 2?
        r = ratio
        s_actual = 0
        while r > 1:
            r /= 2
            s_actual += 1
        if r == 1:
            print(f"  d={d}, s_claimed={s_claimed}: table={tv}, base={base}, "
                  f"table/base = 2^{s_actual}  (so s_actual={s_actual})")
        else:
            print(f"  d={d}, s_claimed={s_claimed}: table={tv}, base={base}, "
                  f"table/base = {ratio} (NOT a power of 2!)")

print("\n" + "=" * 70)
print("CORRECT n=9 amplitude table:")
print("=" * 70)
print(f"\n| d  | s=0       | s=1       | s=2       | s=3       |")
print(f"|----|-----------|-----------|-----------|-----------| ")
for d in [1, 3, 5, 7]:
    row = f"| {d}  |"
    for s in range(4):
        if 9 - 2 - d >= 0:
            val = correct_M(9, d, s)
            row += f" {str(val):>9} |"
        else:
            row += f"     —    |"
    print(row)

print("\n" + "=" * 70)
print("THE DIMENSION LADDER (complete)")
print("=" * 70)
print()
print("H amplitude at even degree d with r components:")
print("  H_amp(n,d,r) = 2^r * (n-d)! / 2^{n-1}")
print()
print("M amplitude at odd degree d with s unrooted components:")
print("  M_amp(n,d,s) = 2^s * (n-2-d)! / 2^{n-2}")
print()
print("Ladder ratio: H_amp(n,d,r) / M_amp(n,d-1,s=r-1) = n-d")
print()
print("This means: going from M at degree d-1 to H at degree d, the amplitude")
print("increases by factor (n-d). This is the SAME factor as going from")
print("(n-2-(d-1))! to (n-d)! in the numerator:")
print("  (n-d)! / (n-d-1)! = n-d")
print("combined with the normalization shift 2^{n-1}/2^{n-2} = 2 and")
print("the component count shift 2^r/2^{r-1} = 2, giving:")
print("  ratio = (n-d) * (2^r / 2^{r-1}) * (2^{n-2} / 2^{n-1}) = (n-d) * 2 * 1/2 = n-d")
print()

print("Ladder descent (complete table for odd n=3..11):")
print()
for n in [3, 5, 7, 9, 11]:
    print(f"  n={n}:")
    max_d = n - 1
    for d_H in range(2, max_d + 1, 2):
        d_M = d_H - 1
        ratio = n - d_H
        H_amp = Fraction(2 * factorial(n - d_H), 2**(n-1))  # r=1
        M_amp = Fraction(factorial(n - 2 - d_M), 2**(n-2))  # s=0
        print(f"    H(d={d_H},r=1)={H_amp} → M(d={d_M},s=0)={M_amp}  [÷{ratio}]", end="")
        if d_M >= 2:
            d_M2 = d_M - 1
            H2_amp = Fraction(2 * factorial(n - d_M), 2**(n-1))
            print(f"  ← H(d={d_M-1},r=1)={H2_amp} [÷{n-d_M+1}]", end="")
        print()
    print()

print("=" * 70)
print("BIT ACCUMULATION in the Walsh spectrum")
print("=" * 70)
print()
print("Each tournament T has M[a,b] expanded in Walsh basis:")
print("  M[a,b](T) = sum_S hat{M}[S] * chi_S(T)")
print()
print("The 2^s factor means each unrooted component contributes ONE BIT")
print("to the amplitude. This parallels OCF: each independent odd cycle")
print("contributes one bit (factor 2) to I(Omega,2).")
print()
print("The dimension ladder ratio (n-d) is always ODD for odd n, even d:")
print("  v_2(n-d) = 0 when n odd and d even")
print("This means the H/M step NEVER introduces new factors of 2.")
print("ALL powers of 2 come from the 2^r (resp 2^s) component count.")
print()
print("Concretely: the 2-adic valuation of the amplitude is:")
print("  v_2(H_amp) = r - (n-1) + v_2((n-d)!)")
print("  v_2(M_amp) = s - (n-2) + v_2((n-2-d)!)")
print()

for n in [5, 7, 9]:
    print(f"  n={n}:")
    for d in range(1, n-1):
        fact = factorial(n - 2 - d) if d <= n-2 else 0
        if fact > 0:
            # count trailing zeros in fact
            v2_fact = 0
            f = fact
            while f % 2 == 0:
                v2_fact += 1
                f //= 2
            v2_base = v2_fact - (n - 2)
            print(f"    d={d}: (n-2-d)!={fact}, v_2((n-2-d)!)={v2_fact}, "
                  f"v_2(base_amp)={v2_fact}-(n-2)={v2_base}")

print("\nDone.")
