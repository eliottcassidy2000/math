#!/usr/bin/env python3
"""
prime_via_derangements.py — kind-pasteur-2026-03-21-S16d

THE DERANGEMENT CONNECTION.

A000255(n) = number of permutations of [n+1] with no substring [k,k+1]
           = !{n+2}/(n+1)  where !k = subfactorial(k) = derangements of [k]
           = D(n) + D(n+1) where D(k) = subfactorial(k)

Our f(n, a) table counts permutations with EXACTLY a "unit ascents" and
0 "unit descents". The total compatible count (sum over a) = A000255(n-1).

E[H^2] = 2^{-2(n-1)} * sum_a f(n,a) * 2^a
       = 2^{-2(n-1)} * n! * sum_a count(a, d=0) * 2^a

The prime p(n) in CV^2(H) comes from the ARITHMETIC of this weighted sum.

KEY FORMULA:
  Var(H) = [weighted compatible count] / 4^{n-1} - (n!/2^{n-1})^2
         = [sum_a n! * c(n,a) * 2^a] / 4^{n-1} - n!^2/4^{n-1}
         = n!/4^{n-1} * [sum_a c(n,a) * 2^a - n!]

where c(n,a) = #{permutations with a unit ascents, 0 unit descents}.

So CV^2(H) = Var(H)/E[H]^2 = [sum c(n,a)*2^a - n!] / n!
           = [sum c(n,a)*2^a] / n! - 1

THE PRIME p(n) IS THE NUMERATOR OF [sum c(n,a)*2^a]/n! - 1.

Let me compute this exactly and factor it.
"""
from math import factorial, gcd
from fractions import Fraction
from itertools import permutations
from collections import Counter
import sys

sys.stdout.reconfigure(line_buffering=True)

def count_unit_ascents_descents(perm):
    asc = desc = 0
    for i in range(len(perm) - 1):
        diff = perm[i+1] - perm[i]
        if diff == 1: asc += 1
        elif diff == -1: desc += 1
    return asc, desc

def factor(n):
    if n <= 1: return {}
    f = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d,0)+1
            n //= d
        d += 1
    if n > 1: f[n] = f.get(n,0)+1
    return f

print("=" * 72)
print("  THE DERANGEMENT-PRIME CONNECTION")
print("  kind-pasteur-2026-03-21-S16d")
print("=" * 72)

# Compute c(n, a) and the weighted sum for each n
print("\n  EXACT COMPUTATION:")
print()

for n in range(3, 10):
    counts = Counter()
    for perm in permutations(range(n)):
        a, d = count_unit_ascents_descents(perm)
        if d == 0:
            counts[a] += 1

    # Weighted sum W = sum c(n,a) * 2^a
    W = sum(cnt * (2**a) for a, cnt in counts.items())

    # CV^2(H) = (W - n!) / n! = W/n! - 1
    nfact = factorial(n)
    CV2 = Fraction(W - nfact, nfact)
    CV2_num = CV2.numerator
    CV2_den = CV2.denominator

    # OCR from master formula
    R2 = Fraction(2*(n-2), n*(n-1)) / CV2 if CV2 > 0 else Fraction(1)
    one_minus_R2 = 1 - R2

    print(f"  n={n}:")
    print(f"    c(n,a): {dict(sorted(counts.items()))}")
    print(f"    W = sum c*2^a = {W}")
    print(f"    n! = {nfact}")
    print(f"    W - n! = {W - nfact}")
    print(f"    (W - n!)/n! = {CV2} = {float(CV2):.10f}")
    print(f"    CV^2(H) numerator = {CV2_num}, factors: {factor(CV2_num)}")
    print(f"    CV^2(H) denominator = {CV2_den}, factors: {factor(CV2_den)}")
    print(f"    R^2 = {R2} = {float(R2):.10f}")
    print(f"    1-R^2 = {one_minus_R2}")
    print()

# ============================================================
# THE STRUCTURE OF W - n!
# ============================================================

print("=" * 72)
print("  THE STRUCTURE OF W - n!")
print("=" * 72)

# W = sum_{a=0}^{n-1} c(n,a) * 2^a
# Note: c(n, n-1) = 1 (only identity permutation has n-1 unit ascents)
# c(n, n-2) = n-1 (remove one ascent = one "skip")

# W = sum c(n,a) * 2^a = 2^{n-1} + (n-1)*2^{n-2} + ...
# n! = product of 1..n

# The key arithmetic:
# W - n! is what produces the prime.

print()
for n in range(3, 10):
    counts = Counter()
    for perm in permutations(range(n)):
        a, d = count_unit_ascents_descents(perm)
        if d == 0:
            counts[a] += 1

    W = sum(cnt * (2**a) for a, cnt in counts.items())
    nfact = factorial(n)
    diff = W - nfact

    print(f"  n={n}: W={W}, n!={nfact}, W-n!={diff}")
    print(f"    W-n! factors: {factor(abs(diff)) if diff != 0 else 'zero'}")

    # The fraction (W-n!)/n!
    frac = Fraction(diff, nfact)
    print(f"    (W-n!)/n! = {frac}")
    print()

# ============================================================
# LOOK FOR PATTERN IN THE PRIMES
# ============================================================

print("=" * 72)
print("  THE PRIME SEQUENCE: LOOKING FOR PATTERNS")
print("=" * 72)

primes_list = [1, 1, 19, 13, 131, 131, 1097]
ns = [3, 4, 5, 6, 7, 8, 9]

print(f"\n  Prime sequence: {primes_list}")
print()

# Check: are these related to subfactorials?
# D(n) = !n = subfactorial
subfactorials = [1]  # D(0) = 1
for k in range(1, 15):
    subfactorials.append(k * subfactorials[-1] + (-1)**k)

print(f"  Subfactorials: {subfactorials[:12]}")
print(f"  !0=1, !1=0, !2=1, !3=2, !4=9, !5=44, !6=265, !7=1854, !8=14833, !9=133496, !10=1334961")

# A000255(n) = !(n+2)/(n+1)
print(f"\n  A000255 values:")
for k in range(12):
    val = subfactorials[k+2] // (k+1)
    print(f"    A000255({k}) = !{k+2}/({k+1}) = {subfactorials[k+2]}/{k+1} = {val}")

# Check: is the prime related to subfactorial residues?
print(f"\n  Subfactorial mod primes:")
for n, p in zip(ns, primes_list):
    if p > 1:
        d_n = subfactorials[n]
        d_n1 = subfactorials[n+1]
        d_n2 = subfactorials[n+2]
        print(f"    n={n}: p={p}, !{n}={d_n}, !{n} mod p = {d_n % p}, "
              f"!{n+1}={d_n1}, !{n+1} mod p = {d_n1 % p}")

# ============================================================
# THE TOPOLOGY: What does the prime COUNT?
# ============================================================

print("\n" + "=" * 72)
print("  THE TOPOLOGY: What is the geometric meaning of p(n)?")
print("=" * 72)

print("""
  The prime p(n) = numerator of CV^2(H) = numerator of (W - n!)/n!

  where W = sum_{a=0}^{n-1} c(n,a) * 2^a

  Since c(n,a) counts permutations with exactly a "unit ascents"
  and no "unit descents", the weighted sum W gives each such permutation
  a weight of 2^a. The identity permutation (a=n-1) gets weight 2^{n-1}.

  W/n! = E[2^a | d=0] * P(d=0)

  This is a MOMENT of the distribution of unit ascents among
  anti-succession-free permutations.

  The prime p(n) measures the IRREDUCIBLE PART of this moment.
  It's the part that cannot be decomposed into products of smaller
  combinatorial quantities.

  GEOMETRICALLY: Consider the "path intersection complex" K(n).
  Vertices = permutations (Hamiltonian paths).
  A 1-simplex connects sigma and tau if they are "compatible" (d=0).
  The weight 2^a on each edge measures the "overlap dimension."

  The Euler characteristic chi(K(n)) involves alternating sums
  of face counts. The prime p(n) may be related to a Betti number
  of this complex.
""")

# ============================================================
# THE OCR SEQUENCE: Full exact values
# ============================================================

print("=" * 72)
print("  COMPLETE EXACT OCR SEQUENCE (n=3..9)")
print("=" * 72)
print()
print("  n  | OCR = R^2    | 1 - R^2     | CV^2(H)     | prime p")
print("  ---+--------------+-------------+-------------+--------")

exact_R2 = {}
for n in range(3, 10):
    counts = Counter()
    for perm in permutations(range(n)):
        a, d = count_unit_ascents_descents(perm)
        if d == 0:
            counts[a] += 1

    W = sum(cnt * (2**a) for a, cnt in counts.items())
    nfact = factorial(n)
    CV2 = Fraction(W - nfact, nfact)

    if CV2 > 0:
        R2 = Fraction(2*(n-2), n*(n-1)) / CV2
    else:
        R2 = Fraction(1)

    one_minus = 1 - R2
    p_num = CV2.numerator if CV2.numerator > 0 else 1

    exact_R2[n] = R2
    print(f"  {n}  | {str(R2):12s} | {str(one_minus):11s} | {str(CV2):11s} | {p_num}")

print()
print("  THE PRIMES: 1, 1, 19, 13, 131, 131, 1097")
print("  Differences: _, _, _, -6, +118, 0, +966")
print("  Ratios: _, _, _, 0.68, 10.1, 1.0, 8.4")
print()
print("  1097 = 7*157 - 2 = ?  NO: 1097 is prime")
print("  1097 - 131 = 966 = 2*3*7*23")
print("  131 - 13 = 118 = 2*59")
print("  19 - 13 = 6 = 2*3")
print()

# Check: does p(n) grow like n^2?
print("  p(n) / n^2:")
for n, p in zip(ns[2:], primes_list[2:]):
    print(f"    n={n}: p={p}, p/n^2 = {p/n**2:.4f}")
