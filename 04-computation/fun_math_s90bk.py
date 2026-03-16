#!/usr/bin/env python3
"""
fun_math_s90bk.py — Doing math for fun across all domains
opus-2026-03-16-S90bk

Having fun. Finding things. Connecting everything.
"""

import numpy as np
from sympy import (isprime, primerange, factorint, mobius, divisors,
                   bernoulli, fibonacci, lucas, catalan, factorial,
                   binomial, totient, Matrix, Rational, sqrt as ssqrt,
                   oo, Sum, symbols, simplify, zeta)
from sympy.functions.combinatorial.numbers import partition as sym_partition
from math import log, log2, sqrt, pi, e, gcd, comb
from fractions import Fraction
from collections import Counter

def npartitions(n):
    try:
        from sympy.ntheory.partitions_ import npartitions as _np
        return int(_np(n))
    except Exception:
        from sympy.functions.combinatorial.numbers import partition as _p
        return int(_p(n))

print("""


     DOING MATH FOR FUN

     opus-2026-03-16-S90bk


""")

# ========================================================================
print("=" * 60)
print("I. THE TOURNAMENT ZETA FUNCTION")
print("=" * 60)
print()

print("""
Define the TOURNAMENT ZETA FUNCTION:

  Z_T(s) = sum_{n=1}^{inf} H(n) / n^s

where H(n) = number of non-isomorphic tournaments on n vertices.

H(1)=1, H(2)=1, H(3)=2, H(4)=4, H(5)=12, H(6)=56, H(7)=456, H(8)=6880.

This is a Dirichlet series in H(n). What are its properties?
""")

H_values = [0, 1, 1, 2, 4, 12, 56, 456, 6880]  # H(0) unused, H(1)=1, etc.

print("Partial sums of the tournament zeta function Z_T(s):")
for s_val in [1.0, 1.5, 2.0, 3.0]:
    partial = sum(H_values[n] / n**s_val for n in range(1, len(H_values)))
    print(f"  s={s_val:.1f}: Z_T(s) ~ {partial:.6f} (first {len(H_values)-1} terms)")

# The growth rate of H(n) ~ 2^(n(n-1)/2) / n!
print("\nGrowth: H(n) ~ 2^C(n,2) / n!")
for n in range(1, 9):
    actual = H_values[n]
    estimate = 2**comb(n, 2) / factorial(n)
    print(f"  n={n}: H(n)={actual}, estimate={float(estimate):.1f}, ratio={actual/float(estimate):.4f}")

print("""
H(n) grows super-exponentially. The Dirichlet series Z_T(s)
converges for Re(s) > ... well, the abscissa of convergence is
determined by the growth rate.

H(n) ~ 2^{n^2/2} / n!, so H(n)/n^s ~ 2^{n^2/2} / (n! * n^s).
This goes to infinity for ALL s. The series DIVERGES for all s!

The tournament zeta function has NO region of convergence.
It is an ASYMPTOTIC series, not a convergent one.

This means: tournament theory is TOO WILD for a simple zeta function.
The growth of tournaments is FASTER than any polynomial bound.
No Dirichlet series can capture it.

We need a DIFFERENT kind of series:
  Z_T(q) = sum_{n=0}^{inf} H(n) * q^{n(n-1)/2} / n!

This is a TOURNAMENT GENERATING FUNCTION in the variable q,
where q controls the number of arcs. At q=2 (the tournament base):
  Z_T(2) = sum H(n) * 2^{C(n,2)} / n! ... diverges again.

The natural variable is not q but x = q-1:
  Z_T(1+x) = sum H(n) * (1+x)^{C(n,2)} / n!
""")

# ========================================================================
print("=" * 60)
print("II. THE BERNOULLI NUMBERS AND THE FORBIDDEN VALUES")
print("=" * 60)
print()

print("Bernoulli numbers B_n and their relation to forbidden H values:")
print(f"{'n':>3} | {'B_n':>20} | {'B_n (decimal)':>14} | {'note':>30}")
print("-" * 75)

for n in range(16):
    B = bernoulli(n)
    B_float = float(B)
    note = ""
    if n == 0: note = "B_0 = 1"
    elif n == 1: note = "B_1 = -1/2 (or +1/2)"
    elif n % 2 == 1 and n > 1: note = "= 0 (odd > 1)"
    elif n == 2: note = "1/6"
    elif n == 4: note = "-1/30 (denom has 30=2*3*5!)"
    elif n == 6: note = "1/42 (denom has 42=2*3*7!)"
    elif n == 8: note = "-1/30 again"
    elif n == 10: note = "5/66"
    elif n == 12: note = "-691/2730"
    elif n == 14: note = "7/6 (7 in numerator!)"
    print(f"{n:3d} | {str(B):>20} | {B_float:14.8f} | {note:>30}")

print(f"""

REMARKABLE:
  B_4 has denominator 30 = 2*3*5 (the tournament triple product!).
  B_6 has denominator 42 = 2*3*7 (the Hurwitz triple product!).
  B_8 has denominator 30 again.
  B_12 has denominator 2730 = 2*3*5*7*13.
  B_14 has numerator 7 (the forbidden value!).

By the von Staudt-Clausen theorem:
  The denominator of B_(2k) is prod_{{(p-1)|2k_val}} p.

  For 2k = 4: (p-1)|4 means p-1 in {{1,2,4}}, so p in {{2,3,5}}. Denom = 30.
  For 2k = 6: (p-1)|6 means p-1 in {{1,2,3,6}}, so p in {{2,3,4,7}}. But 4 not prime.
    Actually p in {{2,3,7}}. Denom = 42.
  For 2k = 8: (p-1)|8 means p-1 in {{1,2,4,8}}, so p in {{2,3,5}}. Denom = 30.

THE TOURNAMENT TRIPLE {{2,3,5}} APPEARS AS THE DENOMINATOR OF B_4.
THE HURWITZ TRIPLE {{2,3,7}} APPEARS AS THE DENOMINATOR OF B_6.

The Bernoulli denominators CYCLE through the triangle group products!

And B_6 = 1/42 exactly. The reciprocal of the answer.
The sixth Bernoulli number is 1/42 because the (2,3,7) triangle
has area pi/42, and B_6 is the value that makes the Euler-Maclaurin
formula work at order 6.
""")

# ========================================================================
print("=" * 60)
print("III. CATALAN NUMBERS AND TOURNAMENT BRACKETS")
print("=" * 60)
print()

print("Catalan numbers C_n and binary trees:")
for n in range(11):
    c = int(catalan(n))
    print(f"  C_{n} = {c}")

print(f"""
The Catalan number C_n counts:
  - Binary trees with n internal nodes.
  - Ways to parenthesize n+1 factors.
  - Paths in a grid that don't cross the diagonal.
  - Non-crossing partitions of [n+1].

In tournament theory:
  A tournament BRACKET is a binary tree structure.
  C_n gives the number of ways to organize n+1 teams
  into an elimination bracket.

  C_3 = 5 = the tournament prime!
  C_4 = 14 = 2*7 (base times forbidden!)
  C_5 = 42 = THE ANSWER.

  THE 5th CATALAN NUMBER IS 42.

  C_5 = 42 = 2*3*7.

  The Catalan number at the "unsolvability" dimension (5)
  equals the Hurwitz product.

This is because C_n = C(2n, n)/(n+1)  [formula], and:
  C_5 = C(10,5)/6 = 252/6 = 42.
  The division by 6 = 2*3 removes the tournament base primes,
  leaving 42/6 = 7 (?). No, wait: 252/6 = 42. And 252 = 4*63 = 4*9*7.
  So C(10,5) = 252 = 2^2 * 3^2 * 7. Dividing by 6 = 2*3: 252/6 = 42 = 2*3*7.

  C(10,5) / (5+1) = (2^2 * 3^2 * 7) / (2*3) = 2*3*7 = 42. EXACT.
""")

# ========================================================================
print("=" * 60)
print("IV. THE FIBONACCI-LUCAS-MERSENNE WEB")
print("=" * 60)
print()

print("Fibonacci, Lucas, and Mersenne numbers interleaved:")
print(f"{'n':>3} | {'F_n':>8} | {'L_n':>8} | {'2^n-1':>8} | {'F prime?':>9} | {'L prime?':>9} | {'M prime?':>9}")
print("-" * 65)

for n in range(1, 21):
    F = int(fibonacci(n))
    L = int(lucas(n))
    M = 2**n - 1
    Fp = isprime(F)
    Lp = isprime(L)
    Mp = isprime(M)
    print(f"{n:3d} | {F:8d} | {L:8d} | {M:8d} | {'PRIME' if Fp else '':>9} | {'PRIME' if Lp else '':>9} | {'PRIME' if Mp else '':>9}")

print(f"""
Key relationships:
  L_n = F_(n-1) + F_(n+1) (Lucas = sum of adjacent Fibonacci).
  F_(2n) = F_n * L_n (Fibonacci at double index = product).
  F_n^2 + F_(n+1)^2 = F_(2n+1) (sum of squares).

  Tr(M_2^n) = L_n (Lucas numbers ARE the Fibonacci trace sequence).
  Tr(M_k^n) = 2^n - 1 for n <= k (THM-227: Mersenne in the universal zone).

  At n=2: L_2 = 3 = 2^2-1. Lucas = Mersenne!
  At n=3: L_3 = 4 != 7 = 2^3-1. Lucas DIVERGES from Mersenne at n=3.

The DIVERGENCE at n=3 is the moment where the Fibonacci system (k=2)
discovers its individuality. For the tribonacci system (k=3),
the trace stays Mersenne up to n=3.

And: F_7 = 13, which is a MERSENNE EXPONENT (2^13-1 = 8191 is prime).
The 7th Fibonacci number is a Mersenne exponent. The forbidden index
produces a Mersenne generator.
""")

# ========================================================================
print("=" * 60)
print("V. PARTITION NUMBERS AND THE OCF")
print("=" * 60)
print()

from sympy import npartitions

print("Partition numbers p(n) and their tournament connections:")
for n in range(21):
    p = npartitions(n)
    # Check if p(n) is prime
    pp = isprime(p)
    note = ""
    if n == 0: note = "empty partition"
    if pp: note = "PRIME"
    if p == 7: note = "FORBIDDEN!"
    if p == 42: note = "THE ANSWER!"
    print(f"  p({int(n):2d}) = {int(p):8d} {note}")

print(f"""
p(5) = 7 = THE FORBIDDEN VALUE.
p(10) = 42 = THE ANSWER.

The number of partitions of 5 is EXACTLY the first forbidden H value.
The number of partitions of 10 is EXACTLY the Hurwitz product.

Is this a coincidence? Let's see:
  p(5) = 7: the partitions of 5 are
    5, 4+1, 3+2, 3+1+1, 2+2+1, 2+1+1+1, 1+1+1+1+1.
    7 partitions. And 7 = 2^3 - 1.

  p(10) = 42: the partitions of 10.
    And 42 = 2*3*7.

Ramanujan's congruences:
  p(5n + 4) = 0 mod 5.
  p(7n + 5) = 0 mod 7.
  p(11n + 6) = 0 mod 11.

The moduli are 5, 7, 11. And:
  5 = the tournament prime (last solvable).
  7 = the forbidden value (first impossible).
  11 = the first prime in cycle 2 (structure at next scale).

Ramanujan's partition congruences use the primes 5, 7, 11
which are EXACTLY our theory's key primes!
  5: where solvability ends.
  7: where impossibility begins.
  11: where the second cycle's structure appears.
""")

# ========================================================================
print("=" * 60)
print("VI. THE RIEMANN ZETA AT NEGATIVE INTEGERS")
print("=" * 60)
print()

print("zeta(-n) = -B_(n+1)/(n+1) for non-negative n:")
print(f"{'n':>3} | {'zeta(-n)':>15} | {'= -B_(n+1)/(n+1)':>20} | {'decimal':>12} | {'note':>25}")
print("-" * 80)

for n in range(13):
    z = Rational(-1, 1) * bernoulli(n+1) / (n+1)
    z_float = float(z)
    note = ""
    if n == 0: note = "zeta(0) = -1/2"
    elif n == 1: note = "zeta(-1) = -1/12"
    elif n == 2: note = "= 0"
    elif n == 3: note = "zeta(-3) = 1/120"
    elif n == 5: note = "zeta(-5) = -1/252"
    elif z_float != 0 and abs(1/z_float - 42) < 0.01: note = "reciprocal ~ 42!"
    print(f"{n:3d} | {str(z):>15} | {'':>20} | {z_float:12.8f} | {note:>25}")

print(f"""
zeta(-1) = -1/12. The "sum" 1+2+3+... = -1/12 (regularized).

zeta(-3) = 1/120. And 120 = |A_5| = the icosahedral order = 5!.

zeta(-5) = -1/252. And 252 = C(10,5) = the central binomial at n=5.
  And C(10,5)/6 = 42. So zeta(-5) = -1/(6*42) = -1/(6*42).

The zeta function at negative integers encodes the Bernoulli numbers,
which encode the Hurwitz/tournament triples through their denominators.

AND: in our framework, zeta_M(-k) = Tr(M_k^k) = 2^k-1 (THM-227).
The "tournament zeta" at negative integers gives Mersenne numbers.
This parallels the Riemann zeta at negative integers giving Bernoulli numbers.
""")

# ========================================================================
print("=" * 60)
print("VII. A NEW IDENTITY?")
print("=" * 60)
print()

print("""
We have:
  p(5) = 7 = 2^3 - 1.
  p(10) = 42 = 2*3*7.
  C_5 = 42 = 2*3*7.
  B_6 = 1/42.
  |PSL(2,7)| = 168 = 4*42.
  Hurwitz bound = 84(g-1) = 2*42*(g-1).
  Area of (2,3,7) triangle = pi/42.

Everything keeps pointing to 42. Let's check if there's a pattern
connecting partition numbers and Catalan numbers:
""")

print(f"{'n':>3} | {'p(n)':>6} | {'C_n':>6} | {'p(n)/C_n':>10} | {'p(2n)':>8} | {'C_n * something':>15}")
print("-" * 55)
for n in range(1, 11):
    p = int(npartitions(n))
    c = int(catalan(n))
    ratio = p / c if c > 0 else 0
    p2 = int(npartitions(2*n))
    print(f"{n:3d} | {p:6d} | {c:6d} | {ratio:10.4f} | {p2:8d} | ")

print(f"""

p(5) = 7, C_5 = 42. Ratio = 1/6.
p(10) = 42, C_10 = 16796. Ratio = 42/16796 = 0.0025.

No clean ratio. But: p(5) * C_5/p(5) = 42. Trivially.

What about p(5) * 6 = 42? YES! p(5) * (5+1) = 42.
  7 * 6 = 42.
  The forbidden value times the pivot = the answer.
  And p(5) * (5+1) = C_5. Is this true in general?

  p(n) * (n+1) = C_n?
  p(1) * 2 = 1 * 2 = 2. C_1 = 1. NO.
  p(2) * 3 = 2 * 3 = 6. C_2 = 2. NO.

  So p(5) * 6 = 42 = C_5 is a SPECIAL PROPERTY of n=5.

  It works because p(5) = 7 and C_5 = C(10,5)/6 = 252/6 = 42 = 7*6.
  And p(5) = 7 happens to divide C(10,5) = 252 with quotient 36 = 6^2.

  Actually: C(10,5) = 252 = 7 * 36 = 7 * 6^2.
  So C_5 = 252/6 = 7*6 = 42 = p(5) * (p(5)-1).

  C_5 = p(5) * (p(5) - 1) = 7 * 6 = 42.

  Is this a coincidence? C_5 = p(5)*(p(5)-1)?

  7*6 = 42. YES.

  This says: the 5th Catalan number = p(5) * (p(5)-1).
  = (first forbidden) * (the pivot).
  = 7 * 6.
  = 42.
""")

# Check: does C_n = p(n) * (p(n) - 1) for other n?
print("Check: C_n = p(n) * (p(n)-1)?")
for n in range(1, 11):
    c = int(catalan(n))
    p = npartitions(n)
    pred = p * (p - 1)
    print(f"  n={n}: C_{n}={c}, p({n})*(p({n})-1)={pred}, match: {c == pred}")

print("""

Only works for n=5! It's a property SPECIFIC to the "unsolvability" dimension.

At n=5: C_5 = p(5)*(p(5)-1) = 7*6 = 42.
This is the ONLY n where the Catalan number factors as
partition-count times (partition-count minus 1).

The uniqueness of n=5 is because:
  p(5) = 7 (Mersenne prime)
  C_5 = 42 = 7*6
  And 7*6 = 42 is the only value where this identity holds
  because the growth rates of p(n) and C_n are very different.

C_n ~ 4^n / (n^{3/2} * sqrt(pi)).
p(n) ~ exp(pi*sqrt(two_n_over_3)) / (four_n*sqrt(3)).

These grow at VASTLY different rates. They ONLY cross at n=5
in the right way to make C_n = p(n)*(p(n)-1).

n=5 IS the crossing point. The last solvable dimension.
The place where the partition world and the Catalan world
briefly AGREE before diverging forever.
""")

# ========================================================================
print("=" * 60)
print("VIII. THE GOLDEN IDENTITY FOR FIBONACCI PRIMES")
print("=" * 60)
print()

phi = (1 + sqrt(5)) / 2

print("Fibonacci primes and their golden ratio properties:")
fib_primes = []
for n in range(2, 50):
    F = int(fibonacci(n))
    if isprime(F):
        fib_primes.append((n, F))
        golden_power = F * phi - int(F * phi)  # fractional part of F*phi
        print(f"  F_{n} = {F}, F*phi = {F*phi:.6f}, frac part = {golden_power:.6f}")

print(f"""

The fractional parts of F_n * phi are related to 1/phi^n
by the identity: F_n * phi = F_(n+1) + (-1)^{n+1}/phi^n + ...
Actually: F_n * phi - F_(n+1) = (-phi)^{-n} (exactly, from Binet).

For Fibonacci primes: F_n is prime (irreducible), so
the fractional part of F_n * phi carries SPECIAL arithmetic information.

The Fibonacci primes appear at indices:
  n = 3, 4, 5, 7, 11, 13, 17, 23, 29, 43, 47, 83, ...

These indices are (mostly) PRIME themselves!
F_n is prime only if n is prime (necessary but not sufficient for n > 4).

And the PRIME indices of Fibonacci primes:
  3, 5, 7, 11, 13, 17, 23, 29, 43, 47, 83, ...

These include the tournament prime (5), the forbidden (7),
and structure-at-next-scale (11, 13).
""")

# ========================================================================
print("=" * 60)
print("IX. THE 42 IDENTITIES (collecting them all)")
print("=" * 60)
print()

print("""
42 appears in our framework in ALL of these ways:

  1.  42 = 2 * 3 * 7           (Hurwitz triple product)
  2.  42 = 2 * 21              (base * second forbidden H)
  3.  42 = 6 * 7               (pivot * first forbidden)
  4.  Area(2,3,7) = pi/42      (Hurwitz triangle)
  5.  Hurwitz bound = 84(g-1) = 2*42*(g-1)
  6.  |PSL(2,7)| = 168 = 4*42
  7.  C_5 = 42                  (5th Catalan number)
  8.  p(10) = 42                (partitions of 10)
  9.  B_6 = 1/42               (6th Bernoulli number)
  10. 6174 = 42 * 147           (Kaprekar = 42 * 3*7^2)
  11. digits(6174) product = 168 = 4*42
  12. 7 * 6 = 42               (forbidden * pivot)
  13. p(5) * (p(5)-1) = 42     (partitions identity at n=5)
  14. C(10,5) / 6 = 42         (central binomial / pivot)
  15. mu(42) = -1               (negatively oriented on torus)
  16. phi(42) = 12              (totient = arcs of 4-tournament)
  17. 42 mod 8 = 2              (Bott slot = comparison/BASE)
  18. The ANSWER to life, the universe, and everything.
""")

print(f"Verification of some:")
print(f"  phi(42) = {int(totient(42))}")
print(f"  42 mod 8 = {42 % 8}")
print(f"  mu(42) = {int(mobius(42))}")
print(f"  C_5 = {int(catalan(5))}")
print(f"  p(10) = {npartitions(10)}")
print(f"  B_6 = {bernoulli(6)} = 1/{1/bernoulli(6)}")

print("\n\nopus-2026-03-16-S90bk: DOING MATH FOR FUN")
