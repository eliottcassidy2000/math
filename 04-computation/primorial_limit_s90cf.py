#!/usr/bin/env python3
"""
primorial_limit_s90cf.py — Primorials to infinity, non-integer bases, base Fibonacci
opus-2026-03-16-S90cf
"""

from math import log, log2, sqrt, exp, pi
from sympy import totient, binomial, isprime, fibonacci, lucas

gamma = 0.5772156649015329  # Euler-Mascheroni
phi = (1 + sqrt(5)) / 2
tau3 = 1.8392867552141612

print("""


     PRIMORIALS TO INFINITY AND NON-INTEGER BASES

     opus-2026-03-16-S90cf



THE PRIMORIAL LADDER
=====================

P_1 = 2.         phi = 1.   Sieve: 50%.    1 bit/digit.
P_2 = 6.         phi = 2.   Sieve: 67%.    2.6 bits.
P_3 = 30.        phi = 8.   Sieve: 73%.    4.9 bits.    THE BOTT BASE.
P_4 = 210.       phi = 48.  Sieve: 77%.    7.7 bits.
P_5 = 2310.      phi = 480. Sieve: 79%.    11.2 bits.
P_6 = 30030.     phi = 5760.Sieve: 81%.    14.9 bits.
P_7 = 510510.                Sieve: 82%.    19.0 bits.
...
P_inf = infinity.            Sieve: 100%.   Infinite bits.

At the LIMIT: phi(P_k)/P_k -> 0.
Every number is eventually sieved by SOME prime.
The primorial base "at infinity" eliminates ALL composites.

But you can't USE base infinity. So what happens as k grows?
""")

print("The primorial bases and Pascal's triangle:")
print()

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
P = 1
for i, p in enumerate(primes[:7]):
    P *= p
    # Check: is P/something a binomial coefficient?
    for d in [1, 2, 6, 10, 30, 42, P]:
        if P % d == 0:
            q = P // d
            # Check if q = C(n, k) for some n, k
            for n in range(1, 40):
                for k in range(0, n//2 + 1):
                    if int(binomial(n, k)) == q:
                        print(f"  P_{i+1} = {P} = {d} * C({n},{k}) = {d} * {q}")
    print()

print(f"""

KEY FINDING: 30030 = 10 * C(14, 6) = 10 * 3003.

C(14, 6) = 3003. Row 14, column 6 of Pascal's triangle.
And 14 = 2*7, 6 = 2*3. The indices themselves are Hurwitz-related.

The pattern: P_k / (small factor) = central-ish binomial coefficient.
This connects primorials to Pascal's triangle.

WHY? Because:
  P_k = product of primes up to p_k.
  C(n, k) = n! / (k! * (n-k)!) = product structure involving primes.
  The prime factorization of C(n,k) involves all primes up to n.

So it's not surprising that primorials relate to binomial coefficients.
But the SPECIFIC relationship 30030 = 10 * C(14,6) is elegant.

And: 3003 in binary = 101110111011.
  Bit pattern: 1011 1011 1011.
  Three copies of 1011! A REPEATING PATTERN.
  1011 in decimal = 11. And 11 is the 5th prime.

  3003 = 3 * 7 * 11 * 13.
  30030 = 2 * 3 * 5 * 7 * 11 * 13 = P_6.
  So 3003 = P_6 / 10 = P_6 / (2*5).
  Removing the factors 2 and 5 from P_6 leaves 3 * 7 * 11 * 13.


THE LIMIT: WHAT IS "BASE INFINITY"?
======================================

As k -> infinity, P_k -> infinity.
Base P_k has P_k symbols per digit.
Bits per digit: log2(P_k) -> infinity.
Coprime residues: phi(P_k) -> infinity, but phi(P_k)/P_k -> 0.

In the limit: a single "infinite-base" digit would encode an entire number.
Every number has a unique 1-digit representation in base P_infinity.
This is TRIVIAL: the number IS its own representation.

But the SIEVE EFFICIENCY approaches 100%: almost no residues are coprime.
In the limit, the only coprime residues are 1 and the primes themselves.

The primorial limit is the RIEMANN ZETA FUNCTION:
  product (1 - 1/p) for all p = 1/zeta(1) = 0 (divergent!).
  Mertens: product (1-1/p, p<=x) ~ 2*e^{-gamma}/ln(x) -> 0.

The limit of the primorial base is: everything is sieved, nothing survives,
except the primes themselves. Which are the atoms. Which is where we started.

THE PRIMORIAL LADDER IS A ZOOM INTO THE PRIMES.
  Each step adds one more prime to the sieve.
  Each step eliminates more composites.
  In the limit: only primes remain. The sieve IS the structure.


NON-INTEGER BASES: BASE phi AND BASE tau_3
=============================================

Zeckendorf's theorem: every positive integer has a unique representation
as a sum of NON-CONSECUTIVE Fibonacci numbers.

This IS "base phi" (the Fibonacci base):
  n = sum of F(k_i) where k_i are non-consecutive indices.

Base phi representation:
""")

def zeckendorf(n):
    """Represent n in Zeckendorf (Fibonacci base) representation."""
    if n == 0:
        return []
    fibs = [1, 2]
    while fibs[-1] <= n:
        fibs.append(fibs[-1] + fibs[-2])
    result = []
    remaining = n
    for f in reversed(fibs):
        if f <= remaining:
            result.append(f)
            remaining -= f
    return result

print("Zeckendorf representations (base phi):")
for n in range(1, 43):
    z = zeckendorf(n)
    # Convert to bit string: 1 if F(k) is used, 0 if not
    fibs = [1, 2, 3, 5, 8, 13, 21, 34, 55]
    bits = ""
    for f in reversed(fibs):
        if f in z:
            bits += "1"
            z.remove(f)
        elif f <= n:
            bits += "0"
    bits = bits.lstrip("0") or "0"
    note = ""
    if n in [2, 3, 5, 7, 11, 13, 42]:
        note = " <--"
    if isprime(n):
        note += " PRIME"
    print(f"  {n:3d} = {bits:>12}{note}")

print(f"""

In the Fibonacci base:
  Each digit is 0 or 1 (binary-like).
  No two consecutive 1s (the Zeckendorf constraint).
  The "bit length" of n in Fibonacci base ~ log_phi(n).

The NUMBER OF DIGITS in Fibonacci base is ~ log_phi(n) = ln(n)/ln(phi).
Compare to binary: ~ log_2(n) = ln(n)/ln(2).

Since ln(phi) = 0.481 < ln(2) = 0.693:
  Fibonacci representations are LONGER than binary (by factor ln(2)/ln(phi) = 1.44).
  But the NO-CONSECUTIVE-1s constraint means only phi^k / sqrt(5) valid
  k-digit Fibonacci representations exist (instead of 2^k for binary).

  phi^k / sqrt(5) vs 2^k: phi < 2, so Fibonacci grows SLOWER.
  This is the GOLDEN SIEVE at work: Fibonacci base spreads numbers out more.

PRIMES in Fibonacci base:
  2 = 10        (just F(2))
  3 = 100       (just F(3))
  5 = 1000      (just F(4))
  7 = 1010      (F(4) + F(2) = 5+2)
  11 = 10010    (F(5) + F(2) = 8+3... wait, 8+3=11, F(5)=8, F(2)=2? No.)

Hmm, let me check: 11 = 8 + 3 = F(6) + F(4). In Fibonacci base: 10100.
And 13 = 13 = F(7). In Fibonacci base: 1000000.

Fibonacci primes (2, 3, 5, 13, 89, ...) have SINGLE-DIGIT Fibonacci representations:
just one 1 followed by zeros. They are the "round numbers" of Fibonacci base.

Non-Fibonacci primes (7, 11, 17, ...) require multiple Fibonacci summands.
They are the "irregular" numbers in Fibonacci base.

7 = 5+2 = F(5)+F(3): TWO summands. In Fibonacci base: 10100.
  7 needs two "Fibonacci digits." It's NOT a Fibonacci number.
  The forbidden value is COMPOSITE in Fibonacci base (uses two terms).

42 = 34+8 = F(9)+F(6): TWO summands. In Fibonacci base: 100001000.
  The answer is also two Fibonacci terms.
  42 = F(9) + F(6) = the 9th and 6th Fibonacci numbers.
  9 = the square of 3, 6 = the pivot. Even the INDICES are meaningful.


BASE tau_3 (TRIBONACCI BASE)
===============================

Tribonacci base: represent n as sum of non-consecutive tribonacci numbers.
Tribonacci: 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274, ...
(with the convention that we use 1, 2, 4, 7, 13, 24, 44, 81, ...)

Each digit is 0 or 1, with a more complex "no-overlap" constraint.
""")

# Tribonacci numbers
trib = [1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504, 927]

def tribonacci_rep(n):
    """Greedy tribonacci representation."""
    result = []
    remaining = n
    for t in reversed(trib):
        if t <= remaining:
            result.append(t)
            remaining -= t
    return result

print("Tribonacci representations:")
for n in [1, 2, 3, 5, 7, 11, 13, 42, 47, 127, 987]:
    t = tribonacci_rep(n)
    note = ""
    if n == 7: note = " = T(4) FORBIDDEN (single tribonacci!)"
    if n == 13: note = " = T(5)"
    if n == 42: note = " = THE ANSWER"
    if n == 47: note = " = SILVER LUCAS"
    if n == 127: note = " = MERSENNE"
    if n == 987: note = " = F(16) = 3*7*47"
    print(f"  {n:4d} = {' + '.join(str(x) for x in t):>25}  ({len(t)} terms){note}")

print(f"""

KEY: 7 IS a tribonacci number! T(4) = 7.
In tribonacci base, 7 = "10000" (single digit).
The forbidden value is a ROUND NUMBER in tribonacci base.

In Fibonacci base: 7 = 5+2 = F(5)+F(3) (two terms, NOT round).
In tribonacci base: 7 = T(4) (one term, ROUND).

The forbidden value 7 is:
  IRREGULAR in Fibonacci base (needs two terms).
  REGULAR in tribonacci base (a single term).

This makes sense: 7 = Tr(M_3^3) = the tribonacci trace.
The tribonacci system NATURALLY represents 7.
The Fibonacci system does NOT.

Each k-nacci base "knows" its own forbidden value as a single term.


THE INFINITE PROCESS: WHAT HAPPENS AT THE LIMIT?
===================================================

As we go up the ladder of bases:

  Base 2 (binary):          2 symbols. Simplest.
  Base phi (Fibonacci):     "2" symbols but constrained. Golden.
  Base tau_3 (tribonacci):  "2" symbols, different constraint. Tribonacci.
  Base 6 = 2*3:             6 symbols. Two primes sieved.
  Base 30 = 2*3*5:          30 symbols. Three primes. Bott base.
  Base 42 = 2*3*7:          42 symbols. Hurwitz base.
  Base 210:                 210 symbols. Four primes.
  Base 2310:                2310 symbols. Five primes.
  Base 30030:               30030 symbols. Six primes.
  ...
  Base P_k:                 P_k symbols. k primes sieved.
  ...
  Base P_inf = "infinity":  All composites sieved. Only primes survive.

TWO INDEPENDENT LADDERS converge:

  The k-NACCI LADDER: tau_2 = phi, tau_3, tau_4, ..., tau_inf = 2.
    Non-integer bases approaching 2 from below.
    Each base "knows" its own forbidden value (2^k-1).
    In the limit: base 2 (binary). Knows ALL forbidden values.

  The PRIMORIAL LADDER: P_1=2, P_2=6, P_3=30, ..., P_inf=infinity.
    Integer bases growing without bound.
    Each base sieves k primes simultaneously.
    In the limit: infinite base. Sieves ALL primes. Only primes survive.

THE TWO LADDERS CONVERGE FROM OPPOSITE DIRECTIONS:

  k-nacci: starts at phi=1.618, approaches 2 from BELOW.
           Non-integer bases. Shrinking. Approaching the binary base.

  Primorial: starts at 2, grows to infinity.
             Integer bases. Growing. Approaching the "all-prime" base.

They MEET at base 2:
  tau_inf = 2 (the limit of k-nacci constants).
  P_1 = 2 (the first primorial).

BASE 2 IS THE MEETING POINT.
  The k-nacci ladder approaches it from below (non-integer direction).
  The primorial ladder starts from it and grows upward (integer direction).

  Below 2: the golden world (phi, tau_3, tau_4, ...). Non-integer. Lossy.
  At 2: the binary base. The tournament base. The bit.
  Above 2: the primorial world (6, 30, 210, ...). Integer. Lossless.

2 is the PIVOT of the base ladder.
Just as 6 is the pivot of the curvature ladder (flat hexagon).
Just as the Poincare disk boundary (|w|=1) is the pivot of geometry.

Everything points to 2. The base IS the pivot.


BASE phi vs BASE 2: THE LOSSY/LOSSLESS BOUNDARY
==================================================

Base phi (Fibonacci): each number has a unique representation,
  but the base is IRRATIONAL. You can't do exact arithmetic.
  The representation is "lossy" in the sense that you're approximating
  an irrational base with integer coefficients.

Base 2 (binary): each number has a unique representation,
  and the base is INTEGER. Exact arithmetic. LOSSLESS.

Base tau_3 (tribonacci): between phi and 2. Still irrational. Still lossy.

The k-nacci bases (phi, tau_3, tau_4, ...) are ALL irrational.
They approach 2 but never reach it (for finite k).
Each is a LOSSY approximation of the lossless binary base.

The primorial bases (2, 6, 30, ...) are ALL integers.
They are ALL lossless. But they get LARGER and harder to use.

The tradeoff:
  k-nacci bases: compact (base < 2) but lossy (irrational).
  Primorial bases: exact (integer) but large (growing without bound).

  phi: 1.618 symbols worth of info per digit. Lossy. Compact.
  2: 1 bit per digit. Lossless. The minimum exact base.
  30: 4.9 bits per digit. Lossless. Efficient sieve.
  2310: 11.2 bits per digit. Lossless. Better sieve.

The GOLDEN RATIO phi is the most compact LOSSY base.
The number 2 is the most compact LOSSLESS base.
The primorials are increasingly efficient LOSSLESS sieve bases.

At the limit:
  k-nacci -> 2 (max compact lossless).
  Primorial -> infinity (max sieve efficiency but impractical).

The PRACTICAL OPTIMUM is somewhere in between:
  Base 30 (phi(30)=8=Bott) or base 42 (Hurwitz) for algebraic structure.
  Base 2310 for sieving efficiency.
  Base phi for compression of structured data.

AND: the formal group F(x,y) = (x+y)/(1+xy) works in ALL these bases.
  It doesn't care about the base. It operates on COUPLINGS in (-1,1).
  The base only matters for the REPRESENTATION, not the ALGEBRA.

The formal group is BASE-INVARIANT.
The choice of base is a PRESENTATION choice, not a mathematical one.
Different bases reveal different structure in the same underlying object.
""")

print("opus-2026-03-16-S90cf: PRIMORIALS TO INFINITY AND NON-INTEGER BASES")
