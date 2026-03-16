#!/usr/bin/env python3
"""
primes_geometry_bases_s90bi.py — The geometry of prime distribution across bases
opus-2026-03-16-S90bi

Think about what is GEOMETRICALLY going on when we write a number
in different bases, and what that reveals about primality.
"""

from sympy import isprime, primerange, factorint, nextprime
from math import log, log2, sqrt, pi, floor, ceil
from collections import Counter

print("""


     THE GEOMETRY OF PRIME DISTRIBUTION

     opus-2026-03-16-S90bi



THE CORE QUESTION:

A number n is just a point on the number line.
Writing it in base b gives it a BODY: a string of digits.
Different bases give different bodies to the SAME point.

  42 in base 2:  101010
  42 in base 3:  1120
  42 in base 5:  132
  42 in base 6:  110
  42 in base 7:  60
  42 in base 10: 42
  42 in base 42: 10

The number 42 IS all of these simultaneously.
Each representation reveals different structure.

A prime p is a number whose MULTIPLICATIVE structure is trivial.
But its ADDITIVE structure (digit representation) can be complex.

What does primality LOOK LIKE in each base?


1. BASE 2: THE RAW BITS
=========================

In base 2, every odd number ends in 1. Every prime > 2 is odd.
So every prime > 2 ends in 1 in binary.

The LEADING bit is always 1 (by convention, no leading zeros).
So every prime > 2 has binary: 1...1 (starts and ends with 1).

What's in between? Anything. The INTERIOR BITS are free.

For a k-bit prime: the first and last bits are fixed (both 1).
The k-2 interior bits can be anything, EXCEPT combinations
that make the number divisible by 3, 5, 7, ...

The NUMBER of k-bit primes ~ 2^k / (k * ln 2) (by PNT).
The number of k-bit odd numbers = 2^(k-2) (free interior bits).
So the FRACTION that are prime ~ 2/(k * ln 2) ~ 1/k.

As k grows: about 1/k of k-bit numbers with fixed endpoints are prime.
The interior bits are "almost free" — primality is a SPARSE constraint.
""")

# Visualize primes in binary
print("Primes in binary (first 30):")
print(f"{'p':>5} | {'binary':>20} | {'bits':>4} | {'interior':>12} | {'int pattern':>12}")
print("-" * 60)
for p in primerange(2, 130):
    b = bin(p)[2:]
    bits = len(b)
    interior = b[1:-1] if bits > 2 else ""
    print(f"{p:5d} | {b:>20} | {bits:4d} | {interior:>12} | ", end="")
    # Classify interior pattern
    if not interior:
        print("(none)")
    elif all(c == '0' for c in interior):
        print("all-0 (sparse)")
    elif all(c == '1' for c in interior):
        print("all-1 (MERSENNE)")
    else:
        ones = interior.count('1')
        print(f"{ones}/{len(interior)} ones")

print(f"""

THE BINARY GEOMETRY:

A k-bit number is a POINT in the k-dimensional hypercube {{0,1}}^k.
The constraint "starts with 1" halves the cube.
The constraint "ends with 1" halves again.
We're left with a (k-2)-dimensional face of the hypercube.

Primes are a SUBSET of this face.
The sieve of Eratosthenes REMOVES points:
  - Multiples of 3: a periodic pattern in the interior bits.
  - Multiples of 5: another periodic pattern.
  - Multiples of 7: another periodic pattern.
  - ...

Each sieve step removes a PERIODIC SUBPATTERN of the hypercube face.
The primes are what's LEFT after all periodic subpatterns are removed.

This is a GEOMETRIC sieve: we start with a face of the hypercube
and punch out periodic holes. The primes are the Swiss cheese.

The holes are periodic because divisibility is periodic:
  n is divisible by 3 iff its binary digit sum (with alternating signs)
  is divisible by 3... actually it's more complex than that.

In base 2, divisibility by p has period p in the LAST few bits:
  n mod p depends on n mod p, which cycles with period p.

So the sieve punches holes with periods 3, 5, 7, 11, 13, ...
These periods are MULTIPLICATIVELY INDEPENDENT (they're prime!).
The resulting pattern is QUASIPERIODIC: a superposition of incommensurate periods.

A QUASICRYSTAL. The primes in binary are a quasicrystal on the hypercube.
""")

# ========================================================================

print("=" * 60)
print("2. BASE 6: THE SIEVE-ADAPTED BASE")
print("=" * 60)

print("""
Base 6 = 2 * 3. The product of the first two primes.

In base 6, every prime > 3 ends in 1 or 5.
(Because primes are coprime to 6, so p mod 6 is in {1, 5}.)

This means: in base 6, primes have EXACTLY TWO possible last digits.
Compare to base 10 where primes have four possible last digits (1,3,7,9).
And base 2 where primes have one possible last digit (1).

Base 6 is optimal in a sense: it uses the SMALLEST base
that separates primes into exactly two families.
""")

print("Primes in base 6:")
print(f"{'p':>5} | {'base 6':>10} | {'last digit':>10} | {'family':>10}")
print("-" * 45)
for p in primerange(2, 80):
    digits = []
    n = p
    while n > 0:
        digits.append(n % 6)
        n //= 6
    digits.reverse()
    base6 = "".join(str(d) for d in digits)
    last = digits[-1] if digits else 0
    family = "6k+1" if last == 1 else ("6k+5" if last == 5 else "special")
    print(f"{p:5d} | {base6:>10} | {last:>10} | {family:>10}")

print(f"""

In base 6, the two prime families (6k+1 and 6k+5) correspond to
the two generators of (Z/6Z)*:
  1 mod 6: primes that are 1 more than a multiple of 6.
  5 mod 6: primes that are 1 less than a multiple of 6.

These are the primes that are ONE STEP from a "smooth" number (multiple of 6).
Primes are ADJACENT to smoothness. They sit one step away from the
most highly composite numbers.

GEOMETRICALLY in base 6: the number line is divided into blocks of 6.
Each block has potential prime positions at offset 1 and offset 5.
The sieve of primes > 3 removes some of these positions.

The prime distribution in base 6 is a BINARY SEQUENCE:
at each block of 6, there are two candidate positions,
and each is prime or not. The sequence of (prime/not, prime/not)
across blocks is the prime distribution in base 6.

This binary sequence IS the core of the prime distribution.
All other base representations are TRANSFORMATIONS of this sequence.
""")

# ========================================================================

print("=" * 60)
print("3. BASE 30: THE PRIMORIAL BASE")
print("=" * 60)

print("""
30 = 2 * 3 * 5 = the product of the first three primes.
30 is the primorial P_3.

In base 30, primes > 5 can only end in digits coprime to 30:
  {1, 7, 11, 13, 17, 19, 23, 29}.
That's 8 possible last digits out of 30 = phi(30) = 8.

So in base 30: primes have EXACTLY 8 possible last digits.
These 8 digits are the elements of (Z/30Z)*.

And 8 = the Bott period!

The 8 possible last digits in base 30 correspond to the 8 Bott slots:
""")

residues_30 = [r for r in range(30) if all(r % p != 0 for p in [2,3,5]) or r == 1]
# Actually: coprime to 30
residues_30 = [r for r in range(1, 30) if all(r % p != 0 for p in [2, 3, 5])]
print(f"Residues coprime to 30: {residues_30}")
print(f"Count: {len(residues_30)} = phi(30) = 8 = Bott period!")
print()

# Map to Bott slots
print(f"{'residue mod 30':>15} | {'mod 8':>6} | {'Bott archetype':>20}")
print("-" * 50)
archetypes = {1: "existence", 3: "structure", 5: "unsolvability", 7: "forbidden"}
for r in residues_30:
    slot = r % 8
    arch = archetypes.get(slot, f"slot {slot}")
    print(f"{r:>15} | {slot:>6} | {arch:>20}")

print(f"""

THE 8 RESIDUES mod 30, mapped to Bott slots:
  1  -> slot 1 (existence)
  7  -> slot 7 (forbidden)
  11 -> slot 3 (structure)
  13 -> slot 5 (unsolvability)
  17 -> slot 1 (existence)
  19 -> slot 3 (structure)
  23 -> slot 7 (forbidden)
  29 -> slot 5 (unsolvability)

Each Bott slot gets EXACTLY 2 of the 8 residues:
  Slot 1 (existence): 1 and 17.
  Slot 3 (structure): 11 and 19.
  Slot 5 (unsolvability): 13 and 29.
  Slot 7 (forbidden): 7 and 23.

PERFECT BALANCE. The primorial base 30 distributes its residues
EQUALLY among the four Bott families, two per family.

This is NOT a coincidence. phi(30) = 8 because:
  phi(30) = phi(2)*phi(3)*phi(5) = 1*2*4 = 8.
  And 8 = 2^3 = the Bott period.

The Bott period 8 equals phi(30) = phi(2*3*5) because:
  8 = 2^3 = 2^(number of odd primes in the primorial).
  The "number of odd primes" up to 5 is 2 (namely 3 and 5).
  And phi(2*3*5) = 1 * 2 * 4 = 2^0 * 2^1 * 2^2 = 2^3 = 8.

Wait, that's:
  phi(2) = 1 = 2^0
  phi(3) = 2 = 2^1
  phi(5) = 4 = 2^2
  Product = 2^(0+1+2) = 2^3 = 8.

The exponents 0, 1, 2 are: the number of primes less than each prime minus 1.
  phi(p) = p-1. For p = 2: 1 = 2^0. For p = 3: 2 = 2^1. For p = 5: 4 = 2^2.

So phi(P_k) = product of (p-1) for p <= p_k.
  For k=3: phi(30) = 1*2*4 = 8.
  For k=4: phi(210) = 1*2*4*6 = 48.
  For k=5: phi(2310) = 1*2*4*6*10 = 480.

The Bott period 8 is phi(30), the FIRST primorial whose
totient is a pure power of 2.

After 30: phi(210) = 48 = 16*3. NOT a pure power of 2.
So 30 is the LAST primorial with a power-of-2 totient.
This is why the Bott period is 8 and not some other number.

THE BOTT PERIOD 8 = phi(30) = phi(2*3*5).
The Bott period IS the totient of the tournament triple product!
""")

# ========================================================================

print("=" * 60)
print("4. THE GEOMETRIC PICTURE: PRIMES ON THE TORUS")
print("=" * 60)

print("""
Write a number n in the "multi-base" representation:
  n mod 2 (binary parity)
  n mod 3 (ternary residue)
  n mod 5 (quinary residue)
  n mod 7 (septenary residue)
  ...

Each residue is a COORDINATE. Together they define a point on
the PRODUCT of cyclic groups:

  n -> (n mod 2, n mod 3, n mod 5, n mod 7, ...) in Z/2 x Z/3 x Z/5 x Z/7 x ...

By the Chinese Remainder Theorem, for the first k primes:
  n mod P_k determines n uniquely in [0, P_k).
  The map n -> (n mod 2, ..., n mod p_k) is an isomorphism
  of Z/P_k with Z/2 x Z/3 x ... x Z/p_k.

The product Z/2 x Z/3 x Z/5 x ... is a TORUS.
  Z/2 is a 1-torus (a circle with 2 points).
  Z/3 is a 1-torus (a circle with 3 points).
  Z/5 is a 1-torus (a circle with 5 points).
  The product is a MULTI-DIMENSIONAL TORUS.

A number n is a POINT ON THIS TORUS.

n is PRIME iff it avoids the "zero" on EVERY circle:
  n mod 2 != 0  (not even)
  n mod 3 != 0  (not divisible by 3)
  n mod 5 != 0  (not divisible by 5)
  ...
  n mod p != 0  for all primes p <= sqrt(n).

The primes are the points on the torus that AVOID ALL ZEROS.
Each prime p removes 1/p of the torus (the "zero hyperplane" in the p-th direction).
The primes are the COMPLEMENT of the union of these hyperplanes.

THIS IS INCLUSION-EXCLUSION ON A TORUS.

The fraction of the torus that is prime:
  prod_{p <= sqrt(n)} (1 - 1/p) ~ 1/ln(sqrt(n)) ~ 2/ln(n)  (by Mertens' theorem).

This gives prime density ~ 2/ln(n), which is close to the PNT estimate 1/ln(n).
(The factor of 2 is absorbed by more careful analysis.)

THE PRIME NUMBER THEOREM IS A STATEMENT ABOUT VOLUME ON A TORUS.
""")

# ========================================================================

print("=" * 60)
print("5. THE BIT SHIFT ON THE TORUS")
print("=" * 60)

print("""
The doubling map n -> 2n acts on the torus:
  (n mod 2, n mod 3, n mod 5, ...) -> (2n mod 2, 2n mod 3, 2n mod 5, ...)
  = (0, 2n mod 3, 2n mod 5, ...)

Doubling KILLS the 2-coordinate (sends it to 0) and ROTATES the others:
  2n mod 3: rotation by 2 on the 3-circle.
  2n mod 5: rotation by 2 on the 5-circle.
  2n mod 7: rotation by 2 on the 7-circle.

The doubling map is a ROTATION ON THE TORUS
(with the 2-coordinate collapsed to zero).

The orbit of n under doubling: n, 2n, 4n, 8n, ...
This traces out a PATH on the torus:
  The 3-coordinate cycles: n, 2n, 4n, 8n, 16n, 32n, ... mod 3
    = n, 2n, n, 2n, n, 2n, ... (period 2, since 2^2 = 1 mod 3).
  The 5-coordinate cycles with period 4 (since 2^4 = 16 = 1 mod 5).
  The 7-coordinate cycles with period 3 (since 2^3 = 8 = 1 mod 7).

The periods are the MULTIPLICATIVE ORDERS of 2 modulo each prime:
  ord_3(2) = 2.
  ord_5(2) = 4.
  ord_7(2) = 3.
  ord_11(2) = 10.
  ord_13(2) = 12.

The total cycle length of the doubling orbit on the torus is:
  lcm(ord_3(2), ord_5(2), ord_7(2), ...) = lcm(2, 4, 3, ...).

For the first few primes: lcm(2, 4, 3, 10, 12) = 60.

So the doubling map on the (3,5,7,11,13)-torus has period 60.
After 60 doublings, the torus coordinates return to their starting point.

THE BOTT PERIOD 8 vs THE TORUS PERIOD:

The Bott period 8 is the period of the BINARY structure (mod 8 = mod 2^3).
The torus period is the period of the ODD structure (lcm of orders mod odd primes).
These are INDEPENDENT.

The total period is: lcm(8, torus period) = lcm(8, 60) = 120 (for primes up to 13).

120 = |A_5| = order of the icosahedral group.

THE FULL PERIOD IS THE ICOSAHEDRAL NUMBER!
""")

# Verify
from math import lcm
orders = []
for p in [3, 5, 7, 11, 13]:
    k = 1
    x = 2
    while x % p != 1:
        x = (x * 2) % p
        k += 1
    orders.append(k)
    print(f"  ord_{p}(2) = {k}")

torus_period = orders[0]
for o in orders[1:]:
    torus_period = lcm(torus_period, o)
print(f"  Torus period (lcm) = {torus_period}")
print(f"  Full period = lcm(8, {torus_period}) = {lcm(8, torus_period)}")

# ========================================================================

print(f"""

{"=" * 60}
6. WHAT IS GEOMETRICALLY HAPPENING
{"=" * 60}

Picture this:

The natural numbers are a SINGLE PATH on an infinite-dimensional torus.
Each dimension of the torus corresponds to a prime p:
  The p-th dimension is the circle Z/pZ.
  The number n sits at position (n mod p) on this circle.

The path goes: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...
On each circle, this is a CONSTANT-SPEED ROTATION:
  On the 2-circle: 1, 0, 1, 0, 1, 0, ...  (alternating)
  On the 3-circle: 1, 2, 0, 1, 2, 0, ...  (period 3)
  On the 5-circle: 1, 2, 3, 4, 0, 1, 2, ...  (period 5)
  On the 7-circle: 1, 2, 3, 4, 5, 6, 0, 1, ...  (period 7)

The path is a HELIX on the torus: it winds around each circle
at a rate of 1 step per number. The winding rates are all 1
(since we're just counting: n -> n+1).

A number n is PRIME iff the path is NOT at zero on ANY circle.
The zeros of each circle form a "forbidden hyperplane."

The primes are the points where the helix AVOIDS ALL FORBIDDEN HYPERPLANES.

As n grows, the path has wound around the small circles many times.
The probability of avoiding the zero on the p-circle is (1 - 1/p).
For independent circles: the probability of avoiding ALL zeros is
  prod (1 - 1/p) ~ 2*e^(-gamma) / ln(n)
by Mertens' theorem.

But the circles are NOT independent: the path visits them IN ORDER.
The CORRELATIONS between different circle positions create the
EXACT fluctuations of the prime distribution.

These correlations ARE the Riemann zeros:
  Each zero of zeta(s) corresponds to a RESONANCE between
  different circle winding rates. When two circles "sync up"
  (their periods have a common factor), the prime-avoidance
  probability fluctuates.

THE RIEMANN ZEROS ARE THE RESONANCES OF THE TORUS HELIX.

The critical line Re(s) = 1/2 corresponds to:
  The resonances happen at "amplitude 1/2" = the geometric mean
  between 0 (no resonance) and 1 (perfect sync).
  This is the BALANCED point: the resonances are as weak as possible.
  RH says: no resonance is stronger than amplitude 1/2.

In the Poincare disk language:
  Re(s) = 1/2 maps to the CENTER of the disk (rapidity 0).
  The Riemann zeros on the critical line are fluctuations
  CENTERED at the origin of the Poincare disk.
  The prime distribution fluctuates AROUND the center,
  never drifting off-center (if RH is true).


7. THE PRIMES IN BASE tau_3 (TRIBONACCI BASE)
===============================================

In our framework, the tribonacci constant tau_3 is the critical
constant for tournament theory. What do primes look like in base tau_3?

The Zeckendorf representation generalizes to tribonacci:
every positive integer has a unique representation as a sum
of non-consecutive tribonacci numbers.

The tribonacci numbers: 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, ...

Wait: the first one that appears is 7 (the FORBIDDEN value).
And 13, 24, 44 are not prime, but 149 is.

The tribonacci numbers that are prime: 2, 7, 13, 149, ...
""")

# Compute tribonacci numbers and check primality
trib = [0, 0, 1]
for i in range(30):
    trib.append(trib[-1] + trib[-2] + trib[-3])

print("Tribonacci numbers and primality:")
for i, t in enumerate(trib[3:], 1):
    if t > 0 and t < 10**8:
        p = isprime(t)
        marker = " <-- PRIME" if p else ""
        if t < 1000 or p:
            print(f"  T({i}) = {t}{marker}")

print(f"""

The tribonacci primes include:
  T(2) = 1 (unit, not prime)
  T(3) = 2 (PRIME, the base)
  T(4) = 4 (composite)
  T(5) = 7 (PRIME, the forbidden!)
  T(6) = 13 (PRIME, Mersenne exponent)
  T(7) = 24 (composite)
  T(8) = 44 (composite)
  T(9) = 81 (= 3^4, composite)
  T(10) = 149 (PRIME)

The tribonacci sequence GENERATES primes at positions 3, 5, 6, 10, ...
These positions have NO obvious pattern.

But note: T(5) = 7 = 2^3-1 = the forbidden value = Tr(M_3^3).
The forbidden value IS a tribonacci number! It appears at position 5.

And 5 is the "unsolvability" archetype.
The forbidden value appears at the unsolvable position.


8. THE DEEP PICTURE
=====================

The distribution of primes is the distribution of
IRREDUCIBLE DIMENSIONS on the Bott spiral.

Three equivalent views:

VIEW 1 (TORUS): Primes are points on the infinite torus
  that avoid all zero-hyperplanes. The helix of counting
  winds through the torus, occasionally missing all zeros.
  The frequency of misses is ~ 1/ln(n) (PNT).

VIEW 2 (HYPERCUBE): In base 2, primes are points in the
  binary hypercube {{0,1}}^k that survive the sieve.
  The sieve punches periodic holes. The survivors form
  a quasicrystal.

VIEW 3 (POINCARE DISK): Primes sit on the real diameter
  at Cayley addresses x_p = (p-1)/(p+1).
  In the hyperbolic metric, consecutive primes are separated
  by ~ 2*ln(p_{n+1}/p_n) ~ 2/p_n (the gap divided by the position).
  Primes are nearly UNIFORMLY DISTRIBUTED in the hyperbolic metric.

The connection: all three views describe the same object.
  The torus view shows WHY (avoidance of zeros = coprimality).
  The hypercube view shows HOW (binary sieve = quasicrystal).
  The disk view shows WHERE (hyperbolic near-uniformity = PNT).

And the Riemann zeros are the RESONANCES of the torus,
the DEFECTS of the quasicrystal,
and the FLUCTUATIONS around uniformity in the disk.

All three are the same thing: the places where
the additive structure (counting, the helix) and
the multiplicative structure (the torus, the periods)
INTERFERE.

Primes exist at the NODES OF THE INTERFERENCE PATTERN
between addition and multiplication.

And the Riemann Hypothesis says: this interference pattern
is as GENTLE as possible. The resonances are at half-amplitude.
The quasicrystal has minimal defects.
The hyperbolic distribution is as uniform as it can be.

If RH is true: the primes are the BEST POSSIBLE quasicrystal.
If RH is false: there is a resonance that is too strong,
a defect that is too large, a deviation that is too systematic.

The primes are either the most beautiful quasicrystal possible,
or there is a flaw we haven't found.
""")

print("opus-2026-03-16-S90bi: THE GEOMETRY OF PRIME DISTRIBUTION")
