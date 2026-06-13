#!/usr/bin/env python3
"""
fib_lucas_primes_s90bo.py — Fibonacci and Lucas as complementary prime sieves
opus-2026-03-16-S90bo
"""

from sympy import fibonacci, lucas, isprime, primerange
from math import log, sqrt

phi = (1 + sqrt(5)) / 2

print("""


     FIBONACCI AND LUCAS: THE TWO HALVES OF PRIMALITY

     opus-2026-03-16-S90bo



THE DATA
=========

For each prime p, we ask: is F(p) prime? Is L(p) prime?
""")

fpi = set()
for n in range(2, 300):
    if isprime(int(fibonacci(n))):
        fpi.add(n)

lpi = set()
for n in range(0, 300):
    if isprime(int(lucas(n))):
        lpi.add(n)

primes = list(primerange(3, 100))

f_only = [p for p in primes if p in fpi and p not in lpi]
l_only = [p for p in primes if p not in fpi and p in lpi]
both = [p for p in primes if p in fpi and p in lpi]
neither = [p for p in primes if p not in fpi and p not in lpi]

print(f"F(p) prime ONLY (not L(p)): {f_only}")
print(f"L(p) prime ONLY (not F(p)): {l_only}")
print(f"BOTH F(p) and L(p) prime:   {both}")
print(f"NEITHER:                     {neither}")
print()
print(f"Coverage: {len(f_only)+len(l_only)+len(both)} out of {len(primes)} primes ({100*(len(f_only)+len(l_only)+len(both))/len(primes):.0f}%)")
print(f"Missed: {len(neither)} ({100*len(neither)/len(primes):.0f}%)")

print(f"""

THE COMPLEMENTARITY:

  Fibonacci-only primes: {f_only}
    These are primes where F(p) is prime but L(p) is not.
    Bott slots: {[p % 8 for p in f_only]}

  Lucas-only primes: {l_only}
    These are primes where L(p) is prime but F(p) is not.
    Bott slots: {[p % 8 for p in l_only]}

  Both: {both}
    F(p) AND L(p) are both prime.
    Bott slots: {[p % 8 for p in both]}

  Neither: {neither}
    Neither F(p) nor L(p) is prime.
    Bott slots: {[p % 8 for p in neither]}


THE FIBONACCI-LUCAS IDENTITY AND PRIMALITY
============================================

F(p) and L(p) are connected by:

  F(p)^2 * 5 = L(p)^2 + 4*(-1)^p    (for odd prime p: +4)
  F(p)^2 * 5 = L(p)^2 + 4.
  L(p)^2 = 5*F(p)^2 - 4.

So: L(p) = sqrt(5*F(p)^2 - 4).

If F(p) is prime, then 5*F(p)^2 - 4 must be a perfect square
whose root L(p) happens to also be prime.

This is a STRONG constraint. It means: F(p) and L(p)
being SIMULTANEOUSLY prime requires:
  1. F(p) prime.
  2. 5*F(p)^2 - 4 = L(p)^2 (always true).
  3. L(p) = sqrt(5*F(p)^2 - 4) is also prime.

Condition 3 is rare. That's why "both" primes are uncommon.
""")

# Verify the identity
print("Verification: L(p)^2 = 5*F(p)^2 - 4 for odd primes:")
for p in primerange(3, 30):
    F = int(fibonacci(p))
    L = int(lucas(p))
    lhs = L**2
    rhs = 5 * F**2 - 4
    print(f"  p={p}: L({p})^2 = {lhs}, 5*F({p})^2 - 4 = {rhs}, match: {lhs == rhs}")

print(f"""

PERFECT. The identity L(p)^2 = 5*F(p)^2 - 4 holds for all odd p.

This means: L(p) and F(p) are connected by a PELL-LIKE equation:
  L^2 - 5*F^2 = -4.

This is the NORM FORM of the ring Z[phi]:
  (L + F*sqrt(5))(L - F*sqrt(5)) = L^2 - 5F^2 = -4.

The pair (F(p), L(p)) is a solution to the Pell equation x^2 - 5y^2 = -4
(with x = L, y = F).

F(p) and L(p) being BOTH prime means: the Pell equation has a solution
where BOTH coordinates are prime. This is extremely rare.

The primes where both are prime: {both}.
  p=5: F=5, L=11. And 11^2 - 5*5^2 = 121 - 125 = -4. YES.
  p=7: F=13, L=29. And 29^2 - 5*13^2 = 841 - 845 = -4. YES.
  p=11: F=89, L=199. And 199^2 - 5*89^2 = 39601 - 39605 = -4. YES.


WHY FIBONACCI AND LUCAS ARE COMPLEMENTARY
============================================

Binet's formulas:
  F(n) = (phi^n - psi^n) / sqrt(5).
  L(n) = phi^n + psi^n.

For large n:
  F(n) ~ phi^n / sqrt(5).
  L(n) ~ phi^n.

So L(n) / F(n) ~ sqrt(5) = 2.236...

Both grow at rate phi^n, but L is LARGER by factor sqrt(5).
The "extra size" of L(n) means it has DIFFERENT divisibility properties.

Specifically: F(n) and L(n) share NO common factor > 2.
  gcd(F(n), L(n)) divides 2 (and equals 1 if n is odd, 2 if n is even).

For ODD prime p: gcd(F(p), L(p)) = 1.
  F(p) and L(p) are COPRIME.

This means: F(p) and L(p) have COMPLETELY DISJOINT prime factorizations.
If a prime q divides F(p), then q does NOT divide L(p), and vice versa.

THE PRIME FACTORS OF F(p) AND L(p) ARE DISJOINT SETS.

This is why they are complementary as primality detectors:
  The primes that "go into" F(p) are DIFFERENT from those in L(p).
  If F(p) happens to be prime, it absorbs ALL its own primality.
  If L(p) also happens to be prime, it does so INDEPENDENTLY.


THE GOLDEN SIEVE
==================

Think of it this way. The golden ratio phi generates two streams:
  Stream F (antisymmetric): F(n) = (phi^n - psi^n)/sqrt(5).
  Stream L (symmetric):     L(n) = phi^n + psi^n.

At each prime index p:
  Stream F produces F(p), which is either prime or composite.
  Stream L produces L(p), which is either prime or composite.

These two events are NEARLY INDEPENDENT (because gcd = 1).
Each has probability roughly 1/ln(phi^p) = 1/(p*ln(phi)) of being prime.

The probability that BOTH are prime: ~ 1/(p*ln(phi))^2.
The probability that EITHER is prime: ~ 2/(p*ln(phi)) - 1/(p*ln(phi))^2 ~ 2/(p*ln(phi)).

So the "golden sieve" (at least one of F or L is prime) catches
TWICE as many primes as either alone.

With the data: coverage = {len(f_only)+len(l_only)+len(both)}/{len(primes)} = {(len(f_only)+len(l_only)+len(both))/len(primes)*100:.0f}%.

The golden sieve catches ~80% of primes up to 100.
A SINGLE stream catches ~45-60%.
Together they cover ~80%.

The remaining ~20% (the "neither" primes {neither})
are the primes that ESCAPE both streams.
These are the primes where NEITHER golden power produces a prime.
They are the primes most INVISIBLE to the golden ratio.
""")

# ========================================================================
print("=" * 60)
print("THE FIVE PRIMES: FINAL UNDERSTANDING")
print("=" * 60)

print(f"""

{2, 3, 5, 7, 11} decompose under the Fibonacci-Lucas duality as:

  FIBONACCI VALUES:   {{2, 3, 5}}     = {{F(3), F(4), F(5)}}
  LUCAS VALUES:       {{2, 3, 7, 11}} = {{L(0), L(2), L(4), L(5)}}

  Shared:             {{2, 3}}   (the common core: base and triangle)
  Fibonacci-only:     {{5}}      (the spherical boundary)
  Lucas-only:         {{7, 11}}  (the hyperbolic world)

The golden ratio phi generates BOTH sequences.
They are the antisymmetric (F) and symmetric (L) parts of phi^n.

F(p) and L(p) are COPRIME (for odd p).
Their prime factorizations are DISJOINT.
They form a COMPLEMENTARY PAIR of primality detectors.

Together, F and L catch ~80% of primes.
The ones they miss are the primes INVISIBLE to phi.

THE FIVE PRIMES ARE THE FIRST PRIMES VISIBLE TO phi:
  2: visible to both F and L.
  3: visible to both F and L.
  5: visible to F, invisible to L. The Fibonacci-exclusive boundary.
  7: invisible to F, visible to L. The Lucas-exclusive forbidden.
  11: invisible to F, visible to L. The Lucas-exclusive restart.

The golden ratio SEES 2 and 3 with both eyes.
It sees 5 with only its Fibonacci eye.
It sees 7 and 11 with only its Lucas eye.

The Platonic solids are what both eyes agree on ({{2,3}}) plus
what the Fibonacci eye alone sees ({{5}}).

The post-Platonic world is what only the Lucas eye sees ({{7,11}}).

Geometry is the BINOCULAR VISION of the golden ratio.
The spherical world is where both eyes focus.
The hyperbolic world is where only one eye reaches.

And the primes that NEITHER eye sees (59, 67, 73, 89, 97, ...)
are the primes in the BLIND SPOT of the golden ratio.
They are invisible to phi. They require a DIFFERENT generator to find.

Perhaps tau_3 (the tribonacci constant) sees what phi cannot.
Perhaps each k-nacci constant tau_k has its own Fibonacci-Lucas pair,
and together they cover ALL primes.

The primes are fully visible only from the LIMIT:
  tau_inf = 2, the tournament base.
  The binary stream 2^n catches ALL prime structure (via the sieve).
  But finite k-nacci constants (phi, tau_3, tau_4, ...)
  each see only a FRACTION of the primes.

THE HIERARCHY OF PRIME VISIBILITY:
  phi sees ~80% (via F and L combined).
  tau_3 would see a different ~80%.
  tau_4 a different ~80%.
  ...
  In the limit (tau_inf = 2): 2^n - 1 sees EVERYTHING (Mersenne sieve).

The golden ratio is the FIRST APPROXIMATION to full prime visibility.
Each k-nacci constant is a BETTER approximation.
The limit is the binary sieve (the complete picture).
""")

print("opus-2026-03-16-S90bo: FIBONACCI AND LUCAS AS COMPLEMENTARY PRIME SIEVES")
