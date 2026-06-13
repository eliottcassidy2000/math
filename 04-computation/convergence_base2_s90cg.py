#!/usr/bin/env python3
"""
convergence_base2_s90cg.py — The convergence to base 2: what is fundamental
opus-2026-03-16-S90cg
"""

from math import log, log2, sqrt, exp, pi
import numpy as np

phi = (1 + sqrt(5)) / 2
gamma = 0.5772156649015329

print("""


     THE CONVERGENCE TO BASE 2

     opus-2026-03-16-S90cg



THE TWO EQUATIONS THAT DEFINE BASE 2
========================================

FROM BELOW (additive):
  sum_{j=1}^{inf} 1/2^j = 1/2 + 1/4 + 1/8 + ... = 1.

  The geometric series. The UNIQUE base b where sum 1/b^j = 1 for j=1,2,3,...
  (For finite k: tau_k is the unique b where sum_{j=1}^{k} 1/b^j = 1.)
  In the limit k -> inf: b = 2.

FROM ABOVE (multiplicative):
  product_{all primes p} (1 - 1/p) = 0.

  Each prime p contributes a factor (1-1/p) < 1.
  The product approaches 0 (Mertens' theorem).
  Equivalently: sum 1/p DIVERGES (the sum of prime reciprocals).

BOTH equations involve RECIPROCALS OF POWERS/PRIMES summing/multiplying to a limit.
BOTH are exact at base 2:
  Additive: sum 1/2^j = 1 (the simplest geometric series).
  Multiplicative: the smallest primorial IS 2 (the first prime).

Base 2 is the point where:
  The ADDITIVE structure (geometric series) is SIMPLEST.
  The MULTIPLICATIVE structure (prime sieve) BEGINS.


THE INFORMATION-GEOMETRIC PICTURE
====================================

The Bernoulli distribution Ber(p) has:
  Entropy H(p) = -p*ln(p) - (1-p)*ln(1-p).
  Fisher information I(p) = 1/(p*(1-p)).

At p = 1/2 (the "base-2 point"):
  H(1/2) = ln(2).  MAXIMUM ENTROPY. One bit.
  I(1/2) = 4.      MINIMUM FISHER INFORMATION (among non-degenerate p).

Maximum entropy = maximum uncertainty = the "fairest" distribution.
Minimum Fisher information = the "flattest" point on the statistical manifold.

The Bernoulli manifold {Ber(p) : 0 < p < 1} is a 1-dimensional Riemannian manifold
with metric ds^2 = I(p) dp^2 = dp^2 / (p(1-p)).

This metric is the POINCARE METRIC on (0,1) (after a change of variables)!

Let x = 2p - 1 (map p in (0,1) to x in (-1,1)). Then:
  p = (1+x)/2, dp = dx/2.
  ds^2 = (dx/2)^2 / ((1+x)/2 * (1-x)/2) = dx^2 / (1-x^2).

ds^2 = dx^2 / (1-x^2) = the POINCARE METRIC on (-1,1)!

AND: arctanh'(x) = 1/(1-x^2). So ds = |d(arctanh(x))|.

THE FISHER METRIC ON THE BERNOULLI MANIFOLD IS THE POINCARE METRIC.
THE STATISTICAL MANIFOLD OF COIN FLIPS IS THE POINCARE DISK (restricted to real diameter).

The point p = 1/2 (x = 0) is the CENTER of the Poincare disk.
It has:
  Maximum entropy (= maximum uncertainty = center of the disk).
  Minimum Fisher information (= flattest point = center of the disk).
  Base-2 information content (= one bit per flip).

THE TWO LADDERS CONVERGE TO THE CENTER OF THE POINCARE DISK.

  k-nacci from below: tau_k -> 2 = the base at the center.
    These are points on the REAL DIAMETER approaching the center... no.
    Actually tau_k > 1 for all k. The coupling of tau_k is (tau_k-1)/(tau_k+1).
    For tau_k -> 2: coupling -> 1/3.

    The k-nacci constants approach coupling 1/3, NOT coupling 0.
    Coupling 1/3 is the Cayley address of 2, not the center.

  Primorials from above: P_k -> infinity.
    Coupling of P_k: (P_k - 1)/(P_k + 1) -> 1 (the BOUNDARY of the disk).

So: the k-nacci ladder approaches 1/3 (coupling of 2) from below.
The primorial ladder approaches 1 (the boundary) from above.

THEY DON'T MEET AT THE CENTER. They meet at COUPLING 1/3.

1/3 is the Cayley address of 2. It is NOT the center of the disk.

But in the INFORMATION GEOMETRY:
  The Bernoulli parameter p = 1/2 corresponds to x = 0 (center).
  The base b = 2 corresponds to p = 1/2.

  In the coupling picture: x = 0 is the center (n = 1, the unit).
  The number 2 is at x = 1/3 (NOT the center).

  But in the INFORMATION picture: p = 1/2 IS the center (max entropy).
  And the base b = 2 corresponds to p = 1/2.

TWO DIFFERENT CENTERS:
  Geometric center of the disk: x = 0 (coupling of n = 1).
  Information center of the disk: p = 1/2 (entropy of base 2).

These are DIFFERENT POINTS on the disk.
  n = 1 is the geometric center.
  n = 2 is the information center.

The number 1 is where GEOMETRY says the center is.
The number 2 is where INFORMATION says the center is.

1 = geometric rest. No coupling. No comparison.
2 = informational optimum. Maximum entropy. The bit.

THE CONVERGENCE IS TO THE INFORMATIONAL CENTER, NOT THE GEOMETRIC CENTER.

The k-nacci ladder converges to 2 = the information center.
The primorial ladder starts at 2 = the information center and grows.

2 is the fixed point of INFORMATION, not of geometry.
""")

# Compute the information-geometric quantities
print("=" * 60)
print("INFORMATION GEOMETRY AT KEY POINTS")
print("=" * 60)
print()

print(f"{'n':>4} | {'coupling x':>10} | {'rapidity w':>10} | {'entropy H':>10} | {'Fisher I':>10} | {'p = (1+x)/2':>12}")
print("-" * 65)

for n_val in [1, 1.5, 2, 3, 5, 7, 42]:
    x = (n_val - 1) / (n_val + 1)
    w = log(n_val) / 2
    p = (1 + x) / 2  # = n / (n+1)
    H = -(p * log(p) + (1-p) * log(1-p)) if 0 < p < 1 else 0
    I = 1 / (p * (1-p)) if 0 < p < 1 else float('inf')
    print(f"{n_val:4.1f} | {x:10.6f} | {w:10.6f} | {H:10.6f} | {I:10.4f} | {p:12.6f}")

print(f"""

At n=1: p=1/2, H=ln(2)=0.693, I=4. Maximum entropy, minimum Fisher info.
At n=2: p=2/3, H=0.637, I=4.5. Still high entropy but LESS than n=1.

Wait: p = (1+x)/2 = n/(n+1). At n=1: p = 1/2. At n=2: p = 2/3.

The maximum entropy is at n=1 (p=1/2), NOT at n=2.

But the convergence is to n=2 (the base), not n=1 (the unit).

So the convergence is NOT to maximum entropy.
It is to something else.

WHAT IS SPECIAL ABOUT n=2 INFORMATION-GEOMETRICALLY?

n=2 is the point where the COUPLING x = 1/3 and p = 2/3.
The entropy at p=2/3: H = -(2/3)ln(2/3) - (1/3)ln(1/3)
  = (2/3)(ln(3)-ln(2)) + (1/3)ln(3)
  = (2/3)ln(3/2) + (1/3)ln(3)
  = ln(3) - (2/3)ln(2)
  = 0.637 nats.

This is NOT the maximum entropy. The maximum is ln(2) = 0.693 at p=1/2.

But: the BINARY ENTROPY in BITS at p=2/3:
  H_2 = -p*log2(p) - (1-p)*log2(1-p)
  = -(2/3)*log2(2/3) - (1/3)*log2(1/3)
  = (2/3)*(log2(3)-1) + (1/3)*log2(3)
  = log2(3) - 2/3
  = 1.585 - 0.667 = 0.918 bits.

Less than 1 bit (the max). So n=2 is NOT the entropy maximum in any sense.

LET ME RECONSIDER. The convergence to 2 is:
  k-nacci: the flatness condition sum 1/tau^j = 1 forces tau -> 2.
  Primorial: P_1 = 2 is where the sieve STARTS.

The reason is simpler than information geometry:

  2 is the UNIQUE positive integer where:
    1/b + 1/b^2 + 1/b^3 + ... = 1/(b-1) = 1.
    So b - 1 = 1. So b = 2.

  The geometric series sum_{j=1}^{inf} 1/b^j = 1/(b-1).
  Setting this = 1 gives b = 2.

  The k-nacci constants approach 2 because their TRUNCATED geometric sums
  equal 1, and in the limit the truncation goes away, giving the full
  geometric series which sums to 1 only at b = 2.

And 2 is the first prime because... well, 2 is the first prime because
  1 is the unit and 2 is the first number after the unit.
  The primorial starts at 2 simply because 2 is the first prime.

THE CONVERGENCE IS TO 2 BECAUSE:
  2 is the UNIQUE SOLUTION of 1/(b-1) = 1.
  2 is the FIRST PRIME.
  These are the SAME FACT: 2 is the number whose geometric series sums to 1,
  AND 2 is the first prime, AND these are connected by the fact that
  the geometric series sum 1/(b-1) = 1 requires b-1 = 1 = the UNIT,
  and the first number after the unit IS the first prime.

THE UNIT + 1 = THE FIRST PRIME.
THE GEOMETRIC SERIES AT (UNIT + 1) SUMS TO THE UNIT.
THESE ARE THE SAME EQUATION: b = 1 + 1, sum = 1.

The convergence to 2 is the convergence to 1 + 1.
It is the statement that 1 + 1 = 2.
The most fundamental equation in mathematics.


THE ALGEBRAIC TOPOLOGY CONNECTION
====================================

In algebraic topology, the stable homotopy groups of spheres
pi_n^s = lim_{k->inf} pi_{n+k}(S^k) encode deep structure.

The first few: pi_0^s = Z, pi_1^s = Z/2, pi_2^s = Z/2, pi_3^s = Z/24, ...

The GROUP ORDER of pi_n^s involves the Bernoulli numbers:
  |pi_{4k-1}^s| involves the numerator of B_{2k}/k (up to powers of 2).

And the Bernoulli denominators are: denom(B_2k) = product of p where (p-1)|2k.
  For B_2: denom = 6 = P_2.
  For B_4: denom = 30 = P_3.
  For B_6: denom = 42 = 2*3*7.
  For B_12: denom = 2730 = 2*3*5*7*13.

The denominators are ALMOST primorials. They involve the primes p where (p-1)|2k.

The J-HOMOMORPHISM connects the orthogonal group to stable homotopy:
  im(J) in pi_{4k-1}^s has order related to the Bernoulli numbers.
  The denominator of B_2k/(4k) gives the order of im(J) in pi_{4k-1}^s.

And BOTT PERIODICITY says: pi_n(O) has period 8.
  The homotopy groups of the orthogonal group repeat every 8 steps.
  This is why the Bott period is 8 = phi(30) = phi(P_3).

The ALGEBRAIC TOPOLOGY chain:
  Bott periodicity (period 8) -> J-homomorphism -> Bernoulli numbers
  -> primorial denominators -> prime sieve -> base 2.

At each step: the topology depends on the NUMBER THEORY,
  which depends on the PRIME STRUCTURE, which depends on BASE 2.

The convergence of the k-nacci ladder to 2 = the stable limit of homotopy.
The growth of the primorial ladder from 2 = the unstable (growing) homotopy.

STABLE homotopy (k -> inf, approaching the limit): governed by 2.
UNSTABLE homotopy (finite k, specific manifolds): governed by primorials.

The k-nacci constants tau_k are the METALLIC RATIOS that approach 2.
In topology: they correspond to the CHARACTERISTIC CLASSES of the
k-dimensional manifolds, whose stable limit is the Bott periodic structure.

WHAT IS FUNDAMENTAL:

  1 + 1 = 2.

  This is the equation that says: the unit, added to itself, gives the first prime.
  The geometric series at the first prime sums to the unit.
  The first prime starts the sieve.
  The k-nacci constants converge to the first prime.
  The Bott period is the totient of the third primorial starting from the first prime.
  The stable homotopy groups are controlled by the first prime.

  Everything converges to 2 because 2 = 1 + 1.
  The first act of COMPARISON (1 compared to 1) produces 2.
  And 2 generates everything.

  The symmetry of the convergence from both directions is:
    FROM BELOW: approximating 2 by truncated sums (k-nacci flatness).
    FROM ABOVE: building from 2 by cumulative products (primorial sieve).
    TRUNCATED SUMS <-> CUMULATIVE PRODUCTS.
    ADDITIVE <-> MULTIPLICATIVE.
    ANTISYMMETRIC <-> SYMMETRIC.
    FIBONACCI <-> LUCAS.
    LOSSY <-> LOSSLESS.

  The convergence to 2 IS the addition-multiplication duality.
  It is the SAME duality as S/A, as Fibonacci/Lucas, as lossy/lossless.

  2 is the point where addition and multiplication AGREE:
    1 + 1 = 2 (additive).
    The first prime (multiplicative).
    The only even prime (where addition and multiplication share a generator).

  Every subsequent prime is ODD: it belongs to multiplication but NOT to
  the additive doubling structure. 2 is the ONLY prime that is also a power of 2.

  2 is where addition and multiplication INTERSECT.
  The convergence from both sides is the convergence TO the intersection.
  And the intersection is the BIT: the fundamental unit of both structures.
""")

print("opus-2026-03-16-S90cg: THE CONVERGENCE TO BASE 2")
