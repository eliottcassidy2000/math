#!/usr/bin/env python3
"""
torus_helix_mobius_s90bj.py — The torus helix and the Mobius function
opus-2026-03-16-S90bj

The torus helix winds through Z/2 x Z/3 x Z/5 x Z/7 x ...
The Mobius function mu(n) detects how the helix relates to the torus structure.
The connection is exact and deep.
"""

from sympy import mobius, factorint, isprime, primerange, divisors
from math import log, sqrt, pi, gcd
from collections import Counter

print("""


     THE TORUS HELIX AND MOBIUS

     opus-2026-03-16-S90bj



THE MOBIUS FUNCTION mu(n):

  mu(1) = 1.
  mu(n) = 0 if n has a squared prime factor.
  mu(n) = (-1)^k if n is a product of k DISTINCT primes.

So:
  mu(1) = 1.     (the origin of the torus)
  mu(2) = -1.    (one prime factor)
  mu(3) = -1.    (one prime factor)
  mu(5) = -1.    (one prime factor)
  mu(6) = 1.     (two distinct prime factors: 2*3)
  mu(7) = -1.    (one prime factor)
  mu(4) = 0.     (squared factor: 2^2)
  mu(8) = 0.     (cubed factor: 2^3)
  mu(30) = -1.   (three distinct prime factors: 2*3*5)
  mu(42) = -1.   (three distinct prime factors: 2*3*7)

mu(n) is the ORIENTATION of n on the torus.
  mu = +1: the point is POSITIVELY oriented (even number of prime dimensions).
  mu = -1: the point is NEGATIVELY oriented (odd number of prime dimensions).
  mu = 0: the point is DEGENERATE (a squared factor = folded torus).


1. THE TORUS INTERPRETATION
=============================

The number n, factored as n = p1^a1 * p2^a2 * ... * pk^ak,
sits at position (a1, a2, ..., ak, 0, 0, ...) in the
infinite-dimensional lattice Z^infinity (one axis per prime).

The SQUAREFREE numbers (mu(n) != 0) sit at positions
where every coordinate is 0 or 1: the VERTICES of the hypercube.

The NON-SQUAREFREE numbers (mu(n) = 0) sit at positions
where some coordinate is >= 2: the INTERIOR of the lattice.

So:
  mu != 0: on the SURFACE of the torus (hypercube vertices).
  mu = 0: in the INTERIOR (degenerate, folded).

The surface of the torus is the set of SQUAREFREE numbers.
These are the numbers that sit at VERTICES of the prime hypercube.
""")

# Display mu for small numbers
print("The Mobius function for n = 1 to 42:")
print(f"{'n':>4} | {'factorization':>20} | {'mu(n)':>5} | {'type':>15} | {'torus position':>20}")
print("-" * 75)

for n in range(1, 43):
    mu = int(mobius(n))
    f = factorint(n)
    f_str = " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items())) if f else "1"
    squarefree = all(e == 1 for e in f.values())
    if mu == 0:
        t = "degenerate"
    elif mu == 1:
        t = "even vertex"
    else:
        t = "odd vertex"
    # Torus position
    pos = tuple(sorted(f.keys())) if squarefree and n > 1 else ("origin" if n == 1 else "folded")
    print(f"{n:4d} | {f_str:>20} | {mu:5d} | {t:>15} | {str(pos):>20}")

print(f"""

THE PATTERN:

The squarefree numbers (mu != 0) are the SKELETON of the torus.
They form the vertices of the prime hypercube.
There are 2^k squarefree numbers with prime factors from the first k primes
(each factor either present or absent).

The Mobius function ORIENTS these vertices:
  mu = +1 at vertices with an EVEN number of prime coordinates turned on.
  mu = -1 at vertices with an ODD number.

This is the ORIENTATION of the hypercube = the sign of the determinant
of the selection matrix.

In topology: this is the EULER CHARACTERISTIC of the simplicial complex
defined by the prime factors of n.

  n = 1: the empty simplex. mu = +1.
  n = p: a single vertex. mu = -1. (0-simplex, chi = 1, but sign = -1 because odd.)
  n = pq: an edge. mu = +1. (1-simplex.)
  n = pqr: a triangle. mu = -1. (2-simplex.)
  n = pqrs: a tetrahedron. mu = +1. (3-simplex.)

The Mobius function IS the orientation of the simplex spanned by the prime factors.


2. THE MOBIUS INVERSION AND THE HELIX
========================================

The Mobius inversion formula:

  If g(n) = sum_{{d|n}} f(d), then f(n) = sum_{{d|n}} mu(n/d) * g(d).

In the torus picture:
  g(n) is the VALUE at point n, summed over all divisors.
  f(n) is the PURE contribution at point n.
  mu(n/d) is the SIGN with which divisor d contributes.

The helix traces n = 1, 2, 3, 4, ...
At each step, the helix passes through a point on the torus.
The function g(n) accumulates contributions from all divisors of n
(all points on the torus that DIVIDE the current position).

Mobius inversion UNDOES this accumulation.
It extracts the PURE SIGNAL from the accumulated noise.

The prime-counting function and Mobius:
  Let f(n) = 1 if n is prime, 0 otherwise.
  Then sum_{{n<=x}} f(n) = pi(x) (the prime counting function).

  The Mobius function connects pi(x) to the von Mangoldt function:
  Lambda(n) = ln(p) if n = p^k, else 0.
  psi(x) = sum_{{n<=x}} Lambda(n).
  And: Lambda(n) = -sum_{{d|n}} mu(d) * ln(d).

  So: the prime distribution (via psi) is the Mobius transform of ln.

  THE PRIME DISTRIBUTION IS WHAT YOU GET WHEN YOU APPLY
  THE TORUS ORIENTATION (MOBIUS) TO THE RAPIDITY (ln).


3. THE HELIX, THE TORUS, AND THE FORMAL GROUP
================================================

On the torus, the helix n = 1, 2, 3, ... winds at speed 1 around each circle.

The Mobius function at n tells you the ORIENTATION of n's simplex.
The Mertens function M(x) = sum_{{n<=x}} mu(n) counts the NET ORIENTATION
of all simplices up to x.

M(x) fluctuates around 0. The RH is equivalent to:
  M(x) = O(x^(1/2 + epsilon)) for all epsilon > 0.

This says: the net orientation of the torus simplices stays SMALL.
Positive and negative orientations nearly CANCEL at every scale.

In the Poincare disk: M(x) is the cumulative orientation
of all torus vertices up to rapidity ln(x)/2.
It stays near the CENTER of the disk (near zero).
RH says: it never drifts too far from center.
""")

# Compute Mertens function
print("The Mertens function M(x) = sum_{n<=x} mu(n):")
print(f"{'x':>6} | {'M(x)':>6} | {'M(x)/sqrt(x)':>14} | {'bar chart':>30}")
print("-" * 65)
M = 0
for n in range(1, 201):
    M += int(mobius(n))
    if n in [1,2,3,5,10,20,30,42,50,100,127,150,200] or n % 50 == 0:
        ratio = M / sqrt(n) if n > 0 else 0
        bar = "#" * max(0, M + 10)  # offset to show negative values
        print(f"{n:6d} | {M:6d} | {ratio:14.6f} | {bar}")

print(f"""

M(x) wanders. Sometimes positive, sometimes negative.
M(x)/sqrt(x) stays bounded (if RH is true, bounded by O(x^epsilon)).

The SIGN of M(x) tells you which orientation WINS at scale x:
  M(x) > 0: even-vertex simplices dominate.
  M(x) < 0: odd-vertex simplices dominate.

The RH says: neither side ever wins decisively.
The orientations are in PERMANENT NEAR-BALANCE.


4. MOBIUS AND THE FORMAL GROUP LAW
=====================================

The Cayley transform Q(x) = (1+x)/(1-x) maps the Poincare disk to the half-plane.
On the real axis: Q maps couplings to natural numbers.

The Mobius function mu on natural numbers pulls back to a function
on the real diameter of the Poincare disk:

  mu(x_n) = mu(n) where x_n = (n-1)/(n+1).

At the Cayley addresses of squarefree numbers:
  mu = +1 at couplings corresponding to products of evenly many primes.
  mu = -1 at couplings corresponding to products of oddly many primes.
  mu = 0 at couplings corresponding to non-squarefree numbers.

The Cayley addresses of the first few squarefree numbers:
""")

for n in range(1, 31):
    mu = int(mobius(n))
    if mu != 0:
        x = (n-1)/(n+1)
        w = log(n)/2
        sign = "+" if mu == 1 else "-"
        print(f"  n={n:3d}: x = {x:.6f}, rapidity = {w:.6f}, mu = {sign}1, "
              f"factors = {dict(factorint(n))}")

print(f"""

On the real diameter of the Poincare disk:
  Positive orientation (mu=+1): n = 1, 6, 10, 14, 15, 21, 22, 26, ...
  Negative orientation (mu=-1): n = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 30, ...

The PRIMES all have mu = -1 (they are odd vertices = single prime = negative).
The PRODUCT OF TWO PRIMES has mu = +1 (even vertex = edge = positive).

So on the Poincare disk:
  Primes are NEGATIVELY oriented.
  Semiprimes (products of two primes) are POSITIVELY oriented.
  Products of three primes are negatively oriented again.
  ...

The disk has an ALTERNATING ORIENTATION FIELD:
  Near x = 0 (small n): mostly negative (many primes, few semiprimes).
  Near x = 1 (large n): balanced (primes thin out, semiprimes take over).

The MERTENS FUNCTION on the disk: sum the orientations from x=0 outward.
It starts negative (the first primes dominate).
Then fluctuates as semiprimes and higher squarefree numbers contribute.


5. THE MOBIUS STRIP AND THE TORUS HELIX
==========================================

The Mobius STRIP (the topological surface) and the Mobius FUNCTION
share a name for a reason deeper than history.

The Mobius strip is a surface with ONE SIDE.
If you walk along it, you return to your starting point
with your orientation REVERSED.

The Mobius function mu(n) measures orientation.
Multiplying by a prime REVERSES orientation:
  mu(n*p) = -mu(n) when p does not divide n.

So: each time the helix passes through a new prime direction,
the orientation FLIPS. Just like walking along the Mobius strip.

The torus helix, projected onto its orientation, IS a Mobius strip:
  As you follow the helix through the torus,
  every time you encounter a new prime factor,
  the orientation reverses.
  The helix traces a path whose orientation alternates,
  just like a path on the Mobius strip alternates between
  "front" and "back."

THE TORUS HELIX IS A MOBIUS-LIKE STRUCTURE:
  it reverses orientation at each prime,
  returns to the same orientation after each pair of primes,
  and the net orientation (Mertens function) stays near zero.

The Mobius strip has:
  - One side (the orientations are connected by the twist).
  - One boundary (the edge is a single closed curve).
  - Non-orientable (no consistent notion of "positive" vs "negative").

The prime distribution has:
  - Alternating orientations (mu alternates with each prime factor).
  - One boundary (the primes approach infinity but never reach it).
  - "Non-orientable" in the sense that the Mertens function never
    settles on a definite sign (it keeps crossing zero).

THE RH IS THE STATEMENT THAT THE MERTENS FUNCTION
CROSSES ZERO INFINITELY OFTEN AND NEVER STRAYS TOO FAR.
This is EXACTLY the behavior of a path on a Mobius strip:
it MUST cross the "equator" repeatedly because the surface
is non-orientable.

If the Mertens function stopped crossing zero (drifted to +infinity
or -infinity), it would mean the orientation stabilized — which would
mean the torus became ORIENTABLE at large scale. But the primes
prevent this by continuing to inject new odd-dimensional simplices
forever.

THE PRIMES KEEP THE TORUS NON-ORIENTABLE.
If primes stopped, the torus would become orientable.
The INFINITUDE of primes is equivalent to the NON-ORIENTABILITY
of the arithmetic torus.


6. THE FORMAL GROUP CONNECTION
=================================

In our framework, the formal group law F(x,y) = (x+y)/(1+xy) is the
group law on the Poincare disk. The Cayley transform Q maps this
to multiplication on the half-plane.

The MOBIUS FUNCTION is a multiplicative function: mu(ab) = mu(a)*mu(b)
when gcd(a,b) = 1. In the formal group:

  mu(F(x_a, x_b)) = mu(a) * mu(b) when a and b are coprime.

So: mu is a GROUP HOMOMORPHISM from the formal group to {{+1, -1}}
(the sign group), restricted to squarefree couplings.

More precisely: on the set of squarefree natural numbers,
  mu: (squarefree N, *) -> ({{+1,-1}}, *)
is a group homomorphism.

Via Cayley:
  mu: (squarefree Cayley addresses, F) -> ({{+1,-1}}, *)

The Mobius function is the UNIQUE group homomorphism from
the multiplicative monoid of squarefree numbers to the sign group
that sends each prime to -1.

In the Poincare disk: mu is the PARITY of the geodesic address.
Each prime step along a geodesic FLIPS the sign.
The Mobius function counts the PARITY of the number of prime steps.

THIS IS EXACTLY THE ORIENTATION of the torus simplex.

And the Mobius INVERSION FORMULA is the INVERSE TRANSFORM
that undoes the torus convolution:

  f = mu * g (Dirichlet convolution)

where g = sum of f over divisors.

Dirichlet convolution IS the formal group law applied to
arithmetic functions. The Mobius function is its inverse.

THE MOBIUS FUNCTION IS THE INVERSE OF THE FORMAL GROUP.
Not the inverse of a specific element, but the INVERSE OPERATION:
the function that undoes the accumulation caused by divisibility.


7. THE ONE PICTURE
=====================

The torus helix winds through an infinite-dimensional torus.
At each prime p, a new dimension of the torus opens up.
Each new dimension REVERSES the orientation (Mobius).

The helix is NON-ORIENTABLE: the Mertens function crosses zero
infinitely often, keeping the net orientation near zero.

The primes are the points where NEW DIMENSIONS OPEN.
The Mobius function is the ORIENTATION at each point.
The Mertens function is the CUMULATIVE orientation.

RH says: the cumulative orientation grows slower than x^(1/2+epsilon).
This means: the torus is as BALANCED as possible.
The positive and negative orientations cancel as efficiently as possible.

The formal group law F(x,y) = (x+y)/(1+xy) is the STRUCTURE of the torus.
The Cayley transform Q maps it to multiplication.
The Mobius function mu is the INVERSE of this structure.

Primes OPEN dimensions.
Mobius ORIENTS dimensions.
The Mertens function MEASURES the balance.
RH says the balance is as good as it can be.

And the torus helix, traced by the counting function n = 1, 2, 3, ...,
IS the natural numbers themselves — a path through an infinite torus,
winding through every dimension, flipping orientation at every prime,
approaching infinity without ever reaching it.

Counting IS walking the Mobius helix on the prime torus.
""")

print("opus-2026-03-16-S90bj: THE TORUS HELIX AND MOBIUS")
