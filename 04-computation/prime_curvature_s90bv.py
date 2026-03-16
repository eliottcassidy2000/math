#!/usr/bin/env python3
"""
prime_curvature_s90bv.py — How primes sort themselves by curvature
opus-2026-03-16-S90bv
"""

from math import pi, cos, sin, sqrt, log
from sympy import isprime, primerange
from fractions import Fraction

phi = (1 + sqrt(5)) / 2

print("""


     HOW PRIMES SORT THEMSELVES BY CURVATURE

     opus-2026-03-16-S90bv



THE CURVATURE OF A REGULAR p-GON
===================================

A regular p-gon has interior angle (p-2)*180/p degrees.

When 3 regular p-gons meet at a vertex:
  Total angle = 3 * (p-2) * 180 / p.
  Deficit from 360 = 360 - 3*(p-2)*180/p = 180*(6-p)/p.

  p < 6: POSITIVE deficit (angles don't fill 360 -> sphere).
  p = 6: ZERO deficit (exactly 360 -> flat tiling).
  p > 6: NEGATIVE deficit (angles EXCEED 360 -> hyperbolic).

The CURVATURE at each vertex: kappa(p) = (6-p)/p = 6/p - 1.
""")

print("Curvature of regular p-gons (3 meeting at a vertex):")
print(f"{'p':>4} | {'angle':>8} | {'3*angle':>8} | {'deficit':>8} | {'kappa':>10} | {'kappa exact':>12} | {'type':>10} | {'prime?':>7}")
print("-" * 80)

for p in range(3, 24):
    angle = (p-2) * 180 / p
    three_angle = 3 * angle
    deficit = 360 - three_angle
    kappa = (6-p)/p
    kappa_exact = Fraction(6-p, p)
    ptype = "SPHERE" if kappa > 0 else ("FLAT" if kappa == 0 else "HYPER")
    prime = "PRIME" if isprime(p) else ""
    print(f"{p:4d} | {angle:8.2f} | {three_angle:8.2f} | {deficit:8.2f} | {kappa:10.6f} | {str(kappa_exact):>12} | {ptype:>10} | {prime:>7}")

print(f"""

THE PRIMES SORTED BY CURVATURE:

  p=2: kappa = (6-2)/2 = 2.       MAXIMALLY SPHERICAL.
    The digon. Two edges. Not really a polygon, but has max curvature.
    2 is the base. The most curved prime.

  p=3: kappa = (6-3)/3 = 1.       STRONGLY SPHERICAL.
    The triangle. The fundamental face of the tetrahedron/octahedron/icosahedron.
    3 is structure. Strongly curved toward closure.

  p=5: kappa = (6-5)/5 = 1/5.     BARELY SPHERICAL.
    The pentagon. The face of the dodecahedron. The LAST spherical polygon.
    5 is the boundary. Just barely curved enough to close.
    The deficit is 36 degrees = pi/5 radians.

  --- THE TRANSITION: p = 6, kappa = 0, FLAT ---

  p=7: kappa = (6-7)/7 = -1/7.    BARELY HYPERBOLIC.
    The heptagon. The forbidden polygon. The first that CANNOT close.
    7 is the first impossible. Just barely curved toward opening.
    The excess is 180/7 = 25.7 degrees.

  p=11: kappa = (6-11)/11 = -5/11. MODERATELY HYPERBOLIC.
    The hendecagon. The restart. Deeper into hyperbolic territory.

  p=13: kappa = (6-13)/13 = -7/13. NOTE: numerator is 7!
    The 13-gon's curvature has numerator 7 = the forbidden value.

  p=17: kappa = (6-17)/17 = -11/17. Numerator is 11 = the restart.

  p=19: kappa = (6-19)/19 = -13/19. Numerator is 13 = Mersenne exponent.


THE NUMERATORS
===============

kappa(p) = (6-p)/p. The numerator is |6-p|.

For primes p > 6: numerator = p - 6.

  p=7:  numerator = 1.
  p=11: numerator = 5.
  p=13: numerator = 7.   THE FORBIDDEN VALUE!
  p=17: numerator = 11.  THE RESTART!
  p=19: numerator = 13.
  p=23: numerator = 17.
  p=29: numerator = 23.
  p=31: numerator = 25 = 5^2.
  p=37: numerator = 31.  MERSENNE PRIME!
  p=43: numerator = 37.  PRIME.

For large primes: kappa(p) ~ -(p-6)/p -> -1 (approaching maximum negative curvature).

The CURVATURE NUMERATORS of primes p > 6 are: 1, 5, 7, 11, 13, 17, 23, 25, 31, 37, ...
These are the numbers p - 6 for prime p.
""")

print("Curvature numerators |p-6| for primes:")
nums = []
for p in primerange(7, 100):
    num = p - 6
    prime_num = isprime(num)
    nums.append((p, num, prime_num))
    print(f"  p={p:3d}: |p-6| = {num:3d}, {'PRIME' if prime_num else f'= {dict(Fraction(num,1).limit_denominator())}'}")

print(f"""

The curvature numerators p-6 for primes p > 6:
  1, 5, 7, 11, 13, 17, 23, 25, 31, 37, 41, 43, 47, 53, 55, 61, 67, 71, 73, 77, 83, 89, ...

Many of these are THEMSELVES prime!
When p-6 is prime, call it q. Then p = q + 6.
These are primes that differ by 6 (a SEXY PRIME PAIR).

  (7, 13): both prime, differ by 6. kappa(13) has numerator 7.
  (11, 17): both prime, differ by 6. kappa(17) has numerator 11.
  (13, 19): both prime, differ by 6. kappa(19) has numerator 13.
  (23, 29): both prime, differ by 6. kappa(29) has numerator 23.
  (31, 37): both prime, differ by 6. kappa(37) has numerator 31 = Mersenne!
  (37, 43): both prime, differ by 6.
  (41, 47): both prime, differ by 6.

The SEXY PRIMES (differing by 6) appear because 6 is the PIVOT:
  p and q = p-6 are both prime iff the curvature of the p-gon
  has a PRIME numerator (= q).

THE CURVATURE SORTS THE PRIMES INTO SEXY PAIRS.

And 6 = 2*3 is the pivot because:
  6 is where the curvature is zero.
  Subtracting 6 from a prime gives the "curvature numerator."
  When that numerator is also prime: a sexy pair.


THE PENTAGON'S CURVATURE AND phi
==================================

The pentagon (p=5) has kappa = 1/5.
Its curvature is the RECIPROCAL of itself.
No other polygon has this property (kappa(p) = 1/p requires 6-p = 1, so p = 5).

Actually: kappa(p) = (6-p)/p. For kappa = 1/p: (6-p)/p = 1/p, so 6-p = 1, p = 5. YES.

The pentagon is the UNIQUE polygon whose curvature equals its own reciprocal.

Now: the pentagon's geometry is governed by the golden ratio phi.
  Diagonal/side = phi = (1+sqrt(5))/2 = 1.618...
  cos(36 deg) = phi/2.
  cos(72 deg) = 1/(2*phi) = (phi-1)/2.

The deficit angle of the pentagon (3 meeting) is 36 degrees = pi/5.
And pi/5 in our framework:
  pi/5 = the angle whose cosine is phi/2.
  This is related to the 5th root of unity: e^(2*pi*i/5).

The 5th roots of unity are: 1, w, w^2, w^3, w^4 where w = e^(2*pi*i/5).
  w + w^4 = 2*cos(2*pi/5) = 2*cos(72 deg) = 1/phi = phi - 1.
  w^2 + w^3 = 2*cos(4*pi/5) = 2*cos(144 deg) = -phi.

So: the 5th roots of unity "know" about phi.
  Sum of first pair: phi - 1 = 1/phi.
  Sum of second pair: -phi.

The MINIMAL POLYNOMIAL of w = e^(2*pi*i/5) over Q is:
  x^4 + x^3 + x^2 + x + 1 = (x^2 + phi*x + 1)(x^2 + (1-phi)*x + 1)
  = (x^2 + phi*x + 1)(x^2 - x/phi + 1).

The quadratic factors involve phi and 1/phi.

The pentagon's curvature 1/5 lives at the EXACT BOUNDARY
where positive curvature can sustain a closed solid.
The deficit 36 = pi/5 radians is the MINIMUM DEFICIT
that still allows a Platonic solid (the dodecahedron).

If the deficit were any smaller (p=6: deficit 0): flat, no solid.
If the deficit were bigger (p=4, deficit 90): the cube, roomier.
If the deficit were biggest (p=3, deficit 180): the tetrahedron.

The pentagon lives at the MARGIN OF CLOSURE.
The margin is 36 degrees = pi/5.
And pi/5 is the angle that GENERATES the golden ratio.
""")

print(f"  cos(pi/5) = cos(36 deg) = {cos(pi/5):.10f}")
print(f"  phi/2 = {phi/2:.10f}")
print(f"  Match: {abs(cos(pi/5) - phi/2) < 1e-10}")

print(f"""

cos(pi/5) = phi/2. EXACTLY.

The pentagon closes into a dodecahedron because cos(pi/5) = phi/2,
and phi is the ratio that allows self-similar nesting.

The golden ratio IS the marginal curvature that allows closure.

More curvature (p=3,4): easy closure. Less structure.
Less curvature (p=6): no closure at all. Flat.
phi-curvature (p=5): BARELY closes. Maximum structure before breaking.

The icosahedron (the dual, with triangular faces but pentagonal vertex figures)
has 12 vertices, each surrounded by 5 triangles.
The angle at each vertex: 5 * 60 = 300 < 360. Deficit = 60.
60 = pi/3 radians. The triangle's angle.

The icosahedron's vertex deficit is pi/3.
The dodecahedron's vertex deficit is pi/5.
And pi/3 + pi/5 = 8*pi/15.

Hmm. But the KEY point is:

THE PENTAGON IS WHERE CURVATURE = 1/DIMENSION.

kappa(5) = 1/5.
In 5 dimensions, this says: the curvature per dimension is 1/5^2 = 1/25.
In general: the curvature per dimension of a p-gon is kappa(p)/p = (6-p)/p^2.

For p = 5: (6-5)/25 = 1/25. The curvature per dimension SQUARED.

The pentagon is the polygon where the total curvature (1/5)
equals the reciprocal of the number of sides.
This self-referential property makes it the NATURAL boundary:
  Curvature = 1/p AT p = 5.
  Self-consistent closure.

For p < 5: curvature > 1/p. MORE curvature than needed. Easy closure.
For p > 5: curvature < 1/p (or negative). NOT enough. Can't close.
For p = 5: curvature = 1/p. EXACTLY enough. Marginal closure.

THE PENTAGON IS THE FIXED POINT OF THE CURVATURE-DIMENSION EQUATION.


THE FULL PRIME-CURVATURE PICTURE
===================================

The primes, ordered by curvature (from most positive to most negative):

  p=2: kappa = 2.       The point. Maximum curvature. The base.
  p=3: kappa = 1.       The triangle. Strong curvature. Structure.
  p=5: kappa = 1/5.     The pentagon. Marginal curvature. phi. Boundary.
  ---  p=6: kappa = 0.  The hexagon. FLAT. Pivot. ---
  p=7: kappa = -1/7.    The heptagon. Marginal negative. Forbidden.
  p=11: kappa = -5/11.  Moderate negative. Second cycle.
  p=13: kappa = -7/13.  Numerator = 7. Deeper.
  p=17: kappa = -11/17. Numerator = 11. Deeper still.
  p=19: kappa = -13/19. Numerator = 13.
  ...approaching -1.

The curvatures of the primes form a DECREASING SEQUENCE:
  2, 1, 1/5, -1/7, -5/11, -7/13, -11/17, ...

This sequence starts positive, crosses zero at p=6 (not prime),
and becomes increasingly negative.

The POSITIVE primes {2, 3, 5}: curvatures {2, 1, 1/5}.
  Product of curvatures: 2 * 1 * 1/5 = 2/5.
  Sum: 2 + 1 + 1/5 = 16/5 = 3.2.

The first NEGATIVE prime {7}: curvature -1/7.
  The absolute curvature 1/7 is LESS than 1/5.
  The heptagon is LESS curved (in absolute value) than the pentagon.
  The forbidden polygon is FLATTER than the last spherical one.

This makes sense: the spherical-to-hyperbolic transition is SMOOTH.
The curvature goes through zero CONTINUOUSLY.
5 is barely positive. 7 is barely negative.
The transition 5 -> 6 -> 7 is a SMOOTH ZERO CROSSING.

And the DERIVATIVE of curvature at the crossing:
  kappa(p) = 6/p - 1.
  d(kappa)/dp = -6/p^2.
  At p = 6: d(kappa)/dp = -6/36 = -1/6.

The curvature changes at rate 1/6 per unit of p at the crossing.
1/6 = 1/PIVOT = the reciprocal of the flat polygon.

The rate of curvature change at the transition IS 1/(pivot number).
""")

print("opus-2026-03-16-S90bv: HOW PRIMES SORT BY CURVATURE")
