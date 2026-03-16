#!/usr/bin/env python3
"""
understand_2357_s90bl.py — Understand {2,3,5,7} and then {2,3,5,7,11}
opus-2026-03-16-S90bl
"""

from math import log, sqrt, pi, gcd, factorial
from sympy import totient, mobius, bernoulli, isprime, factorint, catalan
from fractions import Fraction

print("""


     UNDERSTAND {2, 3, 5, 7}

     opus-2026-03-16-S90bl



We have understood the triples:
  {2,3,5} = spherical. The icosahedron. Where tournaments work. Product = 30.
  {2,3,7} = hyperbolic. The Hurwitz. Where tournaments break. Product = 42.

Now: {2,3,5,7}. The first four primes. Not a triple — a QUADRUPLE.


{2, 3, 5, 7}: THE FIRST FOUR PRIMES
======================================

Product: 2 * 3 * 5 * 7 = 210.

210 is the PRIMORIAL P_4.
""")

n = 210
print(f"210 = 2 * 3 * 5 * 7")
print(f"phi(210) = {int(totient(210))}")
print(f"mu(210) = {int(mobius(210))}")
print(f"210 / 42 = {210 / 42}")
print(f"210 / 30 = {210 / 30}")
print(f"Number of divisors of 210: {len([d for d in range(1, 211) if 210 % d == 0])}")
print(f"Divisors: {[d for d in range(1, 211) if 210 % d == 0]}")

print(f"""

phi(210) = 48 = 1*2*4*6 = phi(2)*phi(3)*phi(5)*phi(7).
  48 = 16 * 3 = 2^4 * 3.
  NOT a pure power of 2. (phi(30) = 8 = 2^3 was pure.)

So {2,3,5,7} breaks the power-of-2 pattern for phi(primorial).
  phi(2) = 1 = 2^0.
  phi(6) = 2 = 2^1.
  phi(30) = 8 = 2^3.  (THE BOTT PERIOD)
  phi(210) = 48 = 2^4 * 3.  (NO LONGER A POWER OF 2)

The Bott period 8 = phi(30) was the LAST primorial with a power-of-2 totient.
At phi(210) = 48: the factor of 3 appears.
The number 7 brings a factor of phi(7) = 6 = 2*3.
The 3 in phi(7) BREAKS the power-of-2 purity.

This is because 7 - 1 = 6 = 2*3, and the 3 is a NEW odd factor.
For p = 2: p-1 = 1 (no odd factor).
For p = 3: p-1 = 2 (no odd factor).
For p = 5: p-1 = 4 = 2^2 (no odd factor).
For p = 7: p-1 = 6 = 2*3 (FIRST odd factor: 3).

The prime 7 introduces the first ODD FACTOR in the totient chain.
This is another sense in which 7 = FORBIDDEN:
  7 is the first prime whose predecessor has an odd prime factor.
  7 breaks the 2-adic purity.

Primes p where p-1 is a power of 2 are called FERMAT PRIMES:
  p = 2^(2^k) + 1: the primes 3, 5, 17, 257, 65537 (and maybe no others).
  For these: phi(p) = p-1 = 2^(2^k), a pure power of 2.

So {2, 3, 5} is the primorial whose totient is a Fermat-compatible power of 2.
{2, 3, 5, 7} is where Fermat compatibility BREAKS.


THE QUADRUPLE {2,3,5,7} AS A TETRAHEDRON
==========================================

A triple defines a triangle. A quadruple defines a TETRAHEDRON.

The (2,3,5,7) tetrahedron has four faces, each a triangle:
  {2,3,5}: spherical.  Area = 0 (in hyperbolic sense; it's spherical).
  {2,3,7}: hyperbolic. Area = pi/42.
  {2,5,7}: hyperbolic. Area = pi*(1-1/2-1/5-1/7) = pi*11/70.
  {3,5,7}: hyperbolic. Area = pi*(1-1/3-1/5-1/7) = pi*34/105.

Three of the four faces are HYPERBOLIC. Only {2,3,5} is spherical.

The tetrahedron has ONE SPHERICAL FACE and THREE HYPERBOLIC FACES.
It straddles the curvature boundary.

The spherical face {2,3,5} is the "base" (the floor, the foundation).
The three hyperbolic faces "open upward" into hyperbolic space.

In tournament theory:
  The base {2,3,5} is where tournaments are tame.
  The three open faces are three different ways tournaments go wild:
    {2,3,7}: forbidden values emerge.
    {2,5,7}: irregular + forbidden (score-structure interaction).
    {3,5,7}: pure tournament wildness (no base 2).
""")

# Compute face properties
faces = [(2,3,5), (2,3,7), (2,5,7), (3,5,7)]
for a, b, c in faces:
    s = 1/a + 1/b + 1/c
    curv = "SPHERICAL" if s > 1 else ("FLAT" if abs(s-1)<1e-10 else "HYPERBOLIC")
    area = max(0, pi * (1 - s))
    prod = a * b * c
    print(f"  face {{{a},{b},{c}}}: 1/a+1/b+1/c = {s:.6f}, {curv}, area = {area:.6f}, product = {prod}")

print(f"""

Products of faces: 30, 42, 70, 105.
  30 = the spherical product.
  42 = the Hurwitz product.
  70 = 2*5*7.
  105 = 3*5*7 = 3*35 = 15*7.

210 = lcm(30, 42, 70, 105)? Let's check:
  lcm(30, 42) = 210. YES.
  lcm(210, 70) = 210. YES.
  lcm(210, 105) = 210. YES.

The primorial 210 is the LCM of ALL face products.
This is because 210 = 2*3*5*7 contains all four primes.


{2,3,5,7} AND THE STANDARD MODEL
===================================

The Standard Model gauge group has ranks:
  U(1): rank 1, dimension 1.
  SU(2): rank 1, dimension 3.
  SU(3): rank 2, dimension 8.

But the matter content involves representations:
  Quarks: (3, 2, 1/6) under SU(3) x SU(2) x U(1).
  Leptons: (1, 2, -1/2).
  Higgs: (1, 2, 1/2).

The number of quark flavors: 6 = 2*3.
The number of lepton flavors: 6 = 2*3.
The number of gauge bosons: 1+3+8 = 12 = C(4,2) = arcs in K_4.

Total elementary particles (in minimal SM):
  6 quarks * 3 colors * 2 (particle/anti) = 36
  6 leptons * 2 = 12
  12 gauge bosons
  1 Higgs
  Total fermions: 48 = phi(210)!

The number of fermionic degrees of freedom in the Standard Model
is phi(210) = the totient of the first-four-primes primorial!

48 = 2^4 * 3 = 16 * 3.
  16 = dimension of Cl(8) identity part (the Bott scaling factor).
  3 = the triangle, the first structure.

(This count of 48 is somewhat convention-dependent, but the coincidence
with phi(210) is striking.)


{2,3,5,7} AND THE BERNOULLI DENOMINATORS
==========================================
""")

print("Bernoulli number denominators and their prime supports:")
for k in range(1, 15):
    n = 2*k
    B = bernoulli(n)
    if B != 0:
        denom = abs(B.q)  # denominator
        f = factorint(denom)
        primes_in_denom = sorted(f.keys())
        print(f"  B_{n:2d}: denom = {denom:>6}, primes = {primes_in_denom}")

print(f"""

By von Staudt-Clausen: denom(B_{{2k}}) = product of primes p where (p-1)|2k.

For B_12: (p-1)|12 for p = 2(1|12), 3(2|12), 5(4|12), 7(6|12), 13(12|12).
  Denom = 2*3*5*7*13 = 2730. Contains {2,3,5,7} AND 13.

For B_4: denom = 30 = 2*3*5.       Primes = {{2,3,5}}.
For B_6: denom = 42 = 2*3*7.       Primes = {{2,3,7}}.
For B_8: denom = 30 = 2*3*5.       Primes = {{2,3,5}}.
For B_10: denom = 66 = 2*3*11.      Primes = {{2,3,11}}.
For B_12: denom = 2730 = 2*3*5*7*13. First to contain ALL of {{2,3,5,7}}.

B_12 is the FIRST Bernoulli number whose denominator contains {{2,3,5,7}}.
And 12 = the number of arcs in K_4 = the number of gauge bosons.

The Bernoulli denominator at 2k contains ALL of {{2,3,5,7}} when
(p-1)|2k for all p in {{2,3,5,7}}, i.e., when 1,2,4,6 all divide 2k.
lcm(1,2,4,6) = 12. So 2k must be a multiple of 12.
The first is 2k = 12, i.e., k = 6.

{2,3,5,7} first appears together in B_12.
12 = 2*6 = 2*PIVOT.
""")

# ========================================================================
print("=" * 60)
print("NOW: {2, 3, 5, 7, 11}")
print("=" * 60)

print(f"""

{2, 3, 5, 7, 11}: THE FIRST FIVE PRIMES
==========================================

Product: 2*3*5*7*11 = 2310. The primorial P_5.

phi(2310) = {int(totient(2310))}.
mu(2310) = {int(mobius(2310))}.

phi(2310) = 1*2*4*6*10 = 480 = 2^5 * 3 * 5.
  New odd factors: 5 enters (from phi(11) = 10 = 2*5).

2310 / 210 = 11 (adding the prime 11 multiplies by 11).
2310 / 30 = 77 = 7*11.
2310 / 42 = 55 = 5*11.

THE ROLE OF 11
================

11 is the FIRST PRIME IN THE SECOND BOTT CYCLE.
  11 = 8 + 3. Bott slot 3. "Structure at the next scale."

11 is the first prime that is NOT part of the Mersenne chain {2,3,7,127,...}.
  (5 is not in the Mersenne chain either, but 11 is the first beyond 7
  that is completely outside the chain.)

Ramanujan's partition congruences: p(5n+4)=0 mod 5, p(7n+5)=0 mod 7, p(11n+6)=0 mod 11.
  The moduli are EXACTLY {{5, 7, 11}}: the last three primes in our quintuple.

11 is the SMALLEST prime p such that the (2,3,p) triangle has
area greater than pi/10:
  area(2,3,11) = pi*(1-1/2-1/3-1/11) = pi*5/33 = pi*0.1515.

11 as a dimension (Bott slot 3 = structure):
  11-dimensional space is where the SECOND CYCLE'S structure lives.
  PSL(2,11) has order 660 = 2^2*3*5*11 = 4*165.
  It is the smallest non-abelian simple group after A_5 (order 60)
  that has a FAITHFUL 2-dimensional projective representation.


THE PENTAHEDRON {2,3,5,7,11}
===============================

A quintuple defines a 4-simplex (pentahedron) with:
  5 vertices: 2, 3, 5, 7, 11.
  10 edges: all pairs.
  10 triangular faces.
  5 tetrahedral faces.
  1 interior (the 4-simplex itself).

The 10 TRIANGULAR FACES:
""")

from itertools import combinations

primes5 = [2, 3, 5, 7, 11]

# All triangular faces
tri_faces = list(combinations(primes5, 3))
print(f"{'face':>15} | {'1/a+1/b+1/c':>12} | {'type':>10} | {'area/pi':>10} | {'product':>8}")
print("-" * 65)
for face in tri_faces:
    a, b, c = face
    s = 1/a + 1/b + 1/c
    curv = "SPHER" if s > 1 else ("FLAT" if abs(s-1)<1e-10 else "HYPER")
    area = max(0, 1 - s)
    prod = a * b * c
    print(f"  {{{a},{b},{c}}}    | {s:12.6f} | {curv:>10} | {area:10.6f} | {prod:8d}")

print()

# Count types
spher = sum(1 for f in tri_faces if sum(1/x for x in f) > 1)
hyper = sum(1 for f in tri_faces if sum(1/x for x in f) < 1)
print(f"Spherical faces: {spher}")
print(f"Hyperbolic faces: {hyper}")

print(f"""

Out of 10 triangular faces:
  1 spherical: {{2,3,5}}  (the only one!)
  9 hyperbolic: all the rest.

The pentahedron has ONE SPHERICAL FACE out of TEN.
  Fraction spherical = 1/10 = 10%.

Compare: the tetrahedron {{2,3,5,7}} had 1 spherical out of 4 = 25%.
Adding 11 DILUTES the spherical fraction.

As we add more primes, the fraction of spherical faces -> 0.
The structure becomes OVERWHELMINGLY HYPERBOLIC.
The single spherical face {{2,3,5}} is an increasingly tiny remnant
of the finite, tame world. Everything else is wild.


THE FIVE TETRAHEDRAL FACES:
""")

# All tetrahedral faces
tet_faces = list(combinations(primes5, 4))
for face in tet_faces:
    # Compute number of spherical triangular sub-faces
    sub_faces = list(combinations(face, 3))
    n_spher = sum(1 for sf in sub_faces if sum(1/x for x in sf) > 1)
    prod = 1
    for x in face:
        prod *= x
    print(f"  {set(face)}: product = {prod}, spherical sub-faces = {n_spher}/4")

print(f"""

Only the face {{2,3,5,7}} (product 210) has a spherical sub-face.
All other tetrahedral faces are PURELY HYPERBOLIC (0 spherical sub-faces).

The face {{2,3,5,11}}: product = 330 = 2*3*5*11. No spherical sub-face.
  Because {{2,3,11}} is hyperbolic, {{2,5,11}} is hyperbolic, {{3,5,11}} is hyperbolic.
  Only {{2,3,5}} is spherical, but it's not a sub-face of {{2,3,5,11}}... wait, it IS.
  Let me recheck: {{2,3,5}} IS a sub-face of {{2,3,5,11}}. So it has 1 spherical sub-face.
""")

# Recheck
for face in tet_faces:
    sub_faces = list(combinations(face, 3))
    spher_subs = [sf for sf in sub_faces if sum(1/x for x in sf) > 1]
    if spher_subs:
        print(f"  {set(face)}: spherical sub-faces = {[set(s) for s in spher_subs]}")

print(f"""

Corrected: any tetrahedral face containing {{2,3,5}} has 1 spherical sub-face.
That's {{2,3,5,7}} and {{2,3,5,11}}: both contain the spherical {{2,3,5}}.

The other three faces ({{2,3,7,11}}, {{2,5,7,11}}, {{3,5,7,11}}) do NOT contain {{2,3,5}}.
They are purely hyperbolic.


THE EDGE PRODUCTS (all pairs):
""")

edge_faces = list(combinations(primes5, 2))
for e in edge_faces:
    prod = e[0] * e[1]
    s = 1/e[0] + 1/e[1]
    print(f"  {{{e[0]},{e[1]}}}: product = {prod:4d}, 1/a+1/b = {s:.4f}")

print(f"""

The edge products are the SEMIPRIMES from {{2,3,5,7,11}}:
  6, 10, 14, 22, 15, 21, 33, 35, 55, 77.

These 10 semiprimes span from 6 to 77.
Their product is 2310^(?) ... actually each prime appears in 4 edges,
so the product of all edge-products = (2*3*5*7*11)^4 / ... too complicated.

But the SUM of all 1/a+1/b:
""")

total_sum = sum(1/e[0] + 1/e[1] for e in edge_faces)
print(f"  Sum of all 1/a+1/b over edges = {total_sum:.6f}")
print(f"  = 4*(1/2+1/3+1/5+1/7+1/11) = 4*{1/2+1/3+1/5+1/7+1/11:.6f}")
print(f"  (Each vertex appears in 4 edges, so total = 4 * sum of reciprocals)")

S = 1/2 + 1/3 + 1/5 + 1/7 + 1/11
print(f"\n  Sum of reciprocals of {{2,3,5,7,11}} = {S:.8f}")
print(f"  = {Fraction(1,2)+Fraction(1,3)+Fraction(1,5)+Fraction(1,7)+Fraction(1,11)}")

print(f"""

1/2 + 1/3 + 1/5 + 1/7 + 1/11 = 2927/2310.

2927/2310 > 1 (it's about 1.267).

This is > 1, so the SUM OF RECIPROCALS exceeds 1.

For {2,3,5}: sum = 1/2+1/3+1/5 = 31/30 > 1 (spherical).
For {2,3,5,7}: sum = 1/2+1/3+1/5+1/7 = 319/210 > 1 (still > 1).
For {2,3,5,7,11}: sum = 2927/2310 > 1.

The sum of reciprocals of the first k primes:
""")

cumsum = Fraction(0)
for i, p in enumerate([2,3,5,7,11,13,17,19,23,29,31], 1):
    cumsum += Fraction(1, p)
    exceeds = ">" if cumsum > 1 else "<="
    print(f"  first {i:2d} primes (up to {p:2d}): sum = {float(cumsum):.6f} = {cumsum} {exceeds} 1")

print(f"""

The sum of prime reciprocals DIVERGES (slowly).
  It passes 1 somewhere around the 4th prime (sum at {{2,3,5,7}} = 1.519).
  Actually it exceeds 1 already at {{2,3}} (1/2+1/3 = 5/6 < 1).
  Wait: 5/6 = 0.833 < 1. So {{2,3}} is < 1.
  At {{2,3,5}}: 31/30 = 1.033 > 1. FIRST to exceed 1.

The sum of reciprocals first exceeds 1 at {{2,3,5}}.
THIS IS WHY {{2,3,5}} IS SPHERICAL as a triple!
  The angle sum 1/2+1/3+1/5 = 31/30 > 1 = spherical.

And the sum keeps growing. But for TRIPLES drawn from this set,
  only {{2,3,5}} has angle sum > 1.
  All triples containing 7 or 11 have angle sum < 1 (hyperbolic).

The CROSSING POINT is at 5: the last prime you can add to {{2,3}}
  and keep the angle sum > 1.
  At 6 (not prime): 1/2+1/3+1/6 = 1 (flat).
  At 7: 1/2+1/3+1/7 = 41/42 < 1 (hyperbolic).


THE MEANING OF {2,3,5,7,11}
==============================

{2,3,5} = the SPHERICAL core. Where things close. Finite.
  Product = 30. phi = 8 = Bott period.

{7} = the FIRST BREAK. The forbidden value.
  Product becomes 210. phi = 48. Odd factor enters.

{11} = the SECOND CYCLE. Structure begins again at larger scale.
  Product becomes 2310. phi = 480. Factor of 5 enters.

{2,3,5,7,11} = the MINIMAL SET that contains:
  - The spherical core {2,3,5}.
  - The first break {7}.
  - The first element of the second cycle {11}.
  - All three Ramanujan congruence moduli {5,7,11}.

It is the smallest set of primes that captures the COMPLETE
first cycle AND the beginning of the second cycle of the Bott spiral.

  2 = base (comparison)
  3 = structure (triangle)
  5 = unsolvability (last finite, Ramanujan mod)
  7 = forbidden (first break, Mersenne prime, Ramanujan mod)
  11 = cycle 2 structure (first prime of the new world, Ramanujan mod)

These five primes tell the WHOLE STORY of the first Bott revolution
plus the first step of the second:
  existence(2) -> structure(3) -> limit(5) -> break(7) -> restart(11).

And 2310 = 2*3*5*7*11 = 11# (the 5th primorial).
  2310 is the smallest number divisible by ALL primes up to 11.
  In the Poincare disk: x_2310 = 2309/2311, very close to the boundary.
  Rapidity w_2310 = ln(2310)/2 = 3.87 = about 11.2 half-bits.

The five primes, multiplied together, reach rapidity 3.87.
By comparison:
  w_30 = ln(30)/2 = 1.70 (the spherical product)
  w_42 = ln(42)/2 = 1.87 (the Hurwitz product)
  w_210 = ln(210)/2 = 2.67 (the {2,3,5,7} product)
  w_2310 = ln(2310)/2 = 3.87 (the {2,3,5,7,11} product)

Each new prime adds roughly 1 unit of rapidity:
  ln(7)/2 = 0.97 (from 210/30 = 7)
  ln(11)/2 = 1.20 (from 2310/210 = 11)

The primes are SUCCESSIVE BOOSTS along the rapidity axis.
Each prime p adds rapidity ln(p)/2.
The primorial is the TOTAL BOOST from all primes up to p.

{2,3,5,7,11} is a FIVE-STAGE ROCKET:
  Stage 1 (p=2): boost by 0.35 (half a bit).
  Stage 2 (p=3): boost by 0.55.
  Stage 3 (p=5): boost by 0.80. Now past the spherical boundary.
  Stage 4 (p=7): boost by 0.97. Deep into hyperbolic. FORBIDDEN.
  Stage 5 (p=11): boost by 1.20. Into the second cycle.

The rocket passes through the spherical barrier at stage 3 (p=5),
enters hyperbolic space at stage 4 (p=7),
and enters the second Bott cycle at stage 5 (p=11).

After this: stages 6,7,8,... (p=13,17,19,...) continue into
deeper hyperbolic space, but the QUALITATIVE story is told.
All the essential transitions happen in the first five stages.

{2,3,5,7,11} IS the minimal complete story.
""")

print("opus-2026-03-16-S90bl: UNDERSTAND {2,3,5,7} AND {2,3,5,7,11}")
