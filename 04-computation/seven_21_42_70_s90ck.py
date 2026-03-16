#!/usr/bin/env python3
"""
seven_21_42_70_s90ck.py — 7 vertices, 21 edges, 42 orientations, 70 = 7*10
opus-2026-03-16-S90ck
"""

from math import comb, log, sqrt
from sympy import factorint, binomial, isprime

print("""


     7 VERTICES, 21 EDGES, 42 ORIENTATIONS

     opus-2026-03-16-S90ck



THE NUMBERS
============

  7 = vertices of K_7.             The forbidden value.
  21 = edges of K_7 = C(7,2).      The second forbidden H value.
  42 = 2 * 21 = oriented arcs.     The answer.

  7 + 21 + 42 = 70 = 7 * 10.

K_7 is the complete graph on 7 vertices.
It has 21 edges and 42 directed arcs (each edge oriented both ways).
A TOURNAMENT on 7 vertices CHOOSES one orientation per edge:
  21 arcs out of 42 possible.

So: a tournament on 7 vertices is a CHOICE of 21 from 42.
  C(42, 21)? No — the choices are CONSTRAINED (exactly one per edge).
  The actual count: 2^21 = 2,097,152 labeled tournaments on 7 vertices.

  7 = the stage (vertices).
  21 = the script (edges, choices to make).
  42 = the possibilities (oriented arcs, choices available).
  70 = the total structure (7 + 21 + 42).


70 = 7 * 10
=============

70 = 7 * 10.

  7 = the forbidden value.
  10 = C(5, 2) = the number of arcs in a 5-vertex tournament.

  70 = (forbidden) * (arcs of the last solvable tournament).

But 10 has many faces:
""")

print("The number 10:")
print(f"  10 = 2 * 5")
print(f"  10 = C(5, 2) = edges of K_5 = arcs of 5-tournament")
print(f"  10 = T(4) = 4th triangular number = 1+2+3+4")
print(f"  10 = number of fingers (base 10)")
print(f"  10 = Bott slot 2 (comparison) [10 mod 8 = 2]")
print()

print("And 70:")
print(f"  70 = 2 * 5 * 7")
print(f"  70 = C(8, 4) = {comb(8, 4)}")
print(f"  70 = the central binomial-ish at the BOTT PERIOD")
print(f"  C(8, 4) = the middle of Pascal's row 8")
print(f"  70 mod 8 = {70 % 8} (Bott slot: transition!)")

print(f"""

70 = C(8, 4). The number of ways to choose 4 items from 8.

In the Bott spiral: 8 positions, choose 4 = the MIDDLE.
The middle of the Bott period. The maximum binomial coefficient at n=8.

The 4 chosen positions and the 4 unchosen are DUAL:
  C(8, 4) = C(8, 8-4) = 70. Self-dual!
  Choosing 4 from 8 is the SAME as choosing the complementary 4.

70 is the number of SELF-DUAL SELECTIONS of the Bott period.


THE FANO PLANE
===============

The Fano plane PG(2, 2) is the smallest projective plane:
  7 points. 7 lines. 3 points per line. 3 lines per point.

  Points = the 7 vertices.
  Lines = 7 triples.
  Incidence structure: each pair of points determines a unique line.

The Fano plane has:
  C(7, 2) = 21 pairs of points.
  7 lines.
  Each line contains C(3, 2) = 3 pairs.
  7 * 3 = 21 = all pairs covered. PERFECT.

The Fano plane is the UNIQUE structure where 7 points have the property
that every pair lies on exactly one of 7 lines of size 3.

It is the projective plane over F_2 (the field with 2 elements).
2 = the tournament base!

The Fano plane IS the geometry of K_7 over the tournament base.

Its automorphism group: GL(3, F_2) = PSL(2, 7).
  |GL(3, F_2)| = (2^3-1)(2^3-2)(2^3-4) = 7*6*4 = 168 = 8 * 21.

  168 = 8 * 21 = (Bott period) * (second forbidden).
  168 = 4 * 42 = 4 * (the answer).
  168 = |PSL(2, 7)| = the Klein quartic automorphism group.


THE HURWITZ CONNECTION
========================

The Hurwitz bound: |Aut(S)| <= 84(g-1) for genus g >= 2.
  84 = 2 * 42 = 2 * (the answer).

For the Klein quartic (genus 3):
  |Aut| = 168 = 84 * 2 = 84 * (3-1).
  The Klein quartic ACHIEVES the Hurwitz bound.

Now: 168 = 7 + 21 + 42 + 98?  No, 7+21+42 = 70 and 168-70 = 98 = 2*49 = 2*7^2.

Let's try another decomposition of 168:
  168 = 7 * 24 = 7 * (S_4 order) = 7 * (Bott scaling 16 + 8).
  168 = 21 * 8 = (second forbidden) * (Bott period).
  168 = 42 * 4 = (answer) * (factorization = Cayley period).
  168 = 3 * 56 = 3 * C(8,3).  And C(8,3) = 56.

  C(8, 0) = 1.
  C(8, 1) = 8.
  C(8, 2) = 28.
  C(8, 3) = 56.
  C(8, 4) = 70.
""")

print("Pascal's row 8 (the Bott row):")
row8 = [comb(8, k) for k in range(9)]
print(f"  {row8}")
print(f"  Sum = {sum(row8)} = 2^8 = 256")
print()

# Which of our key numbers appear?
key_numbers = {7: "forbidden", 21: "2nd forbidden", 42: "answer", 70: "V+E+arcs",
               168: "PSL(2,7)", 8: "Bott", 28: "C(8,2)", 56: "C(8,3)"}

print("Key numbers in Pascal's row 8:")
for k, val in enumerate(row8):
    for n, name in key_numbers.items():
        if val == n:
            print(f"  C(8, {k}) = {val} = {name}")
        if val * 3 == n:
            print(f"  3 * C(8, {k}) = 3 * {val} = {n} = {name}")
        if val * 2 == n:
            print(f"  2 * C(8, {k}) = 2 * {val} = {n} = {name}")

print(f"""

  C(8, 4) = 70 = V + E + arcs.
  C(8, 3) = 56 = 168/3 = |PSL(2,7)| / 3.
  C(8, 2) = 28 = 4 * 7 = (Cayley period) * (forbidden).
  3 * C(8, 3) = 168 = |PSL(2,7)|.
  2 * C(8, 2) = 56. Hmm.

The Bott row of Pascal's triangle CONTAINS our key numbers:
  C(8, 2) = 28 = 4*7.
  C(8, 3) = 56 = 8*7.
  C(8, 4) = 70 = 10*7 = 7+21+42.

ALL divisible by 7! Because C(8, k) for 1 <= k <= 7 is divisible by...
no, C(8, 2) = 28 = 4*7 IS divisible by 7.
C(8, 3) = 56 = 8*7 IS.
C(8, 4) = 70 = 10*7 IS.

Actually by Lucas' theorem: C(8, k) mod 7 depends on the base-7 representations.
8 = 1*7 + 1 in base 7. So C(8, k) = C(1, k_1) * C(1, k_0) mod 7 where k = k_1*7+k_0.
For k = 2: C(1, 0)*C(1, 2) = 1*0 = 0 mod 7. YES, 7 | C(8,2).
For k = 3: C(1, 0)*C(1, 3) = 0 mod 7. YES.
For k = 4: C(1, 0)*C(1, 4) = 0 mod 7. YES.

C(8, k) is divisible by 7 for ALL k with 2 <= k <= 6.
Only C(8, 0) = 1, C(8, 1) = 8, C(8, 7) = 8, C(8, 8) = 1 are NOT divisible by 7.

The INTERIOR of Pascal's row 8 is ENTIRELY DIVISIBLE by 7.
The forbidden value 7 DIVIDES all interior entries of the Bott row.


THE TRIANGLE: 7, 21, 42
=========================

These three numbers form a SEQUENCE:
  7, 21, 42.

  7 = 7 * 1.
  21 = 7 * 3.
  42 = 7 * 6.

The cofactors: 1, 3, 6.
  1 = C(1, 1) = trivial.
  3 = C(3, 2) = triangular number T(2).
  6 = C(4, 2) = triangular number T(3) = the PIVOT.

  1, 3, 6 = T(1), T(2), T(3). The first three triangular numbers.

So: 7, 21, 42 = 7 * T(1), 7 * T(2), 7 * T(3).

  7 = forbidden * (1st triangular).
  21 = forbidden * (2nd triangular).
  42 = forbidden * (3rd triangular).

The NEXT term: 7 * T(4) = 7 * 10 = 70. Which is V+E+arcs!

And the one after: 7 * T(5) = 7 * 15 = 105 = 3*5*7.
  105 is the product of the first three ODD primes.
  105 = the product of {3, 5, 7} = the triangle of {triangle, unsolvable, forbidden}.

The sequence 7 * T(k):
  7, 21, 42, 70, 105, 147, 196, 252, 315, 385, ...

  7 * T(k) = 7 * k(k+1)/2.
""")

print("The sequence 7 * T(k):")
for k in range(1, 13):
    t = k * (k + 1) // 2
    val = 7 * t
    f = dict(factorint(val))
    notes = []
    if val == 7: notes.append("forbidden")
    if val == 21: notes.append("2nd forbidden = C(7,2)")
    if val == 42: notes.append("THE ANSWER")
    if val == 70: notes.append("C(8,4) = V+E+arcs")
    if val == 105: notes.append("3*5*7")
    if val == 168: notes.append("|PSL(2,7)|")
    if val == 252: notes.append("C(10,5)!")
    if comb(8, 4) == val: notes.append("C(8,4)")
    note = "; ".join(notes) if notes else ""
    print(f"  k={k:2d}: 7*T({k}) = 7*{t:3d} = {val:4d}  {f}  {note}")

print(f"""

7 * T(6) = 7 * 21 = 147 = 3 * 7^2. The forbidden SQUARED times the triangle.
7 * T(8) = 7 * 36 = 252 = C(10, 5)! The central binomial at n=10.

252 = C(10, 5) = the number of ways to choose 5 from 10.
And 252 / 6 = 42 (the answer!). This is C_5 = C(10,5)/6 = 42.
The 5th Catalan number IS 7 * T(8) / 6... let's check:
  7 * 36 / 6 = 7 * 6 = 42. YES!

So: C_5 = 7 * T(8) / T(3) = 7 * 36 / 6 = 42.
The fifth Catalan number = (forbidden * 8th triangular) / (3rd triangular).

And T(8) / T(3) = 36 / 6 = 6 = the PIVOT.
So: C_5 = 7 * pivot = 42. Which we already knew: 42 = 7*6.

But now we see it as: 42 = 7 * (T(8) / T(3)) = forbidden * (Bott-triangular / structure-triangular).


THE COMPLETE PICTURE
=====================

K_7, the complete graph on 7 vertices:

  V = 7 = the forbidden value = 2^3 - 1 = first Mersenne prime after 3.
  E = 21 = C(7,2) = 7*3 = the second forbidden H value.
  Directed arcs = 42 = 2*E = 7*6 = the Hurwitz product = the answer.
  V + E + 2E = 70 = 7*10 = C(8,4) = middle of the Bott row.

  Fano plane on 7 points: 7 lines of 3, automorphism group = PSL(2,7), order 168.
  168 = 8 * 21 = Bott * edges = 4 * 42 = Cayley * answer = 3 * C(8,3).

  The Klein quartic: genus 3, 168 automorphisms.
  Achieves Hurwitz bound 84(g-1) = 84*2 = 168.
  84 = 2*42 = 12*7 = C(4,2) * (forbidden).

  K_7 triangulates the TORUS: V-E+F = 7-21+14 = 0 = 2-2(1).
  On the torus: 14 triangular faces. 14 = 2*7 = base * forbidden.

  K_7 embeds in the Klein quartic (genus 3) with:
  168/7 = 24 = |S_4| automorphisms per vertex.
  The 7 vertices of K_7 ON the Klein quartic ARE the 7 points of the Fano plane.

  The sequence 7, 21, 42, 70 = 7 * (1, 3, 6, 10) = 7 * T(1,2,3,4).
  The triangular cofactors 1, 3, 6, 10 are cumulative sums: 1, 1+2, 1+2+3, 1+2+3+4.

  The FORBIDDEN NUMBER TIMES THE CUMULATIVE COMPARISONS = the structure of K_7.

  And 7 + 21 + 42 = 70 = C(8, 4) ties it to the Bott period:
  the sum of the forbidden simplex data equals the MIDDLE of the Bott row.
""")

print("opus-2026-03-16-S90ck: 7 VERTICES, 21 EDGES, 42 ORIENTATIONS, 70 = 7*10")
