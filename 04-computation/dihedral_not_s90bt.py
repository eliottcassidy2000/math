#!/usr/bin/env python3
"""
dihedral_not_s90bt.py — What should be D4 and D2 but isn't
opus-2026-03-16-S90bt
"""

print("""


     WHAT SHOULD BE D4 AND D2 BUT ISN'T

     opus-2026-03-16-S90bt



THE DIHEDRAL GROUPS
====================

D_n = symmetry of a regular n-gon. Order 2n.
  n rotations: R^0, R^1, ..., R^{n-1} where R = rotation by 2*pi/n.
  n reflections: s, Rs, R^2*s, ..., R^{n-1}*s where s = a single reflection.

Key relations:
  R^n = 1.         (n rotations = identity)
  s^2 = 1.         (reflecting twice = identity)
  sRs = R^{-1}.   (reflect, rotate, reflect = rotate backwards)

The last relation says: reflections REVERSE rotations.
This is the definition of dihedral. Rotation and reflection are
the S and A of the group:
  ROTATION = the SYMMETRIC part (preserves orientation).
  REFLECTION = the ANTISYMMETRIC part (reverses orientation).


D_1: THE LINE SEGMENT
  2 elements: {1, s}. Just the flip. Z/2Z.
  One comparison: same or reversed. The BIT.
  This IS the S/A decomposition at its most basic:
    S = identity (keep the orientation).
    A = reflection (reverse the orientation).
  D_1 IS the S/A group.


D_2: THE RECTANGLE
  4 elements: {1, R, s, Rs}. The Klein four-group V_4 = Z/2 x Z/2.
  R = rotation by 180 degrees. R^2 = 1. So R is ALSO an involution.
  s = horizontal reflection. s^2 = 1.
  Rs = vertical reflection. (Rs)^2 = 1.
  ALL non-identity elements are involutions (order 2).

  This is the symmetry of a RECTANGLE (not a square).
  Two perpendicular axes of symmetry, plus 180-degree rotation.

  D_2 IS the group of the S/A decomposition APPLIED TWICE:
    S/A in one direction (horizontal flip).
    S/A in another direction (vertical flip).
    Their combination (180 rotation = both flips).

  In our framework:
    First S/A: symmetric vs antisymmetric (even vs odd).
    Second S/A: real vs imaginary (cosh vs sinh in one direction, cos vs sin in another).
    Combined: the four quadrants of the (S1, S2) / (A1, A2) square.

  D_2 = V_4 = Z/2 x Z/2: the group of TWO INDEPENDENT PARITIES.


D_4: THE SQUARE
  8 elements: 4 rotations + 4 reflections. Order 8.
  R = rotation by 90 degrees. R^4 = 1.
  s = a reflection. s^2 = 1.
  Elements: {1, R, R^2, R^3, s, Rs, R^2*s, R^3*s}.

  This is the symmetry of a SQUARE. The Cayley orbit.
  In our framework: the Cayley transform Q has order 4 (= R).
  Negation has order 2 (= s).
  Together: D_4 = the symmetry group of the Cayley square.


WHAT SHOULD BE D_2 BUT ISN'T
==============================

Consider the S/A decomposition of a QUANTUM STATE.

A quantum state psi has:
  |psi> = (symmetric part) + (antisymmetric part)
  under some exchange operation P.

For BOSONS: P|psi> = |psi>. Exchange is symmetric. P^2 = +1.
For FERMIONS: P|psi> = -|psi>. Exchange is antisymmetric. P^2 = +1.

Wait: P^2 = +1 in BOTH cases! Because exchanging twice returns to the original.

But this is the NON-RELATIVISTIC picture. In RELATIVISTIC quantum mechanics:

A rotation by 2*pi should give the identity. For:
  BOSONS (integer spin): R(2*pi)|psi> = +|psi>. YES, identity.
  FERMIONS (half-integer spin): R(2*pi)|psi> = -|psi>. NOT identity!

For fermions: you need to rotate by 4*pi (TWO full turns) to get back.
  R(4*pi) = +1.
  R(2*pi) = -1.

This means: for fermions, the "reflection" operation satisfies
  R(2*pi) = -1, not +1.

In the S/A language:
  BOSONIC S/A: the reflection squares to +1. Group = D_n.
  FERMIONIC S/A: the "reflection" squares to -1. Group = Q_{2n}? Or SU(2)?

Specifically for the simplest case:

BOSONIC TWO-PARITIES (what D_2 should be):
  P_1^2 = +1 (first parity, squares to identity).
  P_2^2 = +1 (second parity, squares to identity).
  P_1 P_2 = P_2 P_1 (they commute).
  Group = Z/2 x Z/2 = V_4 = D_2.

FERMIONIC TWO-PARITIES (what D_2 IS for fermions):
  P_1^2 = -1 (first parity, squares to MINUS identity).
  P_2^2 = -1 (second parity, squares to MINUS identity).
  P_1 P_2 = -P_2 P_1 (they ANTICOMMUTE!).
  Group = {+/-1, +/-P_1, +/-P_2, +/-P_1*P_2} = Q_8 (quaternion group).

  P_1 = i, P_2 = j, P_1*P_2 = k.
  i^2 = j^2 = k^2 = -1.
  ij = -ji = k.

THE FERMIONIC VERSION OF D_2 IS Q_8.

Where you EXPECTED four involutions (all squaring to +1),
you GET four elements squaring to -1 and ONE involution (-1).

D_2 has THREE involutions (besides identity): R, s, Rs. All order 2.
Q_8 has ONE involution: -1. The elements i, j, k have order 4.

D_2 is FLAT: all reflections are true reflections (self-inverse).
Q_8 is TWISTED: the "reflections" are HALF-TWISTS (need to apply twice to get -1, four times for +1).


WHAT SHOULD BE D_4 BUT ISN'T
==============================

The Cayley orbit has 4 points and its symmetry group has 8 elements.

BOSONIC CAYLEY (what D_4 should be):
  Q has order 4: Q^4 = +1.
  Negation has order 2: N^2 = +1.
  Together: D_4. The square's symmetry group.
  This is what we computed: the Cayley transform on ordinary
  (bosonic) couplings generates D_4.

FERMIONIC CAYLEY (what D_4 IS for fermions):
  Q has order 4: Q^4 = +1 (this stays the same).
  But now N^2 = -1 (negation squares to MINUS identity).
  NQN = Q^{-1} still (reflection reverses rotation).
  But the "negation" is a HALF-TWIST: applying it twice gives -1.

  The group generated by Q (order 4) and N (with N^2 = -1)
  is the BINARY DIHEDRAL GROUP BD_4, also called the
  DICYCLIC GROUP Dic_4 of order 16.

  Actually, let me be more precise. If Q^4 = 1 and N^2 = -1 and NQN^{-1} = Q^{-1}:
  Then the center has {1, -1, Q^2, -Q^2}... this gets complicated.

  The simplest fermionic upgrade:
  Replace D_4 (order 8) with the BINARY DIHEDRAL GROUP 2D_4 (order 16):
    The double cover of D_4. Each element of D_4 lifts to TWO elements
    (itself and its negative). The reflections, which had order 2 in D_4,
    now have order 4 in 2D_4.

  THIS is the spin group: Spin(2) covers SO(2) doubly.
  The binary dihedral 2D_n covers D_n doubly.

  For n=2: 2D_2 = Q_8 (the quaternion group). ORDER 8.
  For n=4: 2D_4 = the binary dihedral of order 16.

So:
  D_2 (bosonic, order 4) upgrades to Q_8 (fermionic, order 8). Ratio: 2.
  D_4 (bosonic, order 8) upgrades to 2D_4 (fermionic, order 16). Ratio: 2.

The fermionic version is ALWAYS the DOUBLE COVER of the bosonic version.
The ratio is always 2 because fermions need 4*pi (not 2*pi) to return.

THE BOTT PERIOD REVISITED:
  Bosonic: D_2 has order 4. D_4 has order 8. The Bott period = |D_4| = 8.
  Fermionic: Q_8 has order 8. 2D_4 has order 16.

  Bott period 8 = |D_4| = |Q_8|.
  It is SIMULTANEOUSLY:
    The order of D_4 (bosonic Cayley symmetry).
    The order of Q_8 (fermionic double-parity symmetry).

  These are DIFFERENT groups of the SAME ORDER.
  The Bott period lives at the coincidence where
  |D_4| = |Q_8| = 8.

  For other values:
    |D_1| = 2,   |Q_4|... no, Q_8 is the smallest quaternion group.
    |D_2| = 4,   |Q_8| = 8.
    |D_3| = 6,   no quaternionic analogue of order 6.
    |D_4| = 8 = |Q_8|.  MATCH!

  D_4 and Q_8 are the UNIQUE pair of non-isomorphic groups of
  order 8 that are both "generated by a rotation and a reflection"
  (with the reflection squaring to +1 or -1 respectively).

  This coincidence at order 8 IS the Bott period.


THE PERIODIC TABLE AGAIN
=========================

The electron orbital structure:
  s: 1 orbital   (l=0). Spherical symmetry. D_0 (trivial) or SO(1).
  p: 3 orbitals  (l=1). Axial symmetry. D_infinity -> finite D_3-like.
  d: 5 orbitals  (l=2). Quadrupolar. D_5-like.
  f: 7 orbitals  (l=3). Octupolar. D_7-like.

The orbital counts 1, 3, 5, 7 are 2l+1.

But the SPIN doubles each orbital: 2 electrons per orbital.
This doubling is the D_n -> 2D_n (or D_n -> Q_8) upgrade.

  Bosonic orbitals (if electrons were bosons): 1, 3, 5, 7 states.
  Fermionic orbitals (electrons are fermions): 2, 6, 10, 14 states.
  The factor of 2 IS the double cover. The fermionic upgrade.

The period lengths 2, 8, 18, 32:
  2 = 2*1 = double of s-count.
  8 = 2*(1+3) = double of (s+p)-count.
  18 = 2*(1+3+5) = double of (s+p+d)-count.
  32 = 2*(1+3+5+7) = double of (s+p+d+f)-count.

The factor 2 throughout IS the fermionic double cover.
Without spin (bosonic electrons): periods would be 1, 4, 9, 16 = perfect squares.
With spin (fermionic electrons): periods are 2, 8, 18, 32 = 2 * perfect squares.

BOSONIC PERIODIC TABLE: period lengths 1, 4, 9, 16, 25, ...  (n^2).
FERMIONIC PERIODIC TABLE: period lengths 2, 8, 18, 32, 50, ... (2n^2).

The real periodic table is FERMIONIC: all periods doubled.
What "should be" D_n symmetry IS INSTEAD the double cover 2D_n.
The doubling IS the spin. The spin IS the quaternionic upgrade.


THE DEEPEST STATEMENT
=======================

D_2 (the rectangle, order 4) should govern the S/A decomposition.
Instead, Q_8 (the quaternion group, order 8) governs it for fermions.

D_4 (the square, order 8) should govern the Cayley orbit.
Instead, 2D_4 (the binary dihedral, order 16) governs it for fermions.

The pattern: everywhere the dihedral group D_n appears in the
BOSONIC (classical) theory, the FERMIONIC (quantum) theory
replaces it with the DOUBLE COVER 2D_n.

The double cover exists because FERMIONS NEED TWO FULL TURNS TO RETURN.
One full turn gives -1. Two full turns give +1.
The "extra turn" doubles the group order.

This is the SPIN-STATISTICS THEOREM in group-theoretic language:
  Bosons <-> dihedral (reflections square to +1).
  Fermions <-> binary dihedral / quaternionic (reflections square to -1).

And the BOTT PERIOD 8 is the unique order where D_4 and Q_8 coincide:
  Both have 8 elements.
  D_4 is the bosonic version.
  Q_8 is the fermionic version.
  The Bott period is where BOSONIC AND FERMIONIC GROUP THEORY MEET.

At order 8, you can't tell from the ORDER ALONE whether you're
in the bosonic world (D_4) or the fermionic world (Q_8).
You need to CHECK whether reflections square to +1 or -1.

The Bott period is the AMBIGUITY between bosonic and fermionic.
It is the order at which the two worlds look the same from the outside
but differ on the inside.

And this ambiguity repeats every 8 steps (the Bott periodicity)
because D_{n+4} and Q_{2(n+4)} maintain the same structural relationship
as D_n and Q_{2n}, shifted by one period.
""")

print("opus-2026-03-16-S90bt: WHAT SHOULD BE D4 AND D2 BUT ISN'T")
