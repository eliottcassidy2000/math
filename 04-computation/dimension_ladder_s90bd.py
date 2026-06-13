#!/usr/bin/env python3
"""
dimension_ladder_s90bd.py — The dimension ladder: from d to d+1
opus-2026-03-16-S90bd

The self-reference hierarchy:
  d=0: a point. The number 1. No comparison.
  d=1: the real line. Addition. The additive group (R,+).
  d=2: the complex plane. Multiplication. The multiplicative group (C*,*).
  d=?: the formal group ((-1,1), F). Where does this sit?

The question: assign a dimension to F, and find the map d -> d+1.
"""

from math import log, sqrt, pi, e, tanh, atanh, cosh, sinh
from cmath import tanh as ctanh, atanh as catanh, exp as cexp, log as clog
import numpy as np

print("""


     THE DIMENSION LADDER: FROM d TO d+1

     opus-2026-03-16-S90bd



THE HIERARCHY AS IT STANDS:

  d=0: A POINT. The number 1. No structure. No comparison.
       The group: {1} (trivial).
       The space: a single point.

  d=1: THE REAL LINE. The additive group (R, +).
       Structure: ordering, distance, addition.
       The rapidity axis. The strip.
       ONE unbounded direction.

  d=2: THE COMPLEX PLANE. The multiplicative group (C*, *).
       Structure: scaling + rotation. Magnitude + phase.
       TWO directions: real (magnitude) and imaginary (phase).
       Or equivalently: the Poincare disk D with F(w1,w2)=(w1+w2)/(1+w1w2).

So the formal group lives at d=2: it IS the Poincare disk,
which IS the complex plane (via arctanh), which has 2 real dimensions.

But the formal group also has a 1-dimensional RESTRICTION:
the real diameter (-1,1) with F restricted to reals.
This is a 1-dimensional formal group = d=1.

So the formal group F operates at d=1 (real restriction)
and d=2 (full complex extension) simultaneously.

THE QUESTION: what is d=3? And how do we get there?


THE LADDER PRINCIPLE
=====================

Going from d to d+1 has always meant the same thing:

  d=0 -> d=1: COMPARE two copies of d=0.
    Two points -> a line between them.
    The act of comparison creates a 1D structure.
    The bit: "which of the two?"

  d=1 -> d=2: COMPARE two copies of d=1.
    Two lines -> a plane between them.
    Real axis vs imaginary axis.
    The complex number: "how much of each direction?"

  d=2 -> d=3: COMPARE two copies of d=2.
    Two planes -> ???

What IS the comparison of two complex planes?

ANSWER: The QUATERNIONS.

  H = C + C*j.
  A quaternion q = a + bi + cj + dk has 4 real components.
  But as a "comparison of two complex planes," it's 2 complex dimensions.
  And as a geometric object, it's a 3-sphere S^3 (unit quaternions).

Wait. Let me think about this more carefully.
""")

# ========================================================================

print("=" * 60)
print("THE d -> d+1 MAP")
print("=" * 60)

print("""
At each level, the structure is:

  d=0: {1} with trivial operation.
       The "group" is just the identity.
       Formal group: F_0(x,y) = 0 (trivial).

  d=1: (-1,1) with F(x,y) = (x+y)/(1+xy).
       The coupling interval. One real dimension.
       Formal group: F_1 = the tanh addition law.
       Linearized by arctanh: coupling -> rapidity.
       The NUMBER 2 lives here at x = 1/3.

  d=2: D (unit disk) with F(w1,w2) = (w1+w2)/(1+w1*w2).
       The SAME formal group, but on complex inputs.
       Two real dimensions: boost + rotation.
       Linearized by arctanh: complex coupling -> complex rapidity.
       The Poincare disk.

Now: the step from d=1 to d=2 was COMPLEXIFICATION.
We took F on the real line and EXTENDED it to the complex plane.
This doubled the dimension (1 real -> 2 real = 1 complex).

The step from d=0 to d=1 was CREATION OF COMPARISON.
We went from "nothing" to "a binary choice."
This created 1 dimension from 0.

What is the PATTERN?

  d=0 -> d=1: CREATION (0 dimensions of comparison -> 1)
  d=1 -> d=2: COMPLEXIFICATION (1 real -> 1 complex = 2 real)

What should d=2 -> d=3 be?

Option A: ANOTHER COMPLEXIFICATION.
  Extend from C to the quaternions H.
  The formal group on H: F(q1,q2) = (q1+q2)/(1+q1*conj(q2))???
  This doesn't work directly because quaternions are noncommutative.

Option B: SELF-REFERENCE.
  The d=2 structure (the Poincare disk) looks at ITSELF.
  The Mersenne chain: 2 -> 3 -> 7 -> 127 -> ...
  Each step takes the output and feeds it back as input.

Option C: FIBRATION.
  Stack copies of the d=2 disk parametrized by a new variable.
  A FIBER BUNDLE: base = d=1 line, fiber = d=1 disk.
  Total space = d=2 in a different sense.

Let me explore Option B more carefully.


THE SELF-REFERENCE LADDER (Option B)
======================================

The Mersenne chain 2 -> 3 -> 7 -> 127 -> ... is:
  Start with the base k=2.
  The k-nacci matrix M_k at power k gives Tr(M_k^k) = 2^k - 1.
  This forbidden value becomes the NEW dimension.
  Repeat.

  Level 0: k = 2.  Base. The bit. F on (-1,1).
  Level 1: k = 3 = 2^2-1. Tribonacci. The transfer matrix. F on D.
  Level 2: k = 7 = 2^3-1. Heptanacci. The first forbidden value.
  Level 3: k = 127 = 2^7-1. The 127-nacci. Deeply hyperbolic.
  Level 4: k = 2^127-1 (Mersenne prime). The ??? -nacci.

Each level is a k-nacci system of dimension k.
The DIMENSION of the system at level n is the Mersenne chain value.

  dim(level 0) = 2
  dim(level 1) = 3
  dim(level 2) = 7
  dim(level 3) = 127
  dim(level n) = 2^{dim(level n-1)} - 1

This is the CATALAN-MERSENNE sequence. Each step:
  The system examines itself (Tr(M_k^k) = 2^k-1).
  The result of self-examination becomes the new dimension.
  The new system examines ITSELF at the new dimension.

THE d -> d+1 MAP IS SELF-EXAMINATION.
  The system of dimension d examines its own trace at power d.
  The result is 2^d - 1 (THM-227).
  This number becomes dimension d' = 2^d - 1.
  But this is NOT d+1 in general. It's d' = 2^d - 1 >> d.

So the Mersenne chain is NOT a d -> d+1 ladder.
It's a d -> 2^d-1 EXPONENTIAL ladder. Each step is enormous.
""")

# ========================================================================

print("=" * 60)
print("THE ACTUAL d -> d+1")
print("=" * 60)

print("""
Let me reconsider. What does d -> d+1 mean GEOMETRICALLY?

The formal group F(x,y) = (x+y)/(1+xy) on (-1,1) is dimension 1.
Its complexification to the Poincare disk is dimension 2.

The formal group at d=1 has k-nacci matrices of ANY size k.
The size k is NOT the same as the dimension of the formal group.
k is the ORDER OF MEMORY: how many predecessors the recurrence uses.

So we have TWO different "dimensions":
  d = the dimension of the formal group (1 real, 2 complex, ...)
  k = the memory depth (2 for Fibonacci, 3 for tribonacci, ...)

The Mersenne chain operates on k (memory depth).
The complexification operates on d (formal group dimension).

These are INDEPENDENT ladders:

  k-ladder: k=2 -> k=3 -> k=7 -> k=127 -> ...
    Self-examination. Exponential growth. Mersenne chain.
    At each step: the system discovers its forbidden value.
    The forbidden value becomes the memory of the next level.

  d-ladder: d=1 -> d=2 -> d=4 -> d=8 -> ...
    Complexification / normed division algebras.
    R -> C -> H -> O.
    At each step: introduce a new "imaginary" direction.
    Doubles the dimension.

And we found that {2, 4, 8} = {d=1 dims, d=2 dims, d=3 dims}
tiles the hyperbolic plane with area pi/8 = one Bott period.

THE BOTT LADDER is the d-ladder:
  d=1: R (reals). Formal group on (-1,1). Tournaments.
  d=2: C (complex). Formal group on D. Poincare disk.
  d=3: H (quaternions). Formal group on ... the unit ball in H?
  d=4: O (octonions). Formal group on ... the unit ball in O?

But H is noncommutative and O is non-associative.
The formal group F(x,y) = (x+y)/(1+xy) uses COMMUTATIVE multiplication.
It doesn't extend to H or O directly.

UNLESS: we modify F for noncommutativity.

For quaternions: F_H(q1, q2) = (q1 + q2)(1 + conj(q1)*q2)^{-1}.
  This IS the Mobius addition on the quaternionic unit ball.
  It's the isometry group of HYPERBOLIC 4-SPACE H^4_R (the quaternionic hyperbolic plane).

For octonions: F_O(o1, o2) = similar but non-associative.
  This is the isometry group of the CAYLEY HYPERBOLIC PLANE (16-dimensional).

So the d-ladder is:

  d=1: F on R^1 -> hyperbolic 1-space H^1 (the real hyperbolic line).
  d=2: F on C^1 -> hyperbolic 2-space H^2 (the Poincare disk).
  d=3: F_H on H^1 -> hyperbolic 4-space H^4_R (quaternionic hyperbolic plane).
  d=4: F_O on O^1 -> hyperbolic 8-space (Cayley hyperbolic plane).

Dimensions: 1, 2, 4, 8. POWERS OF 2. The Bott ladder!

The formal group extends to:
  1D: R with tanh.       Tournaments (binary comparisons).
  2D: C with tanh.       Poincare disk (boost + rotation).
  4D: H with tanh_H.     Quaternionic hyperbolic plane (3-rotations + 1 boost).
  8D: O with tanh_O.     Cayley hyperbolic plane (7 rotations + 1 boost).

And then it STOPS: the octonions are the last normed division algebra.
There is no d=5. The Bott ladder has 4 rungs: R, C, H, O.
(Then it repeats with period 8: Cl(n+8) = Cl(n) tensor M(16,R).)
""")

# ========================================================================

print("=" * 60)
print("THE TWO LADDERS COMBINED")
print("=" * 60)

print("""
We have TWO independent ladders:

  d-LADDER (Bott):     d = 1, 2, 4, 8, [period 8]
    Dimension of the formal group.
    Controlled by normed division algebras R, C, H, O.
    The formal group on R^d: F_d(x,y) = Mobius addition on hyperbolic d-space.

  k-LADDER (Mersenne): k = 2, 3, 7, 127, 2^127-1, ...
    Memory depth of the k-nacci system.
    Controlled by self-examination: Tr(M_k^k) = 2^k-1.
    The forbidden value of level n becomes the dimension of level n+1.

These are PERPENDICULAR AXES of a 2D structure:

         k=2    k=3    k=7    k=127
  d=1:  (1,2)  (1,3)  (1,7)  (1,127)     <-- real formal group at each memory
  d=2:  (2,2)  (2,3)  (2,7)  (2,127)     <-- complex formal group at each memory
  d=4:  (4,2)  (4,3)  (4,7)  (4,127)     <-- quaternionic formal group at each memory
  d=8:  (8,2)  (8,3)  (8,7)  (8,127)     <-- octonionic formal group at each memory

Each cell (d, k) is a FORMAL GROUP of dimension d with memory depth k.

  (1, 3) = the tournament theory we've been studying.
    Real formal group F on (-1,1), tribonacci transfer matrix M_3.
    Tr(M_3^3) = 7 = 2^3-1. The forbidden value.

  (2, 3) = the complexified tournament theory.
    Complex formal group F on D (Poincare disk), complex tribonacci.
    Same trace: Tr(M_3^3) = 7. (Trace doesn't depend on d.)

  (4, 3) = quaternionic tournament theory.
    Quaternionic formal group on the unit ball in H.
    Transfer matrix M_3 with quaternionic entries??
    Trace is STILL 7 (because trace is always over the base field).

  (8, 3) = octonionic tournament theory.
    Octonionic formal group on the unit ball in O.
    Non-associative! But the trace might still be 7.

THE KEY INSIGHT:

THM-227 (Tr(M_k^n) = 2^n-1 for n <= k) does NOT depend on d.
The trace is purely an algebraic identity about Newton's identities
and the all-(-1) coefficients. It works over ANY base field.

So the MERSENNE STRUCTURE is universal across all d-values.
The k-ladder is INDEPENDENT of the d-ladder.

But the d-ladder affects the GEOMETRY:
  d=1: the formal group is a line (real hyperbolic space).
  d=2: the formal group is a disk (complex hyperbolic space).
  d=4: the formal group is a 4-ball (quaternionic hyperbolic space).
  d=8: the formal group is an 8-ball (octonionic hyperbolic space).

The d -> d+1 map (actually d -> 2d) is:
  INTRODUCE A NEW IMAGINARY DIRECTION.
  The formal group acquires a new "rotation axis."
  The hyperbolic space doubles in dimension.

  d=1 -> d=2: introduce i. One rotation axis. The Poincare disk.
  d=2 -> d=4: introduce j. Three rotation axes (i,j,k). The quaternionic ball.
  d=4 -> d=8: introduce the octonion imaginaries. Seven rotation axes.

Each step DOUBLES the number of "comparison directions."

At d=1: a comparison is real: A > B or B > A. ONE BIT.
At d=2: a comparison is complex: A > B by amount x + iy. ONE COMPLEX BIT.
At d=4: a comparison is quaternionic: A > B with 3 rotational degrees.
At d=8: a comparison is octonionic: A > B with 7 rotational degrees.

THE d -> d+1 MAP IS: enrich the comparison with a new rotational degree of freedom.
""")

# ========================================================================

print("=" * 60)
print("THE FORMAL GROUP AT EACH d")
print("=" * 60)

# At d=1: F(x,y) = (x+y)/(1+xy) on (-1,1)
# At d=2: F(w1,w2) = (w1+w2)/(1+conj(w1)*w2) on D
# Wait, the correct Mobius addition on the disk is:
# F(w1,w2) = (w1+w2)/(1+conj(w1)*w2) -- this is the Poincare model formula
# But we've been using F(w1,w2) = (w1+w2)/(1+w1*w2) which works for w on real axis

print("""
The formal group at each d:

  d=1: F_1(x,y) = (x+y)/(1+xy)              on (-1,1) subset R
       Linearization: arctanh
       Metric: ds^2 = dx^2/(1-x^2)^2
       Isometry group: PSL(2,R) restricted to real axis = R (boosts)

  d=2: F_2(w1,w2) = (w1+w2)/(1+w1_bar*w2)   on D subset C
       Linearization: arctanh (complex)
       Metric: ds^2 = 4|dw|^2/(1-|w|^2)^2
       Isometry group: PSU(1,1) = PSL(2,R) = SL(2,R)/{+/-I}

  d=4: F_4(q1,q2) = (q1+q2)(1+q1_bar*q2)^{-1}  on B^4 subset H
       Linearization: arctanh_H (quaternionic)
       Metric: ds^2 = 4|dq|^2/(1-|q|^2)^2
       Isometry group: Sp(1,1)

  d=8: F_8(o1,o2) = (o1+o2)(1+o1_bar*o2)^{-1}  on B^8 subset O
       Linearization: arctanh_O (octonionic)
       Metric: ds^2 = 4|do|^2/(1-|o|^2)^2
       Isometry group: F_4(-20) (exceptional!)

At d=8, the isometry group is an EXCEPTIONAL LIE GROUP (F_4(-20)).
This is why the ladder stops: the octonions are non-associative,
and the Cayley hyperbolic plane is the last hyperbolic space
built from a normed division algebra.

THE FORMAL GROUP DIMENSION SEQUENCE: 1, 2, 4, 8.
These are the dimensions of R, C, H, O.
They are the powers of 2 up to 8.
They tile the hyperbolic plane as the {2,4,8} triangle (area = pi/8).
""")

# ========================================================================

print("=" * 60)
print("COMBINING THE TWO LADDERS: THE GRID")
print("=" * 60)

print("""
The complete structure is a GRID with:
  Horizontal axis: k (memory depth) = 2, 3, 7, 127, ...  (Mersenne chain)
  Vertical axis: d (formal group dim) = 1, 2, 4, 8      (Bott ladder)

At position (d, k):
  - A formal group F_d on the d-dimensional unit ball.
  - A k-nacci transfer matrix M_k (k x k, entries in R^d).
  - Trace Tr(M_k^n) = 2^n - 1 for n <= k (universal, independent of d).
  - The forbidden value 2^k - 1 (universal, independent of d).
  - The geometry is d-dimensional hyperbolic space H^d.

TOURNAMENT THEORY = the cell (d=1, k=3):
  Real formal group, tribonacci transfer matrix.
  H^1 = the real line.
  Forbidden value 7.

THE d -> d+1 MAP (actually d -> 2d):
  At each d, introduce a new imaginary unit.
  The formal group F_d extends to F_{2d}.
  The k-nacci matrix gets entries in a larger algebra.
  BUT THE TRACE DOESN'T CHANGE (it's an algebraic identity).

  So the d-ladder adds GEOMETRIC RICHNESS without changing ALGEBRAIC CONTENT.
  The Mersenne structure (forbidden values, universal zone) is constant.
  Only the SPACE in which the system lives changes.

THE k -> k' MAP (k -> 2^k - 1):
  At each k, the system examines its own trace.
  The result 2^k - 1 becomes the new memory depth.
  The formal group dimension d stays the same.
  The ALGEBRA becomes more complex (larger matrices, more eigenvalues).

  So the k-ladder adds ALGEBRAIC COMPLEXITY without changing GEOMETRIC DIMENSION.
  The formal group (the underlying space) is constant.
  Only the SYSTEM living in that space changes.

THE TWO LADDERS ARE PERPENDICULAR:
  d-ladder: changes GEOMETRY, preserves ALGEBRA.
  k-ladder: changes ALGEBRA, preserves GEOMETRY.

Together they span the full space of "comparison systems":
  (d, k) = a system with d geometric dimensions and k algebraic memory.

THE d -> d+1 MAP (really d -> 2d):

  STEP 1: Take the formal group F_d on the d-ball.
  STEP 2: Extend the base algebra: R -> C -> H -> O.
  STEP 3: The new formal group F_{2d} acts on the 2d-ball.
  STEP 4: A comparison now has d new rotational degrees of freedom.

THE k -> k' MAP (k -> 2^k - 1):

  STEP 1: Take the k-nacci transfer matrix M_k.
  STEP 2: Compute Tr(M_k^k) = 2^k - 1. The self-examination.
  STEP 3: This forbidden value becomes the new dimension k' = 2^k - 1.
  STEP 4: Build M_{k'}. The new system has exponentially more memory.

THE FORMAL GROUP'S "d NUMBER":

  d = 1:  the formal group F on (-1,1).  Area = 0 (it's 1D, no area).
  d = 2:  the formal group F on D.       The Poincare disk. Curvature = -1.
  d = 4:  the formal group F on B^4.     Quaternionic hyperbolic 4-space.
  d = 8:  the formal group F on B^8.     Octonionic Cayley hyperbolic plane.

So the formal group's "d number" is:
  d = dim_R(base algebra) = 1 for R, 2 for C, 4 for H, 8 for O.

And the map d -> d+1 is really d -> 2d (DOUBLING), which is:
  multiplication by 2 in the d-space.

But wait — multiplication by 2 in the k-space gives 2k, not 2^k-1.
While self-examination gives 2^k-1.

The d-ladder is MULTIPLICATIVE (d -> 2d).
The k-ladder is EXPONENTIAL (k -> 2^k - 1).

These are two different growth rates applied to two different axes.
The d-ladder is gentle (doubling: 1, 2, 4, 8).
The k-ladder is explosive (Mersenne: 2, 3, 7, 127, 2^127-1).

And the d-ladder TERMINATES (after O, it repeats with Bott period 8).
The k-ladder NEVER TERMINATES (the Mersenne chain goes on forever,
  as long as 2^k-1 keeps being prime — the Catalan-Mersenne conjecture).

THE FORMAL GROUP SITS AT d=2 (the Poincare disk) in its natural form.
To go from d to d+1:
  d=1 -> d=2: COMPLEXIFY (add i). This is what we did: extend F to the disk.
  d=2 -> d=4: QUATERNIONIFY (add j, k). Extend F to the 4-ball.
  d=4 -> d=8: OCTONIONIFY (add e_1,...,e_7). Extend F to the 8-ball.
  d=8 -> d=1 (mod 8): BOTT PERIODICITY. The cycle repeats.

The sequence of "new dimensions" at each step:
  d=1: 1 new direction (the real line itself).
  d=2: 1 new direction (the imaginary axis).
  d=4: 2 new directions (j and k).
  d=8: 4 new directions (the 4 new octonionic imaginaries).

New directions added: 1, 1, 2, 4. Powers of 2 starting from 1.
Cumulative: 1, 2, 4, 8. The d-values themselves.

Each step DOUBLES the number of new directions.
Each step introduces AS MANY new directions as already existed.
The comparison system DOUBLES its rotational freedom at each step.
""")

print()
print("opus-2026-03-16-S90bd: THE DIMENSION LADDER")
