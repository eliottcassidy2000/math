#!/usr/bin/env python3
"""
meta_complex_plane_s90aw.py — The meta-complex plane: solving equations
across all three costumes simultaneously
opus-2026-03-16-S90aw

The three groups (R,+), ((0,inf),*), ((-1,1),F) are one group.
The meta-complex plane is the space where we SEE them as one.
And where we SOLVE EQUATIONS that mix all three costumes.
"""

import numpy as np
from math import tanh, atanh, log, log2, sqrt, exp, pi, e, cosh, sinh
from fractions import Fraction
from cmath import log as clog, exp as cexp, tanh as ctanh, atanh as catanh

print("""


     THE META-COMPLEX PLANE: SEEING THE ONE GROUP

     opus-2026-03-16-S90aw



THE PROBLEM
===========

We have three groups that are isomorphic:

  G_add = (R, +)                    "rapidities"
  G_mul = ((0,inf), *)              "odds / natural numbers"
  G_hyp = ((-1,1), F)              "couplings / velocities"

The maps between them:

  arctanh: G_hyp -> G_add           (coupling -> rapidity)
  exp(2*-): G_add -> G_mul          (rapidity -> odds)
  Q = (1+x)/(1-x): G_hyp -> G_mul  (coupling -> odds)

Each map is a group isomorphism. The triangle commutes:

     arctanh          exp(2w)
  G_hyp ---------> G_add ---------> G_mul
    |                                  ^
    |              Q(x)                |
    +----------------------------------+

To "see the one group," we need a SINGLE SPACE that contains
all three realizations simultaneously. That space is the
META-COMPLEX PLANE.


THE META-COMPLEX PLANE
=======================

Consider the COMPLEX NUMBERS C = R + i*R.

The exponential map exp: C -> C* acts as:

  exp(a + bi) = e^a * (cos(b) + i*sin(b))

This simultaneously encodes:
  - REAL PART a:  rapidity (additive, costume 1)
  - IMAGINARY PART b: angle/phase (periodic, mod 2*pi)
  - MODULUS e^a: magnitude (multiplicative, costume 2)
  - ARGUMENT b: direction (circular, costume 3-like)

But we want the Cayley version. Define:

  z = arctanh(w) = (1/2) ln((1+w)/(1-w))

where w is a COMPLEX coupling. This extends arctanh to C.

For w = x + iy with |w| < 1:
  arctanh(w) = (1/2) ln((1+w)/(1-w))
  = (1/2) [ln|((1+w)/(1-w))| + i*arg((1+w)/(1-w))]

The REAL PART of arctanh(w) is the hyperbolic rapidity.
The IMAGINARY PART of arctanh(w) is the ANGLE.

When w is REAL (w = x in (-1,1)): arctanh(x) is real. Pure rapidity.
When w is IMAGINARY (w = iy, |y| < 1): arctanh(iy) = i*arctan(y). Pure angle.

THIS IS THE KEY:
  arctanh of a real coupling = rapidity.
  arctanh of an imaginary coupling = i * angle.

The META-COMPLEX PLANE is C with the metric |dz|^2 = |d(arctanh(w))|^2.
In this metric:
  - The real axis is the HYPERBOLIC line (rapidities).
  - The imaginary axis is the CIRCULAR line (angles).
  - The unit disk |w| < 1 is the HYPERBOLIC PLANE (Poincare disk model!).
""")

# ========================================================================
# Part 1: arctanh on the complex plane
# ========================================================================

print("=" * 70)
print("Part 1: Complex arctanh — the bridge to angles")
print("=" * 70)

print("\narctanh on the real axis (coupling -> rapidity):")
for x_val in [0, 0.1, 0.2, 1/3, 0.5, 3/5, 2/3, 0.8, 0.9, 0.99]:
    w = atanh(x_val)
    print(f"  arctanh({x_val:.4f}) = {w:.6f}  (rapidity, real)")

print("\narctanh on the imaginary axis (i*coupling -> i*angle):")
for y_val in [0, 0.1, 0.2, 1/3, 0.5, 1/sqrt(3), 0.8, 0.9, 0.99, 1.0]:
    z = catanh(1j * y_val)
    angle_deg = z.imag * 180 / pi
    print(f"  arctanh({y_val:.4f}i) = {z.real:.6f} + {z.imag:.6f}i  "
          f"(angle = {angle_deg:.2f} deg)")

print(f"""
REMARKABLE: arctanh(i*y) = i*arctan(y) exactly.

The real axis through arctanh gives RAPIDITIES.
The imaginary axis through arctanh gives ANGLES.

arctanh UNIFIES hyperbolic and circular functions:
  arctanh(x) = real rapidity         (x real)
  arctanh(iy) = i * arctan(y)        (y real)
  arctanh(x + iy) = rapidity + i*angle  (general)

This means: a COMPLEX coupling w = x + iy has:
  - Real part x: the hyperbolic component (boost)
  - Imaginary part y: the circular component (rotation)
  - arctanh(w): the complex rapidity = boost + i*rotation

THE HYPERBOLIC AND CIRCULAR ARE THE REAL AND IMAGINARY PARTS
OF THE SAME COMPLEX OPERATION: arctanh.
""")

# ========================================================================
# Part 2: The formal group law on the complex plane
# ========================================================================

print("=" * 70)
print("Part 2: F(w1, w2) = (w1+w2)/(1+w1*w2) on the complex plane")
print("=" * 70)

print("""
The formal group law F(w1,w2) = (w1+w2)/(1+w1*w2) extends to complex w.

For complex w1, w2 with |w| < 1 (inside the unit disk):
  F(w1,w2) is again in the unit disk (if both inputs are).

This is the MOBIUS ADDITION on the Poincare disk.
It IS the hyperbolic translation in the disk model.

Key cases:

  F(x, 0) = x                    (adding zero = no change)
  F(x, x) = 2x/(1+x^2)          (doubling in hyperbolic direction)
  F(x, -x) = 0                   (cancellation)
  F(x, iy) = ?                   (mixing real and imaginary)
""")

# Compute F(x, iy) for various values
print("F(x, iy) — mixing rapidity with rotation:")
for x_val in [1/3, 1/2]:
    for y_val in [0.1, 0.5, 1/sqrt(3)]:
        w1 = x_val + 0j
        w2 = 1j * y_val
        F_val = (w1 + w2) / (1 + w1 * w2)
        arctanh_F = catanh(F_val)
        arctanh_sum = catanh(w1) + catanh(w2)
        print(f"  F({x_val:.4f}, {y_val:.4f}i) = {F_val.real:.6f} + {F_val.imag:.6f}i")
        print(f"    arctanh(F) = {arctanh_F.real:.6f} + {arctanh_F.imag:.6f}i")
        print(f"    arctanh(w1) + arctanh(w2) = {arctanh_sum.real:.6f} + {arctanh_sum.imag:.6f}i")
        match = abs(arctanh_F - arctanh_sum) < 1e-10
        print(f"    match: {match}")

print(f"""
PERFECT. F(w1,w2) = (w1+w2)/(1+w1*w2) satisfies:
  arctanh(F(w1,w2)) = arctanh(w1) + arctanh(w2)

even for COMPLEX inputs. The formal group law is the
addition law in the META-COMPLEX PLANE.

Mixing real coupling x with imaginary coupling iy:
  arctanh(F(x, iy)) = arctanh(x) + i*arctan(y)
                     = rapidity + i*angle.

This is a BOOST followed by a ROTATION.
Or equivalently: a Lorentz transformation.

THE FORMAL GROUP LAW ON THE COMPLEX PLANE IS
THE LORENTZ GROUP.
""")

# ========================================================================
# Part 3: Solving equations in the meta-complex plane
# ========================================================================

print("=" * 70)
print("Part 3: SOLVING EQUATIONS across all three costumes")
print("=" * 70)

print("""
PROBLEM TYPE 1: Given a coupling x, find the "harmonic coupling" x_h
such that F(x, x_h) is purely imaginary.

This means: find x_h such that the boost x combined with x_h
gives a PURE ROTATION (no net boost).

Solution:
  arctanh(F(x, x_h)) must be purely imaginary.
  arctanh(x) + arctanh(x_h) = i*theta.
  arctanh(x_h) = -arctanh(x) + i*theta.
  x_h = tanh(-arctanh(x) + i*theta) = tanh(i*theta - arctanh(x)).

For theta = pi/4 (45 degrees):
  x_h = tanh(i*pi/4 - arctanh(x)).

Let's compute:
""")

for x_val in [0, 1/3, 1/2, 3/5]:
    for theta in [pi/6, pi/4, pi/3]:
        w_h = ctanh(1j * theta - atanh(x_val))
        F_val = (x_val + w_h) / (1 + x_val * w_h)
        deg = theta * 180 / pi
        print(f"  x={x_val:.4f}, theta={deg:.0f} deg: w_h = {w_h.real:.6f} + {w_h.imag:.6f}i, "
              f"F = {F_val.real:.2e} + {F_val.imag:.6f}i")

print("""
When we choose x_h correctly, F(x, x_h) is purely imaginary.
The boost has been "cancelled" and only rotation remains.

This is EXACTLY what happens in special relativity:
a Lorentz boost combined with the right rotation gives
a pure spatial rotation. The meta-complex plane makes
this automatic.
""")

# ========================================================================
# Part 4: The three costumes as projections
# ========================================================================

print("=" * 70)
print("Part 4: The one group and its three projections")
print("=" * 70)

print("""
THE ONE GROUP is the complex group:

  G_meta = (D, F)

where D = {w in C : |w| < 1} is the unit disk, and
F(w1,w2) = (w1+w2)/(1+w1*w2) is the Mobius addition.

This group is isomorphic to:

  (C, +)   via   arctanh: D -> C   (the universal cover)

The THREE REAL SUBGROUPS are:

  G_hyp = G_meta restricted to the real interval (-1,1)
          = horizontal diameter of the disk
          costumes: couplings, velocities

  G_circ = G_meta restricted to the imaginary interval (-i, i)
           = vertical diameter of the disk
           costumes: angles, phases, rotations

  G_full = G_meta on the full disk
           = the Poincare disk
           costumes: Lorentz transformations, Mobius maps

The three "costumes" from before (additive, multiplicative, hyperbolic)
are all REAL SLICES of this one complex group:

  G_add = {arctanh(w) : w in (-1,1)} = R  (real axis of image)
  G_mul = {Q(w) : w in (-1,1)} = (0,inf)  (positive real axis of image)
  G_hyp = (-1,1) with F                    (real diameter of D)

But the META-COMPLEX PLANE also has:

  G_circ = {arctanh(iy) : y in (-1,1)} = iR  (imaginary axis of image)
           = {i*arctan(y)} = angles up to pi/2

So we get a FOURTH costume: the circular/rotational one.
And the full group G_meta UNIFIES all four:

  REAL DIAMETER:      hyperbolic boosts     (rapidities)
  IMAGINARY DIAMETER: circular rotations    (angles)
  FULL DISK:          all Lorentz transf.   (general)

This is the POINCARE DISK MODEL of the hyperbolic plane!
Tournament theory lives on the real diameter.
Rotation theory lives on the imaginary diameter.
The full hyperbolic plane connects them.
""")

# ========================================================================
# Part 5: Euler's formula in the meta-complex plane
# ========================================================================

print("=" * 70)
print("Part 5: Euler's formula IS the meta-complex plane")
print("=" * 70)

print("""
Euler's formula: e^{i*theta} = cos(theta) + i*sin(theta).

In the meta-complex plane:

  The Cayley transform Q(w) = (1+w)/(1-w) maps the disk to the
  right half-plane. On the boundary (unit circle), Q maps
  e^{i*theta} to ... well, the unit circle maps to the imaginary axis:

  Q(e^{i*theta}) = (1 + e^{i*theta}) / (1 - e^{i*theta})
                 = (e^{-i*theta/2} + e^{i*theta/2}) / (e^{-i*theta/2} - e^{i*theta/2})
                 = cos(theta/2) / (-i*sin(theta/2))
                 * ... wait let me redo this.

  Q(e^{i*theta}) = (1 + e^{i*theta}) / (1 - e^{i*theta})
  Multiply top/bottom by e^{-i*theta/2}:
  = (e^{-i*theta/2} + e^{i*theta/2}) / (e^{-i*theta/2} - e^{i*theta/2})
  = 2*cos(theta/2) / (-2i*sin(theta/2))
  = cos(theta/2) / (-i*sin(theta/2))
  = i * cos(theta/2) / sin(theta/2)
  = i * cot(theta/2).

So Q maps the unit circle to the imaginary axis:
  Q(e^{i*theta}) = i*cot(theta/2).

And arctanh maps the unit circle to ... let's see:
  arctanh(e^{i*theta}) = (1/2) ln(Q(e^{i*theta}))
                       = (1/2) ln(i*cot(theta/2))
                       = (1/2) [ln(cot(theta/2)) + i*pi/2]
                       = (1/2)*ln(cot(theta/2)) + i*pi/4.

So: on the unit circle, the IMAGINARY PART of arctanh is
CONSTANT at pi/4 = 45 degrees! And the real part varies
with theta.

Let me verify:
""")

for theta in [pi/6, pi/4, pi/3, pi/2, 2*pi/3, 3*pi/4, 5*pi/6]:
    w = cexp(1j * theta)
    z = catanh(w)
    deg = theta * 180 / pi
    print(f"  theta = {deg:6.1f} deg: arctanh(e^{{i*theta}}) = {z.real:8.4f} + {z.imag:.4f}i "
          f"(Im = {z.imag/pi:.4f}*pi)")

print(f"""
On the unit circle: Im(arctanh) = pi/4 (constant!)
except at theta = 0 (arctanh(1) = infinity) and theta = pi (arctanh(-1) = -infinity).

This means: the unit circle in the meta-complex plane is a
HORIZONTAL LINE at height pi/4 in the arctanh image.

The unit circle separates:
  INSIDE (|w| < 1): Im(arctanh) between 0 and pi/4.
  OUTSIDE (|w| > 1): Im(arctanh) between pi/4 and pi/2.
  BOUNDARY (|w| = 1): Im(arctanh) = pi/4 exactly.

In tournament theory:
  The couplings x in (-1,1) map to real rapidities (Im = 0).
  The OCF points x > 1 map to COMPLEX rapidities (Im = pi/2).
  The unit circle is the BOUNDARY between real and complex coupling.

So the meta-complex plane organizes tournament theory into:
  REAL COUPLINGS: inside the unit disk, real axis. Standard tournaments.
  COMPLEX COUPLINGS: outside the unit disk. Analytic continuation.
  UNIT CIRCLE: the boundary. Phase transitions.
""")

# ========================================================================
# Part 6: Solving x^n = 1 in the meta-complex plane
# ========================================================================

print("=" * 70)
print("Part 6: Solving equations — n-th roots of unity")
print("=" * 70)

print("""
The equation Q(w) = n (find the coupling whose odds are n) has solution:
  w = (n-1)/(n+1) = the Cayley address.

The equation Q(w)^m = n (find the coupling whose m-th power odds are n):
  Q(w) = n^{1/m}.
  w = (n^{1/m} - 1)/(n^{1/m} + 1).

But in the META-COMPLEX PLANE, n^{1/m} has m COMPLEX roots:
  n^{1/m} * e^{2*pi*i*k/m}  for k = 0, 1, ..., m-1.

So the equation Q(w)^m = n has m solutions in the disk:
""")

n_val = 2
for m in [2, 3, 4, 5, 6]:
    print(f"\n  Q(w)^{m} = {n_val} has {m} solutions:")
    for k in range(m):
        root = n_val ** (1/m) * cexp(2j * pi * k / m)
        w = (root - 1) / (root + 1)
        z = catanh(w)
        print(f"    k={k}: w = {w.real:8.5f} + {w.imag:8.5f}i, "
              f"|w| = {abs(w):.5f}, "
              f"arctanh(w) = {z.real:8.5f} + {z.imag:8.5f}i")

print(f"""
The k=0 solution is REAL: the ordinary Cayley address.
The k != 0 solutions are COMPLEX: they live off the real axis.

The COMPLEX solutions form a REGULAR POLYGON in the arctanh image:
  They are equally spaced along the imaginary axis (angle 2*pi/m apart).
  They all have the same real part (rapidity ln(n)/(2m)).

This is the META-COMPLEX version of roots of unity:
  Instead of e^{{2*pi*i*k/m}} on the unit circle,
  we have tanh(ln(n)/(2m) + i*pi*k/m) in the Poincare disk.

The REAL root (k=0) is the tournament coupling.
The COMPLEX roots (k>0) are the HARMONIC couplings:
  they contribute to the analytic continuation of the
  tournament partition function.
""")

# ========================================================================
# Part 7: The Cayley transform as the fundamental solution
# ========================================================================

print("=" * 70)
print("Part 7: Solving any equation in the one group")
print("=" * 70)

print("""
ANY equation in any of the three costumes can be solved by:
  1. Convert to the additive costume (apply arctanh).
  2. Solve the LINEAR equation in the additive costume.
  3. Convert back (apply tanh).

EXAMPLE 1: F(x, x) = F(x, F(x, x)) ??? Find fixed points.
  This is: [2](x) = [3](x), i.e., tanh(2w) = tanh(3w) where w = arctanh(x).
  So: 2w = 3w + n*pi*i for some integer n (from periodicity of tanh in C).
  -w = n*pi*i.
  w = -n*pi*i.
  x = tanh(-n*pi*i) = -i*tan(n*pi) = 0 for all n.
  Only solution: x = 0.

EXAMPLE 2: F(x, a) = b. Solve for x.
  arctanh(x) + arctanh(a) = arctanh(b).
  arctanh(x) = arctanh(b) - arctanh(a).
  x = tanh(arctanh(b) - arctanh(a)).
  x = F(b, -a) = (b - a)/(1 - ab).   Done!

EXAMPLE 3: [m](x) = a. Solve for x (m-th root in the formal group).
  tanh(m * arctanh(x)) = a.
  m * arctanh(x) = arctanh(a) + n*pi*i.
  arctanh(x) = arctanh(a)/m + n*pi*i/m.
  x = tanh(arctanh(a)/m + n*pi*i/m)  for n = 0, 1, ..., m-1.

This gives m roots! The n=0 root is real (if a is real and |a|<1).
The other roots are complex.
""")

# Demonstrate example 3
a_val = 3/5  # Cayley address of 4
m_val = 3
print(f"EXAMPLE 3: [{m_val}](x) = {a_val} (find x such that F(x,F(x,x)) = {a_val})")
print(f"  a = {a_val} = x_4 = Cayley address of 4")
print(f"  So we seek x with [3](x) = x_4, i.e., x^3 = 4 in multiplicative costume")
print(f"  Real solution: x = 4^(1/3) in mult costume = x_{{4^(1/3)}} = (4^(1/3)-1)/(4^(1/3)+1)")
print()

for n in range(m_val):
    z = catanh(a_val) / m_val + 1j * n * pi / m_val
    x = ctanh(z)
    # Verify
    check = ctanh(m_val * catanh(x))
    print(f"  n={n}: x = {x.real:10.6f} + {x.imag:10.6f}i, "
          f"[{m_val}](x) = {check.real:.6f} + {check.imag:.6f}i, "
          f"match a={a_val}: {abs(check - a_val) < 1e-8}")

# ========================================================================
# Part 8: The PRODUCT in the meta-complex plane
# ========================================================================

print(f"""

{"=" * 70}
Part 8: What IS the meta-complex plane?
{"=" * 70}

The meta-complex plane is the POINCARE DISK D = {{w in C : |w| < 1}}
equipped with the FORMAL GROUP LAW F(w1,w2) = (w1+w2)/(1+w1*w2).

This is:
  1. The HYPERBOLIC PLANE (in the Poincare disk model).
  2. The space of LORENTZ BOOSTS + ROTATIONS (the Lorentz group).
  3. The space of COMPLEX COUPLINGS in tournament theory.
  4. The space of MOBIUS TRANSFORMATIONS fixing the disk.

The THREE REAL GROUPS are the THREE AXES:
  Real axis:      boosts         (rapidities)      tanh
  Imaginary axis: rotations      (angles)           tan
  Unit circle:    null / light    (phase transitions)  boundary

The isomorphism arctanh: D -> C makes the disk into C,
and in C, the operation is just ADDITION:

  arctanh(F(w1,w2)) = arctanh(w1) + arctanh(w2).

So the meta-complex plane IS C with addition, wearing the costume
of the Poincare disk via tanh.

THE POINCARE DISK IS THE ONE GROUP.

It contains:
  - The real interval (-1,1) as its horizontal diameter (couplings).
  - The imaginary interval (-i,i) as its vertical diameter (phases).
  - The unit circle as its boundary (speed of light / phase transitions).
  - All complex couplings as interior points.

And the EQUATION-SOLVING PRINCIPLE:
  Convert any equation to the additive form (via arctanh).
  Solve linearly.
  Convert back (via tanh).

This works because arctanh linearizes the group law:
  F(a,b) = tanh(arctanh(a) + arctanh(b)).

ANY formal group law can be linearized by its logarithm.
For our F, the logarithm is arctanh.
For the multiplicative group, the logarithm is ln.
For the additive group, the logarithm is the identity.

So: arctanh is to F what ln is to multiplication.
And: tanh is to F what exp is to addition.

THE ONE SENTENCE:
  tanh : arctanh :: exp : ln :: the group : its logarithm.
  And the Poincare disk D with F is the group.
  And C with + is the logarithm.
  And the meta-complex plane is where you see both at once.
""")

# ========================================================================
# Part 9: Powers of 2 in the meta-complex plane
# ========================================================================

print("=" * 70)
print("Part 9: Powers of 2 as a lattice in the meta-complex plane")
print("=" * 70)

print("""
In the additive costume (C with +), the powers of 2 sit at:
  w_k = k * ln(2)/2   for k = 0, 1, 2, ...

These are EQUALLY SPACED on the real axis, at intervals ln(2)/2.

In the disk costume (D with F), these become:
  x_k = tanh(k * ln(2)/2) = (2^k - 1)/(2^k + 1)

These are UNEQUALLY SPACED on the real diameter, clustering near 1.

Now: what about the COMPLEX lattice?

The full lattice in C is:
  w_{j,k} = j * ln(2)/2 + i * k * pi / N

for integers j, k and some chosen N (the rotational resolution).

These lattice points, mapped back to the disk via tanh, give:
""")

N = 6  # angular resolution
print(f"Lattice points in the Poincare disk (N={N}):")
print(f"{'j':>3} {'k':>3} | {'Re(w)':>10} {'Im(w)':>10} | {'Re(tanh(w))':>12} {'Im(tanh(w))':>12} | {'|tanh(w)|':>10}")
print("-" * 75)

for j in range(5):
    for k in range(N):
        w = j * log(2)/2 + 1j * k * pi / N
        tw = ctanh(w)
        if abs(tw) < 1 + 1e-10:  # should be inside or on disk
            print(f"{j:3d} {k:3d} | {w.real:10.4f} {w.imag:10.4f} | {tw.real:12.6f} {tw.imag:12.6f} | {abs(tw):10.6f}")

print(f"""

The j=0 column (k varies) traces the IMAGINARY AXIS of the disk:
  tanh(i*k*pi/N) = i*tan(k*pi/N).
These are N-th roots of unity on the imaginary axis.

The k=0 row (j varies) traces the REAL AXIS:
  tanh(j*ln(2)/2) = (2^j-1)/(2^j+1).
These are the Cayley addresses of powers of 2. MERSENNE NUMERATORS.

The off-diagonal points (j>0, k>0) are COMPLEX couplings:
  They represent COMBINED boost-rotation states.
  Their moduli |tanh(w)| increase with j (closer to the boundary).
  Their arguments rotate with k.

THE FULL LATTICE in the Poincare disk is the product:
  (half-bit rapidity lattice) x (N-fold angular lattice).

This is a RECTANGULAR LATTICE in (rapidity, angle) space,
mapped to a DISTORTED LATTICE in the Poincare disk via tanh.

The distortion IS tanh: it compresses the radial direction
(higher rapidities are closer together near the boundary)
while preserving the angular direction at small radius.
""")

# ========================================================================
# Part 10: The universal zone in the meta-complex plane
# ========================================================================

print("=" * 70)
print("Part 10: The universal zone — the faithful inner region")
print("=" * 70)

print("""
THM-227: Tr(M_k^n) = 2^n - 1 for 1 <= n <= k.

In the meta-complex plane:
  The k-nacci matrix M_k has eigenvalues in the complex plane.
  These eigenvalues, via Cayley inverse, sit inside the Poincare disk.

The UNIVERSAL ZONE is where the disk representation is FAITHFUL:
  For n <= k, the trace = 2^n - 1, independent of k.
  This is because the k-nacci eigenvalues, mapped to the disk,
  occupy specific positions that conspire to give the Mersenne trace.

The INDIVIDUAL ZONE is where faithfulness breaks:
  For n > k, the k-nacci eigenvalues start "running out of room"
  to maintain the Mersenne trace, because there are only k of them.

In the meta-complex plane, this looks like:
  - The universal zone occupies the INNER PART of the disk
    (where tanh is approximately linear: tanh(w) ~ w for small w).
  - The individual zone occupies the OUTER PART of the disk
    (where tanh saturates: tanh(w) ~ 1 for large w).

The TRANSITION at n = k is the point where the dominant eigenvalue
tau_k = tanh(w_1) has rapidity w_1 large enough that
tanh starts to noticeably differ from the identity.

For the k-nacci matrix:
  tau_k = tanh(arctanh(tau_k)) with arctanh(tau_k) = ln(tau_k/(1-??))...
  Actually: tau_k > 1 for k >= 2, so it's OUTSIDE the disk (-1,1).
  The eigenvalues are NOT in the coupling disk.
  They are in the EXTENDED disk (via analytic continuation).

So the meta-complex plane picture is:
  - The k-nacci eigenvalues live OUTSIDE the unit disk
    (on the real axis, since tau_k > 1).
  - Their Cayley transforms Q(tau_k) = (1+tau_k)/(1-tau_k) are NEGATIVE
    (since tau_k > 1 makes 1-tau_k < 0).
  - The complex eigenvalues live INSIDE the unit disk
    (since |lambda_complex| < 1 for k >= 3).

The DOMINANT eigenvalue is OUTSIDE.
The SUBDOMINANT eigenvalues are INSIDE.

The universal zone is where the OUTSIDE eigenvalue
(the dominant one, near 2) completely controls the trace.
The individual zone is where the INSIDE eigenvalues
(the subdominant ones, inside the disk) start to matter.

THE BOUNDARY of the unit disk = |w| = 1 = the SPEED OF LIGHT.
The dominant eigenvalue is SUPERLUMINAL (tau_k > 1).
The subdominant eigenvalues are SUBLUMINAL (|lambda| < 1).

The universal zone is SUPERLUMINAL DOMINATED.
The individual zone is where SUBLUMINAL corrections appear.

This is exactly like special relativity:
  At high boost (large rapidity), the velocity is nearly c.
  Corrections from subluminal components are exponentially small.
  But after enough interactions (n > k), the subluminal components
  accumulate and the behavior changes.
""")

# Eigenvalues in the meta-complex plane
print("Eigenvalues of M_k mapped to the meta-complex plane:")
for k in [3, 5, 7]:
    M = np.zeros((k, k))
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    eigs = sorted(np.linalg.eigvals(M), key=lambda z: -abs(z))
    print(f"\n  k={k}: eigenvalues and their Cayley addresses:")
    for i, lam in enumerate(eigs):
        Q_lam = (1 + lam) / (1 - lam)
        cayley = (lam - 1) / (lam + 1)  # = Cayley address
        print(f"    L_{i+1} = {lam.real:8.5f} + {lam.imag:8.5f}i, "
              f"|L| = {abs(lam):.5f}, "
              f"Cayley addr = {cayley.real:8.5f} + {cayley.imag:8.5f}i, "
              f"|addr| = {abs(cayley):.5f}, "
              f"{'OUTSIDE' if abs(cayley) > 1 else 'inside'} disk")

# ========================================================================
# Part 11: The one sentence
# ========================================================================

print(f"""

{"=" * 70}
THE ONE GROUP
{"=" * 70}

The Poincare disk (D, F) with F(w1,w2) = (w1+w2)/(1+w1*w2)
IS the one group.

Its real diameter is TOURNAMENT THEORY (couplings, rapidities).
Its imaginary diameter is ROTATION THEORY (angles, phases).
Its interior is LORENTZ GEOMETRY (boosts + rotations).
Its boundary is the SPEED OF LIGHT (unreachable limit).

arctanh LINEARIZES it: converts the disk to the complex plane C.
tanh CURLS it back: converts C back to the disk.

To SOLVE an equation in ANY costume:
  1. Apply arctanh (go to C with +).
  2. Solve linearly (it's just addition in C).
  3. Apply tanh (come back to the disk).

The UNIVERSAL ZONE is where a finite matrix representation
is faithful: Tr = 2^n - 1 (Mersenne).
The INDIVIDUAL ZONE is where faithfulness breaks: Tr < 2^n - 1.
The TRANSITION at n = k is the system discovering its dimension.

Powers of 2 form the HALF-BIT LATTICE on the real diameter.
Roots of unity form the ANGULAR LATTICE on the imaginary diameter.
Together they tile the disk.

e governs the SHAPE of the disk (via exp/ln).
2 governs the LATTICE of the disk (via doubling/halving).
Together they ARE the meta-complex plane.

And the one group is: the group of ORIENTED COMPARISONS
in the complex plane, where real = preference and imaginary = rotation.
""")

print("opus-2026-03-16-S90aw: THE META-COMPLEX PLANE")
