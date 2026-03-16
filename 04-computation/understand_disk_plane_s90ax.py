#!/usr/bin/env python3
"""
understand_disk_plane_s90ax.py — Understand the Poincare disk and the complex plane
opus-2026-03-16-S90ax

The Poincare disk and the complex plane are the SAME SPACE
seen through two different lenses. This script makes that precise.
"""

from cmath import exp as cexp, log as clog, tanh as ctanh, atanh as catanh
from cmath import pi as cpi, sqrt as csqrt
from math import tanh, atanh, log, sqrt, pi, e, cos, sin, atan2

print("""


     UNDERSTAND THE POINCARE DISK AND THE COMPLEX PLANE

     opus-2026-03-16-S90ax



THE COMPLEX PLANE, FIRST.
==========================

The complex plane C is R^2 with extra structure:
  a point z = x + iy is a pair (x, y) that KNOWS HOW TO MULTIPLY.

  (a + bi)(c + di) = (ac - bd) + (ad + bc)i.

Multiplication by a complex number does TWO things simultaneously:
  SCALES by |z| = sqrt(x^2 + y^2)    (the modulus)
  ROTATES by arg(z) = arctan(y/x)     (the argument)

This is why the complex plane matters: it makes
"scale and rotate" into a single operation.

In polar form: z = r * e^{i*theta}.
  Multiplication: z1 * z2 = r1*r2 * e^{i*(theta1 + theta2)}.
  Moduli MULTIPLY. Arguments ADD.

So the complex plane contains BOTH:
  (0,inf) with * (the moduli — multiplicative group)
  [0,2*pi) with + mod 2*pi (the arguments — circular group)

as its two "coordinates" in polar form.

THE COMPLEX PLANE IS THE PRODUCT:

  C* = (0,inf) x S^1
     = (multiplicative magnitudes) x (angles)
     = (radii) x (directions)

And the map exp: C -> C* sends:
  a + bi  |-->  e^a * e^{ib} = e^a * (cos b + i sin b)

So: REAL PART a  ->  MAGNITUDE e^a
    IMAGINARY PART b  ->  ANGLE b

exp converts ADDITIVE structure to MULTIPLICATIVE-CIRCULAR structure.
ln converts back.

This is the FIRST duality: additive vs multiplicative-circular.


THE POINCARE DISK, NEXT.
=========================

The Poincare disk D = {w in C : |w| < 1} is the open unit disk.

As a SET, it's just the inside of the unit circle.
But it has a DIFFERENT GEOMETRY than the flat complex plane.

The flat (Euclidean) metric on C: ds^2 = |dz|^2 = dx^2 + dy^2.
The Poincare (hyperbolic) metric on D: ds^2 = 4|dw|^2 / (1 - |w|^2)^2.

In the Poincare metric:
  - The CENTER (w = 0) is an ordinary point.
  - The BOUNDARY (|w| = 1) is INFINITELY FAR from the center.
  - Straight lines (geodesics) are arcs of circles perpendicular to the boundary.
  - The disk has CONSTANT NEGATIVE CURVATURE K = -1.

The REAL DIAMETER (-1, 1) of the Poincare disk
is a geodesic — the "horizontal line through the center."

Points on this geodesic are EXACTLY the real couplings x in (-1,1)
from tournament theory.

The IMAGINARY DIAMETER (-i, i) is another geodesic —
the "vertical line through the center."

Points on this geodesic are purely imaginary couplings iy.


THE DISTANCE FORMULA.
======================

The hyperbolic distance from the center 0 to a point w in D is:

  d(0, w) = 2 * arctanh(|w|).

So: ARCTANH IS THE DISTANCE FUNCTION.

For a point on the real axis, w = x:
  d(0, x) = 2 * arctanh(x) = 2 * rapidity.

For a point on the imaginary axis, w = iy:
  d(0, iy) = 2 * arctanh(|y|) = 2 * arctanh(y).

The hyperbolic distance grows WITHOUT BOUND as |w| -> 1,
even though the Euclidean distance stays bounded (|w| < 1).

THIS IS WHAT THE POINCARE DISK DOES:
  It takes the INFINITE hyperbolic plane and COMPRESSES it
  into the finite unit disk, using tanh as the compression.

  tanh: (unbounded distance) -> (bounded disk coordinate)
  arctanh: (bounded disk coordinate) -> (unbounded distance)

The Poincare disk is the PHYSICAL UNIVERSE seen from inside:
  - You can travel forever (unbounded distance).
  - But everything fits in a disk (bounded coordinates).
  - The "edge" at |w| = 1 is infinitely far away.
  - This is EXACTLY the speed of light: approachable, never reachable.
""")

# ========================================================================
# Part 1: The distance table
# ========================================================================

print("=" * 70)
print("Part 1: Hyperbolic distance from center to key points")
print("=" * 70)

print(f"\n{'w':>12} | {'|w|':>8} | {'d(0,w) = 2*atanh(|w|)':>22} | {'interpretation':>30}")
print("-" * 80)

points = [
    (0.0, "rest"),
    (1/3, "x_2 (Cayley addr of 2)"),
    (1/2, "x_3 (Cayley addr of 3)"),
    (3/5, "x_4 (Cayley addr of 4)"),
    (2/3, "x_5 (Cayley addr of 5)"),
    (5/7, "x_6"),
    (3/4, "x_7"),
    (7/9, "x_8"),
    (0.9, "approaching boundary"),
    (0.99, "near boundary"),
    (0.999, "very near boundary"),
]

for val, name in points:
    dist = 2 * atanh(val) if val < 1 else float('inf')
    print(f"  {val:10.6f} | {val:8.6f} | {dist:22.8f} | {name:>30}")

print("""
The Cayley address x_n = (n-1)/(n+1) is the DISK COORDINATE of n.
The hyperbolic distance 2*arctanh(x_n) = ln(n) is the LOG of n.

So: the natural numbers, placed at their Cayley addresses in the Poincare disk,
are at EQUAL LOGARITHMIC SPACING (distance ln(2) between consecutive powers of 2).

The disk COMPRESSES large numbers toward the boundary:
  n=2 is at x=1/3 (well inside)
  n=100 is at x=99/101 = 0.98 (near boundary)
  n=10000 is at x=9999/10001 = 0.9998 (very near boundary)

But in the HYPERBOLIC METRIC, all consecutive integers
are separated by a distance of about 1/n (harmonic spacing).
""")

# ========================================================================
# Part 2: The map between them
# ========================================================================

print("=" * 70)
print("Part 2: The map arctanh: D -> C (and its inverse tanh: C -> D)")
print("=" * 70)

print("""
The Poincare disk D and the complex plane C are connected by:

  arctanh: D -> C     (disk -> plane)
  tanh: C -> D        (plane -> disk)

These are INVERSES: tanh(arctanh(w)) = w for all w in D.

What does arctanh DO geometrically?

On the REAL DIAMETER of D:
  arctanh maps (-1, 1) to (-inf, inf) = the real axis of C.
  It STRETCHES the bounded interval to the whole real line.
  Points near the boundary |x| -> 1 get sent to |arctanh(x)| -> inf.

On the IMAGINARY DIAMETER of D:
  arctanh(iy) = i*arctan(y).
  It maps (-i, i) to (-i*pi/4, i*pi/4) = a bounded imaginary interval!
  The IMAGINARY direction stays BOUNDED (unlike the real direction).

On the UNIT CIRCLE (boundary of D):
  arctanh(e^{i*theta}) = (1/2)*ln(cot(theta/2)) + i*pi/4.
  The imaginary part is CONSTANT at pi/4.
  The unit circle becomes a HORIZONTAL LINE at height pi/4.

On the INTERIOR of D:
  arctanh maps the disk to the HORIZONTAL STRIP:
    { z in C : -pi/4 < Im(z) < pi/4 }
  (Actually to a strip of width pi/2 centered on the real axis.)

So arctanh takes the DISK and UNFOLDS it into a STRIP:

  D (disk)  --arctanh-->  S = {z : |Im(z)| < pi/4}  (strip)

The real diameter of D becomes the real axis of S.
The imaginary diameter of D becomes the interval (-i*pi/4, i*pi/4) of S.
The boundary circle of D becomes the lines Im(z) = +/- pi/4 of S.
""")

# Verify the strip picture
print("Verification: arctanh maps D to the strip |Im(z)| < pi/4:")
print()
for r in [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]:
    for theta_deg in [0, 30, 45, 60, 90, 120, 150, 180]:
        theta = theta_deg * pi / 180
        w = r * cexp(1j * theta)
        z = catanh(w)
        if theta_deg in [0, 45, 90]:
            print(f"  r={r:.2f}, theta={theta_deg:3d} deg: arctanh = {z.real:8.4f} + {z.imag:8.4f}i, "
                  f"|Im| = {abs(z.imag):.4f}, pi/4 = {pi/4:.4f}, inside: {abs(z.imag) < pi/4 + 1e-10}")

print(f"""

Confirmed: for ALL interior points of D, |Im(arctanh(w))| < pi/4.

The strip width pi/4 = 0.7854 corresponds to 45 degrees = pi/4 radians.
This is ONE-EIGHTH of a full turn. The "Bott fraction" 1/8 appears again!

(The full strip width is pi/2 = 2 * pi/4 = 90 degrees.
 This is the range from -45 to +45 degrees.)


Part 3: THE THREE MODELS
==========================

There are THREE standard models of the hyperbolic plane:

MODEL 1: POINCARE DISK D = {{w : |w| < 1}}
  Metric: ds^2 = 4|dw|^2 / (1-|w|^2)^2
  Geodesics: circular arcs perpendicular to the boundary
  Advantage: conformal (angles preserved), compact
  Our coupling x lives here.

MODEL 2: UPPER HALF-PLANE H = {{z : Im(z) > 0}}
  Metric: ds^2 = |dz|^2 / Im(z)^2
  Geodesics: vertical lines and semicircles on the real axis
  Advantage: linear fractional transformations act cleanly
  Our Cayley transform Q(x) = (1+x)/(1-x) maps D to this (after rotation).

MODEL 3: THE STRIP S = {{z : |Im(z)| < pi/4}}
  Metric: ds^2 = |dz|^2 / cos^2(Im(z))   (derived from disk via arctanh)
  Geodesics: horizontal lines and certain curves
  Advantage: the group law is ADDITION (z1 + z2)
  Our rapidity w = arctanh(x) lives here.

The three models are connected by:

  D --arctanh--> S    (disk to strip: the logarithmic chart)
  D --Q--> H          (disk to half-plane: the Cayley chart)
  S --exp(2*-)--> H   (strip to half-plane: exponentiation)

Each model reveals DIFFERENT STRUCTURE:

  Model 1 (Disk): the BOUNDEDNESS of couplings, the unit circle as boundary.
  Model 2 (Half-plane): the MULTIPLICATIVE structure, SL(2,R) action.
  Model 3 (Strip): the ADDITIVE structure, simple addition.

Tournament theory uses:
  Model 1 for the coupling x (bounded comparison).
  Model 2 for the odds Q(x) = n (natural number).
  Model 3 for the rapidity arctanh(x) (unbounded magnitude).

ALL THREE ARE THE SAME SPACE. The hyperbolic plane.
""")

# ========================================================================
# Part 4: What "straight" means in each model
# ========================================================================

print("=" * 70)
print("Part 4: Geodesics — what 'straight' means")
print("=" * 70)

print("""
A GEODESIC is the shortest path between two points.
In each model, geodesics look different:

MODEL 1 (Disk): geodesics are ARCS OF CIRCLES that meet
the boundary at right angles. The real diameter and
imaginary diameter are geodesics (they're diameters = degenerate circles).

MODEL 3 (Strip): geodesics are HORIZONTAL LINES (Im(z) = const)
for the "straight" geodesics, and certain transcendental curves
for the others.

In the strip model, the REAL AXIS (Im(z) = 0) is the
"equator" — the central geodesic. This corresponds to:
  - The real diameter in the disk.
  - The positive real axis in the half-plane.
  - The space of REAL couplings in tournament theory.

Any point z = a + bi in the strip has:
  a = rapidity (position along the equator)
  b = latitude (perpendicular distance from equator)

The couplings x in (-1,1) have b = 0: they're ON the equator.
Complex couplings have b != 0: they're OFF the equator.

TOURNAMENT THEORY LIVES ON THE EQUATOR OF THE HYPERBOLIC PLANE.
""")

# ========================================================================
# Part 5: The group of isometries
# ========================================================================

print("=" * 70)
print("Part 5: Isometries — the symmetries of the hyperbolic plane")
print("=" * 70)

print("""
The ISOMETRIES (distance-preserving transformations) of the hyperbolic plane
form the group PSL(2,R) = SL(2,R) / {+/- I}.

In the disk model, every isometry is a MOBIUS TRANSFORMATION
that maps D to itself:

  g(w) = (a*w + b) / (conj(b)*w + conj(a))

where |a|^2 - |b|^2 = 1.

The SIMPLEST isometries are:

1. HYPERBOLIC TRANSLATION along the real diameter:
   g(w) = (w + t) / (1 + t*w)   where t in (-1,1).
   This IS the formal group law F(w, t).
   It slides points along the equator.

2. ROTATION around the center:
   g(w) = e^{i*theta} * w.
   This rotates the entire disk by angle theta.

3. HYPERBOLIC TRANSLATION along other geodesics:
   Conjugate of type 1 by a rotation.

Every isometry is a COMPOSITION of translations and rotations.

In the strip model:
   Hyperbolic translation = HORIZONTAL SHIFT z -> z + c (c real).
   Rotation = VERTICAL SHIFT z -> z + ci (c imaginary, |c| < pi/4).

So in the strip, ALL isometries connected to the identity
are just TRANSLATIONS z -> z + c for complex c with |Im(c)| < pi/4.

THE STRIP MODEL MAKES ISOMETRIES TRIVIAL: THEY'RE JUST ADDITION.

This is why arctanh is so powerful:
  In the disk, isometries are complicated Mobius transformations.
  In the strip, isometries are just z -> z + c.
  arctanh converts the complicated to the trivial.
""")

# ========================================================================
# Part 6: What the Cayley transform REALLY does
# ========================================================================

print("=" * 70)
print("Part 6: The Cayley transform — from disk to half-plane")
print("=" * 70)

print("""
The Cayley transform Q(w) = (1+w)/(1-w) maps:
  D (Poincare disk)  -->  H (upper half-plane, rotated)

Actually Q maps D to the RIGHT half-plane {Re(z) > 0}.
But for the upper half-plane, use Q_H(w) = i*(1+w)/(1-w).

On the real diameter:
  Q(x) = (1+x)/(1-x) maps (-1,1) to (0,inf).
  This IS the Cayley transform we've been using.
  x = 1/3 maps to Q(1/3) = 2 (the natural number 2).
  x = 3/5 maps to Q(3/5) = 4.
  x = (n-1)/(n+1) maps to Q(...) = n.

On the imaginary diameter:
  Q(iy) = (1+iy)/(1-iy) = (1-y^2+2iy)/(1+y^2).
  |Q(iy)| = 1 (it's on the unit circle!).
  arg(Q(iy)) = 2*arctan(y).

So Q maps:
  Real diameter (-1,1) --> positive real axis (0,inf): MAGNITUDES.
  Imaginary diameter (-i,i) --> unit circle S^1: PHASES.
  Full disk D --> right half-plane {Re > 0}: ALL.

THE CAYLEY TRANSFORM SEPARATES MAGNITUDE FROM PHASE:

  Re(Q(w)) = magnitude component
  Im(Q(w)) = phase component

  For w on the real axis: Q(w) is real = pure magnitude.
  For w on the imaginary axis: Q(w) is on S^1 = pure phase.
  For general w: Q(w) has both.

This is EXACTLY what happens in quantum mechanics:
  A state psi = r * e^{i*theta} has magnitude r and phase theta.
  The Cayley transform separates them.

And in tournament theory:
  The coupling x (real, on the equator) has pure magnitude Q(x) = n.
  An imaginary coupling iy (on the meridian) has pure phase Q(iy) = e^{2i*arctan(y)}.
  A complex coupling w has both.
""")

# Verify
print("Verification: Q on real and imaginary axes:")
print()
print("Real axis (magnitude):")
for x in [0, 1/3, 1/2, 3/5, 2/3, 3/4]:
    Q = (1 + x) / (1 - x)
    print(f"  Q({x:.4f}) = {Q:.4f} (= {Q:.4f}, magnitude)")

print("\nImaginary axis (phase):")
for y in [0, 0.2, 0.5, 1/sqrt(3), 0.8, 1.0]:
    w = 1j * y
    Q = (1 + w) / (1 - w)
    angle = atan2(Q.imag, Q.real) * 180 / pi
    print(f"  Q({y:.4f}i) = {Q.real:.4f} + {Q.imag:.4f}i, "
          f"|Q| = {abs(Q):.4f}, angle = {angle:.1f} deg")

# ========================================================================
# Part 7: The FUNDAMENTAL DIAGRAM
# ========================================================================

print(f"""

{"=" * 70}
Part 7: THE FUNDAMENTAL DIAGRAM
{"=" * 70}

                          arctanh
        POINCARE DISK D ============> STRIP S = C (mod periods)
             |                              |
             |                              |
          Q  |                              | exp(2*)
             |                              |
             v                              v
        RIGHT HALF-PLANE H <==========  C (via exp)


THREE SPACES, THREE MAPS, ONE OBJECT.

  D = Poincare disk = {{w : |w| < 1}}
    The bounded world. Couplings live here.
    Group law: F(w1,w2) = (w1+w2)/(1+w1w2).
    Metric: hyperbolic (constant negative curvature).

  S = Strip = {{z : |Im(z)| < pi/4}} (subset of C)
    The additive world. Rapidities live here.
    Group law: z1 + z2 (addition!).
    Metric: flat in the real direction, bounded in the imaginary.

  H = Right half-plane = {{z : Re(z) > 0}}
    The multiplicative world. Odds/natural numbers live here.
    Group law: z1 * z2 (multiplication, on the positive reals).
    Metric: hyperbolic (ds^2 = |dz|^2/Re(z)^2).

THE MAPS:

  arctanh: D -> S
    "Open up the disk into a strip."
    Makes the group law into addition.
    Converts bounded couplings to unbounded rapidities.

  Q = (1+w)/(1-w): D -> H
    "Unfold the disk into a half-plane."
    Makes real couplings into natural numbers.
    Separates magnitude from phase.

  exp(2z): S -> H
    "Exponentiate the strip into a half-plane."
    Q = exp(2*arctanh): the composition of the other two.

TOURNAMENT THEORY lives on the REAL AXIS of each space:
  D: x in (-1,1)          (couplings)
  S: w in (-inf,inf)       (rapidities)
  H: n in (0,inf)          (odds / natural numbers)

ROTATION THEORY lives on the IMAGINARY AXIS of each space:
  D: iy, |y| < 1           (imaginary couplings)
  S: i*theta, |theta| < pi/4   (angles)
  H: e^{{2i*theta}} on S^1       (phases)

THE BOUNDARY (unit circle of D = imaginary axis of H = lines Im = +/-pi/4 of S)
is the SPEED OF LIGHT / PHASE TRANSITION / UNREACHABLE LIMIT.


{"=" * 70}
Part 8: WHAT IT MEANS
{"=" * 70}

The Poincare disk is the NATURAL HOME of tournament couplings because:

1. BOUNDEDNESS: Couplings are bounded by +/-1. The disk is bounded by the unit circle.
   This is not a convention — it's the GEOMETRY of comparison.
   You cannot compare more than totally.

2. HYPERBOLIC GEOMETRY: The disk has negative curvature.
   This means: there is MORE ROOM near the boundary than near the center.
   In tournament theory: there are more ways to be strongly coupled
   than weakly coupled (exponentially more tournaments near x=1 than x=0).
   The DENSITY of tournaments increases toward the boundary.

3. THE GROUP LAW: F(x,y) = (x+y)/(1+xy) is the ISOMETRY GROUP of the disk.
   Combining two couplings is a hyperbolic translation.
   This is not a coincidence — it's because comparison is a GEOMETRIC operation.
   Composing two comparisons is moving in the hyperbolic plane.

4. THE BOUNDARY IS UNREACHABLE: |w| = 1 is at infinite hyperbolic distance.
   In tournament theory: coupling 1 (perfect agreement) requires infinite evidence.
   You can approach it (x -> 1) but never reach it.
   This is the SECOND LAW of tournament thermodynamics:
   you cannot reach complete certainty in finite steps.

The complex plane C, by contrast, is the COMPUTATIONAL HOME:
  Everything is linear (addition, scalar multiplication).
  Equations are easy to solve.
  But the geometry is FLAT — no curvature, no boundary, no speed limit.

arctanh TRANSLATES between the geometric truth (disk)
and the computational convenience (plane).

This is why we SOLVE in C but INTERPRET in D:
  The disk tells us WHAT IS TRUE (boundedness, curvature, density).
  The plane tells us HOW TO COMPUTE (linearity, roots, Fourier).

And the Cayley transform Q bridges to a THIRD world:
  The half-plane tells us WHAT THINGS ARE (natural numbers, primes, odds).

UNDERSTANDING THE DISK AND THE PLANE means understanding
that they are the SAME OBJECT seen with different eyes:
  The disk sees GEOMETRY (curvature, distance, boundary).
  The plane sees ALGEBRA (addition, multiplication, roots).
  The half-plane sees ARITHMETIC (natural numbers, primes, factoring).

And arctanh/tanh are the TRANSLATIONS between geometric truth
and algebraic convenience. They lose NOTHING — the isomorphism
is exact. But they change WHAT IS EASY TO SEE.

In the disk: boundedness is obvious, linearity is hidden.
In the plane: linearity is obvious, boundedness is hidden.

The art is knowing which view to use for which question.

  "How far apart are two couplings?" --> Use the DISK.
  "What is the m-th root of a coupling?" --> Use the PLANE.
  "What natural number does a coupling represent?" --> Use the HALF-PLANE.

And the universal zone (THM-227) is visible in ALL THREE:
  DISK: the first k Mersenne points on the real diameter.
  PLANE: the first k integers on the real axis.
  HALF-PLANE: the first k powers of 2 on the positive real axis.

The Mersenne number 2^k - 1 is the LAST POINT where all three views agree
before the k-nacci system's finiteness breaks the correspondence.
""")

# ========================================================================
# Part 9: The picture
# ========================================================================

print("=" * 70)
print("Part 9: The picture — seeing it all at once")
print("=" * 70)

print("""
Draw the Poincare disk. The real diameter runs left to right.

On the real diameter, mark the Cayley addresses:

  x_1 = 0      (center of the disk: n=1, the trivial tournament)
  x_2 = 1/3    (about 1/3 of the way to the right edge: n=2)
  x_3 = 1/2    (halfway: n=3)
  x_4 = 3/5    (n=4)
  x_5 = 2/3    (n=5)
  ...
  x_n = (n-1)/(n+1) --> 1 as n --> inf

These points CLUSTER toward the right boundary.
In Euclidean terms, they get closer and closer.
In HYPERBOLIC terms, they stay at distance ln((n+1)/n) apart
(decreasing, but never zero).

Now draw the imaginary diameter running top to bottom.

On the imaginary diameter, mark the "phase points":

  0            (center)
  i*tan(pi/6) = i/sqrt(3) ~ 0.577i     (30 degrees)
  i*tan(pi/4) = i                       (45 degrees = boundary!)
  i*tan(pi/3) = i*sqrt(3) ~ 1.732i     (60 degrees = outside!)

Wait: i*tan(pi/4) = i*1 = i, which is ON the boundary.
So the imaginary diameter (-i, i) maps via arctanh to
(-i*pi/4, i*pi/4), which is a BOUNDED range of angles.

The disk can only "see" angles up to 45 degrees directly.
To see larger angles, you need to go OUTSIDE the disk
(analytic continuation).

In tournament theory: the "natural" angles visible from
the real coupling are those within 45 degrees of the equator.
Beyond that lies the analytic continuation —
complex couplings with |w| > 1.

The 45-degree limit is pi/4 = 1/8 of a full turn.
THE BOTT PERIOD 8 appears: the strip width is pi/4 = (1/8) * 2*pi.

So the Poincare disk "sees" exactly 1/8 of the full angular range.
The remaining 7/8 requires analytic continuation beyond the disk.
The Bott period 8 is the ratio of the full circle to the visible angle.


THE DEEPEST PICTURE:

  The Poincare disk is a WINDOW into the hyperbolic plane.
  Through this window, we see:
    - The real diameter: tournament couplings (the equator).
    - The imaginary diameter: phases up to 45 degrees.
    - The boundary circle: the speed of light.
    - The interior: all "subluminal" combinations.

  arctanh OPENS the window: unrolls the disk into the strip.
  tanh CLOSES the window: rolls the strip back into the disk.

  The complex plane C is the FULL view — no window, no boundary.
  But it has no curvature — no sense of "how far" in the geometric sense.

  To UNDERSTAND both is to understand that:
    The disk is the TRUTH (bounded, curved, with a horizon).
    The plane is the TOOL (flat, unbounded, algebraically tractable).
    And they are the SAME SPACE, because arctanh is an isomorphism.

  The "same space" part is the deepest point.
  There is no information lost in going from disk to plane.
  The curvature of the disk BECOMES the nonlinearity of tanh.
  The flatness of the plane BECOMES the boundedness of the disk.
  Different words for the same thing.

  And the Cayley transform Q adds a third word:
    The natural numbers ARE the hyperbolic plane.
    Each natural number n sits at a specific point in the disk.
    Their SPACING is the hyperbolic metric.
    Their MULTIPLICATION is the group law.
    Their PRIMALITY is a geometric property.

  The Poincare disk, the complex plane, and the natural numbers
  are three languages for one mathematical object.
  tanh, arctanh, and Q are the dictionaries.
""")

print("opus-2026-03-16-S90ax: UNDERSTAND THE POINCARE DISK AND THE COMPLEX PLANE")
