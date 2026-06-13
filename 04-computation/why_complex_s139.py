#!/usr/bin/env python3
"""
why_complex_s139.py — Why Complex Numbers Solve Equations
opus-2026-03-21-S139

The user asks: WHY are complex numbers needed? What are the rotations?
What's happening with shapes?

These are foundational questions. Let me answer them honestly.
"""

from math import comb, sqrt, pi, cos, sin, atan2
import cmath

print("=" * 72)
print("  WHY COMPLEX NUMBERS SOLVE EQUATIONS")
print("  opus-2026-03-21-S139")
print("=" * 72)
print()

# ==========================================================================
# PART 1: THE REAL LINE ISN'T ENOUGH
# ==========================================================================

print("PART 1: THE REAL LINE ISN'T ENOUGH")
print("-" * 60)
print()

print("""  Consider x² = -1.

  On the real line, x² is always ≥ 0.
  So x² = -1 has no solution.

  But this equation SHOULD have a solution. It's degree 2,
  and degree-2 equations "should" have 2 roots.

  WHY should it? Because of SHAPES.

  Plot y = x² + 1. It's a parabola sitting entirely above the x-axis.
  It never crosses zero. No real root exists.

  But now think about it differently. Don't plot y = f(x) on a line.
  Instead, think of the FUNCTION f(z) = z² + 1 acting on the PLANE.

  Every point z in the plane gets mapped to f(z) = z² + 1.
  The question "where does f = 0?" becomes:
  "which points in the plane map to the origin?"

  And the answer is z = i and z = -i, two points on the
  vertical axis, at distance 1 above and below the real line.

  THE REAL LINE IS A 1D SLICE OF THE 2D STORY.
  The parabola doesn't cross zero because zero isn't on the
  1D slice. But it IS in the 2D plane.
""")

# ==========================================================================
# PART 2: WHAT THE PLANE ADDS — ROTATION
# ==========================================================================

print()
print("PART 2: WHAT THE PLANE ADDS — ROTATION")
print("-" * 60)
print()

print("""  On the real line, multiplication does two things:
    SCALE (make bigger or smaller)
    FLIP (change sign, ×(-1))

  That's all. Scale and flip. One dimension of operation.

  On the complex plane, multiplication does three things:
    SCALE (change distance from origin)
    FLIP (still available: ×(-1) = rotate 180°)
    ROTATE (turn by ANY angle)

  Multiplication by z = r·e^{iθ} = r(cos θ + i sin θ):
    SCALES by r
    ROTATES by θ

  The real line only sees θ = 0 (positive) and θ = π (negative).
  The complex plane sees ALL angles θ ∈ [0, 2π).

  THIS IS WHY EQUATIONS NEED THE PLANE:
  The roots of z^n = 1 are at angles 2πk/n for k = 0, 1, ..., n-1.
  They form a REGULAR n-GON.

  On the real line, z^n = 1 has at most 2 solutions (±1).
  On the complex plane, z^n = 1 has EXACTLY n solutions,
  equally spaced around the unit circle.
""")

# Show the n-th roots of unity as regular polygons
print("  ROOTS OF UNITY = REGULAR POLYGONS:")
print()
for n in range(2, 8):
    roots = [cmath.exp(2j * cmath.pi * k / n) for k in range(n)]
    shape_name = {2: "segment", 3: "TRIANGLE", 4: "square",
                  5: "pentagon", 6: "hexagon", 7: "heptagon"}[n]
    root_strs = []
    for r in roots:
        if abs(r.imag) < 1e-10:
            root_strs.append(f"{r.real:+.3f}")
        else:
            root_strs.append(f"{r.real:+.3f}{r.imag:+.3f}i")
    print(f"  z^{n} = 1:  {shape_name:9s}  roots: {', '.join(root_strs)}")

print()

# ==========================================================================
# PART 3: THE FUNDAMENTAL THEOREM — TOPOLOGY, NOT ALGEBRA
# ==========================================================================

print()
print("PART 3: WHY EVERY POLYNOMIAL HAS ROOTS — IT'S TOPOLOGY")
print("-" * 60)
print()

print("""  The fundamental theorem of algebra says:
  Every polynomial of degree n has exactly n roots in C.

  But WHY? The proof is TOPOLOGICAL, not algebraic:

  1. Take a polynomial p(z) = z^n + lower terms.
  2. Draw a big circle |z| = R in the complex plane.
  3. As z goes around this circle, p(z) ALSO traces a curve.
  4. For large R, p(z) ≈ z^n, so the image curve winds
     around the origin EXACTLY n TIMES.
  5. As you shrink R → 0, the image curve shrinks to p(0).
  6. The winding number can only change when the curve
     CROSSES the origin, i.e., when p(z) = 0 for some z on the circle.
  7. The winding number starts at n (large R) and must reach 0 (R=0).
  8. So p(z) = 0 for at least one z.
  9. Factor out (z - root), reduce degree by 1, repeat.

  THE DEEP REASON: the complex plane has the right TOPOLOGY.
  It's 2-dimensional, simply connected, and multiplication
  acts by rotation. A function that winds n times around the
  origin MUST pass through the origin at least once.

  The real line doesn't have this property because it's 1D.
  A function can "jump over" zero on a line.
  It can't "jump over" zero in a plane — there's nowhere to jump.

  SHAPES: the polynomial p(z) acts on the plane like a
  RUBBER SHEET that wraps around the origin n times.
  The points where the sheet touches the origin are the roots.
  There MUST be at least n such points because of the wrapping.
""")

# ==========================================================================
# PART 4: THE SHAPES — WHAT ROOTS LOOK LIKE
# ==========================================================================

print()
print("PART 4: THE SHAPES — ROOT GEOMETRY")
print("-" * 60)
print()

# The tribonacci equation: x³ - x² - x - 1 = 0
# This is the gap function g(x) = 0
# It has one real root (τ) and two complex conjugate roots

coeffs = [1, -1, -1, -1]  # x³ - x² - x - 1
# Roots: use the cubic formula or numerical methods
import numpy as np
roots = np.roots(coeffs)

print(f"  THE GAP FUNCTION g(x) = x³ - x² - x - 1 = 0:")
print()
for i, r in enumerate(roots):
    if abs(r.imag) < 1e-10:
        print(f"    root {i+1}: {r.real:.10f}  (REAL = τ, the tribonacci constant)")
    else:
        mod = abs(r)
        arg = cmath.phase(r)
        print(f"    root {i+1}: {r.real:.6f} {r.imag:+.6f}i  "
              f"(|z|={mod:.6f}, arg={arg:.6f} rad = {arg*180/pi:.2f}°)")

print()

# The real root is τ. The complex roots are at |z| and angle ±θ.
tau = max(r.real for r in roots if abs(r.imag) < 1e-10)
complex_root = [r for r in roots if abs(r.imag) > 0.01][0]
mod_c = abs(complex_root)
arg_c = abs(cmath.phase(complex_root))

print(f"  THE THREE ROOTS FORM A SHAPE:")
print(f"    τ = {tau:.6f} (real, on the positive x-axis)")
print(f"    Two complex roots at distance {mod_c:.6f} from origin,")
print(f"    at angles ±{arg_c*180/pi:.2f}° from the real axis.")
print()
print(f"    The three roots form a (distorted) TRIANGLE in the complex plane!")
print(f"    It's not equilateral (the real root is farther from origin)")
print(f"    but it IS a triangle — the 3-fold structure is geometric.")
print()

# Product of roots = constant term (up to sign) = 1
# Sum of roots = coefficient of x² = 1
print(f"  BY VIETA'S FORMULAS:")
print(f"    Sum of roots = {sum(roots):.6f}  (should be 1)")
print(f"    Product of roots = {np.prod(roots):.6f}  (should be 1)")
print(f"    Sum of pairwise products = {sum(roots[i]*roots[j] for i in range(3) for j in range(i+1,3)):.6f}  (should be -1)")
print()

# The complex roots carry the OSCILLATORY behavior
# When you evaluate x³ - x² - x - 1 at real x:
#   x < τ: dominated by the complex roots (oscillation)
#   x > τ: dominated by the real root (growth)
# The tribonacci SEQUENCE T(n) satisfies T(n) = τ^n/... + oscillatory terms

print(f"  THE COMPLEX ROOTS CARRY THE OSCILLATION:")
print(f"    The tribonacci sequence T(n) ~ A·τⁿ + B·z₂ⁿ + C·z₃ⁿ")
print(f"    where z₂, z₃ are the complex roots.")
print(f"    |z₂| = |z₃| = {mod_c:.6f} < 1, so these terms DECAY.")
print(f"    The oscillatory part dies out, leaving pure τⁿ growth.")
print()
print(f"    At x = 2 (tournament fugacity):")
print(f"    g(2) = 8 - 4 - 2 - 1 = 1")
print(f"    We're BEYOND τ = {tau:.6f}, in the growth regime.")
print(f"    The complex roots are invisible at x = 2:")
print(f"    they've already decayed away.")
print()

# ==========================================================================
# PART 5: WHY THREE — THE TRIANGLE AND THE CUBIC
# ==========================================================================

print()
print("PART 5: WHY THREE — THE CUBIC AND THE TRIANGLE")
print("-" * 60)
print()

print("""  The simplest equation that NEEDS complex numbers is x² = -1.
  But the most INTERESTING is the cubic: x³ + px + q = 0.

  HISTORICAL FACT (Cardano, 1545):
  The cubic formula SOMETIMES requires complex numbers
  even when ALL THREE ROOTS ARE REAL.

  This is called "casus irreducibilis" — the irreducible case.
  When the discriminant is negative (three real roots),
  Cardano's formula produces expressions like ∛(a + bi),
  which are complex even though the final answer is real.

  The complex numbers are NEEDED as intermediate steps.
  You must PASS THROUGH the complex plane to reach real roots.

  WHY? Because the three real roots form a triangle
  (when viewed as roots of unity of a related equation),
  and you can't describe a triangle using only the real line.
  You need the plane.

  THE TRIANGLE IS THE FIRST SHAPE THAT REQUIRES THE PLANE.
    - 2 points: a line segment (1D suffices)
    - 3 points: a triangle (NEEDS 2D)
    - 4 points: a tetrahedron (needs 3D... or does it?)

  Equations of degree 2 can be solved with real square roots.
  Equations of degree 3 FORCE you into the complex plane.
  The 3-cycle in tournament theory appears at n=3.
  This is the SAME transition.
""")

# Demonstrate casus irreducibilis
# x³ - 7x + 6 = 0 has roots 1, 2, -3 (all real)
# But Cardano's formula gives complex intermediate values
p_cubic = -7
q_cubic = 6
disc = -(4*p_cubic**3 + 27*q_cubic**2)
print(f"  EXAMPLE: x³ - 7x + 6 = 0")
print(f"  Roots: 1, 2, -3 (ALL real)")
print(f"  Discriminant: {disc} > 0 (three real roots)")
print(f"  But Cardano's formula requires ∛(complex)!")
print()

# Cardano: x = ∛(-q/2 + √(q²/4 + p³/27)) + ∛(-q/2 - √(q²/4 + p³/27))
inner = q_cubic**2/4 + p_cubic**3/27
print(f"  Inner discriminant q²/4 + p³/27 = {inner:.4f} < 0!")
print(f"  √({inner:.4f}) = {cmath.sqrt(inner)}")
print(f"  The formula REQUIRES passing through complex numbers")
print(f"  even though the answer is real.")
print()

# ==========================================================================
# PART 6: THE SHAPE CONNECTION — POLYGONS AND POLYNOMIALS
# ==========================================================================

print()
print("PART 6: POLYGONS AND POLYNOMIALS")
print("-" * 60)
print()

print("""  Every polynomial of degree n defines a SHAPE:
  its n roots form a POLYGON in the complex plane.

  z^n - 1 = 0:  regular n-gon (most symmetric)
  z^n + 1 = 0:  regular n-gon rotated by π/n
  General:      distorted n-gon

  THE SHAPE ENCODES THE POLYNOMIAL.
  By Vieta's formulas:
    sum of roots = -a_{n-1}/a_n          (centroid)
    product of roots = (-1)^n a_0/a_n    (signed area, roughly)
    sum of pairwise products = a_{n-2}/a_n (second moment)

  The polynomial IS its root polygon.
  The coefficients ARE the polygon's moments.

  NOW THE TOURNAMENT CONNECTION BECOMES CLEAR:

  The n-th roots of unity form a regular n-gon.
  The THIRD roots of unity form an equilateral TRIANGLE.
  The third roots are 1, ω, ω² where ω = e^{2πi/3}.

  The QUATERNION imaginary units i, j, k satisfy:
    ij = k (cyclic product)
  Their products trace the same 3-cycle as the cube roots of unity.
  The 3-fold symmetry is the SAME 3-fold symmetry.

  The TOURNAMENT 3-CYCLE a → b → c → a is the SAME shape:
  three elements arranged cyclically, with the direction
  encoding the "sign" (the ±1 of the imaginary unit product).
""")

# Cube roots of unity
omega = cmath.exp(2j * pi / 3)
print(f"  CUBE ROOTS OF UNITY:")
print(f"    1   = {1.0:+.4f} + {0.0:+.4f}i  (angle 0°)")
print(f"    ω   = {omega.real:+.4f} + {omega.imag:+.4f}i  (angle 120°)")
print(f"    ω²  = {(omega**2).real:+.4f} + {(omega**2).imag:+.4f}i  (angle 240°)")
print(f"    These form an equilateral TRIANGLE on the unit circle.")
print()
print(f"  Key property: 1 + ω + ω² = {(1 + omega + omega**2).real:.6f}")
print(f"  The roots SUM TO ZERO. The triangle is centered at the origin.")
print()

# ==========================================================================
# PART 7: THE DEEP ANSWER — WHY THE PLANE, WHY ROTATION
# ==========================================================================

print()
print("PART 7: THE DEEP ANSWER")
print("-" * 60)
print()

print("""  WHY are complex numbers needed to solve equations?

  Because MULTIPLICATION IS ROTATION, and you need the PLANE
  to have enough room for all the rotations.

  On the real line:
    Multiplication by +1: identity (0° rotation)
    Multiplication by -1: flip (180° rotation)
    That's it. Two rotations. Not enough to solve x^n = 1 for n > 2.

  On the complex plane:
    Multiplication by e^{iθ}: rotation by ANY angle θ
    Now x^n = 1 has n solutions: rotations by 360°/n, 2·360°/n, ..., n·360°/n
    These are the vertices of a regular n-gon.

  WHY rotations? Because multiplication SHOULD be:
    (1) continuous (small changes in input → small changes in output)
    (2) invertible (every nonzero number has a reciprocal)
    (3) commutative (ab = ba)
    (4) associative ((ab)c = a(bc))

  The ONLY 2D multiplication satisfying (1)-(4) is
  complex multiplication = rotation + scaling.
  This is a THEOREM (the classification of real division algebras).

  So the complex plane isn't a CHOICE. It's the UNIQUE
  2D system where multiplication works properly.

  WHY shapes? Because the roots of a polynomial ARE the
  vertices of a shape. The polynomial IS a description of
  the shape via its moments (coefficients = symmetric functions
  of roots). Solving the equation = FINDING the shape.

  THE PROGRESSION:
    Degree 1: 1 root = 1 point. No shape needed.
    Degree 2: 2 roots = 2 points = line segment. The real line suffices
              IF both roots are real. If not, you need the plane.
    Degree 3: 3 roots = triangle. You ALWAYS need the plane
              (casus irreducibilis: even all-real roots need complex
              intermediate steps). The triangle is INHERENTLY 2D.
    Degree n: n roots = n-gon. The plane is always sufficient.

  DEGREE 3 IS THE CRITICAL STEP.
  It's where the plane becomes NECESSARY, not just convenient.
  And degree 3 is where:
    - Tournament 3-cycles first appear (n=3)
    - The tribonacci equation lives (x³ - x² - x - 1 = 0)
    - The Cayley-Dickson tower passes through quaternions (3 imaginaries)
    - The {3,∞} tessellation has its face (triangle)

  THE TRIANGLE IS THE REASON FOR EVERYTHING:
  It's the first shape that can't live on a line.
  It forces you into the plane.
  The plane gives you rotation.
  Rotation gives you imaginary units.
  Imaginary units give you the Cayley-Dickson tower.
  The tower gives you octonions, Fano, forbidden values.

  It all starts with: THREE POINTS DON'T FIT ON A LINE.
""")

# ==========================================================================
# PART 8: THE TOURNAMENT VERSION
# ==========================================================================

print()
print("PART 8: THE TOURNAMENT VERSION")
print("-" * 60)
print()

print("""  In tournament theory, the same progression:

  n = 2: Two players. One beats the other. Total order.
         Lives on a LINE. No complex structure needed.
         H = 1 always. One shape: a segment.

  n = 3: Three players. A can beat B, B beats C, C beats A.
         This is a CYCLE. It can't live on a line (a line
         would require a total order: A > B > C or some permutation).
         The cycle IS the triangle. It forces you OFF the line.
         H ∈ {1, 3}. Two shapes: segment or triangle.

  n = 4: Four players. Multiple cycles can coexist.
         H ∈ {1, 3, 5}. The shapes become 3D.
         (But 3D is not needed yet — all H values are achievable.)

  n = 5: Five players. H ∈ {1,3,5,9,11,13,15}. Gap at H = 7.
         The gap IS the "casus irreducibilis": a value that
         SHOULD exist (it's a valid odd number in range)
         but CAN'T be reached without passing through
         impossible intermediate configurations.

  The forbidden H = 7 is to tournaments what
  "casus irreducibilis" is to cubics:
  THE POINT WHERE THE EXTRA DIMENSION BECOMES INESCAPABLE.

  In cubics: you can't avoid complex numbers even for real roots.
  In tournaments: you can't avoid the octonionic structure
  even for integer H values.

  The complex plane IS the extra dimension.
  The tournament gap IS the evidence that the extra dimension exists.
  The 3-cycle IS the triangle that forced us into the plane.
""")

# ==========================================================================
# PART 9: e^{iπ} + 1 = 0 AND H = 1 + 2|Ω|
# ==========================================================================

print()
print("PART 9: TWO EQUATIONS THAT SAY THE SAME THING")
print("-" * 60)
print()

print(f"""  EULER: e^{{iπ}} + 1 = 0

  This says: a half-turn rotation (e^{{iπ}} = -1) plus the identity (+1)
  equals nothing (0). The rotation CANCELS the identity.

  TOURNAMENT: H = 1 + 2|Ω|  (at n=5, where Ω is complete)

  This says: the path count = identity (1 transitive path) plus
  the cycle contribution (2 × number of directed odd cycles).

  Both equations decompose a quantity into:
    IDENTITY PART + ROTATION PART = TOTAL

  In Euler: the total is 0 (exact cancellation).
  In tournaments: the total is H (the Hamiltonian path count).

  The "+1" in both is the SAME "+1": the identity element,
  the vacuum, the transitive path, the ground state.

  The "rotation part" is:
    In Euler: e^{{iπ}} = -1 (a single half-turn, fully destructive)
    In tournaments: 2|Ω| (many quarter/third/... turns, constructive)

  Euler's equation works because ONE rotation exactly cancels.
  Tournaments work because MANY rotations partially add.

  The difference: Euler evaluates at the MOST SPECIAL point (θ = π).
  Tournaments evaluate at a GENERIC point (x = 2).
  Special points cancel; generic points accumulate.
""")

print()
print("opus-2026-03-21-S139")
