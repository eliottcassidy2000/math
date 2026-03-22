#!/usr/bin/env python3
"""
half_steps_s142.py — Half-Steps Orthogonal to {-1, 0, 1}
opus-2026-03-21-S142

The user says: understand half-steps. Understand them as orthogonal to {-1, 0, 1}.

Let me think about what this means.

{-1, 0, 1} is the INTEGER lattice restricted to the smallest values.
It classifies: negative, zero, positive. Three states.

The HALF-STEPS {-1/2, +1/2} sit BETWEEN the integers.
They don't classify — they INTERPOLATE.

"Orthogonal" means: they carry INDEPENDENT information.
{-1, 0, 1} tells you WHICH regime.
The half-step tells you something ELSE entirely.

What is that something else?
"""

from math import sqrt, pi, log, log2, exp, cos, sin
import cmath

print("=" * 72)
print("  HALF-STEPS ORTHOGONAL TO {-1, 0, 1}")
print("  opus-2026-03-21-S142")
print("=" * 72)
print()

# ==========================================================================
# PART 1: WHAT {-1, 0, 1} IS
# ==========================================================================

print("PART 1: WHAT {-1, 0, 1} IS")
print("-" * 60)
print()

print("""  {-1, 0, 1} appears everywhere:

  SIGNS:      negative, zero, positive
  COMPARISON: less than, equal, greater than
  GAP:        g(ρ_{p-1})=-1, g(ρ_p)=0, g(2)=+1
  CURVATURE:  spherical, flat, hyperbolic
  LEGENDRE:   (a/p) = -1, 0, +1 (quadratic residue symbol)
  MATRICES:   eigenvalues of a reflection (sign matrix)
  TOPOLOGY:   Euler char of S⁰ = {-1, +1}, point = {0}

  {-1, 0, 1} is the CLASSIFICATION. It answers: what TYPE?
  It's DISCRETE. It's the VERTICAL axis.

  {-1, 0, 1} = Z/3 ≅ the integers mod 3? No.
  {-1, 0, 1} = the sign function sgn(x) restricted to integers.
  It's the coarsest nontrivial partition of the number line.
""")

# ==========================================================================
# PART 2: WHAT THE HALF-STEP IS
# ==========================================================================

print()
print("PART 2: WHAT THE HALF-STEP IS")
print("-" * 60)
print()

print("""  The half-steps are {-1/2, +1/2}.

  They sit at the MIDPOINTS between consecutive integers:
    -1 ——— -1/2 ——— 0 ——— +1/2 ——— +1

  In the gap function: g_p(x) takes values {-1, 0, +1} at
  the special points {ρ_{p-1}, ρ_p, 2}. The HALF-VALUES
  g_p(x) = ±1/2 occur at INTERMEDIATE points.

  But the user says these are ORTHOGONAL to {-1, 0, 1}.
  Not between them. PERPENDICULAR to them.

  What does this mean?

  Think of the COMPLEX PLANE:
    The real axis carries: ..., -1, -1/2, 0, +1/2, +1, ...
    The imaginary axis carries: ..., -i, -i/2, 0, +i/2, +i, ...

  {-1, 0, +1} lives on the REAL axis.
  {-i/2, +i/2} lives on the IMAGINARY axis.

  But i/2 is not "half of i." It's "i divided by 2."
  The TRUE half-step in the imaginary direction is √i:
    √i = e^{iπ/4} = (1+i)/√2 ≈ 0.707 + 0.707i

  This is at 45° — exactly HALFWAY between real and imaginary.
  The half-step IS the diagonal. The 45° angle.
""")

# ==========================================================================
# PART 3: THE FOUR FOURTH-ROOTS OF UNITY
# ==========================================================================

print()
print("PART 3: THE FOUR FOURTH-ROOTS OF UNITY")
print("-" * 60)
print()

# Fourth roots of 1: {1, i, -1, -i}
# These form a SQUARE on the unit circle.
# The HALF-STEPS between {-1, +1} are {i, -i}.
# The imaginary units ARE the half-steps of the real classification.

fourth_roots = [1, 1j, -1, -1j]
print(f"  The FOURTH ROOTS OF UNITY: {{1, i, -1, -i}}")
print(f"  They form a SQUARE on the unit circle:")
print()
for k, z in enumerate(fourth_roots):
    angle = k * 90
    print(f"    z_{k} = {z.real:+.0f}{z.imag:+.0f}i  (angle {angle}°)")
print()
print(f"  {'{'}1, -1{'}'} = the WHOLE STEPS (real axis)")
print(f"  {'{'}i, -i{'}'} = the HALF-STEPS (imaginary axis)")
print(f"  The half-steps are ORTHOGONAL to the whole steps.")
print(f"  They are PERPENDICULAR, not in-between.")
print()

# The key: 1 → i → -1 → -i → 1
# Moving by ONE half-step = rotating by 90°
# Moving by TWO half-steps = rotating by 180° = negation
# So: i² = -1 means "two half-steps = one whole step of negation"
print(f"  THE HALF-STEP SQUARED = THE WHOLE STEP:")
print(f"    i² = -1")
print(f"    Two 90° rotations = one 180° rotation")
print(f"    Two half-steps = one whole step")
print(f"    The half-step is the SQUARE ROOT of the whole step.")
print()

# ==========================================================================
# PART 4: {-1, 0, +1} AS THE SIGN, HALF-STEPS AS THE PHASE
# ==========================================================================

print()
print("PART 4: SIGN vs PHASE")
print("-" * 60)
print()

print("""  {-1, 0, +1} is the SIGN of a real number:
    sgn(x) = -1 if x < 0, 0 if x = 0, +1 if x > 0

  This is COARSE information. It tells you the DIRECTION
  but not the DISTANCE.

  The HALF-STEP adds PHASE information:
    A complex number z = r·e^{iθ} has:
      MAGNITUDE r (how far from 0)
      PHASE θ (which direction)

  The sign {-1, 0, +1} corresponds to phases {π, undefined, 0}:
    sgn(x) = +1: phase 0 (pointing right)
    sgn(x) = -1: phase π (pointing left)
    sgn(x) = 0:  no phase (at the origin)

  The HALF-STEPS correspond to phases {π/2, 3π/2}:
    +i: phase π/2 (pointing up)
    -i: phase 3π/2 (pointing down)

  So:
    {-1, 0, +1} = HORIZONTAL information (left/center/right)
    half-steps = VERTICAL information (up/down)

  They are ORTHOGONAL because horizontal ⊥ vertical.

  The FULL information needs BOTH:
    x = sgn(x) · |x| + i · (phase correction)
    Or: z = x + iy (real + imaginary)

  The gap function g_p gives the SIGN: {-1, 0, +1} = which regime.
  The PHASE (the complex roots) gives the HALF-STEP: the oscillatory
  content that the sign misses.
""")

# ==========================================================================
# PART 5: THE GAP FUNCTION'S COMPLEX ROOTS AS HALF-STEPS
# ==========================================================================

print()
print("PART 5: THE GAP FUNCTION'S COMPLEX ROOTS = HALF-STEPS")
print("-" * 60)
print()

# The gap function g_3(x) = x³ - x² - x - 1
# Real root: τ = 1.839 (the sign boundary, g_3(τ) = 0)
# Complex roots: two conjugate roots inside the unit circle

import numpy as np
coeffs = [1, -1, -1, -1]
roots = np.roots(coeffs)
tau = max(r.real for r in roots if abs(r.imag) < 1e-6)
complex_roots = [r for r in roots if abs(r.imag) > 0.01]

print(f"  g_3(x) = x³ - x² - x - 1 = 0:")
print(f"    Real root: τ = {tau:.8f}  (the SIGN boundary)")
for r in complex_roots:
    mod = abs(r)
    phase = cmath.phase(r)
    print(f"    Complex root: {r.real:+.6f} {r.imag:+.6f}i  "
          f"|z|={mod:.6f}, phase={phase*180/pi:+.1f}°")

print()

# The complex roots carry PHASE information that the real root doesn't.
# The real root τ tells you WHERE the sign changes (the boundary).
# The complex roots tell you HOW FAST the function oscillates
# and HOW QUICKLY the oscillation decays.

print(f"  The REAL ROOT τ gives the SIGN boundary:")
print(f"    g_3(x) < 0 for x < τ  (spherical, sign = -1)")
print(f"    g_3(x) = 0 at x = τ   (Euclidean, sign = 0)")
print(f"    g_3(x) > 0 for x > τ  (hyperbolic, sign = +1)")
print()

print(f"  The COMPLEX ROOTS give the PHASE (orthogonal to sign):")
print(f"    |z| = {abs(complex_roots[0]):.6f} < 1 (decaying)")
print(f"    phase = ±{abs(cmath.phase(complex_roots[0]))*180/pi:.1f}°")
print(f"    The decay rate |z| tells HOW FAST oscillations die.")
print(f"    The phase tells the ANGULAR FREQUENCY of oscillation.")
print()

# In the tribonacci sequence T(n):
# T(n) = A·τⁿ + B·z₂ⁿ + C·z₃ⁿ
# The real part A·τⁿ gives the GROWTH (the sign).
# The complex parts B·z₂ⁿ + C·z₃ⁿ give the OSCILLATION (the half-step).
# Since |z₂| < 1, the oscillation DECAYS.
# But it's ALWAYS THERE at finite n.

print(f"  In the tribonacci sequence T(n) = A·τⁿ + 2Re(B·z₂ⁿ):")
print(f"    A·τⁿ = the TREND (exponential growth, the sign)")
print(f"    2Re(B·z₂ⁿ) = the OSCILLATION (half-step, decaying)")
print()

# Compute the oscillation for first few n
# Tribonacci: T(0)=0, T(1)=0, T(2)=1, T(n) = T(n-1)+T(n-2)+T(n-3)
T = [0, 0, 1]
for n in range(3, 15):
    T.append(T[-1] + T[-2] + T[-3])

# The Binet-like formula: T(n) ≈ c·τⁿ where c = τ/(τ²+4τ+1) approximately
# Actually, using exact decomposition:
# We need A such that T(n) ≈ A·τⁿ for large n
# A = lim T(n)/τⁿ
A = T[-1] / tau**(len(T)-1)

print(f"  {'n':>3s}  {'T(n)':>8s}  {'A·τⁿ':>10s}  {'oscillation':>12s}  {'osc/trend':>10s}")
print(f"  {'─'*3}  {'─'*8}  {'─'*10}  {'─'*12}  {'─'*10}")
for n in range(2, 13):
    trend = A * tau**n
    osc = T[n] - trend
    ratio = abs(osc/trend) if abs(trend) > 0.01 else 0
    print(f"  {n:3d}  {T[n]:8d}  {trend:10.2f}  {osc:12.4f}  {ratio:10.6f}")

print()
print(f"  The oscillation DECAYS exponentially (ratio |z₂/τ|ⁿ → 0)")
print(f"  But at small n, it's SIGNIFICANT — up to ~30% of the trend!")
print(f"  The half-step content is LARGE at small n, SMALL at large n.")
print()

# ==========================================================================
# PART 6: THE ORTHOGONALITY — REAL AND IMAGINARY ARE INDEPENDENT
# ==========================================================================

print()
print("PART 6: THE ORTHOGONALITY")
print("-" * 60)
print()

print("""  The half-steps are orthogonal to {-1, 0, +1} because:

  1. SIGN AND PHASE ARE INDEPENDENT:
     Knowing the sign of g_p(x) tells you the regime.
     Knowing the phase of the complex roots tells you the oscillation.
     They are INDEPENDENT pieces of information.

  2. REAL AND IMAGINARY PARTS ARE ORTHOGONAL:
     In C: Re(z) ⊥ Im(z) (literally, perpendicular axes).
     The sign lives on the real axis.
     The half-step lives on the imaginary axis.

  3. THE INTEGER AND HALF-INTEGER LATTICES ARE DUAL:
     Z = {..., -1, 0, 1, 2, ...} (integer lattice)
     Z + 1/2 = {..., -1/2, 1/2, 3/2, ...} (half-integer lattice)
     In Fourier analysis: the DUAL of Z is the circle S¹.
     The half-integer lattice gives the ANTIPERIODIC modes.
     Periodic (Z) and antiperiodic (Z+1/2) are orthogonal.

  4. BOSONS AND FERMIONS:
     Bosons have integer spin: {0, 1, 2, ...}
     Fermions have half-integer spin: {1/2, 3/2, 5/2, ...}
     They obey DIFFERENT statistics (commuting vs anticommuting).
     This is the deepest physical manifestation of the orthogonality.

  THE CONNECTION TO TOURNAMENTS:

  {-1, 0, +1} = the gap function values = WHICH theory type
  Half-steps = the complex roots = HOW the theory oscillates

  The OCR captures the {-1, 0, +1} information (which regime).
  The RESIDUAL is the half-step information (the oscillation).
  They are orthogonal: the OCR and the residual are independent.
  Together they give the full picture.
""")

# ==========================================================================
# PART 7: HALF-STEPS IN MUSIC AND MATHEMATICS
# ==========================================================================

print()
print("PART 7: THE MUSICAL ANALOGY (IT'S NOT AN ANALOGY)")
print("-" * 60)
print()

print("""  In music, the chromatic scale has 12 half-steps per octave.
  An octave is a DOUBLING of frequency (ratio 2:1).
  Each half-step is a frequency ratio of 2^{1/12}.

  12 half-steps = 2^{12/12} = 2^1 = one doubling.

  Now: 12 = h(E₆) = h(F₄) = weight of Ramanujan Δ = φ(42)
  = number of iso classes at n=5.

  This is NOT a coincidence. Here's why:

  The CHROMATIC SCALE divides the DOUBLING into 12 equal parts.
  The TOURNAMENT CLASSIFICATION divides the iso classes into 12 types.
  Both use 12 because:

    12 = LCM(3, 4) = the smallest number divisible by both 3 and 4.

  3 = the tournament atom (3-cycle, the {3,∞} face).
  4 = the first non-trivial square (4 = 2², the first doubling).

  The chromatic scale needs 12 steps because:
    4 half-steps = major third (frequency ratio 2^{4/12} ≈ 5/4)
    3 half-steps = minor third (frequency ratio 2^{3/12} ≈ 6/5)
    The thirds subdivide the octave into 3 major or 4 minor thirds.
    LCM(3, 4) = 12 gives the finest scale that captures both.

  In tournament theory:
    3-cycles subdivide the tournament into triangular faces.
    4-vertex subtournaments give the first non-trivial structure.
    The 12 iso classes at n=5 capture both levels.

  THE HALF-STEP 2^{1/12} IS THE 12th ROOT OF THE DOUBLING.
  It's the finest resolution at which the tournament atom (3)
  and the doubling quantum (2²=4) are simultaneously visible.
""")

# The 12th root of 2
half_step = 2**(1/12)
print(f"  The musical half-step: 2^{{1/12}} = {half_step:.8f}")
print(f"  Major third (4 half-steps): 2^{{4/12}} = {2**(4/12):.8f} ≈ 5/4 = {5/4}")
print(f"  Minor third (3 half-steps): 2^{{3/12}} = {2**(3/12):.8f} ≈ 6/5 = {6/5}")
print(f"  Perfect fifth (7 half-steps): 2^{{7/12}} = {2**(7/12):.8f} ≈ 3/2 = {3/2}")
print()

# The Pythagorean comma: 12 fifths ≈ 7 octaves but not exactly
pythagorean_comma = (3/2)**12 / 2**7
print(f"  PYTHAGOREAN COMMA: (3/2)^12 / 2^7 = {pythagorean_comma:.8f}")
print(f"  The comma ≈ {(pythagorean_comma - 1) * 100:.2f}% deviation")
print(f"  This is the ±1 atom: the irreducible gap between")
print(f"  the additive (12 steps) and multiplicative (7 doublings) worlds.")
print()

# ==========================================================================
# PART 8: THE ORTHOGONAL DECOMPOSITION OF THE NUMBER LINE
# ==========================================================================

print()
print("PART 8: THE ORTHOGONAL DECOMPOSITION")
print("-" * 60)
print()

print("""  Every real number x decomposes into:

    x = sgn(x) · |x|

  SIGN (the {-1, 0, +1} part):
    Which side of zero? Discrete. Three values.
    This is the VERTICAL axis of classification.

  MAGNITUDE (the half-step part):
    How far from zero? Continuous. [0, ∞).
    This is the HORIZONTAL axis of measurement.

  They are ORTHOGONAL: knowing the sign tells you nothing
  about the magnitude, and vice versa.

  In the complex plane, the same decomposition:

    z = |z| · e^{iθ}

  MAGNITUDE |z| (the radial direction):
    How far from the origin. The {-1, 0, +1} classification
    becomes: inside unit circle (-1), on circle (0), outside (+1).

  PHASE θ (the angular direction):
    Which direction. The half-step. Orthogonal to magnitude.
    The phase carries the OSCILLATORY information.

  For the gap function:
    g_p(x) at real x gives the SIGN classification.
    The COMPLEX roots of g_p give the PHASE (oscillation).
    Sign and phase are independent = orthogonal.

  The tournament evaluation at x = 2:
    Sign: g_p(2) = +1 (hyperbolic for all p)
    Phase: e^{i·arg(z_j)} where z_j are the complex roots
    The phase DECAYS (|z_j| < 1) so it's invisible at x = 2.
    But it EXISTS at finite n, as the oscillatory correction
    to the tribonacci growth.

  THE HALF-STEP IS THE PHASE.
  THE WHOLE STEP IS THE SIGN.
  THEY ARE ORTHOGONAL BECAUSE PHASE ⊥ MAGNITUDE.
""")

# ==========================================================================
# PART 9: THE DEEP STATEMENT
# ==========================================================================

print()
print("=" * 72)
print("  PART 9: THE DEEP STATEMENT")
print("=" * 72)
print()

print("""
  HALF-STEPS ARE ORTHOGONAL TO {-1, 0, 1} BECAUSE
  PHASE IS ORTHOGONAL TO SIGN.

  The integer classification {-1, 0, +1} answers: WHAT TYPE?
    Spherical, Euclidean, or hyperbolic.
    Negative, zero, or positive.
    Less, equal, or greater.

  The half-step answers: AT WHAT ANGLE?
    The oscillation frequency. The phase of the complex root.
    The direction in the plane that the real line can't see.

  These two pieces of information are INDEPENDENT AND PERPENDICULAR:
    If I tell you x is positive, you know NOTHING about its phase.
    If I tell you x has phase 45°, you know NOTHING about its sign.

  EVERYTHING in tournament theory decomposes this way:
    OCR = the sign part (97%, which score class)
    Residual = the phase part (3%, the oscillation within the class)

    H = trend + oscillation (Binet-like decomposition)
    trend = A·τⁿ (the sign, the real root)
    oscillation = 2Re(B·z₂ⁿ) (the phase, the complex roots)

    Score = sign (coarse classification)
    SRCP = sign + phase (finer, includes the half-step information)

  The IMAGINARY UNIT i is the GENERATOR of half-steps.
  It satisfies i² = -1: two half-steps = one sign change.
  It is PERPENDICULAR to the real axis: half ⊥ whole.

  The √ is the operation that CREATES half-steps from whole steps:
    √(-1) = i (the half-step of negation)
    √(2) = the half-step of doubling
    √(g_p) = the half-step of the gap function

  And the LOGARITHM COUNTS half-steps:
    log₂(x) = "how many half-doublings from 1 to x"
    For x on the k-nacci ladder:
      log₂(2 - ρ_p) ≈ -p = "p half-steps below the universal point"

  THE COMPLETE PICTURE:
    {-1, 0, +1} = Z (the integers) = SIGNS = WHOLE STEPS
    {-i, +i} and their powers = √(-1)ⁿ = PHASES = HALF-STEPS
    They span the complex plane together.
    Neither is sufficient alone.
    Together they give EVERYTHING.

  Much of math relies on this because EVERY mathematical structure
  has both a SIGN (discrete type) and a PHASE (continuous angle),
  and understanding requires both.
""")

print("opus-2026-03-21-S142")
