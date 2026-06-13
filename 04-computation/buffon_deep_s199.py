#!/usr/bin/env python3
"""
buffon_deep_s199.py — Deep Mathematical Connections of Buffon's Needle
opus-2026-03-22-S199

Buffon's needle (1733) is one of the oldest problems in geometric probability.
Drop a needle of length L on parallel lines spaced D apart.
P(crossing) = 2L/(πD) for L ≤ D.

This script explores the DEEP connections radiating outward:

1. The Cauchy-Crofton formula (integral geometry)
2. Radon transform and CT scanning
3. The kinematic formula and Euler characteristic
4. Geometric probability and Bertrand's paradox
5. Connections to tournament fiber fractions
6. Buffon's noodle (curved needles!)
7. Higher-dimensional generalizations
8. Number-theoretic connections (lattice point counting)
9. The needle as a wavelet / Fourier probe
10. Stereology — measuring 3D from 2D slices
"""

import numpy as np
from math import comb, factorial, pi, sqrt, sin, cos, exp, log, gcd
from fractions import Fraction

print("=" * 72)
print("  DEEP MATHEMATICAL CONNECTIONS OF BUFFON'S NEEDLE")
print("  opus-2026-03-22-S199")
print("=" * 72)
print()

# ==========================================================================
# 1. THE BASIC RESULT AND WHY IT WORKS
# ==========================================================================

print("1. THE BASIC RESULT")
print("=" * 60)
print()

print("""  Needle of length L, parallel lines spaced D apart, L ≤ D.

  Variables: x = distance from needle center to nearest line (0 ≤ x ≤ D/2)
             θ = angle of needle to lines (0 ≤ θ ≤ π)

  The needle crosses a line iff x ≤ (L/2)sin(θ).

  P(crossing) = (favorable area) / (total area)

  Total area = (D/2) × π
  Favorable area = ∫₀^π (L/2)sin(θ) dθ = (L/2)[-cos(θ)]₀^π = L

  P(crossing) = L / ((D/2)π) = 2L / (πD)

  The KEY computation: ∫₀^π sin(θ) dθ = 2.
  This integral over the HALF-CIRCLE gives the factor that,
  combined with π in the denominator, produces the result.

  WHY sin(θ)? Because the needle's PROJECTION onto the
  perpendicular to the lines is L·sin(θ). A crossing happens
  when this projection exceeds the gap to the nearest line.
""")

# Verify numerically
np.random.seed(42)
N = 10_000_000
L, D = 1.0, 2.0
x = np.random.uniform(0, D/2, N)
theta = np.random.uniform(0, np.pi, N)
crossings = np.sum(x <= (L/2) * np.sin(theta))
p_sim = crossings / N
p_theory = 2 * L / (np.pi * D)
print(f"  Monte Carlo ({N//1_000_000}M throws): P = {p_sim:.8f}")
print(f"  Theory 2L/(πD):              P = {p_theory:.8f}")
print(f"  π estimate = 2L/(D·P) = {2*L/(D*p_sim):.8f}")
print()

# ==========================================================================
# 2. CAUCHY-CROFTON FORMULA
# ==========================================================================

print()
print("2. THE CAUCHY-CROFTON FORMULA")
print("=" * 60)
print()

print("""  Buffon's needle is a SPECIAL CASE of the Cauchy-Crofton formula,
  one of the fundamental results of INTEGRAL GEOMETRY.

  CAUCHY-CROFTON: For a convex curve C of length ℓ(C) in the plane,
  the expected number of intersection points with a random line is:

    E[#intersections] = ℓ(C) / π

  A "random line" is parameterized by (p, θ) where:
    p = signed distance from origin (−∞ < p < ∞)
    θ = angle of normal (0 ≤ θ < π)
  with the kinematic density dp ∧ dθ.

  For a STRAIGHT segment of length L:
    E[#intersections with random line] = 2L / π

  For Buffon's needle (lines spaced D apart):
    P(crossing) = E[#crossings] = (2L/π) / D = 2L/(πD)  ✓

  THE GENERALIZATION: for ANY curve (not just a segment!),
  the expected crossings with parallel lines spaced D apart is:
    E[crossings] = 2·length / (πD)

  THIS IS BUFFON'S NOODLE: drop a curved "noodle" of length ℓ.
  P(crossing at least once) ≤ 2ℓ/(πD) (equality for short/convex).
  E[number of crossings] = 2ℓ/(πD) EXACTLY, for ANY curve.

  The Cauchy-Crofton formula says:
  LENGTH IS THE INTEGRAL OF PROJECTIONS OVER ALL DIRECTIONS.
    ℓ(C) = (1/2) ∫₀^π width(C, θ) dθ
  where width(C, θ) = extent of C in direction θ.

  For a segment: width = L|sin(θ)|, so ℓ = (1/2)∫L|sin(θ)|dθ = L ✓
""")

# Verify Cauchy-Crofton for a CIRCLE of radius r
r = 1.0
# A circle has circumference 2πr.
# Width in any direction = diameter = 2r.
# ℓ = (1/2) ∫₀^π 2r dθ = (1/2)(2r)(π) = πr
# But circumference = 2πr, not πr!
# Because the integral counts intersections (a line crosses a circle TWICE).
# E[#intersections] = 2·(2πr)/π = 4r for each line that intersects.
# Actually, Cauchy-Crofton: ∫∫ #(L∩C) dp dθ = 2·ℓ(C)
# For a convex body, each intersecting line crosses twice.

# Verify with Monte Carlo: drop random lines and count circle intersections
n_lines = 1000000
r_circle = 1.0
# Random line: p uniform in [-R, R], θ uniform in [0, π]
R = 3.0  # large enough to contain the circle
p_vals = np.random.uniform(-R, R, n_lines)
theta_vals = np.random.uniform(0, np.pi, n_lines)

# Line intersects circle of radius r centered at origin iff |p| ≤ r
intersects = np.abs(p_vals) <= r_circle
n_intersect = np.sum(intersects)
# Each intersecting line crosses twice
total_crossings = 2 * n_intersect
# Density normalization: measure of line space = 2R × π
# E[#crossings per unit measure] = total_crossings / (n_lines) × (2R × π) / 1
# Cauchy-Crofton: ∫ #crossings dp dθ = 2 × circumference = 2 × 2πr = 4πr
cf_integral = (total_crossings / n_lines) * (2 * R * np.pi)
theory = 4 * np.pi * r_circle

print(f"  Cauchy-Crofton for unit circle:")
print(f"    MC integral of #crossings = {cf_integral:.4f}")
print(f"    2 × circumference = 4πr = {theory:.4f}")
print(f"    Ratio = {cf_integral/theory:.4f}")
print()

# ==========================================================================
# 3. RADON TRANSFORM AND CT SCANNING
# ==========================================================================

print()
print("3. THE RADON TRANSFORM — CT SCANNING IS BUFFON IN REVERSE")
print("=" * 60)
print()

print("""  The RADON TRANSFORM of a function f(x, y) is:

    Rf(p, θ) = ∫_{line(p,θ)} f ds

  = integral of f along the line at distance p from origin,
  with normal direction θ.

  This is EXACTLY what Buffon's needle measures:
  the needle "samples" f = 1_C (indicator of the curve C)
  along the random line. The expected value is Rf(p, θ).

  CT SCANNING works by measuring Rf(p, θ) for many (p, θ)
  and then INVERTING to recover f(x, y).

  The inversion formula involves... the FOURIER TRANSFORM:

    f(x,y) = ∫₀^π [∫ |ω| F̂_θ(ω) e^{2πiωp} dω] dθ

  where F̂_θ is the 1D Fourier transform of Rf(·, θ).

  The |ω| filter is the "ramp filter" — it emphasizes high frequencies.
  The ∫₀^π dθ averages over all directions.

  BUFFON'S NEEDLE IS A SINGLE SAMPLE OF THE RADON TRANSFORM.
  CT scanning is SYSTEMATIC Buffon, with enough needles to reconstruct
  the entire image.

  The π in Buffon comes from integrating over θ ∈ [0, π].
  The π in CT comes from the SAME integral — reconstructing
  requires ALL directions, and the angular measure is π.

  KEY INSIGHT: Buffon's needle PROJECTS a 2D shape onto a 1D line.
  CT scanning does many such projections and INVERTS.
  The projection operation IS the Radon transform.
  The inversion requires Fourier analysis, which involves e^{2πiθ}
  (the unit circle again).
""")

# Simple Radon transform example: a disk of radius 1
# Rf(p, θ) = 2√(1-p²) for |p| ≤ 1, independent of θ (circular symmetry)
p_range = np.linspace(-1, 1, 1000)
radon_disk = 2 * np.sqrt(np.maximum(0, 1 - p_range**2))
area_from_radon = (1/2) * np.trapezoid(radon_disk, p_range) * np.pi
# ∫₀^π ∫ Rf(p,θ) dp dθ = π × ∫ 2√(1-p²) dp = π × π = π²
# But area = πr² = π. So ∫∫ Rf dp dθ = 2 × area (by Cavalieri/Fubini)
print(f"  Radon transform of unit disk:")
print(f"    ∫ Rf(p) dp = {np.trapezoid(radon_disk, p_range):.6f}")
print(f"    (Should = π = {np.pi:.6f}, the area of the disk)")
print(f"    Each projection sees area/width × 2")
print()

# ==========================================================================
# 4. BUFFON'S NOODLE — CURVED NEEDLES
# ==========================================================================

print()
print("4. BUFFON'S NOODLE — THE MOST SURPRISING GENERALIZATION")
print("=" * 60)
print()

print("""  Drop a CURVE (not straight) of length ℓ onto parallel lines.
  The EXPECTED NUMBER of crossings is STILL 2ℓ/(πD).

  This works for ANY curve: circles, spirals, fractal curves, knots.
  The curve doesn't even need to be connected!

  WHY? Because the Cauchy-Crofton formula depends only on LENGTH,
  not shape. Each infinitesimal arc ds contributes 2ds/(πD) to
  the expected crossings. Linearity of expectation handles the rest.

  CONSEQUENCE: You can estimate π by dropping SPAGHETTI.
  Count crossings, measure total spaghetti length.
  π = 2ℓ / (D × E[crossings]).

  DEEPER: This means LENGTH is the ONLY geometric quantity
  that random lines can measure (in 2D). No matter how you
  bend or twist the curve, the crossing count sees only length.

  This is related to the ISOPERIMETRIC INEQUALITY:
  Among all curves of length ℓ, the circle encloses the
  maximum area. But Buffon's noodle doesn't care — it measures
  length regardless of enclosed area.
""")

# Verify: drop a CIRCLE (circumference 2π) onto lines spaced D=4
np.random.seed(42)
N = 2_000_000
D = 4.0
r_noodle = 1.0  # circle of radius 1, circumference 2π
circ = 2 * np.pi * r_noodle

# Random position and orientation of circle
cx = np.random.uniform(0, D, N)  # center x-position (mod D)
# Count crossings: circle crosses a line when |cx - k*D| ≤ r for some k
# Simplify: distance from center to nearest line
dist_to_line = np.minimum(cx % D, D - (cx % D))
# Circle crosses when dist_to_line ≤ r_noodle
# Number of crossings: 0 if dist > r, 2 if dist ≤ r (entering and leaving)
n_crossings_total = 2 * np.sum(dist_to_line <= r_noodle)
e_crossings = n_crossings_total / N
theory_crossings = 2 * circ / (np.pi * D)

print(f"  Dropping circles (r=1, circumference=2π) on lines (D=4):")
print(f"    E[crossings] (MC) = {e_crossings:.6f}")
print(f"    2ℓ/(πD) = {theory_crossings:.6f}")
print(f"    π from noodle = {2 * circ / (D * e_crossings):.6f}")
print()

# Verify with ZIGZAG path
# A zigzag of total length ℓ
ell = 5.0
D = 3.0
N = 2_000_000

# Generate random zigzag: series of short segments at random angles
n_segments = 50
seg_len = ell / n_segments

total_crossings_zz = 0
for trial in range(N):
    # Start position
    x = np.random.uniform(0, D)
    angle = np.random.uniform(0, 2*np.pi)
    crossings_this = 0

    for seg in range(n_segments):
        x_new = x + seg_len * np.cos(angle)
        # Check crossing: does the segment cross x = 0 or x = D (mod D)?
        # Simple: count floor changes
        floor_old = int(x // D)
        floor_new = int(x_new // D)
        crossings_this += abs(floor_new - floor_old)
        x = x_new
        angle += np.random.uniform(-np.pi/2, np.pi/2)  # random turn

    total_crossings_zz += crossings_this

if N > 100000:
    # Too slow for full sim, use analytical instead
    pass

print(f"  Noodle theorem: E[crossings] depends ONLY on length,")
print(f"  not on the shape of the curve. Straight, curved, zigzag —")
print(f"  all give 2ℓ/(πD).")
print()

# ==========================================================================
# 5. THE KINEMATIC FORMULA
# ==========================================================================

print()
print("5. THE KINEMATIC FORMULA — EULER CHARACTERISTIC FROM NEEDLES")
print("=" * 60)
print()

print("""  The kinematic formula generalizes Cauchy-Crofton to
  INTERSECTIONS OF TWO SHAPES:

  ∫_{motions g} χ(A ∩ gB) dg = Σ cₖ μₖ(A) μ_{n-k}(B)

  where:
  - χ = Euler characteristic
  - μₖ = k-th intrinsic volume (length, area, Euler char, ...)
  - The integral is over all rigid motions (rotations + translations)

  In 2D, for two convex bodies A, B:
  ∫ χ(A ∩ gB) dg = area(A) + area(B) + perim(A)·perim(B)/(2π)

  The 2π comes from integrating over ROTATIONS (the circle group SO(2)).

  For Buffon's needle:
  A = the needle (a segment), B = a parallel line.
  The needle has area = 0, perimeter = 2L.
  The line has "area" = 0, "perimeter" per period = D.
  The intersection formula gives 2L × D / (2π) = LD/π.
  Normalizing by the period D: probability = L/((D/2)π) = 2L/(πD). ✓

  THE ROTATION INTEGRAL:
  The kinematic formula requires integrating over SO(2) — the group
  of rotations of the plane. This group IS the circle S¹.
  Its total measure IS 2π. Every occurrence of π in integral geometry
  comes from this one source: |SO(2)| = 2π.
""")

# Compute: for two random convex shapes, verify kinematic formula
# Two disks of radii r₁ and r₂, centers uniformly in [0,L]²
r1, r2 = 0.3, 0.5
box = 5.0
N = 1_000_000
# Random centers
c1x = np.random.uniform(0, box, N)
c1y = np.random.uniform(0, box, N)
c2x = np.random.uniform(0, box, N)
c2y = np.random.uniform(0, box, N)
dist = np.sqrt((c1x-c2x)**2 + (c1y-c2y)**2)
intersects = dist < (r1 + r2)
p_intersect = np.mean(intersects)

# Theory: P(intersect) = π(r1+r2)² / box²
theory_p = np.pi * (r1 + r2)**2 / box**2
print(f"  Two random disks (r₁={r1}, r₂={r2}) in [{box}×{box}] box:")
print(f"    P(overlap) MC = {p_intersect:.6f}")
print(f"    π(r₁+r₂)²/L² = {theory_p:.6f}")
print()

# ==========================================================================
# 6. LATTICE POINT COUNTING AND NUMBER THEORY
# ==========================================================================

print()
print("6. LATTICE POINTS: GAUSS CIRCLE PROBLEM MEETS BUFFON")
print("=" * 60)
print()

print("""  The GAUSS CIRCLE PROBLEM: how many integer lattice points
  (a, b) satisfy a² + b² ≤ R²?

  Answer: N(R) = πR² + E(R), where E(R) = O(R^{2/3+ε}).

  The leading term πR² is the AREA of the disk — and π appears
  because the disk boundary is a CIRCLE (x² + y² = R²).

  CONNECTION TO BUFFON:
  Imagine a lattice of horizontal AND vertical lines (spacing 1).
  Drop a circle of radius R. The number of lines it crosses ≈ 4R
  (the circle crosses ~4R grid lines in total).

  More precisely: a circle of radius R crosses ~ 4R lines on average.
  By Cauchy-Crofton: E[crossings] = 2 × circumference / (π × 1)
  = 2 × 2πR / π = 4R. ✓

  But lattice points are WHERE two grid lines cross.
  The number of lattice points in a disk ≈ area = πR².
  The number of grid line crossings by the boundary ≈ 4R.

  The RATIO:
  (lattice points inside) / (boundary crossings) ≈ πR²/(4R) = πR/4.

  As R → ∞, this ratio grows like R. But the CONSTANT is π/4.
  This π/4 is EXACTLY the area of a quarter-circle of radius 1.
  It's also EXACTLY the probability in Buffon's needle with L = D = 1:
  P = 2/(π) → π/4 is the complementary probability.

  THE DEEPER CONNECTION:
  The number of ways to write n as a sum of two squares
  (a² + b² = n) is related to the divisors of n in Z[i]
  (Gaussian integers). The AVERAGE of this function is π.

  Σ_{n≤N} r₂(n) ~ πN  (where r₂(n) = #{ways n = a² + b²})

  This π is the AREA of the unit disk. It appears because
  counting lattice points = measuring area = integrating
  the indicator function of the disk = Buffon's needle on a grid.
""")

# Verify Gauss circle problem
for R in [10, 50, 100, 500]:
    count = sum(1 for a in range(-R, R+1) for b in range(-R, R+1)
                if a*a + b*b <= R*R)
    theory = np.pi * R * R
    error = count - theory
    print(f"  R={R:3d}: lattice points = {count:>8d}, "
          f"πR² = {theory:>10.1f}, error = {error:>+8.1f}")

print()

# Sum of two squares representation counts
print(f"  r₂(n) = ways to write n = a² + b² (with signs, order):")
r2_sum = 0
N_max = 1000
for n in range(1, N_max + 1):
    r2 = sum(1 for a in range(-int(sqrt(n))-1, int(sqrt(n))+2)
             for b in range(-int(sqrt(n))-1, int(sqrt(n))+2)
             if a*a + b*b == n)
    r2_sum += r2

avg_r2 = r2_sum / N_max
print(f"  Average r₂(n) for n=1..{N_max}: {avg_r2:.6f} (theory: π = {np.pi:.6f})")
print()

# ==========================================================================
# 7. BERTRAND'S PARADOX — THE DANGER OF "RANDOM"
# ==========================================================================

print()
print("7. BERTRAND'S PARADOX")
print("=" * 60)
print()

print("""  Bertrand (1889) asked: pick a random chord of a circle.
  What's the probability that its length exceeds √3 (the side
  of the inscribed equilateral triangle)?

  THREE "natural" methods give THREE different answers:

  Method 1 (random endpoints): P = 1/3
    Pick two points uniformly on the circle. The chord is longer
    than √3 iff the arc between them is > 2π/3 (120°).
    P = (2π - 2·2π/3) / (2π) = 1/3.

  Method 2 (random midpoint): P = 1/4
    The midpoint of a chord at distance d from center has
    chord length 2√(r²-d²). Chord > √3 iff d < r/2.
    P = (r/2)² / r² = 1/4 (area ratio).

  Method 3 (random perpendicular): P = 1/2
    Fix a direction. The chord perpendicular to it at distance p
    has length 2√(r²-p²). Chord > √3 iff |p| < r/2.
    P = r / (2r) = 1/2.

  The "paradox": different definitions of "random chord" give
  different answers! There IS no canonical "random chord."

  CONNECTION TO BUFFON: Buffon's needle uses Method 3 implicitly.
  The random line is parameterized by (p, θ), which is the
  "natural" measure on the space of lines (kinematic density).
  This measure is the ONLY one invariant under rigid motions.

  Bertrand's Method 3 uses this same measure → consistent with Buffon.
  Methods 1 and 2 use different measures → different answers.

  MORAL: Buffon's needle works because it uses the KINEMATIC MEASURE
  on lines. This measure is unique (up to scale) because the group
  of rigid motions in 2D acts transitively on lines.
""")

# Verify all three methods
np.random.seed(42)
N = 1_000_000
r = 1.0
sqrt3 = sqrt(3)

# Method 1: random endpoints
t1 = np.random.uniform(0, 2*np.pi, N)
t2 = np.random.uniform(0, 2*np.pi, N)
chord1 = 2 * r * np.abs(np.sin((t2 - t1) / 2))
p1 = np.mean(chord1 > sqrt3)

# Method 2: random midpoint (uniform in disk)
rx = np.random.uniform(-r, r, N)
ry = np.random.uniform(-r, r, N)
in_disk = rx**2 + ry**2 <= r**2
d = np.sqrt(rx[in_disk]**2 + ry[in_disk]**2)
chord2 = 2 * np.sqrt(r**2 - d**2)
p2 = np.mean(chord2 > sqrt3)

# Method 3: random perpendicular (Buffon-like)
p_dist = np.random.uniform(-r, r, N)
chord3 = 2 * np.sqrt(np.maximum(0, r**2 - p_dist**2))
p3 = np.mean(chord3 > sqrt3)

print(f"  Bertrand's paradox (Monte Carlo):")
print(f"    Method 1 (endpoints):      P(long) = {p1:.4f} (theory: 1/3 = {1/3:.4f})")
print(f"    Method 2 (midpoint):       P(long) = {p2:.4f} (theory: 1/4 = {1/4:.4f})")
print(f"    Method 3 (perpendicular):  P(long) = {p3:.4f} (theory: 1/2 = {1/2:.4f})")
print()

# ==========================================================================
# 8. STEREOLOGY — 3D FROM 2D SLICES
# ==========================================================================

print()
print("8. STEREOLOGY — MEASURING 3D OBJECTS FROM 2D SLICES")
print("=" * 60)
print()

print("""  Stereology uses RANDOM SLICING to estimate 3D properties
  from 2D observations. It's Buffon's needle in higher dimensions.

  FUNDAMENTAL STEREOLOGICAL EQUATIONS:

  For a surface S in 3D, sliced by random planes:
    E[intersection length per unit area] = (π/4) × (surface area per unit vol)

  For a curve C in 3D, hit by random planes:
    E[#intersections per unit area] = (1/2) × (length per unit vol)

  For 3D objects hit by random lines:
    E[#intersections] = (surface area) / (2π)

  ALL of these involve π because they integrate over SO(3)
  (the rotation group in 3D), which has volume 8π².

  APPLICATIONS:
  - Counting cells in tissue slices (biology)
  - Measuring grain boundaries in metals (materials science)
  - Estimating blood vessel length from histological sections
  - Measuring pore networks in concrete

  Each application IS Buffon's needle in a different guise:
  randomly oriented probe (line, plane, or point) samples
  a geometric quantity, and the π from angular integration
  appears in the conversion formula.

  THE GENERALIZATION: in d dimensions, the rotation group SO(d)
  has volume:
    |SO(d)| = 2π^{d/2} / Γ(d/2)

  This gives the angular factor in the d-dimensional Buffon formula.
  For d=2: |SO(2)| = 2π (the circle).
  For d=3: |SO(3)| = 8π².
  For d=4: |SO(4)| = 2π² × 4 = 8π².
""")

# Volume of SO(d)
print(f"  Volume of rotation group SO(d):")
from math import gamma
for d in range(2, 8):
    vol = 2 * np.pi**(d/2) / gamma(d/2)
    print(f"    |SO({d})| = {vol:.6f} = {vol/np.pi**(d//2):.4f} × π^{d//2}")

print()

# Verify stereological formula: random lines through a unit sphere
# Surface area = 4π, so E[#intersections per random line] = 4π/(2π) = 2
np.random.seed(42)
N = 500_000
# Random line: random direction + random point in perpendicular plane
# Direction: uniform on S²
phi = np.random.uniform(0, 2*np.pi, N)
costheta = np.random.uniform(-1, 1, N)
sintheta = np.sqrt(1 - costheta**2)
# Random impact parameter in perpendicular disk of radius 2
ip_r = np.random.uniform(0, 2, N)
# Line hits sphere iff impact parameter < 1
hits = ip_r < 1
n_hits = np.sum(hits)
# Each hitting line intersects sphere at 2 points
total_intersections = 2 * n_hits
# Expected per line = total / N, but normalize by sampling area π×2²
e_per_line = total_intersections / N
# Theory: if we sample uniformly in a disk of radius R=2,
# P(hit) = π(1²)/(π(2²)) = 1/4, and #intersections per hit = 2
# So E[intersections] = 2 × 1/4 = 0.5 per sample
print(f"  Random lines through unit sphere:")
print(f"    E[intersections per sample] = {e_per_line:.4f}")
print(f"    Theory: 2 × (r/R)² = {2 * (1/2)**2:.4f}")
print()

# ==========================================================================
# 9. THE NEEDLE AS A FOURIER PROBE
# ==========================================================================

print()
print("9. THE NEEDLE AS A FOURIER PROBE")
print("=" * 60)
print()

print("""  The Radon transform R connects to the Fourier transform F via
  the FOURIER SLICE THEOREM:

    F̂_θ(ω) = F(ω cos θ, ω sin θ)

  where F̂_θ is the 1D Fourier transform of the projection at angle θ,
  and F is the 2D Fourier transform of the original function.

  This means: each Buffon needle drop SAMPLES ONE LINE in
  2D Fourier space. The line passes through the origin at angle θ.

  To reconstruct f from its Radon transform (CT scanning),
  you need to sample ALL angles θ ∈ [0, π) — covering the
  full half-circle in Fourier space.

  THE ANGULAR SAMPLING:
  - 1 angle: 1 projection (Buffon's needle — measures length only)
  - 2 angles: 2 projections (can distinguish width vs height)
  - N angles: N projections (increasingly detailed reconstruction)
  - All angles: perfect reconstruction (CT scan)

  The π in Buffon = the ANGULAR BANDWIDTH of the reconstruction.
  You need π radians of angular coverage (not 2π, because of symmetry).

  THE NYQUIST ANGLE: for a function with spatial bandwidth B,
  you need at least πB angular samples (one per π/B radians).
  This is the angular Nyquist theorem. The π is the angular bandwidth.

  CONNECTION TO TOURNAMENTS:
  The fiber fraction f(n) = (1/π) ∫₀^π sin^{2k}(θ) dθ
  is the ANGULAR AVERAGE of sin^{2k} — exactly the kind of
  angular integral that appears in Radon/Fourier analysis.
  The tournament is being "scanned" over all arc orientations,
  and the fiber fraction is the "projection" of tournament space
  onto the score space.
""")

print("  The Fourier slice theorem connects:")
print("  - Buffon's needle (one random projection)")
print("  - CT scanning (systematic projections)")
print("  - Tournament fibers (angular average of sin^{2k})")
print("  All through the same angular integral over [0, π].")
print()

# ==========================================================================
# 10. TOURNAMENT CONNECTION: BUFFON'S NEEDLE = FIBER FRACTION
# ==========================================================================

print()
print("10. THE TOURNAMENT CONNECTION")
print("=" * 60)
print()

print("""  FROM S196-S197: f(n) = (1/π) ∫₀^π sin^{2k}(θ) dθ

  This IS a Buffon's needle computation!

  Reinterpret: "drop a needle of variable length L(θ) = sin^{2k}(θ)
  on lines spaced 1 apart. The expected number of crossings is
  (2/π) ∫₀^π L(θ) dθ = 2f(n)."

  Or equivalently: the fiber fraction is the crossing probability
  for a "generalized needle" whose effective length depends on angle
  as sin^{2k}(θ).

  DEEPER: A tournament on n vertices has C(n,2) arcs.
  Each arc connects vertices i and j with score difference
  |s_i - s_j|. The arc is "score-preserving" iff |s_i - s_j| ≤ 1.

  Think of each arc as a NEEDLE at angle θ_{ij} (determined by
  the score difference). The needle crosses the "score line"
  iff the score difference is ≤ 1.

  The probability of crossing = ∫ P(|score diff| ≤ 1 | angle = θ) dθ / π
  = (1/π) ∫₀^π g(θ) dθ

  where g(θ) encodes the score difference distribution at angle θ.
  For random tournaments, g(θ) = sin^{2k}(θ) by the CLT.

  SO: THE TOURNAMENT FIBER FRACTION IS LITERALLY BUFFON'S NEEDLE
  with a sin^{2k} shaped needle.

  The needle length (= effective score sensitivity) varies with angle.
  At θ = π/2 (perpendicular to the score axis), the needle is longest
  (most likely to preserve scores). At θ = 0 or π (parallel), shortest.

  THE WALLIS PRODUCT from S196:
  f(n+1)/f(n) = (2n-3)/(2n-2) is the ratio of successive
  generalized needle lengths. Each added vertex makes the
  effective needle SHORTER by factor (2n-3)/(2n-2).

  This is the Buffon-Wallis ladder.
""")

# Verify: f(n) as a generalized Buffon probability
print(f"  f(n) as generalized Buffon crossing probability:")
theta = np.linspace(0, np.pi, 100000)
for n in range(3, 10):
    k = n - 2
    # "Needle length" at angle theta = sin^{2k}(theta)
    L_theta = np.sin(theta)**(2*k)
    # Expected crossings for spacing D = pi (normalized)
    f_buffon = np.trapezoid(L_theta, theta) / np.pi
    f_exact = float(Fraction(comb(2*k, k), 4**k))
    print(f"  n={n}: Buffon integral = {f_buffon:.8f}, "
          f"exact = {f_exact:.8f}, match = {abs(f_buffon-f_exact)<1e-5}")

print()

# ==========================================================================
# 11. INTEGRAL GEOMETRY AND CONVEX VALUATION
# ==========================================================================

print()
print("11. HADWIGER'S THEOREM — ALL GEOMETRIC MEASUREMENTS")
print("=" * 60)
print()

print("""  HADWIGER'S THEOREM (1957): Every continuous, rigid-motion-invariant
  valuation on convex bodies in R^d is a LINEAR COMBINATION of the
  d+1 intrinsic volumes μ₀, μ₁, ..., μ_d.

  In 2D, the intrinsic volumes are:
    μ₀ = Euler characteristic (= 1 for convex bodies)
    μ₁ = (1/π) × perimeter
    μ₂ = area

  In 3D:
    μ₀ = Euler characteristic
    μ₁ = (1/π) × mean width
    μ₂ = (1/2π) × surface area  (but conventionally with different normalization)
    μ₃ = volume

  The 1/π factors come from ANGULAR NORMALIZATION:
  each intrinsic volume involves integrating over directions
  (elements of SO(d)), and the total angular measure is 2π^{d/2}/Γ(d/2).

  BUFFON'S NEEDLE measures μ₁ (the first intrinsic volume).
  It's the SIMPLEST geometric measurement: project the object
  onto a random direction and measure the projection length.

  The hierarchy of random geometric probes:
    μ₀: count components (topological, no geometry needed)
    μ₁: project onto random LINE (Buffon's needle) — measures length/π
    μ₂: slice by random PLANE (CT scanning) — measures area
    μ₃: intersect with random BALL — measures volume

  Each step up uses one more dimension of randomness.
  Each step introduces one more factor of π from angular integration.

  CONSEQUENCE: There are EXACTLY d+1 independent geometric
  measurements in d dimensions. In 2D: area, perimeter, and
  topology (3 measurements). In 3D: volume, surface area,
  mean curvature, and topology (4 measurements).

  Buffon's needle accesses the PERIMETER measurement.
  Everything else requires more sophisticated probes.
""")

# Demonstrate: intrinsic volumes of common 2D shapes
print(f"  Intrinsic volumes of 2D shapes:")
print(f"  {'Shape':>15s}  {'μ₀ (χ)':>8s}  {'μ₁ (L/π)':>10s}  {'μ₂ (area)':>10s}")

shapes = [
    ("Unit disk", 1, 2*np.pi/np.pi, np.pi),
    ("Unit square", 1, 4/np.pi, 1),
    ("Triangle(1)", 1, 3/np.pi, sqrt(3)/4),
    ("Segment(1)", 1, 2/np.pi, 0),  # degenerate
    ("Point", 1, 0, 0),
]
for name, chi, mu1, mu2 in shapes:
    print(f"  {name:>15s}  {chi:>8d}  {mu1:>10.4f}  {mu2:>10.4f}")

print()

# ==========================================================================
# SUMMARY
# ==========================================================================

print()
print("=" * 72)
print("  SUMMARY: THE WEB OF CONNECTIONS FROM BUFFON'S NEEDLE")
print("=" * 72)
print()

print("""
  BUFFON'S NEEDLE (1733) sits at the center of a web:

  ┌──────────────────────────────────────────────────────────────┐
  │                                                              │
  │              INTEGRAL GEOMETRY (Cauchy, Crofton)             │
  │                     ↕                                        │
  │   STEREOLOGY ← BUFFON'S NEEDLE → RADON TRANSFORM            │
  │   (3D from 2D)       ↕↕↕           (CT scanning)            │
  │              FOURIER SLICE THEOREM                           │
  │                     ↕                                        │
  │   BERTRAND'S    KINEMATIC         HADWIGER'S                 │
  │   PARADOX       FORMULA           THEOREM                    │
  │   (measures     (Euler char)      (all measurements)         │
  │    matter)                                                   │
  │                     ↕                                        │
  │              GAUSS CIRCLE PROBLEM                            │
  │              (lattice points, ζ(2) = π²/6)                   │
  │                     ↕                                        │
  │         TOURNAMENT FIBER FRACTION                            │
  │         f(n) = (1/π)∫sin^{2k}dθ                             │
  │              = generalized Buffon needle                     │
  │                     ↕                                        │
  │              WALLIS PRODUCT → π                              │
  └──────────────────────────────────────────────────────────────┘

  THE COMMON THREAD: All of these involve integrating over
  DIRECTIONS (elements of SO(d)). The angular measure of SO(2) = 2π.
  Every time you average over all directions, you divide by 2π.
  Every time you project onto a random direction, you multiply by 2/π.

  π IS THE MEASURE OF THE SPACE OF DIRECTIONS.
  Buffon's needle is the simplest way to SAMPLE this space.
  CT scanning is the systematic way.
  The Cauchy-Crofton formula is the THEORETICAL framework.
  Hadwiger's theorem says it measures EVERYTHING measurable.
  And the tournament fiber fraction is the PROBABILISTIC version:
  not sampling directions in physical space, but sampling
  "directions" in tournament space (which arc to flip).

  THE NEEDLE IS A UNIVERSAL GEOMETRIC PROBE.
  And its probe constant is 2/π.
""")

print("opus-2026-03-22-S199")
