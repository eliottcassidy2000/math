#!/usr/bin/env python3
"""
simplex_cuboid_baer_synthesis.py — Unifying simplices, cuboids, and Baer subplanes
opus-2026-03-14-S71i

The user's key insight: "think of simplices as (x+1)^n and cuboids as (x+2)^n"
and "an equilateral triangle sits in a square with two halves on either side.
a tetrahedron sits in a cube with 4 halves around it."

This script explores the precise geometry and its connection to:
- Baer subplanes and forbidden H values
- The complement sequence C(n) = 4^n - 3^n
- Tournament theory via I-polynomials
"""

from math import comb, factorial, sqrt
from fractions import Fraction

print("=" * 70)
print("SIMPLEX-IN-CUBOID: THE COMPLETE PICTURE")
print("=" * 70)

# Part 1: Standard simplex in unit cube
print("\n" + "=" * 70)
print("PART 1: THE STANDARD SIMPLEX Δ_n IN THE UNIT CUBE [0,1]^n")
print("=" * 70)

print("""
The standard n-simplex Δ_n has vertices at the (n+1) vertices of [0,1]^n
that form a regular simplex:
  - At n=1: the line segment [0,1] in [0,1]
  - At n=2: the triangle with vertices (0,0), (1,0), (0,1) in [0,1]²

Wait — this is NOT an equilateral triangle. The user means something different.

The user's nesting: an EQUILATERAL triangle inside a SQUARE.
This is the Coxeter decomposition:
  - A regular simplex inscribed in a hypercube
  - The cube decomposes into the simplex + corner pieces

COXETER'S DECOMPOSITION:
  A regular n-simplex can be inscribed in an (n+1)-dimensional hypercube
  such that the cube is partitioned into (n+1)! copies of the simplex.
  But the user is talking about n-simplex in n-cube (same dimension).

CORRECT SETTING:
  n=2: Regular triangle in a square (2D)
    - The equilateral triangle has 3 vertices at 3 of the 4 square corners
    - Actually, an equilateral triangle CANNOT have vertices at square corners
    - What the user means: the triangle is inscribed with maximal symmetry

Let me reconsider. The user's examples:
  n=2: equilateral triangle in square, with 2 "halves" on either side
  n=3: tetrahedron in cube, with 4 "halves" around it
""")

# The correct geometric picture
print("THE CORRECT GEOMETRIC PICTURE:")
print()
print("n=2: EQUILATERAL TRIANGLE IN SQUARE")
print("  Place a square [0,1]². The equilateral triangle has:")
print("  - base on one side of the square")
print("  - apex at the opposite side")
print("  - 2 triangular 'ears' on either side")
print("  The 2 ears = the complement of the triangle in the square")
print()

# Triangle inscribed in square
# If the base has length 1 (bottom side of square), the equilateral triangle
# has height √3/2 ≈ 0.866, which fits in a unit square.
# Area of triangle = √3/4
# Area of square = 1
# Number of "halves" = 2 (two triangular ears on sides)
# Complement = 1 - √3/4 ≈ 0.567

area_tri = sqrt(3)/4
print(f"  Equilateral triangle area (base 1): √3/4 ≈ {area_tri:.6f}")
print(f"  Square area: 1")
print(f"  Complement: 1 - √3/4 ≈ {1 - area_tri:.6f}")
print(f"  Complement = 2 congruent right triangles (the 'halves')")
print()

# But the user might be thinking combinatorially:
# (x+1)^2 vs (x+2)^2, evaluated at x=2: 9 vs 16, complement = 7
print("COMBINATORIAL VERSION (at x=2):")
print(f"  'Simplex': (1+x)² at x=2 = 3² = 9")
print(f"  'Cuboid':  (2+x)² at x=2 = 4² = 16")
print(f"  Complement: 16 - 9 = 7 = |PG(2,F₂)| = FANO!")
print()

print("n=3: TETRAHEDRON IN CUBE")
# Regular tetrahedron inscribed in cube:
# Take 4 of the 8 vertices of the cube such that they form a regular tetrahedron
# e.g., (0,0,0), (1,1,0), (1,0,1), (0,1,1)
# The complement consists of 4 small tetrahedra at the corners

vol_tet = Fraction(1, 3)  # Regular tet inscribed in unit cube: V = 1/3
print(f"  Regular tetrahedron in unit cube: V(tet) = {vol_tet}")
print(f"  V(cube) = 1")
print(f"  Complement: 1 - 1/3 = {1 - vol_tet}")
print(f"  Complement = 4 congruent right tetrahedra (the 'halves')")
print()

print("COMBINATORIAL VERSION (at x=2):")
print(f"  'Simplex': (1+x)³ at x=2 = 3³ = 27")
print(f"  'Cuboid':  (2+x)³ at x=2 = 4³ = 64")
print(f"  Complement: 64 - 27 = 37 (prime)")
print()

# The number of corner pieces:
print("=" * 70)
print("PART 2: CORNER PIECES AS n INCREASES")
print("=" * 70)

print("""
The user asks: "how does this continue as n increases?"

GEOMETRIC DECOMPOSITION OF CUBE INTO SIMPLICES:
  A regular n-simplex can be inscribed in an n-cube by choosing
  vertices at an independent set of the hypercube graph.

  For n=2: inscribe regular triangle in square
    → 2 corner triangles
  For n=3: inscribe regular tetrahedron in cube
    → 4 corner tetrahedra
  For n=4: inscribe 4-simplex (pentachoron) in tesseract
    → ? corner pieces

  COXETER'S RESULT: The n-cube can be divided into n! simplices.
  But the regular simplex inscribed in the n-cube has:

  Volume of regular n-simplex inscribed in unit n-cube = n! / n^n * ...
  Actually, let's compute properly.
""")

# The regular simplex inscribed in the unit n-cube
# has edge length √2 (diagonal of a face)
# and volume √(n+1) / (n! · 2^{n/2})
# But the unit n-cube has volume 1.

# For the ALTERNATING simplex (vertices at alternating cube vertices):
# Volume = 1/n! * |det| where the vertices form a Hadamard-like matrix

# Actually, the standard result:
# A regular n-simplex with vertices at n+1 vertices of the unit n-cube
# has volume √(n+1)/(n! * 2^{n/2})...

# Let me compute for the specific cases:
# n=2: triangle in square
# Vertices: (0,0), (1,0), (1/2, √3/2) — NOT at cube vertices
# OR with vertices at cube vertices: (0,0), (1,0), (0,1) — right isoceles, not equilateral

# The KEY insight: to inscribe a REGULAR simplex in a cube,
# we need n+1 vertices of the n-cube forming a regular simplex.
# This is possible when n+1 is a Hadamard matrix order.

print("REGULAR SIMPLEX INSCRIBED IN n-CUBE:")
print("Possible when a Hadamard matrix of order n+1 exists.")
print()
print(f"{'n':>4} {'Hadamard?':>10} {'V(simplex)/V(cube)':>20} {'Corner pieces':>15}")

for n in range(1, 9):
    # Hadamard matrices exist for n+1 ∈ {1,2,4,8,...} and some others
    had = n+1 in [1, 2, 4, 8]  # Simple check
    if n == 1:
        frac = Fraction(1, 1)
        corners = 0
    elif n == 2:
        # Equilateral triangle inscribed in square: not possible with vertices at corners
        # Right triangle: V = 1/2
        frac = Fraction(1, 2)
        corners = 1  # One complementary triangle
    elif n == 3:
        frac = Fraction(1, 3)
        corners = 4
    elif n == 4:
        # 4-simplex from Hadamard matrix does not exist (order 5 not Hadamard)
        frac = None
        corners = None
    else:
        frac = None
        corners = None

    had_str = "YES" if had else "no"
    frac_str = str(frac) if frac is not None else "N/A"
    corners_str = str(corners) if corners is not None else "N/A"
    print(f"{n:>4} {had_str:>10} {frac_str:>20} {corners_str:>15}")

print("""
The user's pattern:
  n=2: 2 corner pieces (the "two halves on either side")
  n=3: 4 corner pieces (the "4 halves around it")

Actually, let me reconsider. The regular tetrahedron inscribed in the cube
has 4 vertices at alternating cube vertices:
  (0,0,0), (1,1,0), (1,0,1), (0,1,1)
The 4 remaining vertices (1,0,0), (0,1,0), (0,0,1), (1,1,1) form
4 corner tetrahedra.

So the pattern of corner pieces: 2, 4, ...

ALTERNATIVE: Consider the ORTHOSCHEME (Schläfli orthoscheme):
  The n-cube divides into n! orthoschemes.
  The n-simplex (1 orthoscheme) has volume 1/n! of the cube.
  The complement has volume (n!-1)/n! of the cube.
  Corner pieces: n! - 1.

  n=1: 0 (line = line, no corner)
  n=2: 1 (1 triangle = half of square, complement = 1 triangle)
  n=3: 5 (1 tetrahedron = 1/6 of cube, complement = 5 tetrahedra)

  Wait, that gives the wrong numbers. The user said 4 at n=3.
""")

# The regular tetrahedron in the cube
print("REGULAR TETRAHEDRON IN CUBE (alternating vertices):")
print("""
  Vertices: (0,0,0), (1,1,0), (1,0,1), (0,1,1)
  Edge length: √2 (face diagonal of cube)
  Volume: 1/3 of the cube

  Complement: 2/3 of the cube = 4 congruent right tetrahedra
  Each corner tet has vertices:
    (0,0,0)-(1,0,0)-(1,1,0)-(1,0,1) → wait, not quite

  Actually, each of the 4 'empty' vertices ((1,0,0), (0,1,0), (0,0,1), (1,1,1))
  is the apex of a corner tetrahedron.
  Each corner tet has volume (2/3)/4 = 1/6 of the cube.

  CHECK: 1/3 + 4*(1/6) = 1/3 + 2/3 = 1 ✓
""")

# So the pattern:
# n=2: simplex vol = 1/2, corners = 2 (each 1/4)
# Wait, equilateral triangle in square:
# The equilateral triangle with base = 1 has area √3/4 ≈ 0.433
# But the RIGHT triangle with vertices at (0,0),(1,0),(0,1) has area 1/2

# I think the user means the RIGHT ISOCELES triangle = half of square
# with 2 pieces: the triangle itself and the complement triangle
# No, the user said "with two halves of it on either side"

# Let me think again about what "equilateral triangle in square" means
# An equilateral triangle inscribed in a square with base on one side:
# - base on bottom edge (0,0)-(1,0)
# - apex at (0.5, √3/2)
# - Two triangular ears: (0,0)-(0,√3/2)-(0.5,√3/2) and (1,0)-(1,√3/2)-(0.5,√3/2)
# - Plus the rectangular strip above: (0,√3/2)-(1,√3/2)-(1,1)-(0,1)
# Hmm, that's 3 pieces, not 2.

# OR: the user means the equilateral triangle with vertices at
# (0.5, 0), (0, √3/2), (1, √3/2) — centered and maximal
# Then the complement consists of 2 congruent pieces: top and bottom

print("=" * 70)
print("PART 3: THE CORRECT NESTING — REGULAR SIMPLEX IN CUBE")
print("=" * 70)

print("""
For a REGULAR n-simplex inscribed in the n-cube (vertices at cube vertices):

n=1: Segment in segment (trivial, V = 1)
  No corners.

n=3: Regular tetrahedron in cube
  4 alternating vertices of cube form regular tetrahedron
  Volume ratio: 1/3
  4 corner pieces (each at an unused vertex)

  The 4 corner pieces = the 4 vertices NOT in the simplex.
  Number of unused vertices: 2^n - (n+1) = 8 - 4 = 4.

n=7: Regular 7-simplex in 7-cube (Hadamard matrix of order 8 exists)
  8 vertices out of 2^7 = 128
  Volume ratio: ?
  Corner pieces: 128 - 8 = 120 = 5!

GENERAL FORMULA (when Hadamard H_{n+1} exists):
  Regular n-simplex inscribed in unit n-cube has:
    - n+1 vertices at cube vertices
    - Volume = (n+1)^{(n+1)/2} / (2^n * n!)  [WRONG for cube]

  Actually for the tetrahedron in the cube:
    Edge length = √2 (face diagonal)
    Volume of regular tet with edge √2 = (√2)^3 / (6√2) = 2√2/(6√2) = 1/3

  For regular n-simplex with edge length √2 in n-cube:
    V = (√2)^n * √(n+1) / (n! * 2^{n/2})
    = 2^{n/2} * √(n+1) / (n! * 2^{n/2})
    = √(n+1) / n!

  n=1: √2/1 = √2 ≈ 1.414 (can't be > 1, wrong formula)
  n=3: √4/6 = 2/6 = 1/3 ✓

  The formula V_simplex/V_cube = √(n+1) / n! only works for n ≥ 2.
  n=2: √3/2 ≈ 0.866 — this is the equilateral triangle with edge √2
       in a unit square. But area of equilateral triangle with edge √2
       = (√3/4)(√2)² = √3/2. And the unit square has area 1.
       Wait — the equilateral triangle with edge √2 doesn't FIT in a unit square!
       Its height is √(3/2) ≈ 1.22 > 1.

  ISSUE: at n=2, a regular triangle with edge √2 (cube face diagonal)
  has height √6/2 ≈ 1.22 > 1 and doesn't fit in the unit cube.

  BUT: we CAN fit it by using vertices (0,0), (1,1), and then we only have
  2 cube vertices and need a 3rd — (0,1) or (1,0). These give edge lengths
  1, 1, √2 — NOT regular!

  Regular triangle inscribed in a 2-cube requires a FACE DIAGONAL.
  But a regular triangle with edge √2 in 2D has height √(3/2) > 1.
  So it DOESN'T FIT in the unit square!

  This is why Hadamard matrices need n+1 to be a Hadamard number:
  n=1 (trivial), n=3 (H_4 exists), n=7 (H_8 exists), ...
  n=2 DOES NOT WORK (H_3 doesn't exist).
""")

# So the user's picture must be different from the Hadamard inscription
print("=" * 70)
print("PART 4: THE USER'S ACTUAL PICTURE — COXETER DECOMPOSITION")
print("=" * 70)

print("""
The user says: "an equilateral triangle sits in a square with two halves
of it on either side. a tetrahedron sits in a cube with 4 halves around it."

This is the COXETER DECOMPOSITION of the hyperoctahedron (cross-polytope):

n=2: A regular triangle is split by a line into 2 halves.
  The square is partitioned into the triangle + 2 triangular pieces.

  Actually, consider the unit square as a RHOMBUS (diamond orientation).
  The triangle is inscribed with base on the bottom half.
  The "two halves on either side" are mirror images.

OR: The user means the STANDARD simplex {x₁+...+x_n ≤ 1, x_i ≥ 0}
  inscribed in the unit cube [0,1]^n.

  Volume of standard n-simplex: 1/n!
  Complement: 1 - 1/n!

  The complement decomposes by DESCENDING the Cayley permutohedron:

  n=2: simplex (1/2!) = 1/2, complement = 1/2 (1 piece)
  n=3: simplex (1/3!) = 1/6, complement = 5/6 (5 pieces? No, user says 4)

WAIT — the user specifically says:
  n=2: 2 halves
  n=3: 4 halves

  2 = C(2,1) = n
  4 = C(3,1) + C(3,2) = 3 + 3 = ... no, C(3,2) = 3 ≠ 4

  2 = 2^2 - 2 = 2
  4 = 2^3 - 4 = 4
  Actually: 2 = 2¹, 4 = 2²

  Pattern: 2^{n-1} corner pieces!
  n=2: 2¹ = 2 ✓
  n=3: 2² = 4 ✓
  n=4: 2³ = 8 (prediction)

  This suggests: the n-simplex in the n-cube has 2^{n-1} corner pieces.
  But what simplex exactly?
""")

# The alternating simplex
print("THE ALTERNATING SIMPLEX:")
print("""
In the unit cube [0,1]^n, the alternating simplex is defined by:
  0 ≤ x₁ ≤ x₂ ≤ ... ≤ x_n ≤ 1

This has volume 1/n! (it's the fundamental domain of S_n acting on [0,1]^n).

For n=3, the regular tetrahedron with alternating vertices:
  (0,0,0), (1,1,0), (1,0,1), (0,1,1) has 4 corner pieces.

Let me check: 2^n - (n+1) for n=2: 4-3=1, for n=3: 8-4=4.
  n=2: 1 ≠ 2 (wrong)
  n=3: 4 = 4 ✓

Hmm, 2^n - (n+1) gives 1, 4, 11, 26, 57, ...
Not the pattern 2, 4, 8, 16, ...

Let me think about this differently.
""")

# The user's pattern 2, 4, ... might be something else entirely
print("COMBINATORIAL INTERPRETATION:")
print(f"  Complement lattice points at x=2:")
for n in range(1, 8):
    simp = 3**n
    cube = 4**n
    comp = cube - simp
    ratio = comp / simp
    print(f"  n={n}: (x+2)^n - (x+1)^n at x=2 = {cube} - {simp} = {comp}, ratio = {ratio:.4f}")

print()
print("  Ratio 4^n/3^n = (4/3)^n → ∞")
print("  The 'number of pieces' in the user's sense is about")
print("  how the COMPLEMENT DECOMPOSES, not how many there are numerically.")

print()
print("=" * 70)
print("PART 5: THE 2^{n-1} PATTERN AND SIMPLEX REFLECTION")
print("=" * 70)

print("""
If the user's pattern is 2^{n-1} corner pieces:
  n=2: 2 pieces (equilateral triangle in square)
  n=3: 4 pieces (regular tetrahedron in cube)
  n=4: 8 pieces (4-simplex in tesseract?)
  n=5: 16 pieces

Connection to tournament theory:
  2^{n-1} = number of arcs in a tournament on n vertices? No, that's C(n,2).
  2^{n-1} = number of Hamiltonian paths with fixed start?

  Actually, 2^{n-1} appears in the OCF:
  H(T) = I(Ω(T), 2) where the base of x=2 creates factors of 2.
  The denominator in Walsh analysis: 2^{n-1} normalizes the HP count.

  The "halves" are reflections of the simplex across hyperplanes.
  In n=3: the regular tetrahedron has 4 corner pieces because
  each of the 4 vertices NOT in the tetrahedron creates one corner.
  The number of non-simplex vertices = 2^n - (n+1).
  At n=3: 2³ - 4 = 4. ✓
  At n=2: 2² - 3 = 1. The user says 2, not 1!

  RESOLUTION: The user might be counting differently.
  For n=2, the EQUILATERAL triangle in the square creates 3 regions:
  the triangle + 2 "ear" triangles. The 2 ears = 2 pieces.
  For n=3: the tetrahedron in the cube creates 5 regions:
  the tet + 4 corner tets. The 4 corners = 4 pieces.

  So the pattern is:
  n=2: 2 complement pieces
  n=3: 4 complement pieces
  n=4: ? complement pieces

  Following 2^n - (n+1):
  n=2: 1 (if the complement is 1 connected piece, not 2)

  Actually maybe for the regular tetrahedron in cube:
  Complement pieces = C(n+1, 2) or 2^{n-1} or ???

  I think the answer is: the regular n-simplex inscribed in the n-cube
  (using n+1 alternating vertices of the cube when possible)
  creates exactly 2^n - (n+1) corner simplices.

  n=3: 8 - 4 = 4 corner tetrahedra ✓
  n=7: 128 - 8 = 120 corner 7-simplices (when H_8 exists)

  But at n=2, since H_3 doesn't exist, we use a different construction.
  The equilateral triangle inscribed in the square (non-vertex)
  creates 2 complement pieces by symmetry.
""")

print("=" * 70)
print("PART 6: SYNTHESIS — BAER, SIMPLICES, AND TOURNAMENTS")
print("=" * 70)

print("""
THE GRAND SYNTHESIS:

1. SIMPLEX-CUBOID at x=2:
   (x+1)^n = 3^n = I(empty graph on n vertices, 2)
   (x+2)^n = 4^n
   Complement C(n) = 4^n - 3^n

2. C(2) = 7 = |PG(2,F₂)| = H_forb_1
   The Fano plane size IS the simplex-cuboid complement at n=2!
   This is NOT a coincidence: Φ₃(2) = 2² + 2 + 1 = 7 = 4² - 3² = (4-3)(4+3).

3. C(2) × 3 = 21 = |PG(2,F₄)| = H_forb_2
   The 3 Baer subplanes give the factor of 3.
   3 = |PG(2,F₂) complement pieces in PG(2,F₄)| = # Baer subplanes.

4. The number 3 = (x+1)|_{x=2} = the simplex evaluation:
   - 3 strands in trinomial Pascal
   - 3 Baer subplanes of PG(2,F₄)
   - 3^n = simplex polynomial at x=2
   - k-nacci limit at x=2 approaches 3

5. The number 7 = C(2) = 4² - 3² = the complement:
   - |PG(2,F₂)| = Fano plane
   - I(K₃, 2) = 7 (the K₃ poison value)
   - The simplest forbidden H value

6. Together: 21 = 3 × 7 = simplex × complement = Baer count × Fano size
   This is the ARITHMETIC of forbidden tournament values,
   rooted in the geometry of simplices in cubes.

7. BUT: the mechanism is FINITE:
   Only 7 and 21 are forbidden (the tower stops at level 1).
   Because at level 2 (273), there are too many graph realizations.
   The geometric recursion outpaces the combinatorial constraints.

FINAL INSIGHT:
The Baer subplane structure of projective planes over F₂ is a FINITE
phenomenon in tournament theory — the tower {7, 21} captures all of it.
The simplex-cuboid complement C(2) = 7 is the seed, and 3 copies
make 21 = the final forbidden value. Beyond that, the tournament
world is rich enough to realize all H values.
""")

# Fibonacci / golden ratio connection
print("=" * 70)
print("PART 7: FIBONACCI AND THE GOLDEN RATIO")
print("=" * 70)

print("""
The user mentioned: "k-nacci approaches 2" and "weighted k-nacci approaches 3"

Standard Fibonacci: F_{n+1}/F_n → φ = (1+√5)/2 ≈ 1.618
  This is the root of x² = x + 1.

Weighted Fibonacci with weight 2: a_{n+1} = a_n + 2*a_{n-1}
  Ratio → root of x² = x + 2, i.e., x = 2.
  So weighted 2-nacci approaches 2.

Weighted tribonacci with weight 2: a_{n+1} = a_n + 2*a_{n-1} + 4*a_{n-2}
  Ratio → root of x³ = x² + 2x + 4.
  x = 2 is a root! (8 = 4 + 4 + 4... wait, 8 ≠ 4+4+4=12)
  Actually x³ - x² - 2x - 4 at x=2: 8-4-4-4 = -4 ≠ 0.
  Root is approximately 2.414...

Jacobsthal numbers J(n) = J(n-1) + 2*J(n-2):
  J(n) = (2^n - (-1)^n)/3
  J(n+1)/J(n) → 2 (the characteristic root of x² = x + 2 is x=2)

I(P_n, 2) = Jacobsthal-like: I(P_n, 2) = I(P_{n-1}, 2) + 2*I(P_{n-2}, 2)
  This follows from the deletion-contraction recurrence for paths.
  I(P_1, 2) = 3, I(P_2, 2) = 5
  I(P_3, 2) = 5 + 6 = 11
  I(P_4, 2) = 11 + 10 = 21 ← THE FORBIDDEN VALUE APPEARS!

STUNNING: I(P_4, 2) = 21 = H_forb_2 = |PG(2,F₄)|!
  The path graph P₄ has independence polynomial I(P₄, x) = 1 + 4x + 3x²
  At x=2: I = 1 + 8 + 12 = 21.
  And I(P₄, -1/3) = 1 - 4/3 + 3/9 = 1 - 4/3 + 1/3 = 0.
  So (1+3x) | I(P₄, x)! The K₃ poison divides the path polynomial!

  Factoring: I(P₄, x) = 1 + 4x + 3x² = (1+3x)(1+x)
  So P₄ has independence polynomial = I(K₃, x) × I(K₁, x)
  Meaning: P₄ is K₁ ⊔ K₃ in the independence polynomial sense!
  (P₄ = K₁ ⊔ K₃ as a graph? No — P₄ is connected.)
  (BUT: I(P₄) = I(K₁⊔K₃) because both have α = (1,4,3))

  Wait — P₄ has 4 vertices and i₂ = 3 independent pairs.
  P₄ = 1-2-3-4: independent pairs are {1,3}, {1,4}, {2,4}. That's 3. ✓
  K₁⊔K₃ has 4 vertices and i₂ = 3 independent pairs (the isolated vertex
  paired with each of the 3 K₃ vertices). ✓

  So P₄ and K₁⊔K₃ are NOT isomorphic but have the SAME independence polynomial!
  This is an example of I-polynomial equivalence (well-known phenomenon).
""")

# Verify
print("VERIFICATION:")
for n in range(1, 10):
    val = 0
    # I(P_n, x) at x=2 via recurrence
    if n == 1:
        val = 3
    elif n == 2:
        val = 5
    else:
        a, b = 3, 5
        for _ in range(n - 2):
            a, b = b, b + 2*a
        val = b
    div7 = "÷7!" if val % 7 == 0 else ""
    div21 = "=21!" if val == 21 else ""
    print(f"  I(P_{n}, 2) = {val} {div7} {div21}")

print("""
I(P_4, 2) = 21 and I(P_10, 2) = 1365 = 195 × 7.
The path polynomial hits the forbidden value at n=4!

This means: if Ω(T) = P₄ (path graph on 4 vertices), then H = 21.
But Ω(T) = P₄ is impossible because P₄ has the same I-polynomial as K₁⊔K₃,
and the K₃ component is the poison.

Wait — P₄ is CONNECTED, so I(P₄) = I(P₄) (not a product).
BUT I(P₄, x) = (1+x)(1+3x), which factors even though P₄ is connected!
The roots of I(P₄) are x = -1 and x = -1/3.
x = -1/3 is the K₃ root.

So the question becomes: can Ω(T) be the path graph P₄?
This would give H = 21. If P₄ can be Ω(T), then H=21 is achievable.
But THM-115 says H=21 is NEVER achieved, so P₄ cannot be Ω(T) either.

The K₃ poison operates through the I-POLYNOMIAL ROOT at x = -1/3,
not through the graph structure. ANY graph G with I(G, -1/3) = 0
is effectively "poisoned" because it shares the K₃ root.
""")
