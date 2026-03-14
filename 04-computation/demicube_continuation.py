#!/usr/bin/env python3
"""
demicube_continuation.py — opus-2026-03-14-S71g

The user asks: triangle in square (2 halves), tetrahedron in cube (4 halves).
How does this continue?

ANSWER: The pattern is the DEMICUBE — the convex hull of cube vertices
with even coordinate sum. At n=3 this is a regular tetrahedron.

Key question: What is Vol(demicube_n) / Vol(cube_n)?
And how many "corner pieces" are there?

Also: deeper investigation of the Fibonacci-Jacobsthal category structure.
"""

import math
import numpy as np
from itertools import product
from fractions import Fraction

print("=" * 70)
print("DEMICUBE CONTINUATION AND FIBONACCI CATEGORY")
print("opus-2026-03-14-S71g")
print("=" * 70)

# ============================================================
# Part 1: Demicube volumes — exact computation
# ============================================================
print("\n" + "=" * 70)
print("PART 1: DEMICUBE VOLUMES")
print("=" * 70)

# The demicube (half-cube) h_n is the convex hull of vertices of [0,1]^n
# with even coordinate sum.
#
# FACT: The volume of h_n in [0,1]^n is EXACTLY 1/2 for all n ≥ 2.
#
# PROOF: The map (x_1,...,x_n) → (x_1,...,x_{n-1}, 1-x_n) is a
# volume-preserving involution that swaps even/odd parity vertices.
# So the even-half and odd-half each have volume 1/2.
#
# Wait — that's the volume of the REGION, not the POLYTOPE.
# The demicube polytope is the convex hull of even-parity vertices.
# But the region {x : Σx_i even} is NOT the same as conv(even vertices).
#
# Example: at n=3, even vertices are (0,0,0),(1,1,0),(1,0,1),(0,1,1).
# Their convex hull is a tetrahedron with volume 1/3, not 1/2.
# But the "even region" (where the rounding to nearest even vertex is closer)
# is different from the convex hull.

# Let me compute the actual convex hull volumes using scipy
try:
    from scipy.spatial import ConvexHull

    print("\nDemicube convex hull volumes:")
    for n in range(2, 9):
        verts = [v for v in product([0, 1], repeat=n) if sum(v) % 2 == 0]
        pts = np.array(verts, dtype=float)

        if len(pts) < n + 1:
            print(f"  n={n}: only {len(pts)} vertices, not enough for full-dim hull")
            continue

        hull = ConvexHull(pts)
        vol = hull.volume
        n_verts = len(verts)
        n_complement = 2**n - n_verts

        # Try to express as simple fraction
        frac = Fraction(vol).limit_denominator(10000)

        print(f"  n={n}: {n_verts} vertices, Vol = {vol:.10f} ≈ {frac}, "
              f"complement verts = {n_complement}")

        # The complement volume
        comp_vol = 1.0 - vol
        comp_frac = Fraction(comp_vol).limit_denominator(10000)
        print(f"         Complement vol = {comp_vol:.10f} ≈ {comp_frac}")

        # Volume per complement vertex (if each defines one piece)
        if n_complement > 0:
            vol_per = comp_vol / n_complement
            vp_frac = Fraction(vol_per).limit_denominator(100000)
            print(f"         Vol per complement vertex = {vol_per:.10f} ≈ {vp_frac}")

except ImportError:
    print("  scipy not available")

# ============================================================
# Part 2: Analytic formula for demicube volume
# ============================================================
print("\n" + "=" * 70)
print("PART 2: ANALYTIC FORMULA")
print("=" * 70)

# The n-dimensional demicube (half-cube) inscribed in [0,1]^n:
# Known formula: Vol(h_n) = 1/2 * (1 - 1/2^{n-1}) for n ≥ 2?
# Let me check:
# n=2: 1/2 * (1 - 1/2) = 1/4? But we computed... let's see.
# At n=2, even vertices = (0,0), (1,1). Convex hull = diagonal segment.
# This has 1D "volume" = √2, but as a 2D object, volume = 0.
# So the convex hull is degenerate at n=2!

# Actually for n=2, the even vertices (0,0) and (1,1) form a line segment,
# which has zero 2D volume. So the formula breaks at n=2.

# For n=3: Vol = 1/3 (from our computation)
# Formula attempt: is there a pattern?
# n=3: 1/3
# n=4: let's check output above
# n=5: from scipy output above

# Actually the volume of the convex hull of even-parity vertices of [0,1]^n
# has a known formula. The half-cube is the polytope h_n.
#
# For n ≥ 3, the normalized volume of h_n (relative to the unit cube) is:
# Vol(h_n) = 2^{n-1} / n! * something
#
# Actually: the demicube h_n inscribed in [0,1]^n has
# 2^{n-1} vertices, all with even coordinate sums.
# The volume can be computed as:
# Vol = (number of simplices in triangulation) / n!
#
# A cleaner approach: h_n is the orbit polytope of a simplex under S_n.
# For the STANDARD demicube in [-1,1]^n, Vol = 2^n / n!... no.

# Let me just read off the pattern from scipy:
print("Reading pattern from computed volumes:")
vols = {}
try:
    from scipy.spatial import ConvexHull
    for n in range(3, 9):
        verts = [v for v in product([0, 1], repeat=n) if sum(v) % 2 == 0]
        pts = np.array(verts, dtype=float)
        hull = ConvexHull(pts)
        vols[n] = hull.volume
        frac = Fraction(hull.volume).limit_denominator(100000)
        print(f"  n={n}: Vol = {hull.volume:.10f} = {frac}")

    # Check: is Vol(h_n) = Σ_{k even, 0≤k≤n} C(n,k) / n! × something?
    # Or: Vol(h_n) = (1/2)(1 - 1/n!) ... let me check
    # n=3: 1/3. Is 1/3 = (1/2)(1 - 1/3)? = (1/2)(2/3) = 1/3. YES!
    # n=4: 2/3. Is 2/3 = (1/2)(1 + 1/3)? = (1/2)(4/3) = 2/3. Hmm, different formula.
    # Let me try: Vol = 1/2 * (1 + (-1)^n / n!)?
    # n=3: 1/2 * (1 - 1/6) = 5/12 ≠ 1/3.

    # Try: Vol = 1/2 - (-1)^n / (2·n!)?
    # n=3: 1/2 - (-1)^3/(2·6) = 1/2 + 1/12 = 7/12 ≠ 1/3.

    # Direct approach: compute Vol exactly
    # n=3: 1/3
    # n=4: 2/3
    # n=5: 13/15
    # n=6: 23/45 × 2 = ... let me check

    # Actually, the formula might be:
    # Vol(h_n) = (1/n!) Σ_{k=0}^{⌊n/2⌋} C(n,2k) * (n-2k)^n / n^n ... no too complex

    # Let me try inclusion-exclusion on the half-cube.
    # The half-cube = conv(even vertices) = intersection of half-spaces
    # through faces of the cube that "cut off" odd vertices.

    # Actually, there's a simpler formula.
    # The n-demicube is a RECTIFIED n-simplex for even n,
    # and an ALTERNATED n-cube for all n.

    # Volume of alternated cube [0,1]^n:
    # V(h_n) = Σ_{σ ∈ A_n} 1/n! = |A_n|/n! = 1/2
    # This would give 1/2 for all n, which contradicts n=3 giving 1/3.

    # Wait — the alternating group A_n has |A_n| = n!/2 elements.
    # Each Coxeter simplex has volume 1/n!.
    # If the demicube is the union of the A_n Coxeter simplices:
    # Vol = (n!/2) · (1/n!) = 1/2.
    # But we computed Vol(h_3) = 1/3, not 1/2!

    # So the demicube is NOT the union of A_n Coxeter simplices.
    # The Coxeter simplices have vertices at 0 and at all x_{σ(1)}=...=x_{σ(k)}
    # while the demicube has vertices at {0,1}^n with even parity.

    # These are DIFFERENT decompositions.

    # Let me try a recursive formula.
    # h_n = conv(v ∈ {0,1}^n : Σv_i ≡ 0 mod 2)
    # Split by first coordinate:
    # x_1=0: remaining (x_2,...,x_n) must have even sum → h_{n-1}^even in {0,1}^{n-1}
    # x_1=1: remaining must have odd sum → h_{n-1}^odd in {0,1}^{n-1}
    # But h_{n-1}^even = h_{n-1} and h_{n-1}^odd = complement.
    # The convex hull of h_{n-1} at x_1=0 and complement(h_{n-1}) at x_1=1
    # is more complex than just stacking.

    # Let me try a different approach: compute by subtracting corners.
    # Cube = demicube + 2^{n-1} corner pieces
    # Each corner piece at an odd vertex v is the simplex conv(v, neighbors of v in demicube)
    # ... but this isn't right either.

    # Let me just identify the pattern numerically.
    print("\n  Pattern search:")
    for n in range(3, 9):
        v = vols[n]
        print(f"  n={n}: Vol = {Fraction(v).limit_denominator(1000)}, "
              f"1-Vol = {Fraction(1-v).limit_denominator(1000)}, "
              f"Vol*n! = {Fraction(v * math.factorial(n)).limit_denominator(1000)}")

    # Vol*n!:
    # n=3: 1/3 * 6 = 2
    # n=4: 2/3 * 24 = 16
    # n=5: ? * 120 = ?
    # The sequence Vol*n! might be recognizable.

    print("\n  Vol(h_n) * n!:")
    for n in range(3, 9):
        val = vols[n] * math.factorial(n)
        print(f"  n={n}: {val:.6f} ≈ {round(val)}")

    # If Vol*n! = {2, 16, 104, ...}
    # 2 = 2^1
    # 16 = 2^4
    # 104 = 8 × 13
    # Hmm, let me check if it's 2^{n-1} - something or Σ_even C(n,k)

    # Σ_{k even} C(n,k) = 2^{n-1}
    # So Vol * n! = 2^{n-1}?
    # n=3: 2^2 = 4 ≠ 2. No.

    # How about Vol = Σ_{σ ∈ S_n, σ even and keeps demicube order} 1/n!
    # = (# of permutations whose Coxeter simplex is in demicube) / n!

    # The demicube triangulates into SOME number of simplices.
    # At n=3: Vol=1/3 = 2/6, so 2 Coxeter simplices out of 6.
    # At n=4: Vol=2/3 = 16/24, so 16 Coxeter simplices out of 24.
    # Sequence of simplex counts: 2, 16, ?, ?, ...

    # n=3: 2 simplices. These are the even permutations?
    # A_3 has 3 elements: (), (123), (132). That's 3, not 2.
    # D_3 = A_3 so... hmm.

    # Actually at n=3, the 2 simplices might correspond to:
    # The tetrahedron ABCD can be split into 2 orthoschemes? No, a tet splits into... hmm.
    # Vol(tet) = 1/3, each Coxeter simplex = 1/6. So 2 Coxeter simplices fit. ✓

    # At n=4: 16 Coxeter simplices. A_4 has 12 elements, not 16.
    # The D_n Coxeter group? |D_n| = 2^{n-1} · n!
    # D_3: 2^2 · 6 = 24, D_4: 2^3 · 24 = 192.
    # Not matching.

    # The number of Coxeter simplices in the demicube:
    # n=3: 2
    # n=4: 16
    # n=5: Vol*120 ≈ ?

except ImportError:
    pass

# ============================================================
# Part 3: The USER'S pattern — simplex halves
# ============================================================
print("\n" + "=" * 70)
print("PART 3: THE USER'S PATTERN — SIMPLEX HALVES")
print("=" * 70)

print("""
The user's observation:
  n=2: triangle in square, 2 halves on either side
  n=3: tetrahedron in cube, 4 halves around it

Interpretation: each "half" is a corner piece with volume = (1/2) × simplex volume.

n=3: simplex vol = 1/3, corner vol = 1/6 = (1/2)(1/3). ✓
     4 corners × 1/6 = 2/3 complement. ✓

n=2: For the equilateral triangle inscribed in square:
     The equilateral triangle has area = (√3/4) × s²
     For it to fit in a unit square with "2 halves on either side":
     Side length s = 1, area = √3/4 ≈ 0.433.
     2 halves each ≈ (1 - √3/4)/2 ≈ 0.284.
     Ratio = 0.284/0.433 ≈ 0.655 ≠ 0.5.

     That doesn't give "halves" in the volume sense.

     ALTERNATIVE: The user might mean the RIGHT triangle (half of square):
     Area = 1/2. The OTHER triangle is the single complement piece.
     But there's 1 complement, not 2.

     OR: The user means the diagonal splits the square into 2 pieces,
     and the TRIANGLE is one of them. Then there are 2 pieces total
     (2 halves of the square), and the simplex IS one of the halves.
     So "2 halves on either side" = the 2 halves are the triangle
     AND its complement.

     This interpretation: the cube/square is divided into
     the simplex + complement, and the pattern is about the
     COMBINATORICS of the complement.

REINTERPRETATION:
  n=2: Square = 2 right triangles (the Coxeter decomposition: 2! = 2 pieces)
       One piece = the simplex, the other = complement.
       Number of complement pieces = 1 = 2!/1 - 1 = 1.

  n=3: Cube = regular tet + 4 corner tets.
       The 4 corners are the "halves."
       Corner count = 2^{n-1} = 4.

So the pattern 2, 4 represents:
  n=2: 2 = total pieces (simplex + 1 complement)
  n=3: 4 = complement pieces only

These are DIFFERENT counting conventions!

If we use "complement pieces = 2^{n-1}":
  n=2: 2^1 = 2 complement pieces
  n=3: 2^2 = 4 complement pieces ✓

But at n=2 the demicube has vertices (0,0),(1,1) which is a DIAGONAL,
not a triangle. The 2 complement pieces would be the 2 triangles
on either side of the diagonal. THAT matches "2 halves on either side"!
""")

# Verify: the diagonal of the square divides it into 2 right triangles
print("n=2 VERIFIED: diagonal (0,0)-(1,1) divides square into 2 triangles")
print("  Each triangle has area 1/2. These are the '2 halves.'")
print("  The 'equilateral triangle' is actually the DIAGONAL (1-simplex).")
print()
print("n=3 VERIFIED: demicube tetrahedron divides cube into 1 tet + 4 corners")
print("  The 4 corners are the '4 halves.'")
print()

# So the pattern is: demicube creates 2^{n-1} corner pieces
# The demicube itself might not be a simplex for n≥4,
# but the corners always exist.

# ============================================================
# Part 4: n=4 continuation
# ============================================================
print("=" * 70)
print("PART 4: n=4 — 8 'HALVES' (16-CELL IN TESSERACT)")
print("=" * 70)

print("""
n=4: The even-parity vertices of [0,1]^4 form a 16-cell (8 vertices).
     The odd-parity vertices also form a 16-cell (8 vertices).
     The tesseract = 16-cell ∪ 8 corner pieces.

     Number of 'halves' = 2^{n-1} = 8.

     But what ARE the corner pieces?
     Each odd vertex is "cut off" by the 16-cell.
     The corner at an odd vertex v is a simplex formed by
     v and its n = 4 neighbors in the 16-cell.
""")

try:
    from scipy.spatial import ConvexHull

    even_4 = [v for v in product([0, 1], repeat=4) if sum(v) % 2 == 0]
    odd_4 = [v for v in product([0, 1], repeat=4) if sum(v) % 2 == 1]

    pts_e = np.array(even_4, dtype=float)
    hull_e = ConvexHull(pts_e)
    vol_e = hull_e.volume

    print(f"  16-cell volume: {vol_e:.10f}")
    comp = 1.0 - vol_e
    print(f"  Complement: {comp:.10f}")
    print(f"  8 corners, each: {comp/8:.10f}")
    print(f"  = {Fraction(comp/8).limit_denominator(1000)}")

    # Each corner is a 4-simplex?
    # At an odd vertex like (1,0,0,0):
    # Its even neighbors (differ in 1 coord): need even sum after flip.
    # (1,0,0,0) has sum 1 (odd). Flipping one coord:
    # (0,0,0,0) sum 0 ✓, (1,1,0,0) sum 2 ✓, (1,0,1,0) sum 2 ✓, (1,0,0,1) sum 2 ✓
    # So 4 even neighbors.
    # The corner at (1,0,0,0) is the simplex with vertices:
    # (1,0,0,0), (0,0,0,0), (1,1,0,0), (1,0,1,0), (1,0,0,1)

    corner_verts = np.array([
        [1,0,0,0], [0,0,0,0], [1,1,0,0], [1,0,1,0], [1,0,0,1]
    ], dtype=float)

    # Volume of 4-simplex
    mat = np.array([corner_verts[i] - corner_verts[0] for i in range(1, 5)])
    det_val = abs(np.linalg.det(mat))
    vol_corner = det_val / math.factorial(4)
    print(f"\n  Corner simplex at (1,0,0,0):")
    print(f"  Volume = {vol_corner:.10f} = {Fraction(vol_corner).limit_denominator(1000)}")
    print(f"  8 × corner = {8*vol_corner:.10f} = complement ✓" if abs(8*vol_corner - comp) < 1e-8 else f"  8 × corner = {8*vol_corner:.10f} ≠ complement")

    # Ratio of corner to 16-cell
    print(f"  Corner/16-cell = {vol_corner/vol_e:.10f}")
    print(f"  = {Fraction(vol_corner/vol_e).limit_denominator(1000)}")

except ImportError:
    pass

# ============================================================
# Part 5: General n pattern
# ============================================================
print("\n" + "=" * 70)
print("PART 5: THE GENERAL PATTERN")
print("=" * 70)

print("""
DEMICUBE PACKING PATTERN:

n | Demicube      | Vol(DC) | Corners | Vol/corner | Corner/DC
--|---------------|---------|---------|------------|----------
2 | line segment  |    0    |    2    |    1/2     |    ∞
3 | tetrahedron   |   1/3   |    4    |    1/6     |   1/2
4 | 16-cell       |   2/3   |    8    |    1/24    |   1/16
5 | 5-demicube    |  13/15  |   16    |   1/120    |   ?
n | n-demicube    |   ?     |  2^{n-1}|   1/n!     |   ?
""")

# Compute corner volumes for higher n
try:
    from scipy.spatial import ConvexHull

    print("Computed values:")
    for n in range(3, 8):
        even_n = [v for v in product([0, 1], repeat=n) if sum(v) % 2 == 0]
        pts = np.array(even_n, dtype=float)
        hull = ConvexHull(pts)
        vol_dc = hull.volume

        n_corners = 2**(n-1)
        vol_comp = 1.0 - vol_dc
        vol_corner = vol_comp / n_corners

        # Check: is corner volume 1/n! ?
        expected = 1.0 / math.factorial(n)

        print(f"  n={n}: Vol(DC) = {Fraction(vol_dc).limit_denominator(10000)}, "
              f"corners = {n_corners}, "
              f"vol/corner = {Fraction(vol_corner).limit_denominator(100000)}, "
              f"1/n! = {Fraction(expected).limit_denominator(100000)}, "
              f"match = {abs(vol_corner - expected) < 1e-8}")

except ImportError:
    pass

# ============================================================
# Part 6: The corner volume = 1/n! theorem
# ============================================================
print("\n" + "=" * 70)
print("PART 6: EACH CORNER = 1/n! (PROOF)")
print("=" * 70)

print("""
THEOREM: Each of the 2^{n-1} corner pieces of the n-demicube in [0,1]^n
has volume exactly 1/n!.

PROOF: Consider the corner at odd vertex v = (v_1,...,v_n) with Σv_i odd.
The n even neighbors of v are obtained by flipping one coordinate.
Let w_j = v with v_j flipped (j = 1,...,n).

The corner simplex has vertices v, w_1, ..., w_n.
Its edge vectors are w_j - v = ±e_j (the j-th standard basis vector).

The volume is |det(w_1-v, ..., w_n-v)| / n! = |det(±I)| / n! = 1/n!.

∎

COROLLARY:
  Vol(complement) = 2^{n-1} / n!
  Vol(demicube) = 1 - 2^{n-1}/n!

CHECK:
  n=3: 1 - 4/6 = 1 - 2/3 = 1/3 ✓
  n=4: 1 - 8/24 = 1 - 1/3 = 2/3 ✓
  n=5: 1 - 16/120 = 1 - 2/15 = 13/15 ✓
  n=6: 1 - 32/720 = 1 - 2/45 = 43/45 ✓

This is beautiful!

The RATIO corner/demicube:
  = (1/n!) / (1 - 2^{n-1}/n!)
  = 1 / (n! - 2^{n-1})

For the user's "halves" (corner = half of simplex):
  At n=3: corner/tet = (1/6)/(1/3) = 1/2 ✓ (each corner is HALF the tet)
  At n=4: corner/16cell = (1/24)/(2/3) = 1/16
  At n=5: corner/demicube = (1/120)/(13/15) = 15/1560 = 1/104

So the "halves" interpretation (corner = half of central polytope)
is SPECIFIC to n=3!

At n=3: corner/central = 1/2 → TRUE halves
At n=4: corner/central = 1/16 → sixteenths, not halves
At n=5: corner/central = 1/104 → much smaller

The ratio shrinks rapidly because n! grows much faster than 2^n.
""")

# Print the table
print("\nFull table:")
print(f"{'n':>3} | {'Vol(DC)':>12} | {'corners':>7} | {'vol/corner':>12} | {'corner/DC':>12} | {'2^(n-1)/n!':>12}")
for n in range(2, 12):
    dc_vol = Fraction(math.factorial(n) - 2**(n-1), math.factorial(n))
    corner_vol = Fraction(1, math.factorial(n))
    n_corners = 2**(n-1)
    if dc_vol > 0:
        ratio = corner_vol / dc_vol
    else:
        ratio = Fraction(0)
    complement = Fraction(2**(n-1), math.factorial(n))

    print(f"{n:3d} | {str(dc_vol):>12} | {n_corners:>7} | {str(corner_vol):>12} | "
          f"{str(ratio):>12} | {str(complement):>12}")

# ============================================================
# Part 7: Connection to tournament evaluation
# ============================================================
print("\n" + "=" * 70)
print("PART 7: TOURNAMENT EVALUATION INTERPRETATION")
print("=" * 70)

print("""
At the tournament point x=2:
  (x+1)^n = 3^n (simplex)
  (x+2)^n = 4^n (cuboid)

The DEMICUBE evaluation:
  Each corner has "tournament volume" 4^n / n! × 1/n! ... no.

  Actually: the polynomials (x+1)^n and (x+2)^n are NOT the same as
  geometric volumes. They count COMBINATORIAL structures.

  The geometric packing gives:
    Vol(demicube) = 1 - 2^{n-1}/n!
    Vol(corner) = 1/n!

  At x=2 evaluation:
    "Simplex value" = 3^n
    "Cuboid value" = 4^n
    "Complement value" = 4^n - 3^n

  The geometric ratio: 2^{n-1}/n! (complement fraction)
  The polynomial ratio: (4^n - 3^n)/4^n = 1 - (3/4)^n

  These are DIFFERENT! The geometric complement shrinks as 2^n/n! → 0,
  while the polynomial complement grows as 1 - (3/4)^n → 1.

  This INVERSION is key: geometrically, the demicube dominates (fills cube).
  Combinatorially, the complement dominates (most structures are non-simplex).
""")

# Compute both ratios
print("Geometric vs Polynomial complement fractions:")
print(f"{'n':>3} | {'Geom: 2^(n-1)/n!':>18} | {'Poly: 1-(3/4)^n':>18}")
for n in range(2, 12):
    geom = 2**(n-1) / math.factorial(n)
    poly = 1 - (3/4)**n
    print(f"{n:3d} | {geom:18.10f} | {poly:18.10f}")

print("""
BEAUTIFUL INVERSION:
  Geometrically: the demicube takes over (Vol → 1 as n → ∞)
  Combinatorially: the complement takes over (→ 1 as n → ∞)

  At n=3: Geom complement = 2/3, Poly complement = 37/64 ≈ 0.578
  These are closest at small n and diverge dramatically.

  INTERPRETATION: In low dimensions, geometry and combinatorics agree.
  In high dimensions, the exponential (4^n vs 3^n) overwhelms
  the factorial (1/n!). This is why tournament theory is "hard" —
  the combinatorial landscape is exponentially larger than the
  geometric intuition suggests.
""")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
