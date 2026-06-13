"""
simplex_cuboid_eisenstein.py -- opus-2026-03-14-S71l
SIMPLEX-IN-CUBOID NESTING: Continuing the user's pattern and connecting to Eisenstein

The user asks:
- "equilateral triangle sits in a square with two halves on either side"
- "tetrahedron sits in a cube with 4 halves around it"
- "how does this continue as n increases?"
- "think of simplices as (x+1)^n and cuboids as (x+2)^n"

KEY CONNECTIONS to establish:
1. The nesting decomposition: (x+2)^n = sum C(n,k)(x+1)^k
2. Coxeter decomposition: n-cube = n! simplices (one per permutation)
3. Corner count = 2^n - (n+1) = nonlinear Walsh monomials
4. Hadamard simplex exists at n+1 = Hadamard order: n = 1, 3, 7, ...
5. At n=7: the Fano plane IS the Hadamard(8) construction
"""

import sys
import numpy as np
from math import factorial, comb, sqrt, log2
from itertools import combinations

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("SIMPLEX-IN-CUBOID: HOW DOES THE NESTING CONTINUE?")
    print("opus-2026-03-14-S71l")
    print("=" * 70)

    # Part 1: The user's examples
    print(f"\n{'='*70}")
    print("PART 1: THE USER'S EXAMPLES")
    print(f"{'='*70}")

    print(f"""
  n=2: Equilateral triangle in a square.
    The triangle occupies the center.
    2 "halves" (triangular wedges) remain on either side.
    complement_pieces = 2

  n=3: Regular tetrahedron in a cube.
    Inscribe using alternating cube vertices:
      (0,0,0), (1,1,0), (1,0,1), (0,1,1)  [even parity]
    Volume of tetrahedron = 1/3 of cube.
    4 right-tetrahedra corners remain.
    complement_pieces = 4

  PATTERN: 2, 4, ... = 2^1, 2^2, ...
  Hypothesis: complement_pieces(n) = 2^{{n-1}}.

  BUT WAIT. The n=3 tetrahedron in cube is special:
  The 8 cube vertices split into 2 groups of 4:
    Even parity: (0,0,0), (1,1,0), (1,0,1), (0,1,1) -> tetrahedron 1
    Odd parity:  (1,0,0), (0,1,0), (0,0,1), (1,1,1) -> tetrahedron 2

  The TWO tetrahedra together fill the cube!
  Each has volume 1/3, and the 4 "corner pieces" are the
  INTERSECTION regions between the two dual tetrahedra.

  Actually: Vol(tet1) + Vol(tet2) = 1/3 + 1/3 = 2/3 < 1.
  The remaining 1/3 is the 4 corner octahedra? No.

  Let me recompute. Regular tetrahedron in unit cube:
  Vertices: A=(0,0,0), B=(1,1,0), C=(1,0,1), D=(0,1,1).
  Edge length: |AB| = sqrt(2), all edges = sqrt(2).
  Volume = (sqrt(2))^3 / (6*sqrt(2)) = 2*sqrt(2)/(6*sqrt(2)) = 1/3.

  The complement = 1 - 1/3 = 2/3.
  This decomposes into 4 equal right tetrahedra:
    Corner at (1,0,0): right tet with legs along 3 edges of cube
    Corner at (0,1,0): same
    Corner at (0,0,1): same
    Corner at (1,1,1): same
  Each has volume = 1*(1/2)*(1)/3 = 1/6.
  Check: 4 * 1/6 = 2/3. YES!

  So: 1 central piece + 4 corner pieces = 5 total pieces.
  Or: 1 simplex + 2^n - (n+1) corners (at n=3: 1 + 4 = 5)
  Wait, 2^3 - (3+1) = 4. Yes!

  The formula 2^n - (n+1) gives:
    n=1: 2-2 = 0 corners (line segment = simplex, no leftovers)
    n=2: 4-3 = 1 corner (?? user said 2)
    n=3: 8-4 = 4 corners (matches)
    n=4: 16-5 = 11 corners
    n=5: 32-6 = 26 corners
    n=7: 128-8 = 120 corners

  Hmm, at n=2 this gives 1 corner, but user said 2.
  The discrepancy: for n=2, you can't inscribe a regular
  simplex (equilateral triangle) in a square using VERTICES.
  No 3 vertices of a square form an equilateral triangle.
  The user is using a DIFFERENT construction for n=2.
""")

    # Part 2: The two constructions
    print(f"\n{'='*70}")
    print("PART 2: TWO CONSTRUCTIONS — VERTEX-INSCRIBED vs CENTERED")
    print(f"{'='*70}")

    print(f"""
  CONSTRUCTION A (vertex-inscribed):
    Use vertices of the cube to form a regular simplex.
    Only works when n+1 is a HADAMARD ORDER (1, 2, 4, 8, ...).
    So n = 0, 1, 3, 7, 15, 31, ...

  CONSTRUCTION B (centered/general):
    Place a regular n-simplex CENTERED in the n-cube,
    not necessarily using cube vertices.
    Works for all n.

  The USER'S CONSTRUCTION (for n=2):
    An equilateral triangle placed centrally in a square.
    The triangle touches 3 sides of the square, dividing
    the complement into 2 roughly equal pieces.

  For TOURNAMENT THEORY, Construction A is more relevant:
    n=3: 4 vertices -> tetrahedron, 4 corner pieces
    n=7: 8 vertices -> 7-simplex, 120 corner pieces

  The CORNER COUNT for Construction A:
    Each non-simplex vertex of the cube generates one corner piece.
    Count = 2^n - (n+1) (when n+1 is a power of 2).

  BEAUTIFULLY, this equals the number of NONLINEAR Walsh monomials!
  Linear Walsh: degree 0 (one: the constant) + degree 1 (n terms) = n+1
  Nonlinear: degree >= 2 = 2^n - (n+1)

  The Walsh decomposition of H corresponds to:
    Linear part = "simplex" = mean + individual arc effects
    Nonlinear part = "corners" = interaction effects
""")

    # Part 3: Volume ratios
    print(f"\n{'='*70}")
    print("PART 3: VOLUME RATIOS — SIMPLEX/CUBE")
    print(f"{'='*70}")

    print(f"  For vertex-inscribed regular simplex (Hadamard construction):")
    for n in [1, 3, 7, 15]:
        # Volume of regular n-simplex with edge sqrt(2) * sqrt(n/2)...
        # Actually edge = sqrt((n+1)/2) for Hadamard construction
        # Regular simplex volume: V = edge^n * sqrt(n+1) / (n! * 2^{n/2})
        # Edge from Hadamard: all pairs at Hamming distance (n+1)/2
        edge = sqrt((n+1)/2)
        V_simplex = edge**n * sqrt(n+1) / (factorial(n) * 2**(n/2))
        V_cube = 1.0
        corners = 2**n - (n+1)
        V_corners = V_cube - V_simplex
        V_each = V_corners / corners if corners > 0 else 0

        print(f"  n={n:2d}: edge={edge:.4f}, V_simplex={V_simplex:.6f}, "
              f"corners={corners}, V_each_corner={V_each:.6f}")
        print(f"        ratio simplex/cube = 1/{1/V_simplex:.1f}, "
              f"corner fraction = {V_corners:.6f}")

    # Part 4: The (x+1)^n inside (x+2)^n decomposition
    print(f"\n{'='*70}")
    print("PART 4: ALGEBRAIC NESTING — (x+2)^n = sum C(n,k)(x+1)^k")
    print(f"{'='*70}")

    print(f"""
  The user says: "think of simplices as (x+1)^n and cuboids as (x+2)^n."

  Algebraically: (x+2)^n = ((x+1) + 1)^n = sum_k C(n,k)(x+1)^k

  This is a DECOMPOSITION of the cuboid into simplicial pieces!

  At x=1 (the OCF evaluation point):
    3^n = sum C(n,k) * 2^k

  This is just 3^n = (1+2)^n, but the interpretation is:
  The "3-adic" tournament (Phi_3 structure) decomposes into
  a sum of "2-adic" simplicial pieces.

  The pieces have sizes C(n,k) * 2^k:
""")

    for n in range(2, 8):
        total = 3**n
        pieces = []
        for k in range(n+1):
            piece = comb(n, k) * 2**k
            pieces.append(piece)
        print(f"  n={n}: 3^{n} = {total} = ", end="")
        print(" + ".join(f"C({n},{k})*2^{k}={p}" for k, p in enumerate(pieces)))

    # Part 5: The Phi_3 connection
    print(f"\n{'='*70}")
    print("PART 5: Phi_3 AND THE SIMPLEX-CUBOID RATIO")
    print(f"{'='*70}")

    print(f"""
  (x+2)^n / (x+1)^n = ((x+2)/(x+1))^n = (1 + 1/(x+1))^n

  At x = 1: 3^n / 2^n = (3/2)^n

  The ratio cuboid/simplex grows as (3/2)^n.
  This is the FUNDAMENTAL RATIO of tournament theory:
  Mean(H) = n!/2^{{n-1}} and the total tournament count is 2^{{C(n,2)}}.

  NOW: what about Phi_3?
  Phi_3(x) = x^2 + x + 1 = (x+1)^2 + (x+1) + 1 ... no, that's not right.
  (x+1)^2 + (x+1) + 1 = x^2 + 2x + 1 + x + 1 + 1 = x^2 + 3x + 3 = Phi_3(x+1)? No.
  x^2 + 3x + 3 evaluated at x=1: 7 = Phi_3(2). YES!

  So Phi_3(x+1) = (x+1)^2 + (x+1) + 1 = x^2 + 3x + 3.

  The forbidden value 7 = Phi_3(2) = Phi_3(1+1) = 1 + 3 + 3 = 7.
  And 21 = Phi_3(4) = Phi_3(3+1) = 9 + 9 + 3 = 21.

  WAIT: Phi_3(x+1) = x^2 + 3x + 3.
  This factors in Z_3: x^2 + 3x + 3 equiv x^2 mod 3.
  So Phi_3(x+1) equiv x^2 mod 3.
  And Phi_3(2) = Phi_3(1+1) = 1 + 3 + 3 = 7 equiv 1 mod 3.

  The simplex-cuboid nesting gives:
  (x+2)^n - (x+1)^n = "corners" = sum_{{k=0}}^{{n-2}} C(n,k)(x+1)^k + (n-1)(x+1)^{{n-1}}
  Wait, let me just expand:
  (x+2)^n - (x+1)^n = sum_{{k=0}}^n C(n,k)(x+1)^k - (x+1)^n
                     = sum_{{k=0}}^{{n-1}} C(n,k)(x+1)^k

  At x=1: 3^n - 2^n = sum_{{k=0}}^{{n-1}} C(n,k) * 2^k = 3^n - 2^n.

  Trivially true, but the POINT is: the "corners" are the
  LOWER-DEGREE simplicial pieces.
""")

    for n in range(2, 8):
        corners = 3**n - 2**n
        simplex = 2**n
        ratio = corners / simplex
        print(f"  n={n}: simplex=2^{n}={simplex}, corners=3^{n}-2^{n}={corners}, "
              f"ratio corners/simplex = {ratio:.4f}")

    # Part 6: The continuation pattern (direct answer)
    print(f"\n{'='*70}")
    print("PART 6: THE CONTINUATION — DIRECT ANSWER")
    print(f"{'='*70}")

    print(f"""
  The user's pattern (GEOMETRIC nesting):

  n=2: Triangle in square.
    Complement: 2 pieces.
    Volume fractions: triangle = sqrt(3)/4 = 0.433
                      2 pieces = 0.567 total

  n=3: Tetrahedron in cube.
    Complement: 4 right-tetrahedra corners.
    Volume fractions: tetrahedron = 1/3 = 0.333
                      4 pieces = 2/3 total, each = 1/6

  n=4: Can't do vertex-inscribed regular simplex in tesseract.
    BUT: can do COXETER decomposition: 24 = 4! simplices
    fill the tesseract.

  n=7: 7-simplex in 7-cube (Hadamard construction).
    Complement: 120 = 5! pieces.
    This is the FANO PLANE geometry!

  The continuation:
  n=1: 0 extra pieces
  n=3: 4 = 2^2 extra pieces (tetrahedra)
  n=7: 120 = 5! extra pieces (7-simplices? No, right simplices)
  n=15: 32752 = 2^15 - 16 extra pieces

  FORMULA: At Hadamard dimensions n = 2^k - 1:
    corner_count = 2^n - (n+1) = 2^{{2^k-1}} - 2^k
    simplex_volume / cube_volume = (n+1)^{{(n+1)/2}} / (n! * 2^n)

  The corner count at the first few Hadamard dimensions:
    n=1:  0 corners, 1 simplex fills line
    n=3:  4 corners = 2^2 = C(4,2) = "pairs of vertices"
    n=7:  120 corners = 5! = C(8,3)*..., related to octonion algebra
    n=15: 32752 corners = 2^15 - 16

  BEAUTIFUL: 120 = dimension of the Lie algebra of SO(16).
  And 120 = C(16,2) = "pairs in 16-element set."
  No wait: C(16,2) = 120. And C(8,3) = 56.
  Actually 120 = 5! = number of HPs of transitive T_5 = 1... no.
  120 = 5! = number of permutations of 5 elements.
  Also 120 = number of long roots in E8? No, that's 240.
  120 = order of S_5 = icosahedral group = 120.

  The sequence 0, 4, 120, 32752, ... = 2^{{2^k-1}} - 2^k.
""")

    print(f"  Hadamard simplex corner counts:")
    for k in range(1, 8):
        n = 2**k - 1
        corners = 2**n - (n+1)
        print(f"    k={k}: n=2^{k}-1={n}, corners = 2^{n} - {n+1} = {corners}")
        if corners > 10**15:
            break

    # Part 7: Tournament meaning
    print(f"\n{'='*70}")
    print("PART 7: TOURNAMENT INTERPRETATION")
    print(f"{'='*70}")

    print(f"""
  The simplex-cuboid nesting has a DIRECT tournament meaning:

  The n-cube represents 2^n Walsh monomials (subsets of [n]).
  The simplex represents the n+1 LINEAR monomials (empty + singletons).
  The corners represent 2^n - (n+1) NONLINEAR monomials.

  LINEAR WALSH COMPONENTS:
    hat{{H}}[empty] = mean(H) = n!/2^{{n-1}}  (the "volume" of the simplex)
    hat{{H}}[{{i}}] = degree-1 terms (individual vertex effects)
    These n+1 terms are the "simplex" of H.

  NONLINEAR WALSH COMPONENTS:
    hat{{H}}[S] for |S| >= 2 (interaction effects)
    These 2^n-(n+1) terms are the "corners."

  VARIANCE DECOMPOSITION:
    Var(H) = sum_{{|S|>=1}} hat{{H}}[S]^2
    Level-2 variance = sum_{{|S|=2}} hat{{H}}[S]^2 (C(n,2) terms)
    Higher-level variance = sum_{{|S|>=3}} hat{{H}}[S]^2 (remaining terms)

  The FRACTION of variance in linear vs nonlinear:
    Linear share = C(n,1) * hat{{H}}_1^2 / Var(H)
    Nonlinear share = 1 - linear share

  From THM-204 (Grand Fourier Level Formula):
    E_2/E_0 = (n-2)/C(n,2) = 2/(n-1)(n)... wait.
    E_{{2k}}/E_0 = 2(n-2k)^k / P(n,2k)

  Level 2 (k=1): E_2/E_0 = 2(n-2)/n(n-1) = (n-2)/C(n,2)... no.
  Let me just compute the level-2 share for small n.
""")

    for n in range(3, 9):
        # From THM-204: Var/Mean^2 = sum_{k>=1} 2*(n-2k)^k / P(n,2k)
        # Level-2k contributes 2*(n-2k)^k / P(n,2k) to Var/Mean^2
        total_var_ratio = 0
        level_shares = []
        for k in range(1, n//2 + 1):
            p_falling = 1
            for j in range(2*k):
                p_falling *= (n - j)
            if p_falling == 0:
                break
            contrib = 2 * (n - 2*k)**k / p_falling
            total_var_ratio += contrib
            level_shares.append((2*k, contrib))

        print(f"  n={n}: Var/Mean^2 = {total_var_ratio:.6f}")
        for level, share in level_shares:
            pct = share / total_var_ratio * 100 if total_var_ratio > 0 else 0
            print(f"    Level {level}: {share:.6f} ({pct:.1f}% of variance)")

    print(f"\n{'='*70}")
    print("DONE — SIMPLEX-IN-CUBOID NESTING")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()
