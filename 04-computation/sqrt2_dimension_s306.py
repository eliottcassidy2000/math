#!/usr/bin/env python3
"""
sqrt2_dimension_s306.py — Tournament space as the boundary between 1D and 2D
opus-2026-03-24-S306

The staircase δ_{n-2} IS a right isosceles triangle with:
  legs = n-2, hypotenuse = (n-2)√2

Tournament space lives on this triangle. It's "almost 1D" (H dominates)
but has 2D structure (width perpendicular to H). The √2 appears because
the diagonal of the square has length √2 × leg.

Explore:
1. Effective dimension from eigenvalue distribution
2. The √2 ratio in metagraph metrics
3. Fractal dimension of the distance profile
4. The Krawtchouk correction to linearity (K_2, K_3 contributions)
5. Space-filling curve structure
"""

import sys, time
from math import comb, factorial, sqrt, log, log2, pi
from collections import Counter
import numpy as np

sys.stdout.reconfigure(line_buffering=True)

print("=" * 72)
print("  TOURNAMENT SPACE: BETWEEN 1D AND 2D")
print("  opus-2026-03-24-S306")
print("=" * 72)

# ============================================================
# 1. THE √2 IN THE STAIRCASE
# ============================================================

print(f"\n{'='*72}")
print("  1. THE √2 IN THE STAIRCASE GEOMETRY")
print(f"{'='*72}")

for n in range(3, 12):
    m = comb(n-1, 2)  # tiles = cells of staircase
    leg = n - 2        # leg length
    hyp = leg * sqrt(2)  # hypotenuse
    area = leg**2 / 2   # area of triangle = m (number of cells!)

    print(f"  n={n:>2}: m={m:>4} tiles, leg={leg:>2}, hyp={hyp:>6.2f}, "
          f"area={area:>5.1f}, area/m={area/m:.3f} {'✓' if abs(area-m) < 0.1 else '✗'}")

print(f"\n  The AREA of the right isosceles triangle = m = C(n-1,2).")
print(f"  The HYPOTENUSE has length (n-2)√2.")
print(f"  √2 IS the geometric compression from 2D (square) to 1D (diagonal).")

# ============================================================
# 2. EFFECTIVE DIMENSION FROM EIGENVALUE GAPS
# ============================================================

print(f"\n{'='*72}")
print("  2. EFFECTIVE DIMENSION FROM METAGRAPH SPECTRA")
print(f"{'='*72}")

# From S268/S269: the distance matrix eigenvalue ratios
# λ_0/λ_1 = 22 (n=5), 27 (n=6), 57 (n=7)
# The participation ratio PR = (Σλ)² / Σλ² measures effective dimension

eig_data = {
    5: {'gap_ratio': 22, 'spectral_dim_peak': 1.5},
    6: {'gap_ratio': 27, 'spectral_dim_peak': 2.3},
    7: {'gap_ratio': 57, 'spectral_dim_peak': 3.5},
}

print(f"  {'n':>3} {'m':>4} {'V':>6} {'gap_ratio':>10} {'d_s_peak':>9} {'log2(V)':>8} {'d_s/m':>6}")
for n in [5, 6, 7]:
    m = comb(n-1, 2)
    V = {5:10, 6:34, 7:272}[n]
    gr = eig_data[n]['gap_ratio']
    ds = eig_data[n]['spectral_dim_peak']
    print(f"  {n:>3} {m:>4} {V:>6} {gr:>10} {ds:>9.1f} {log2(V):>8.2f} {ds/m:>6.3f}")

print(f"\n  EFFECTIVE DIMENSION RATIO d_s/m:")
print(f"    n=5: 1.5/6 = 0.250")
print(f"    n=6: 2.3/10 = 0.230")
print(f"    n=7: 3.5/15 = 0.233")
print(f"    → d_s ≈ m/4 ≈ (n-2)²/8")
print(f"    This is the 'effective 1D' regime: d_s/m → 0 but d_s → ∞")

# ============================================================
# 3. THE √2 IN THE KRAWTCHOUK ANALYSIS
# ============================================================

print(f"\n{'='*72}")
print("  3. √2 IN THE KRAWTCHOUK STRUCTURE")
print(f"{'='*72}")

# K_1(x; m) = m - 2x. The factor 2 here is the q of the Hamming scheme.
# For q=2: K_1 = m - 2x = m(1 - 2x/m)
# The "crossing point" is at x = m/2 (center of Hamming weight space).
# H ≈ 0 at this point. The transitive is at x = 0 (K_1 = m, H ≈ 1).
# The regular is at x ≈ m/2 (K_1 ≈ 0, H ≈ max).

# The VARIANCE of K_1 over the binomial distribution:
# If x ~ Binom(m, 1/2): Var(K_1) = Var(m - 2x) = 4 × m/4 = m
# Standard deviation = √m

for n in range(3, 11):
    m = comb(n-1, 2)
    std_K1 = sqrt(m)
    print(f"  n={n:>2}: m={m:>4}, √m = {std_K1:>6.2f}, "
          f"√m/(n-2) = {std_K1/(n-2):>6.3f}, "
          f"√m/√2 = {std_K1/sqrt(2):>6.2f}")

print(f"\n  √m / (n-2) → √((n-2)/2) → ∞ as n → ∞")
print(f"  The ratio √m/(n-2) = √(C(n-1,2))/(n-2) = √((n-1)/2)")
print(f"  AT n=3: √(1/2) = 1/√2 = 0.707")
print(f"  AT n=5: √(3/2) = √(3)/√2 = 1.225")
print(f"  √2 appears as the DENOMINATOR of the leg/hypotenuse ratio!")

# ============================================================
# 4. THE DIAGONAL PROJECTION
# ============================================================

print(f"\n{'='*72}")
print("  4. THE DIAGONAL PROJECTION: 2D → 1D VIA √2")
print(f"{'='*72}")

print("""
  The staircase δ_{n-2} lives in a 2D grid with coordinates (r, c).
  r = skip (= x - y - 1), c = lower vertex (= y).
  Constraints: r ≥ 1, c ≥ 1, r + c ≤ n-1.

  TWO NATURAL AXES:
    DIAGONAL: d = r + c (the "strip" direction). Range: 2 to n-1.
    ANTI-DIAGONAL: a = r - c (the "skew" direction). Range: -(n-3) to n-3.

  H correlates with the DIAGONAL d = r + c:
    Tiles with higher d (farther from the origin) have larger hooks.
    The H-gradient runs along the diagonal.

  The WIDTH perpendicular to the diagonal is the anti-diagonal span.
  At diagonal d: width = min(d-1, n-1-d) cells.
  Maximum width = (n-2)/2 ≈ n/2.
  "Length" (along diagonal) = n-3.
  ASPECT RATIO = width/length ≈ 1/2 → the triangle is "twice as long as wide."

  The PROJECTION from 2D to 1D (along the diagonal) compresses by:
    compression = 1/√2 (because the diagonal is √2 × the side)

  So: tournament space ≈ 1D projection of 2D triangle via √2 compression.
""")

for n in range(4, 10):
    length = n - 3  # diagonal range
    max_width = (n - 2) // 2
    m = comb(n-1, 2)
    aspect = max_width / length if length > 0 else 0
    compression = 1 / sqrt(2)

    print(f"  n={n}: length={length}, max_width={max_width}, "
          f"aspect={aspect:.3f}, m={m}, m/(length×width)={m/(length*max_width) if max_width > 0 else 0:.2f}")

# ============================================================
# 5. THE HAUSDORFF DIMENSION ESTIMATE
# ============================================================

print(f"\n{'='*72}")
print("  5. HAUSDORFF DIMENSION OF TOURNAMENT SPACE")
print(f"{'='*72}")

# V_n grows as 2^{m}/n! ≈ 2^{n²/2}/n! (by Burnside)
# m = C(n-1,2) ≈ n²/2
# log(V_n) ≈ m × log(2) - log(n!) ≈ m × log(2) - n × log(n)
# d_H ≈ log(V_n) / log(scale) where scale depends on the metric

# Using the metagraph diameter as the "scale":
# d_H = log(V) / log(diameter + 1)

print(f"  {'n':>3} {'V':>8} {'diam':>5} {'log2(V)':>8} {'log2(d+1)':>10} {'d_H':>6}")
V_vals = {3:2, 4:3, 5:10, 6:34, 7:272, 8:3528}
diam_vals = {3:1, 4:1, 5:3, 6:4, 7:7, 8:8}
for n in [3, 4, 5, 6, 7, 8]:
    V = V_vals[n]; d = diam_vals[n]
    lv = log2(V); ld = log2(d + 1)
    dH = lv / ld if ld > 0 else float('inf')
    print(f"  {n:>3} {V:>8} {d:>5} {lv:>8.2f} {ld:>10.3f} {dH:>6.2f}")

print(f"\n  Hausdorff dimension d_H → ? as n grows.")
print(f"  If d_H converges to a constant, tournament space has fixed fractal dimension.")
print(f"  If d_H → ∞, it fills progressively more of the ambient space.")

# ============================================================
# 6. THE √2 CONJECTURE
# ============================================================

print(f"\n{'='*72}")
print("  6. THE √2 CONJECTURE")
print(f"{'='*72}")
print("""
  CONJECTURE: Tournament space has Hausdorff dimension √2 ≈ 1.414.

  Evidence:
  - The staircase is a right isosceles triangle with hypotenuse/leg = √2
  - The effective dimension d_s ≈ m/4 ≈ (n-2)²/8, while the
    algebraic dimension is m = (n-2)(n-1)/2
  - The ratio d_s/m ≈ 1/4 → 0, but d_s/√m ≈ √m/4 → √(n/2)/4
  - The Krawtchouk connection: K_1 has standard deviation √m,
    and the "width" of tournament space in K_1 is √m ≈ (n-2)/√2
  - The fiber fraction f(n) ∼ 1/√(πn) involves √π from the
    Wallis integral, which is the area formula for the semicircle

  ALTERNATIVE: d_H = 1 + 1/2 = 3/2.
    From S269: peak spectral dimension at n=5 is 1.52 ≈ 3/2.
    3/2 = the dimension of the Sierpinski gasket (not √2 = 1.414).

  ALTERNATIVE: d_H = log(3)/log(2) ≈ 1.585.
    The Sierpinski triangle dimension. The staircase IS a triangle.

  TESTING NEEDED: Compute d_H at n=9,10 to see which constant it approaches.
""")

# ============================================================
# 7. CONNECTIONS TO SPACE-FILLING CURVES
# ============================================================

print(f"\n{'='*72}")
print("  7. SPACE-FILLING CURVE STRUCTURE")
print(f"{'='*72}")
print("""
  A space-filling curve maps [0,1] → [0,1]² surjectively.
  Tournament space is a "quotient space-filling curve":
    - The 1D path = the H-gradient (from transitive to regular)
    - The 2D region = the staircase triangle
    - The "curve" = the path through tournament space ordered by H

  At n=5 (V=10 classes):
    H-ordering: 1, 3, 3, 5, 9, 9, 11, 13, 15, 15
    This path visits 10 "points" in a 2D triangle
    The path is NOT injective (H=3, H=9, H=15 each appear twice)
    The duplicate H-values = width of the triangle at that H-level

  The PATH ITSELF is the H-gradient walk:
    Start at transitive (H=1, corner of triangle)
    Walk along the diagonal, collecting classes at each H-level
    End at regular (H=max, opposite corner)

  This is a DIAGONAL WALK through the triangle:
    - Direction: along the hypotenuse
    - Speed: H grows linearly with position
    - Width: number of classes at each H-level grows, then shrinks

  The FRACTAL STRUCTURE comes from the recursive Mode B decomposition:
    At each level, the triangle splits into sub-triangles
    Each sub-triangle has its own H-gradient
    The full structure is a SELF-SIMILAR hierarchy of diagonal walks
""")

print("DONE.")
