#!/usr/bin/env python3
"""
cartan_packing_unified.py — opus-2026-03-14-S71e

UNIFIED EXPLORATION: Cartan matrices, simplex-cuboid packing, and α₁≥α₂

Key discoveries to build on:
1. det(A_{n-1}) = n — the Cartan determinant gives the tournament size
2. Keys 2,3 = det(A₁), det(A₂) — the recurrence roots
3. α₁ ≥ α₂ proved for n ≤ 9 via Cauchy-Schwarz (S75)
4. The breakdown at n=10: can we fix it using longer cycles?

CENTRAL INSIGHT: The Cauchy-Schwarz proof only uses 3-cycles.
Longer odd cycles use MORE vertices per cycle, making the clique
structure in CG(T) DENSER. This should extend the bound.
"""

import sys
import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PART 1: CARTAN DETERMINANT SEQUENCE AND SIMPLEX PACKING")
print("=" * 70)

# The type A_{n-1} Cartan matrix has det = n
# This means the root lattice has INDEX n in the weight lattice
# Geometrically: n! orthoschemes tile the cube, each of volume 1/n!
# The permutohedron has n! vertices

print("\n  The Cartan matrix A_{n-1} encodes the type A root system.")
print("  Its determinant is the index [weight lattice : root lattice] = n.")
print()

for n in range(2, 13):
    # Build Cartan matrix of type A_{n-1}
    m = n - 1
    if m == 0:
        det_val = 1
        eigs = [1]
    else:
        A = np.zeros((m, m))
        for i in range(m):
            A[i, i] = 2
            if i > 0:
                A[i, i-1] = -1
            if i < m-1:
                A[i, i+1] = -1
        det_val = round(np.linalg.det(A))
        eigs = sorted(np.linalg.eigvalsh(A))

    # Volume of orthoscheme = 1/n!
    vol = 1.0 / factorial(n)
    # Number of orthoschemes = n!
    n_orth = factorial(n)
    # Simplex (regular, Hadamard) volume ratio to cube
    # For Hadamard dimension n: simplex_vol/cube_vol = sqrt(n+1)/(2^n * n!)... actually complex

    print(f"  n={n:2d}: det(A_{n-1}) = {det_val:4d}, "
          f"n! = {n_orth:8d}, "
          f"vol(orthoscheme) = 1/{n_orth}, "
          f"eigenvalue range [{eigs[0]:.3f}, {eigs[-1]:.3f}]")

print()
print("  KEY: det(A₁) = 2, det(A₂) = 3 — these are the tournament keys!")
print("  det(A_{n-1}) = n means n is fundamentally a 'Cartan index'.")
print()

# The eigenvalues of A_{n-1} are 2 - 2cos(kπ/n) for k=1,...,n-1
# Product = det = n (Chebyshev identity)
# This means: n = ∏_{k=1}^{n-1} (2 - 2cos(kπ/n))

print("  CHEBYSHEV PRODUCT IDENTITY:")
print("  n = ∏_{k=1}^{n-1} (2 - 2cos(kπ/n))")
print()
for n in range(2, 10):
    product = 1.0
    for k in range(1, n):
        product *= (2 - 2 * np.cos(k * np.pi / n))
    print(f"  n={n}: ∏ = {product:.6f} (should be {n})")

print()
print("  This connects to cyclotomic polynomials:")
print("  Φ_n(2) = ∏_{gcd(k,n)=1} (2 - 2cos(2πk/n))")
print()

# Compute Φ_n(2) for small n
def cyclotomic_at_2(n):
    """Compute Φ_n(2) using the definition."""
    from math import gcd
    result = 1.0
    for k in range(1, n+1):
        if gcd(k, n) == 1:
            result *= (2 - np.exp(2j * np.pi * k / n))
    return round(result.real)

print("  Φ_n(2) values:")
for n in range(1, 16):
    phi_val = cyclotomic_at_2(n)
    print(f"    Φ_{n:2d}(2) = {phi_val}")

print()

print("=" * 70)
print("PART 2: EXTENDING α₁ ≥ α₂ BEYOND n=9")
print("=" * 70)

# The S75 proof:
# Each 3-cycle uses 3 vertices → Σ d_v = 3α₁
# Cauchy-Schwarz: Σ C(d_v,2) ≥ n·C(3α₁/n, 2) = (9α₁²-3nα₁)/(2n)
# Need this ≥ α₁(α₁-3)/2 → 9 ≥ n
#
# KEY EXTENSION: include ALL odd cycles, not just 3-cycles
# A cycle of length ℓ uses ℓ vertices
# If we have cycles with lengths ℓ₁,...,ℓ_m, then Σ d_v = Σ ℓᵢ ≥ 3m = 3α₁
# Average cycle length L = (Σ ℓᵢ)/α₁ ≥ 3
# Cauchy-Schwarz: Σ C(d_v,2) ≥ n·C(Lα₁/n, 2) = (L²α₁²-Lnα₁)/(2n)
# Need: (L²α₁-Ln)/n ≥ α₁-3
# i.e., L²α₁ - Ln ≥ nα₁ - 3n
# i.e., (L²-n)α₁ ≥ n(L-3)
# For L=3: (9-n)α₁ ≥ 0 → n ≤ 9
# For L=5: (25-n)α₁ ≥ 2n → n ≤ 25 (if α₁ ≥ 2n/15 ≈ n/7.5)
# For L=7: (49-n)α₁ ≥ 4n → n ≤ 49 (if α₁ ≥ 4n/42 ≈ n/10.5)

print("\n  S75 proof bounds (3-cycles only): works for n ≤ 9")
print()
print("  EXTENSION: if average cycle length L > 3, bound extends!")
print("  Need: (L²-n)α₁ ≥ n(L-3) [for α₁ > 0]")
print()

for L in [3, 5, 7, 9, 11]:
    # Maximum n for which the bound works unconditionally
    # (L²-n)α₁ ≥ n(L-3) → true for all α₁ > 0 iff L²-n ≥ 0 AND n(L-3) ≤ 0
    # The second is only true for L=3. For L>3, need α₁ ≥ n(L-3)/(L²-n)
    max_n_unconditional = L**2
    print(f"  L={L}: unconditional for n ≤ {max_n_unconditional}")
    if L > 3:
        for n in [10, 15, 20, 25, 30, 50]:
            if n < L**2:
                min_a1 = n * (L - 3) / (L**2 - n)
                print(f"    n={n:2d}: need α₁ ≥ {min_a1:.1f}")

print()

# But what is the average cycle length at large n?
# The number of 3-cycles grows as O(n³), 5-cycles as O(n⁵), etc.
# For random tournaments:
# E[dc_ℓ] = C(n,ℓ) · ℓ! / (2·ℓ) = C(n,ℓ) · (ℓ-1)! / 2

print("  Expected cycle counts in random tournaments:")
for n in [7, 9, 11, 15, 20]:
    total_cycles = 0
    total_length = 0
    print(f"\n  n={n}:")
    for ell in range(3, min(n+1, 16), 2):  # odd cycles only
        # E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / 2^ℓ
        # Wait, more precisely:
        # Number of ℓ-cycles = (1/ℓ) · Σ_{ℓ-subsets} (number of directed cycles on subset)
        # For random tournament on ℓ vertices: prob of specific directed ℓ-cycle = 1/2^ℓ
        # Number of directed ℓ-cycles on ℓ labeled vertices = (ℓ-1)!
        # E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / (2^ℓ)  ... no
        # Actually: E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / 2^ℓ   (not dividing by ℓ again)
        # Wait: a directed ℓ-cycle on vertices {v_1,...,v_ℓ} is a cyclic permutation.
        # There are (ℓ-1)! such cycles. Each has probability 1/2^ℓ of all arcs matching.
        # But we're double-counting because each undirected cycle corresponds to 2 directed ones.
        # So: E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / (2 · 2^ℓ) ... no
        # Let me be more careful.
        # An ℓ-cycle in the UNDIRECTED sense is a set {v1,...,vℓ} plus a cyclic order.
        # Number of distinct cyclic orders on ℓ elements = (ℓ-1)!/2 (unordered cycles)
        # Each has 2 directed orientations, probability 1/2^(C(ℓ,2)) for a specific orientation.
        # Wait: for a TOURNAMENT, the arcs within the ℓ-subset are determined.
        # An ℓ-cycle exists on {v1,...,vℓ} if the tournament on that subset is a "cyclic tournament".
        # Actually: dc_ℓ = number of distinct directed ℓ-cycles / ℓ
        # (dividing by ℓ for the ℓ rotations of the same cycle).
        # Prob that a specific ordered sequence v1→v2→...→vℓ→v1 is a directed cycle = (1/2)^ℓ
        # Number of ordered cycles on ℓ labeled vertices from n = P(n,ℓ) / ℓ
        # Wait, I'm getting confused. Let me just use the formula:
        # E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / 2^ℓ  (ℓ-cycles on ℓ vertices, directed, div by ℓ rotations)
        # No: the number of directed ℓ-cycles on ℓ vertices is (ℓ-1)! (fix first vertex, permute rest cyclically)
        # Each uses ℓ arcs, each arc has prob 1/2 of matching direction
        # E[directed ℓ-cycles on specific ℓ-subset] = (ℓ-1)! / 2^ℓ
        # Sum over C(n,ℓ) subsets, but each cycle counted ℓ times (once per starting vertex)?
        # No: (ℓ-1)! already accounts for fixing one vertex.
        # So E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / 2^ℓ
        # But this counts both CW and CCW... for directed: both are distinct cycles.
        # Hmm, actually for DIRECTED cycles, rotation equivalence: (v1→v2→v3) = (v2→v3→v1) = (v3→v1→v2)
        # But (v1→v3→v2) is a DIFFERENT directed cycle.
        # So on ℓ vertices: (ℓ-1)! directed cycles / ... no.
        # On ℓ labeled vertices, the directed cycles are cyclic permutations.
        # Total cyclic permutations = (ℓ-1)! (fixing first element).
        # These include both orientations.
        # For UNDIRECTED cycles: (ℓ-1)!/2.
        # For directed: (ℓ-1)!.
        # Each directed cycle uses ℓ specific arcs, each correct with prob 1/2.
        # E[# directed ℓ-cycles on subset] = (ℓ-1)! · (1/2)^ℓ
        # E[dc_ℓ (directed)] = C(n,ℓ) · (ℓ-1)! · (1/2)^ℓ
        # But we count each undirected cycle TWICE (CW and CCW).
        # The standard dc_ℓ counts UNDIRECTED: dc_ℓ = C(n,ℓ)·(ℓ-1)!/(2·2^ℓ)
        # Hmm, but for 3-cycles at n=7:
        # C(7,3)·2!/(2·2³) = 35·2/16 = 70/16 = 4.375
        # We know E[dc₃] at n=7 for random tournament = C(7,3)·(1/4) = 35/4 = 8.75
        # So the formula is: dc₃ = C(n,3)/4 (for random)
        # = C(n,3) · 1/2² = 35/4
        # More generally: E[dc_ℓ] = C(n,ℓ) · (ℓ-1)! / 2^ℓ (directed cycles)
        # For ℓ=3: C(n,3)·2/8 = C(n,3)/4. ✓
        # For ℓ=5: C(n,5)·24/32 = C(n,5)·3/4.
        expected = comb(n, ell) * factorial(ell - 1) / (2 ** ell)
        total_cycles += expected
        total_length += expected * ell
        if expected > 0.01:
            print(f"    ℓ={ell:2d}: E[dc_ℓ] = {expected:12.1f}")

    if total_cycles > 0:
        avg_length = total_length / total_cycles
        print(f"    Total expected cycles: {total_cycles:.1f}, avg length: {avg_length:.2f}")

print()

print("=" * 70)
print("PART 3: THE SIMPLEX-CUBOID PACKING RATIO AT EACH n")
print("=" * 70)

# (x+2)^n / (x+1)^n evaluated at integer points
# At x=1: 3^n / 2^n — the fundamental ratio
# At x=0: 2^n / 1 — the simplex size
# At x=-1: 1^n / 0 — degenerate

print("\n  The ratio (x+2)^n / (x+1)^n = ((x+2)/(x+1))^n")
print("  At x=1: (3/2)^n — grows exponentially")
print()

for n in range(1, 13):
    ratio = (3/2) ** n
    # Number of orthoschemes in cube = n!
    # Each orthoscheme has volume 1/n!
    # Regular simplex (Hadamard) has volume sqrt(n+1) / (2^n · n!)
    # Wait, that's the volume in R^n. Let me use the Hadamard formula.
    # For a regular simplex inscribed in {0,1}^n (Hadamard), the volume is:
    # V = sqrt(n+1) / (n! * 2^{n/2}) ... not sure of exact formula

    # Number of halves:
    # At n=2: 2 halves
    # At n=3: 4 halves
    # At n=k (Hadamard): 2^n - (n+1) corner points minus simplex...
    # Actually: the unit cube has 2^n vertices.
    # Regular simplex uses n+1 of them.
    # The "complementary" region has 2^n - (n+1) extra vertices.
    # Number of "corner pieces" = 2^n - (n+1)... but only for Hadamard n.

    has_hadamard = (n + 1) % 4 == 0 or n <= 1
    corner_count = 2**n - (n + 1) if has_hadamard else None

    print(f"  n={n:2d}: (3/2)^n = {ratio:10.3f}, "
          f"n! = {factorial(n):8d}, "
          f"2^n = {2**n:5d}, "
          f"corners = {corner_count if corner_count is not None else 'N/A (no Hadamard)':>8}")

print()

# The user's packing question: the simplex sits in the cube, and the
# remaining pieces are the "corners". At n=2: 2 corners. At n=3: 4 corners.
# At n=7: 2^7 - 8 = 120 corner vertices. But corners aren't determined by vertices alone.

# The key question: what is vol(cube)/vol(simplex) for Hadamard simplices?
print("  Volume ratios for Hadamard simplices:")
print("  (Regular simplex inscribed in unit cube)")
print()

for n in [1, 2, 3, 7]:
    # For Hadamard simplex in [0,1]^n:
    # Vertices are rows of normalized Hadamard matrix
    # Volume = sqrt(n+1) / (2^n · n!)... let me compute directly
    # Actually, volume of regular simplex with edge length a in R^n:
    # V = a^n * sqrt(n+1) / (n! * 2^(n/2))
    # For Hadamard in {0,1}^n: edge length a = sqrt((n+1)/2) (from our earlier computation)
    # But the simplex is embedded in n-dim affine subspace of R^n (sum of coords = (n+1)/2)
    # so we need the n-dimensional volume in that subspace.

    # Easiest: direct computation using the Cayley-Menger determinant
    # For a regular simplex with n+1 vertices in R^n, all edges = a:
    # n-vol = a^n * sqrt(n+1) / (n! * 2^{n/2})

    if n == 1:
        # Simplex = line segment [0,1], vol = 1 = cube vol
        ratio_vc = 1
    elif n == 2:
        # Equilateral triangle in unit square
        # Edge = 1 (if we use the standard embedding)
        # But for Hadamard: vertices at (0,0), (1,0), (0,1) → NOT equilateral
        # For equilateral in square: edge = sqrt(2), area = sqrt(3)/2*sqrt(2)^2/4... complex
        # Let me use: equilateral triangle with one side as base of unit square
        # Base = 1, height = sqrt(3)/2, area = sqrt(3)/4
        # Square area = 1*sqrt(3)/2 = sqrt(3)/2 (bounding rect)
        # Actually the user says "equilateral triangle sits in a square"
        # If triangle side = s, minimal bounding square has side = s * sin(75°)... complex
        # Simpler: triangle inscribed in unit square with vertices at
        # (0,0), (1,0), (1/2, sqrt(3)/2)
        # Bounding rectangle: [0,1] x [0, sqrt(3)/2]
        # Area of rectangle = sqrt(3)/2 ≈ 0.866
        # Area of triangle = sqrt(3)/4 ≈ 0.433
        # Ratio = 2. Two "halves" of equal area = sqrt(3)/8 each? No.
        # Remaining area = sqrt(3)/2 - sqrt(3)/4 = sqrt(3)/4
        # Split into 2 corner pieces: each sqrt(3)/8 = triangle_area/2. ✓!
        ratio_vc = 2
    elif n == 3:
        # Regular tetrahedron in unit cube
        # From our computation: cube/tet = 3, so 4 corner pieces total vol = 2/3
        # Each corner = 1/6 = (1/3)/2 = tet_vol/2
        ratio_vc = 3
    elif n == 7:
        # From simplex_cube_hadamard.py: ratio = 157.5
        ratio_vc = 157.5

    n_corners_calc = 2**n - (n + 1)
    if ratio_vc > 0:
        corner_each = (ratio_vc - 1) / n_corners_calc if n_corners_calc > 0 else 0
        print(f"  n={n}: cube/simplex = {ratio_vc}, "
              f"corners = {n_corners_calc}, "
              f"each corner/simplex = {corner_each:.4f}, "
              f"{'= 1/2 ✓' if abs(corner_each - 0.5) < 0.01 else ''}")

print()
print("  The 'halves' property (each corner = simplex/2) holds ONLY at n=2,3!")
print("  At n=7: each corner = 1.29× the simplex — they're LARGER than the simplex.")
print()

print("=" * 70)
print("PART 4: THE (x+1)/(x+2) FRAMEWORK FOR TOURNAMENTS")
print("=" * 70)

# For each tournament T at n=5, compute I(Ω,x) for x = 1,2,3
# and the simplex/cuboid evaluations
n = 5
edges = [(i, j) for i in range(n) for j in range(i+1, n)]
num_edges = len(edges)

print(f"\n  Computing I(Ω, x) for all n={n} tournaments:")
print()

# Collect unique (α₁, H, I₃) triples
from collections import Counter
triples = Counter()

for bits in range(2**num_edges):
    adj = [[False]*n for _ in range(n)]
    for idx, (i, j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i][j] = True
        else:
            adj[j][i] = True

    # Count 3-cycles
    dc3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    dc3 += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    dc3 += 1

    # Count 5-cycles
    dc5 = 0
    for verts in combinations(range(n), 5):
        for perm in permutations(verts):
            if all(adj[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                dc5 += 1
    dc5 //= 5

    alpha1 = dc3 + dc5
    # At n=5, no disjoint pairs possible (need 6 vertices for two 3-cycles)
    alpha2 = 0

    # I(x) = 1 + α₁x + α₂x²
    I2 = 1 + 2*alpha1 + 4*alpha2  # = H
    I3 = 1 + 3*alpha1 + 9*alpha2
    Im1 = 1 - alpha1 + alpha2     # = I(-1)

    triples[(alpha1, I2, I3, Im1)] += 1

print(f"  {'α₁':>4s} {'H=I(2)':>7s} {'I(3)':>5s} {'I(-1)':>5s} {'I(3)/I(2)':>9s} {'count':>6s}")
for (a1, i2, i3, im1), count in sorted(triples.items()):
    ratio = i3 / i2 if i2 > 0 else float('inf')
    print(f"  {a1:4d} {i2:7d} {i3:5d} {im1:5d} {ratio:9.4f} {count:6d}")

print()

# Now the key: the simplex and cuboid evaluations
print("  SIMPLEX-CUBOID DECOMPOSITION of I(Ω, x):")
print("  I(Ω, x) = 1 + α₁x + α₂x² + ...")
print("  In the (x+1)^k basis: simplex decomposition")
print("  In the (x+2)^k basis: cuboid decomposition")
print()
print("  Basis change: x = (x+1) - 1 = (x+2) - 2")
print("  So: I = 1 + α₁((x+1)-1) = (1-α₁) + α₁(x+1)")
print("  At n=5 (α₂=0): I = (1-α₁) + α₁·(x+1)")
print("  This is LINEAR in (x+1)!")
print()
print("  Simplex coefficients [c₀, c₁]: I = c₀ + c₁(x+1)")
print("  c₀ = I(x=-1) = I(-1) = 1-α₁")
print("  c₁ = α₁")
print()
print("  Cuboid coefficients [d₀, d₁]: I = d₀ + d₁(x+2)")
print("  d₀ = I(x=-2) = I(-2) = 1-2α₁")
print("  d₁ = α₁")
print()
print("  At x=1 (key 2): I(1) = (1-α₁) + 2α₁ = 1+α₁ (= H for 3-cycles only)")
print("  Wait, I(2) = H. Let me recompute.")
print("  I(x) = 1 + α₁x: I(2) = 1+2α₁ = H. ✓")
print("  In simplex basis: I = (1-α₁) + α₁(x+1)")
print("    At x=2: (1-α₁) + 3α₁ = 1+2α₁ = H. ✓")
print()

# The key identity: H = I(-1) + 3α₁ (at n=5 where α₂=0)
# More generally: H = I(-1) + 3α₁ + 12α₂ + ...
# Since I(-1) = 1-α₁+α₂-..., and H = 1+2α₁+4α₂+...:
# H - I(-1) = 3α₁ + 3α₂ + ... = 3(α₁+α₂+...) ? No.
# H - I(-1) = (1+2α₁+4α₂) - (1-α₁+α₂) = 3α₁+3α₂ = 3(α₁+α₂)
# For higher: H - I(-1) = 3α₁ + 3α₂ + ... hmm
# Actually: I(2) - I(-1) = Σ αₖ(2^k-(-1)^k) = Σ αₖ(2^k+1 for odd k, 2^k-1 for even k)
# For k=1: 3, k=2: 3, k=3: 9, k=4: 15, ...

print("  H - I(-1) = Σ αₖ(2^k - (-1)^k)")
print("  For k=1: coefficient 3 = 2+1")
print("  For k=2: coefficient 3 = 4-1")
print("  For k=3: coefficient 9 = 8+1")
print("  For k=4: coefficient 15 = 16-1")
print()
print("  These are 2^k + (-1)^{k+1} = 3, 3, 9, 15, 33, 63, ...")
print("  For k=1,2: both equal 3! The first two levels contribute equally.")
print("  This is why α₁ ≥ α₂ is 'natural': levels 1 and 2 have equal weight in H-I(-1).")
print()

print("=" * 70)
print("PART 5: THE α₁≥α₂ EXTENSION — AVERAGE CYCLE LENGTH BOUND")
print("=" * 70)

# For n ≥ 10, we need average cycle length > 3.
# Let's check: at n=10, what is the minimum average cycle length?
# A tournament on 10 vertices can have at most floor(10/3) = 3 disjoint 3-cycles
# using 9 vertices. The remaining vertex can't form a cycle alone.
# Maximum number of 3-cycles: C(10,3)/4 = 30 (for random, max is higher)
# But these 30 cycles certainly have average length 3.

# The question is: for tournaments where α₁ is LARGE enough to threaten α₁ < α₂,
# does the average cycle length increase?

# At n=10, the maximum α₁ (3-cycles only) is about 40 (for near-regular tournaments).
# We need the clique density to be high enough.

# Actually, let me look at this differently.
# The key formula: α₂ = C(α₁,2) - e(CG)
# We need α₁ ≥ α₂ = C(α₁,2) - e(CG)
# i.e., e(CG) ≥ C(α₁,2) - α₁ = α₁(α₁-3)/2

# The S75 proof uses: e(CG) ≥ Σ_v C(d_v,2) - s₂
# where Σ d_v = 3α₁, and bounds Σ C(d_v,2) by Cauchy-Schwarz.

# But we can also use a STRONGER clique property:
# not only do cycles through vertex v form a clique,
# but cycles through EDGE (u,v) also interact!

# Two 3-cycles sharing an edge must share 2 vertices.
# In a tournament, if arc u→v exists, then 3-cycles through edge u→v
# are: u→v→w→u for each w where v→w and w→u.
# Number of such w = |N_out(v) ∩ N_in(u)| = |{w : v→w, w→u}|

print("\n  For n=10, need a sharper bound than Cauchy-Schwarz on vertices.")
print("  Alternative: use the EDGE clique structure.")
print()
print("  For each arc u→v, the 3-cycles through this arc form a clique in CG.")
print("  Number of such 3-cycles = |N_out(v) ∩ N_in(u)| = common out-neighbors")
print()

# Let's verify this computationally at n=7 and n=8
# and check if the edge-based bound gives a stronger result

for n in [5, 6, 7]:
    if n > 7:
        break  # too slow
    print(f"  n={n}: checking edge-clique structure")
    edges_n = [(i, j) for i in range(n) for j in range(i+1, n)]
    num_edges_n = len(edges_n)

    max_ratio_a2_a1 = 0
    worst_tournament = None

    num_checked = 0
    for bits in range(2**num_edges_n):
        adj = [[False]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges_n):
            if bits & (1 << idx):
                adj[i][j] = True
            else:
                adj[j][i] = True

        # Count 3-cycles
        cycles = set()
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if adj[i][j] and adj[j][k] and adj[k][i]:
                        cycles.add(frozenset([i, j, k]))
                    elif adj[j][i] and adj[i][k] and adj[k][j]:
                        cycles.add(frozenset([i, j, k]))

        cycles = list(cycles)
        a1 = len(cycles)
        if a1 < 2:
            continue

        # Count disjoint pairs
        a2 = sum(1 for i in range(len(cycles)) for j in range(i+1, len(cycles))
                 if len(cycles[i] & cycles[j]) == 0)

        if a1 > 0:
            ratio = a2 / a1
            if ratio > max_ratio_a2_a1:
                max_ratio_a2_a1 = ratio
                worst_tournament = (a1, a2, bits)

        num_checked += 1

    if worst_tournament:
        a1, a2, bits = worst_tournament
        print(f"    max α₂/α₁ = {max_ratio_a2_a1:.4f} at α₁={a1}, α₂={a2}")
    print(f"    checked {num_checked} tournaments")

print()

print("=" * 70)
print("PART 6: THE PACKING PATTERN — HOW IT CONTINUES")
print("=" * 70)

# The user's question: equilateral triangle → 2 halves, tetrahedron → 4 halves
# What's the pattern?

# n=2: 2 corners = 2^2 - (2+1) = 1? No, 2^2=4 vertices, simplex uses 3.
# Wait: in R^2, the equilateral triangle has 3 vertices.
# Bounding rectangle (not cube) has 4 vertices.
# Corners = 4-3 = 1 extra vertex region → 2 corner triangles
# Each corner = simplex_vol/2

# n=3: in R^3, regular tet has 4 vertices from {0,1}^3.
# Cube has 8 vertices. Extra = 4 vertices → 4 corner tets
# Each corner = tet_vol/2

# The pattern: at n=2, number of corners = 2
# At n=3, number of corners = 4
# At n=k: number of corners = 2^k - (k+1) for Hadamard k

# The "halves" property:
# Total corner volume = cube - simplex
# If each corner = simplex/2, then total corner = (number of corners) * simplex/2
# So: cube - simplex = corners * simplex/2
# cube = simplex * (1 + corners/2)
# cube/simplex = 1 + corners/2 = 1 + (2^n-(n+1))/2

# At n=2: 1 + (4-3)/2 = 1.5. But cube/simplex = 2. Contradiction!
# Hmm, the "2 corners each = simplex/2" gives total corner = simplex.
# So cube = 2 * simplex. cube/simplex = 2. ✓

# At n=3: "4 corners each = simplex/2" gives total corner = 2*simplex.
# So cube = 3 * simplex. cube/simplex = 3. ✓

# General: if each corner = simplex/2, total corner = (2^n-n-1)*simplex/2
# cube/simplex = 1 + (2^n-n-1)/2 = (2^n-n+1)/2

# At n=2: (4-2+1)/2 = 3/2. But actual = 2. ❌
# The "corners" at n=2 aren't all the same size as at n=3!

# Let me think about this more carefully.
# At n=2: bounding RECTANGLE (not square!) has area = base × height = 1 × √3/2
# Triangle area = √3/4
# Ratio = 2. Two corners, each = √3/8 = (√3/4)/2 = simplex/2. ✓
# But wait: the "cube" in n=2 is the bounding rectangle of the triangle.
# This isn't the unit square — it's sized to exactly contain the equilateral triangle.

# The right framing: the SIMPLEX determines the cube/rectangle it sits in.
# At n=2: equilateral triangle → bounding rect = 2×triangle
# At n=3: regular tet → bounding cube = 3×tet (when inscribed in {0,1}^3 × Hadamard)

# THE PATTERN: cube/simplex = n at Hadamard dimensions?
# n=2: 2 ✓, n=3: 3 ✓, n=7: 157.5 ❌

# So the pattern BREAKS at n=7. Why?
# Because at n=7, the "cube" is [0,1]^7 and the simplex is tiny relative to it.
# The ratio depends on the specific embedding.

# Actually, I think the right comparison is:
# The ORTHOSCHEME has volume 1/n!
# The cube [0,1]^n has volume 1
# Ratio = n!
# At n=2: 2! = 2 ✓
# At n=3: 3! = 6, but we said ratio = 3 ❌

# Wait, the orthoscheme is NOT the regular simplex.
# The orthoscheme Δ = {1 ≥ x₁ ≥ x₂ ≥ ... ≥ x_n ≥ 0} is a RECTANGULAR simplex
# with volume 1/n!.
# The REGULAR simplex (equilateral) inscribed in {0,1}^n via Hadamard
# has a DIFFERENT volume.

# For the Hadamard simplex at n=3:
# Vertices: {(0,0,0), (1,1,0), (1,0,1), (0,1,1)} (from H_4 Sylvester)
# Edge length: sqrt(2) for all edges ✓ (distance between any two = sqrt(2))
# Volume of regular tet with edge a: V = a³/(6√2) = 2√2/(6√2) = 2/6 = 1/3
# Cube volume = 1
# Ratio = 3 ✓ (not 6)

# The Hadamard simplex uses n+1 vertices of {0,1}^n, each with
# exactly (n+1)/2 coordinates equal to 1.
# Its volume in R^n is sqrt(n+1)/(2^{n/2} · n!) × (edge length)^n... complex

print("\n  THE PACKING PATTERN:")
print()
print("  n=2: equilateral triangle in bounding rectangle")
print("    Vol ratio = 2. Two corners, each = triangle/2.")
print("    'Complementary piece' = 1 triangle")
print()
print("  n=3: regular tetrahedron in unit cube (Hadamard)")
print("    Vol ratio = 3. Four corners, each = tet/2.")
print("    'Complementary piece' = 2 tets")
print()
print("  THE RATIO cube/simplex = n+1 - 1 = n at n=2,3")
print("  Does this continue?")
print()

# Compute Hadamard simplex volumes
import numpy as np

def hadamard_simplex_volume(n):
    """Volume of regular simplex inscribed in {0,1}^n via Hadamard."""
    if n == 1:
        # Segment [0,1], volume = 1
        return 1.0
    if (n + 1) % 4 != 0 and n > 1:
        return None  # No Hadamard matrix

    # Build Sylvester Hadamard matrix of order n+1 (if n+1 is power of 2)
    # Otherwise use Paley or other construction
    m = n + 1
    # Check if m is power of 2
    if m & (m - 1) == 0:
        # Sylvester construction
        H = np.array([[1]])
        while H.shape[0] < m:
            H = np.block([[H, H], [H, -H]])
        # Normalize: map {-1,1} to {0,1}
        V = (1 + H) // 2
        # Take n+1 rows as vertices in R^n... but H is (n+1)×(n+1)
        # We need n+1 vertices in R^n
        # Drop the first column (all 1s after normalization)
        vertices = V[:, 1:].astype(float)
    else:
        return None  # Only handle Sylvester for now

    # Compute volume using the Cayley-Menger determinant
    # Or: volume = |det(M)|/n! where M has rows v_i - v_0
    M = vertices[1:] - vertices[0]
    vol = abs(np.linalg.det(M)) / factorial(n)
    return vol

for n in [1, 3, 7, 15]:
    vol = hadamard_simplex_volume(n)
    if vol is not None:
        cube_vol = 1.0
        ratio = cube_vol / vol
        n_corners = 2**n - (n + 1)
        corner_each_ratio = (ratio - 1) / n_corners if n_corners > 0 else 0
        print(f"  n={n:2d}: simplex vol = {vol:.6e}, "
              f"cube/simplex = {ratio:.1f}, "
              f"corners = {n_corners}, "
              f"each/simplex = {corner_each_ratio:.4f}")

print()

# Check: is cube/simplex = n always for Hadamard?
# n=1: ratio = 1/1 = 1, should be 1. ✓ (but n=1 is degenerate)
# n=3: ratio = 3. ✓
# n=7: let's see...

print("  CHECKING: is cube/simplex always = n for Hadamard?")
# n=3: vol = 1/3, ratio = 3 ✓
# n=7: vol = ?, ratio = ?
# The Hadamard matrix for n+1=8 (Sylvester):
if True:
    n = 7
    vol = hadamard_simplex_volume(n)
    if vol:
        print(f"  n=7: vol = {vol:.10e}, 1/vol = {1/vol:.1f}")
        # If ratio = 7: vol should be 1/7 ≈ 0.1429
        # If ratio = 7!: vol should be 1/5040

print()

print("=" * 70)
print("PART 7: THE n+1 DECOMPOSITION — WHY cube = n·simplex AT n=3")
print("=" * 70)

print()
print("  At n=3, the regular tet uses vertices {0000, 1100, 1010, 0110}")
print("  Wait: {(0,0,0), (1,1,0), (1,0,1), (0,1,1)} in R³")
print("  These are the EVEN-WEIGHT vertices of {0,1}³.")
print("  The ODD-WEIGHT vertices {(1,0,0), (0,1,0), (0,0,1), (1,1,1)}")
print("  form ANOTHER regular tet — the COMPLEMENT!")
print()
print("  The cube is decomposed into:")
print("  - Tet₁ (even weight vertices)")
print("  - Tet₂ (odd weight vertices)")
print("  - 4 corner orthoschemes (joining vertices of Tet₁ and Tet₂)")
print()
print("  Vol(Tet₁) = Vol(Tet₂) = 1/3")
print("  Vol(4 corners) = 1 - 2/3 = 1/3")
print("  Each corner = 1/12 = (1/3)/4 = Tet/4? No...")
print("  Wait: from our earlier result, each corner = Tet/2 = 1/6")
print("  4 × 1/6 = 2/3. But 1 - 1/3 = 2/3. ✓")
print("  (We're only considering ONE tet, not both)")
print()
print("  With both tets:")
print("  Cube = Tet₁ + Tet₂ + 4 edges × ... no, they overlap")
print("  Actually: Cube = 6 orthoschemes (= 3! chambers)")
print("  = n! = 6 chambers of the type A₂ arrangement")
print("  Tet₁ uses 2 chambers (even permutations)")
print("  Tet₂ uses 2 chambers (odd permutations)")
print("  Remaining: 2 chambers form the 'twist' region")
print()
print("  TOURNAMENT INTERPRETATION:")
print("  Transitive tournament: H=1 (1 chamber)")
print("  Regular tournament: H=3 (3 chambers)")
print("  Total chambers = 6 = 3!")
print("  H ranges from 1 to 3 at n=3.")
print()
print("  H=3 means the tournament 'fills' 3/6 = 1/2 of the cube")
print("  = volume of ONE regular tet = 1/3... no")
print("  Volume = H × vol(orthoscheme) = 3 × 1/6 = 1/2")
print("  At n=3: vol(union of 3 chambers) = 3/6 = 1/2")
print("  This equals (1/3 + 1/6) = the tet + 1 corner")
print()

# At n=3, the 3 Hamiltonian path chambers for the regular tournament:
# Path 0→1→2: chamber x₀ ≥ x₁ ≥ x₂
# Path 1→2→0: chamber x₁ ≥ x₂ ≥ x₀
# Path 2→0→1: chamber x₂ ≥ x₀ ≥ x₁
# These are the 3 EVEN permutations of (0,1,2)!
# Even permutations correspond to one of the two regular tets.

print("  REMARKABLE: For the regular 3-tournament,")
print("  the 3 HP chambers = 3 even permutations = one regular tet!")
print("  vol(3 chambers) = 1/2 ≠ vol(tet) = 1/3")
print("  Because the orthoscheme (1/6) ≠ tet slice (1/3 per tet)")
print()
print("  The orthoscheme {x₀≥x₁≥x₂} is a RIGHT simplex (1/6 of cube)")
print("  The regular tet (1/3 of cube) consists of 2 orthoschemes")
print("  So tet = 2 chambers. But H=3, and the 3 chambers aren't 1.5 tets...")
print()

print("=" * 70)
print("PART 8: THE MASTER IDENTITY — det(A_{n-1})=n AND packing")
print("=" * 70)

print()
print("  SYNTHESIS: det(A_{n-1}) = n encodes the fundamental packing ratio")
print()
print("  The Cartan matrix A_{n-1} of type A has:")
print("  - Diagonal entries: 2 (= KEY₁)")
print("  - Off-diagonal: -1 (adjacent) or 0")
print("  - Determinant: n")
print("  - Eigenvalues: 2 - 2cos(kπ/n), k=1,...,n-1")
print()
print("  The ratio cube/Hadamard-simplex depends on the NORMALIZATION.")
print("  In the standard {0,1}^n embedding:")
print("    n=3: ratio = 3 = det(A₂)")
print("    n=7: ratio ≈ 315/2 ≠ 7")
print()
print("  But in the WEIGHT LATTICE normalization:")
print("    The weight lattice Λ_w has index n = det(A_{n-1}) in the root lattice Λ_r")
print("    Vol(Λ_w)/Vol(Λ_r) = n (up to normalization)")
print("    This IS the fundamental packing ratio!")
print()
print("  CONNECTION TO OCF:")
print("  H(T) = I(Ω, 2) = I(Ω, det(A₁))")
print("  I(Ω, 3) = I(Ω, det(A₂))")
print("  The evaluation points are CARTAN DETERMINANTS!")
print()
print("  I(Ω, n) = I(Ω, det(A_{n-1})) for small n")
print("  At n=5: I(Ω, 5) = I(Ω, det(A₄))")
print("  This gives the PENTAGONAL evaluation of the independence polynomial.")
print()

# Compute I(Ω, k) for k = det(A_{k-1}) = k at n=5
print("  I(Ω, k) for all n=5 tournaments:")
print(f"  {'α₁':>4s} {'I(2)=H':>7s} {'I(3)':>6s} {'I(4)':>6s} {'I(5)':>6s}")
for (a1, i2, i3, im1), count in sorted(triples.items()):
    i4 = 1 + 4*a1  # α₂=0 at n=5
    i5 = 1 + 5*a1
    print(f"  {a1:4d} {i2:7d} {i3:6d} {i4:6d} {i5:6d}  (×{count})")

print()
print("  RATIOS I(k+1)/I(k):")
print(f"  {'α₁':>4s} {'I(3)/I(2)':>9s} {'I(4)/I(3)':>9s} {'I(5)/I(4)':>9s}")
for (a1, i2, i3, im1), count in sorted(triples.items()):
    i4 = 1 + 4*a1
    i5 = 1 + 5*a1
    r32 = i3/i2 if i2 > 0 else 0
    r43 = i4/i3 if i3 > 0 else 0
    r54 = i5/i4 if i4 > 0 else 0
    print(f"  {a1:4d} {r32:9.4f} {r43:9.4f} {r54:9.4f}")

print()
print("  As α₁→∞: I(k)/I(k-1) → k/(k-1)")
print("  I(3)/I(2) → 3/2, I(4)/I(3) → 4/3, I(5)/I(4) → 5/4")
print("  These are consecutive Cartan determinant RATIOS!")
print()
print("  CONJECTURE (HYP-1009):")
print("  For any tournament T with α₁ cycles:")
print("  I(k+1)/I(k) → (k+1)/k as α₁→∞")
print("  This is the ratio of consecutive Cartan determinants.")
print("  The 'simplex' tournament (H=1) has ratio = 1 (degenerate)")
print("  The 'maximally cyclic' tournament approaches (k+1)/k")

print("\nDone.")
