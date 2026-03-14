#!/usr/bin/env python3
"""
simplex_cuboid_nesting.py — How simplices nest inside hypercubes

When an n-simplex is inscribed in an n-cube (using vertices of the cube),
the complement of the simplex inside the cube decomposes into pieces.

Key examples:
  n=2: right triangle in square → 1 complementary triangle
  n=3: regular tetrahedron in cube → 4 corner right-tetrahedra

This script investigates:
  1. The standard simplex-in-cube embedding (Coxeter decomposition)
  2. The regular Hadamard simplex (n=3,7,11,...)
  3. Volume ratios and complement piece counts
  4. Connection to tournament combinatorics
"""

import numpy as np
from itertools import combinations, permutations
from math import factorial, comb, sqrt
from fractions import Fraction

def simplex_volume(vertices):
    """Volume of simplex with given vertices (as numpy array)."""
    n = len(vertices) - 1
    if n == 0:
        return 1
    M = np.array([vertices[i] - vertices[0] for i in range(1, n+1)], dtype=float)
    return abs(np.linalg.det(M)) / factorial(n)

# ═══════════════════════════════════════════════════════════════════════
print("=" * 72)
print("SIMPLEX-IN-CUBE NESTING: How n-simplices sit inside n-cubes")
print("=" * 72)

# ═══════════════════════════════════════════════════════════════════════
# PART 1: The Coxeter Decomposition (universal, all n)
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 1: COXETER DECOMPOSITION — The n-cube as n! simplices
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The unit n-cube [0,1]^n can be decomposed into n! congruent simplices,
one for each permutation σ ∈ S_n:

  Δ_σ = { x ∈ [0,1]^n : 0 ≤ x_{σ(1)} ≤ x_{σ(2)} ≤ ··· ≤ x_{σ(n)} ≤ 1 }

Each simplex has volume 1/n!. They tile the cube perfectly.
The simplex Δ_id (identity permutation) has vertices:
  v_0 = (0,...,0), v_1 = (1,0,...,0), v_2 = (1,1,0,...,0), ..., v_n = (1,...,1)

This is a RIGHT-ANGLED simplex (orthoscheme), not regular.
Its complement in the cube consists of the other n! - 1 simplices.
""")

print(f"{'n':>3} │ {'n! simplices':>12} │ {'Vol each':>12} │ {'complement':>12}")
print("─" * 48)
for n in range(1, 11):
    nf = factorial(n)
    print(f"{n:>3} │ {nf:>12} │ {'1/'+str(nf):>12} │ {nf-1:>12}")

# ═══════════════════════════════════════════════════════════════════════
# PART 2: n=3 — Regular Tetrahedron in Cube (the classic example)
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 2: n=3 — REGULAR TETRAHEDRON IN CUBE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The 8 vertices of the unit cube split by coordinate-sum parity:

  Even parity (sum ≡ 0 mod 2): (0,0,0), (1,1,0), (1,0,1), (0,1,1)
  Odd  parity (sum ≡ 1 mod 2): (1,0,0), (0,1,0), (0,0,1), (1,1,1)

The 4 even-parity vertices form a REGULAR tetrahedron:
  All pairwise Euclidean distances = √2.
""")

tet = np.array([[0,0,0],[1,1,0],[1,0,1],[0,1,1]], dtype=float)
vol_tet = simplex_volume(tet)
print(f"  Vol(regular tetrahedron) = {Fraction(vol_tet).limit_denominator(100)}")

# Verify all distances
dists = set()
for i in range(4):
    for j in range(i+1, 4):
        d = np.linalg.norm(tet[i] - tet[j])
        dists.add(round(d, 6))
print(f"  Pairwise distances: {dists}  →  REGULAR ✓")

# Corner pieces
print(f"""
The 4 unused (odd-parity) vertices each generate a corner piece:
each unused vertex + its 3 cube-adjacent tetrahedron vertices
form a right tetrahedron.
""")

corners_3 = {
    (1,0,0): [(0,0,0), (1,1,0), (1,0,1)],
    (0,1,0): [(0,0,0), (1,1,0), (0,1,1)],
    (0,0,1): [(0,0,0), (1,0,1), (0,1,1)],
    (1,1,1): [(1,1,0), (1,0,1), (0,1,1)],
}

total_vol = vol_tet
for uv, adj in corners_3.items():
    verts = np.array([list(uv)] + [list(a) for a in adj], dtype=float)
    v = simplex_volume(verts)
    total_vol += v
    print(f"  Corner at {uv}: adjacent to {adj}, vol = {Fraction(v).limit_denominator(100)}")

print(f"\n  Total: {Fraction(vol_tet).limit_denominator(100)} + 4×(1/6) = "
      f"{Fraction(total_vol).limit_denominator(100)} ✓")
print(f"\n  RESULT: 1 tetrahedron + 4 corner pieces = 5 pieces total")
print(f"  Complement = 4 pieces")

# ═══════════════════════════════════════════════════════════════════════
# PART 3: n=2 — Triangle in Square
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 3: n=2 — TRIANGLE IN SQUARE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The 4 vertices of the unit square: (0,0), (1,0), (0,1), (1,1).
We cannot inscribe an EQUILATERAL triangle (no 3 square vertices are
equidistant — distances are always 1,1,√2).

Instead, take any 3 vertices, e.g., (0,0), (1,0), (0,1).
This right triangle has area = 1/2.
The complement is the other right triangle: area = 1/2.

The diagonal of the square splits it into 2 congruent right triangles.
  → 2 pieces total (= 2! = n!)
  → 1 complement piece (= n! - 1 = 1)

This IS the Coxeter decomposition for n=2:
  Δ_{(1,2)} = { x_1 ≤ x_2 } (below diagonal)
  Δ_{(2,1)} = { x_2 ≤ x_1 } (above diagonal)
""")

# ═══════════════════════════════════════════════════════════════════════
# PART 4: Max-volume inscribed simplices for n=2,3,4,5,6
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 4: MAX-VOLUME SIMPLEX INSCRIBED IN n-CUBE (vertices at {0,1}^n)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

For each n, we find the maximum volume of an n-simplex whose n+1
vertices are chosen from the 2^n vertices of {0,1}^n.

This is related to the HADAMARD MAXIMAL DETERMINANT PROBLEM:
  max |det(M)| over (n+1)×(n+1) matrices with entries in {0,1}
  = max simplex volume × n!
""")

max_vol_data = {}

for n in range(1, 7):
    nf = factorial(n)
    num_verts = 2**n

    # Generate all cube vertices
    all_v = np.array([[(i >> bit) & 1 for bit in range(n)]
                      for i in range(num_verts)], dtype=float)

    # Brute force for small n, smarter for larger
    if comb(num_verts, n+1) <= 2_000_000:
        max_v = 0
        best = None
        count = 0
        for combo in combinations(range(num_verts), n+1):
            verts = all_v[list(combo)]
            M = np.array([verts[i] - verts[0] for i in range(1, n+1)], dtype=float)
            det_abs = abs(np.linalg.det(M))
            vol = det_abs / nf
            if vol > max_v + 1e-12:
                max_v = vol
                best = verts.copy()
                count = 1
            elif abs(vol - max_v) < 1e-12:
                count += 1

        max_vol_data[n] = (max_v, best, count)
        k = round(max_v * nf)

        # Check regularity
        dist_set = set()
        for i in range(n+1):
            for j in range(i+1, n+1):
                d2 = np.sum((best[i] - best[j])**2)
                dist_set.add(round(d2, 6))
        regular = len(dist_set) == 1

        unused = num_verts - (n + 1)
        vol_frac = Fraction(max_v).limit_denominator(100000)

        print(f"  n={n}: Vol = {str(vol_frac):>8}, "
              f"Coxeter cells = {k:>3}/{nf:<5}, "
              f"unused verts = {unused:>3}, "
              f"regular = {'YES' if regular else 'no':>3}, "
              f"#{count} max simplices")
        print(f"        Vertices: {[tuple(v.astype(int)) for v in best]}")
        print(f"        Dist² set: {sorted(dist_set)}")
        print()
    else:
        print(f"  n={n}: C({num_verts},{n+1}) = {comb(num_verts,n+1):,} — "
              f"too many for brute force, using known results")
        # Known results from Hadamard maximal determinant problem
        # n=6: max |det| of 7×7 {0,1} matrix
        # Best known: not regular, max vol = ?
        # Will compute via random sampling below
        max_vol_data[n] = None
        print()

# ═══════════════════════════════════════════════════════════════════════
# PART 5: Hadamard construction for n=7
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 5: HADAMARD SIMPLEX — Regular simplex in n-cube
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

A REGULAR n-simplex can be inscribed in {0,1}^n IFF a Hadamard matrix
of order n+1 exists. This requires n+1 ≡ 0 (mod 4) [for n≥3], so:

  n = 1, 3, 7, 11, 15, 19, 23, 27, 31, ...

For n=3: the regular tetrahedron (analyzed above).
For n=7: use the Sylvester-Hadamard construction H_8.
""")

# Build H_8 via Kronecker products
H2 = np.array([[1,1],[1,-1]])
H4 = np.kron(H2, H2)
H8 = np.kron(H2, H4)

# Map ±1 to {0,1}: entry = (1+h)/2. Drop first column (all 1s).
cube01 = ((1 + H8) // 2).astype(int)
simp7 = cube01[:, 1:]  # 8 vertices in {0,1}^7

print("n=7 Hadamard simplex vertices:")
for i, v in enumerate(simp7):
    print(f"  V{i} = {tuple(v)}  (sum = {sum(v)})")

# Verify regularity
d2_set = set()
for i in range(8):
    for j in range(i+1, 8):
        d2_set.add(int(np.sum((simp7[i] - simp7[j])**2)))
print(f"\nPairwise squared distances: {sorted(d2_set)}")
print(f"Regular: {'YES' if len(d2_set)==1 else 'NO'}")

# Volume
M7 = np.array([simp7[i] - simp7[0] for i in range(1, 8)], dtype=float)
vol7 = abs(np.linalg.det(M7)) / factorial(7)
vol7_frac = Fraction(vol7).limit_denominator(100000)
k7 = round(vol7 * factorial(7))

print(f"\nVolume = {vol7_frac} ≈ {vol7:.8f}")
print(f"Coxeter cells: {k7} of {factorial(7)}")
print(f"Complement Coxeter cells: {factorial(7) - k7}")
print(f"Unused cube vertices: {2**7 - 8} = {128-8}")

# ═══════════════════════════════════════════════════════════════════════
# PART 6: The corner piece structure
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 6: CORNER PIECE ANALYSIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

At n=3, each of the 4 unused cube vertices gives rise to exactly 1
corner simplex (right tetrahedron), and these 4 corner pieces exactly
tile the complement. This works because:

  1. Each unused vertex v is adjacent (Hamming distance 1) to exactly
     n simplex vertices.
  2. The convex hull of {v} ∪ {adjacent simplex vertices} is a right
     simplex of volume 1/n!.
  3. The 2^n - (n+1) such right simplices tile the complement.

Does this hold for general n? Let's check adjacency counts.
""")

for n_check in [2, 3, 4, 5]:
    nf = factorial(n_check)
    if n_check not in max_vol_data or max_vol_data[n_check] is None:
        continue
    max_v, best, _ = max_vol_data[n_check]
    if best is None:
        continue

    best_set = set(tuple(v.astype(int)) for v in best)
    all_v_set = set()
    for i in range(2**n_check):
        v = tuple((i >> bit) & 1 for bit in range(n_check))
        all_v_set.add(v)
    unused = sorted(all_v_set - best_set)

    print(f"\nn={n_check}: {len(unused)} unused vertices")
    adj_counts = {}
    for uv in unused:
        adj = [sv for sv in best_set if sum(a!=b for a,b in zip(uv,sv))==1]
        adj_count = len(adj)
        adj_counts[adj_count] = adj_counts.get(adj_count, 0) + 1
        if len(unused) <= 12:
            print(f"  {uv}: {adj_count} adjacent simplex vertices")

    print(f"  Adjacency distribution: {dict(sorted(adj_counts.items()))}")

    # Check if corner pieces tile the complement
    total_corner_vol = 0
    for uv in unused:
        adj = [sv for sv in best_set if sum(a!=b for a,b in zip(uv,sv))==1]
        if len(adj) >= n_check:
            # Can form an n-simplex
            corner_verts = np.array([list(uv)] + [list(a) for a in adj[:n_check+1]], dtype=float)
            try:
                cv = simplex_volume(corner_verts[:n_check+1])
                total_corner_vol += cv
            except:
                pass

    comp_vol = 1 - max_v
    print(f"  Complement volume: {Fraction(comp_vol).limit_denominator(10000)}")
    print(f"  Sum of corner volumes: {Fraction(total_corner_vol).limit_denominator(10000)}")
    if abs(total_corner_vol - comp_vol) < 1e-10:
        print(f"  Corner pieces EXACTLY tile complement ✓")
    else:
        print(f"  Corner pieces do NOT tile complement (gap = {comp_vol - total_corner_vol:.6f})")

# ═══════════════════════════════════════════════════════════════════════
# PART 7: The sequence of complement Coxeter cells
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 7: HOW MANY COXETER CELLS IN THE MAX INSCRIBED SIMPLEX?
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The max-volume inscribed simplex occupies k Coxeter cells of volume 1/n!.
So k = Vol(simplex) × n!
""")

print(f"{'n':>3} │ {'Vol':>10} │ {'k cells':>8} │ {'n!':>8} │ {'comp cells':>11} │ {'2^n-(n+1)':>10}")
print("─" * 60)
for n in range(1, 7):
    if n in max_vol_data and max_vol_data[n] is not None:
        max_v, _, _ = max_vol_data[n]
        nf = factorial(n)
        k = round(max_v * nf)
        comp = nf - k
        unused = 2**n - (n+1)
        vol_frac = Fraction(max_v).limit_denominator(100000)
        print(f"{n:>3} │ {str(vol_frac):>10} │ {k:>8} │ {nf:>8} │ {comp:>11} │ {unused:>10}")

# Add n=7 Hadamard
n = 7
nf = factorial(7)
k = k7
comp = nf - k
unused = 2**7 - 8
print(f"  7 │ {str(vol7_frac):>10} │ {k7:>8} │ {nf:>8} │ {comp:>11} │ {unused:>10}")

print("""
Coxeter cells in simplex: 1, 1, 2, 3, ?, ?, 32, ...
""")

# ═══════════════════════════════════════════════════════════════════════
# PART 8: Volume ratios
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 8: VOLUME RATIOS — How much of the cube does the simplex fill?
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Three different simplices in the n-cube:
  A) Coxeter orthoscheme: Vol = 1/n!  (right-angled, always exists)
  B) Max-volume inscribed: Vol = max over (n+1)-subsets of {0,1}^n
  C) Hadamard (regular): Vol = √(n+1)·((n+1)/4)^{n/2}/n!  (when exists)
""")

print(f"{'n':>3} │ {'Coxeter 1/n!':>14} │ {'Max inscribed':>14} │ {'Regular':>14}")
print("─" * 54)

for n in range(1, 11):
    nf = factorial(n)
    cox = f"1/{nf}"

    # Max inscribed from our data
    if n in max_vol_data and max_vol_data[n] is not None:
        max_v = max_vol_data[n][0]
        max_str = str(Fraction(max_v).limit_denominator(100000))
    elif n == 7:
        max_str = str(vol7_frac)
    else:
        max_str = "?"

    # Hadamard regular (exists for n+1 ≡ 0 mod 4 and n≥1)
    if (n+1) % 4 == 0 or n <= 1:
        # Volume formula: sqrt(n+1) * ((n+1)/4)^{n/2} / n!
        if n == 1:
            had_str = "1"
        elif n == 3:
            had_str = "1/3"
        elif n == 7:
            had_str = str(vol7_frac)
        else:
            had_vol = sqrt(n+1) * ((n+1)/4)**(n/2) / nf
            had_str = f"~{had_vol:.6f}"
    else:
        had_str = "---"

    print(f"{n:>3} │ {cox:>14} │ {max_str:>14} │ {had_str:>14}")

# ═══════════════════════════════════════════════════════════════════════
# PART 9: Connection to Tournaments
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 9: CONNECTION TO TOURNAMENTS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The Coxeter decomposition connects cubes, simplices, and tournaments:

1. Each of the n! Coxeter simplices corresponds to a PERMUTATION σ ∈ S_n.
   Each permutation defines a TRANSITIVE TOURNAMENT on [n]:
     vertex i beats vertex j iff σ(i) < σ(j).

2. So: n! Coxeter simplices ↔ n! transitive tournaments ↔ n! linear orders.

3. The n-cube [0,1]^n tiles into these n! simplices the same way the
   set of all orderings decomposes into n! total orders.

4. The user's binomial connection:
     (x+1)^n at x=2: 3^n
     (x+2)^n at x=2: 4^n
     Difference: 4^n - 3^n
   This is the number of {0,1,2,3}^n points NOT in {0,1,2}^n.
   In a scaled cube, this counts lattice points in the "outer shell."
""")

print(f"{'n':>3} │ {'n!':>8} │ {'2^C(n,2)':>10} │ {'n!/2^C(n,2)':>12} │ {'4^n-3^n':>10}")
print("─" * 52)
for n in range(2, 10):
    nf = factorial(n)
    t = 2**comb(n,2)
    ratio = Fraction(nf, t)
    diff = 4**n - 3**n
    print(f"{n:>3} │ {nf:>8} │ {t:>10} │ {str(ratio):>12} │ {diff:>10}")

print("""
Note: n!/2^{C(n,2)} = fraction of tournaments that are transitive.
This goes to 0 exponentially — most tournaments have cycles.
The Coxeter simplex (vol = 1/n!) captures the "transitive fraction"
of the cube's geometry.
""")

# ═══════════════════════════════════════════════════════════════════════
# PART 10: Resolving the user's pattern 2, 4, ?, ...
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 10: RESOLVING THE PATTERN — "2 halves, 4 halves, ?, ..."
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The user observes:
  n=2: "equilateral triangle in square, 2 halves on either side"
  n=3: "tetrahedron in cube, 4 halves around it"

INTERPRETATION: "halves" = complement pieces. The user asks how this
pattern continues.

For the REGULAR (Hadamard) simplex construction:
  n=3: 4 unused vertices → 4 corner simplices → 4 complement pieces

For the COXETER (right-angled) simplex:
  n=2: 2 total pieces (diagonal cuts square in half)
  n=3: 6 total pieces (but complement = 5, not 4)

RECONCILIATION: The user's n=2 case ("2 halves") = 2 total pieces
in the Coxeter decomposition. The n=3 case ("4 halves around it")
= the 4 corner pieces of the regular tetrahedron construction.

These are DIFFERENT constructions! For n=2 there is no regular simplex
in the square, so the Coxeter view (2 total) applies.

THE TWO NATURAL SEQUENCES:
""")

print("Sequence A (Coxeter): complement pieces = n! - 1")
print(f"  n:      ", end="")
for n in range(1, 9):
    print(f"{n:>6}", end="")
print(f"\n  pieces: ", end="")
for n in range(1, 9):
    print(f"{factorial(n)-1:>6}", end="")
print(f"\n  (0, 1, 5, 23, 119, 719, 5039, 40319)")

print(f"\nSequence B (Hadamard corner pieces): 2^n - (n+1)")
print(f"  n:      ", end="")
for n in range(1, 12):
    print(f"{n:>6}", end="")
print(f"\n  pieces: ", end="")
for n in range(1, 12):
    print(f"{2**n-(n+1):>6}", end="")
print(f"\n  (0, 1, 4, 11, 26, 57, 120, 247, 502, 1013, 2036)")
print(f"  = OEIS A000325")

print(f"""
For the Hadamard case (n = 3, 7, 11, 15, 19, 23, ...):
  The complement of the regular simplex has AT MOST 2^n - (n+1)
  convex corner pieces, one per unused vertex.

  At n=3: 4 corners, each with vol 1/6, total 4/6 = 2/3 ✓
  At n=7: 120 unused vertices, but corner pieces may overlap or
  have lower dimension — need careful analysis.
""")

# Check: at n=7, how many simplex vertices is each unused vertex adjacent to?
print("n=7: adjacency of unused vertices to simplex vertices:")
simp7_set = set(tuple(v) for v in simp7)
all7 = set()
for i in range(128):
    v = tuple((i >> bit) & 1 for bit in range(7))
    all7.add(v)
unused7 = sorted(all7 - simp7_set)

adj_dist_7 = {}
for uv in unused7:
    adj = sum(1 for sv in simp7_set if sum(a!=b for a,b in zip(uv,sv))==1)
    adj_dist_7[adj] = adj_dist_7.get(adj, 0) + 1

print(f"  {len(unused7)} unused vertices")
print(f"  Adjacency distribution (# adjacent simplex verts → count):")
for k in sorted(adj_dist_7):
    print(f"    {k} adjacent: {adj_dist_7[k]} vertices")

# ═══════════════════════════════════════════════════════════════════════
# PART 11: Master Summary
# ═══════════════════════════════════════════════════════════════════════

print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 PART 11: MASTER SUMMARY TABLE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")

print(f"{'n':>3} │ {'n!':>6} │ {'Coxeter':>8} │ {'Max Vol':>8} │ {'Max cells':>9} │ "
      f"{'Comp cells':>10} │ {'Unused':>6} │ {'Regular?':>8}")
print("─" * 74)

for n in range(1, 8):
    nf = factorial(n)
    cox_comp = nf - 1
    unused = 2**n - (n+1)

    if n in max_vol_data and max_vol_data[n] is not None:
        mv, _, _ = max_vol_data[n]
        k = round(mv * nf)
        comp = nf - k
        vf = str(Fraction(mv).limit_denominator(100000))
        # Check regularity
        _, best, _ = max_vol_data[n]
        d2s = set()
        for i in range(n+1):
            for j in range(i+1, n+1):
                d2s.add(round(np.sum((best[i]-best[j])**2), 4))
        reg = "YES" if len(d2s)==1 else "no"
    elif n == 7:
        k = k7
        comp = nf - k7
        vf = str(vol7_frac)
        reg = "YES"
    else:
        k = "?"
        comp = "?"
        vf = "?"
        reg = "?"

    print(f"{n:>3} │ {nf:>6} │ {cox_comp:>8} │ {vf:>8} │ {str(k):>9} │ "
          f"{str(comp):>10} │ {unused:>6} │ {reg:>8}")

print("""
ANSWERS TO THE USER'S QUESTIONS:

1. PIECES when n-simplex inscribed in n-cube:

   COXETER construction (right-angled simplex, all n):
     Total pieces = n!. Complement pieces = n! - 1.
     n=2: 2 total (matches "2 halves")
     n=3: 6 total, 5 complement
     n=4: 24 total, 23 complement

   HADAMARD construction (regular simplex, n=3,7,11,...):
     n=3: 5 total (1 tetra + 4 corners)
     n=7: 1 simplex + up to 120 corner pieces

2. VOLUME RATIOS:
   Coxeter: Vol(simplex)/Vol(cube) = 1/n! always.
   Max inscribed: 1/2, 1/3, 1/8, ?, ?, 2/315, ... (Hadamard bound)

3. The user's pattern "2, 4, ?, ?" is BEST explained as:
   n=2: 2 total Coxeter pieces (the diagonal cuts the square in 2)
   n=3: 4 corner pieces (regular tetrahedron complement)
   n=4: 23 Coxeter complement pieces (no regular simplex exists)
   n=5: 119 Coxeter complement pieces

   The n=2 and n=3 cases use DIFFERENT constructions, which is why
   a single clean sequence is hard to find.

   If forced to one sequence: n! - 1 = 1, 5, 23, 119, 719, ...
   (Coxeter complement) grows as n! and has no simpler form.

4. TOURNAMENT CONNECTION:
   The n! Coxeter simplices ↔ n! transitive tournaments on n vertices.
   The simplex-in-cube geometry IS the geometry of linear orders
   inside the space of all orderings.
""")
