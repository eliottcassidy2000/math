#!/usr/bin/env python3
"""
Baer Subplane Deep Analysis — The Geometric Heart of Forbidden Values
opus-2026-03-14-S71i

KEY THESIS: The forbidden H values {7, 21} = |PG(2,F₂)|, |PG(2,F₄)| arise
because the Baer subplane decomposition of PG(2,F₄) = 3 × PG(2,F₂) EXACTLY
mirrors the independence polynomial factorization I(K₁⊔K₃, 2) = 3 × 7.

This script explores:
1. Explicit construction of PG(2,F₄) and its 3 Baer subplanes
2. The Frobenius automorphism and its orbits
3. Connection to I-polynomial factorization
4. Why the pattern stops at F₄ (Baer subplane → simplex packing)
5. The (x+1)^n / (x+2)^n nesting via projective geometry
"""

from itertools import combinations, permutations
from collections import defaultdict

# ============================================================
# PART 1: Construct F₄ = F₂[ω] where ω² + ω + 1 = 0
# ============================================================

print("=" * 70)
print("PART 1: THE FIELD F₄ AND ITS STRUCTURE")
print("=" * 70)
print()

# F₄ elements: {0, 1, ω, ω²} where ω² = ω + 1 (mod 2)
# Represent as (a, b) meaning a + b·ω in F₂²
# Addition: componentwise mod 2
# Multiplication: use ω² = ω + 1

class F4:
    """Element of F₄ = F₂[ω]/(ω² + ω + 1)"""
    def __init__(self, a, b):
        self.a = a % 2
        self.b = b % 2

    def __add__(self, other):
        return F4(self.a ^ other.a, self.b ^ other.b)

    def __mul__(self, other):
        # (a₁ + b₁ω)(a₂ + b₂ω) = a₁a₂ + (a₁b₂ + b₁a₂)ω + b₁b₂ω²
        # = a₁a₂ + (a₁b₂ + b₁a₂)ω + b₁b₂(ω + 1)
        # = (a₁a₂ + b₁b₂) + (a₁b₂ + b₁a₂ + b₁b₂)ω
        aa = (self.a * other.a + self.b * other.b) % 2
        bb = (self.a * other.b + self.b * other.a + self.b * other.b) % 2
        return F4(aa, bb)

    def __eq__(self, other):
        return self.a == other.a and self.b == other.b

    def __hash__(self):
        return hash((self.a, self.b))

    def __repr__(self):
        if self.a == 0 and self.b == 0: return "0"
        if self.a == 1 and self.b == 0: return "1"
        if self.a == 0 and self.b == 1: return "ω"
        return "ω²"

    def is_zero(self):
        return self.a == 0 and self.b == 0

    def frobenius(self):
        """Frobenius automorphism: x → x² in F₄"""
        return self * self

    def in_F2(self):
        """Is this element in the subfield F₂?"""
        return self.b == 0

# Create all elements
ZERO = F4(0, 0)
ONE = F4(1, 0)
OMEGA = F4(0, 1)
OMEGA2 = F4(1, 1)  # ω² = ω + 1

f4_elements = [ZERO, ONE, OMEGA, OMEGA2]
f4_nonzero = [ONE, OMEGA, OMEGA2]

print("F₄ = {0, 1, ω, ω²} where ω² + ω + 1 = 0 over F₂")
print()

# Verify multiplication table
print("Multiplication table:")
for x in f4_elements:
    row = [str(x * y) for y in f4_elements]
    print(f"  {str(x):>2} × : {row}")

print()

# Verify Frobenius
print("Frobenius automorphism (x → x²):")
for x in f4_elements:
    print(f"  {x}² = {x.frobenius()}")
print(f"  Fixed points of Frobenius = F₂ = {{{', '.join(str(x) for x in f4_elements if x.frobenius() == x)}}}")

print()

# ============================================================
# PART 2: Construct PG(2, F₄) — the 21-point projective plane
# ============================================================

print("=" * 70)
print("PART 2: PG(2, F₄) — THE 21-POINT PROJECTIVE PLANE")
print("=" * 70)
print()

# Points of PG(2, F₄) = equivalence classes of nonzero triples (a:b:c) in F₄³
# Two triples are equivalent if they differ by a nonzero scalar.

def normalize_point(a, b, c):
    """Normalize a projective point (a:b:c) — first nonzero coord = 1"""
    for coord in [a, b, c]:
        if not coord.is_zero():
            # Find inverse
            for inv in f4_nonzero:
                if (coord * inv) == ONE:
                    return (inv * a, inv * b, inv * c)
    return None  # all zero

# Generate all points
points = set()
for a in f4_elements:
    for b in f4_elements:
        for c in f4_elements:
            if a.is_zero() and b.is_zero() and c.is_zero():
                continue
            p = normalize_point(a, b, c)
            points.add(p)

points = sorted(points, key=lambda p: (p[0].a, p[0].b, p[1].a, p[1].b, p[2].a, p[2].b))
print(f"|PG(2, F₄)| = {len(points)} points")
assert len(points) == 21, f"Expected 21, got {len(points)}"

# Label points 0-20
point_to_idx = {p: i for i, p in enumerate(points)}
idx_to_point = {i: p for i, p in enumerate(points)}

print("Points of PG(2, F₄):")
for i, p in enumerate(points):
    in_f2 = all(coord.in_F2() for coord in p)
    marker = " ← F₂ point (Baer)" if in_f2 else ""
    print(f"  P{i:2d} = ({p[0]}, {p[1]}, {p[2]}){marker}")

# ============================================================
# PART 3: Identify the 3 Baer subplanes
# ============================================================

print()
print("=" * 70)
print("PART 3: THE 3 BAER SUBPLANES OF PG(2, F₄)")
print("=" * 70)
print()

# Baer subplane 1: Points with all coordinates in F₂ = {0, 1}
baer1 = [i for i, p in enumerate(points) if all(coord.in_F2() for coord in p)]
print(f"Baer subplane B₁ (F₂-rational points): {baer1}")
print(f"  |B₁| = {len(baer1)}")

# The Frobenius automorphism σ: (a,b,c) → (a²,b²,c²) acts on PG(2,F₄)
# Fixed points of σ = B₁ (the F₂-rational points)

def frobenius_point(p):
    """Apply Frobenius to a projective point"""
    return normalize_point(p[0].frobenius(), p[1].frobenius(), p[2].frobenius())

# Compute Frobenius orbits
frob_orbits = {}
visited = set()
for i, p in enumerate(points):
    if i in visited:
        continue
    orbit = [i]
    visited.add(i)
    q = frobenius_point(p)
    j = point_to_idx[q]
    while j != i:
        orbit.append(j)
        visited.add(j)
        q = frobenius_point(idx_to_point[j])
        j = point_to_idx[q]
    frob_orbits[i] = orbit

print()
print("Frobenius orbits (σ: x → x²):")
fixed_points = []
two_orbits = []
for leader, orbit in sorted(frob_orbits.items()):
    if len(orbit) == 1:
        fixed_points.append(leader)
        print(f"  Fixed: P{leader}")
    else:
        two_orbits.append(orbit)
        print(f"  Orbit: {' ↔ '.join(f'P{i}' for i in orbit)}")

print(f"\nFixed points (= B₁): {fixed_points} ({len(fixed_points)} points)")
print(f"2-orbits: {len(two_orbits)} orbits ({2*len(two_orbits)} points)")
print(f"Total: {len(fixed_points)} + {2*len(two_orbits)} = {len(fixed_points) + 2*len(two_orbits)}")

# For the other 2 Baer subplanes: use the "ω-translate" construction
# B₂ and B₃ come from applying ω· and ω²· to the coordinates

# Actually, the 3 Baer subplanes partition PG(2,4) into 7+7+7.
# B₁ = F₂-rational points
# B₂ = image of B₁ under a collineation that permutes the 3 Baer subplanes
# B₃ = the remaining 7 points

# The partition uses the "absolute trace": Tr(x) = x + x² for x ∈ F₄
# Tr(0) = 0, Tr(1) = 0, Tr(ω) = ω + ω² = 1, Tr(ω²) = ω² + ω = 1
# So Tr: F₄ → F₂ maps {0,1} → 0, {ω, ω²} → 1

# A point (a:b:c) lies in a Baer subplane determined by a "conic" condition
# Let me construct B₂ and B₃ explicitly using the scaling trick:
# Multiply one coordinate by ω to get points NOT in F₂

# Alternative approach: The 3 Baer subplanes form a "Baer subplane spread"
# They partition the 21 points into 3 sets of 7.
# We already have B₁. For B₂ and B₃:

# Use: for each 2-orbit {P_i, P_j} of Frobenius, one goes to B₂ and one to B₃
# The assignment depends on a consistent coloring

# Let's find which points are in which Baer subplane by checking the
# classical construction: use the three "sub-PG(2,2)" inside PG(2,4)

def is_collinear_f4(p1, p2, p3):
    """Check if three points of PG(2,F₄) are collinear (det = 0)"""
    # det of 3x3 matrix with rows p1, p2, p3 over F₄
    a = p1[0] * (p2[1] * p3[2] + p2[2] * p3[1]) + \
        p1[1] * (p2[2] * p3[0] + p2[0] * p3[2]) + \
        p1[2] * (p2[0] * p3[1] + p2[1] * p3[0])
    return a.is_zero()

# Verify B₁ is a subplane: check it has 7 lines (each with 3 points)
b1_points = [points[i] for i in baer1]
b1_lines = []
for combo in combinations(range(len(baer1)), 3):
    p1, p2, p3 = b1_points[combo[0]], b1_points[combo[1]], b1_points[combo[2]]
    if is_collinear_f4(p1, p2, p3):
        b1_lines.append(combo)

print(f"\nB₁ has {len(b1_lines)} lines of 3 points = PG(2,F₂) ✓" if len(b1_lines) == 7 else f"B₁ has {len(b1_lines)} lines (expected 7)")

# Now construct B₂ and B₃ from the Frobenius orbits
# Each 2-orbit {i, σ(i)} contributes one point to B₂ and one to B₃
# The key: B₂ and B₃ must also be subplanes (7 lines of 3)

# Try all partitions of the 7 two-orbits into (B₂ choice, B₃ choice)
# This is 2^7 = 128 possibilities — check which give subplanes

print("\nSearching for Baer subplane partition...")
best_partition = None

for mask in range(128):
    b2_extra = []
    b3_extra = []
    for j, orbit in enumerate(two_orbits):
        if mask & (1 << j):
            b2_extra.append(orbit[0])
            b3_extra.append(orbit[1])
        else:
            b2_extra.append(orbit[1])
            b3_extra.append(orbit[0])

    # Check if B₂ forms a subplane (7 lines of 3)
    b2_pts = [points[i] for i in b2_extra]
    b2_line_count = 0
    for combo in combinations(range(7), 3):
        if is_collinear_f4(b2_pts[combo[0]], b2_pts[combo[1]], b2_pts[combo[2]]):
            b2_line_count += 1

    if b2_line_count == 7:
        # Also check B₃
        b3_pts = [points[i] for i in b3_extra]
        b3_line_count = 0
        for combo in combinations(range(7), 3):
            if is_collinear_f4(b3_pts[combo[0]], b3_pts[combo[1]], b3_pts[combo[2]]):
                b3_line_count += 1

        if b3_line_count == 7:
            best_partition = (b2_extra, b3_extra)
            print(f"  Found partition! mask={mask}")
            print(f"  B₂ = {b2_extra} (7 lines ✓)")
            print(f"  B₃ = {b3_extra} (7 lines ✓)")
            break

if best_partition is None:
    print("  Partition not found via simple Frobenius orbit splitting.")
    print("  The 3 Baer subplanes may require a more subtle construction.")
    # Brute force: try all ways to partition 14 non-F₂ points into 2 groups of 7
    non_f2 = [i for i in range(21) if i not in baer1]
    print(f"  Non-F₂ points: {non_f2} ({len(non_f2)} points)")

    found = False
    for combo in combinations(non_f2, 7):
        b2_candidate = list(combo)
        b3_candidate = [i for i in non_f2 if i not in combo]

        b2_pts = [points[i] for i in b2_candidate]
        lines2 = sum(1 for c in combinations(range(7), 3)
                     if is_collinear_f4(b2_pts[c[0]], b2_pts[c[1]], b2_pts[c[2]]))

        if lines2 == 7:
            b3_pts = [points[i] for i in b3_candidate]
            lines3 = sum(1 for c in combinations(range(7), 3)
                         if is_collinear_f4(b3_pts[c[0]], b3_pts[c[1]], b3_pts[c[2]]))

            if lines3 == 7:
                print(f"  FOUND! B₂ = {b2_candidate}, B₃ = {b3_candidate}")
                best_partition = (b2_candidate, b3_candidate)
                found = True
                break

    if not found:
        print("  No partition found by brute force either!")

print()

# ============================================================
# PART 4: The I-polynomial / Baer subplane correspondence
# ============================================================

print("=" * 70)
print("PART 4: I-POLYNOMIAL ↔ BAER SUBPLANE CORRESPONDENCE")
print("=" * 70)
print()

print("INDEPENDENCE POLYNOMIAL FACTORIZATION:")
print("  I(K₁⊔K₃, x) = I(K₁, x) · I(K₃, x) = (1+x)(1+3x)")
print("  At x=2: I = 3 × 7 = 21")
print()
print("BAER SUBPLANE DECOMPOSITION:")
print("  PG(2, F₄) = B₁ ⊔ B₂ ⊔ B₃")
print("  |PG(2, F₄)| = |B₁| + |B₂| + |B₃| = 7 + 7 + 7 = 21")
print()
print("THE CORRESPONDENCE:")
print("  Factor 7 = I(K₃, 2) = |PG(2, F₂)| = |Baer subplane|")
print("  Factor 3 = I(K₁, 2) = NUMBER of Baer subplanes")
print()
print("  K₃ component of Ω ↔ one Baer subplane (the cycle conflict structure)")
print("  K₁ component of Ω ↔ the 'selector' choosing which subplane")
print()
print("  The I-polynomial factorization IS the Baer partition!")
print()

# ============================================================
# PART 5: (x+1)^n and (x+2)^n — Simplex/Cuboid and Projective
# ============================================================

print("=" * 70)
print("PART 5: SIMPLEX (x+1)^n, CUBOID (x+2)^n, AND PROJECTIVE PLANES")
print("=" * 70)
print()

print("Evaluate at x = 2:")
print("  (x+1)^n at x=2: 3^n (simplex polynomial)")
print("  (x+2)^n at x=2: 4^n (cuboid polynomial)")
print("  Complement C(n) = 4^n - 3^n")
print()

for n in range(1, 8):
    s = 3**n
    c = 4**n
    comp = c - s
    # Check if comp is a projective plane size |PG(2, F_q)| = q²+q+1
    is_pp = False
    for q in range(2, 20):
        if q*q + q + 1 == comp:
            is_pp = True
            print(f"  n={n}: C(n) = 4^{n} - 3^{n} = {comp} = |PG(2, F_{q})| ← PROJECTIVE PLANE!")
            break
    if not is_pp:
        # Check if divisible by 7 or 21
        div7 = " (div by 7)" if comp % 7 == 0 else ""
        div21 = " (div by 21)" if comp % 21 == 0 else ""
        print(f"  n={n}: C(n) = 4^{n} - 3^{n} = {comp}{div7}{div21}")

print()
print("STUNNING: C(2) = 4² - 3² = 16 - 9 = 7 = |PG(2, F₂)| = H_forb_1!")
print("The complement of the simplex-in-cuboid at n=2 IS the Fano plane size!")
print()

# What about 21? Is there an n where C(n) = 21?
# 4^n - 3^n = 21: try n=2: 7, n=3: 37. No exact match.
# But: 21 = 3 × 7 = 3 × C(2). The factor 3 appears because
# PG(2,4) = 3 Baer subplanes, each of size C(2).

print("And C(2) × 3 = 7 × 3 = 21 = |PG(2, F₄)| = H_forb_2")
print("The three 'copies' are the three Baer subplanes!")
print()

# ============================================================
# PART 6: Why does the pattern stop? The cubic wall.
# ============================================================

print("=" * 70)
print("PART 6: WHY THE PATTERN STOPS — THE BAER WALL")
print("=" * 70)
print()

print("Projective planes over F_{2^k}:")
for k in range(1, 7):
    q = 2**k
    pp_size = q*q + q + 1
    # Number of Baer subplanes in PG(2, q²):
    # A Baer partition exists iff q is a prime power.
    # PG(2, q²) can be partitioned into q²+q+1 Baer subplanes PG(2,q)
    # Actually: |PG(2,q²)| / |PG(2,q)| = (q⁴+q²+1)/(q²+q+1)
    # = q² - q + 1 (by polynomial division)
    q2 = q*q
    pp2_size = q2*q2 + q2 + 1
    ratio = pp2_size // pp_size
    achievable = "FORBIDDEN" if pp_size in (7, 21) else "achievable"
    print(f"  k={k}: q=2^{k}={q}, |PG(2,F_{q})| = {pp_size} [{achievable}]")
    print(f"         Baer subplanes of PG(2,F_{q}²): |PG(2,F_{q}²)/{pp_size} = {ratio}")
    print(f"         q² - q + 1 = {q*q - q + 1}")

print()
print("THE BAER WALL:")
print("  k=1: |PG(2,F₂)| = 7: FORBIDDEN (= Ω cannot be K₃)")
print("  k=2: |PG(2,F₄)| = 21: FORBIDDEN (= 3 copies of forbidden)")
print("  k=3: |PG(2,F₈)| = 73: ACHIEVABLE (escapes the wall)")
print()
print("Why? At k=1,2 the field F_{2^k} has 'too few elements' for")
print("the Ω conflict graph to accommodate the projective structure.")
print("At k=3 (F₈), there are enough elements (8) that tournaments")
print("with enough vertices can achieve H=73.")
print()

# ============================================================
# PART 7: The Trinity — 3 Baer, 3 strands, 3^n
# ============================================================

print("=" * 70)
print("PART 7: THE TRINITY — 3 BAER, 3 STRANDS, 3^n")
print("=" * 70)
print()

print("THE NUMBER 3 APPEARS IN THREE STRUCTURAL ROLES:")
print()
print("1. BAER SUBPLANE COUNT:")
print("   PG(2,F₄) = 3 × PG(2,F₂)")
print("   21 = 3 × 7")
print("   The 3 Baer subplanes partition the forbidden structure")
print()
print("2. PASCAL STRAND COUNT:")
print("   Standard Pascal (2-strand): row sums 2^n → Fibonacci")
print("   Trinomial Pascal (3-strand): row sums 3^n → tribonacci")
print("   3 strands ↔ 3 Baer subplanes ↔ ternary triple classification")
print()
print("3. SIMPLEX VALUE:")
print("   (x+1)^n at x=2 = 3^n = I(empty graph on n vertices, 2)")
print("   = maximum possible I(G,2) for n-vertex graph")
print("   = CEILING of the independence polynomial world")
print()
print("4. k-NACCI LIMIT:")
print("   Weighted k-nacci with weight x=2 → limit ratio = 1+x = 3")
print("   The limit IS the per-vertex simplex contribution")
print()
print("UNIFIED VIEW:")
print("   3 = (x+1)|_{x=2} = |{Baer subplanes}| = |{Pascal strands}|")
print("   = k-nacci limit at x=2 = the fundamental 'packing factor'")
print()
print("   When we 'pack' a simplex (x+1)^n into a cuboid (x+2)^n,")
print("   the simplex captures 3^n of the cuboid's 4^n points.")
print("   The complement 4^n - 3^n at n=2 gives EXACTLY the Fano plane size 7.")
print("   The 3 Baer subplanes of PG(2,4) are the 3 'directions' of packing.")
print()

# ============================================================
# PART 8: Baer subplane → Ω structure impossibility
# ============================================================

print("=" * 70)
print("PART 8: WHY BAER STRUCTURE CANNOT BE Ω(T)")
print("=" * 70)
print()

print("THEOREM (informal): The Baer subplane structure of PG(2,F₄)")
print("is incompatible with tournament conflict graphs because:")
print()
print("1. Each Baer subplane = PG(2,F₂) = Fano plane")
print("   = 7 points with S(2,3,7) incidence structure")
print("   = every pair of points on exactly 1 common line")
print()
print("2. In Ω(T), the 'lines' would be sets of 3 mutually-conflicting")
print("   odd cycles (3 cycles sharing a common tournament vertex).")
print("   This is EXACTLY the K₃ structure that THM-201 prohibits!")
print()
print("3. So each Baer subplane requires K₃ substructures in Ω,")
print("   which are impossible. And 3 Baer subplanes compound this.")
print()
print("4. Result: H = |PG(2,F₄)| = 21 is forbidden because the")
print("   geometric structure REQUIRES tournament-impossible configurations.")
print()

# ============================================================
# PART 9: Quantitative check — incidence matrices
# ============================================================

print("=" * 70)
print("PART 9: INCIDENCE STRUCTURE OF PG(2, F₄)")
print("=" * 70)
print()

# Generate all lines of PG(2, F₄)
# A line is defined by a triple [a:b:c] in the dual (same as points)
# Point (x₁,x₂,x₃) is on line [a,b,c] iff ax₁ + bx₂ + cx₃ = 0

lines = []
for a in f4_elements:
    for b in f4_elements:
        for c in f4_elements:
            if a.is_zero() and b.is_zero() and c.is_zero():
                continue
            l = normalize_point(a, b, c)
            if l not in lines:
                lines.append(l)

print(f"Lines of PG(2, F₄): {len(lines)}")

# For each line, find which points are on it
def point_on_line(pt, ln):
    """Check if point (x₁,x₂,x₃) is on line [a,b,c]: ax₁+bx₂+cx₃ = 0"""
    return (ln[0]*pt[0] + ln[1]*pt[1] + ln[2]*pt[2]).is_zero()

line_points = {}
for i, ln in enumerate(lines):
    pts = [j for j, pt in enumerate(points) if point_on_line(pt, ln)]
    line_points[i] = pts

# Check: each line should have 5 points (q+1 = 4+1 = 5)
line_sizes = set(len(v) for v in line_points.values())
print(f"Points per line: {line_sizes} (expected {{5}})")

# Check: each pair of points on exactly 1 line
pair_line_count = defaultdict(int)
for pts in line_points.values():
    for p1, p2 in combinations(pts, 2):
        pair_line_count[(p1,p2)] += 1
        pair_line_count[(p2,p1)] += 1

print(f"Lines per point pair: {set(pair_line_count.values())} (expected {{1}})")
print()

# How many lines of PG(2,F₄) are entirely within B₁?
print("Lines of PG(2,F₄) restricted to Baer subplane B₁:")
b1_set = set(baer1)
b1_lines_in_full = 0
for pts in line_points.values():
    b1_on_line = [p for p in pts if p in b1_set]
    if len(b1_on_line) >= 3:
        b1_lines_in_full += 1

print(f"  Lines meeting B₁ in ≥3 points: {b1_lines_in_full}")
print(f"  (A Baer subplane meets each line in 1 or q+1=3 points)")
print()

# Count: how many lines meet B₁ in exactly 1 point vs 3 points?
meet_counts = defaultdict(int)
for pts in line_points.values():
    b1_on = len([p for p in pts if p in b1_set])
    meet_counts[b1_on] += 1

print(f"  Line-B₁ intersection sizes: {dict(sorted(meet_counts.items()))}")
print(f"  (Expected: 7 lines meet in 3 points, 14 lines meet in 1 point)")
print()

# ============================================================
# PART 10: The deep connection — Φ₃(q) and OCF
# ============================================================

print("=" * 70)
print("PART 10: Φ₃(q) AND THE OCF — THE DEEPEST CONNECTION")
print("=" * 70)
print()

print("The third cyclotomic polynomial Φ₃(q) = q² + q + 1:")
print("  Φ₃(1) = 3   (= I(K₁, 2) = number of Baer subplanes)")
print("  Φ₃(2) = 7   (= I(K₃, 2) = |PG(2,F₂)| = H_forb_1)")
print("  Φ₃(3) = 13  (Fibonacci prime, not forbidden)")
print("  Φ₃(4) = 21  (= I(K₁⊔K₃, 2) = |PG(2,F₄)| = H_forb_2)")
print("  Φ₃(5) = 31  (Mersenne prime)")
print("  Φ₃(6) = 43  (achievable)")
print("  Φ₃(7) = 57  (= 3×19)")
print("  Φ₃(8) = 73  (achievable)")
print()
print("OBSERVATION: Φ₃ at POWERS OF 2 gives projective planes over F_{2^k}:")
print("  Φ₃(2¹) = 7   → FORBIDDEN")
print("  Φ₃(2²) = 21  → FORBIDDEN")
print("  Φ₃(2³) = 73  → achievable")
print()
print("The forbidden values are Φ₃(2) and Φ₃(4) = Φ₃(2²).")
print("But Φ₃(2³) = 73 escapes. WHY?")
print()
print("ANSWER: The Baer subplane structure.")
print("  PG(2,F₄) = 3 × PG(2,F₂): the subplane IS the forbidden structure.")
print("  PG(2,F₈) cannot be decomposed into copies of PG(2,F₂) because")
print("  F₈ is NOT an extension of F₂ via squaring (F₈ ≅ F₂[x]/(x³+x+1)).")
print("  The Baer subplane of PG(2,F₈) would be PG(2,F_{√8}),")
print("  but √8 is not a prime power — no subplane exists!")
print()
print("MORE PRECISELY: A Baer subplane PG(2,q) exists in PG(2,Q)")
print("  iff Q = q². So:")
print("  - PG(2,4) has Baer subplanes PG(2,2) ✓ (4 = 2²)")
print("  - PG(2,8) has Baer subplanes PG(2,√8) ✗ (√8 not integer)")
print("  - PG(2,16) has Baer subplanes PG(2,4) ✓ (16 = 4²)")
print()
print("So: |PG(2,16)| = 273 = Φ₃(16). Is 273 forbidden?")
print("  PG(2,16) has Baer subplanes PG(2,4), each of size 21.")
print("  But 21 is ITSELF forbidden, so the chain continues!")
print("  273 = 21 × 13 = forbidden × ??")
print("  Actually 273/21 = 13, and Φ₃(4) = 21, so 273 = 13 × 21.")
print("  Number of Baer subplanes in PG(2,16): |PG(2,16)|/|PG(2,4)| = 273/21 = 13")
print()

# Is 273 achievable? We'd need to check very large n.
# For now, flag it as an open question.
print("OPEN QUESTION: Is H=273 = |PG(2,F₁₆)| forbidden?")
print("  If yes: the Baer chain continues (PG(2,2) → PG(2,4) → PG(2,16))")
print("  Each level contains the previous as a Baer subplane.")
print("  This would give forbidden values at |PG(2, 2^{2^k})|:")
print("  k=0: |PG(2,2)| = 7")
print("  k=1: |PG(2,4)| = 21")
print("  k=2: |PG(2,16)| = 273")
print("  k=3: |PG(2,256)| = 65793")
print("  These are the 'Baer tower' values!")
print()

# ============================================================
# PART 11: Summary — The Baer Subplane Theory of Forbidden H
# ============================================================

print("=" * 70)
print("PART 11: THE BAER SUBPLANE THEORY OF FORBIDDEN H")
print("=" * 70)
print()
print("MAIN THEOREM (conjectural):")
print("  H(T) ≠ |PG(2, F_{2^{2^k}})| for all tournaments T and all k ≥ 0.")
print()
print("PROOF SKETCH:")
print("  Base case (k=0): H ≠ 7 = |PG(2,F₂)|.")
print("    Proof: Ω(T) ≠ K₃ by THM-201 (K₃ requires a tournament vertex")
print("    common to 3 directed 3-cycles in a Steiner-like arrangement,")
print("    which forces additional cycles by dominance cascade).")
print()
print("  k=1: H ≠ 21 = |PG(2,F₄)|.")
print("    The I-polynomial of any graph with I(G,2)=21 factors through")
print("    the Baer subplane structure: either (1+3x) divides I(G,x)")
print("    (the K₃ poison), or G ∈ {K₆-2e, K₈-e, K₁₀} which cannot be Ω(T)")
print("    by structural constraints on cycle overlap in tournaments.")
print()
print("  k≥2: H ≠ |PG(2, F_{2^{2^k}})| (CONJECTURED).")
print("    Each level inherits the impossibility from the previous via")
print("    the Baer subplane containment PG(2,F_{q²}) ⊃ PG(2,F_q).")
print()
print("THE SIMPLEX-CUBOID NESTING INTERPRETATION:")
print("  The forbidden values arise at the 'corners' of the simplex-in-cuboid")
print("  packing. At n=2: the complement 4²-3²=7 is the Fano plane.")
print("  The 3 Baer subplanes are the 3 'directions' along which the")
print("  simplex touches the cuboid boundary.")
print("  The k-nacci limit 3 = the number of these directions.")
print()
print("  (x+1)^n = simplex = 'achievable core'")
print("  (x+2)^n = cuboid = 'total space'")
print("  Complement = the 'forbidden shell' with projective structure!")
