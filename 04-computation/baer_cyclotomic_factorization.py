#!/usr/bin/env python3
"""
baer_cyclotomic_factorization.py
opus-2026-03-14-S71k

THE CYCLOTOMIC-BAER-TOURNAMENT DEEP STRUCTURE

Key insight: The Baer subplane decomposition of PG(2,q^2) = Phi_3(q^2) points
factors as Phi_3(q^2) = Phi_3(q) * Phi_6(q), where:
  - Phi_3(q) = q^2+q+1 = size of each Baer subplane PG(2,q)
  - Phi_6(q) = q^2-q+1 = NUMBER of Baer subplanes in the partition

At q=2: 21 = 7 * 3 exactly mirrors 3-cycle generator * Fano plane.

This script explores:
1. The cyclotomic factorization Phi_3(q^2) = Phi_3(q)*Phi_6(q)
2. Why period 6 = LCM(2,3) creates EXACTLY the Baer structure
3. The Frobenius involution parallel: F_4/F_2 ↔ tournament complement
4. Baer subplane intersection patterns and tournament conflict graphs
5. The tower Phi_3(2), Phi_3(4), Phi_3(16), ... and forbidden value termination
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter
from math import gcd

print("=" * 70)
print("PART 1: THE CYCLOTOMIC BAER FACTORIZATION")
print("=" * 70)
print()

def phi3(x):
    return x**2 + x + 1

def phi6(x):
    return x**2 - x + 1

def phi1(x):
    return x - 1

def phi2(x):
    return x + 1

print("The master identity: Phi_3(x^2) = Phi_3(x) * Phi_6(x)")
print()
for q in range(2, 10):
    lhs = phi3(q**2)
    rhs = phi3(q) * phi6(q)
    print(f"  q={q}: Phi_3({q}^2) = Phi_3({q**2}) = {lhs}")
    print(f"       = Phi_3({q}) * Phi_6({q}) = {phi3(q)} * {phi6(q)} = {rhs}  {'✓' if lhs == rhs else '✗'}")
    print(f"       Baer decomposition: PG(2,F_{q**2}) = {phi6(q)} copies of PG(2,F_{q})")
    print()

print("THE TOURNAMENT INSTANCE (q=2):")
print(f"  Phi_3(4) = 21 = Phi_3(2) * Phi_6(2) = 7 * 3")
print(f"  PG(2,F_4) = 3 copies of PG(2,F_2) (Fano planes)")
print(f"  Phi_3(2) = 7 = Fano plane = K₃ poison = FORBIDDEN H-value 1")
print(f"  Phi_6(2) = 3 = cycle generator = NUMBER of Baer subplanes")
print(f"  Their product 21 = FORBIDDEN H-value 2")
print()

print("=" * 70)
print("PART 2: THE SIXTH CYCLOTOMIC STRUCTURE")
print("=" * 70)
print()

print("x^6 - 1 = Phi_1(x) * Phi_2(x) * Phi_3(x) * Phi_6(x)")
print()
for x in range(2, 8):
    p1, p2, p3, p6 = phi1(x), phi2(x), phi3(x), phi6(x)
    product = p1 * p2 * p3 * p6
    actual = x**6 - 1
    print(f"  x={x}: {x}^6-1 = {actual}")
    print(f"       = Phi_1({x}) * Phi_2({x}) * Phi_3({x}) * Phi_6({x})")
    print(f"       = {p1} * {p2} * {p3} * {p6} = {product}  {'✓' if product == actual else '✗'}")

print()
print("At x=2:")
print(f"  2^6 - 1 = 63 = 1 * 3 * 7 * 3")
print(f"  But 1*3*7*3 = 63 ✓")
print(f"  Phi_3(2) * Phi_6(2) = 7 * 3 = 21 = (2^6-1)/(2^2-1) = 63/3 = 21 ✓")
print()

print("DEEPER: x^4 + x^2 + 1 = Phi_3(x) * Phi_6(x)")
print("This polynomial is the 'Baer polynomial' — it generates the subplane partition.")
for x in range(2, 8):
    lhs = x**4 + x**2 + 1
    rhs = phi3(x) * phi6(x)
    print(f"  x={x}: {x}^4+{x}^2+1 = {lhs} = {phi3(x)}*{phi6(x)} = {rhs}  {'✓' if lhs==rhs else '✗'}")

print()
print("=" * 70)
print("PART 3: THE FROBENIUS-COMPLEMENT PARALLEL")
print("=" * 70)
print()

print("Two order-2 involutions in parallel algebraic structures:")
print()
print("  PROJECTIVE GEOMETRY            TOURNAMENT THEORY")
print("  ─────────────────────          ──────────────────")
print("  F_4 / F_2                      Tournament space T")
print("  Frobenius σ: x → x^2           Complement τ: T → T^op")
print("  Fixed field: F_2               Fixed set: SC tournaments")
print("  σ has order 2                  τ has order 2")
print("  F_4 = F_2[α] / (α²+α+1)      |T| = 2^C(n,2) points")
print("  PG(2,F_4) = 21 pts            H-spectrum = {odd integers}")
print("  Fixed pts = PG(2,F_2) = 7     SC: H(T) = H(T^op)")
print()

# Build F_4 and its Frobenius
print("F_4 = {0, 1, α, α+1} where α² = α+1 (so α²+α+1=0 in F_2)")
print("Frobenius σ(x) = x²:")
f4_elements = ['0', '1', 'α', 'α+1']
# In F_4: 0²=0, 1²=1, α²=α+1, (α+1)²=α²+1=α+1+1=α
f4_frobenius = {'0': '0', '1': '1', 'α': 'α+1', 'α+1': 'α'}
for x in f4_elements:
    print(f"  σ({x}) = {f4_frobenius[x]}  {'(fixed)' if x == f4_frobenius[x] else '(moved)'}")

print()
print("Fixed points of Frobenius = F_2 = {0, 1}")
print("Orbits of non-fixed points: {α, α+1}")
print()
print("On PG(2,F_4):")
print(f"  21 points total")
print(f"  7 fixed points = PG(2,F_2) (one Baer subplane)")
print(f"  14 moved points = 7 orbits of size 2")
print(f"  These 14 points form 2 more Baer subplanes (each of size 7)")
print()

# Self-complementary tournaments parallel
print("On tournament space (n=5):")
print(f"  C(5,2) = 10 arcs, 2^10 = 1024 tournaments")
print(f"  Complement τ(T) reverses all 10 arcs")
print(f"  Fixed points: 0 (no tournament can equal its complement)")
print(f"  Instead: SC tournaments satisfy T ≅ T^op (isomorphic, not equal)")
print(f"  8 SC isomorphism classes, 2 NSC pairs = 12 total classes")
print()

print("THE PARALLEL DEEPENS:")
print(f"  Frobenius: PG(2,F_4) → 3 Baer subplanes (partition)")
print(f"  Complement: tournament iso classes → SC + NSC pairs")
print(f"  At n=5: 12 classes = 8 + 2*2 = 8 SC + 4 NSC (in 2 pairs)")
print(f"  At n=7: 456 classes = 88 SC + 184*2 NSC (in 184 pairs)")
print()

print("=" * 70)
print("PART 4: BAER INTERSECTION PATTERNS → CONFLICT GRAPHS")
print("=" * 70)
print()

print("In PG(2,F_4), the 3 Baer subplanes B_0, B_1, B_2 partition the 21 points.")
print("But the LINES of PG(2,F_4) intersect the Baer subplanes in structured ways.")
print()

# Construct PG(2,F_4) explicitly
# F_4 = {0, 1, a, a+1} where a^2 = a+1
# Elements as integers: 0, 1, 2(=a), 3(=a+1)
class F4:
    """Finite field with 4 elements."""
    def __init__(self, val):
        self.val = val % 4 if isinstance(val, int) else val
        # But F_4 addition is XOR, multiplication is polynomial mod x^2+x+1

# Actually, let's use a lookup table approach
# F_4 = GF(4) with elements {0, 1, a, a^2} where a^2+a+1=0
# Represent as integers 0,1,2,3 where 2=a, 3=a+1=a^2
add_table = [
    [0, 1, 2, 3],
    [1, 0, 3, 2],
    [2, 3, 0, 1],
    [3, 2, 1, 0]
]
mul_table = [
    [0, 0, 0, 0],
    [0, 1, 2, 3],
    [0, 2, 3, 1],
    [0, 3, 1, 2]
]

def f4_add(a, b):
    return add_table[a][b]

def f4_mul(a, b):
    return mul_table[a][b]

# Build PG(2,F_4): points are equivalence classes of (x,y,z) ≠ (0,0,0)
# Normalize: first nonzero coordinate = 1
def normalize_pg2(point):
    """Normalize a point in PG(2,F_4)."""
    x, y, z = point
    for coord in [x, y, z]:
        if coord != 0:
            # Divide all by this coordinate
            inv = {1: 1, 2: 3, 3: 2}[coord]  # multiplicative inverses in F_4
            return (f4_mul(x, inv), f4_mul(y, inv), f4_mul(z, inv))
    return None  # zero vector

# Generate all 21 points
points_pg2_f4 = set()
for x in range(4):
    for y in range(4):
        for z in range(4):
            if (x, y, z) != (0, 0, 0):
                p = normalize_pg2((x, y, z))
                if p:
                    points_pg2_f4.add(p)

points_pg2_f4 = sorted(points_pg2_f4)
print(f"PG(2,F_4) has {len(points_pg2_f4)} points ✓" if len(points_pg2_f4) == 21 else f"ERROR: {len(points_pg2_f4)} points")

# Identify the Baer subplane PG(2,F_2) = points with coordinates in {0,1}
baer_0 = sorted([p for p in points_pg2_f4 if all(c in [0, 1] for c in p)])
print(f"Baer subplane B_0 = PG(2,F_2): {len(baer_0)} points")

# Frobenius: (x,y,z) → (x^2, y^2, z^2)
def frobenius(point):
    return normalize_pg2(tuple(f4_mul(c, c) for c in point))

# Find orbits under Frobenius
fixed_pts = [p for p in points_pg2_f4 if frobenius(p) == p]
moved_pts = [p for p in points_pg2_f4 if frobenius(p) != p]
print(f"Fixed points of Frobenius: {len(fixed_pts)}")
print(f"Moved points: {len(moved_pts)} = {len(moved_pts)//2} orbits of size 2")

# The fixed points should be exactly B_0
assert set(fixed_pts) == set(baer_0), "Fixed points ≠ Baer subplane!"
print(f"Fixed points = B_0 ✓")
print()

# Construct the other two Baer subplanes
# They consist of the moved points, partitioned by Frobenius orbits
# Each orbit {p, σ(p)} contributes one point to B_1 and one to B_2
orbits = []
seen = set()
for p in moved_pts:
    if p not in seen:
        fp = frobenius(p)
        orbits.append((p, fp))
        seen.add(p)
        seen.add(fp)

print(f"Frobenius orbits of moved points: {len(orbits)}")

# The two non-trivial Baer subplanes: each selects one element from each orbit
# But which selection gives a subplane? This requires checking the subplane axiom.
# Actually, the Baer subplanes come from the cosets of the Frobenius fixed field.
# B_1 = {(x, y, z) : x^2=y, y^2=z, z^2=x equivalent to Frobenius orbit structure}

# Alternative: use the F_4 automorphism group structure
# The 3 Baer subplanes correspond to the 3 "spread lines" of a specific spread

# Let me find them by checking which 7-element subsets are subplanes
def is_line(p1, p2, p3):
    """Check if three points in PG(2,F_4) are collinear."""
    # Three points are collinear if det = 0
    # det = x1(y2*z3 - y3*z2) - y1(x2*z3 - x3*z2) + z1(x2*y3 - x3*y2)
    x1,y1,z1 = p1
    x2,y2,z2 = p2
    x3,y3,z3 = p3

    t1 = f4_mul(y2, z3)
    t2 = f4_mul(y3, z2)
    term1 = f4_mul(x1, f4_add(t1, t2))

    t3 = f4_mul(x2, z3)
    t4 = f4_mul(x3, z2)
    term2 = f4_mul(y1, f4_add(t3, t4))

    t5 = f4_mul(x2, y3)
    t6 = f4_mul(x3, y2)
    term3 = f4_mul(z1, f4_add(t5, t6))

    det = f4_add(f4_add(term1, term2), term3)
    return det == 0

# A Baer subplane is a set of 7 points such that:
# - Every line of PG(2,F_4) meets it in 1 or 3 points (q+1 = 3 for q=2)
# - It contains 7 lines of 3 points each

# Find lines of PG(2,F_4)
lines_pg2_f4 = []
for i, j in combinations(range(21), 2):
    line_pts = [points_pg2_f4[i], points_pg2_f4[j]]
    for k in range(j+1, 21):
        if is_line(points_pg2_f4[i], points_pg2_f4[j], points_pg2_f4[k]):
            line_pts.append(points_pg2_f4[k])
    if len(line_pts) == 5:  # lines in PG(2,4) have 5 points
        line_set = frozenset([points_pg2_f4.index(p) for p in line_pts])
        if line_set not in [frozenset(l) for l in lines_pg2_f4]:
            lines_pg2_f4.append(sorted(line_set))

# Actually, each line has q+1 = 5 points. Let me find lines properly.
lines = []
for i in range(21):
    for j in range(i+1, 21):
        # Find all points on the line through points i and j
        line = {i, j}
        for k in range(21):
            if k != i and k != j:
                if is_line(points_pg2_f4[i], points_pg2_f4[j], points_pg2_f4[k]):
                    line.add(k)
        line = frozenset(line)
        if line not in [frozenset(l) for l in lines]:
            lines.append(sorted(line))

print(f"Lines of PG(2,F_4): {len(lines)}")
print(f"Points per line: {set(len(l) for l in lines)}")

# Check intersection of lines with B_0
baer_0_indices = set(points_pg2_f4.index(p) for p in baer_0)
print(f"\nB_0 indices: {sorted(baer_0_indices)}")

intersection_sizes = Counter()
for line in lines:
    inter = len(set(line) & baer_0_indices)
    intersection_sizes[inter] += 1

print(f"Line-Baer intersection sizes:")
for size, count in sorted(intersection_sizes.items()):
    print(f"  |line ∩ B_0| = {size}: {count} lines")

# For a Baer subplane, each line meets it in 1 or q+1 = 3 points
# Lines meeting in 3 points are "Baer lines" (internal to the subplane)
# Lines meeting in 1 point are "external" lines

print()
print("BAER SUBPLANE INTERSECTION THEOREM:")
print("  Every line of PG(2,q^2) meets a Baer subplane in either 1 or q+1 points.")
print(f"  At q=2: every line meets B_0 in 1 or 3 points")
print(f"  Lines with 3-point intersection: {intersection_sizes.get(3,0)} = 7 = lines of PG(2,F_2)")
print(f"  Lines with 1-point intersection: {intersection_sizes.get(1,0)} = 14 = remaining lines")
print()

# Find ALL 3 Baer subplanes
print("Finding all Baer subplane partitions of PG(2,F_4)...")

# A Baer subplane has 7 points, every pair of its points determines a unique line
# that meets the subplane in exactly 3 points (for q+1=3 lines within the subplane)

def check_baer(point_indices, all_lines):
    """Check if a set of 7 points forms a Baer subplane."""
    if len(point_indices) != 7:
        return False
    pset = set(point_indices)

    # Check: every line meets it in 1 or 3 points
    for line in all_lines:
        inter = len(set(line) & pset)
        if inter not in [1, 3]:
            return False

    # Check: C(7,2) = 21 pairs, each pair on a unique line
    # and there should be 7 internal lines (meeting in 3)
    internal_lines = 0
    for line in all_lines:
        if len(set(line) & pset) == 3:
            internal_lines += 1

    return internal_lines == 7

# We already know B_0. Find complementary Baer subplanes.
remaining = set(range(21)) - baer_0_indices

# Try to find a Baer subplane among the remaining 14 points
from itertools import combinations as combs
baer_subplanes = [sorted(baer_0_indices)]
found_second = False

for subset in combs(sorted(remaining), 7):
    if check_baer(set(subset), lines):
        baer_subplanes.append(sorted(subset))
        complement = sorted(remaining - set(subset))
        if check_baer(set(complement), lines):
            baer_subplanes.append(complement)
        found_second = True
        break

if found_second:
    print(f"Found {len(baer_subplanes)} Baer subplanes forming a partition!")
    for i, bp in enumerate(baer_subplanes):
        print(f"  B_{i}: {bp}")
        pts = [points_pg2_f4[j] for j in bp]
        f2_pts = all(all(c in [0,1] for c in p) for p in pts)
        if f2_pts:
            print(f"       = PG(2,F_2) (coordinates in F_2)")
        else:
            print(f"       (coordinates involve α or α+1)")
else:
    print("Partition not found by brute force; trying alternative construction...")

print()
print("=" * 70)
print("PART 5: THE TOURNAMENT-BAER CONFLICT ANALOGY")
print("=" * 70)
print()

print("CONFLICT GRAPH Ω(T) has:")
print("  Vertices = odd cycles of T")
print("  Edges = cycles sharing a tournament vertex")
print()
print("BAER INCIDENCE has:")
print("  'Vertices' = Baer subplanes B_i (i=0,1,2)")
print("  'Edges' = lines connecting points in different subplanes")
print()
print("KEY ANALOGY:")
print("  In Ω(T): I(Ω(T), 2) = H(T) counts Hamiltonian paths")
print("  In Baer: 21 = 7 * 3 = I(K₃, 2) * (number of subplanes)")
print()
print("  The K₃ structure in Ω(T) IS the Baer partition structure!")
print("  Three 'subplanes' (K₃ vertices) each of size 7 (independence count)")
print("  K₃ is the UNIQUE graph with I(G,2) = 7")
print("  And 3 copies of K₃ IS K₃ itself (K₃ is self-replicating)")
print()

# Compute I(K₃, x) and the Baer connection
print("I(K₃, x) = 1 + 3x + 3x(x-1)/2 + ... let me compute properly")
print("K₃ = triangle. Independent sets: ∅, {1}, {2}, {3}, plus NO edges/triples")
print("Wait: K₃ has edges between ALL pairs, so max independent set size = 1")
print("I(K₃, x) = 1 + 3x")
print(f"I(K₃, 2) = 1 + 6 = 7 ✓")
print()

print("But the Baer partition 21 = 3 * 7:")
print("  Think of 21 as I(3K₁ + K₃-structure, 2)")
print("  Actually: I(3*K₁, x) = (1+x)^3 = 1+3x+3x²+x³")
print(f"  I(3*K₁, 2) = 3^3 = 27 ≠ 21")
print()
print("  Better: I(P₄, x) = (1+x)(1+3x)")
print(f"  I(P₄, 2) = 3 * 7 = 21 ✓")
print(f"  P₄ (path on 4 vertices) factors as: (1+x)(1+3x)")
print(f"  Factor (1+x) at x=2: 3 = Φ₆(2) = number of Baer subplanes")
print(f"  Factor (1+3x) at x=2: 7 = Φ₃(2) = Fano plane size")
print()

# The I-polynomial factorization mirrors the Baer partition!
print("THE I-POLYNOMIAL FACTORIZATION = BAER PARTITION:")
print("  I(P₄, x) = (1+x)(1+3x)")
print("  At x=2: 3 * 7 = 21")
print("  (1+x) = Φ₂(x) evaluated at x, = x+1")
print("  (1+3x) at x=2 gives 7 = Φ₃(2)")
print()
print("  But 1+3x = 3(x + 1/3) = 3Φ₃^{-1}(x)")
print("  The root x = -1/3 of (1+3x) is the K₃ POISON ROOT!")
print("  This root is 1/Φ₃(1) = 1/3 = the CONE RATIO!")
print()

print("=" * 70)
print("PART 6: TOWER OF BAER PLANES AND FORBIDDEN VALUES")
print("=" * 70)
print()

print("The Baer tower: PG(2,F_{2^{2^k}}) for k=0,1,2,...")
print()
for k in range(5):
    q = 2**(2**k)
    n_pts = phi3(q)
    q_sub = int(q**0.5) if k > 0 else None
    n_baer = phi6(q_sub) if k > 0 else "N/A"
    subplane_size = phi3(q_sub) if k > 0 else "N/A"
    forbidden = "YES" if k <= 1 else "NO (too many graph types)"
    print(f"  k={k}: q=2^{2**k}={q}")
    print(f"       |PG(2,F_{q})| = Φ₃({q}) = {n_pts}")
    if k > 0:
        print(f"       Baer partition: {n_baer} copies of PG(2,F_{q_sub})")
        print(f"       Each subplane has {subplane_size} points")
        print(f"       {n_pts} = {subplane_size} × {n_baer}")
    print(f"       Forbidden H value? {forbidden}")
    print()

print("WHY THE TOWER STOPS AT k=1:")
print(f"  k=0: Φ₃(2) = 7.    Only 1 graph type with I(G,2)=7 (K₃). BLOCKABLE.")
print(f"  k=1: Φ₃(4) = 21.   Only 6 graph types with I(G,2)=21. ALL BLOCKABLE.")
print(f"  k=2: Φ₃(16) = 273. Hundreds of types. NOT ALL BLOCKABLE → H=273 achievable.")
print()
print("  The blocking mechanism is the K₃ AVOIDANCE theorem (THM-029):")
print("  Ω(T) never contains K₃ as a component.")
print("  At I(G,2)=7: the ONLY graph is K₃ → always blocked")
print("  At I(G,2)=21: all 6 types contain K₃ substructure → blocked")
print("  At I(G,2)=273: some types avoid K₃ entirely → not blocked")

print()
print("=" * 70)
print("PART 7: THE CYCLOTOMIC CASCADE")
print("=" * 70)
print()

print("All cyclotomic polynomials evaluated at x=2:")
# Compute cyclotomic polynomials at x=2
# Phi_n(2) = product over d|n of (2^d - 1)^{mu(n/d)}
# Simpler: use the fact that x^n - 1 = product_{d|n} Phi_d(x)

def euler_phi(n):
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

def mobius(n):
    if n == 1:
        return 1
    factors = []
    temp = n
    p = 2
    while p * p <= temp:
        if temp % p == 0:
            count = 0
            while temp % p == 0:
                temp //= p
                count += 1
            if count > 1:
                return 0
            factors.append(p)
        p += 1
    if temp > 1:
        factors.append(temp)
    return (-1)**len(factors)

def cyclotomic_at_2(n):
    """Compute Phi_n(2) using the product formula."""
    # Phi_n(x) = product_{d|n} (x^d - 1)^{mu(n/d)}
    result_num = 1
    result_den = 1
    for d in range(1, n+1):
        if n % d == 0:
            m = mobius(n // d)
            val = 2**d - 1
            if m == 1:
                result_num *= val
            elif m == -1:
                result_den *= val
    return result_num // result_den

print(f"  n  | Phi_n(2) | deg(Phi_n) | Prime? | H-forbidden?")
print(f"  ---+----------+------------+--------+-------------")
forbidden_set = {7, 21}
for n in range(1, 25):
    val = cyclotomic_at_2(n)
    deg = euler_phi(n)
    is_prime = all(val % p != 0 for p in range(2, int(val**0.5)+1)) and val > 1
    is_forbidden = val in forbidden_set
    marker = " ← FORBIDDEN" if is_forbidden else ""
    marker += " ← CONE DENOM" if val == 3 else ""
    print(f"  {n:2d} | {val:8d} | {deg:10d} | {'YES' if is_prime else 'no ':3s}    | {marker}")

print()
print("OBSERVATIONS:")
print("  Φ₁(2) = 1: trivial")
print("  Φ₂(2) = 3: the cycle generator / simplex factor / Baer partition count")
print("  Φ₃(2) = 7: FORBIDDEN (Fano plane)")
print("  Φ₄(2) = 5: NOT forbidden (achievable)")
print("  Φ₅(2) = 31: NOT forbidden")
print("  Φ₆(2) = 3: same as Φ₂(2)! Both give 3.")
print()
print("  WHY Φ₃ IS SPECIAL:")
print("  - Φ₃(2) = 7 is the ONLY cyclotomic value that's forbidden")
print("  - Φ₃ is the minimal polynomial of F_4 over F_2")
print("  - The K₃ avoidance theorem makes I(K₃,2) = 7 unreachable")
print("  - K₃ is the UNIQUE graph with I(G,2) = Φ₃(2)")
print()

print("=" * 70)
print("PART 8: THE FIBONACCI-BAER BRIDGE VIA CYCLOTOMIC NORMS")
print("=" * 70)
print()

print("Recall from MEMORY: F_8 = 21 = F_4 × L_4 = 3 × 7")
print("And from Eisenstein norms (kind-pasteur S105c):")
print("  Φ₃(b) = |1 - b·ω²|² (norm of Eisenstein integer)")
print("  At b=2: Φ₃(2) = 7 = N(1 - 2ω²)")
print("  At b=4: Φ₃(4) = 21 = N(1 - 4ω²)")
print()

# Fibonacci and Baer
print("The Fibonacci-Lucas identity at index 4:")
print(f"  F_8 = F_4 × L_4 = 3 × 7 = 21")
print(f"  F_4 = 3 = Φ₆(2) = Φ₂(2) = number of Baer subplanes")
print(f"  L_4 = 7 = Φ₃(2) = Fano plane size")
print(f"  Index 4 = |F₄| = the field that creates the Baer structure!")
print()

print("THE TRIPLE FACTORIZATION OF 21:")
print(f"  21 = 3 × 7         (arithmetic)")
print(f"     = Φ₂(2) × Φ₃(2) (cyclotomic)")
print(f"     = F_4 × L_4     (Fibonacci-Lucas)")
print(f"     = Φ₃(4)         (projective plane)")
print(f"     = |PG(2,F_4)|   (geometry)")
print(f"     = F_8            (Fibonacci)")
print(f"     = I(P_4, 2)     (independence polynomial)")
print(f"     = H_forb_2      (tournament theory)")
print()

print("ALL SEVEN are the SAME number, seen from different mathematical windows.")
print()

# The Fibonacci-Eisenstein connection
print("DEEPER: The Fibonacci sequence mod 7 has period 16:")
fibs = [0, 1]
for i in range(50):
    fibs.append(fibs[-1] + fibs[-2])

fib_mod7 = [f % 7 for f in fibs[:20]]
print(f"  F_n mod 7: {fib_mod7}")
# Find period
for p in range(1, 50):
    if all(fibs[i] % 7 == fibs[i+p] % 7 for i in range(20)):
        print(f"  Period of Fibonacci mod 7 (Pisano period π(7)): {p}")
        break

print()
fib_mod3 = [f % 3 for f in fibs[:15]]
print(f"  F_n mod 3: {fib_mod3}")
for p in range(1, 50):
    if all(fibs[i] % 3 == fibs[i+p] % 3 for i in range(15)):
        print(f"  Period of Fibonacci mod 3 (Pisano period π(3)): {p}")
        break

print()
print("  π(7) = 16 = 2^4 = |F₄|²")
print("  π(3) = 8 = 2^3 = |F₄| * 2")
print("  LCM(π(3), π(7)) = LCM(8, 16) = 16")
print("  The Pisano periods at the Baer primes are powers of 2!")
print()

# Check: at which Fibonacci numbers does 21 divide?
print("Fibonacci numbers divisible by 21:")
for i in range(30):
    if fibs[i] % 21 == 0:
        print(f"  F_{i} = {fibs[i]}, 21 | F_{i}")
# Pisano period of 21
for p in range(1, 100):
    if all(fibs[i] % 21 == fibs[i+p] % 21 for i in range(30)):
        print(f"  Pisano period π(21) = {p}")
        break

print()
print("=" * 70)
print("PART 9: THE SIMPLEX-CUBOID-BAER NESTING")
print("=" * 70)
print()

print("Recall the user's observation:")
print("  'An equilateral triangle sits in a square with two halves on either side.'")
print("  'A tetrahedron sits in a cube with 4 halves around it.'")
print()
print("The volume ratios:")
for n in range(1, 8):
    simplex_vol = 1  # relative: n-simplex in n-cube
    # Volume of regular n-simplex inscribed in n-cube: sqrt(n+1)/(n! * 2^{n/2})
    # Ratio simplex/cube varies by dimension

    # Actually: for the standard simplex with vertices at unit vectors,
    # volume = 1/n! and the cube is [0,1]^n with volume 1
    # So ratio = 1/n!

    # But for regular simplex inscribed in [-1,1]^n cube:
    # At n=2: equilateral triangle in unit square, ratio = sqrt(3)/4 ≈ 0.433
    # At n=3: regular tetrahedron in unit cube, ratio = 1/3

    # The "complementary pieces" count:
    # n=2: triangle in square → 2 extra pieces (3 total pieces)
    # n=3: tetrahedron in cube → 4 extra pieces (5 total pieces)
    # Pattern: n+1 extra pieces?

    # Actually the user says:
    # n=2: 2 halves on either side = 2 pieces
    # n=3: 4 halves around it = 4 pieces

    extra = 2**(n) - (n+1) if n <= 3 else None
    print(f"  n={n}: {n}-simplex in {n}-cube")
    if n == 1:
        print(f"       Line segment in line segment: no extra pieces")
    elif n == 2:
        print(f"       Triangle in square: 2 extra pieces")
    elif n == 3:
        print(f"       Tetrahedron in cube: 4 extra pieces")
        print(f"       Ratio = 1/3 = 1/Φ₃(1) = CONE RATIO!")

print()

# Compute 2^n - (n+1) for various n
print("Extra pieces when inscribing n-simplex in n-cube: 2^n - (n+1)")
for n in range(1, 10):
    extra = 2**n - (n + 1)
    proj = phi3(n-1) if n >= 2 else "N/A"
    print(f"  n={n}: 2^{n} - {n+1} = {extra:5d}    |PG(2,F_{n-1})| = {proj}")

print()
print("  At n=3: extra = 4 = |F₄|")
print("  At n=7: extra = 120 = |S_5| = |Aut(icosahedron rotations)| = 5!")
print("  The simplex-cuboid gap 2^n-(n+1) connects to Baer at n=3")
print("  because |F₄| = 4 is the field that creates the Baer structure!")
print()

# The deeper connection: 2^n - (n+1) = sum_{k=2}^n C(n,k)
print("Note: 2^n - (n+1) = Σ C(n,k) for k=2..n")
print("  = number of subsets of size ≥ 2")
print("  = number of 'interaction terms' beyond individual elements")
print()
print("  In the Walsh decomposition, these correspond to")
print("  the Fourier coefficients at degree ≥ 2.")
print("  H depends on degree 2 and higher — the 'extra pieces'!")
print()

print("=" * 70)
print("PART 10: THE 1/3 GRAND UNIFICATION")
print("=" * 70)
print()

print("The number 1/3 appears in AT LEAST 8 independent contexts:")
print()
print("  1. Φ₃(1) = 3, so 1/3 = 1/Φ₃(1)")
print("     (reciprocal of cone polynomial)")
print()
print("  2. K₃ poison root: I(P₄, -1/3) = 0")
print("     (the root that kills path independence poly)")
print()
print("  3. Var(H)/Mean(H)² = 1/3 at n=3,4")
print("     (tournament variance ratio)")
print()
print("  4. Volume ratio: tetrahedron/cube = 1/3")
print("     (simplex-cuboid nesting at n=3)")
print()
print("  5. Φ₆(2) = 3 Baer subplanes in PG(2,F₄)")
print("     (Baer partition count)")
print()
print("  6. Parseval: 1/Σ|c_k|² = 1/3 for Φ₃(e^{iθ})")
print("     (unit circle energy of cyclotomic)")
print()
print("  7. Cuboid-simplex excess: (2^n - (1+x)^n) / 2^n → 1 - (3/4)^n")
print("     per-dimension factor: 1/(x+1) = 1/3 at x=2")
print()
print("  8. Eisenstein norm: Φ₃(1) = N(1-ω²) = 3")
print("     (smallest non-trivial Eisenstein prime norm)")
print()
print("  ALL EIGHT are manifestations of Φ₃(1) = 3.")
print("  The THIRD cyclotomic polynomial evaluated at 1 gives the")
print("  fundamental constant of tournament theory.")
print()

# Final synthesis
print("=" * 70)
print("GRAND SYNTHESIS: THE BAER-CYCLOTOMIC-TOURNAMENT TRIANGLE")
print("=" * 70)
print()
print("                    Φ₃(x) = x² + x + 1")
print("                    /         |         \\")
print("                   /          |          \\")
print("          Φ₃(1)=3       Φ₃(2)=7       Φ₃(4)=21")
print("         /    \\          |    |         |    \\")
print("    Var/μ²  3 Baer   Fano   K₃     PG(2,4)  F₈")
print("    =1/3    planes   plane  poison  Baer    Fibonacci")
print("                      |")
print("                   I(K₃,2)")
print("                      |")
print("                  THM-029")
print("                      |")
print("                  H ≠ 7")
print()
print("The single polynomial Φ₃ generates:")
print("  - The variance structure (x=1)")
print("  - The first forbidden value (x=2)")
print("  - The second forbidden value (x=4)")
print("  - The field F₄ = F₂[x]/Φ₃(x)")
print("  - The Baer partition via Φ₃(q²) = Φ₃(q)·Φ₆(q)")
print("  - The Eisenstein integers Z[ω] where ω is a root of Φ₃")
print()
print("=" * 70)
print("DONE — BAER CYCLOTOMIC FACTORIZATION")
print("=" * 70)
