#!/usr/bin/env python3
"""
exceptional_platonic_tournament.py — opus-2026-03-14-S77

The 5 exceptional simple Lie groups and the 5 Platonic solids:
How do they relate to tournament theory and the keys 2, 3?

PLATONIC SOLIDS:
  Tetrahedron:   4 vertices, 6 edges, 4 faces  (self-dual)
  Cube:          8 vertices, 12 edges, 6 faces  (dual: octahedron)
  Octahedron:    6 vertices, 12 edges, 8 faces  (dual: cube)
  Dodecahedron: 20 vertices, 30 edges, 12 faces (dual: icosahedron)
  Icosahedron:  12 vertices, 30 edges, 20 faces (dual: dodecahedron)

EXCEPTIONAL LIE GROUPS (5 of them):
  G₂: dim 14, rank 2   (automorphisms of octonions)
  F₄: dim 52, rank 4   (automorphisms of exceptional Jordan algebra)
  E₆: dim 78, rank 6
  E₇: dim 133, rank 7
  E₈: dim 248, rank 8  (the biggest exceptional)

THE CONNECTIONS:
1. Both "fivefold": exactly 5 Platonic solids, exactly 5 exceptional groups
2. The Platonic solids correspond to finite subgroups of SO(3)
3. The ADE classification connects Dynkin diagrams to both
4. McKay correspondence: finite subgroups of SU(2) ↔ ADE Dynkin diagrams

HOW DOES THIS RELATE TO TOURNAMENTS?
- Tournament theory lives in the type A root system (Weyl group = S_n)
- The keys 2, 3 are Cartan determinants of A₁, A₂
- The exceptional groups have different keys...

THIS SCRIPT EXPLORES:
1. Cartan determinants of all simple Lie algebras
2. Independence polynomial analogs for exceptional root systems
3. The ADE-tournament correspondence
4. Platonic solid symmetry groups as tournament constraints
"""

import math
from itertools import combinations, permutations

# ====================================================================
print("=" * 70)
print("PART 1: CARTAN DETERMINANTS — THE KEYS FOR EACH LIE ALGEBRA")
print("=" * 70)
print()
print("For type A_n: det = n+1")
print("  A₁: det=2, A₂: det=3, A₃: det=4, ..., A₈: det=9")
print("  THESE ARE OUR KEYS: 2 and 3 from A₁ and A₂")
print()

# Cartan matrices and determinants
def cartan_det_A(n):
    return n + 1

def cartan_det_B(n):
    return 2

def cartan_det_C(n):
    return 2

def cartan_det_D(n):
    return 4

# For exceptional:
exceptional_dets = {
    'G₂': 1,  # det of 2x2 Cartan matrix [[2,-1],[-3,2]] = 4-3 = 1
    'F₄': 1,
    'E₆': 3,
    'E₇': 2,
    'E₈': 1,
}

print("Cartan determinants of all simple Lie algebras:")
print()
print("Classical series:")
for n in range(1, 10):
    print(f"  A_{n}: det = {cartan_det_A(n)}")
print()
for n in range(2, 8):
    print(f"  B_{n}: det = {cartan_det_B(n)}")
print()
for n in range(3, 8):
    print(f"  C_{n}: det = {cartan_det_C(n)}")
print()
for n in range(4, 9):
    print(f"  D_{n}: det = {cartan_det_D(n)}")
print()
print("Exceptional:")
for name, det in exceptional_dets.items():
    print(f"  {name}: det = {det}")

print()
print("OBSERVATION: The Cartan determinant = |weight lattice / root lattice|")
print("  = index of root lattice in weight lattice")
print("  = order of the center of the simply-connected group")
print()
print("For A_n: det = n+1 = |Z(SU(n+1))| = |Z_{n+1}|")
print("  A₁: Z₂ (2 = KEY₁)")
print("  A₂: Z₃ (3 = KEY₂)")
print()
print("For the exceptionals:")
print("  G₂: trivial center (det=1)")
print("  F₄: trivial center (det=1)")
print("  E₆: Z₃ center (det=3 = KEY₂!)")
print("  E₇: Z₂ center (det=2 = KEY₁!)")
print("  E₈: trivial center (det=1)")
print()
print("REMARKABLE: E₆ has center Z₃ and E₇ has center Z₂!")
print("The two LARGEST non-trivial exceptionals carry our KEYS!")
print("E₆ ↔ KEY₂ = 3, E₇ ↔ KEY₁ = 2")

# ====================================================================
print()
print("=" * 70)
print("PART 2: PLATONIC SOLIDS AND THEIR SYMMETRY NUMBERS")
print("=" * 70)
print()

platonic = [
    ("Tetrahedron", 4, 6, 4, 12, "A₄", "S₄"),
    ("Cube",        8, 12, 6, 24, "S₄", "S₄"),
    ("Octahedron",  6, 12, 8, 24, "S₄", "S₄"),
    ("Dodecahedron",20, 30, 12, 60, "A₅", "A₅"),
    ("Icosahedron", 12, 30, 20, 60, "A₅", "A₅"),
]

print(f"{'Solid':<15} {'V':<4} {'E':<4} {'F':<4} {'|Rot|':<6} {'χ=V-E+F':<8} {'V·F':<6} {'Schläfli'}")
for name, V, E, F, rot, _, _ in platonic:
    chi = V - E + F
    schlafli = ""
    if name == "Tetrahedron": schlafli = "{3,3}"
    elif name == "Cube": schlafli = "{4,3}"
    elif name == "Octahedron": schlafli = "{3,4}"
    elif name == "Dodecahedron": schlafli = "{5,3}"
    elif name == "Icosahedron": schlafli = "{3,5}"
    print(f"{name:<15} {V:<4} {E:<4} {F:<4} {rot:<6} {chi:<8} {V*F:<6} {schlafli}")

print()
print("ALL have Euler characteristic χ = V - E + F = 2 (genus 0)")
print("This is the SPHERE relation — every Platonic solid is a 2-sphere.")
print()
print("KEY NUMBERS in Platonic solids:")
print(f"  Rotation groups: |A₄|=12, |S₄|=24, |A₅|=60")
print(f"  12 = 2²·3, 24 = 2³·3, 60 = 2²·3·5")
print(f"  ALL are divisible by 6 = 2·3 = KEY₁·KEY₂!")
print()
print("  Vertex counts: 4, 6, 8, 12, 20")
print(f"  4 = 2², 6 = 2·3, 8 = 2³, 12 = 2²·3, 20 = 2²·5")
print(f"  8 and 6 are our KEY POWERS: 8 = KEY₁³, 6 = KEY₁·KEY₂")
print()
print("  Edge counts: 6, 12, 30")
print(f"  6 = KEY₁·KEY₂, 12 = 2·KEY₁·KEY₂, 30 = 2·3·5")

# ====================================================================
print()
print("=" * 70)
print("PART 3: THE ADE CLASSIFICATION")
print("=" * 70)
print()
print("The ADE classification connects:")
print("  1. Simply-laced Dynkin diagrams (A_n, D_n, E_6, E_7, E_8)")
print("  2. Finite subgroups of SU(2)")
print("  3. Simple singularities")
print("  4. Regular polyhedra")
print()
print("The McKay correspondence:")
print("  A_n ↔ Z_{n+1} (cyclic group)")
print("  D_n ↔ Dic_{n-2} (binary dihedral, order 4(n-2))")
print("  E_6 ↔ BT (binary tetrahedral, order 24)")
print("  E_7 ↔ BO (binary octahedral, order 48)")
print("  E_8 ↔ BI (binary icosahedral, order 120)")
print()
print("TOURNAMENT THEORY LIVES IN TYPE A:")
print("  A_{n-1} ↔ Z_n ↔ S_n (Weyl group)")
print("  Tournament on n vertices ↔ orientation of A_{n-1} root system")
print()
print("BUT: the Platonic solids enter via E₆, E₇, E₈!")
print("  E₆ (tetrahedron): det = 3 = KEY₂")
print("  E₇ (cube/octahedron): det = 2 = KEY₁")
print("  E₈ (dodecahedron/icosahedron): det = 1")
print()
print("THE KEYS ARE ENCODED IN THE EXCEPTIONAL-PLATONIC CORRESPONDENCE:")
print("  KEY₂ = 3 = det(E₆) ↔ tetrahedron (the SELF-DUAL solid)")
print("  KEY₁ = 2 = det(E₇) ↔ cube/octahedron (the DUAL PAIR)")
print("  1 = det(E₈) ↔ dodecahedron/icosahedron (the GOLDEN RATIO pair)")

# ====================================================================
print()
print("=" * 70)
print("PART 4: TOURNAMENTS ON PLATONIC SOLID VERTICES")
print("=" * 70)
print()

# Build tournaments on Platonic solid vertex sets
# Tetrahedron: K₄ (complete graph on 4 vertices)
# Cube: 8 vertices
# Octahedron: 6 vertices

# For the tetrahedron (4 vertices), there are 2^C(4,2) = 64 tournaments
# The rotation group A₄ (order 12) acts on these
# Number of non-isomorphic tournaments on 4 vertices = 4

print("TETRAHEDRON (4 vertices):")
print("  Total tournaments: 2^6 = 64")
print("  Non-isomorphic: 4")
print("  Symmetry group: A₄ (order 12) — alternating on 4 elements")
print("  # orbits = 64/12 = 5.33... (not uniform orbits)")
print()

# Count H values for all 64 tournaments on 4 vertices
n = 4
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
h_dist = {}
for bits in range(2**len(edges)):
    adj = [0]*n
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)
    H = 0
    for perm in permutations(range(n)):
        ok = True
        for idx in range(n-1):
            if not (adj[perm[idx]] & (1 << perm[idx+1])):
                ok = False
                break
        if ok:
            H += 1
    h_dist[H] = h_dist.get(H, 0) + 1

print(f"  H distribution: {dict(sorted(h_dist.items()))}")
print(f"  Possible H values: {sorted(h_dist.keys())}")
print(f"  Note: H ∈ {{1, 3, 5}} at n=4")
print()

# Octahedron (6 vertices)
print("OCTAHEDRON (6 vertices):")
print("  Total tournaments: 2^15 = 32768")
print("  Symmetry group: S₄ (order 24, acting on 3 pairs of antipodal vertices)")
print()

# Key question: can we orient octahedron edges to make a tournament?
# The octahedron has 12 edges. A tournament on 6 vertices has C(6,2)=15 edges.
# The octahedron is NOT the complete graph — it's missing 3 edges
# (the 3 diameters connecting antipodal vertices).
# So we can't directly "orient the octahedron" to get a tournament.
# But we CAN consider tournaments on the vertex set {0,...,5}.

# HOWEVER: the octahedron graph structure constrains which 3-cycles exist.
# A 3-cycle in the tournament that is also a triangle in the octahedron
# is a "geometric 3-cycle."

# The octahedron has 8 triangular faces.
# Each face is a triple of vertices. In a tournament, each triple
# is either a 3-cycle or transitive.

import random
random.seed(42)
n = 6
nsamples = 1000
geo_cycle_counts = {}

# Octahedron faces (vertices 0-5, with 0-5 as opposite pairs: {0,5}, {1,4}, {2,3})
oct_faces = [
    (0,1,2), (0,2,4), (0,4,1), (0,1,3), (0,3,4), (0,2,3),
    # Actually, the octahedron with antipodal pairs {0,5},{1,4},{2,3}:
    # has faces: triangles sharing no antipodal pair
]
# Standard octahedron: vertices at ±e_1, ±e_2, ±e_3
# Label: 0=+e₁, 1=+e₂, 2=+e₃, 3=-e₁, 4=-e₂, 5=-e₃
# Edges: all except antipodal (0-3, 1-4, 2-5)
# Faces: 8 triangles
oct_faces = [
    (0,1,2), (0,2,4), (0,4,1),  # top cap (around +e₁)
    (3,1,2), (3,2,4), (3,4,1),  # ... wait, need to think about this properly
]
# Vertices adjacent to 0 (+e₁): 1,2,4,5 (all except 3=-e₁)
# Faces containing 0: (0,1,2), (0,2,5), (0,5,4), (0,4,1) — 4 faces
# Similarly for each vertex.
# Total faces = 6·4/3 = 8 ✓

oct_faces = [
    (0,1,2), (0,2,4), (0,4,5), (0,5,1),  # 4 faces around vertex 0
    (3,1,5), (3,5,4), (3,4,2), (3,2,1),  # 4 faces around vertex 3 (opposite of 0)
]

print(f"  Octahedron faces: {len(oct_faces)} triangles")
for trial in range(5):
    adj = [0]*n
    for i in range(n):
        for j in range(i+1,n):
            if random.random() < 0.5:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)

    # Count 3-cycles among octahedron faces
    geo_3cyc = 0
    for face in oct_faces:
        a, b, c = face
        if (adj[a]&(1<<b)) and (adj[b]&(1<<c)) and (adj[c]&(1<<a)):
            geo_3cyc += 1
        elif (adj[a]&(1<<c)) and (adj[c]&(1<<b)) and (adj[b]&(1<<a)):
            geo_3cyc += 1

    # Count all 3-cycles
    all_3cyc = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (adj[i]&(1<<j)) and (adj[j]&(1<<k)) and (adj[k]&(1<<i)):
                    all_3cyc += 1
                elif (adj[i]&(1<<k)) and (adj[k]&(1<<j)) and (adj[j]&(1<<i)):
                    all_3cyc += 1

    H = 0
    for perm in permutations(range(n)):
        ok = True
        for idx in range(n-1):
            if not (adj[perm[idx]] & (1 << perm[idx+1])):
                ok = False
                break
        if ok:
            H += 1
    print(f"  T{trial}: geo_3cyc={geo_3cyc}/8, all_3cyc={all_3cyc}/C(6,3)=20, H={H}")

# ====================================================================
print()
print("=" * 70)
print("PART 5: THE COXETER NUMBER AND THE TOURNAMENT POLYNOMIAL")
print("=" * 70)
print()
print("Coxeter numbers of simple Lie algebras:")
print("  A_n: h = n+1")
print("  B_n: h = 2n")
print("  C_n: h = 2n")
print("  D_n: h = 2(n-1)")
print("  G₂: h = 6")
print("  F₄: h = 12")
print("  E₆: h = 12")
print("  E₇: h = 18")
print("  E₈: h = 30")
print()
print("COXETER NUMBERS OF THE EXCEPTIONALS:")
print(f"  G₂: h=6 = 2·3 = KEY₁·KEY₂ = product of keys")
print(f"  F₄: h=12 = 2²·3 = KEY₁²·KEY₂")
print(f"  E₆: h=12 = 2²·3 = KEY₁²·KEY₂ (same as F₄!)")
print(f"  E₇: h=18 = 2·3² = KEY₁·KEY₂² (9=KEY₂² appears!)")
print(f"  E₈: h=30 = 2·3·5")
print()
print("THE 5 APPEARS ONLY IN E₈!")
print("E₈ is the only exceptional group whose Coxeter number involves")
print("a prime other than 2 and 3. The prime 5 enters through the")
print("icosahedral symmetry (dodecahedron has C₅ symmetry).")
print()
print("TOURNAMENT POLYNOMIAL z² - 5z + 6 = 0:")
print("  5 = Coxeter number of A₄ (rank 4)")
print("  6 = Coxeter number of G₂ (rank 2)")
print("  So the tournament polynomial has:")
print("  coefficient 5 = h(A₄)")
print("  coefficient 6 = h(G₂)")
print("  roots 2 = det(A₁) = h(A₁)")
print("  roots 3 = det(A₂) = h(A₂)")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 6: THE 5-FOLD WAY")
print("=" * 70)
print()
print("WHY FIVE? The number 5 appears as:")
print()
print("  5 Platonic solids")
print("  5 exceptional Lie groups")
print("  5 = 2 + 3 (sum of keys)")
print("  5 = Φ₄(2) (fourth cyclotomic at 2)")
print("  5 = h(A₄) (Coxeter number of A₄)")
print("  5 = F(5) (Fibonacci fixed point)")
print()
print("THE DEEP REASON:")
print("  5 = 2 + 3 is the diameter of the key pair.")
print("  The 'classification gap' between keys 2 and 3.")
print()
print("  In ADE: A, D are INFINITE series. E₆, E₇, E₈ are FINITE exceptions.")
print("  The E series ENDS at rank 8 (E₈). Why?")
print("  Because E₉ would require a Coxeter number of 30+...")
print("  Actually: the E-series stops because the Cartan matrix must")
print("  be positive definite. E₈ is the last.")
print()
print("  In Platonic solids: they stop at 5 because higher-dimensional")
print("  analogs would require dihedral angles impossible on a 2-sphere.")
print("  The angle constraint is: π/(Schläfli number) must satisfy")
print("  a positivity condition — same as Cartan matrix positive-definiteness!")
print()
print("  THE SAME POSITIVITY CONSTRAINT stops both sequences at 5!")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 7: THE DUAL PAIRS AND THE KEYS")
print("=" * 70)
print()
print("Platonic dual pairs:")
print("  Tetrahedron ↔ Tetrahedron (self-dual)")
print("  Cube ↔ Octahedron")
print("  Dodecahedron ↔ Icosahedron")
print()
print("Tournament duality: T ↔ T̄ (reversal)")
print("  H(T) = H(T̄) (proved in S75)")
print("  The reversal is an INVOLUTION, like Platonic duality.")
print()
print("THE SELF-DUAL SOLID (tetrahedron) has:")
print("  4 vertices, 6 edges, 4 faces")
print("  Symmetry group |A₄| = 12 = 2²·3")
print("  V = F = 4 = 2² = KEY₁²")
print("  E = 6 = 2·3 = KEY₁·KEY₂")
print()
print("  In tournament theory, a 'self-dual' tournament is one with")
print("  T ≅ T̄ (has an anti-automorphism).")
print("  These are called SELF-COMPLEMENTARY tournaments.")
print("  They exist only for n ≡ 3 (mod 4): n = 3, 7, 11, 15, ...")
print()
print("  At n=3: the 3-cycle C₃ is self-complementary (and is a")
print("  tournament on the faces of a tetrahedron!)")
print("  The tetrahedron has 4 vertices but its FACE graph is K₄,")
print("  and orientations of K₃ (= face triples) give tournaments.")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 8: THE ICOSAHEDRAL CONNECTION — WHERE 5 MEETS 2 AND 3")
print("=" * 70)
print()
print("The icosahedron/dodecahedron pair is special:")
print("  It introduces the GOLDEN RATIO φ = (1+√5)/2")
print("  φ² = φ + 1 (the Fibonacci recurrence!)")
print("  φ = 1.618... (the 2-nacci root!)")
print()
print("  The k-nacci roots:")
print("  k=2: φ ≈ 1.618 (golden ratio, icosahedral)")
print("  k=3: ≈ 1.839 (tribonacci, n=9 regime)")
print("  k→∞: → 2 (KEY₁)")
print()
print("  The GOLDEN RATIO is the FIRST k-nacci root.")
print("  It lives in Q(√5), the simplest quadratic extension.")
print("  And √5 comes from the DISCRIMINANT of z²-5z+6:")
print(f"  Δ = 25 - 24 = 1 → √Δ = 1 (trivial!)")
print()
print("  BUT: φ satisfies z² - z - 1 = 0, with Δ = 5.")
print("  And 5 = the sum of keys = 2 + 3.")
print("  So √5 = √(KEY₁ + KEY₂).")
print()
print("  The golden ratio φ = (1 + √(KEY₁+KEY₂)) / KEY₁")
print()
print("  ICOSAHEDRAL SYMMETRY meets tournament theory:")
print("  The binary icosahedral group (BI) has order 120 = 5!")
print("  It corresponds to E₈ in the McKay correspondence.")
print("  E₈ is the LARGEST exceptional group.")
print("  Its Coxeter number h(E₈) = 30 = 2·3·5.")
print("  This is the PRODUCT of the three primes {2, 3, 5}.")
print()
print("  The HIERARCHY:")
print("  Level 1: KEY₁ = 2 (binary, A₁, cycle orientations)")
print("  Level 2: KEY₂ = 3 (ternary, A₂, cycle length)")
print("  Level 3: 5 = KEY₁+KEY₂ (quintic, golden ratio, icosahedral)")
print()
print("  2, 3, 5 are the first three primes.")
print("  They generate: 6=2·3, 10=2·5, 15=3·5, 30=2·3·5")
print("  These are the Coxeter numbers of the ADE system!")
print("  h(G₂)=6, h(A₉)=10, h(A₁₄)=15, h(E₈)=30")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 9: THE TOURNAMENT ON THE ICOSAHEDRON")
print("=" * 70)
print()
print("The icosahedron has 12 vertices, 30 edges, 20 triangular faces.")
print("A tournament on 12 vertices has C(12,2) = 66 edges.")
print("The icosahedron graph has 30 edges — less than half of 66.")
print()
print("But the FACE graph of the icosahedron has 20 faces,")
print("each a triple of vertices. In a tournament, each triple")
print("is either a 3-cycle (prob 3/4) or transitive (prob 1/4).")
print()
print("Expected # of 3-cycle faces: 20 · 3/4 = 15")
print("Expected α₁ at n=12: C(12,3)/4 = 55")
print("Face 3-cycles / total 3-cycles ≈ 15/55 ≈ 27%")
print()
print("The 20 icosahedral faces decompose into:")
print("  10 pairs of opposite faces (antipodal)")
print("  5 'rings' of 4 faces each (around each edge pair)")
print()
print("If we orient the tournament to MAXIMIZE 3-cycles among faces:")
print("  Want all 20 face triples to be 3-cycles")
print("  This constrains 20 · 3 = 60 edge orientations (with overlaps)")
print("  The icosahedron has 30 edges, each shared by 2 faces")
print("  So 40 constraints on 30 edges — OVERDETERMINED")
print("  Can we achieve all 20? Probably not.")
print()

# ====================================================================
print()
print("=" * 70)
print("PART 10: E₈ AND THE TOURNAMENT AT n=8")
print("=" * 70)
print()
print("E₈ has rank 8, dim 248.")
print("A tournament on 8 vertices lives in type A₇ (rank 7).")
print()
print("BUT: E₈ has deep connections to tournaments:")
print("  dim(E₈) = 248 = 8·31 = 8·(2⁵-1) = 8·Mersenne(5)")
print("  The number of positive roots: 120 = 5!")
print("  The Coxeter number: 30 = 2·3·5")
print("  The exponents: 1, 7, 11, 13, 17, 19, 23, 29 (primes < 30!)")
print()
print("  E₈ exponents are EXACTLY the integers coprime to 30 in [1,29]:")
print("  These are the totatives of 30 = 2·3·5.")
print("  φ(30) = 30·(1-1/2)·(1-1/3)·(1-1/5) = 8")
print("  So E₈ has rank = φ(h) = φ(30) = 8!")
print()
print("  THIS IS THE KEY RELATIONSHIP:")
print("  rank(E₈) = φ(h(E₈)) = φ(2·3·5) = 8")
print("  The rank of E₈ is the EULER TOTIENT of its Coxeter number!")
print()
print("  For type A_{n-1}: rank = n-1, h = n")
print("  φ(h) = φ(n) — not equal to rank in general")
print("  But: φ(p) = p-1 = rank for prime n=p!")
print("  So A_{p-1} has rank = φ(h) when h=p is prime.")
print()
print("  h(A₁) = 2 (prime): rank = 1 = φ(2) ✓")
print("  h(A₂) = 3 (prime): rank = 2 = φ(3) ✓")
print("  h(A₄) = 5 (prime): rank = 4 = φ(5) ✓")
print("  h(A₆) = 7 (prime): rank = 6 = φ(7) ✓")
print()
print("  The KEY primes 2, 3, 5, 7 all give rank = φ(h)!")
print("  These are the primes dividing 30 = h(E₈)... plus 7.")

# ====================================================================
print()
print("=" * 70)
print("PART 11: THE RECURRENCE VIEW — (2,3,5) HIERARCHY")
print("=" * 70)
print()
print("  LEVEL 1 (KEY₁ = 2):")
print("  z² - z - 1 = 0 → φ = (1+√5)/2 ≈ 1.618 (Fibonacci)")
print("  z - 2 = 0 → x = 2 (evaluation point for H)")
print("  Binary choice: each cycle has 2 orientations")
print()
print("  LEVEL 2 (KEY₂ = 3):")
print("  z² - 5z + 6 = 0 → roots 2, 3 (tournament polynomial)")
print("  Ternary structure: each cycle uses 3 vertices")
print("  n=9 = 3² is the CS boundary")
print()
print("  LEVEL 3 (5 = KEY₁ + KEY₂):")
print("  z² - 5z + 5 = 0 → roots (5±√5)/2 ≈ 3.618, 1.382")
print("  These are φ² = φ+1 and 1/φ² = ... no.")
print("  Actually: the roots are φ³ = 2+φ ≈ 3.618 and φ⁻³ ≈ ... no.")
print()
print("  Better: 5 appears as the DISCRIMINANT of the Fibonacci polynomial.")
print("  Δ(z²-z-1) = 1+4 = 5")
print("  And Δ(z²-5z+6) = 25-24 = 1 (trivial!)")
print()
print("  The tournament polynomial has DISCRIMINANT 1 (perfect square).")
print("  The Fibonacci polynomial has DISCRIMINANT 5 (= sum of keys).")
print()
print("  TOGETHER:")
print("  Fibonacci (k=2): gives φ → approaches KEY₁=2 as k→∞")
print("  Tournament poly: gives KEY₁=2, KEY₂=3 directly")
print("  The bridge: 5 = discriminant of Fibonacci = sum of keys")
print()
print("  THE THREE POLYNOMIALS:")
print("  z - 2 = 0            (KEY₁)")
print("  z² - z - 1 = 0       (Fibonacci, disc = 5 = KEY₁+KEY₂)")
print("  z² - 5z + 6 = 0      (tournament, disc = 1, roots = keys)")
print()
print("  They form a HIERARCHY:")
print("  1st poly degree 1: gives 2")
print("  2nd poly degree 2: gives φ (approaches 2)")
print("  3rd poly degree 2: gives 2 and 3")

# ====================================================================
print()
print("=" * 70)
print("PART 12: SYNTHESIS — THE FIVE-FOLD TOURNAMENT UNIVERSE")
print("=" * 70)
print()
print("  THE FIVE EXCEPTIONAL GROUPS ↔ FIVE TOURNAMENT STRUCTURES:")
print()
print("  G₂ (det=1, h=6=2·3):")
print("    The MINIMAL exceptional. 'Octonion structure.'")
print("    Tournament analog: the 6-vertex tournament (first α₂>0)")
print("    h = 6 = product of keys = first non-trivial independence level")
print()
print("  F₄ (det=1, h=12):")
print("    'Exceptional Jordan algebra.' 52-dimensional.")
print("    Tournament analog: n=12 (first α₄>0)")
print("    h = 12 = 4·3 = threshold for 4 disjoint 3-cycles")
print()
print("  E₆ (det=3, h=12):")
print("    Center = Z₃ = KEY₂!")
print("    McKay: binary tetrahedral (order 24)")
print("    Tournament analog: the KEY₂ structure")
print("    Tetrahedron = self-dual solid = self-complementary tournament")
print()
print("  E₇ (det=2, h=18):")
print("    Center = Z₂ = KEY₁!")
print("    McKay: binary octahedral (order 48)")
print("    Tournament analog: the KEY₁ structure (binary orientation)")
print("    h = 18 = 2·9 = 2·KEY₂² (the CS boundary!)")
print()
print("  E₈ (det=1, h=30):")
print("    Center = trivial. The LARGEST exceptional.")
print("    McKay: binary icosahedral (order 120 = 5!)")
print("    Tournament analog: the COMPLETE structure")
print("    h = 30 = 2·3·5 = product of first 3 primes")
print("    rank = 8 = φ(30) = 2³ = KEY₁³")
print()
print("  THE MASTER TABLE:")
print(f"  {'Group':<6} {'det':<5} {'h':<5} {'McKay order':<14} {'Key':<8} {'Tournament n'}")
print(f"  {'G₂':<6} {1:<5} {6:<5} {12:<14} {'2·3':<8} {'n=6'}")
print(f"  {'F₄':<6} {1:<5} {12:<5} {24:<14} {'2²·3':<8} {'n=12'}")
print(f"  {'E₆':<6} {3:<5} {12:<5} {24:<14} {'KEY₂':<8} {'n=12 (self-dual)'}")
print(f"  {'E₇':<6} {2:<5} {18:<5} {48:<14} {'KEY₁':<8} {'n=18=2·9=2·3²'}")
print(f"  {'E₈':<6} {1:<5} {30:<5} {120:<14} {'2·3·5':<8} {'n=30'}")
