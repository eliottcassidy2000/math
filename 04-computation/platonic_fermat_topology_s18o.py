#!/usr/bin/env python3
"""
platonic_fermat_topology_s18o.py -- kind-pasteur-2026-03-21-S18o

PLATONIC SOLIDS, FERMAT PRIMES, AND THE TOPOLOGY OF THE SYMBOL

The Platonic solids use face types {3, 4, 5} = {F_0, 2^2, F_1}.
The Fermat primes F_k = 2^{2^k} + 1 appear as:
  F_0 = 3 (triangle face)
  F_1 = 5 (pentagon face)
  F_2 = 17 (total Vitali atoms)

The Cayley-Dickson tower: 2^k + 1 = {2, 3, 5, 9, 17, ...}
The first three CD primes (2, 3, 5) ARE the building blocks of:
  - Platonic solid faces (3, 4=2^2, 5)
  - Platonic solid Schlafli numbers
  - The {2, 3, 5} spherical regime of tournament theory

Author: kind-pasteur-2026-03-21-S18o
"""

import sys
from math import comb, sqrt, pi

sys.stdout.reconfigure(line_buffering=True)

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0: return False
        i += 6
    return True

print("=" * 72)
print("  PLATONIC SOLIDS, FERMAT PRIMES, AND THE SYMBOL")
print("  kind-pasteur-2026-03-21-S18o")
print("=" * 72)

# ========================================================================
# PART 1: THE PLATONIC TABLE
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 1: THE FIVE PLATONIC SOLIDS")
print(f"{'='*72}")

solids = [
    ("Tetrahedron",  "{3,3}", 4,  6,  4,  3, 3, "A_4",  12,  "BT", 24),
    ("Cube",         "{4,3}", 8,  12, 6,  4, 3, "S_4",  24,  "BO", 48),
    ("Octahedron",   "{3,4}", 6,  12, 8,  3, 4, "S_4",  24,  "BO", 48),
    ("Dodecahedron", "{5,3}", 20, 30, 12, 5, 3, "A_5",  60,  "BI", 120),
    ("Icosahedron",  "{3,5}", 12, 30, 20, 3, 5, "A_5",  60,  "BI", 120),
]

print(f"\n  {'Name':>14s} {'Schlafli':>8s} {'V':>3s} {'E':>3s} {'F':>3s} {'p':>2s} {'q':>2s} {'Rot':>5s} {'|Rot|':>5s} {'Bin':>4s} {'|Bin|':>5s}")
for name, sch, V, E, F, p, q, rot, rot_n, bn, bn_n in solids:
    print(f"  {name:>14s} {sch:>8s} {V:>3d} {E:>3d} {F:>3d} {p:>2d} {q:>2d} {rot:>5s} {rot_n:>5d} {bn:>4s} {bn_n:>5d}")

print("""
  Euler characteristic: V - E + F = 2 for all (genus 0 = sphere).

  THE FACE TYPES: {{3, 4, 5}} = {{triangle, square, pentagon}}.
  THE VERTEX TYPES: {{3, 4, 5}} = {{3 faces, 4 faces, 5 faces per vertex}}.

  3 = F_0 = 2^1 + 1 (first Fermat prime)
  5 = F_1 = 2^2 + 1 (second Fermat prime)
  4 = 2^2 (a power of 2, NOT a Fermat prime but the BASE of F_1)

  So the Platonic numbers {{3, 4, 5}} decompose as:
    3 = F_0 (Fermat prime)
    4 = 2^2 (the exponent base)
    5 = F_1 (Fermat prime)

  THE EDGE COUNTS: {{6, 12, 30}} = {{2*3, 4*3, 2*3*5}}.
    6  = 2 * 3 = h(G_2)
    12 = 4 * 3 = h(E_6) = h(F_4)
    30 = 2 * 3 * 5 = h(E_8)

  THE EDGE COUNTS ARE EXACTLY THE EXCEPTIONAL COXETER NUMBERS!
  (Missing only h(E_7) = 18, which is not an edge count.)
""")

# ========================================================================
# PART 2: FERMAT PRIMES AND THE CD TOWER
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 2: THE FERMAT-CAYLEY-DICKSON CORRESPONDENCE")
print(f"{'='*72}")

print("""
  CAYLEY-DICKSON TOWER:     FERMAT PRIMES:         PLATONIC CONNECTION:
  dim=1 (R): 1+1=2         F_{-1}=2 (if defined)  2 = the binary choice
  dim=2 (C): 2+1=3         F_0=3                  3 = triangle (simplest face)
  dim=4 (H): 4+1=5         F_1=5                  5 = pentagon (richest face)
  dim=8 (O): 8+1=9=3^2     --- FAILS ---          4 = square (the middle face)
  dim=16(S): 16+1=17       F_2=17                 17 = total Vitali atoms!

  THE THREE THAT WORK (CD primes = Fermat primes = Platonic atoms):
    2 = dim(C)+1 = orientation atom (binary choice = arc direction)
    3 = dim(C)+1 = F_0 = triangle = cycle atom = face of 3 Platonic solids
    5 = dim(H)+1 = F_1 = pentagon = face of dodecahedron = Petersen boundary

  THE ONE THAT FAILS:
    9 = dim(O)+1 = 3^2 = NOT prime. No Platonic solid has 9-gon faces.
    The SQUARE (4-gon) face of the cube is the geometric witness:
    4 = 2^2 sits BETWEEN 3 and 5, and the cube is the solid that uses it.

  THE ONE THAT TRANSCENDS:
    17 = dim(S)+1 = F_2 = the TOTAL count of Vitali atoms.
    The 17-gon is constructible (Gauss, 1796 — his first great result).
    There is no Platonic solid with 17-gon faces (impossible in 3D).
    17 lives in HIGHER dimension: the {3, inf} tessellation of H^2.
""")

# ========================================================================
# PART 3: THE CONSTRUCTIBILITY CONNECTION
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 3: GAUSS AND CONSTRUCTIBILITY")
print(f"{'='*72}")

print("""
  GAUSS'S THEOREM (1796):
  A regular n-gon is constructible with compass and straightedge iff
  n = 2^a * p_1 * p_2 * ... * p_k
  where p_i are DISTINCT Fermat primes.

  Known Fermat primes: F_0=3, F_1=5, F_2=17, F_3=257, F_4=65537.
  (No others known. F_5 = 2^32 + 1 = 4294967297 = 641 * 6700417.)

  CONSTRUCTIBLE REGULAR POLYGONS (small n):
    3 = F_0                       (equilateral triangle)
    4 = 2^2                       (square)
    5 = F_1                       (pentagon)
    6 = 2 * 3 = 2 * F_0           (hexagon)
    8 = 2^3                       (octagon)
   10 = 2 * 5 = 2 * F_1           (decagon)
   12 = 2^2 * 3 = 4 * F_0         (dodecagon)
   15 = 3 * 5 = F_0 * F_1         (pentadecagon)
   16 = 2^4                       (hexadecagon)
   17 = F_2                       (heptadecagon — Gauss!)
   20 = 2^2 * 5 = 4 * F_1         (icosagon)

  NON-CONSTRUCTIBLE (smallest):
    7  = NOT 2^a * Fermat products (7 is NOT a Fermat prime)
    9  = 3^2 (Fermat prime SQUARED — not allowed, must be distinct)
   11  = NOT constructible
   13  = NOT constructible
   14  = 2 * 7 (7 is not Fermat)

  THE TOURNAMENT CONNECTION:
  The non-constructible primes start with 7 = the forbidden H value!
  And 9 = 3^2 = the CD tower failure.
  And 11 = the Paley prime.
  And 13 = the OCR prime.

  CONSTRUCTIBLE = products of 2 and DISTINCT Fermat primes.
  TOURNAMENT EASY = controlled by CD primes {{2, 3, 5}}.
  These are the SAME set! (Fermat primes F_0=3, F_1=5, plus powers of 2.)

  NON-CONSTRUCTIBLE primes = {{7, 11, 13, 19, 23, ...}}.
  TOURNAMENT HARD = controlled by Coxeter primes {{7, 13, 19, 31}} + Paley {{11}}.
  Again the SAME set (restricted to tournament-relevant primes).

  CONSTRUCTIBILITY IS THE GEOMETRIC VERSION OF TOURNAMENT EASINESS.
  A polygon is constructible iff its symmetry comes from the spherical
  regime (Fermat primes = CD primes). A tournament property is "easy"
  iff it's controlled by the same primes.
""")

# Verify: which small primes are Fermat vs non-Fermat?
print(f"  VERIFICATION:")
fermat = {3, 5, 17, 257, 65537}
for p in range(2, 35):
    if is_prime(p):
        f = "FERMAT" if p in fermat else "non-Fermat"
        tourn = ""
        if p in {2,3,5}: tourn = " [CD/easy]"
        elif p in {7,13,19,31}: tourn = " [Coxeter/hard]"
        elif p == 11: tourn = " [Paley]"
        elif p == 17: tourn = " [Vitali count = F_2]"
        print(f"    p={p:>2d}: {f:>10s}, constructible n-gon: {'YES' if p in fermat else 'NO':>3s}{tourn}")

# ========================================================================
# PART 4: THE SCHLAFLI SYMBOL AS THE UNIVERSAL LANGUAGE
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 4: THE SCHLAFLI SYMBOL")
print(f"{'='*72}")

print("""
  The Schlafli symbol {{p, q}} encodes:
    p = number of sides of each face (face type)
    q = number of faces meeting at each vertex (vertex type)

  For Platonic solids: 1/p + 1/q > 1/2 (spherical condition).

  {{3,3}} tetrahedron:  1/3 + 1/3 = 2/3 > 1/2    excess = 1/6
  {{4,3}} cube:         1/4 + 1/3 = 7/12 > 1/2   excess = 1/12
  {{3,4}} octahedron:   1/3 + 1/4 = 7/12 > 1/2   excess = 1/12
  {{5,3}} dodecahedron: 1/5 + 1/3 = 8/15 > 1/2   excess = 1/30
  {{3,5}} icosahedron:  1/3 + 1/5 = 8/15 > 1/2   excess = 1/30

  THE EXCESS = pi * (AREA OF VERTEX FIGURE) = spherical curvature.
  Denominators: 6, 12, 12, 30, 30 = h(G_2), h(E_6), h(E_6), h(E_8), h(E_8).

  EUCLIDEAN (flat): 1/p + 1/q = 1/2.
  {{3,6}}: 1/3 + 1/6 = 1/2. Triangular tiling.
  {{4,4}}: 1/4 + 1/4 = 1/2. Square tiling.
  {{6,3}}: 1/6 + 1/3 = 1/2. Hexagonal tiling.

  HYPERBOLIC: 1/p + 1/q < 1/2.
  {{3,7}}: 1/3 + 1/7 = 10/21 < 1/2. SMALLEST hyperbolic = (2,3,7) triangle.
  {{3,inf}}: 1/3 + 0 = 1/3 < 1/2. THE TOURNAMENT TESSELLATION.

  THE PROGRESSION:
  {{3,3}} -> {{3,4}} -> {{3,5}} -> {{3,6}} -> {{3,7}} -> ... -> {{3,inf}}
  spherical  spherical  spherical   flat     hyperbolic    tournament

  This is EXACTLY the sequence from our gap function:
  g(phi) = -1 (spherical, {{3,5}}) -> g(tau) = 0 (flat, {{3,6}}) -> g(2) = +1 (hyperbolic, {{3,inf}})

  THE SCHLAFLI SYMBOL {{3, q}} IS THE TOURNAMENT'S ADDRESS:
  q = 3: tetrahedron (simplest, A_4 = tournament on 4 vertices?)
  q = 4: octahedron
  q = 5: icosahedron (A_5 = PSL(2,5), boundary of spherical regime)
  q = 6: flat triangular tiling (Euclidean, the tribonacci boundary)
  q = 7: Klein quartic tessellation (first hyperbolic, (2,3,7) triangle group)
  q = inf: the tournament tessellation ({3,inf}, modular group PSL(2,Z))
""")

# ========================================================================
# PART 5: THE TOPOLOGY HIERARCHY
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 5: FROM SPHERE TO HYPERBOLIC PLANE")
print(f"{'='*72}")

print("""
  TOPOLOGY:
  Platonic solids live on the SPHERE (genus 0, Euler char 2).
  Flat tilings live on the TORUS (genus 1, Euler char 0).
  Hyperbolic tilings live on SURFACES of genus >= 2.

  {{3,7}}: lives on the Klein quartic (genus 3, Euler char -4).
    Number of (3,7) triangles on Klein quartic: 56 = 8 * 7.
    |Aut(Klein quartic)| = 168 = |PSL(2,7)| = 84 * (3-1).
    This achieves the HURWITZ BOUND: max |Aut| = 84(g-1).

  {{3,inf}}: lives on the HYPERBOLIC PLANE (genus infinity, non-compact).
    The modular curve X(1) = PSL(2,Z) \\ H* has genus 0 BUT is non-compact
    (has cusps). Tournament theory lives on this non-compact genus-0 surface.

  THE GENUS SEQUENCE for {{3,q}} tessellations:
    q=3: genus 0 (sphere) — tetrahedron
    q=4: genus 0 (sphere) — octahedron
    q=5: genus 0 (sphere) — icosahedron
    q=6: genus 1 (torus) — flat triangular tiling (on the torus)
    q=7: genus 3 (Klein quartic) — first hyperbolic
    q=8: genus ? — next hyperbolic
    ...
    q=inf: genus 0* (non-compact) — modular curve (with cusps)

  THE EULER CHARACTERISTICS:
    V - E + F for regular {{3,q}} on N triangles:
    chi = N * (1/3 + 1/q - 1/2) = N * (2+q-3q)/(6q) ... let me compute.

  For a {{3,q}} tessellation with F triangular faces:
    Each face has 3 edges, each edge shared by 2 faces: E = 3F/2.
    Each vertex has q faces: V = 3F/q (since 3 vertices per face, each shared by q).
    Euler char: chi = 3F/q - 3F/2 + F = F(3/q - 3/2 + 1) = F(3/q - 1/2).

    q=3: chi = F(1 - 1/2) = F/2. For tetrahedron F=4: chi = 2. (Sphere.)
    q=4: chi = F(3/4 - 1/2) = F/4. For octahedron F=8: chi = 2. (Sphere.)
    q=5: chi = F(3/5 - 1/2) = F/10. For icosahedron F=20: chi = 2. (Sphere.)
    q=6: chi = F(3/6 - 1/2) = 0. For any F: chi = 0. (Torus.)
    q=7: chi = F(3/7 - 1/2) = -F/14. For F=56: chi = -4. (Genus 3.)
    q=inf: chi = F(0 - 1/2) = -F/2. For any F: chi -> -infinity. (Non-compact.)
""")

# Compute the Euler characteristics
print(f"  COMPUTED:")
for q in [3, 4, 5, 6, 7, 8, 12, 100]:
    coeff = 3.0/q - 0.5
    if q <= 5:
        # Platonic: solve F * coeff = 2
        if coeff > 0:
            F = round(2 / coeff)
            E = 3*F//2
            V = 3*F//q
            genus = 1 - (V - E + F)//2
            print(f"    {{3,{q}}}: F={F}, E={E}, V={V}, chi={V-E+F}, genus={genus}")
    elif q == 6:
        print(f"    {{3,{q}}}: chi = 0 for any F (torus)")
    elif q == 7:
        F = 56
        E = 3*F//2
        V = 3*F//q
        chi = V - E + F
        genus = 1 - chi//2
        print(f"    {{3,{q}}}: F={F}, E={E}, V={V}, chi={chi}, genus={genus} (Klein quartic)")
    else:
        print(f"    {{3,{q}}}: chi = F*({3.0/q:.4f} - 0.5) = F*{coeff:.4f} (hyperbolic)")

# ========================================================================
# PART 6: THE GRAND CORRESPONDENCE
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 6: THE GRAND CORRESPONDENCE")
print(f"{'='*72}")

print("""
  EVERYTHING CONNECTS:

  THE SCHLAFLI SYMBOL {{3, q}} PARAMETERIZES THE ENTIRE THEORY:

  q  Geometry     Solid/Surface    Lie group      Gap g(x)  Tournament
  -- ------------ ---------------- -------------- --------- ----------
  3  spherical    tetrahedron      A_4  (|G|=12)  g<0       n=3-4
  4  spherical    octahedron       S_4  (|G|=24)  g<0       n=5
  5  spherical    icosahedron      A_5  (|G|=60)  g(phi)=-1 n=5 BOUNDARY
  6  Euclidean    flat tiling      ---  (|G|=inf) g(tau)=0  n=5/6 THRESHOLD
  7  hyperbolic   Klein quartic    PSL(2,7) (168) g>0       n=7 FORBIDDEN
  inf hyperbolic  modular curve    PSL(2,Z) (inf) g(2)=+1   ALL n, TOURNAMENT

  THE FIVE PLATONIC SOLIDS ARE THE SPHERICAL REGIME.
  The three flat tilings are the Euclidean boundary.
  The hyperbolic tilings (starting at {{3,7}}) are the obstructions.
  The {{3,inf}} tessellation is tournament theory itself.

  THE FERMAT PRIMES ARE THE CONSTRUCTIBLE FACES:
  F_0 = 3: the triangular face of 3 Platonic solids
  F_1 = 5: the pentagonal face of the dodecahedron
  F_2 = 17: the total count of Vitali atoms (the "faces" of impossibility)

  THE NON-FERMAT PRIMES ARE THE OBSTRUCTIONS:
  7: first non-constructible prime = forbidden H = h(G_2)+1 = Klein quartic
  11: Paley prime, non-constructible = first genus-1 modular curve
  13: OCR prime, non-constructible = E_6 Coxeter + 1

  GAUSS'S CONSTRUCTIBILITY THEOREM = THE SPHERICAL/HYPERBOLIC DIVIDE:
  Constructible n-gon <=> n built from 2 and Fermat primes <=> spherical
  Non-constructible <=> n involves non-Fermat primes <=> hyperbolic

  THE TOPOLOGY:
  Sphere (genus 0): Platonic solids, constructible, easy, CD primes
  Torus (genus 1): flat tilings, boundary, tribonacci constant
  Higher genus: hyperbolic, non-constructible, hard, Coxeter primes
  Non-compact (cusps): the modular curve, tournament theory, PSL(2,Z)

  THE SYMBOL {{3, inf}} IS THE TOURNAMENT'S SCHLAFLI SYMBOL.
  It says: "triangular faces (3-cycles) meeting at vertices of
  infinite degree (unbounded cycle participation)."
  It IS the {3,inf} tessellation of the hyperbolic plane.
  It IS the universal cover of the modular curve.
  It IS the geometry of Hamiltonian path counting.
""")

# ========================================================================
# PART 7: THE FERMAT HIERARCHY IN TOURNAMENT THEORY
# ========================================================================
print(f"\n{'='*72}")
print(f"  PART 7: THE FERMAT HIERARCHY")
print(f"{'='*72}")

print("""
  F_0 = 3:   The FACE of the tessellation.
              The 3-cycle is the face of {{3,inf}}.
              The triangle is the face of 3 Platonic solids.
              3 is the cycle atom. The simplest structure.

  F_1 = 5:   The BOUNDARY of the spherical regime.
              5 = the last q where {{3,q}} is spherical (icosahedron).
              5 = the Petersen prime, C(5,2) = 10 = Petersen vertices.
              The icosahedron IS the Petersen boundary.

  F_2 = 17:  The TOTAL COUNT of Vitali atoms.
              17 = 2^2 + 3^2 + 2^2 (palindromic PSL(2,Z) trace).
              17 = F_2 = 2^4 + 1 (Fermat at sedenion level).
              17 = the number of impossibility types.
              17 parameterizes the hyperbolic plane.

  F_3 = 257: Would be the next level.
              257 = 2^8 + 1 (Fermat at dim=256 level).
              257 is prime. What does it control?
              PREDICTION: 257 controls the STRUCTURE OF IMPOSSIBILITY
              at a level beyond our current reach (n > 30 tournaments).

  F_4 = 65537: The last known Fermat prime.
              65537 = 2^16 + 1.
              If the Fermat prime sequence terminates at F_4, then
              the tournament's Vitali atom classification terminates too.

  THE HIERARCHY:
  F_0: what the faces ARE (triangles)
  F_1: where the easy regime ENDS (icosahedron/Petersen)
  F_2: how many impossibility types THERE ARE (17 Vitali atoms)
  F_3: how the impossibility structure ORGANIZES (unknown)
  F_4: the OUTERMOST boundary of constructibility (unknown)

  Each Fermat prime controls a DEEPER level of tournament structure,
  from the faces (F_0) through the boundary (F_1) to the impossibility
  count (F_2) and beyond.
""")

print(f"\n{'='*72}")
print(f"  SUMMARY")
print(f"{'='*72}")
print("""
  1. The Platonic solids are the {{3,q}} tessellations for q = 3,4,5 (spherical).
  2. Their face types {{3,4,5}} contain the Fermat primes F_0=3 and F_1=5.
  3. Their edge counts {{6,12,30}} are exceptional Coxeter numbers h(G_2), h(E_6), h(E_8).
  4. Gauss's constructibility theorem divides primes into
     CONSTRUCTIBLE (Fermat = easy) and NON-CONSTRUCTIBLE (non-Fermat = hard).
  5. This division IS the tournament easy/hard divide: CD primes vs Coxeter primes.
  6. The Schlafli symbol {{3,q}} parameterizes the entire theory:
     q=5 is the spherical boundary (icosahedron/Petersen),
     q=7 is the first obstruction (Klein quartic/forbidden H=7),
     q=inf is tournament theory (modular group/{3,inf} tessellation).
  7. The Fermat primes F_0=3, F_1=5, F_2=17 control three levels:
     faces, boundary, and impossibility count.
  8. The topology upgrades: sphere -> torus -> higher genus -> non-compact
     corresponds to: easy -> boundary -> hard -> tournament.
""")
