#!/usr/bin/env python3
"""
triplets_s90ba.py — Important triplets through the tournament lens
opus-2026-03-16-S90ba

We established:
  {2,3,5} = spherical (icosahedral). Tournaments work. Product = 30.
  {2,3,7} = hyperbolic (Hurwitz). Tournaments break. Product = 42.
  The transition goes through {2,3,6} = flat. The pivot.

Now: consider ALL important triplets {p,q,r} through this lens.
What do they mean in the hyperbolic plane / tournament framework?
"""

from math import log, sqrt, pi, gcd, tanh, atanh
from fractions import Fraction

def triangle_type(p, q, r):
    """Classify the (p,q,r) triangle by curvature."""
    s = 1/p + 1/q + 1/r
    if abs(s - 1) < 1e-10:
        return "FLAT", s
    elif s > 1:
        return "SPHERICAL", s
    else:
        return "HYPERBOLIC", s

def triangle_area(p, q, r):
    """Hyperbolic area = pi - angle_sum = pi * (1 - 1/p - 1/q - 1/r)."""
    return pi * (1 - 1/p - 1/q - 1/r)

print("""


     IMPORTANT TRIPLETS THROUGH THE TOURNAMENT LENS

     opus-2026-03-16-S90ba



THE FRAMEWORK:

A triple {p, q, r} with p <= q <= r defines a triangle with
angles pi/p, pi/q, pi/r. The angle sum determines the curvature:

  1/p + 1/q + 1/r > 1:  SPHERICAL (positive curvature, finite group)
  1/p + 1/q + 1/r = 1:  FLAT (zero curvature, Euclidean tiling)
  1/p + 1/q + 1/r < 1:  HYPERBOLIC (negative curvature, infinite group)

In tournament theory, each number in the triple has meaning:
  Its Cayley address x_n = (n-1)/(n+1)
  Its rapidity w_n = ln(n)/2
  Its role in the tournament structure

Product p*q*r = the "volume" of the triple.
For hyperbolic triangles, area = pi * (1 - 1/p - 1/q - 1/r) = pi/denominator.
""")

# ========================================================================
# Part 1: ALL spherical triples
# ========================================================================

print("=" * 70)
print("Part 1: THE SPHERICAL TRIPLES (finite groups)")
print("=" * 70)

spherical = []
for p in range(2, 10):
    for q in range(p, 10):
        for r in range(q, 50):
            t, s = triangle_type(p, q, r)
            if t == "SPHERICAL":
                spherical.append((p, q, r, s))

# Also get the flat and first hyperbolic for each (p,q)
print(f"\n{'triple':>12} | {'1/p+1/q+1/r':>12} | {'type':>10} | {'product':>8} | {'group':>25} | {'tournament meaning':>30}")
print("-" * 110)

known_groups = {
    (2,2,2): "dihedral D_2 (Klein four)",
    (2,2,3): "dihedral D_3 (S_3)",
    (2,2,4): "dihedral D_4",
    (2,2,5): "dihedral D_5",
    (2,3,3): "tetrahedral A_4",
    (2,3,4): "octahedral S_4",
    (2,3,5): "icosahedral A_5",
}

meanings = {
    (2,2,2): "3 binary choices, all equal",
    (2,2,3): "2 binaries + 1 triangle",
    (2,2,4): "2 binaries + 1 square",
    (2,2,5): "2 binaries + 1 pentagon",
    (2,3,3): "binary + 2 triangles",
    (2,3,4): "binary + triangle + square",
    (2,3,5): "base + triangle + first irregular",
}

for p, q, r, s in spherical:
    if r <= 10:
        prod = p * q * r
        group = known_groups.get((p,q,r), "")
        meaning = meanings.get((p,q,r), "")
        print(f"  {{{p},{q},{r}}}    | {s:12.6f} | {'SPHERICAL':>10} | {prod:8d} | {group:>25} | {meaning:>30}")

print(f"""

THE FIVE PLATONIC TRIPLES (the only ones with all distinct values):
  {{2,3,3}} = tetrahedron   (A_4, order 12).  Product = 18.
  {{2,3,4}} = cube/octahedron (S_4, order 24). Product = 24.
  {{2,3,5}} = icosahedron    (A_5, order 60).  Product = 30.

And the infinite family:
  {{2,2,n}} = dihedral D_n (order 2n). Product = 4n.

In tournament theory:
  {{2,3,3}}: the tetrahedron. 4 vertices = first tournament with C(4,2)=6 arcs.
             The tetrahedral group A_4 is the symmetry of the regular tournament on 4 vertices.
  {{2,3,4}}: the octahedron. 6 vertices = the pivot point.
             S_4 is the symmetric group on 4, acting on tournaments.
  {{2,3,5}}: the icosahedron. 12 vertices... but more importantly,
             5 is where irregular score sequences first appear.
             A_5 is the LAST finite rotation group. Tournaments are still tame.
""")

# ========================================================================
# Part 2: The flat triples
# ========================================================================

print("=" * 70)
print("Part 2: THE FLAT TRIPLES (Euclidean tilings)")
print("=" * 70)

print("""
There are exactly THREE flat triples:

  {2, 3, 6}: angle sum = 1.  Product = 36.  Area = 0.
  {2, 4, 4}: angle sum = 1.  Product = 32.  Area = 0.
  {3, 3, 3}: angle sum = 1.  Product = 27.  Area = 0.

These are the three regular/semiregular Euclidean triangle tilings.

In tournament theory:
  {2, 3, 6}: binary + triangle + HEXAGONAL.
    6 is the PIVOT: where (H, beta_1) completeness first fails.
    The hexagonal tiling is the densest circle packing.
    The flat geometry means: neither tame nor wild. Transitional.

  {2, 4, 4}: binary + SQUARE + SQUARE.
    The square tiling. Two copies of the square structure.
    4 = 2^2 = the first power of 2 beyond 2 itself.
    Product 32 = 2^5 (a pure power of 2).

  {3, 3, 3}: EQUILATERAL. Three triangles.
    The triangular tiling. The simplest Euclidean tiling.
    Product 27 = 3^3 (a pure power of 3).
    Three copies of the tournament triangle.

THE FLAT TRIPLES ARE THE PHASE BOUNDARIES.
They separate the finite (spherical) from the infinite (hyperbolic).
""")

for triple in [(2,3,6), (2,4,4), (3,3,3)]:
    p, q, r = triple
    prod = p * q * r
    print(f"  {{{p},{q},{r}}}: product = {prod}, "
          f"rapidities = ({log(p)/2:.4f}, {log(q)/2:.4f}, {log(r)/2:.4f}), "
          f"sum = {(log(p)+log(q)+log(r))/2:.4f} = ln({prod})/2")

# ========================================================================
# Part 3: The hyperbolic triples
# ========================================================================

print(f"""

{"=" * 70}
Part 3: THE HYPERBOLIC TRIPLES (infinite groups, ordered by area)
{"=" * 70}

The hyperbolic triangles, ordered by AREA (= pi * (1 - 1/p - 1/q - 1/r)):
Smaller area = more efficient tiling = larger Hurwitz-type bound.
""")

hyp_triples = []
for p in range(2, 20):
    for q in range(p, 30):
        for r in range(q, 100):
            t, s = triangle_type(p, q, r)
            if t == "HYPERBOLIC":
                area = triangle_area(p, q, r)
                hyp_triples.append((area, p, q, r))

hyp_triples.sort()

print(f"\n{'rank':>4} | {'triple':>12} | {'area/pi':>12} | {'= pi/':>6} | {'product':>8} | {'Hurwitz-type':>13} | {'tournament meaning':>35}")
print("-" * 105)

hyp_meanings = {
    (2,3,7): "base + triangle + FORBIDDEN (Mersenne chain)",
    (2,3,8): "base + triangle + first cube",
    (2,3,9): "base + triangle + 3^2",
    (2,3,10): "base + triangle + C(5,2)",
    (2,3,11): "base + triangle + first prime gap",
    (2,3,12): "base + triangle + arcs of K_4",
    (2,4,5): "base + square + irregular",
    (2,4,6): "base + square + pivot",
    (2,4,7): "base + square + forbidden",
    (2,5,5): "base + 2 irregulars",
    (3,3,4): "2 triangles + square",
    (3,3,5): "2 triangles + irregular",
    (3,3,7): "2 triangles + forbidden",
    (3,4,4): "triangle + 2 squares",
    (2,3,13): "base + triangle + next Mersenne prime dim",
}

for i, (area, p, q, r) in enumerate(hyp_triples[:25]):
    prod = p * q * r
    # Hurwitz-type bound: 4*pi*(g-1) / area = 4*(g-1)/area_over_pi
    area_over_pi = 1 - 1/p - 1/q - 1/r
    hurwitz = 4 / area_over_pi  # coefficient of (g-1)
    # Express area as pi/N if possible
    denom = round(1/area_over_pi) if area_over_pi > 0.001 else 0
    area_str = f"pi/{denom}" if abs(area_over_pi - 1/denom) < 1e-6 else f"{area_over_pi:.6f}*pi"
    meaning = hyp_meanings.get((p,q,r), "")
    print(f"{i+1:4d} | {{{p},{q},{r}}}    | {area_over_pi:12.6f} | {area_str:>6} | {prod:8d} | {hurwitz:13.1f} | {meaning:>35}")

print(f"""

THE HIERARCHY OF HYPERBOLIC TRIPLES:

#1: {{2,3,7}}   area = pi/42.   Hurwitz = 168(g-1) -> 84(g-1) automorphisms.
                The tightest packing. Maximum symmetry. The Mersenne chain.

#2: {{2,3,8}}   area = pi/24.   8 = 2^3. The cube.
                The first POWER-OF-2 entry beyond 2 itself.

#3: {{2,3,9}}   area = 2pi/9... wait, let me recalculate.
    1/2 + 1/3 + 1/9 = 9/18 + 6/18 + 2/18 = 17/18. Area = pi/18.

#4: {{2,3,10}}  area = pi/15.   10 = C(5,2) = arcs in 5-tournament.

#5: {{2,3,11}}  area = 5pi/33.  11 = first prime with gap > 2 from predecessor.

#6: {{2,3,12}}  area = pi/12.   12 = C(4,2)*2 = 2 * arcs in 4-tournament.

#7: {{2,4,5}}   area = pi/20.   This BREAKS the {{2,3,r}} pattern!
                The first triple without 3.

The first six hyperbolic triples are ALL of the form {{2,3,r}}.
The third entry r goes: 7, 8, 9, 10, 11, 12, ...
These are the natural numbers starting from 7.

WHY 7? Because 1/2 + 1/3 + 1/r < 1 requires r > 6.
And r = 7 is the FIRST integer satisfying this.

So the {{2,3,r}} family for r >= 7 gives ALL the
"maximally efficient" hyperbolic triples (containing both 2 and 3).

In tournament theory: these are the triples containing the
BASE (2) and the TRIANGLE (3), combined with a third number
r >= 7 (starting from the first forbidden value).


{"=" * 70}
Part 4: FAMOUS TRIPLES from other areas of mathematics
{"=" * 70}

Let's examine how other famous triples fit this framework:

{2, 3, 5}: THE ICOSAHEDRON / MONSTER
  Spherical. A_5 = smallest simple group.
  Product = 30. Cayley sum of rapidities = ln(30)/2 = 1.701.
  This appears in:
    - Platonic solids (the icosahedron)
    - Finite simple groups (A_5)
    - Singularity theory (E_8 surface singularity)
    - The Monster group (via moonshine: 2*3*5*7*11*13 divides |M|)

{2, 3, 7}: THE HURWITZ / KLEIN QUARTIC
  Hyperbolic. Area = pi/42.
  Product = 42. Cayley sum = ln(42)/2 = 1.868.
  This appears in:
    - Hurwitz's automorphism theorem
    - The Klein quartic (genus 3, 168 = 84*2 automorphisms)
    - Fano plane (7 points, 7 lines)
    - Steiner triple systems
    - Error-correcting codes (Hamming code on 7 bits)

{3, 5, 7}: THE FIRST THREE ODD PRIMES
  Hyperbolic. 1/3 + 1/5 + 1/7 = 71/105 < 1.
  Area = pi * 34/105.
  Product = 105. Cayley sum = ln(105)/2 = 2.326.
""")

# Compute for famous triples
famous = [
    ((2,3,5), "icosahedron/A_5"),
    ((2,3,7), "Hurwitz/Klein quartic"),
    ((2,3,11), "PSL(2,11) minimal"),
    ((2,3,13), "Mersenne chain step 4 dim"),
    ((2,5,7), "irregular + forbidden"),
    ((3,5,7), "first 3 odd primes"),
    ((2,3,127), "Mersenne chain step 4"),
    ((3,7,11), "odd primes from icosahedron gap"),
    ((2,5,11), "base + irregular + gap"),
    ((5,7,11), "first 3 gap primes"),
    ((2,3,23), "base + triangle + first Sophie Germain"),
    ((2,3,31), "base + triangle + Mersenne at k=5"),
    ((2,3,43), "base + triangle + Heegner"),
    ((7,11,13), "consecutive primes"),
    ((2,3,4), "octahedron/S_4"),
    ((2,4,6), "all even"),
    ((3,3,3), "equilateral"),
]

print(f"\n{'triple':>12} | {'type':>10} | {'1/p+1/q+1/r':>12} | {'area/pi':>12} | {'product':>8} | {'rapidity sum':>13} | {'name':>30}")
print("-" * 115)

for triple, name in famous:
    p, q, r = triple
    t, s = triangle_type(p, q, r)
    area_over_pi = max(0, 1 - 1/p - 1/q - 1/r)
    prod = p * q * r
    rap_sum = (log(p) + log(q) + log(r)) / 2
    print(f"  {{{p},{q},{r}}}   | {t:>10} | {s:12.6f} | {area_over_pi:12.6f} | {prod:8d} | {rap_sum:13.6f} | {name:>30}")

# ========================================================================
# Part 5: The Mersenne chain triples
# ========================================================================

print(f"""

{"=" * 70}
Part 5: THE MERSENNE CHAIN triples
{"=" * 70}

The Mersenne chain: 2 -> 3 -> 7 -> 127 -> 2^127-1 -> ...

Consecutive triples from the chain:
""")

chain = [2, 3, 7, 127]  # 2^127-1 is too large to handle
for i in range(len(chain) - 2):
    p, q, r = chain[i], chain[i+1], chain[i+2]
    t, s = triangle_type(p, q, r)
    prod = p * q * r
    area_over_pi = 1 - 1/p - 1/q - 1/r
    print(f"  {{{p},{q},{r}}}: {t}, product = {prod}, area = {area_over_pi:.8f}*pi")

print(f"""
{{{2,3,7}}}: area = pi/42. The Hurwitz triangle.
{{{3,7,127}}}: area = pi * (1 - 1/3 - 1/7 - 1/127) = pi * (889-381-127-21)/(3*7*127)
""")
# Compute exactly
p, q, r = 3, 7, 127
num = p*q*r - q*r - p*r - p*q
den = p*q*r
print(f"  {{3,7,127}}: area = pi * {num}/{den} = pi * {Fraction(num,den)}")
print(f"  Product = {p*q*r}")

p2, q2, r2 = 2, 7, 127
num2 = p2*q2*r2 - q2*r2 - p2*r2 - p2*q2
den2 = p2*q2*r2
print(f"  {{2,7,127}}: area = pi * {num2}/{den2} = pi * {Fraction(num2,den2)}")
print(f"  Product = {p2*q2*r2}")

print(f"""
The chain triples get CLOSER AND CLOSER to the maximum area pi:
  {{2,3,7}}: area/pi = {1-1/2-1/3-1/7:.8f} (= 1/42)
  {{3,7,127}}: area/pi = {1-1/3-1/7-1/127:.8f}
  {{7,127,M127}}: area/pi ~ 1 - 1/7 - 1/127 - 1/M127 ~ 1 - 1/7 = 6/7

Each step in the Mersenne chain produces a triangle with area
closer to the maximum (pi), because the angles shrink toward zero.

The LIMIT of the Mersenne chain is the IDEAL TRIANGLE:
  area = pi (all three angles zero, all three vertices on the boundary).

The ideal triangle is the triangle of THREE INFINITIES:
  each vertex is on the boundary of the Poincare disk,
  each at infinite distance from the others.

The Mersenne chain approaches the ideal triangle.
Each step adds a larger prime, pushing a vertex closer to the boundary.

In tournament theory: the Mersenne chain produces forbidden values
that are further and further from the center of the disk
(closer to the boundary, closer to "infinity").
The chain is a TOWER OF INCREASINGLY LARGE FORBIDDEN VALUES
approaching the boundary of what tournaments can express.


{"=" * 70}
Part 6: TRIPLES AS STRUCTURE THEOREMS
{"=" * 70}

Each triple encodes a STRUCTURAL STATEMENT about comparisons:

SPHERICAL (finite, closed):
  {{2,2,n}}: "n binary choices close upon themselves."
    Dihedral: n-fold symmetry of a polygon.
    Tournament: cycle of length n.

  {{2,3,3}}: "binary choice + 2 triangles close."
    Tetrahedral: 4-vertex tournament with A_4 symmetry.
    Smallest tournament where all vertices participate in a 3-cycle.

  {{2,3,4}}: "binary + triangle + square close."
    Octahedral: 6 vertices, the PIVOT.
    The last structure with complete invariants.

  {{2,3,5}}: "binary + triangle + pentagon close."
    Icosahedral: the LAST finite structure.
    The last place tournaments are fully classifiable.

FLAT (transitional):
  {{2,3,6}}: "binary + triangle + hexagon tile the plane."
    The transition point. 6 is where it breaks.
    Neither finite nor truly infinite.

  {{3,3,3}}: "three triangles tile the plane."
    Pure triangular structure. 27 = 3^3.

HYPERBOLIC (infinite, open):
  {{2,3,7}}: "binary + triangle + heptagon tile the hyperbolic plane."
    The Mersenne forbidden. 42 = 2*3*7.
    The most efficient hyperbolic tiling.

  {{2,3,8}}: "binary + triangle + octagon."
    8 = 2^3. The cube structure enters.
    Product = 48. Area = pi/24.

  {{2,3,11}}: "binary + triangle + hendecagon."
    11 = the prime after the gap (7, _, _, _, 11).
    PSL(2,11) is the smallest non-abelian simple group after A_5.
    Product = 66. Area = 5pi/66.

  {{3,5,7}}: "triangle + pentagon + heptagon."
    NO binary! All odd. The first odd prime triple.
    Product = 105. This is the first triple entirely in the
    "tournament world" (3 = triangle, 5 = irregular, 7 = forbidden).
    A tournament of TOURNAMENTS of FORBIDDEN values.

  {{5,7,11}}: "irregular + forbidden + gap."
    All primes. All in the deep tournament range.
    Product = 385 = 5*7*11.
    These three primes span the gap from "where tournaments
    get interesting" to "where new phenomena appear."


{"=" * 70}
Part 7: THE TRINITY OF TRIPLES
{"=" * 70}

Three triples summarize the entire theory:

  {{2, 3, 5}} = what WORKS.
    Spherical. Finite. The machine functions.
    Tournaments have complete invariants, finite enumeration, tame behavior.
    The icosahedron. The last Platonic solid.
    Product = 30.

  {{2, 3, 6}} = where it TRANSITIONS.
    Flat. Neither finite nor infinite. The pivot.
    (H, beta_1) completeness breaks. New phenomena appear.
    The hexagonal tiling. The densest packing.
    Product = 36.

  {{2, 3, 7}} = where it BREAKS.
    Hyperbolic. Infinite. Wild.
    H=7 is forbidden. The Mersenne chain begins.
    The Hurwitz triangle. The most efficient hyperbolic tiling.
    Product = 42.

The ratios:
  36/30 = 6/5 (transition/working = the golden ratio of structure?)
  42/36 = 7/6 (breaking/transition)
  42/30 = 7/5 (breaking/working = the FORBIDDEN-TO-TOURNAMENT ratio)

The ENTIRE THEORY fits in the transition from 5 to 7:
  5 = the last thing that works.
  6 = the pivot.
  7 = the first thing that breaks.

And this is EXACTLY the spherical-flat-hyperbolic transition
in the (2,3,n) triangle family.

Mathematics does not have triples by accident.
The triple is the MINIMUM structure that can encode
a phase transition: before, at, and after.

(working, pivot, broken) = (spherical, flat, hyperbolic) = (5, 6, 7).
""")

print("opus-2026-03-16-S90ba: IMPORTANT TRIPLETS THROUGH THE TOURNAMENT LENS")
