#!/usr/bin/env python3
"""
polygon_simplex_dihedral_s90ab.py — Polygons, Simplices, Dihedral Groups on the Unit Circle
opus-2026-03-16-S90ab

The n vertices of a regular n-gon sit at ω^k = e^{2πik/n} on the unit circle.
The dihedral group D_n (order 2n) is the symmetry group.
The n-simplex is the complete graph K_{n+1} — which IS a tournament space.
An alkane C_nH_{2n+2} has its carbon backbone as a Hamiltonian path of the simplex.

EVERYTHING CONVERGES on the unit circle.
"""

import numpy as np
from math import sqrt, log, log2, pi, cos, sin, factorial

phi = (1 + sqrt(5)) / 2

# ======================================================================
# PART 1: ELEMENTS ON THE UNIT CIRCLE
# ======================================================================
print("=" * 70)
print("PART 1: THE PERIODIC TABLE ON THE UNIT CIRCLE")
print("=" * 70)

print("""
Each element Z has Cayley address x_Z = (Z-1)/(Z+1) in [0,1).
Map this to the unit circle: place element Z at angle θ_Z = 2π·x_Z.
The position on the unit circle: e^{iθ_Z}.
""")

elements = [(1,"H"), (2,"He"), (3,"Li"), (4,"Be"), (5,"B"), (6,"C"),
            (7,"N"), (8,"O"), (9,"F"), (10,"Ne"), (11,"Na"), (12,"Mg"),
            (13,"Al"), (14,"Si"), (17,"Cl"), (19,"K"), (26,"Fe"), (79,"Au")]

print(f"{'Z':>3} {'Sym':>3} | {'x_Z':>8} | {'angle/π':>8} | {'Re':>8} {'Im':>8} | note")
print("-" * 65)

for Z, sym in elements:
    xZ = (Z-1)/(Z+1)
    theta = 2*pi*xZ
    re = cos(theta)
    im = sin(theta)
    note = ""
    if Z == 1: note = "at 1 (real, positive)"
    elif Z == 3: note = "at -1 (ANTIPODAL to H!)"
    elif Z == 6: note = "near -i (PIVOT)"
    elif Z == 7: note = "at e^{3πi/2} = -i (FIXED POINT!)"
    elif Z == 11: note = "at e^{5πi/3}"
    print(f"{Z:3d} {sym:>3} | {xZ:8.4f} | {theta/pi:8.4f} | {re:8.4f} {im:8.4f} | {note}")

print(f"""
NITROGEN (Z=7) sits at angle 3π/2 = 270° = the position of -i!
And -i is a FIXED POINT of the Cayley transform Q!

  Q(-i) = (1+(-i))/(1-(-i)) = (1-i)/(1+i) = (1-i)²/2 = -2i/2 = -i ✓

The FORBIDDEN element sits at the FIXED POINT of the monad!

LITHIUM (Z=3) sits at angle π = 180° = at -1 on the unit circle.
  -1 is where Q has its OTHER pole: Q(-1) → 0.
  Li is ANTIPODAL to H on the unit circle.

The LiH bond: from 1 (H at angle 0) to -1 (Li at angle π).
  This is a DIAMETER of the unit circle!
  The LiH bond spans the FULL DIAMETER of the chemical circle.
""")

# ======================================================================
# PART 2: THE ALKANE BACKBONE AS A SIMPLEX PATH
# ======================================================================
print("=" * 70)
print("PART 2: ALKANES = HAMILTONIAN PATHS THROUGH THE SIMPLEX")
print("=" * 70)

print("""
The n-SIMPLEX Δ_{n-1} has n vertices, all pairwise connected = K_n.
A tournament on K_n is an ORIENTATION of the complete graph.
A Hamiltonian path through the tournament IS an ordering of the vertices.

ALKANE ISOMERS as tournament paths:
  The n-carbon backbone of an alkane can be arranged in different ways.
  Each arrangement = a TREE on n vertices.
  The number of labeled trees on n vertices = n^{n-2} (Cayley's formula!).
  The number of Hamiltonian paths = H(T) for the tournament T.

  For the TRANSITIVE tournament (total order): H = 1.
  This corresponds to the LINEAR alkane (n-pentane, etc.).

  For OTHER tournaments: H > 1.
  These correspond to BRANCHED isomers (isopentane, neopentane, etc.).
""")

# The simplex and its tournament paths
for n in range(2, 8):
    num_tournaments = 2**(n*(n-1)//2)
    cayley_trees = n**(n-2) if n >= 2 else 1
    print(f"  n={n}: C_{n}H_{2*n+2} has {cayley_trees} labeled trees, "
          f"{num_tournaments} tournaments, {factorial(n)} Hamiltonian paths in K_{n}")

print(f"""
THE SIMPLEX-TOURNAMENT BRIDGE:
  The (n-1)-simplex K_n has C(n,2) edges.
  Each edge can be oriented 2 ways → 2^{{C(n,2)}} tournaments.
  Each tournament has H(T) Hamiltonian paths.
  The AVERAGE H over all tournaments = n!/2^{{n-1}}.

  The alkane C_nH_{{2n+2}} uses ONE specific Hamiltonian path.
  The isomer count is related to the NUMBER OF DISTINCT PATHS
  modulo the symmetry group (automorphisms of the tree).
""")

# ======================================================================
# PART 3: DIHEDRAL GROUPS AND CYCLOALKANES
# ======================================================================
print("=" * 70)
print("PART 3: D_n = SYMMETRY OF CYCLO-C_nH_{2n}")
print("=" * 70)

print("""
The DIHEDRAL GROUP D_n (order 2n) is the symmetry group of the
regular n-gon inscribed in the unit circle.
  n rotations: r^k for k=0,...,n-1 (multiply by ω^k = e^{2πik/n}).
  n reflections: sr^k (complex conjugation composed with rotation).

CYCLOALKANE C_nH_{2n} has carbon backbone forming a regular n-gon.
The symmetry group of the cycloalkane IS D_n!
""")

print(f"{'n':>3} | {'D_n order':>9} | {'molecule':>14} | {'strain':>12} | polygon | our numbers")
print("-" * 80)

dihedral_data = [
    (3, "Cyclopropane", "115 kJ/mol", "triangle", "3 = dim(M)"),
    (4, "Cyclobutane", "110 kJ/mol", "square", "4 = Q order"),
    (5, "Cyclopentane", "26 kJ/mol", "pentagon", "5 = golden √5"),
    (6, "Cyclohexane", "0 kJ/mol", "HEXAGON", "6 = PIVOT = 2·3"),
    (7, "Cycloheptane", "26 kJ/mol", "heptagon", "7 = FORBIDDEN"),
    (8, "Cyclooctane", "40 kJ/mol", "octagon", "8 = |D₄| = Bott"),
    (10, "Cyclodecane", "50 kJ/mol", "decagon", "10 = C(5,2)"),
    (11, "Cycloundecane", "~50 kJ/mol", "hendecagon", "11 = DISC"),
    (12, "Cyclododecane", "16 kJ/mol", "dodecagon", "12 = icosa verts"),
]

for n, mol, strain, poly, note in dihedral_data:
    print(f"{n:3d} | {2*n:9d} | {mol:>14} | {strain:>12} | {poly:>9} | {note}")

print(f"""
THE DIHEDRAL HIERARCHY:
  D_3: cyclopropane = the TRIANGLE = the 3-cycle = β₁=1 (STRAINED)
  D_4: cyclobutane = the SQUARE = Q⁴=Id (but strained!)
  D_5: cyclopentane = the PENTAGON = golden ratio
  D_6: cyclohexane = the HEXAGON = UNSTRAINED = THE PIVOT
  D_7: cycloheptane = the HEPTAGON = first forbidden
  D_8: cyclooctane = the OCTAGON = Bott period

OUR D₄ (the Cayley group) corresponds to CYCLOBUTANE (n=4)!
  D₄ has 8 elements = 4 rotations + 4 reflections.
  The 4 rotations of D₄ ARE our Q⁰, Q¹, Q², Q³.
  The 4 reflections include negation (x → -x).

  But cyclobutane is STRAINED (110 kJ/mol).
  Our D₄ theory has "strain" too: the complex eigenvalues
  (the helix) represent the "bending" of the flat structure.
""")

# ======================================================================
# PART 4: THE REGULAR n-GON ON THE UNIT CIRCLE
# ======================================================================
print("=" * 70)
print("PART 4: REGULAR POLYGONS AND THEIR CAYLEY TRANSFORMS")
print("=" * 70)

print("""
Place the regular n-gon on the unit circle with vertices at ω^k = e^{2πik/n}.
Apply the Cayley transform Q to each vertex: Q(ω^k) = (1+ω^k)/(1-ω^k).
""")

for n in [3, 4, 5, 6, 7, 8, 11]:
    print(f"\nn-gon, n={n}: vertices at the {n}th roots of unity")
    Q_vals = []
    for k in range(n):
        wk = np.exp(2j*pi*k/n)
        if abs(1 - wk) < 1e-10:
            print(f"  k={k}: ω^{k} = 1.000, Q = ∞ (POLE)")
            Q_vals.append(None)
        else:
            Qk = (1 + wk) / (1 - wk)
            Q_vals.append(Qk)
            if n <= 7:
                print(f"  k={k}: ω^{k} = e^{{2πi·{k}/{n}}}, Q = {Qk.real:+.4f}{Qk.imag:+.4f}i, |Q|={abs(Qk):.4f}")

    # The Q-images (excluding pole) form a new polygon
    valid_Q = [q for q in Q_vals if q is not None]
    if len(valid_Q) >= 2:
        # Are the Q-images themselves on a circle?
        mods = [abs(q) for q in valid_Q]
        args = [np.angle(q) for q in valid_Q]
        print(f"  Q-image moduli: min={min(mods):.4f}, max={max(mods):.4f}, "
              f"ratio={max(mods)/min(mods):.4f}")
        # Check if all on the imaginary axis (arg = π/2 or -π/2)
        on_imag = all(abs(abs(a) - pi/2) < 0.01 for a in args)
        if on_imag:
            print(f"  ALL Q-images lie on the IMAGINARY AXIS!")

print(f"""
KEY OBSERVATION: Q maps the unit circle to the IMAGINARY AXIS!

This is because for |z| = 1:
  Q(z) = (1+z)/(1-z)
  |Q(z)|² = |1+z|²/|1-z|²
  For z = e^{{iθ}}: |1+z|² = 2+2cos(θ), |1-z|² = 2-2cos(θ).
  Q(z) is PURELY IMAGINARY when z is on the unit circle (except z=±1).

Actually: Q(e^{{iθ}}) = (1+e^{{iθ}})/(1-e^{{iθ}}) = -i·cot(θ/2).

So Q maps the unit circle to the IMAGINARY AXIS via -i·cot(θ/2)!

This is the CAYLEY TRANSFORM IN OPERATOR THEORY:
  The Cayley transform maps UNITARY operators (unit circle)
  to SELF-ADJOINT operators (real/imaginary axis).
  Our Q maps the circle of phases to the line of measurements.

For the n-gon: the vertices at ω^k map to Q(ω^k) = -i·cot(πk/n).
These are n-1 EQUALLY SPACED points on the imaginary axis (in cot-space)!
""")

# ======================================================================
# PART 5: LiH ON THE UNIT CIRCLE
# ======================================================================
print("=" * 70)
print("PART 5: THE LiH BOND ON THE UNIT CIRCLE")
print("=" * 70)

# H at angle 0 (position 1+0i), Li at angle π (position -1+0i)
print("""
H  (Z=1): angle 0, position = 1   on the unit circle.
Li (Z=3): angle π, position = -1  on the unit circle.

The LiH bond = the DIAMETER from 1 to -1.

Q(1) = ∞ (the pole — hydrogen is at the boundary of the Cayley domain)
Q(-1) = 0 (the zero — lithium maps to the origin)

The Cayley transform maps:
  H → ∞ (hydrogen is "infinitely fast" = at the speed of light)
  Li → 0 (lithium is "at rest" = zero rapidity)

This seems backwards! But it makes sense:
  H (Z=1) is the LIGHTEST element = highest speed/energy per mass.
  Li (Z=3) is the LIGHTEST metal = the "resting" reference frame.

The LiH bond is a CONNECTION between:
  ∞ (the boundary, the pole, the light cone) and 0 (the origin, rest).
  It spans the ENTIRE Cayley line!

In terms of the n-gon: LiH is a DIAMETER of the unit circle.
The dihedral group D_2 (= Z/2 × Z/2, order 4, the Klein four-group)
is the symmetry group of a diameter. And D_2 is a SUBGROUP of our D_4.
""")

# ======================================================================
# PART 6: THE METHANE TETRAHEDRON
# ======================================================================
print("=" * 70)
print("PART 6: METHANE CH₄ = THE 3-SIMPLEX ON THE UNIT SPHERE")
print("=" * 70)

print("""
METHANE: one carbon + four hydrogens in a TETRAHEDRON.

The tetrahedron is the 3-SIMPLEX = K_4 = the complete graph on 4 vertices.
It has 6 edges = C(4,2) = the number of arcs in a tournament on 4 vertices!

In our theory: a tournament on 4 vertices has 2^6 = 64 possible orientations.
METHANE corresponds to the TRANSITIVE tournament (all H atoms equivalent,
one canonical ordering by position).

THE TETRAHEDRON IN THE UNIT SPHERE:
  Place the 4 vertices of a regular tetrahedron on the unit sphere.
  Standard coordinates: (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1) (normalized).

  These 4 points project to 4 points on the complex unit circle:
  via stereographic projection from the north pole.

  The DIHEDRAL GROUP of the tetrahedron:
    Full symmetry: S_4 (order 24) = rotations + reflections.
    Rotation subgroup: A_4 (order 12) = even permutations.
    This is NOT a dihedral group! The tetrahedron breaks dihedral symmetry.

  But the face of the tetrahedron IS a triangle = D_3 symmetry.
  And cyclopropane C_3H_6 (the smallest cycloalkane) has D_3 symmetry.

METHANE ↔ TOURNAMENT ON 4:
  4 H atoms form a tournament under the "comparison" of their positions.
  The tournament is transitive (H=1) because all H atoms are equivalent
  in the tetrahedral symmetry → no 3-cycle possible.

  But if we BREAK the symmetry (e.g., CHFClBr), the tournament becomes
  non-trivial, and H > 1 (multiple orderings of the substituents).
  The NUMBER OF ORDERINGS = H(T) = the tournament invariant!
""")

# Compute: tournament on 4 vertices = isomers of substituted methane
print("Tournaments on 4 vertices (= C-X₁X₂X₃X₄ configurations):")
from itertools import permutations
from collections import Counter

def count_H_4(T):
    count = 0
    for perm in permutations(range(4)):
        ok = all(T[perm[p]][perm[p+1]] for p in range(3))
        if ok:
            count += 1
    return count

h_dist = Counter()
for bits in range(64):
    T = [[0]*4 for _ in range(4)]
    idx = 0
    for i in range(4):
        for j in range(i+1, 4):
            if (bits >> idx) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    h = count_H_4(T)
    h_dist[h] += 1

print(f"  H distribution: {dict(sorted(h_dist.items()))}")
print(f"  H=1 (transitive): {h_dist[1]} tournaments = {h_dist[1]}/64 = {h_dist[1]/64:.1%}")
print(f"  H=3 (one 3-cycle): {h_dist[3]} tournaments")
print(f"  H=5 (two 3-cycles): {h_dist[5]} tournaments = the NEAR-TRANSITIVE!")

# ======================================================================
# PART 7: THE UNIT CIRCLE UNIFIES EVERYTHING
# ======================================================================
print()
print("=" * 70)
print("PART 7: THE UNIT CIRCLE UNIFIES EVERYTHING")
print("=" * 70)

print(f"""
On the UNIT CIRCLE e^{{iθ}} for θ ∈ [0, 2π):

  ELEMENTS sit at θ_Z = 2π·(Z-1)/(Z+1):
    H at 0, Li at π, N at 3π/2 = -i (FIXED POINT!)

  REGULAR n-GONS have vertices at e^{{2πik/n}}:
    Triangle (n=3): cyclopropane, D₃, first cycle, β₁=1
    Square (n=4): cyclobutane, D₄ = our Cayley group
    Pentagon (n=5): cyclopentane, golden ratio, icosahedral
    Hexagon (n=6): cyclohexane, PIVOT, zero strain, benzene
    Heptagon (n=7): cycloheptane, FORBIDDEN, first hyperbolic

  THE CAYLEY TRANSFORM maps the circle to the imaginary axis:
    Q(e^{{iθ}}) = -i·cot(θ/2)
    Unitary (circle) → self-adjoint (line)
    Phases → measurements
    Quantum states → classical observables

  SIMPLICES: the n-simplex K_{{n+1}} is the tournament space.
    Each simplex face is a sub-tournament.
    A Hamiltonian path = a total order = an alkane backbone.
    H(T) = the number of distinct molecular orderings.

  THE LiH BOND = the DIAMETER (from 1 to -1):
    Spans the entire circle.
    Q maps it from ∞ to 0 (from boundary to center).
    The simplest bond is the longest chord.

  ALKANES C_nH_{{2n+2}}:
    Linear: transitive tournament (sim_H=1, one ordering)
    Branched: non-trivial tournament (H > 1, multiple orderings)
    Cyclic: ring tournament (β₁=1, Möbius twist, strain energy)

  THE DIHEDRAL HIERARCHY D_n:
    D₃: triangle = 3-cycle = tournament primitive
    D₄: square = Cayley monad = our symmetry group
    D₅: pentagon = golden ratio = icosahedral bridge
    D₆: hexagon = PIVOT = flat = zero strain = benzene
    D₇: heptagon = FORBIDDEN = hyperbolic threshold

  Ring strain energy is the CURVATURE of D_n:
    Strain(5) = 26 = Strain(7) = MIRROR IMAGES across D₆!

Everything meets on the unit circle:
  elements, molecules, symmetry groups, tournaments, and the Cayley transform.
  The circle is the GROUND STATE of the theory.
  Departures from the circle = curvature = strain = the helix.
""")

print("opus-2026-03-16-S90ab: POLYGONS, SIMPLICES, DIHEDRAL GROUPS ON THE UNIT CIRCLE")
