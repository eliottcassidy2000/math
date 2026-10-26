# Session71: dissection norms, boundary memory, and exact bridges

2026-10-26. First verified checkpoint: arithmetic and inherited-geometry lane.
The matching/Penrose synthesis is in progress. New elementary deductions are
separated from CITED external results and finite replay.

## A source correction before using the optimality theorem

The requested paper is Demaine–Kamata–Uehara,
[Dudeney's Dissection is Optimal, v3](https://arxiv.org/html/2412.03865v3).
The [official v4 notice](https://arxiv.org/abs/2412.03865v4) says v3 accidentally
compiled an outdated file and v4 corrects the cut-graph analysis. Our version
comparison is recorded separately; the served artifacts did not establish the
expected substantive delta. We do not silently replace versions or claim an
independent certification of the entire case analysis.

The external result is CITED: four pieces are optimal for square/equilateral
triangle dissections using polygonal pieces and no flipping. The paper uses
matching diagrams, topology, angles and lengths. Fermat's two-square theorem
is not a stated dependency. Flips and curved-piece extensions are not settled
by that theorem.

## The actual arithmetic bridge

[The norm and boundary report](session71-norms/REPORT.md) recovers the repo's
Gaussian-square/Farey/Pythagorean maps and the earlier rigid-ear obstruction
PS-D2 in `prime_shells_20260921_dissections.md`.

- Fermat's criterion distinguishes 5=1²+2² from 44, whose prime factor 11
  occurs to odd exponent. Gaussian squaring gives the exact Pythagorean map
  (m,n) ↦ (m²−n²,2mn,m²+n²). Norm-only records lose angular/determinant data.
- Normalize the equilateral triangle to side 2. Its equal-area square has
  side θ=3^(1/4), outside Q(√3). For regular-closed polygonal pieces, each square side is a sum of whole
  exposed polygon-edge lengths. Thus each side needs an exposed edge length outside
  Q(√3). This is an edge-length module obstruction, not a count of pieces:
  four exposed edge occurrences need not come from four distinct pieces.
- No square of that area has all real x/y coordinates in Q(√3), even tilted:
  its edge would require a²+b²=√3, and the conjugation √3↦−√3 contradicts
  a sum of real squares. Adjoining θ removes this target-coordinate
  obstruction. Coordinate fields and edge-length fields are distinct.

These facts explain an arithmetic component of geometric compatibility,
NOT the three-versus-four threshold. Indeed scaling the triangle/square pair
to area 5 or area 44 leaves its optimum unchanged, while square-to-itself at
either area uses one piece. Representability of an area cannot determine the
minimum dissection size without the shapes and permitted motions.

## Two sharper exact interfaces

[The independent norm/geometry audit](session71-norm-audit/REPORT.md) proves:

- Simultaneous target vertices require a real Cartesian coordinate field of
  degree at least four, attained by Q(θ). For complex vertex coordinates the
  minimum is eight, attained by Q(θ,i). No abelian complex number field works;
  every degree-eight Galois example has dihedral group D4. The obstruction is
  conjugate positivity of z·conjugate(z), not the number of pieces. The theorem
  concerns the target vertices, not the coordinates of Dudeney's actual cuts.
- For rational matrices and a Q-independent boundary basis β, the unconstrained
  length equations Ax=BCβ are solvable iff rank[A|BC]=rank A. This retains
  labelled relations, not just their number. Positivity, topology, orientations,
  nonoverlap, and a piece budget remain additional geometric requirements.

## What 44 and 5 count

Several recovered 44s and 5s concern genuinely different objects: strong-atom
rooted-complement targets, an order-six pure-black tournament count, !5=44
permutations, and the paper's displayed cut-type inventory. Equal cardinality
alone supplies no incidence-, symmetry-, or operation-preserving isomorphism.
The report pins the individual sources and their proved/cited/finite scopes.

## Verification and next boundary

Native sources, portable exact checks, normal/optimized receipts and reviewed
hashes accompany the arithmetic report. The full session will add the
matching-diagram transport carrier, precise Penrose placement tests and an
independent geometric audit. None of these implies a graph-to-geometry
isomorphism without a specified map and retained gluing data.
