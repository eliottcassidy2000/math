# Independent audit: boundary counts, seven modes, and the prescribed unit interface

2026-10-25. **PASS**, after the two prepublication repairs recorded below,
for [cube_boundary/RESULT.md](../cube_boundary/RESULT.md). The geometric
examples are prescribed six-palette interfaces, not ordinary six-color
obstructions. The independent script imports no producer implementation;
it reads the published coordinates and recomputes their consequences.

## Proof review

**Feasibility shadows.** A two-coordinate projection fiber on the binary
three-cube is exactly an edge. Thus a missing projected assignment means
an entire cube edge lies in the complement. This proves the independent-
set criterion. The independent-set counts `1,8,16,8,2` in sizes zero
through four are correct. In particular an independent set containing
opposite-parity points must contain an antipodal pair, which blocks every
third vertex. The claimed 35 full-shadow relations and their size table
follow without extrapolation from a computational sample.

**Numerical marginals.** Coordinate sums kill exactly those Walsh
characters using an omitted coordinate. All proper marginals therefore
leave a one-dimensional kernel spanned by the full-support character.
For integer tables, integrality of the coefficient follows from any
single entry difference. Nonnegativity gives the stated lower and upper
bounds with their signs as printed. For Boolean indicators any nonzero
difference must be exactly plus or minus that character, forcing the
complementary parity tetrahedra. This is stronger than feasibility
projection; the six-state NAE table and full cube do *not* have identical
numerical marginals. The document keeps that distinction throughout.

**Gluing.** Under the repaired hypothesis of literal graph union with
disjoint private interiors and no undecomposed cross-constraints, extension
choices multiply at each fixed boundary assignment. Fourier inversion
then gives the `1/8` convolution and total-count formula. Added shared-port
edges can be separate indicators; interior cross-constraints require more
joint boundary data. The minimal linear dimension assertion is correctly
scoped to the full eight-dimensional counting-table interface with
prescribed point-pin tests. It is not a lower bound for every specialized
fragment family or an unprescribed gadget-size theorem.

**Cube spectrum.** The stated matrix is `-2I+2A(Q3)`, so a weight-k Walsh
character has eigenvalue `-2+2(3-2k)=4(1-k)`. In particular the singleton
modes are a genuine kernel. The two nonnegative tables `1+chi_j` and
`1-chi_j` have equal matrix images and unequal pinned responses. This
proves exactly the claimed failure of that divisor-relation quotient
as a general coloring-table compression. No transfer of complex surface
geometry to colorability is inferred.

## Prescribed planar unit geometry

The six original points have exactly the triangle and three spokes under
the exact norm `a^2+ab+b^2` in the basis `1,alpha`. The rational circle
parametrization in the source has norm one identically. Each desired leaf
can be placed while avoiding finitely many bad circle points; its integer
parameter sequence has infinitely many distinct positions, so the greedy
selection terminates. The saved coordinates provide explicit witnesses.

The independent verifier recomputes all 496 pair distances of graph F
and all 171 pair distances of graph G using exact fractions. There are
exactly the declared 31 and 15 unit edges; all vertices are distinct and
every prescribed added vertex is a leaf of degree one. No hidden unit
edge between attached leaves, old ports, or the separate star survives.

For each of all eight prescribed binary port assignments, the independent
checker exhausts all colors on the remaining unpinned vertices (at most
`6^4` assignments), rather than importing the producer's recursive menu
search. It recovers

```
w_F = (0,4,4,4,4,4,4,0),
w_G = (3,3,3,3,3,3,3,3).
```

The fixed leaves forbid precisely the stated colors. There is no claim
that the underlying graph forces these leaf colors. Graph F is a
triangle with trees and a separate star; graph G is a forest. Thus the
examples cannot be mistaken for unprescribed six-uncolorable unit graphs.

Their shared mass is 24, all singleton and cubic Walsh coefficients
vanish, and both completion parities have mass 12. Their pair coefficients
differ by `-8` in each of the three pair channels. Their existential
pair shadows are all full. Pinning all ports to zero produces zero versus
three extensions, exactly as claimed.

## Repairs made before promotion

1. The draft prose said graph F had 28 edges. Its construction and checker
   gave `6+9+12+4=31`; the prose now says 31. This was a local arithmetic
   typo, not a changed graph or surviving mathematical claim.
2. The product formula originally said cross-edges were declared without
   explicitly excluding private-interior cross-constraints. The final
   text requires literal graph union with no undecomposed constraints,
   explains shared-port indicators, and requires boundary enlargement for
   interior cross-edges. The reviewed formula is correct under that scope.

No false promoted theorem or inherited-canon correction was needed.

## Independent finite checks and reproduction

[check_boundary_independent.py](check_boundary_independent.py) verifies all
256 relations, all 65,536 ordered gluing pairs, all 6,561 integer tables in
`{0,1,2}^8`, the complete exact coordinate graphs, and every port assignment.
It uses tuple coordinates rather than the producer's bit-mask ordering
for its Fourier and projection calculations. Normal and optimized
independent outputs agree byte for byte. The producer's normal and `-O`
runs were also replayed successfully with identical receipt hashes.

The [independent receipt](boundary_independent.json) pins the reviewed
proof, producer checker, producer receipt, and independent source hashes.
After a purely editorial status promotion, rerun the audit checker to
refresh its proof-byte pin.

```
python -B 05-knowledge/results/hadwiger_six_seven_eight_20261025/audits/check_boundary_independent.py 05-knowledge/results/hadwiger_six_seven_eight_20261025/audits/boundary_independent.json
python -B -O 05-knowledge/results/hadwiger_six_seven_eight_20261025/audits/check_boundary_independent.py 05-knowledge/results/hadwiger_six_seven_eight_20261025/audits/boundary_independent_optimized.json
```

No independent claim about the chromatic lower bound of the whole plane,
unprescribed exterior realizability, or a new K3 lattice is made here.
