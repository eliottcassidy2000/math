# Independent audit of the plane cube obstruction and relation-space bridge

2026-10-25. Reviewer: the counting/recurrence agent, independent of the
geometry proof's author. **PASS:** Theorem U, the square-flat coloop
formula, and the typed K3 relation-space bridge are proved as stated in
[README.md](README.md). This review does not claim a new result about
six-colour forcing, the whole K3 action, or tournament connectivity.

## Proof audit

1. For a four-cycle with four distinct images in the Euclidean plane, the
   two opposite vertices are the distinct intersections of two unit
   circles. Their midpoint is the midpoint of the other opposite pair.
   Thus the parallelogram equation is valid even without an assumption
   about the drawing's crossings. Distinctness and ambient dimension two
   are both necessary hypotheses and are explicitly retained.
2. For distinct nonzero generators a,b of an elementary abelian 2-group,
   the four vertices of their square are distinct. The parallelogram
   equation makes Delta_a invariant under every other generator b.
   A dependence expressing a by the others then makes Delta_a invariant
   under a, while the involution makes it change sign. It vanishes,
   contradicting unit edge length. This proves the necessary linear
   independence with no assumption that the graph action acts by plane
   isometries.
3. For independent directions, generic sums of unit plane vectors give
   the required hypercubes. The inductive exceptional sets are finite:
   collisions exclude finitely many unit vectors, and each nonzero old
   difference d excludes at most two solutions of |u|=|u+d|=1. The old
   internal distances persist under translation. This also proves the
   stronger induced realization and permits rational unit vectors.
   Empty S and multiple connected components satisfy the stated result.
4. Walsh characters diagonalize every square equation. A character is
   allowed precisely when it is negative on at most one selected
   generator. A unique negative direction exists exactly for a coloop;
   the all-positive character is constant on each component. The dimension
   is therefore 2^(r-t)(1+b), including S empty. Every noncoloop direction
   is collapsed by every such map. There is no converse error from a
   dimension count alone: the argument identifies the actual functions.
5. For Q3 the three nonconstant permitted characters have support size
   one. The inherited matrix M=-2I+2A(Q3) has eigenvalue 4-4|J| on support
   J, so these are exactly ker M. Centered coordinate functions of every
   plane unit cube drawing lie in that kernel; the additional vector
   norms 1/2 and the stated noncollision condition are exactly what turn
   such functions into an injective unit drawing. The result identifies
   coefficient-vector spaces. It does not identify a Euclidean metric
   with the K3 intersection pairing.

## Inherited geometric inputs

The review checked the specific dependencies in
[the Fano-code result, section 6a](../../six_seven_geometry_20261025/fano_codes/RESULT.md)
and [the root-deletion result, section 4](../../six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md).
For the seven distinct nonzero labels, d1+d2+d3=c is nonzero. With all
d_j nonzero, absence of a bad triple is equivalent to their pairwise
distinctness, which together with their nonzero sum is exactly rank three.
The translation map then acts regularly on the eight component labels.
The two half-cube divisors give the three actual linear relations. The
other four nonconstant characters have negative intersection eigenvalues
and are not silently promoted to null divisor relations. Freeness on the
resolved hourglass is an inherited audited geometric input; freeness on
the entire K3 is neither used nor inferred here.

## Independent finite controls

An independent ephemeral standard-library Python audit ran successfully
with both `python` and `python -O`, using explicit exceptions throughout:

- All 128 subsets S of the seven nonzero vectors in F2^3 were checked by
  exact rational row reduction of the literal square-equation matrix.
  The matrix contains the rows f(x)+f(x+a+b)-f(x+a)-f(x+b), for every x
  and every distinct pair a,b in S. Its nullity agreed with
  2^(3-dim span(S))(1+number of coloops) in every case.
- Unit vectors (1,0), (3/5,4/5), (5/13,12/13) gave eight distinct centered
  signed half-sums. All 28 squared distances were checked exactly:
  distance one occurred iff the cube labels differed in one bit. Both
  coordinate columns were annihilated by -2I+2A(Q3).
- For S={1,2,4,7}, the noninjective map x -> (parity(x),0) gives all
  edges length one but violates a square equation. This is the explicit
  K4,4 hostile showing why injectivity cannot be dropped.
- The four vectors (1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1) have every
  squared mutual distance eight but violate the parallelogram equation.
  Rescaling yields the regular unit tetrahedron hostile in R3.

The maintained `check.py` supplements these review controls. The finite
checks support, but do not substitute for, the general proofs above.

No mathematical repair was required. The useful failure boundary is
precise: a regular elementary-2 action can preserve the abstract component
set while its dependent generator squares obstruct an injective plane
unit-distance realization.
