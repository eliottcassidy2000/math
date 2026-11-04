# Independent proof review: realizable ears and complete boundary inversion

**Verdict: G1--G4 pass within their stated scopes.** This is an independent
mathematical and source-code review of
[`../geometry/PROOF.md`](../geometry/PROOF.md) and `../geometry/verify.py`.
The root worker performs separate direct colour counts and final replays;
those executions are not claimed as performed by this review.

## The generic placement argument is valid

For distinct u,v at distance strictly between0 and3 there is an open arc
on the unit circle about u for which the distance to v lies strictly
between0 and2. The circle intersections defining y are transverse there.
Thus local branches exist and can be restricted to a smaller arc.

For an old point w other than u, the simultaneous conditions
`|x-u|=|x-w|=1` have at most two solutions. Hence the only branch along
which y could be identically an old point is y=u, possible precisely when
`|u-v|=1`. The proof discards that branch. The other branch is nonconstant
and remains available on an open arc. Every other forbidden old-point
collision, and every unintended unit contact of x, removes finitely many
x positions. An unwanted y-w unit contact first restricts y to the at-most-
two intersections of circles centred at v and w; for each such y!=u,
the circle constraints on x again leave at most two positions. This
justifies the finite-exception argument, including the distance-one case.

The proof distinguishes intended contacts: x-u and y-v are allowed;
x-v and y-u are forbidden unless they are already one of the intended
edges, which distinctness prevents. The x-y unit contact itself prevents
their collision. All earlier private vertices can be treated as old
points during the next placement, so simultaneous probes are proved,
rather than assumed independent. There is no requirement that the old
induced unit graph be planar in the graph-theoretic sense; the assertion
concerns distinct points and Euclidean distances in the plane.

The boundary treatment is correct. At endpoint distance3, triangle
inequality equality fixes the two private locations on the straight
segment; an old point there can prevent privacy. At distance greater
than3 no three-unit-edge path is possible. The strict bound is therefore
the correct uniform guarantee with the old placement fixed.

## The full partition inverse uses actual available contexts

For a forest F_pi joining the smallest member of every block to all other
members, every edge subset S is again exactly the canonical star forest
of its component partition pi(S). Removing star edges isolates leaves
and preserves the smallest vertex as the centre of the remaining star.
This detail is essential: the inversion never calls for an unmeasured
forest context.

For each old colouring, a private ear contributes `K-delta_ij` with
K=q^2-3q+3. Multiplying over a forest gives the actual joint count because
its private vertices and edges have no additional interactions. Expanding
the product of deltas as the product of `K-(K-delta)` gives equation4
exactly. The resulting equality moment requires each block of pi to be
monochromatic, which is equivalent to its actual equality partition
coarsening pi. Subtraction over strict coarsenings then gives equation5.

Thus the Bell(m) query transformation is invertible on the entire formal
space of partition masses, not only on the masses realized by a graph.
Indeed it is unimodular for every integer K: ordering partitions by
refinement makes the forest expansion triangular with diagonal
`(-1)^(m-|pi|)`, and the coarsening zeta matrix is unit triangular. This
observation is already implied by the displayed integer inverse; it
does not add a new physical construction or a minimal-query theorem.

All original equality partitions remain present. For more than q blocks
the actual mass is zero; division by(q)_k is used only when k<=q.
Palette symmetry justifies the recovery of the weight of each fixed
labelled assignment within its equality pattern. The all-pairs distance
condition allows one common arrangement with an ear for every boundary
pair; induced subsets of its private vertices realize every forest query.

## The forcing obstruction is preserved

For q>=3, both endpoint extension counts are strictly positive. Therefore
every proper base colouring extends independently across every ear, and
the full boundary support is unchanged. Exact signed inversion of
measurements cannot be substituted for graph gluing: actual counts
multiply positive factors, while the inverse subtracts observations.
At q=2 the equal-endpoint extension count is zero, and the supplied
two-edge path plus three-edge ear is a valid odd-cycle hostile. The scope
excludes extra private edges and identifications, exactly where the
positivity product argument would need to be reconsidered.

## Arithmetic/source inspection

The implementation certifies the intended ear edges by exact algebraic
identities: the rational circle parameter has norm one, the displacement
and its perpendicular have dot product zero, and
`d2*(D+1/4)=1`. The interval containment check for a required edge is
therefore not being used as its sole proof of exact unit length.

The square-root bounds use integer square roots of rational endpoints
scaled by2^160, with outward rounding. Interval products handle all four
endpoint combinations, and squaring handles intervals crossing zero.
Strict positive squared-distance lower bounds exclude collisions; every
nonedge interval excludes1. The common sqrt3 and other radical dependencies
can widen these enclosures but cannot make them unsound. The coordinate
artifact correctly records the larger radical field, with radicands in
Q(sqrt3); it does not claim all private coordinates lie in Q(sqrt3).

The complete boundary table and the three-colouring certificate agree
with the mathematical scope: recovery of a zero orbit does not imply an
empty relation, and a three-colourable37-point example provides no new
six-chromatic witness. No defect was found in G2's existence proof or
G3's full forest tomography.
