# An infinite native family separating unindexed arcs from indexed path flags

Status: PROVED elementary construction for every prime q=3 mod 4, q>=7;
FINITE-EXACT checks at q=3,7,11,19,23,31 (the q=3 boundary is different).

## Construction

Let P_q have arc x->y exactly when y-x is a nonzero square modulo q. Let
D(T) have left and right copies of T, a new center c, left internal T,
right internal converse T, center->left, right->center, and left x->right y
exactly when x=y or x->y in T. Take A=D(D(P_q)), of order 4q+3.

Write L=q+1. Its four q-vertex copies have labels iL+x (i=0,1,2,3),
and its centers are q,2q+1,3q+2. Consider the ordered path

(3q+2), followed for x=0,1,...,q-1 by
(x, L+x, 3L+x, 2L+x), followed by (2q+1,q).

Every displayed step is an arc: the three within-coordinate transitions
are diagonal doubling edges or their internal converse, the next-coordinate
step 2L+x -> x+1 follows from 1 being a quadratic residue, and the initial
and final center steps follow directly from the doubling convention.
Thus this is an actual Hamiltonian path P with actual endpoints 3q+2 and q.

For each nonzero quadratic residue a and each b modulo q, define g_(a,b)
to fix all centers and send iL+x to iL+(ax+b mod q). Multiplication by a
preserves the Paley arcs, and doubling extends every source automorphism,
so all these q(q-1)/2 maps are source automorphisms. Their images of P are
pairwise distinct because equality of ordered paths forces equality on every
vertex. Denote this explicit family of paths by O.

## Exact incidence multiplicities

Put h=(q-1)/2. For an arc of the original P:

1. Each of the 3q within-coordinate arcs, and the final center-to-center arc,
   occurs in every one of the qh paths. The affine map merely permutes their
   coordinate parameter.
2. Each of the q-1 transitions 2L+x -> x+1 occurs in exactly q-1 paths.
   If the preimage transition has coordinate y, equality of the coordinate
   difference forces a=1, and b=x-y. The q-1 allowed y are all residues
   except q-1, giving q-1 distinct translations.
3. The initial center-to-copy arc occurs in h paths, those with b=0.
   The final copy-to-center arc occurs in h paths, those with
   a(q-1)+b=q-1.

These are all 4q+2 path arcs. For q>=7 the least multiplicity is h>=3.
Consequently NO single unindexed positive arc-occurrence observation
identifies P within O. This provides native counterexamples of unbounded
order, not merely an abstract inconsistent attachment model.

The initial and final copy-center arcs together impose b=0 and a=1.
Exactly two unindexed positive arc observations therefore suffice, and two
are necessary. This is relative to the declared affine orbit O; no claim
is made that these two arcs distinguish P among every Hamiltonian path of A.

## One indexed flag still suffices

At zero-based path positions 4 and 5, the ordered vertices are 2L and 1.
Requiring these actual owners at those actual positions forces b=0 and
then a=1. Thus one indexed path flag distinguishes P within O. An indexed
flag includes the path-position correspondence and both vertex identities;
it must not be replaced by a statement that the same arc occurs somewhere.

For q=3, h=1, so an initial arc already isolates P in its three-path family.
The exception q=3 and first obstruction q=7 follow from the exact stabilizer
sizes, not from a numerical analogy with forbidden Hamiltonian counts.

## Controls and scope

check_arc_family.py directly constructs the full tournaments and verifies
all candidate automorphisms, paths, incidence counts, the two-arc certificate,
and the indexed flag at q=3,7,11,19,23,31. Normal and optimized Python agree.
This proves no assertion about Hamiltonian-count gaps, complete metagraph
wiring, or arithmetic Collatz dynamics. Its contribution is an exact
all-size information distinction useful to an attachment-signature design.

Independent proof review by the parent confirmed the doubling edge directions,
all three incidence categories, and the two distinct observation certificates.
No full automorphism-group classification is needed for the counterexample.
