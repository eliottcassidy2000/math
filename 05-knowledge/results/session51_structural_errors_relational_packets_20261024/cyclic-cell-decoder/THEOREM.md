# Prime cyclic cells: edge rigidity and cycle frustration

Status: PROVED conditional on the explicitly verified cell domain; computational
certificate in check.py. No novelty priority claim. Not an all-tournament
classification, metagraph wiring theorem, or native routing theorem.

## Inheritance

Session39's affine triangle connections and Session50's signed-cell error
theorem already give the triangle kernel and error 6 per frozen triangle.
Session40 PRIME_CELL_BURNSIDE_FORMULAS.md classifies two matching-coupled
prime cells; regular-blocks/REGULAR_BLOCK_COUNTS.md counts arbitrary regular
crossings between two cyclic prime cells. We retain those facts as inherited.
The extension here is their complete multi-cell compatibility system, including
nonmatching rigid crossings, mixed prime sizes, and an exact root-error decoder.

Let stable cells mean the individual cells of directed equitable refinement
started from one color. They are invariant individually, not an arbitrary
block system. Assume each cell is a singleton or an induced R_p, for an odd
prime p, where vertices of R_p are F_p and x points to x+d for
1 <= d <= (p-1)/2. The condition is algorithmically verified, not inferred
from regularity or prime cardinality. In particular arbitrary regular seven
cells are outside the domain.

## 1. Internal and crossing actions

Write m=(p-1)/2. In R_p the arc x to x+d has m-d common outneighbors.
Consequently the arcs with m-1 common outneighbors recover the directed
successor cycle. Every automorphism preserves it, and all rotations preserve
R_p. Therefore Aut(R_p)=C_p. This argument is inherited from Session40.

Cells of coprime sizes have uniform crossings: counting the crossing in two
ways gives a*r=b*s, and coprimality forces r=0 or b. Thus singleton crossings
and crossings between distinct prime sizes impose no rotation constraint.

For two p-cells let B(x,y) be their binary crossing matrix, with constant
row and column degree d. Define its full translation stabilizer

    L(B) = {(a,b) in F_p^2: B(x+a,y+b)=B(x,y) for all x,y}.

It is a vector subspace because it is a subgroup of F_p^2. If d=0 or p it
is the full space. Otherwise it contains neither a nonzero horizontal nor a
nonzero vertical vector: invariance under all row shifts, for example, would
make every column constant and its degree either 0 or p, contrary to
0<d<p. Hence precisely two possibilities remain:

1. L(B)={0}, a rigid crossing;
2. L(B)={(a,lambda*a):a in F_p}, for one lambda in F_p^*, a gain crossing.

This classification uses the actual full crossing matrix. Row degrees alone
do not determine the case or the gain. Complementing B leaves L(B) unchanged.
Reversing the ordered cell pair exchanges coordinates and inverts lambda.

## 2. Complete decoder

For each p form a graph whose vertices are the p-cells. Put an oriented
gain edge i to j, labelled lambda, for every gain crossing, with relation
k_j=lambda*k_i. A rigid crossing pins BOTH endpoints to zero. Uniform
crossings impose nothing and give no edges. Pins then propagate throughout
each connected component of the gain graph.

A gain component is balanced if every oriented cycle has gain product 1;
inverse gains are used when an edge is traversed backwards. A component is
free exactly when balanced and containing no pinned vertex. Let b_p be the
number of free components and u_p the number of p-cells not in free components.

THEOREM. The complete ordinary automorphism permutation group is

    Aut(T) isomorphic to product_p C_p^(b_p).

Each free p-cell is a single vertex orbit, each nonfree p-cell consists of
p singleton orbits, and each singleton cell is fixed. In particular the
refinement error defined in Session50 is exactly

    e(T) = sum_p p*(p-1)*u_p.

PROOF. Every automorphism fixes each stable cell, so its restrictions must
be independent internal rotations k_i. Crossing preservation is exactly
membership of (k_i,k_j) in each L(B), with no additional arc constraint left.
Distinct primes decouple. Along a connected gain component all values are
determined by a root value. A pin forces that value to zero. An inconsistent
cycle gives (g-1)k=0 with g!=1, forcing zero over F_p. Otherwise every root
value is permitted. A nonzero root rotation acts transitively in each cell
of a free component. Fixed p-cells contribute p*(p-1) to the sum of
cell size minus actual orbit size over vertex orbits. This proves every claim.

The checker constructs a generating rotation for each free component,
including its actual vertex permutation, so the conclusion retains the
permutation action rather than merely its abstract group order.

## 3. What a cycle-only model misses

At p=3 every nonuniform biregular crossing is a permutation matrix or its
complement, and every permutation is affine. Rigid crossings do not occur;
only inconsistent loop gains can freeze a component. At p=5, rigid crossings
already occur between just two cells, a graph with no cycle. A nonaffine
matching is an explicit witness. Thus cycle balance alone is insufficient
outside the triangle domain. The repair is a gain graph WITH rigid-edge pins.

This does not say an unbalanced component is impossible. It says its internal
rotational symmetry is trivial. Symmetry loss, state existence, charge coverage,
and executable no-reuse routing are different predicates.

## 4. Scope and representation loss

The source is the entire tournament plus its verified cyclic stable-cell
coordinates; the target is gain graphs with rigid-edge pins. The map preserves
the complete automorphism action, vertex orbits, and refinement error. It loses
crossing translations and many nontranslation matrix features, and therefore
does not reconstruct the tournament isomorphism class or complete metagraph.
Full crossing matrices must remain available for those questions. It is
invariant under changes of origins and cell relabellings; converse is handled
by transforming the complete tournament, not by dropping its coordinates.

The algorithm tests each crossing's p^2 translations against p^2 entries,
then solves gain components by propagation. It is polynomial in the number
of input vertices and avoids a full automorphism search in this domain.
Recognition of R_p is by its successor cycle and complete arc verification.
No claim is made that constructing all parent tournaments is polynomial,
or that this domain exhausts the parent atlas at arbitrary order.

## 5. Relation to indivisible routing changes

Rigid crossing pins and directed implication constraints both show why
acyclic interaction graphs can still impose severe restrictions. The formal
analogy is loss of freedom on a tree. No predicate-preserving map from these
automorphism kernels to native no-reuse routes is established here. In
particular a fixed vertex in this theorem need not be a frozen owner record
in a routing problem. That distinction is a required sidecar, not a caveat
to omit when applying the decoder.
