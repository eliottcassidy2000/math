# Even graphs, tree trijections, and odd symmetry descent

2026-10-29. A research synthesis with explicit proofs, independently tested
counterexamples, and reusable code. The source snapshot is
`3eb49fa63d7ff4e58742b1ae84944f832e87a118`. Classical results and inherited
repository mechanisms are distinguished from this session's extensions.
No external-priority claim or major-open-problem solution is made.

## The main change of viewpoint

The useful common object is **a representation together with the operations and
symmetries it must preserve**. Four different successes must not be conflated:

1. An invertible object encoding preserves the actual structure. Under the
   same acting group, an equivariant bijection also preserves stabilizers;
   enlarged graph encodings must separately exclude extra automorphisms.
2. A bijection of isomorphism classes can change stabilizers and sampling weights.
3. A Fourier change of basis preserves invariant functions and counts, while
   using signed combinations rather than matching individual objects.
4. A sufficient state for one counting observer can forget information that
   later attachment operations still need.

These are familiar mathematical distinctions, not a newly claimed universal
theory. The contributions below make them meet on the same tournament examples,
with explicit maps, exact losses, and concrete repairs.

## 1. Even graphs extend to paired tournaments, with a joint sign test

[Royle et al.'s theorem](https://arxiv.org/abs/2204.01947) equates tournament
classes with **automorphism-even** graph classes. Here the sign counts reversed
edge endpoints. It is neither vertex-degree parity nor the parity of a
permutation of the edge set. K2 separates these definitions immediately.

The [decorated Fourier theorem](even/DECORATED_EVEN_GRAPH_DUALITY.md) extends
the inherited Walsh mechanism to any finite decoration set. For consecutive
pairs, write their two aligned cores as T and T XOR R, with root mask d. The
partner objects are triples (X,R,d) whose **joint stabilizer** has trivial
Royle sign on X. This gives counts

    1, 2, 8, 88, 2752

through n=4, hence the corresponding odd self-converse class counts through
order9 via Session89's lift. Requiring X alone to be even incorrectly gives
1,2,6,60,1092. Attaching distinguishable ports really changes which signs survive.

This is a linear correspondence, not a natural object bijection. At n=5 an
explicit Sylow3/5 Hall witness has520 source classes and only516 possible
partner classes for an equivariant class-bijective map. The obstruction extends
to every n>=5 for decorations with a global fixed point.

The [integral refinement](even/FOURIER_LATTICE.md) is stronger than count equality:
the change-of-basis cokernel is killed by2^binom(n,2). Every odd prime is
unexceptional for this operator, even if it divides the acting group order.
At n=3 the exact88-by-88 matrix has determinant magnitude2^116 and cokernel
exponent8. At n=4 the stabilizer formula gives2^7940; the large matrix was not
independently eliminated. Thus odd-prime cancellation phenomena elsewhere in
the repository need a different operator, not just the same symmetry vocabulary.

## 2. A genuine tree-based trijection, and an enlarged sparse graph codec

The [tree-boundary construction](trijection/TREE_COMPRESSION.md)
uses a simple exact map: a nonedge {u,v} of a spanning tree corresponds to the
first and last edges of its unique u-v path. The inverse takes the outer
endpoints of those two edges. Ordered versions preserve reversal and the
actual tree automorphism action.

Consequently a tournament on n vertices is equivalently an orientation of a
chosen spanning tree plus a tournament on that tree's n-1 edges. This extends
the earlier path/star examples to every tree. It does not make H or tournament
transitivity equal to the corresponding edge-tournament statistic, nor make
forgetting the marked tree free.

A [second construction](trijection/SPARSE_EVEN_CARRIERS.md) encodes any tournament faithfully into a larger ordinary
simple graph with exactly one fewer edge than vertices, preserving its full
automorphism group. The resulting graph is also Royle-even. A parallel generic
graph codec preserves the actual Royle character, even for source groups of
even order. This joins both requested graph themes by an explicit map, while
respecting the size change.

The condition e=v-1 alone is not tree-like: it says total cycle rank equals
number of components minus one. An arbitrary hard cyclic component can be
balanced by tree components. The generic construction gives a polynomial,
automorphism-exact reduction of graph isomorphism to this disconnected sparse
class. This diagnoses the weakness of a global edge budget; it does not solve GI.

[Capacity and padding bounds](capacity/CAPACITY_AND_EULER_REPAIR.md) quantify
the price. Already at n=3 there is one unmarked two-edge graph class but two
tournament classes, so a fixed-length tuple of such same-size graphs cannot
encode both. Asymptotically a fixed k-tuple with N vertices per sparse graph
needs N >= (1/(4k)+o(1))n²/log2(n). For a retained connected core of cycle rank
beta and odd symmetry, disjoint exact-symmetry repair needs Omega(beta log beta)
extra vertices; the explicit rigid-tree padding attains this order.

## 3. Reversal markings have exact transports, sometimes with nontrivial loops

Given two anti-involutions s,t of a self-converse tournament, the permutation

    g(t,s) = sqrt(t s)

conjugates s to t. Its odd square root is computed cyclewise in linear time,
assuming the maps are supplied. This repairs an already found isomorphism
without searching the automorphism group; finding that isomorphism remains a
separate problem.

[Transport and holonomy](carriers/TRANSPORT_AND_HOLONOMY.md) proves an exact
flatness criterion: all triangular compositions of this pair-local rule agree
with direct transport iff the inverted twisted subgroup is abelian.
The first odd-order counterexample is C3[C3,C3,C3] on9 vertices:

- full automorphism group order81;
-27 anti-involutions, each with centralizer order3;
-11664 of19683 ordered transport triangles have nonidentity error;
- the associated classical twisted-subgroup construction is a nonassociative
  Bruck loop of order27, with its identities checked directly.

A shared global base mark gives flat transports even here. Thus the obstruction
belongs to the particular pair-local convention, not to every conceivable
choice of coordinates. If that base is discarded, centralizer corrections must
be retained whenever later operations can see named ports.

The [iterated C3 family](carriers/ITERATED_TRIANGLE_FAMILY.md) turns this into
an all-order procedural hostile. For N=3^k,

    |Aut| = N |C|²,       number of anti-involutions = N |C|.

Even after naming the fixed vertex there remain |C| possible reversal markings:
243 at order27. Specifying the exact mark requires asymptotically
(N/4)log2(3) extra bits in this family. This is a statement about the exact-mark
observation task, not about every decision on the tournament.

The [mass calculation](carriers/MARKED_MASS.md) explains another loss: the lift
preserves classes but replaces full Aut by an anti-involution centralizer.
Uniform matrix inputs therefore bias classes by the number of anti-involutions.

## 4. Oddness is a descent hypothesis, not a numerical resemblance

[Fixed-orbit descent](oddness/FIXED_ORBIT_DESCENT.md) recovers the classical
coprime-action mechanism: for a p-group A acting compatibly with a p'-group G,

    X^A / C_G(A)  ->  (X/G)^A

is a bijection. The elementary proof works on both an orbit and its transporter
coset. It preserves classes, not the old stabilizers. C2 and C4 hostiles show
that existence and uniqueness can fail separately without coprimality.
A C3 two-port example shows why applying descent to separate marginals loses
relative alignment.

The Monster/Baby Monster work contributes real adjacent mechanisms: the
Monster7A local quotient7:6 contains the Paley7 symmetry7:3, and marking a2A
involution retains2.B. But the Monster has even order and multiple involution
classes, so tournament odd-kernel uniqueness does not transfer to it. The
older level-two supersingular Frobenius defect similarly distinguishes a
fixed base object from its nonfixed extra structure. These are scoped bridges,
not a moonshine or {7,21} spectrum proof.

In particular H, H/|Aut|, reversal-fixed path count h_tau, and its centralizer
quotient are four different statistics. The folded observer actually attains7;
it does not inherit the all-order7/21 exclusions for H.

## 5. One observer really can discard more information

The recovered pair-mask path DP computes h_tau in O(n²2^n) elementary transitions
on a2n+1-vertex lift. The sign choices of already used pairs need not be retained:
only the current endpoint can affect continuation. The recurrence was already
in the October20 corpus and is credited accordingly.

The [centralizer-orbit refinement](oddness/CENTRALIZER_ORBIT_PATH_DP.md) compresses
this observer exactly when the actual centralizer action is supplied. It uses
orbit totals and weighted transitions, including nonfree intermediate states.
On the nine-vertex wreath example,64 states become24 orbits and96 arcs become32
weighted entries; h_tau=33 and the blue quotient is11. Incorrect unit
initialization gives19. All531 small labelled lifts agree with the uncompressed
recurrence, but their aggregate reduction is modest; no general runtime win is
claimed.

Flipping the three internal mirror-pair arcs of Paley7 gives eight tournaments
with the same h_tau=9 but different full groups and blue counts. This is a
precise example of information that is unnecessary for one count yet necessary
for symmetry-aware reconstruction.

## Verification and remaining frontier

[Independent audit](audit/INDEPENDENT_AUDIT.md) records separate action conventions,
ordinary-bank small-order checks, fullS9 symmetry checks, every labelled tree
through7, rigid-tree recognition, integer Fourier matrices, and hostile tests.
All released Python checkers use active runtime checks under normal execution
and `python3 -O`. They are standard-library-only and generate receipts locally.
Their exact command scopes are documented beside each theorem.

The main open goals remain natural unlabelled object correspondences beyond
these obstructions, efficient general tournament classification, and minimal
closed carriers for larger operation languages. No new A000568 count formula,
new all-order prime atlas, polynomial GI/TI algorithm, LRC closure, or major
open conjecture is claimed. The gains are exact structural maps, a reusable
transport algorithm, an infinite hostile family, and a verified observer
compression rather than a numerical analogy.
