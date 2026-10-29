# Moving selectors, ordered runs, and the minimum appropriate to each observer

2026-10-28. Local results C1--C6. **PROVED** means the arguments in this
file, not a formalization or an independent-agent review. Enumerated numerical
claims are **FINITE-EXACT**. There is no claim of external priority.

## C1. An exact carrier for a single read-only moving selector

Start with the [singleton language M1](../minimal_blue_selectors_20261028/THEORY.md).
X_n consists of tournaments with the identity Hamilton path. Tau flips the
non-path arcs. J_L/J_R insert a universal source/sink. The old blue observer
has minimal residual B(T), with4n-8 states for n>=5.

Add one cursor i. It can move to the preceding or next path position, with
observable boundary failure. It retains its vertex under tau and insertion.
The extra observation is the isomorphism type

    R_i(T) = [T; first, last, cursor=i].

Marks are ordered and may coincide. Fresh re-marking chooses a uniformly
random Hamilton path and retains the cursor vertex, now at its new position;
it forgets the old endpoint roles. Counted re-marking is also allowed.
The current index is observed. This is a read-only selector: it does not
deposit new markers at previously visited vertices.

**Theorem C1 (PROVED).** The exact minimal behavioral carrier is

    K_scan(T,i) = (i, (R_j(T),R_j(tau T))_(0<=j<n), B(T)).

Necessity: move the cursor to every j and observe before/after tau. Every
entry of the displayed profile is therefore observable. Distinct B states
are separated by the old pure-singleton tests. Sufficiency: tau exchanges
the two arrays. For J_L, each rooted source type gains a universal source;
each rooted complemented type gains a vertex beating only the old first
port. All required old ports are among first,last,j. The new cursor's own
entry uses just the old endpoint type. J_R is dual. B updates by M1.
For reset, the whole next-state distribution depends only on [T;cursor],
obtained by forgetting the endpoint marks from the observed R_i. Rooted
isomorphisms biject the eligible paths and the new cursor positions. Thus
the laws, including every later adaptive continuation, descend exactly.

This need not store n copies of a relation. Given E(T)=[T;first,last], each
R_j is just the orbit of vertex j under the endpoint stabilizer. Store one
rooted model and the ordered list of its vertex-orbit IDs; do the same for
tau T. The two lists cost at most2n ceil(log2 n) bits in addition to the
rooted relations and B. These lists preserve marginal orbit information,
not a common transporter for all positions.

If either endpoint stabilizer is trivial, the corresponding list assigns
every position a unique vertex in its rooted model, reconstructing the full
presentation. Whether K_scan reconstructs X_n for every n is **OPEN here**.
The ordered-run counterexample in C2 is a different input interface and
must not be misreported as a counterexample to this Hamilton-presentation claim.

**The old nine-vertex obstruction is separated.** Pad the1072/4176 seven-pair
with one tau J_L tau and one tau J_R tau. Positions0 and1 have identical
cursor observations on both complement sheets. Position2 separates after
tau. Thus from the first vertex the cursor needs two advances, or one from
the formerly fixed successor. More generally padding k vertices on each
side gives first separation at position k+1. Before then all queried vertices
are exterior or the original first endpoint; the old maps fix them. At k+1
the complemented seven-map moves the old second vertex. This proves the
claim for every k, with direct checks for k=0,...,4.

FINITE-EXACT: the39 Hamilton presentations of C3[C3,C3,C3] also all have
distinct K_scan profiles. This is one nine-vertex source class, not a census
of all191536 unlabelled nine-vertex tournaments.

## C2. A moving selector on child runs still misses a joint bit

Now use the general composition interface. A child comes with an ordered
directed spanning path cover, possibly consisting entirely of singleton runs.
A read-only cursor scans the ordered run list. At each position it observes
the whole child relation rooted at the two outer ports and the current run
port; it may also observe the complemented child relation with these marks.
It does not retain a previous interior port as another marked vertex.

Let A=C3[C3,1,1], with vertices0,1,2 in the cyclic child and singletons3,4:

    A rows = (10,12,9,16,7),
    Aut(A) = <(0 1 2)>, fixing3 and4.

Take the all-singleton covers in orders

    U=(3,0,1,2,4),       V=(3,0,2,1,4).

All five singleton runs are legal. Their forest complement is simply A^op,
whose automorphism group is the same. Both outer ports are3,4. At each
cursor index the selected middle vertices are in the same orbit; therefore
every finite adaptive read-only scan and complement sequence has the same
observations. Yet no one automorphism maps the whole U to V: the required
permutation is the transposition(1 2).

**An actual composition detects the loss.** Use prime quotient C3 and three
strong children A, D, D, where D has rows(2,4,9,3). Label the latter children
5,...,8 and9,...,12. For W=U or V, assemble the legal Hamilton path

    (W_0,5,9,W_1,6,10,W_2,7,11,W_3,8,12,W_4).

Both paths have quotient walk (A,B,C)^4 A, where B and C are the two distinct
quotient vertices carrying copies of D. Every run
is a singleton. Complementing along these paths gives strong tournaments
with respectively

    H = 2,899,675 and 2,906,091.                        (FINITE-EXACT)

Two independent forward/backward subset recurrences agree on both integers.
Thus the complete uncoloured metagraph targets differ, not just an invented
port statistic. Source and target rows and paths are in probe-receipt.json.
Replacing the two D children by transitive4 gives the separate positive
control2,483,993 versus2,491,877; a canonicalizer also separates that pair.

**Exactly one extra joint bit repairs this fiber.** With outer ports3,4 fixed,
the six orders of0,1,2 form two orbits under C3: even versus odd permutation
relative to(0,1,2). This relative parity is a complete discriminator here.
It can equivalently be read as the direction between the first two interior
run ports. All six orders are checked: each parity yields its corresponding
H value above. A retained bookmark at the first interior port also repairs
the fiber, because it kills the residual C3 stabilizer. Calling that bookmark
"free" without retaining it through later operations would lose the same bit.

Five is the least strong child size for this particular obstruction with
two distinct fixed outer vertices and all-singleton runs. At3 the cyclic
group is regular. A nontrivial automorphism of a strong4-tournament would
have a3-cycle and one fixed vertex; that vertex must see the cycle uniformly,
making it a source or sink, contrary to strongness. Thus the endpoint
stabilizers at smaller strong orders are trivial. This does not claim13 is
the least possible total size for every composition interface.

The lesson is not that a cursor cannot gather information. A cursor that
leaves persistent attachment gadgets has written memory. C2 refutes
composition from the list of isolated cursor responses after their relative
vertex identities have been forgotten.

## C3. The exact order fiber and the minimum anchors for marginal probes

Retain the inherited aligned relation carrier: A and its forest-complement
B on the same vertex set, modulo simultaneous isomorphism. Their agreeing
arcs recover a directed spanning linear forest F. The directed paths of F
recover every run internally. Only the ordering of its k components remains
to be specified. Set G=Aut(A,B), acting on those components.

**Theorem C3a (PROVED).** G acts faithfully on the components, has odd order,
and the number of ordered-cover classes over the given aligned pair is

    k! / |G|.                                         (C3a)

An automorphism fixing every component setwise fixes each directed path
pointwise, hence every vertex. G is a subgroup of a tournament automorphism
group, so has odd order. Its action on complete component orderings is free,
giving the formula. The minimal fixed-length discriminator for faithful
cover classification over this pair has ceil(log2(k!/|G|)) bits.

All component permutations in G are even. Consequently relative ordering
parity is well-defined on the quotient; for k>=2 the two parity sheets have
equal cardinality. Swapping two positions pairs the sheets. Parity depends
on a reference ordering: intrinsically it is a two-element choice with no
preferred zero. It need not classify the remaining k!/(2|G|) possibilities.

This refines the inherited aligned-pair/ordered-endpoint carrier
[C1](../strong_atom_continuation_vitali_20261026/README.md). It removes
redundant internal path words and endpoint lists because F already determines
them, and replaces the list by its component-order orbit. Complement swaps
A,B and preserves F and G. Re-marking changes the pair and requires a new
component action. Assembly expands the chosen ordered runs and adds the
quotient hinges with multiplicity. These are the inherited exact operations.

This is a sufficient closed carrier and an exact minimum for **faithful
ordered-cover classification over a fixed aligned pair**. For contexts
restricted to quotient-walk assembly of supplied runs and complement, we
have not proved that every distinct such class is behaviorally distinguished.
That restricted minimum remains OPEN. C2 shows a concrete fiber where both
classes are observable. This is not an openness claim for every richer
operation menu: incoming arbitrary-position C3 substitutions force complete
Hamilton-presentation memory, as audited in CONNECTIONS.md.

**Theorem C3b (PROVED).** Suppose some component anchors are retained, with
pointwise stabilizer H<=G. The sequence of individual anchored-component
orbit responses classifies every complete ordering of the remaining
components if and only if H=1.

If H=1, each individual response identifies its component. Conversely, a
nontrivial H has an orbit of odd size at least3. Swap two components in that
orbit in a complete ordering. Every marginal response is unchanged, but
the full orderings cannot be H-equivalent: their unique required permutation
is a transposition, of order2, and cannot belong to H. Thus the least number
of retained anchors for universally sufficient individual probes is the
permutation-group base size, relative to any anchors already fixed.

This is a scope-specific necessity statement, not a claim that storing a
base is the fewest-bit encoding of the entire order orbit. C2 needs one bit
for its two possibilities, even though a literal vertex bookmark has a larger
raw address. Bases and stabilizer chains are classical; this application
quantifies the missing component-order information.

FINITE-EXACT:593 labelled forests on all strong sources through5 (including
the singleton) satisfy C3a, faithful component action and the two parity
sheets. Minimum vertex bases for all-singleton forests are1 for A,2 for
Paley7,3 for wreath9, and3 for each81/27 hostile. The latter all-singleton
order fibers have5,913,600 and17,740,800 classes, respectively. The original
binary-attachment176/512 and transporter27/0 controls are also rerun.

## C4. Growing memory, Frobenius phases, and LCA

For W_l=C3[W_(l-1),W_(l-1),W_(l-1)] on N=3^l vertices, use the full group
proved in the incoming [wreath-family theorem](../session90_even_graphs_odd_descent_20261029/carriers/ITERATED_TRIANGLE_FAMILY.md).
Its vertex-action base size is exactly N/3 for l>=1. Each bottom triangle
has an independent rotation supported there, so a base must meet all N/3
triangles. Conversely choosing one vertex in each bottom triangle kills
their rotations; the induced upper block rotations are then killed
recursively. This proves both directions without enumerating the large groups.

On all-singleton forests this is also the component-anchor requirement in
C3b. After one named vertex is supplied, N/3-1 further suitably chosen
vertices suffice and are necessary. Compare three explicitly different tasks:

* The singleton blue verdict needs only4N-8 residual states.
* Choosing an actual reversal mark at a named pole requires distinguishing
  3^((N-2l-1)/4) possibilities in the incoming family:243 at N=27.
* Classifying a complete component order over (W_l,W_l^op) requires
  N!/3^((N-1)/2) orbit values. One-point queries need a base to recover these
  joint values from marginal ones.

The orders and information tasks differ. A small folded-path or blue
observer does not compress every frame or component ordering.

On free frames, the inherited LCA rule retains g_i^(-1)g_j on tree edges.
Multiplication along the unique path recovers every relative alignment.
C3b specifies when a component-anchor tuple is a free frame, hence when this
shortcut is justified. Without that condition, pairwise orbit data can
select incompatible intermediate witnesses. Replacing a group by its order
does not repair the problem, as81/27 continues to demonstrate.

The [incoming square-root rule](../session90_even_graphs_odd_descent_20261029/carriers/TRANSPORT_AND_HOLONOMY.md)
g(t,s)=sqrt(ts) transports reversal marks exactly. Its triangle error
h=g(u,s)^(-1)g(u,t)g(t,s) lies in a centralizer and may be nonidentity.
Our replayed wreath9 example gives

    h=(0,1,2,4,5,3,8,6,7),

two3-cycles. Its permutation parity is even. Thus the cheap parity bit from
C3 cannot detect this phase error; actual transport remains necessary.
Permutation composition is associative. The Bruck loop of order27 arises
from the projected square-root operation on reversal-frame coordinates,
not from a failure of associativity of actual graph substitution.

After choosing a normal basis, Frobenius is a cyclic permutation of basis
ports. Such an action on distinct tournament ports can preserve orientation
only when every cycle is odd. Therefore for degree6 the possible shifts
are0,2,4, and for degree8 only0. We recheck these permutation guards; the
actual F64/F256 field models are inherited from
[N3](../strong_attachment_tail_compiler_20261028/THEORY.md). Allowed phases
are even permutations, so again relative order parity cannot recover them.
The safe compression is to use an actual free frame and tree transport,
not to identify a Frobenius phase with the blue polarity or order-parity bit.

There is also a literal finite-field instance of the lost joint alignment.
Construct Paley27 on F3[a]/(a^3+2a+1), orienting u->v when v-u is a nonzero
square. Frobenius x->x^3 is an order3 automorphism fixing0,1,2 and cycling
a,a+2,a+1. In a complete component order with outer ports0,1, interchange
a and a+2. Every individual rooted-port query still agrees, witnessed by
an appropriate Frobenius power for that one port. No single automorphism
implements the whole interchange, because the required map is a transposition.
All singleton runs are legal and their complemented relation is the converse,
so both relation sheets have the same marginal ambiguity. Once a is also
retained, the **displayed Frobenius subgroup** has trivial stabilizer. The
code checks the field, orientations, strongness and maps. We do not assume
this subgroup is the full automorphism group, or claim a minimal graph gadget
for this27-vertex example. This Paley27 is not the iterated W_3 tournament.

## C5. The folded-path observer has a smaller exact minimum

Use the inherited folded DAG for an odd tournament W with a supplied
anti-involution. A state (S,v) records used mirror pairs and the last vertex.
It moves to an unused signed vertex w when v->w; after all pairs are used,
acceptance means v->pole. Retain level |S| and this **final** acceptance flag.
No interior vertex identity is observed unless explicitly added.

Among partitions preserving these observations and weighted one-step
transition counts into each cell, there is a unique coarsest one. Start with
equal level/acceptance observations, then repeatedly split by the vector
of outgoing multiplicities into current cells. Stability is exactly
equitability. Any admissible partition refines the initial one and, by
induction, every refinement step; therefore the stable partition is minimal.
This is the standard equitable-partition mechanism, not a new general
minimization theorem. On this acyclic DAG an independent reverse-level
signature calculation obtains the same partition directly.

For wreath9 the exact chain is

    64 labelled states -> 24 centralizer orbits -> 16 equitable cells.

The16 cells have3,7,4,2 cells at levels1,2,3,4 and22 nonzero weighted
transition entries. Initialize each level-one cell with its number of
labelled states, then propagate its total mass using its representative's
outgoing multiplicities. This preserves the count33, and division by the
full centralizer order3 still gives11 blue presentation orbits. This final
division is not a division of every intermediate cell size.

The16-state quotient is minimal in the declared partition category, not
among arbitrary arithmetic circuits or all linear representations computing
the final scalar. Adding the Boolean probe [v=3] at every state raises the
minimum to35 cells. Allowing a selector to test [v=a] for every a, equivalently
observing v itself, restores all64 cells in this example. These are exact
guard costs, rather than extrapolation from automorphism orbit counts.

FINITE-EXACT: both refinement implementations agree on all531 normalized
lifts with core size0..3 and on wreath9. Aggregate labelled, group-orbit and
minimal equitable states are12420,12292,5396. All folded totals agree with
the inherited independently accumulated labelled recurrence. The wrong
unit initialization of the24 orbit cells still gives19 instead of33.

## C6. Tree trijection conjugates complement, with a different operation menu

Fix an oriented spanning tree F. The inherited end-edge bijection sends each
non-tree vertex pair to the first and last edges of its unique F-path. It
identifies the free completion with a tournament on E(F), equivariantly for
the marked tree's automorphisms. Define C_F to reverse all non-tree arcs.
Every free coordinate flips, so

    tree_codec(C_F T) = converse(tree_codec(T)).

This is an exact operation correspondence. Re-marking the tree requires
changing the chart and retaining its induced orientation; forgetting that
chart does not preserve Hamilton counts or attachment data.

For a deliberately larger language, allow C_F for **every underlying spanning
tree**, oriented as induced by the current tournament. No Hamilton-path
legality constraint is imposed. On the labelled m=C(n,2) edge bits these are
translations by the characteristic vectors of E(K_n) minus F.

**Theorem C6 (PROVED).** For n>=2 their span is all F2^m when n=0 or3 modulo4,
and the even-weight subspace when n=1 or2 modulo4. Consequently the labelled
operation graph has respectively one or two components; its unlabelled
quotient has one component.

For n>=3 any two edges e,f of K_n lie on a common cycle. Complete that cycle
minus f to a spanning tree containing e; exchanging e for f produces two
tree-complement vectors whose difference is e+f. These differences span the
even-weight subspace. Every generator has weight C(n-1,2), odd precisely at
n=0,3 modulo4, proving the span statement. The n=2 generator is zero and
gives the stated even subspace. When there are two labelled parity sheets,
an odd relabeling exchanges them: the parity change of the forward-edge
bit sum is the sign of the relabeling. Hence the unlabelled quotient is
connected. The singleton is trivial separately.

FINITE-EXACT: all18248 labelled trees at orders2..7 give the predicted ranks;
all1050 tournament/tree pairs through4 verify the codec/complement identity.
This theorem is **not** a connectivity theorem for the actual Hamilton
metagraph: there the path must be legal, and the translation is unavailable
when its realization guard fails. The larger tree language identifies
exactly which simplification would discard those guards.
