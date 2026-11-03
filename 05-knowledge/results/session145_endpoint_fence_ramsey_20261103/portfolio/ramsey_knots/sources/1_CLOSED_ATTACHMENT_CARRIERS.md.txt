# Closed attachment carriers: joint relations, actions, and exact loss

Session86, 2026-10-28. **PROVED elementary statements and scoped deductions;
FINITE-EXACT controls.** No priority claim, universal affine model, or new
graceful-tree/square-sum existence theorem. The inherited source snapshot is
`5ed1835235b36a57b8cd11eb705aa259b6bc0f08`.

## 1. Findings first

1. The closed object for ordinary labelled attachment counting is a **typed
   weighted joint relation**, not a list of endpoint marginals. Its composition
   is a compatible join followed by summation over the shared assignment.
2. There is no absolute smallest carrier. The exact minimum for a specified
   grammar is contextual equivalence: two states agree iff every allowed
   future context gives the same observation, including failure. If the grammar
   allows every named singleton boundary probe, the full weighted relation is
   already minimal. Without those probes, that conclusion does not follow.
3. A boundary table together with the action on its ports is still insufficient
   for stabilizers or unlabelled counts. Retain the **action on the actual
   witness fibers**, or its complete subgroup fixed-fiber tables. A three-point
   trivial C3-set and a regular C3-set have the same total mass and identical
   one-point boundary action, but three and one orbits respectively.
4. Frobenius and LCA can recover correlations when the supplied equations
   actually determine them. They cannot supply an unspecified phase, offset,
   branch decoration, label map, or nonlinear validity constraint. Session85's
   precise missing-data quotient is `W / sum_B (W intersect V_B)`.
5. The mandatory tournament hostile survives the stronger test: the 81/27 pair
   has the **same labelled cyclic-triple hypergraph**, not merely the same
   number 67, and H=462915 on both sides. Actual child decorations determine
   the different stabilizers and Hamilton-word orbit counts 5715/17145.

## 2. Declare the attachment language before minimizing

Each port has a name and a finite type/domain. For a finite named interface B,
let `D_B` be its joint assignment set, including resource coordinates needed
by guards. An object consists of a finite set X of actual internal witnesses,
a restriction map `b:X -> D_B`, and weights in a commutative semiring K. Its
ordinary labelled response is

    R_X(d) = sum_(x in X, b(x)=d) weight(x).

Use the Boolean semiring for existence, natural numbers for multiplicity,
or a specified field when discussing linear dimension. Signed field weights
can cancel, so a zero field coefficient must not be read as nonexistence.

The primitive grammar is: named relabelling/type-preserving transport;
disjoint product; identification of designated common ports; an explicit
guard depending on exposed coordinates; and hiding/summing an interface.
For disjoint interior witnesses, two pieces with interfaces A+B and B+C compose
by

    (R ⋆ S)(a,c) = sum_(b in D_B) R(a,b) S(b,c).                 (1)

A cross-guard g contributes its explicit factor `g(a,b,c)`. A failed Boolean
guard yields the empty relation. Shared variables denote the **same** choice,
not fresh independent samples. If two pieces also share an unexposed witness
or resource, (1) is not justified until that identity is exposed. Counting a
shared choice once and assigning each multiplicative factor once are part of
the model, not optional corrections after multiplication.

**Closure theorem.** Finite weighted relations are closed under these
operations, and any parenthesization of a fixed compatible gluing gives the
same table.

**Proof.** A composite witness is a tuple of interior witnesses with equal
shared values and satisfied guards. Group the sum of the product of its
weights by the boundary assignment. This gives (1). Reassociating the finite
sum/product cannot change the result. Named transports simply permute its
coordinates. QED.

This is an exact sufficient carrier; it is not automatically compact. A table
can have size exponential in the number of ports. Nor does a table at one
fixed arity automatically support later exposure of a forgotten interior
vertex. Such an operation must already be in the grammar, with its necessary
state. In particular, complement/re-marking is not included for free merely
because ordinary attachment is closed.

## 3. Exact minimality, with its quantifiers visible

Fix the permitted contexts, a terminal observation O, and failure semantics.
For states of one type put

    x ~ y  iff  O(C[x]) = O(C[y]) for every permitted context C.   (2)

Contexts include the empty context and are closed under inserting any allowed
operation before a continuation.

**Context theorem.** Relation (2) is the coarsest operation-stable equivalence
that preserves O. Every exact deterministic carrier factors through a
refinement of (2); the equivalence classes themselves form an exact closed
carrier.

**Proof.** If x~y, appending an operation cannot distinguish them because any
subsequent context is also an original context. Conversely, any stable
O-preserving equivalence preserves O along every context by induction, and
therefore cannot identify a pair distinguished by (2). QED.

For weighted tables, suppose a context may pin every exposed coordinate to
any specified `d in D_B` and return the resulting K-weight. Its output is
exactly R(d). Thus distinct tables are separated by a single pin context.
This proves minimality among deterministic exact table summaries in this
**unrestricted named-probe grammar**. For Boolean relations it likewise
separates any two supports by membership of a tuple in their symmetric
difference. It is not a lower bound for a restricted arithmetic grammar that
cannot realize the probe, nor does it require storing a table literally.

There is a separate linear dimension statement. Over a field F, terminal
one-hole contexts in which the input table occurs exactly once define linear
functionals on `F^(D_B)`. If L is their span, the
smallest linear response representation has dimension dim L, namely the
quotient by their common annihilator. This is not the number of deterministic
states, not a bit count without a finite input universe, and not the number
of named guards. With all singleton probes, dim L=|D_B|. If the coefficient
domain has an explicit finite-string encoding (for example Boolean values,
natural numbers or rational numbers), one unbounded integer can encode a
whole finite table. This does not apply to arbitrary real-valued tables.
Counting scalar storage slots without a declared encoding model would
therefore be meaningless.
If a grammar duplicates the same input table, its observations can be
polynomial rather than linear; the stated linear dimension theorem does not
silently include such contexts.

The inherited Session70 results are different specializations of (2):

* Fixed named Hamilton path plus all named tournament re-marking guards:
  the complete free-arc presentation, `2^((n-1)(n-2)/2)` valid states.
  This uses its cited arc-traceability reconstruction input, not a new proof
  in this packet; an explicit failure sink is extra.
* r mandatory disjoint nontrivial paths with 2r distinct named ends and
  arbitrary fresh connector paths: exactly the endpoint perfect matching,
  with `(2r-1)!!` states.
* Its factorial graceful family, with named edge-overwrites, no-ops allowed,
  and only validity observed: exactly the color-to-component permutation;
  the independent left/right signs can genuinely be forgotten.

Changing the observer or forbidding no-ops changes the experiment. These
three answers cannot be substituted for one another.

## 4. Symmetry needs a witness action, not just a port action

In this section witnesses are unweighted finite sets (counting multiplicity
means retaining actual distinct witnesses). Suppose a specified finite group
G acts on the actual witness set X and on its boundary D, with the restriction
map equivariant. For every subgroup
H≤G and every H-fixed boundary assignment d, retain

    R_H(d) = |{x in X : b(x)=d and h x=x for all h in H}|.       (3)

For ordinary labelled counts use H={1}. For orbit counts, Burnside uses
the cyclic-subgroup cases; full stabilizer distributions need more than
elementwise fixed counts in general.

In fact, for the fixed-G equivariant grammar below and **only final orbit
counts**, the smaller family of element-fixed boundary tables is already
closed and sufficient. Full subgroup marks must not be called the minimum
for that weaker observer. They supply the stronger classification claimed
here.

**Equivariant closure theorem.** For gluing two such G-sets over the same
equivariant interface, using G-equivariant transports and G-invariant guards,
the H-fixed boundary tables compose by the same
formula (1), summing over **H-fixed** intermediate assignments. Hiding ports
likewise sums over fixed assignments. Thus the full family (3) is closed.

An arbitrary named pin is not G-invariant: its surviving symmetry is the
stabilizer of the pinned data, not automatically the original G. Restricting
the action to that subgroup is valid and is covered by the complete subgroup
tables. If numerical weights matter, invariant weight values may be retained
as extra colors, with separate fixed tables for each color. Weighted fixed
sums alone need not classify the underlying weighted action.

**Proof.** A pair in a fiber product is fixed by H iff both its coordinates
are fixed by H. The common boundary of an H-fixed witness is H-fixed by
equivariance. This gives the assertion by a bijection, not an independence
assumption about orbits. QED.

Moreover (3) determines the finite G-set over D up to equivariant isomorphism.
Choose one d in each G-orbit. Its fiber is a `G_d`-set. For subgroups H≤G_d,
its fixed counts are triangular sums of the orbit multiplicities
`G_d/K`: a term can occur only when H is conjugate into K, and the diagonal
at H=K is `|N_(G_d)(K)|/|K|`, nonzero. Descending through subgroup orders
therefore recovers all orbit multiplicities. Transport recovers the fibers
over the rest of D. This is the elementary marks argument; computation may
be expensive. Under the observer that permits each fixed-subgroup pin count,
different tables are distinguishable and the representation is minimal up
to equivariant isomorphism. A weaker observer may compress further.

**Small hostile.** Let D be one point and G=C3. Let X be three witnesses on
which G acts trivially, and Y its three-element regular orbit. Both ordinary
tables are `[3]`; the boundary actions are the same. Their element fixed-count
vectors are `(3,3,3)` and `(3,0,0)`, with respectively three and one orbits.
Their witness stabilizers are respectively C3 and the trivial group.
Therefore neither boundary mass nor boundary action recovers hidden symmetry.

**Element-count hostile.** Let G=C3×C3. Take X to be one regular nine-point
orbit together with three globally fixed points. Take Y to be the disjoint
union of G/H over the four order-three subgroups H. Both sets have twelve
points; every nonidentity g fixes exactly three points in each. Hence their
element fixed-count functions agree and both have four G-orbits. But X has
three globally fixed points and Y has none; their stabilizers have orders
1 or 9 versus always 3. Thus element counts suffice for Burnside without
classifying the action or its full stabilizer distribution. The independent
review supplied this odd-order control, and `check_attachment_carriers.py`
verifies it explicitly.

Even knowing all point stabilizers need not recover transport between points.
Let H=C3^3 rotate the triangles inside three four-vertex blocks, fixing each
block's singleton; let G=H semidirect C3 also rotate the three blocks. Every
element of G outside H moves every vertex to a different block. Consequently
`G_v=H_v` as literal permutation subgroups for every named vertex v, and the
same equality holds for pointwise stabilizers of every nonempty ordered port
tuple. Nevertheless G has transporters between corresponding ports in
different blocks and H does not. Thus their port orbit partitions differ.
This is exactly the 81/27 pair's action mechanism, not merely a different
abstract toy group. It concerns ordered roles / pointwise stabilizers;
unordered setwise stabilizers need not agree. A list of local isotropy groups
is not a replacement for the actual common permutation action or transporter
relation.

This theorem assumes the **actual common group and transports are given**.
It does not recover an unknown parent automorphism group from child group
orders. If pieces have different groups, the allowed common-frame action,
embeddings and parent transport must first be specified. Retaining the full
witness groupoid is always sufficient, but calling it a smaller numerical
invariant would conceal the data it contains.

## 5. Precisely when Frobenius constraints give joint data for free

Let `C=x0+K` be a nonempty affine family of explicitly identified binary
coordinates, and `W=K^perp`. For retained coordinate scopes B define

    W_local = sum_B (W intersect V_B),
    J = join_B projection_B(C).

Session85 proves, with right-hand sides inherited from x0,

    J = x0 + W_local^perp,
    |J|/|C| = 2^(dim W - dim W_local).                         (4)

Thus the missing joint data are the quotient `W/W_local`, not the absence of
an LCA label. A supplied family of valid additional parity equations with
row span F repairs the projected join exactly iff

    W_local + F = W.                                         (5)

Here F⊆W and every right-hand side must be retained; inconsistent equations
give the empty set, not a negative free-bit count. More generally, under
observed linear forms O a target linear value q is determined exactly when
`q in W+O`. Joint targets Q retain
`dim((W+O+span Q)/(W+O))` unknown binary dimensions on each nonempty uniform
affine fiber. These assertions are representation-relative linear algebra,
not an entropy model for arbitrary graphs or tournaments.

**What Frobenius can actually select.** The inherited marked field statement
retains `E=F_(q^2)`, `K=F_q`, the spread of K-lines, and absolute Frobenius.
Alternating binary forms isotropic on every K-line are exactly

    B_c(x,y)=Tr_(E/F2)(c x y^q), c in K.

Frobenius invariance gives c=sqrt(c), so c∈F2; nondegeneracy excludes zero.
Hence the unique form is B_1. This is genuine reconstruction from the
retained structures. If an attachment law is **defined** to be B_1(u,v)=0
and its ports are already named elements u,v of E, its entire joint
compatibility matrix is indeed free to compute: no new choice bit is needed.

But determining B_1 is not determining which vectors an unnamed pair of
ports denotes, an arbitrary subset of the B_1-compatible pairs, the witness
multiplicities, or the action on that fiber. It does not select an absolute
normal-basis phase or an arithmetic prime labeling. The earlier proof says
this explicitly. For F16, Frobenius alone admits two nondegenerate forms;
the F4 spread alone admits three; their marked intersection has one.

**Distinct-port Frobenius-phase hostile.** Work in
`F16=F2[alpha]/(alpha^4+alpha+1)`, with binary polynomial labels. The orbit
of alpha is `(2,4,3,5)`. Thus pairs `(2,3)` and `(2,5)` have the same separate
port-orbit labels and distinct nonzero ports, but

    B(2,3)=Tr(alpha^2)=0,   B(2,5)=Tr(alpha^3)=1.

Only the second port has changed phase, since `F(3)=5`. Simultaneous
Frobenius does preserve B, for example `B(4,5)=B(2,3)=0`. This is an explicit
failure of replacing a joint orbit by separate port orbits even **after**
the form has been uniquely reconstructed. Indeed
`(alpha+1)^4=alpha` and `(alpha^2+1)^4=alpha^2`; the displayed traces follow
by repeatedly substituting `alpha^4=alpha+1`. Independent polynomial
convolution/long-division arithmetic confirms both values.

**Cyclic symmetry hostile.** On three bits the even and odd parity relations
are both invariant under cyclic coordinate rotation. They have identical
one- and two-port marginals, identical LCA rank invoices, and the same
homogeneous parity row `(1,1,1)`. Yet their joint supports are disjoint: the
right-hand side is the missing bit. Symmetry of the family or knowledge of
the equation's form does not choose that offset. All proper marginals of
even/odd parity agree in every arity at least two.

**LCA's exact role.** With actual leaf labels, quotient orientations or
transport maps retained, LCA locates the interface at which two branches
meet. For linear transports it combines stored maps along
`u -> LCA(u,v) -> v`, with inverses and order. It does not manufacture missing
map coefficients or certify that a discarded cross-branch constraint held.
Assigning each parity row to its support's LCA invoices only new rank, and
does not make an arbitrary guard linear. Session56's running-intersection
theorem assembles the retained constraints exactly; it does not restore a
constraint that was first projected away.

## 6. The three applications, with native guards retained

### Square-sum paths

For Q_n the vertices are the actual integers 1,...,n and an edge ab requires
`a != b` and a+b to be an integer square. A selected partial Hamilton object
is a linear forest. With all vertices fixed and 2r distinct exposed endpoints
of nontrivial components, its exact unrestricted Boolean attachment state
is the endpoint matching M. Joining a to b is legal iff they are in distinct
matched pairs; replace `aa'` and `bb'` by `a'b'`. A connector between partners
closes a cycle. Different matchings are separated by choosing a pair present
in one and not the other. This proves the inherited `(2r-1)!!` minimum.

Arithmetic attachments must additionally carry the **integer endpoint names**
and square-sum legality. Singleton components require their degree-zero /
coincident-end convention, not two fictitiously distinct vertices. A family
of covers needs joint multiplicities. Degree counts, component sizes, or
independent endpoint histograms lose which ends belong together.

For the restricted birth operation adding the prescribed next integer, all
perfect matchings and all distinguishing connectors need not be realizable.
The matching table is sufficient, but the unrestricted minimum is not a
proved minimum for that arithmetic language. Nor does a graph cycle-space
or incidence rank test enforce the final single-component Hamilton condition.

### Graceful trees

A prescribed m-edge tree T is gracefully labelled by a bijection
`x:V(T)->{0,...,m}` whose edge differences `|x_u-x_v|` are exactly 1,...,m.
Thus joint boundary data must support label equality at glued vertices,
noncollision of labels elsewhere, nonreuse of difference colors, and the
specified shape/forest condition. Partial-labeling tables may expose used
label/color sets as finite resource coordinates, rather than pretending
all-different constraints are linear. If tree shape is fixed, retain it;
if edge selections are allowed to change shape, forest component incidence
is a separate state coordinate.

For the inherited rainbow-tree formulation with fixed labelled forest F,
contract F's components but retain actual colored candidate contacts,
including loops and parallel edges. Its exact graphic closure is the
component partition, not just the number of components. The joint choice
of one edge per color is still a rainbow constraint. The factorial-family
validity quotient is an example where signs disappear under its exact
overwrite grammar; it is not a proof that label signs can always be erased.
The actual 13-vertex frozen-color example has five completions split 3+2
under width-two moves; this already refutes the unrestricted frozen-color
connectivity claim, not the Graceful Tree Conjecture.

### Tournament substitution and re-marking

For `Q[T_1,...,T_k]`, arcs inside children are their actual arcs, and every
arc between children is determined by Q. Hamilton paths can revisit a child
in several runs. An ordinary composition response must therefore retain
joint ordered local path-cover roles and their multiplicities, then sum
over compatible quotient-walk words. A one-entry/one-exit scalar or one
visit to each quotient vertex is insufficient. Session69 already proves
the all-arity uncoloured response recursion; it is inherited, not claimed
as a new closure result here.

If the language adds complement and named re-marking, retain their actual
aligned source-target transport or the full named tournament presentation.
Independent endpoint-rooted types discard relative alignment; resetting
ports by fresh sampling changes the language. Session70 shows that all
named legality guards recover the entire tournament presentation. A weaker
counted-multiset or normalized-distribution observer has different quotients.

The mandatory pair is

    T+ = C3[A+,A+,A+],  T- = C3[A+,A+,A-],
    A+ = C3 -> 1,       A- = (A+)^op.

Both have the same 67 labelled cyclic triples and H=462915. Their actual
ordinary groups have orders 81 and 27: three identical decorated children
admit the outer C3 rotation; one odd child breaks it, while the three inner
C3 rotations remain. Automorphisms act freely on ordered Hamilton words,
so the word-orbit counts are 5715 and 17145. Keeping children only up to
converse, a form, total H, or the cyclic-triple hypergraph loses that choice.
The full H-fiber does determine the possible choices and their group action
(Session81), but not which realization was selected.

**Label convention.** `check_attachment_carriers.py` uses local vertices
0,1,2 for the directed triangle and local vertex 3 for A+'s singleton sink;
the third block is reversed in T-. The separate `81_27_TRANSPORTERS.md` convention
instead starts with the singleton source at local vertex 0 and triangle
vertices 1,2,3. These are different **converse and label conventions**, not
merely different vertex names: the all-source and all-sink tournaments have
different degree multisets. No named-port statistic is transferred between
them here. In particular, the named-port-8 Hamilton-count split in
`81_27_TRANSPORTERS.md` belongs to its presentation, not silently to this checker's
port 8. This checker certifies its unmarked 81/27 pair and does not assert
that named-marker count for its different converse convention.

## 7. The signed-incidence bridge is exact, but not an equivalence of problems

For a directed edge use

    x_v = s_e x_u + w_e,  s_e in {+1,-1}.                     (6)

Signed differences have s_e=+1; prescribed square sums have s_e=-1 and
`w_e=k_e^2`. Tree propagation gives `x_v=rho_v z+t_v`; each chord gives

    (rho_v-s_e rho_u)z = w_e-t_v+s_e t_u.                      (7)

The coefficient is zero or ±2. Balanced cycles test consistency; an
unbalanced cycle pins z, with an integer parity condition and consistency
against every other chord. For all-negative sums these are even/odd cycles.
Reducing modulo two merges + and - and can erase the factor-two obstruction.

On a bipartite graph with parts U,V put y_u=x_u and y_v=-x_v. For an edge
u∈U,v∈V its sum x_u+x_v becomes the oriented difference y_u-y_v. This is an
exact change of linear coordinates. It also transforms every side guard:
positivity becomes y_U>0,y_V<0, injectivity becomes signed-coordinate
noncollision, and prescribed square values remain prescribed squares. It
does not become a graceful permutation of 0,...,m. Conversely, gracefulness
requires absolute differences; oriented values -1,+1 are distinct but have
the same forbidden repeated absolute color 1. Retaining signs/mirror factors
is essential. No rank or Smith-form statement proves the nonlinear label
guards, Hamiltonicity, or existence for every abstract tree.

## 8. Relations becoming objects does not erase their interfaces

The earlier 7→15 construction takes the four-dimensional Hamming7 kernel,
chooses an actual basis, and makes its 15 nonzero relations into the 15
points of F2^4, with the declared fifteen-prime realization. Their 35
three-term laws become the next incidence relations. This is a genuine
relation-to-object map, not merely a matching count.

For attachment purposes retain the basis/transport, coefficient ring,
boundary coordinate identifications, and any further form or decoration.
The new Hamming15 kernel has dimension 11, not 4; it is not the old object
under an unnamed numerical renaming. Integer subset-module types and
presented maps distinguish modular relations arising from actual integer
kernels from those arising from cokernel torsion. Keeping only rank,
nullity, a full cokernel, or an abstract 15-element set does not retain the
composition law on particular named ports.

## 9. Reproduction and stopping boundary

    python3 -B check_attachment_carriers.py > attachment-normal.json
    python3 -B -O check_attachment_carriers.py > attachment-optimized.json
    cmp attachment-normal.json attachment-optimized.json

The two runs agree byte for byte: **12,775 explicit checks, PASS**. The
self-contained checker uses only the Python standard library, requires no
pre-existing JSON inputs, and uses no removable asserts.
It checks all 4,096 binary 2-by-2 table composition triples over integer
counting, all singleton probes on ternary-valued tables, complete parity
hostiles for arities 3,4,5, all endpoint matchings through eight ports and
their legal updates, the hidden-C3-action and C3-squared character hostiles,
and the literal 12-vertex
81/27 pair using separate subset Hamilton DP and unrestricted arc-preserving
permutation backtracking. These finite checks support the displayed general
proofs; they do not prove missing global conjectures.

Closed here: exact relational/equivariant closure, observer-relative
minimality criteria, and the precise linear test for whether extra equations
repair lost joint data. Not closed: optimal carriers for restricted arithmetic
birth languages, a compact universal complement/re-marking aggregate algebra,
or arbitrary-tournament affine coordinates. Those require a specified new
operation grammar or additional mathematical input.

## Source ledger (all repository-relative, inherited)

* `05-knowledge/results/relation_gluing_7_15_20261028/THEOREMS.md`, T1–T4:
  kernel/image gluing, relation-to-object map, coefficient change.
* `05-knowledge/results/free_bits_signed_incidence_20261028/THEORY.md`,
  P3–P4: marked subfield/Frobenius form reconstruction, phase boundary,
  affine signed-incidence/LCA transport.
* `05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/LINEAR_GLUING_INFORMATION.md`,
  §§2–5: exact scope loss, separator invoice, affine messages, LCA row rank.
* `05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/SIGNED_INCIDENCE_RELATIONS.md`:
  signed parameters, torsion, typed factor-node laws versus bare graph cycles.
* `05-knowledge/results/session56_rmq_lca_connections_20261025/correlated-gluing/THEOREM.md`
  and `tournament-hierarchy/REPORT.md`: compatible identities, parity hostile,
  actual quotient data at LCAs.
* `05-knowledge/results/session70_attachment_guard_alignment_20261026/session70-root/SYNTHESIS.md`,
  `session70-graph-labeling/REPORT.md`, `session70-guards/REPORT.md`:
  three distinct minimality statements and native graph guards.
* `05-knowledge/results/fresh_remarking_quotient_and_divisor_arrangements_20261026.md`,
  §§1–4: named versus fresh protocols, 81/27 phase, Vitali observer history.
* `05-knowledge/results/session81_counterexample_families_20261028/FIBER_ACTION.md`:
  actual signed realization-fiber action and stabilizers.

The Vitali contribution used here is its finite observer-fiber/continuation
principle and native correlation hostiles, not measure-theoretic
nonmeasurability. The note does not identify every historical object called
a Vitali atom with one mathematical construction.
