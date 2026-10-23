# A marked overlap must detect its loop group

23 October 2026. **PROVED, scoped:** the base-marked overlap criterion,
the explicit three-cell hostile, and information bounds for the inherited
zero-slot recolouring family. **FINITE-EXACT:** the declared small group,
colour, and shared-boundary controls. **OPEN:** a uniformly small observer
for arbitrary tournament wiring. These are elementary mechanisms and new
connections within this research packet, not claims of literature novelty.

The useful connection is between two kinds of middle data. The tournament
seven-term cover needs enough overlap information to synchronize its local
charts. The Fibonacci rows need enough occurrence information to synchronize
their shared coloured diagonals. In both cases a nonempty overlap can retain
the visible object while failing to retain the choice that produced it.

## 1. Inheritance and exact recursion types

Closest proved mechanism: [THM-801, three-face Cech descent](../../../../01-canon/theorems/THM-801-mobius-cech-descent-and-pure-cubic-colour-law.md),
especially its nonempty-triple-overlap phase argument. The canonical hostile
is its order-five cover: 64 pairwise-compatible bare triples, only32 of
which have zero holonomy. The corrected near miss is calling all three
faces deletions of fixed original tournament vertices. The least-used
sidecar is the pointwise stabilizer of the actual overlap owners.

The current [partial-chart gluing rule](../../odd_cell_gluing_20261023/README.md)
already handles compatible functional attachments at all sizes. This note
does not replace that rule. It identifies exactly when a smaller collection
of overlap owners suffices to certify the chart consistency it needs.

The full staircase identity is

    F_n = 3 F_(n-1) - 3 F_(n-2) + F_(n-3),
    F_n = binom(n-1,2).

Its three size-(n-1) faces are the high, low, and middle gap-window faces.
The middle face is an honest tournament operation in path-gap coordinates,
but not induced deletion of one fixed vertex of the original tournament.
The exact object statement is agreement and gluing of binary assignments
on those overlapping faces. The cardinality identity alone does not glue
arbitrary isomorphism classes of three tournaments.

The different odd half-chart has sizes

    A,B: n-1;  C,D: n-2;  E,F: n-3;  G: n-4,
    h_n = A+B-C+D-E-F+G,
    h_n = floor((n-1)^2/4), n odd.

The equal sizes of C and D cancel numerically; their actual address slots
do not cancel. See [THM-550, half-tiling parity recurrence](../../../../01-canon/theorems/THM-550-tournament-half-tiling-parity-recurrence.md)
and [THM-553, tile-address recurrence](../../../../01-canon/theorems/THM-553-tile-address-recurrence.md).
The two seven-letter sign strings must therefore retain their own slot
types, not be substituted into each other as a numerical recurrence.

Our concept board is: (i) base-marked overlap, (ii) jointly chosen charts,
(iii) residual-cycle potentials, (iv) zero-length index slots, and
(v) arbitrary future queries at a shared port. Every result below states
which of these data it preserves.

## 2. Base-marked overlap criterion

Let a connected graph have bijective transports between copies of a fixed
finite set X. Reverse edges carry inverse transports. Choose a root and a
spanning tree. Each non-tree edge determines a based fundamental-cycle
permutation h_e of X. Suppose every h_e belongs to a known subgroup
K <= Sym(X). A set B of actual named root owners is a **base** for K if
the only element of K fixing every member of B is the identity.

**Theorem.** If B is a base for K and every fundamental-cycle transport
fixes B pointwise, all transports admit a common system of vertex charts.
Equivalently, all closed-walk transports are identity. Conversely, if B
is not a base for K, a three-vertex triangle with transports in K can pass
every B-owner loop check and still have no common chart system.

**Proof.** The hypotheses force each h_e to be the identity. Tree transport
then defines a chart at every vertex, and the non-tree edge equations hold
precisely because their h_e are identity. All closed walks consequently
have identity transport. Conversely choose nonidentity h in K fixing B,
put identity on two sides of a triangle and h on the third. It passes the
pointwise checks on B but its loop product is h. QED.

This criterion is sharp for the specified possible-defect group K. If the
actual set of possible defects is smaller than K, fewer marks might suffice.
The tree, the base owners, and each transport's actual domain must remain
part of a certificate. Partial maps require the inherited domain checks;
an undefined loop is not an identity loop. The result concerns compatible
charts, not the existence of a single parallel owner section.

### Why the old nonempty triple overlap worked

In THM-801 the only chart change is simultaneous binary complement. Its
C2 action on the value at any one actual tile is free: a nonidentity
complement fixes no bit. Thus one tile in the triple overlap is already a
base. The general replacement for "G is nonempty" is not simply "there is
one shared owner". It is "the retained owners detect every possible loop
defect". Several triangles also require a cycle-generating collection of
overlap constraints; arbitrary graph loops cannot be discarded merely
because some triangular checks pass.

For a cyclic regular-five cell the rotation group acts freely, so one
owner is a base. For the Paley-seven automorphism group

    H = { x -> a*x+b over F7 : a in {1,2,4} }, |H|=21,

one owner leaves three maps; two distinct owners leave only the identity.
Indeed fixing x and y gives (a-1)(x-y)=0, hence a=1 and b=0. Fixing x alone
allows b=(1-a)x for all three slopes. The nonidentity map x->2x fixes0.
This is exactly the group type surviving in the six high-symmetry cases
of the inherited order-31 host, not an identification of every seven-cell
attachment with this affine group.

## 3. A literal tournament hostile with shared owners

Take three named copies X_0,X_1,X_2 of P7. Within each copy orient u->v
when v-u is in {1,2,4} modulo7. For i<j put u in X_i -> v in X_j exactly
when v=t_ij(u); all other cross arcs point backwards. The maps t_ij are
P7 automorphisms. This produces a tournament whose cells have distinct
outdegrees5,10,15, so any isomorphism must preserve the three actual cells.

Define source S by t_01=t_12=t_02=id. Define target T by

    t'_01=t'_12=id,  t'_02(x)=2x.

An internal chart triple (g_0,g_1,g_2) must satisfy

    g_1=g_0,  g_2=g_1,  g_2=(x->2x) g_0.

No triple does. Each pair of cells separately is isomorphic, however:

| Actual owners pinned in each cell | Charts on pair01 | pair02 | pair12 | Joint triples |
|---|---:|---:|---:|---:|
| none |21|21|21|0|
|0|3|3|3|0|
|0 and1|1|0|1|0|

The source and target even have the same supplied owner section (0,0,0)
through all three transports. This section does not extend to a chart
isomorphism. The second distinct owner makes the offending two-cell block
fail locally. The witness is a21-vertex construction; it is not claimed
to be a new order-31 native Hamiltonian-port census or a minimum-order
counterexample among all tournaments.

The checker independently enumerates all7! permutations to recover exactly
the21 P7 automorphisms, then verifies the actual cross arcs. It exhausts
all group-chart choices, with0,1,2 retained owner marks. Thus the obstruction
is not merely an abstract permutation diagram.

### Flatness and marking cost are different questions

The defect group in Section2 certifies compatibility of the specified
relative transitions. Uniqueness of a surviving chart instead depends on
the **residual automorphism group**, which centralizes the retained loop
transports after tree synchronization. These groups must not be exchanged.

The same literal example makes the distinction sharp. For automorphisms
of the flat source S, a diagonal copy of all H21 survives. Its loop defect
is already identity, but recovering a unique chart still requires two
distinct root owners. For automorphisms of the twisted target T, the
surviving group is the three multiplications x->a*x, a in{1,2,4}:

| System | No root owner | root0 fixed | root1 fixed | roots0,1 fixed |
|---|---:|---:|---:|---:|
|flat S|21|3|3|1|
|twisted T|3|3|1|1|

Here tree compatibility determines the other cells' charts from the root,
so these really are root-owner costs. A nontrivial loop can reduce the
remaining marking cost; choosing an owner that its centralizer moves is
essential. No marking creates an S-to-T isomorphism: those comparison
charts were inconsistent before any marks were added. The second owner
in the earlier table simply exposes the incompatibility locally.

## 4. Residual potentials belong after compatible chart selection

The [projected transport theorem](../../session38_incidence_ports_electrical_20261023/theory/PROJECTED_TRANSPORT_UNIQUENESS.md)
tests an exact *abelian* statement: a specified vector observation is
constant on feasible integer transport tables iff every residual cycle
has zero observed cost. Its potential equation is additive. It remains
correct with its actual port channels, capacities, and chart choices fixed.

There is a different nonabelian equation at the chart level:

    additive wiring: cost(u,v)=p(v)-p(u),
    chart gluing:    t(u,v)=g(v) g(u)^(-1).

The second cannot in general be certified by abelianizing the first.
Inside H above, put a(x)=2x and t(x)=x+1. Direct composition gives

    a t a^(-1) t^(-1) = t != id.

Every homomorphism from H to an abelian group kills this commutator,
although the actual loop translates every owner and has no fixed point.
The slope register alone also vanishes on it. This is a cheap exact hostile
to replacing ordered chart transport by additive character totals.

The resulting protocol is precise: first retain a jointly compatible chart
choice, or a sound representation of all such choices; then build and test
the corresponding residual transport channels. A potential on the union of
independently allowed channels can mix incompatible chart choices. The
theorem does not claim that all additive observations are powerless: an
observation retaining full actual channel identity is a different object
from a homomorphic abelian quotient of H.

## 5. The zero slot permits positive information density

Use the incoming [zero-slot theorem](../../odd_cell_gluing_20261023/ZERO-SLOT.md):

    h=(1_A,0,1_B,1_C,2,3,5,...),
    W_4=CA, W_5=CAB, W_k=W_(k-1)W_(k-2), k>=6.

Its preferred limit begins CABCACABCABCACAB, exactly the current sixteen
letters supplied by the user. Every preferred occurrence of C can be
changed independently to B while keeping all finite rows compatible and
legally decomposable. The choice belongs to the actual occurrence, not
just the numerical weight of a Fibonacci atom. Keeping the same set of
chosen positions makes all shared diagonals agree.

The adjacency change is essential. In the unspaced alphabet1_A,1_B,1_C,
the atoms A and B are consecutive. In the zero-slot alphabet their indices
are0 and2, so B+A is legal. It can replace C+A inside a block of length2.
Erasing a zero-length interval geometrically does not authorize erasing
its index from the nonconsecutivity rule. The old two-or-three-fibre
theorem remains distinct from the new three-or-four-fibre theorem; see
[the earlier three-unit lift](../../reset_20260926_colours.md).

Let c(n) be the number of C occurrences in the preferred length-n prefix.
For k>=4 its full-block value is

    c(h_k)=F_(k-3),  h_k=F_(k-1),  F_0=0,F_1=1.

If Z(n) is the ordinary Zeckendorf support in indices3,4,..., then

    c(n)=1_(3 in Z(n)) + sum_(i in Z(n), i>=4) F_(i-3).       (1)

This follows from the incoming exact prefix concatenation, and gives an
integer computation with no floating-point floor decisions. Put
alpha=phi^(-2) and psi=-phi^(-1). Binet's identity gives, for i>=4,

    F_(i-3)-alpha F_(i-1) = -alpha psi^(i-3).

Summing the absolute errors over all possible indices bounds their sum by
1/phi. The exceptional unit term contributes at most1-alpha=1/phi. Hence

    |c(n)-alpha n| <= 2/phi,
    c(n)/n -> phi^(-2).                                  (2)

**Corollary (exact family-relative information cost).** There are exactly
2^c(n) distinct length-n prefixes in the independent recolouring family.
A lossless fixed-length binary code for these prefixes needs at least
c(n) bits. The prefix growth rate is phi^(-2) log2 per unit length.

Proof: the set of changed C positions is an arbitrary subset of c(n)
distinct places, recovered by inspecting the output at those places.
Thus it is a bijection with the Boolean cube of that dimension. Equation(2)
gives the rate. This is a statement about this declared family, not an
entropy calculation for every admissible Fibonacci expansion.

### Charge and an endpoint product still omit occurrence information

The zero-slot register says that each C-to-B replacement adds one unit
to the hidden refinement debt. If r replacements are made, remembering
that scalar debt distinguishes r but not the chosen positions. Assign
each letter an actual permutation of some finite group G and remember
the exact whole-word product as well. At fixed n the combined summary
has at most(c(n)+1)|G| values. Some fibre therefore has size at least

    ceil( 2^c(n) / ((c(n)+1)|G|) ).                       (3)

Even without the product, the exact-r fibre has binom(c(n),r) histories.
Any future query allowed to inspect an arbitrary actual position can
separate distinct histories. Thus neither finite phase nor total debt
can replace the occurrence sidecar for arbitrary future shared-port use.

At n=21 there are8 free positions and256 compatible histories. With the
three distinct S3 letter maps

    A=(0->1->2->0), B=(1 2), C=(0 1),

the history changing position6 and the history changing position9 have
the same replacement count1 and the same final product(0 1):

    CABCABABCABCACABCACAB,
    CABCACABBABCACABCACAB.

The largest combined debt/product fibre among all256 histories has36
members. These are complete-word product collisions, not failed grammar
checks; every one of their shorter rows is admissible and compatible.

## 6. Fixed shared ports permit exact compression

There is a positive counterpart to the information lower bound. Suppose
an actual bijection word is a path, and the positions that may be shared
with other pieces are known. Split the path at those positions and retain
the ordered product on every resulting interval, with its actual endpoint
sets and integer length. This exactly preserves all boundary-state and
pin queries at the declared shared positions.

For total bijections, those interval products are also necessary: fixing
an arbitrary input at the left end of one interval and asking for its
right value recovers the whole interval map. Two different interval maps
therefore cannot be merged by a summary that must answer every such query.
The checker tests all6^3 S3 words of length3 with a shared port after the
first step: the36 interval-product pairs give exactly36 different boundary
relations. A single whole-word product would give only6.

Partial maps require the actual composite domain, including every internal
guard. An interior state later shared with another attachment becomes a
new boundary: the old whole-interval summary may no longer suffice.
This recovers the purpose of the inherited private-path theorem rather
than asserting a new free compression of all future incidence data.

For a Fibonacci macro the recursive product still costs logarithmic
construction depth when no interior ports are shared. The price of adding
shared ports is to refine the interval partition, not to abandon the
grammar. The precise all-size design rule is therefore **compress between
actual shared ports; retain base-detecting data around actual cycles**.

## 7. Map ledger and next decisive target

| Source -> target | Map | Preserved predicate | Lost data / required sidecar |
|---|---|---|---|
|Three-face tiling cover -> scalar seven terms|count actual slots|carrier size|overlap addresses and complement phase|
|Odd-cell chart loop -> overlap owner tests|evaluate permutation at named owners|pointwise agreement|pointwise stabilizer; choose a base|
|Chart loop -> abelian character|group homomorphism|character product|commutators; retain ordered transport|
|Recoloured Fibonacci row -> length, debt, product|scalar and group summaries|visible length, total debt, endpoint map|chosen occurrence positions|
|Bijection path -> shared-boundary products|compose separately between named ports|all boundary pin queries|interior queries require new boundary splits|

The next useful target is not another numerical match between three
colours and three cells. It is to retain enough overlap owners to detect
incompatible loop defects, then choose an economical base for the remaining
automorphism group of each compatible component. These are separate tests.
They can use far fewer marks than independently marking every cell, but
must use the same actual owners and chart choices around every shared cycle.
The parent packet develops those joint odd-cell rules. No Collatz descent
claim, universal metagraph classification, or physical interpretation of
the hidden units follows from these lemmas.

## Reproduction

    python 05-knowledge/results/odd_cell_completion_20261023/fibonacci/check.py
    python -O 05-knowledge/results/odd_cell_completion_20261023/fibonacci/check.py

The frozen receipt covers all7! candidate P7 automorphisms; every H21
pair product; every triple of charts for the three-cell witness with each
declared pin set; all128 owner subsets for each of its10 subgroups;
the prefix formula through2000; all256 recolourings through row21 and
their shorter shared rows; and all216 S3 three-step boundary words.
The checker uses explicit failures, not removable assertions. It imports
the inherited zero-slot grammar to replay legality, but derives the P7
automorphism group independently by permutation enumeration. The general
theorems are established above; these finite controls test their boundaries.
