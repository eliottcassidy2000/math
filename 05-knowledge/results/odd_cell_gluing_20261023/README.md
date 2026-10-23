# Compatible odd-cell attachments, zero-slot Fibonacci rows, and joint charts

23 October 2026. **PROVED:** the scoped partial-bijection gluing and
chart-conditioned transport rules, zero-slot row family, and prime-power
affine clock identities. **FINITE-EXACT:** the complete ordered two-mark tests on
P31 and D(D(P7)), including cross-target chart changes. **OPEN:** uniformly
small all-order tournament observers, unrestricted metagraph reconstruction,
black-core connectivity, and universal Collatz entry.

The direct result separates stalled refinement from real ambiguity.
Of the 426 nondiscrete two-mark partitions of D(D(P7)), 168 have no
nonidentity mark-fixing automorphism. Keeping the same chart choice through
several cell attachments detects this rigidity. A small cell-relative
four-vertex flag reduces 126 of the 132 seven-cell cases to cyclic triples;
six cases genuinely retain a group of order21 and are not forced into that
model.

The user's zero-length Fibonacci slot supplies a second concrete repair:
it changes adjacency of atom indices, permitting A+B=2. The previously
blocked B at diagonal35 now extends to row36. The precise result is a
whole family of compatible rows, not a unique continuation chosen by a
finite prefix. See [ZERO-SLOT.md](ZERO-SLOT.md).

**Late incoming reconciliation.** Commits14df2c7f01 andd45865c332 arrived
while this packet was being tested. The [owned endpoint packet](../endpoint_overlap_ports_20261023/README.md)
independently closes the same order31 census and minimum extra-owner counts.
Our `cross_check.py` compares every marked case and every literal initial
partition. Its subset-equality selector and our binary cyclic-triple flag
produce the SAME repaired partition in all132 seven-cell cases. Credit
these counts as convergent results, not two separate discoveries. The
additional contributions here are the partial-bijection state-graph rule,
chart-conditioned residual-potential test, alternative bit observer,
zero-slot compatible-row family, and ordered Fibonacci clocks.

Final incoming commit295d14b301 sharpens the next target further: the
[prime-cell displacement necklace](../session40_marked_ports_square_sum_20261023/PRIME_CELL_BURNSIDE_FORMULAS.md)
is already a complete all-size carrier for its matching-coupled cyclic
cells. Its scalar contractions fail from prime7 onward. Thus the next
step is to glue its retained rotation choices across several attachments,
or handle more general regular blocks, rather than seek another scalar
replacement for that complete word.

## Inheritance, current truth, and concept board

Started from clean a1f4dcda5c, fetched and fast-forwarded through69648c025d.
Read the incoming rooted-five/marked-phase frontier, Session39 affine
connections and girth hostiles, and Session40 shared-incidence/Hopf packet.
These strengthen the earlier session; the affine connection theorem below
is inherited, not re-proved as a new discovery.

- Closest proved mechanism: [Session39 affine connections](../session39_signed_phases_and_prime_cylinders_20261023/marked-cells/AFFINE_CONNECTIONS.md),
  complete for named equitable singleton/triple partitions; and
  [Session38 residual-cycle potentials](../session38_incidence_ports_electrical_20261023/theory/PROJECTED_TRANSPORT_UNIQUENESS.md),
  exact for a specified integer transport relaxation.
- Canonical hostile: [the order31 mark tests](../incidence_wiring_20261023/MARKED-CELLS.md).
  One mark leaves every DRT at1+15+15; two marks can leave seven-cells.
- Corrected near miss: [girth-phase hostiles](../session39_signed_phases_and_prime_cylinders_20261023/marked-cells/GIRTH_PHASE_HOSTILES.md)
  refute fixed local incidence completeness. No matching owner section is
  required for a tournament to exist. The chart section constructed here
  is instead an isomorphism between two specified attachment systems.
- Least-used sidecars: the projected stabilizer of a complete two-cell
  block, and the zero slot's effect on nonconsecutive atom indices.
- Recovered older niche: [three-unit rows](../reset_20260926_colours.md),
  [the historical supplied stripe](../crossroads_crossing_20260926_colour.md),
  and [the four-state Fibonacci charge](../duck_zeckendorf_20260925.md).
  The latter charge is not the user's three-colour word.

| Lane / concept | Operation and preserved target | Lost coordinate / next test |
|---|---|---|
| Anchor: several odd cells | join actual local chart relations | separately chosen chart witnesses may conflict |
| Port margins | residual cycles after fixing joint charts | a union of allowed cells can introduce phantom switches |
| Seven-cells | cyclic-triple flag in a named partner cell | unmarked motif totals forget the partner |
| Niche: zero-slot rows | indexed nonconsecutive decompositions | deleting the empty slot changes admissibility |
| Wildcard: Fibonacci transport | ordered macro products with exact lengths | colour totals and periodic phase forget order and clock |

The method card used is **Retain the product and shared labels before
 taking scalar totals**, in `00-navigation/META-PATTERNS.md`. This packet
supplies a new application, not a duplicate generic method card.

## 1. An all-size exact rule for functional attachment constraints

Let G be a connected finite simple graph with named vertices i. Each i
owns a finite chart domain C_i. For each oriented edge i->j retain an
actual partial bijection R_ij:C_i --> C_j; reverse transport is its inverse.
Empty domains and empty relations are allowed. A compatible section is
one s_i in each C_i satisfying R_ij(s_i)=s_j on EVERY edge, with every
required map defined. Retained source/port marks can restrict C_i.

**Theorem (partial-chart gluing).** Choose any root r and spanning tree.
For each s_r in C_r, propagate along the tree, rejecting undefined maps.
Then check all remaining edges, including their domains. This enumerates
all compatible sections exactly once. In particular there are at most
min_i |C_i| sections; the bound does not multiply over cells.

Proof. A section's tree values are forced by its root value. Propagation
recovers them and no genuine section is rejected. Conversely every
accepted assignment satisfies all edges, hence is a section. Using any
vertex as root proves the minimum bound. This is finite-domain constraint
propagation, with no claim of general novelty or uniformly bounded domains.

There is an equivalent certificate avoiding a choice of spanning tree.
Make the **state graph** with vertices(i,s), s in C_i, and edges joining
(i,s) to(j,R_ij(s)) whenever defined. A section is exactly a connected
component satisfying BOTH:

1. projection to G has one state above every base vertex;
2. each state has an edge above every incident base edge.

Proof. A section contains all neighbors of its states: functionality and
injectivity prevent a second neighbor above any base edge. Its states
therefore form a whole connected component. Conversely the two conditions
supply all the required equations. Coverage without condition2 is false:
a lifted spanning path can cover all three vertices while a triangle's
closing relation is undefined. The checker exhausts all34^3=39,304 systems
of partial bijections on three states over a triangle, comparing the
component decoder against all27 assignments.

After the relations are supplied, component testing is linear in the
explicit lifted graph size. Enumerating the chart domains or constructing
the relations can be expensive. This is not a polynomial-time algorithm
for arbitrary tournament isomorphism. For disconnected G, apply the rule
per component; any extra constraints between components must still be
checked jointly.

### Actual tournament blocks provide partial bijections

Let T,T' have corresponding NAMED cells D_i,D'_i and retained marks.
Take C_i to be all internal tournament isomorphisms D_i->D'_i preserving
those marks. For a pair i,j, define compatibility by

    T(u,v) = T'(f_i(u),f_j(v)),  u in D_i, v in D_j.       (1)

If the target block has pairwise distinct columns, fixing f_i determines
at most one f_j: each source column, reindexed by f_i, must match one
specific target column. Distinct rows give the reverse uniqueness.
Thus every row-and-column-distinguishing block gives an actual partial
bijection between chart domains. This statement allows unequal cell sizes.

If a chosen connected graph of such blocks spans the nontrivial cells,
root propagation gives at most |C_r| full candidate maps. Check every
remaining block, singleton attachment, and retained path flag on those
candidates. The surviving maps are exactly the cell/mark-preserving
isomorphisms. Uniform blocks impose no transport and cannot be treated as
functional: two C3 cells joined uniformly allow all3x3 chart pairs.
Other nondistinguishing blocks must remain explicit additional constraints.

This is an all-size **verification/decoding rule for supplied blocks**.
It does not manufacture missing blocks from their scalar signatures and
does not identify every relaxed port table with a realizable tournament.

### Port transport remains literal

An accepted f retains cell names and the selected source marks. If P is
an actual Hamiltonian path, transport that same ordered path to f(P).
For the path-complement operation K_P, entrywise arc inspection gives

    f(K_P(T)) = K_(f(P))(f(T)).                            (2)

Path arcs are kept on both sides and exactly the other arcs are reversed.
This is the incoming [path-flag transport](../metagraph_flag_transport_20261023/README.md)
mechanism. It preserves endpoints and their actual incidence. Independently
renaming an unresolved cell for each incident block does not give an f
and cannot be used in (2). The census below tests all ordered pairs of
marks, without asserting that each pair is realized as Hamiltonian endpoints.

## 2. Complete order31 audit: ambiguity is smaller than the cells

Use the earlier explicit doubling construction, reproduced by import from
`incidence_wiring_20261023/marked_cells.py`. Vertex and cell gauges and the
source adjacency hashes are frozen in cells.json. For each of the930
ordered distinct mark pairs, run outgoing-neighbor refinement with the
marks distinguished, retain every final cell, and solve (1) on the exact
blocks. Every nondiscrete case here has a complete informative block graph.
Every accepted full map is checked against all31^2 source arcs and both marks.

P31 becomes discrete for all930 pairs, as inherited. For D(D(P7)):

| Largest stable cell | Exact mark-fixing automorphism group size | Ordered mark pairs |
|---:|---:|---:|
|1|1|504|
|3|1|126|
|3|3|168|
|7|1|42|
|7|3|84|
|7|21|6|

The rigid total is672, of which168 were hidden by ordinary refinement.
The seven-cell blocks do not all have the same internal isotype: some
have only three internal automorphisms. The worst root chart domain has
size21, not21^4. The triple cases need only three root candidates, not3^8.
All426 nondiscrete partitions were also tested against independently
changed local target charts; the accepted maps were exactly the transported
source maps. No source/target chart identity was silently reused.

There are two independent controls beyond propagation: complete product
searches of21^4=194,481 and3^8=6,561 local chart tuples for the declared
seven-cell and triple examples, and comparison at EVERY mark pair with
explicit automorphisms lifted from

    x -> a*x+b on F7,  a in{1,2,4}, b in F7.              (3)

These21 permutations preserve P7 and its two doubles, entrywise.
They give exactly the maps found by the decoder at every pair.

**FINITE-EXACT consequence:** the full automorphism group of this D(D(P7))
has order21. Here is a short completeness proof from the certified counts.
A tournament automorphism group has odd order: an involution swapping two
vertices would reverse their arc, and inverse-pairing in any even-order
finite group supplies a nonidentity involution. Its action is free on the
672 ordered pairs whose stabilizers are trivial, so its order divides
672=32*21. It contains the explicit group(3), so its odd order is exactly21.
This argument concerns this tournament, not all doubled tournaments.

### A small flag repairs most seven-cells

For ordered seven-cells D,E, their block degree in the audited bank is3 or4.
Keep the three selected neighbors when the degree is3; when it is4 keep
the three selected nonneighbors. For each u in D define

    b_(D,E)(u) = 1 iff those three vertices induce a directed C3.    (4)

This is a cell-relative rooted four-vertex flag, computable directly from
the block and E's internal arcs. It is covariant under the same joint
cell isomorphism. It is not an unmarked four-census.

For each D, choose the first partner E in retained cell order for which
(4) is nonconstant; keep the cell identity and this bit, then refine.
This reduces126 of132 seven-cell cases to seven singletons plus eight
cyclic triples. The remaining six pairs, precisely pairs among the three
connector vertices7,15,23 in this source gauge, have constant flags for
every partner and retain their group of order21.

Within the126 repaired cases,84 have balanced signs and a remaining C3
rotation;42 have12 negative fundamental triangle signs and no nonidentity
rotation. The original294 triple cases likewise split into168 balanced
and126 unbalanced cases. These are actual shared-block tests, not inferred
from cell sizes. The flag is small and sufficient on this bank; no global
minimality claim is made.

### The exact inherited triple rule closes the joint rotation step

Normalize every nonuniform3x3 block to its selected matching, and use
cyclic charts x in F3. It has form x->epsilon*x+delta. In a self-isomorphism,
cell rotations a_i satisfy

    a_j=epsilon_ij*a_i.                                 (5)

A negative sign cycle forces a_r=-a_r, hence a_r=0 in F3 and every
rotation vanishes. A balanced connected support admits exactly three
rotations. This is Session39's all-size theorem; the new work applies it
to the complete order31 bank and to the seven-cell repair.

For comparing DIFFERENT attachment systems, retain the phase differences
as well: a_j-epsilon_ij*a_i=delta'_ij-delta_ij. The affine part cannot be
dropped. A three-cell hostile in the checker gives three compatible chart
pairs on every edge and three tree survivors, but no globally compatible
map: the closing cycle rejects all three. Pairwise feasibility alone is
not enough.

### A source-owned base, rather than independently fixing every cell

Suppose the functional block graph is connected and G is the complete
cell/mark-preserving automorphism group. Restriction of G to any root cell
is faithful: identity there propagates to identity everywhere. If G is
nontrivial, choose a moved vertex in that root cell and retain its actual
owner. Its stabilizer has index at least3, because G is an odd-order
tournament automorphism group. Repeating needs at most floor(log_3|G|)
additional root-cell owners before the group is trivial. At each step the
new stabilizer remains faithful on that cell. This is an all-size bound
relative to the supplied G, not a uniform numerical bound on all graphs.

For the declared order31 source the exact extra-root-owner counts are
672 cases needing zero,252 needing one, and6 needing two. The checker
orders candidate owners by a retained actual Hamiltonian path, chooses
the first moved root-cell vertex, and compares its count with exhaustive
root-cell base searches. A group of order21 acting on its seven-cell
has no free point, so the two-owner lower bound is real. These counts
independently agree with the incoming endpoint-overlap work; its additional
932 distinct actual-path controls are inherited, not claimed as rerun here.
Pointwise naming costs these owners; classification up to the remaining
symmetry is a different target and need not pay the same cost.

### Change the represented question when the old boundary becomes empty

All420 authentic eight-triple systems reached here have ZERO parallel
vertex-owner sections, independently confirming the incoming result.
Nevertheless each has at least the identity chart-isomorphism to itself;
252 have three such maps and168 have one. These are different predicates:

- a vertex-owner section asks for x_j=pi_ij(x_i) on every selected matching;
- a chart-isomorphism asks whether one whole attachment system can be
  relabelled into another while preserving every owned block.

Thus the 64-state parallel-boundary observer sends all420 whole systems
to its empty state and cannot reconstruct them. The relative chart-domain
construction retains an identity even when the old section problem is
inconsistent. This is a precise reroute around the obstruction; it does
not declare the empty systems invalid tournaments, and does not claim
that comparing against a complete target inventory is free.

## 3. Residual potentials must be conditioned on compatible charts

The state graph above is not the residual flow graph. They have different
objects and different conservation laws. There is nevertheless an exact
composition of the two certificates.

Suppose a finite set S of compatible shared-chart choices has been
certified. For each s, let T_s be a bounded nonnegative integer transport
model with specified allowed cells, row/column margins and reflection
handling. Discard a branch only with an infeasibility proof; otherwise
retain a feasible table x_s. Fix the same linear observed wiring map P.
Then P is constant on the union of the remaining T_s iff:

1. in every branch, each directed residual cycle at x_s has zero P-cost;
2. the baseline vectors P(x_s) agree across branches.                    (6)

Proof. The inherited residual-cycle theorem makes condition1 equivalent
to constancy on each branch. Condition2 is exactly constancy between those
branches. Both directions follow without averaging or assuming every
relaxed table has a native tournament realization. Empty total union does
not authenticate an actual source; its existence must be retained.

A small hostile shows why one must not first union the allowed cells.
Take two2x2 permutation blocks. A shared chart bit permits either both
identity matchings or both swapped matchings. Let P be the number of
off-diagonal entries in block1 minus that in block2. Both legitimate
branches have P=0. Forgetting that the SAME bit controls both blocks admits
independent swaps and creates P=+2 and P=-2. Each is a real table in the
larger relaxation but a phantom of the shared-chart model. Conversely,
checking uniqueness in each branch without comparing baselines can miss
actual differences between branches.

The exact rule therefore uses both certificates: chart compatibility
first, residual potentials within each compatible branch, and common
projected baselines. This does not prove a small bound on the number of
branches for arbitrary nondistinguishing blocks.

## 4. Source-preserving refinement: private paths versus shared ports

A path whose internal vertices have no other incidences or observations
may be replaced by the composition of its actual partial bijections.
The composite is a partial bijection; its domain contains precisely the
starting states for which every intermediate map is defined. Each accepted
endpoint pair has one and only one internal lift. Thus this contraction
preserves all boundary sections, and expanding it preserves the source.
The checker compares all34^2 pairs of partial bijections on three states.

If an internal vertex is a shared port, suppressing it independently in
two paths is invalid. Retain its state or the joint relation between all
exposed endpoints. The same principle governs recursive Fibonacci pieces:
compress the ordered transport and its exact length, not separate colour
counts or unrelated copies of a boundary state. The new zero-slot row
family retains occurrence addresses and compatible choices across all rows.

### Incoming aligned-pair theorem sharpens the target

The [new aligned refinement proof](../aligned_port_refinement_20261023/PROOF.md)
shows that an actual source A and its ALIGNED partner B jointly refine to
singletons within floor(n/2) rounds: A+B-(J-I) is the signed Hamiltonian
path matrix. Their common nonnegative doubly stochastic centralizer is
also trivial. This already proves full-data rigidity; it does not find
an unknown B or its alignment. Our new gluing certificate must therefore
be evaluated on compressed candidate alignment records, not sold as a
stronger theorem about an already supplied actual pair. The incoming
linear-mark bounded-core rule is another all-size fallback with its own
budget; a fixed tiny observer remains open.

## 5. What remains open

This session closes the joint-chart classification for the declared
order31 controls and supplies an exact all-size rule whenever enough
actual blocks give functional chart relations. It does not derive those
blocks from a fixed small motif bank, supply a full order31 Hamiltonian-port
census, or improve unrestricted A000568 complexity.

The next sharp target is **nonfunctional block relations**: retain their
projected stabilizer cosets and determine when a small shared separator
makes the remaining relations functional. A related optimization is to
choose one source-owned mark that shrinks all surviving chart domains,
rather than independently fixing a vertex in every odd cell. The six
order31 order21 cases are the required positive-symmetry control, not
examples to be falsely declared rigid.

The new prime-cell necklace provides a concrete positive control: its
ordered differences reconstruct a matching up to the two cell rotations,
and the same rotation must be used on every block incident to that cell.
Do not replace these permutations by affine maps for primes at least5;
the incoming exact orbit count grows as p!/p^2. By contrast, our affine
Fibonacci clock is an explicitly restricted group model. Its eight-step
translation return lifts to exact periods8p^r, and at modulus7^k the
whole pair map has order16*7^k. These statements compress ordered group
products; they do not classify the larger matching carrier.

Two incoming hostiles retain the missing ownership coordinate. At order5,
the SAME complete source and actual endpoints allow complement targets
with Hamiltonian counts15 and11, depending on the internal path order
([source-only obstruction](../session40_marked_ports_square_sum_20261023/ENDPOINTS_DO_NOT_DETERMINE_WIRING.md)).
For the separately defined square-sum birth colouring, actual endpoints
and the complete ordered colour word first fail to determine a path at36
([finite census](../session40_marked_ports_square_sum_20261023/square-sum/COLOR_FIRST_FAILURE.md)).
The zero-slot word therefore certifies its specified recursive row family,
not arbitrary path ownership. The shared lesson is exact: retain the
position-to-owner map required by the target predicate.

## Reproduction and audit

    python3 -B 05-knowledge/results/odd_cell_gluing_20261023/check_cells.py
    python3 -O -B 05-knowledge/results/odd_cell_gluing_20261023/check_cells.py
    python3 -B 05-knowledge/results/odd_cell_gluing_20261023/check_gluing.py
    python3 -O -B 05-knowledge/results/odd_cell_gluing_20261023/check_gluing.py
    python3 -B 05-knowledge/results/odd_cell_gluing_20261023/check_fibonacci.py
    python3 -O -B 05-knowledge/results/odd_cell_gluing_20261023/check_fibonacci.py
    python3 -B 05-knowledge/results/odd_cell_gluing_20261023/cross_check.py
    python3 -O -B 05-knowledge/results/odd_cell_gluing_20261023/cross_check.py
    python3 -B 05-knowledge/results/odd_cell_gluing_20261023/check_clock_lifts.py
    python3 -O -B 05-knowledge/results/odd_cell_gluing_20261023/check_clock_lifts.py

Explicit failure checks survive -O. `--write` intentionally regenerates
receipts. [AUDIT.md](AUDIT.md) records proof boundaries, rejected trial
claims, independent computational paths and the final hash manifest.
