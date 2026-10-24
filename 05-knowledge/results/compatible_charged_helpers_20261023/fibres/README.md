# Binary boundary obstructions and occupied charged contexts

23 October 2026. **FINITE-EXACT:** complete fixed-margin/full-paired-charge
fibres in the authenticated 28-source two-target menu. **PROVED:** the
reduction to binary occupancy, the minimal-margin obstruction, complement
transport, and the stated native realization of the finite certificates.
No general all-order wiring theorem is asserted.

## 1. Inheritance and the precise question

The incoming [six-trade lattice theorem](../../balanced_native_20261023/SIX_TRADE_LATTICE.md)
proves that 27 four-slot trades and 1,917 six-slot trades generate the entire
integer charge kernel, using an index-2,048 determinant certificate. It
explicitly leaves binary-fibre connectivity open. The earlier
[charged workspace packet](../../native_packet_decomposition_20261023/native/README.md)
proved connectivity for one particular margin-eight signature and repaired
its sixteen-source endpoints with genuine compensators.

We test the missing implication in the **complete inherited supply**, not
a chosen two-path sublist: does additive generation supply a feasible route
without exceeding any binary source capacity? The concept board is integral
generation / binary boundary / actual normalized ports / charged background
/ helper-choice hypergraph. The least-used sidecar is whether a supplied
source is already occupied at target V. Neutral capacity and charged
occupancy are different coordinates.

The source is a full native packet record. Its target representation is a
28-bit vector b, with b_s=1 when the actual source-s path lands at V=2234,
and b_s=0 when it lands at U=1965. The map preserves source frames, fixed
owner clocks, actual inverse/converse resources, and the effective paired
charge

    Sigma(s)=(X(s) mod128,Y1(s) mod2,...,Y4(s) mod2)
             in Q=Z/128 x (Z/2)^4.

The inverse map is the inherited literal path menu. No source kernel is
discarded: all 28 source automorphism groups are trivial. All binary
mixtures are resource-compatible because the 56 actual paired choices
have disjoint four-word inverse/converse resource orbits.

The binary domain requires the target-margin restriction as well as pins:
raw common-clock pins allow one extra source1923 path to target1559; its
required zero occupancy excludes it. This is the previously audited
Session46 boundary, not a new claim of pin-only completeness.

At fixed margin k and full charge q, a feasible state satisfies

    sum b_s=k,   sum b_s Sigma(s)=q.

Every pair of such states is a valid simultaneous native packet replacement.
Its changed-slot cost is their Hamming distance. A route of budget six
means every consecutive replacement changes at most six slots. Integer
span alone does not ensure its intermediate b coordinates remain in {0,1}.

## 2. Complete small-margin census

Every subset at each listed margin is enumerated, with no selected-charge
filter. All 2,048 charge values occur at every listed margin.

|Margin|States|Connected fibres at budget6|Disconnected fibres|Components in the disconnected fibres|
|---|---:|---:|---:|---|
|4|20,475|1,851|197|175 have2;21 have3;1 has4|
|5|98,280|2,037|11|Each has one isolated state and one component of size28–55|
|6|376,740|2,048|0|—|
|7|1,184,040|2,048|0|—|
|8|3,108,105|2,048|0|—|

Budget eight connects **every** margin-four and margin-five fibre; because
some are disconnected at six, eight is their sharp uniform bound. Fixed
margin makes all Hamming distances even, so budget seven adds no edges.

Margins zero through three are automatically connected: any pair of
k-subsets differs in at most 2k<=6 coordinates. Thus four is the smallest
possible failing margin, even before the native census.

Global bit complementation sends (k,q) to
(28-k, sum_s Sigma(s)-q), preserving Hamming distance, native resources
and pins. Hence the census also covers margins20–24, and automatic
connectivity covers25–28. Margins9–19 have **not** been exhaustively tested
here. The complete six-trade integer lattice does not settle them.

## 3. A minimal two-state obstruction inside the full supply

The charge q=(15,3), where the second coordinate packs four binary rows,
has exactly two margin-four states:

    A={176,1425,1549,2087},
    B={1394,1923,2837,3369}.

They are disjoint and both have exact integer X-sum527. All other20
sources land at U in each endpoint. There is no third allowed state in the
**entire28-source supply**, so neither endpoint admits any nontrivial
six-slot move. Their direct eight-slot replacement is feasible. This is
minimal both in the possible failing margin and in the size of a
nontrivial disconnected fibre. It is not a claim that eight is a globally
minimal source universe among all imaginable native menus.

Their difference lies in the incoming six-trade integer lattice. Its
integer expansion therefore necessarily crosses a non-binary state if
interpreted as sequential applications starting at either endpoint. The
first failed implication is “generated as an integer vector” implies
“applicable through the capacity-one region”. No additional additive
character can fix that implication: the incoming index proof already
accounts for the whole charge kernel.

## 4. Every margin-four obstruction has an occupied-context repair

For every pair A,B lying in different six-move components of the same
margin-four fibre, consider a source h outside A union B and the new
endpoints A union {h}, B union {h}. These belong to margin five and charge
q+Sigma(h). This **changes the occupied background**. It is not a route
between the original margin-four packets.

The exact census tests all1,300 cross-component state pairs, not merely
component representatives. Every pair has a six-move connection after
such stabilization. The numbers of eligible choices h are:

|Eligible occupied source choices|Pairs|
|---|---:|
|18|11|
|19|127|
|20|1,162|

A fixed **choice bank** of four source types suffices:

    {176,1394,1477,2071}.

For each hostile pair, choose a bank member outside both endpoints whose
stabilized margin-five states connect. All3,276 possible three-source banks
have an explicit failing pair, so four is the exact minimum bank size in
this supplied universe. This is a choice among four possible occupied
contexts, not a statement that all four sources are occupied together.
The resulting route may use other source slots in the complete28 inventory.

The literal finite relation “h repairs pair(A,B)” is the useful extra
structure here. It makes helper selection a finite hitting-set problem;
the number four follows from its authenticated relation, not a numerical
analogy with other four-object constructions.

## 5. A complete ten-source native repair with exact integer current

For the minimal hostile in section3, a single added source does not repair
the eight-source core at budget six. If it initially lands at U, the full
margin-four obstruction already rules out a route. If it initially lands
at V, all20 possible nine-source restrictions still have exactly the two
endpoints. We test every supplied candidate.

Two source types suffice: add source1480 initially at V and source3000
initially at U. The ten-source supply has **exactly three** feasible
margin-five/full-charge states, all with exact integer X-sum599:

    {176,1425,1480,1549,2087}
       -> {1394,1425,2087,3000,3369}
       -> {1394,1480,1923,2837,3369}.

Each arrow changes six slots; their actual ordinary owner-clock work is
12 and15, respectively. Budget four has three isolated states, so six is
sharp for this complete three-state family. Source1480 moves from V to U
and back; source3000 moves from U to V and back. Thus the helper context
is restored at the end.

Among all380 ordered occupied/unoccupied helper pairs outside the core,
18 permit a route;12 permit one retaining the exact integer source sum.
The selected pair is a complete native witness, not an abstract codeword:
the checker reconstructs its source adjacency matrices, verifies every
Hamiltonian path and every common owner-clock pin, and authenticates all
ordinary, inverse, converse and inverse-converse normalized owners. Its20
available paired resources have80 distinct closed endpoint words; each
state uses ten distinct paired resources. An unchanged background must
avoid that full reservoir. No full-atlas rewiring is claimed.

The target owners have X(U)=X(V)=146. At a complete event, target margins
are fixed and the full source X sum is unchanged, so the target physical
integer current changes vanish. Source rows vanish because the targets
have equal X values. For a partial subgroup one must instead use

    (DeltaL X)_U = 2 sum_s w_s (X_s-X_U),
    (DeltaL X)_V = -2 sum_s w_s (X_s-X_V),

including the target-degree diagonal; no subgroup source-charge sum is
mislabelled as physical current. Equality of full paired Q sums pays every
original dyadic aggregate current at both actual target owners. The chosen
original per-slot X phase is retained. Universal per-slot phase retention
is not asserted; the inherited 4eU character separates the target ports.

## 6. Evidence, reproduction and remaining target

`census.py` enumerates all subsets at a margin and computes exact components
within every charge fibre. Its union-find scan may stop adding edges once
all already processed vertices are connected: further edges between them
cannot change their partition, and each later vertex is still tested.
This optimization never skips a possible disconnected component.

`check.py` independently uses recursive subset enumeration and intersection-
based BFS for complete margins4/5. It rechecks every stabilization pair,
every minimum-bank obstruction, all380 helper contexts, and the literal
ten-source native records. The larger positive margins6/7/8 have exhaustive
producer replays in normal and optimized Python; these are finite receipts,
not an extrapolation from the small margins.

Regenerate each census with:

    python 05-knowledge/results/compatible_charged_helpers_20261023/fibres/census.py 4 --write

and likewise for margins5,6,7,8. Then:

    python 05-knowledge/results/compatible_charged_helpers_20261023/fibres/stabilization.py --write
    python 05-knowledge/results/compatible_charged_helpers_20261023/fibres/check.py
    python -O 05-knowledge/results/compatible_charged_helpers_20261023/fibres/check.py

Omit `--write` to compare regenerated censuses and stabilization data with
their frozen receipts. `reproduce.py` runs all normal/optimized replays.
Input text hashes normalize line endings; compressed atlas hashes retain
raw bytes. No discovery cache is needed by the checks.

The open next question is a feasible-context supply theorem: which actual
occupied backgrounds force enough safe small-trade entry? The new boundary
shows why integer kernel generation cannot be that theorem. The helper
relation and its minimum bank provide a finite constructive substitute
for the margin-four cases, while native all-size wiring and the other
major conjectures remain open.
