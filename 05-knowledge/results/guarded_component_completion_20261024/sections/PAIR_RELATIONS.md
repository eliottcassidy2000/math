# Pairwise coverage, a native three-cell obstruction, and its relative-block repair

Status: **FINITE-EXACT** for the stated order-eight supply; **PROVED** for the general implications of the finite certificates. No all-size tournament classification is asserted.

## 1. Exact carrier and inheritance

Use the authenticated53 helper sources in
`../../compatible_charged_helpers_20261023/supply/selected_supply.json`, with
the original full path records, indexed owner pins, inverse/converse resource
inventory, and fixed targets U=1965, V=2234. A bit1 selects V. The charge group
is Q=Z/128 x F2^4 and the source label is

    lambda(s) = (X_s mod128, Y1_s mod2, ..., Y4_s mod2).

All guards are the actual original native electrical guards. The complete
paired-code theorem in
`../../injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md`
converts fixed target margin plus equal total Q-charge into preservation of
all original dyadic **aggregate** currents. This does not preserve every
ordinary per-port harmonic phase. The chosen original X phase is retained
per slot because X_U=X_V=146. Exact integer X sums are not imposed here.

Closest mechanism: Session47's12 disjoint quadruple witnesses per charge.
Canonical hostile: one V together with one U can remove the last witness for
32 directed pairs. Corrected near miss: a full charge section is not
automatically closed under conjunction of independently safe predicates.
The least-used sidecar below is a directed consecutive owner block whose
position may move, rather than an indexed owner value.

## 2. Every single two-bit relation of size at least two is safe

For every unordered pair of distinct helpers i,j and every q in Q, at least
three of the four patterns00,01,10,11 occur in an occupancy-four state of
charge q. Therefore **every relation R subset {0,1}^2 with |R|>=2 leaves a
nonempty charge-q fibre**, simultaneously for every q and every helper pair.

Proof of implication: a set of at least three patterns cannot be disjoint
from a set of at least two patterns in a four-element space. This handles
equality, XOR, both implication directions, OR, NAND, and one forced bit.

The finite statement is exhaustive over all292825 quadruples,1378 pairs,
and2048 charges. With pattern bits00=1,01=2,10=4,11=8, the availability
histogram is:

| Available mask | Number of pair/charge cases |
|---|---:|
|15|1,340,057|
|7|1,482,055|
|13|16|
|11|16|

Singleton00 is safe for every pair; each oriented mixed singleton is safe
for1362 of1378 pairs; singleton11 is safe for no pair. The latter does not
contradict the theorem, which requires at least two allowed patterns.
`pair_sections.py` uses witness bitsets; `check.py` independently uses
single-occurrence and pair-occurrence counts, reconstructing the four
nonnegative pattern counts by inclusion-exclusion at every charge.

## 3. Two individually safe factors can fail together

The satisfiable three-node XOR tree

    x368 != x776,       x368 != x3139

has exactly two local states: (1,0,0) and (0,1,1). On the full53/k4 bank its
image is Q minus {1509}, where q is encoded as16 times its Z/128 coordinate
plus the four parity bits. This is minimal in number of such pair factors:
one factor is always safe by Section2. No claim of a complete catalogue of
three-node trees, or of minimality among arbitrary larger-arity predicates,
is made.

The two failed branches explain the obstruction exactly:

* If368 is V, the remaining triple must have charge305 and avoid776,3139.
  The complete charge305 triple list is
  {776,1182,2638}, {1391,1948,3139}, {1400,3111,3139},
  {2439,2728,3139}. Every triple violates that avoidance.
* If368 is U, the two leaves776,3139 are V. The remaining pair must have
  charge239. No such pair occurs in the53 bank. In the complete72 actual
  U/V-source candidate alphabet the only such pair is {1596,2125}.

The companion native packet `../native/XOR_TREE.md` authenticates this very
tree as an actual common-owner construction, using equality of three actual
indexed owner coordinates and checking both source and local partner arcs.
That additional host theorem is necessary: a Boolean relation alone would
not certify an actual shared-owner tournament embedding.

## 4. A new guard grammar supplies the missing helper

Source2125 has401 Hamiltonian paths and exactly two paths whose targets lie
in {U,V}, one per target:

    U: 5,6,2,0,7,3,4,1
    V: 7,1,0,3,4,5,6,2.

It was excluded from the53 bank solely because the old selection rule
required a common indexed owner clock. These two paths have none. Their
relative clock permutation is the single cycle
0,5,3,2,7,1,6,4, so even a nonempty proper unordered clock-set pin cannot be
retained. This is a genuine failure of that pin grammar, not missing native
paths or missing electrical charge.

Both paths retain the directed consecutive blocks5->6->2 and3->4. Requiring
these **moving owner blocks** gives19 raw paths; imposing the inherited
U/V target-margin filter leaves precisely the displayed two. No old53 pin
has been weakened or removed. The new slot has a different, explicitly
declared full-record guard. Its actual owner-clock change costs8.

`native_2125.py` authenticates the literal source matrix, all401 paths,
automorphism order1, unique full lift of each selected normalized resource,
pure-black capacity2 at each target, literal K partner, converse paths with
reversed blocks, and all source/inverse/converse words. Its8 reserved words
are disjoint from the previous648, making656 in the complete reservation.
The inverse target words give the target-port injection without an arbitrary
target isomorphism. All mixed choices remain injective because the entire
reserved menu union is disjoint. The construction is a selected native
subnetwork with an untouched disjoint background, not a replacement of the
full order-eight atlas.

The enlarged54/k4 bank with the same XOR tree covers every q in Q. Previously
represented charges keep2125 at U. Charge1509 now has exactly two witnesses:

    {368,676,1770,2125},   {776,1596,3139,2125}.

The second witness is forced by the unique full72 pair239; the first shows
that the newly supplied helper can also repair the other tree phase. This
repair therefore preserves occupancy4. The companion host may embed2125
privately, adding no shared-owner conflicts.

A stronger boundary statement is available. With center368 fixed V and
both leaves fixed U, the original bank covered2047 charges, missing only
1509; the enlarged bank covers all2048. Its51 remaining private helpers
therefore form a complete occupancy-three section after subtracting the
fixed control charge. In the companion embedding this keeps the common
owner value literally3. The opposite fixed phase covers only946 charges
before and972 after: the theorem does not grant arbitrary externally forced
phases. `phase_sections.py` records a full2048-entry fixed-phase section.

## 5. Actual reset cost and useful large-component survivors

Two feasible constrained54/k4 states differ in at most8 helper slots.
Within the center-V phase their control bits agree and at most6 free slots
change; within center-U at most4 free slots change. Across phases all three
controls change and at most5 free slots change. Their actual costs are
6,5,6 for368,776,3139; the free list begins8,7,7,7,7,7,... . Therefore the
three actual clock-occurrence bounds are53,43,29, and the uniform bound is
**53 helper clock occurrences**. This retains repeated occurrences of shared
owners; it does not replace clock cost by number of independent variables.
A simultaneous core move adds its actual support and clock cost. For a
core two-slot move of clock cost at most8 the totals are at most10 slots and
61 clock occurrences. This is a source-specific reset section, not a claim
that every predicate-preserving one-bit move is legal.

Once the complete center-V boundary section is chosen, all three controls
remain fixed. Subsequent resets need at most6 free helper changes and43
clock occurrences. Initial or prescribed endpoints in the opposite phase
still require the cross-phase bound; they cannot silently be replaced by
convenient same-phase records.

Component size alone is not the obstruction. On the original53 bank an
arbitrary outward implication star h->j is safe by fixing h=U; an arbitrary
inward implication star j->h is safe by fixing h=V. Each such forced-bit
fibre covers Q by Section2, regardless of the number of leaves. More
generally any collection of implications with at most11 distinct
antecedents is safe: set all those antecedents to U, and choose one of the
12 disjoint quadruple witnesses avoiding them. This is a theorem about
actual record predicates; embedding a chosen implication system into a
common-owner host still requires the partner-gluing and resource checks.

## 6. Reproduction and map ledger

Run `pair_sections.py`, `native_2125.py`, `phase_sections.py`, and `check.py` from this directory
(also with `python -O`). Receipts are frozen; `--write` regenerates them.
The pair producer and independent count checker use different algorithms.
The native helper enumerates8! permutations directly, using the original
atlas only for class ownership and conductance.

Source: full actual path records and their common-owner compatibility.
Target: bit predicates plus Q-charge section. Map: choose U/V target and
sum authenticated source labels. Preserved: fixed margins, all original
dyadic aggregate current equations, chosen X per-slot phase, supplied pins,
literal normalized resources, inverse/converse and actual arrival clocks.
Destroyed by the bare bit/charge quotient: the full owner order, kernels,
local partner arcs and capacity conflicts. Needed sidecar: the retained
full record tables and the actual host/partner agreement proof. The cheapest
decisive test was a three-node XOR tree, followed by the complete401-path
source2125 census; neither required an invented electrical formula.
