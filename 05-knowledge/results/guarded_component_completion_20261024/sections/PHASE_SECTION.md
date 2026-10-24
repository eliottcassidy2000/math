# A fixed owner phase gives an eight-slot global route

Status: **PROVED finite-native connectivity theorem**, using the exhaustive
boundary section in `phase_sections_receipt.json` and the actual host in
`../native/XOR_TREE.md`. This sharpens that packet's valid ten-slot/61-clock
controller bound without changing any of its records or guards.

## Exact state space

Use all28 old core cells and all54 authenticated helper cells, each with its
fixed full U/V record and original pins. The only additional shared-owner
constraints are the actual XOR-tree equalities among helpers368,776,3139.
All old core cells and all other helper cells have private owners in the
652-owner host. Retain fixed core V count k, helper V count4, and fixed total
Q-charge q. Actual inverse/converse resources, capacities, clocks and local
partner matrices are those already authenticated. An optional immutable
background avoids the entire reserved alternative-word inventory.

This is the complete binary-menu state space under these declared factors,
not a search through one chosen path subset. Its helper branch with368=V
and776=3139=U is called the hub phase. It literally pins the shared owner
to3. The complete51-free-cell occupancy-three section covers every residual
charge in this phase. The other phase covers only972 charges; no arbitrary
external boundary value is assumed harmless.

## Connectivity theorem

For every k in {0,...,28} and every q in Q, any two feasible states in this
space are connected by simultaneous actual full-record replacements that
change at most **eight slots and51 owner-clock occurrences per event**.
Prescribed feasible initial and final helper records are retained. At most
min(k,28-k)+2 events are needed.

Proof. Let A,B be the initial and final core V sets.

1. Keep the actual initial core fixed. If its helper state is already in
   the hub phase, retain it. Otherwise use its source-specific nearest
   same-charge hub state. The complete admission certificate enumerates
   all1275 off-hub states and all20825 hub states: every off-hub state has
   an actual same-charge hub witness of clock cost at most49. Both states
   have occupancy4 and satisfy the actual tree equalities, changing at
   most8 helper slots. The worst minimum49 is attained by exactly two
   off-hub states; the generic53 diameter bound was not needed.
2. Move A to B by a Johnson path of |A\B| exchanges, replacing one actual
   core V cell by one U cell at each step. At each resulting core state,
   choose a hub-phase section at q minus its actual core charge. At the
   last core arrival, choose the nearest hub state to the prescribed final
   helper state (or that final state itself if it is already in the hub).
   The tree
   controls stay fixed; only three free helper V cells can vary, hence at
   most6 helper slots and43 helper clock occurrences change. Each core
   exchange changes2 slots, costing at most8 clocks. Thus these events
   change at most8 slots and cost at most51 clocks.
3. At the final core, restore the prescribed final helper records. Symmetry
   of actual Hamming clock cost gives the same bound8 slots/49 clocks.

If A=B, no Johnson exchange is needed. Same-phase helper endpoints differ
in at most43 clocks in the hub phase or29 in the other phase, so one direct
reset suffices. For different phases, take the nearest hub state to the
off-hub endpoint, then the prescribed hub endpoint, using at most two
events of costs49 and43. This treats the zero-exchange case without adding
an uncounted intermediate event.

The hub section exists at every core state because all core cells are
private and the section covers Q at the fixed shared-owner value3. Every
arrival keeps both exact target margins and total Q charge. The native
paired-code theorem therefore preserves all original dyadic aggregate
currents at every physical row, with the correct target-degree terms.
Resource injection and all local owner predicates follow from the actual
full-menu inventory and host theorem. This proves the route, with
|A\B|<=min(k,28-k). Zero-change events may be omitted. QED.

This is an upper bound, not a minimal budget claim. Events are coupled
replacements in the selected native subnetwork, not asserted single native
metagraph microedges. Exact integer-X neutrality, odd-primary characters,
every ordinary per-port phase, and a Hamiltonian observer of the entire
652-owner host are not claimed.

## Why boundary freedom matters

An extra external pin forcing the shared owner to4 would exclude the hub
phase. Because the remaining phase misses1076 charges, the proof would then
fail at its first admission. Similarly, the twelve-control construction in
`../native/NATIVE_COMPONENTS.md` has control phases forced by core owners;
one may not fix all those phases freely. It uses its genuinely shielded40
reservoir instead. These are different compatible state spaces, not two
gauges of the same unmarked charge table.

The mechanism is a complete section at one *actual fixed boundary value*.
Completeness at the union of several boundary values would not by itself
give the smaller route: the controls might have to change along every core
exchange. The moving-block source2125 repairs precisely the missing value
of the fixed-boundary section.

## Exact controls

`phase_sections.py` enumerates all fixed-phase subsets in both old and new
supplies, provides2048 literal hub witnesses, and derives the43-clock bound
from the actual path-change costs. `admission.py` exhausts every one of the
1275 actual off-hub configurations, retaining a minimum-cost witness and
the number of minima in its complete same-charge target fibre. Its two
hardest states certify the exact admission bound49 for this fixed phase.
It additionally checks365 same-charge off-hub pairs, all having a common
hub at peak cost49; this extra control is not needed for the theorem.
`phase_route.py` realizes a route whose
two endpoints use the opposite phase, checking each arrival's literal path
records, shared-owner equations, all normalized source/inverse/converse
resources, target margin, Q-charge and original physical generator currents.
Both normal and optimized Python runs reproduce the frozen receipts.
