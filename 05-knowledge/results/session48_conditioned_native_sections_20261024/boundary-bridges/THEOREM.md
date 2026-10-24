# Sharp33-clock routing with two-point native boundary sections

Status: **PROVED** in the declared finite-native scope, using the explicitly
**CITED** authenticated native-menu and complete-old-core results below. The
new arithmetic claims are **FINITE-EXACT**, independently reproduced. The clock bound33 is sharp
in this declared inventory. Eight slots is an upper bound, not a claim of
sharp support. No all-order or arbitrary-host assertion is made.

## 1. Exact universe and inherited dependencies

Freeze the82 source slots, both choices per source, source652-owner host,
all indexed pins, retained source2125 blocks, full kernel lifts, actual
converse/partner operations, and normalized capacity-one resources from
`external_pin_sections_20261024/native/CORE_PAYS_BOUNDARY.md`, at repository
commit `c09b642a629b7d12b529b88c2697abaee78e73e7`. Impose the same actual equality

    path368[1] = embedded_path776[4] = embedded_path3139[7] = 4.

Thus368 is always U,776 and3139 are always V. Exactly two of the remaining
51 free helpers choose V. The old28 core sources have fixed V population k.
No other source, option, target, pin, extra owner identification or helper
population is admitted. A binary record is its actual chosen native path,
not an independently chosen charge label or normalized kernel lift.

Let Q=Z/128 x (Z/2)^4, packed as cyclic coordinate times16 plus the parity
numeral formed from original Y1,Y2,Y3,Y4. All arithmetic below uses this exact
coordinate convention. Write lambda(A) for the sum of source labels.

The following results are inherited, not new discoveries here:

- **CITED native completion:** `guarded_component_completion_20261024/native/XOR_TREE.md`
  authenticates the literal host, the complete controlled kernels and menus.
- **CITED current equivalence:** `injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md`
  authenticates the saturated paired code: at fixed integer margins, Q equality
  is the original complete aggregate nodewise dyadic current constraint.
- **CITED full-old-core connectivity:** `sixslot570_20261024/FIBRES.md` proves
  every nonempty full28 core fibre connected with at most eight changed slots
  and22 primary clock occurrences per event. This incoming Session47 result
  is explicitly not claimed anew. It applies with helpers fixed, but not
  with further core records frozen. Our constructions never make that mistake.

The inherited42-clock construction is replaced, not merely replayed. Exact
upstream source files are included in `upstream/`; archival theorem snapshots
are in `archival_markdown_snapshots.txt`. No Git or outside checkout is needed
to run the new checks.

## 2. The sheet-access principle

Fix a total charge t and helper pair H. Its core sheet is

    F_H = { C : |C|=k, lambda(C)=t-lambda(H) }.

Each nonempty sheet is connected at joint(8,22), by the cited old-core theorem.
To connect sheets H and H', it suffices to exhibit TWO actual core records
C and C' with their separately prescribed charges, such that the simultaneous
transition (C,H)->(C',H') fits the desired slot and clock budgets.

The route from an arbitrary actual endpoint first moves within its own sheet
to C, executes this actual bridge, then moves within the new sheet to the
next bridge record. A spanning graph of sheets therefore connects the full
fibre. Initial and final within-sheet paths retain and recover the entire
prescribed endpoint records. If the endpoint helper pairs coincide, the
old-core route alone fixes every helper throughout. Otherwise the final
helper assignment is exactly the prescribed target assignment; helpers are
not left at a convenient auxiliary assignment.

This is an existence-and-construction statement: finite BFS in the cited
core move graph and in the verified sheet graphs produces the route. The
algorithm never chooses one minimum at the first charge and an unrelated
minimum at the second while pretending their difference has a certified cost.

## 3. Genuine two-point sections at grades5 and6

Let R be the20 cost-two core sources and E the remaining eight cost-three/four
sources, exactly as in the inherited20+8 partition. A chart consists of fixed
selected coordinates T and three variable selected coordinates in R\T.
Its charge image I_T retains an actual full subset for every charge it contains.
Every two charges in I_T therefore have actual representatives agreeing on T,
with at most six changed cheap slots and at most12 clock occurrences.

The complete chart calculations are:

| Core grade | Fixed selected chart coordinates | Charts | Charge-image entries |
|---|---|---:|---:|
|5|two selected sources in R|190|132,658|
|5|one source in E and one source in R|160|128,392|
|6|three selected sources in R|1,140|682,678|

The grade5 cheap-only charts miss exactly1,008 ORDERED charge pairs; (2,333)
is one explicit hostile. Adding the160 charts with one fixed E source covers
ALL4,194,304 ordered pairs. The grade6 charts likewise cover all ordered pairs.

The checker forms, for every q, the exact bitset union

    union of I_T over every chart T with q in I_T.

Every row equals all Q precisely when every prescribed pair(q,q') has a
single compatible chart. This is a two-charge statement, not difference
coverage, an average, independent minima or a sampled bridge claim.
`two_charge.py` checks every943,728 actual image-entry witness and independently
reconstructs the images by enumerating full records first, then their fixed
coordinates. Both constructions agree exactly.

For every virtual grade g in6..14, fix any g-6 sources of E and use the grade6
charts on R; this padding exists because g-6<=8. The two required R charges
are translated by the SAME padding charge. For g=5 use the two grade5 chart
families. No claim of22-clock connectivity with the padding frozen is used:
the padding is fixed only during the one cross-charge bridge.

Any helper Johnson exchange changes two actual helper slots of combined
cost at most15. Choose a two-point core bridge for the two required charges.
The combined event costs at most12+15=27 and changes at most6+2=8 slots.
The full Johnson graph of two-subsets of51 helpers is connected. Hence
all core grades5..23 have joint(8,27) connectivity, including endpoint access
at(8,22). This improves the old graded-section normalization without invoking
its20-slot/42-clock intermediate events.

## 4. Core-only complementation and the low boundary grades

For k>14 use a virtual core C of grade g=28-k and sign sigma=-1; otherwise
g=k and sigma=+1. Helpers are NOT complemented. With B=the full core charge
for sigma=-1 and B=0 otherwise, the physical charge is

    B + sigma*lambda(C) + lambda(H) + lambda({776,3139}).

Accordingly each helper sheet in virtual total charge q requests exactly

    lambda(C) = sigma*(q-lambda(H)).

Complementation preserves actual changed slots and clock work. Every signed
sheet check retains this formula, both integer margins and unchanged helper
labels. Complementing the entire total charge would be wrong here.

For virtual grades2,3,4, enumerate ALL core subsets and choose the actual
minimum(selected clock weight, numeric binary mask) representative r_g(a)
for each attained core charge a. There are respectively378,3276,20475 core
states and351,1646,2048 attained charges.

For each sign and each of2048 virtual total charges, the graph vertices are
exactly the helper pairs H for which r_g(sigma*(q-lambda(H))) exists. Join
H,H' only when they differ by one helper Johnson exchange, their actual
chosen core representatives differ in at most six sources, and their combined
literal clock work is at most33. Thus every recorded edge is an actual event
with at most eight changed slots, not an abstract quotient debt.

All12,288 graphs are nonempty and connected. `sheets.cpp` checks the complete
graphs twice: recursive core enumeration plus BFS, and numeric Gosper masks
plus union-find. Both verify the same edge fingerprints and exact totals:

| Grade | Sheet vertices per sign | Accepted sheet edges per sign | Largest accepted event work |
|---:|---:|---:|---:|
|2|447,525|3,750,803|29|
|3|2,098,650|82,618,312|33|
|4|2,611,200|81,965,534|31|

These are finite exact intermediate bounds, not sharpness claims for those
grades. The sheet-access principle gives joint(8,33) for each full fibre.

For virtual grades0 and1, `low_grade.py` directly enumerates every complete
joint state(C,H), separates exact signed charges, and computes Kruskal
bottlenecks in the complete actual weighted fibre graph. There are1,275
states/972 fibres at grade0 and35,700 states/2048 fibres at grade1, per sign.
Their exact uniform bottlenecks are28 and33, respectively. Every possible
edge changes at most four or six actual slots, so support8 is automatic.
This calculation is independent of the helper-sheet BFS/union-find route.

## 5. Exact sharp33-clock lower cut

At physical core grade1, virtual charge1156 has exactly the following nine
states; fixed helpers776 and3139 remain selected in every row:

    core1457; helpers753,2135
    core1532; helpers1650,2983
    core1841; helpers2397,2677
    core1923; helpers1596,2798
    core1923; helpers1948,2616
    core2087; helpers811,1624
    core2176; helpers811,1596
    core2835; helpers862,2397
    core3360; helpers2218,2483

After adding the fixed helper charge, their physical packed total charge is398.
The state core1532/helpers1650,2983 is isolated at clock budget32; the other
eight form the other component. Its minimum crossing work is33. One crossing
is to core1457/helpers753,2135, with changed-source cost ledger

    1532:4, 1457:2, 1650:7, 2983:7, 753:7, 2135:6.

This is the ENTIRE state space of the fibre, not a selected subset or a
restricted move grammar. Hence no route at budget32 exists, regardless of
support allowance or number of events. Budget33 connects it. The independent
complete-state audit also finds a13-state opposite-sign hostile at virtual
charge700, physical grade27 and total charge265; its isolated state is the
complement of core2553 with helpers1624,2125. All states, both32-cuts and
minimum33 crossings are recorded in `low_grade.normal.json`.

## 6. The resulting theorem and native sidecars

For EVERY k in0..28 and EVERY NONEMPTY fibre of the declared externally pinned
82-slot menu with core population k, helper population four and original full
Q charge, any two prescribed full records are connected by actual simultaneous
replacement events satisfying

    at most8 changed source slots and33 primary slot-clock occurrences.

Thirty-three is the sharp uniform clock budget in this state space. No claim
that eight is the minimum possible uniform support is made. The whole packet
still has164 alternative normalized resources/656 reserved closed words,
occupies82 resources/328 closed words at each arrival, and has capacity one.
No kernel lift, opposite/converse copy, or metagraph multiplicity adds capacity.
Any unchanged outside remainder must avoid the whole reserved option reservoir
and supply the chosen baseline current signature. Clock work counts primary
source-slot occurrences; separately charged reflected/inverse copies require
their own extra ledger.

`native_check.py` reconstructs and authenticates all164 actual source paths,
normalized classes, partner/converse squares, inherited owner pins, source2125
consecutive blocks,652-owner host compatibility and the immutable external
owner4 equality. It checks the exact input closure and retained kernel
certificate, rather than inventing a new host or silently dropping a lift.

Every accepted event preserves both integer margins and Q. For the literal
paired current, each source-row change is killed by2(f_U-f_V)=0 at its modulus;
fixed population cancels the target diagonal contribution. Q equality kills
the X target sum modulo128 and Y1..Y4 sums modulo2; the checked identity
Y0=X+Y4 modulo2 supplies the omitted parity coordinate. Thus all original
physical source AND target rows, with factor-two black weight and diagonals,
are preserved. The cited saturated-code theorem supplies every genuine
dyadic level. The checker also literally replays every state of both sharp
fibres and two actual27-clock interior bridges. These examples supplement,
but do not substitute for, the universal argument and complete finite checks.

The endpoint records are exact. No resource is simultaneously reused.
Released options MAY be reoccupied later. Therefore strict temporal nonreturn
is a separate open routing restriction, not a consequence of this theorem.
Atomic accepted events need not decompose into permitted single-metagraph-edge
moves. No universal ordinary per-port phases, odd-primary guards, unchanged
integer X current, all-order native supply, global Hamiltonian observer,
occupied-atlas theorem, LRC exclusion or Collatz result is claimed.
