# Twelve moving owner pins paid by a protected native charge reservoir

Status: **PROVED construction and lifting mechanism; FINITE-EXACT native
record authentication.** Two actual shared-owner components coexist with the
original28-core packet and53 new helpers. Every declared core margin and every
paired charge has compatible coverage. This is a selected order-eight packet
family, not an all-order native supply theorem.

## 1. Inheritance and the missing coordinate

The closest proved mechanism is the state-dependent section theorem in
`compatible_charged_helpers_20261023/controller/COMPLETE_CHARGE_BANK.md` §5.
The closest actual realization is `../supply/shared_owner.md` in that same
packet: one core/helper indexed-owner equation, with source and partner host
gluing. The current native data is its immutable
`supply/selected_supply.json`, including complete pin domains, full kernels,
black coefficient two and disjoint normalized alternative resources.

The canonical hostile is the occupancy-four bank with two forced-V helpers:
only C(51,2)=1275 choices remain, fewer than2048 charges. Thus merely adding
independent pin constraints does not inherit the unrestricted charge section.
The newer Session47 mixed-pair census leaves32 bad directed one-V/one-U
restrictions. Our construction uses a different mechanism: preserve a
complete residual bank outside every forced component.

The active coordinates are (i) owner-pair arc orientation, (ii) adjacency of
that pair along each actual path, (iii) exact control occupancy, and (iv) a
disjoint native charge reservoir. Arc orientation glues source cells; path
adjacency is the additional coordinate needed to glue their K-partners.
Neither charge nor source-frame agreement alone contains that information.

## 2. Exact shared-pair gluing criterion

Let source cells C,H each have two supplied Hamiltonian paths, indexed by
bits0,1. At a changing indexed clock in C let the owners be a,b, and at a
changing clock in H let them be c,d. Embed c,d into a,b in the unique order
that makes their directed source edge agree. Embed every other H owner as a
fresh vertex. Write epsilon=0 if c maps to a, and epsilon=1 if c maps to b.
Then the actual indexed-owner equality is exactly

    helper bit = core bit XOR epsilon.                       (1)

Let alpha_i indicate that a,b are consecutive in core path i; let beta_j
indicate that c,d are consecutive in helper path j. The two local K-partner
matrices glue on their shared pair for both allowed bit states iff

    alpha_i = beta_(i XOR epsilon),  i=0,1.                  (2)

Proof. The source cells share exactly this edge and agree by construction.
K retains a consecutive path edge and reverses every other edge. Therefore
the two K outcomes agree exactly when they both retain or both reverse the
shared edge. All remaining arcs have at most one prescribed value. Complete
unprescribed cross arcs arbitrarily. This proves sufficiency and necessity
for this particular two-owner overlap. QED.

For a star of many helpers on one core pair, pairwise agreement on that
SAME pair is also sufficient for the entire source and partner gluing:
every pair of helper cells intersects in precisely the shared pair, and
their other owners are private. Several such stars can be given disjoint
owner sets. The full intersection census in the checker rules out unintended
cross-component overlaps.

## 3. The ten-pin star

The two paths of inherited source176 are

    U: (5,7,2,3,1,4,6,0),
    V: (5,7,6,3,1,4,2,0).

Its clocks2 and6 interchange owners2 and6. That pair is nonconsecutive in
both paths, so its adjacency profile is(0,0). The core source has6→2 and
both local K partners have2→6.

The original53 helpers split into the complete40-helper triple reservoir and
the thirteen remaining sources

    368,556,753,1369,1952,2066,2135,2397,2439,2483,2638,2948,3187.

All except1952 and2638 have a changing clock whose two owner values are
nonconsecutive in BOTH supplied paths. Choose ten of the eleven eligible
helpers. The exact control rule is

    helper bit = bit176 for 368,556,753,1369,2066;
    helper bit = 1-bit176 for 2135,2397,2439,2483,3187.

Use the first eligible helper clock in the frozen receipt. Align its pair
by source-edge orientation, then choose core clock2 or6 to obtain the desired
sign in(1). Choosing between these two clocks changes the pin equation; it
does not silently relabel a fixed source matrix. Every full owner embedding
and both indexed clocks are recorded.

Exactly five controls use V in either phase. Keep helpers1952,2638,2948 at
their actual U records. Let exactly three of the protected forty use V.
The total helper V occupancy is therefore eight at every arrival. All53
helpers remain occupied, including the fixed-U slots.

These actual source cells embed in a628-owner host:81×8 local owners minus
two identifications per attached helper. The ten controlled cells and their
core form one68-owner component; the remaining cells have private owners.
Every original pin is retained through its fixed owner injection. Local
partner cells also admit coherent host completions for both phases and for
every independent reservoir choice.

## 4. Two exceptional helpers give a second native component

The exclusions above have a structural explanation: at every changing clock
of either1952 or2638, the corresponding pair is adjacent in at least one
supplied path. Hence neither can join the always-reversed176 seam while
keeping both full local K-partner choices. Matching source arcs alone would
produce a false positive.

Inherited core1532 repairs that exact obstruction. At clock2 its owners are
6 in the U path and4 in the V path. The pair is adjacent in BOTH paths, so
its profile is(1,1). Attach

    helper1952 clock3, owners6/5, with bit1952=bit1532;
    helper2638 clock4, owners1/0, with bit2638=1-bit1532.

Both helper pairs also have profile(1,1). The orientation-consistent fixed
embeddings satisfy(1) and(2), so their source and partner matrices glue. In
this component the shared edge is retained by K; in the first it is reversed.
This is a concrete use of two different local adjacency types.

The second component has one V control in either phase. Together the twelve
controls have exactly six V records, for all four independent core-bit
patterns. Keep helper2948 at U and the protected reservoir at occupancy three.
Total helper occupancy is now nine. The two components use disjoint global
owners; all81 cells embed in a624-owner source tournament. The first component
has68 owners and the second20. Their local partner matrices glue for every
allowed root pattern and every independent reservoir assignment.

This is a new packet family relative to the ten-pin family: its helper margin
is nine rather than eight, and its immutable embeddings are different.
Switching between the two families is not advertised as a free move.

## 5. Why the blocks pay together

Here is the general protected-reservoir mechanism. Let I range over a
connected authenticated core graph. Let F(I) be the exact forced-V set in a
fixed set of controlled helper slots. Assume its cardinality r0 is constant,
and every resulting full assignment satisfies the declared owner and kernel
compatibility. Let R be a disjoint reservoir with an authenticated complete
charge section at fixed occupancy k. All other slots use fixed full records.

At total charge q choose

    J(I) = section_R(q - lambda(I) - lambda(F(I)) - fixed_charge).

The complete helper state is F(I) union J(I), of cardinality r0+k. It retains
all capacities and all declared shared pins. For a core edge I→I', a
simultaneous lifted event needs at most

    |I symmetric_difference I'|
    + |F(I) symmetric_difference F(I')| + 2min(k,|R|-k)       (3)

changed slots. Its actual clock cost is the corresponding sum of the
literal per-slot path distances. A feasible prescribed initial or final
reservoir record can replace the chosen section at that endpoint, since it
already has the necessary residual charge. QED.

The crucial preserved coordinate is the protected reservoir's actual
capacity. The changing controls are not allowed to delete or occupy its
alternative normalized units. Their charges can be large or unfavorable;
the section pays their combined residual charge at each arrival. This does
not claim that any of the individual forced constraints leaves an arbitrary
unpartitioned occupancy-four bank complete.

In this native application all original28 core slots remain independent
apart from their own retained local pins. At any core margin their Johnson
graph is connected. The controlled helpers provide exactly the stated
dependence on roots176 and1532, and the reservoir provides all2048 residual
charges. Thus EVERY feasible full packet state of a given core margin and
total Q-charge lies in one connected lifted fibre, with prescribed feasible
helper endpoints retained. This includes all margins0..28; it is not an
extrapolation from a sample route.

## 6. Exact boundary costs and returning routes

The ten control paths have total clock distance61. The two exceptional paths
add11. The first component therefore has exact forced boundary cost
(10 slots,61 occurrences) whenever bit176 changes; the second contributes
(2 slots,11 occurrences) whenever bit1532 changes. If both bits change the
forced boundary is their sum. These costs cannot be erased by viewing all
constraints as observations of only two physical owner names: the work
metric counts literal slot-clock occurrences.

The reservoir needs at most six changed slots and at most42 clock occurrences.
Using the actual old core path costs gives these uniform bounds:

| Packet family | Helper V occupancy | Changed slots per event | Primary clock occurrences per event |
|---|---:|---:|---:|
|Ten-pin single component|8|18|109|
|Twelve pins, two components|9|20|120|

For the ten-star bound, changing176 costs two clocks and the other changed
core costs at most four:2+4+61+42=109. For the two-component bound, the
largest case changes176 and1532 together, costing2+4+72+42=120. Changes of
neither root or only the second root are cheaper. The checker enumerates all
old-core exchange pairs and saves the maximizing pair. These are certified
upper bounds on complete events, not claims that every charge requires the
largest reservoir replacement.

The explicit two-component core route is

    {1394,2185} → {176,2185} → {1532,2185}
      → {176,1532} → {176,1841} → {1477,1841}.

Core margin is TWO throughout; balancing slots are retained, rather than
silently flipping a single core root. The two root patterns are
00,10,01,11,10,00. The endpoint core charges agree. Begin the reservoir at
{650,676,760}, use the exact residual section at interior arrivals, and return
to that SAME reservoir at the final arrival. All controls also return to
their original bits, so every helper returns to its exact initial full record.

The five saved event(support,clock) pairs are

    (18,101), (20,114), (16,89), (8,40), (16,89).

The ten-star's shorter returning route uses three events
(18,101),(8,39),(16,89). Every arrival is checked against the original physical
Laplacian rows, all source pins, the shared-owner equations and all normalized
resources. The phase/group calculation is not used as a substitute for those
literal arrival checks.

## 7. Scope, capacity and reproducibility

Both constructions retain81 actual order-eight source slots:28 old plus53
new. The complete option reservoir has162 distinct normalized units and648
closed endpoint words; each arrival occupies81 units and324 words. The
immutable external background must avoid the entire option reservoir. Later
reuse of a released alternative is permitted, but concurrent resource reuse
or counting kernel lifts as extra capacity is not.

The full canonical-source lifts of every controlled normalized resource are
retained and independently reenumerated, including the self-converse source2439.
Every inherited common pin and complete U/V target-domain restriction remains
in place. Source/partner words and the chosen original X phase are unchanged
by the fixed global owner injections. Converse reversal transports every
shared indexed equation to its reflected clocks.

The electrical graph remains the ORIGINAL ORDER-EIGHT graph. Its full paired
code preserves every original dyadic aggregate harmonic current, as proved in
`injective_native_blocks_20261023/electrical/ELECTRICAL_CODES.md` §2. The host
tournaments above witness coherent sharing of local owners and arcs. Their
local eight-owner paths are not asserted to be Hamiltonian paths on the whole
host, and no order624/628 electrical theorem is imported. Exact integer-X,
odd-primary current and every ordinary individual-port phase remain outside
the claimed guard.

Reproduction from the repository root:

    python -B 05-knowledge/results/guarded_component_completion_20261024/native/check.py
    python -O -B 05-knowledge/results/guarded_component_completion_20261024/native/check.py

`receipt.json` contains all fixed embeddings, exact indexed equations, full
special-cell kernels, original input hashes, the complete boundary-cost
tables and every route state/event. `--write` regenerates it. Text hashes
normalize CRLF to LF; optimized execution retains all explicit checks.
