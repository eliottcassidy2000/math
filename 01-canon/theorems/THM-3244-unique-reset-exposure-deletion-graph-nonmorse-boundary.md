---
id: THM-3244
title: "Unique-reset Rips routing and one-pole deletion-flow boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the complete
  4,319-state THM-3238 physical bank, the exact exposed-reset
  functional has 32 nonreset one-pole local maxima, so global exposure does
  not induce a local deletion flow.  Its one-pole superlevel merge tree has
  33 births, 32 finite bars, and as many as 19 simultaneous components; the
  second-highest state persists as a separate component for 2,096 rank steps.
  Nevertheless, rows 2 and 10 of the lawful THM-3238 response bank form a
  state-dependent one-pole atlas: at every nonreset state at least one row has
  a strict reset-distance-monotone ascent.  Iteration reaches the reset in
  exactly the edit distance.  Among the 22 rows there are exactly 31 covering
  pairs, while no single row covers all states.  The abstract two-chart
  set-cover nerve is contractible, but 3,453 of the 3,911 overlap states have
  no common ascent direction.  Minimal chart-switch depth has histogram
  716/3,600/2 at depths 0/1/2, so two switches are sharply necessary.  No
  fixed positive blend of rows 2 and 10 removes the switching: two explicit
  states have complementary rational trap intervals covering every positive
  blend ratio, and two states are minimal.
  The reset-distance-monotone escape radius is at most 10, with a unique sharp
  radius-10 trap; at radius 10 the reset is the unique directed sink.  More
  strongly, two explicit hubs route every state to the reset in at most two
  strict, reset-distance-monotone radius-10 jumps, and no one-hub two-level
  atlas can do so.  This is a finite directed-graph theorem, not a simplicial
  collapse or a result for another support/bank.
source: root/multiscale-newton-flag/2026-08-03
audit: >
  The exact companion reconstructs the canonical 4,319 physical
  submultisets, pins the immutable THM-3238 script and output, verifies an
  independently extracted exact rank certificate for all THM-3238 integer
  coordinates, and exhausts every one-pole neighbor, every
  reset-distance-monotone upward pair, every possible radius-10 hub, and the
  complete elder-rule
  union-find merge tree of the one-pole superlevel filtration.  A second
  independently extracted bit certificate verifies all 22 lawful one-pole
  row covers, the distinguished row-(2,10) atlas, the 31-pair census, and the
  no-single-row boundary.  A direction-bit certificate additionally checks
  the overlap graph, common-direction histogram, and exact switch-depth
  dynamic program.  Exact edge differences verify the sharp two-state
  constant-blend obstruction.  Two independent audits compared the rank,
  cover, direction, and edge-difference certificates pointwise against the
  full exact response cache.  They confirmed pairwise distinctness, the unique
  maximum, and zero rank-certificate mismatches; rederived the merge tree,
  induced graphs, dynamic program, rational trap intervals, and scope
  implications; and matched every cover, direction, and displayed edge
  difference.  One reconstruction of all 22 unscaled response rows on all
  4,319 states gave exact row-bank digest
  366905601b960854d249e1f12ce02edeb55af17a638bb56673bd49c6dbba26e9.
  Fresh normal and optimized runs byte-match the stored output and declared
  hashes.  The audits repaired the folded reset-distance-preorder, local-flow,
  one-row-universe, set-cover-nerve, route-count, and scalar-gauge boundaries.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3216-depth-nine-degree-fourteen-unique-reset-face-and-omega-cone-boundary
script: 04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py
output: 05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out
script_sha256: 3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba
output_sha256: 27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9
hash_basis: LF-normalized bytes
---

# THM-3244 -- unique-reset Rips routing and one-pole deletion-flow boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3238 exposes the reset

```text
Q=(1,3,3,4,5,6,7,8)                                    (1)
```

as the unique maximum of an exact lawful functional `H` on all 4,319
nonempty physical submultisets of

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                  (2)
```

Exposure is global.  This theorem asks the different, local question: does
the order defined by `H` flow to `Q` through small pole edits?  The answer is
negative at one-pole scale but sharply positive at edit radius 10.  There is,
however, a different one-pole repair: retain the individual lawful response
coordinates and allow a state-dependent choice between two charts.

## 1. Physical edit metric and directed Rips graphs

Write `n_j(sigma)` for the multiplicity of `j` in a physical state `sigma`.
Thus

```text
0 <= (n_1,...,n_8) <= (4,3,2,2,2,1,1,1),               (3)
```

coordinatewise, with the zero vector deleted.  Put

```text
d(sigma,tau)=sum_(j=1)^8 |n_j(sigma)-n_j(tau)|.          (4)
```

The one-pole graph joins states at distance one.  Say that `tau` is
reset-distance-monotone from `sigma`, and write `tau <=_Q sigma`, when

```text
|n_j(tau)-n_j(Q)| <= |n_j(sigma)-n_j(Q)|                (5)
```

for every `j`.  This is a folded coordinatewise preorder on absolute distance
from `Q`; it is not a single orthant and may cross `Q`.  For an integer
`r>=1`, define the directed reset-Rips graph
`R_r^+(H,Q)` by

```text
sigma --> tau
 iff tau <=_Q sigma, d(sigma,tau)<=r, and H(tau)>H(sigma). (6)
```

The name records only this threshold graph.  No simplicial Vietoris--Rips
complex is being asserted to collapse.

THM-3238's 4,319 exact integer values are pairwise distinct.  Their canonical
state order has digest

```text
c060e22f900232b608ad3bce1f6b24cae51b6eb45c138b679f4c698fbad2c6a2, (7)
```

and the compressed rank permutation used by the companion has raw digest

```text
ac7aeaf94958c04034137dbbffaf18f494485446421ad9006b6657a15b3487d7. (8)
```

Rank zero is the largest coordinate, namely `Q`.  The certificate is an exact
ordering of the full THM-3238 integer vector, not a floating approximation;
the companion pins the dependency script and output hashes before using it.

## 2. One-pole exposure does not induce a local deletion flow

Exhausting every legal insertion and deletion in the nonempty physical bank
gives exactly

```text
32                                                           (9)
```

nonreset states whose `H`-value is strictly larger than the value at every
one-pole neighbor.  Hence the unique global maximum `Q` does not by itself
supply an `H`-increasing one-pole edge from every nonreset state and therefore
does not induce a local deletion flow.  No discrete-Morse conclusion is drawn.

The reset-distance-monotone obstruction is larger still.  For `sigma!=Q`, define its
escape radius

```text
epsilon(sigma)=min {d(sigma,tau):
                    tau<=_Q sigma and H(tau)>H(sigma)}.  (10)
```

The minimum exists because `tau=Q` is always admissible.  The exact histogram
of `(10)` is

```text
radius       1    2   3   4  5  6  7  8  9 10
states    4231   45  23  10  1  2  2  2  1  1.          (11)
```

Thus 87 nonreset states have no radius-one reset-distance-monotone ascent,
even though only 32 of them are local maxima in the full one-pole graph.  This
separates failure under the folded reset-distance preorder from failure of
every local direction.

### The order-persistent obstruction

Filter the one-pole graph by adding states in decreasing `H` order and use
the elder rule when two connected components merge.  Equivalently, rank zero
is `Q`, rank one is the next largest exact coordinate, and the filtration at
rank `k` is induced by the first `k+1` states.  The exact zero-dimensional
merge tree has

```text
33 component births, 32 finite bars,
maximum simultaneous component count 19.               (11a)
```

The longest finite rank bar is

```text
birth state (1), birth rank 1,
merge state (1,1,1,1,2,2,2,3,3,4,5,8), merge rank 2097,
rank persistence 2096.                                  (11b)
```

Thus the second-highest global state is not adjacent to the reset basin in
the superlevel graph: its component remains separate until nearly half the
bank has entered.  Two further named bars anticipate the coarse atlas below:

```text
(3): birth 16, merge at (1,3) in rank 597, persistence 581;
T:   birth 12, merge in rank 496,              persistence 484. (11c)
```

The first state in `(11c)` is the lower hub of Section 4 and `T` is the sharp
radius-10 trap.  The upper hub has rank 36 but is not a component birth.  This
rank persistence uses only the exact order of `H`, so it is invariant under
strictly increasing reparameterization; it is not a claim about gaps between
the enormous integer coordinates.

## 3. Sharp radius-10 sink filtration

Since `H` strictly increases along every directed edge, `R_r^+(H,Q)` is
acyclic.  Its number of sinks for `r=1,...,10` is

```text
(88,43,20,10,9,7,5,3,2,1).                              (12)
```

The unique state with escape radius 10 is

```text
T=(1,1,1,1,2,2,2,3,3,4).                               (13)
```

Consequently `Q` and `T` are both sinks at radius 9, while `Q` is the unique
sink at radius 10.  Every state therefore has a directed path to `Q` in
`R_10^+(H,Q)`.  Radius 10 is sharp for this fixed metric, order, and folded
reset-distance preorder.

## 4. A two-hub, two-jump atlas

The unique-sink conclusion does not bound path length.  There is nevertheless
an explicit stronger routing.  Every state with `d(sigma,Q)<=10` jumps
directly to `Q`.  Exactly 133 states lie farther away.  Put

```text
A=(3),
B=(1,1,2,2,3,3,4,5,6,7,8).                             (14)
```

Their distances to `Q` are 7 and 3.  Among the 133 far states:

```text
A is a strict reset-distance-monotone radius-10 ascent from  36 states,
B is a strict reset-distance-monotone radius-10 ascent from 127 states,
both are available from                              30 states. (15)
```

The two coverage sets therefore have union size

```text
36+127-30=133.                                          (16)
```

After the first jump, either hub jumps directly to `Q`.  Hence every one of
the 4,319 states reaches `Q` by at most two edges of `R_10^+(H,Q)`.

This two-hub claim is minimal within such a two-level radius-10 atlas.  Any
hub which itself jumps to `Q` must lie within radius 10 of `Q`.  Exhausting all
such states, one hub covers at most 127 of the 133 far states.  Exactly two
single hubs attain that bound:

```text
(1,1,2,3,4,5,6,7,8),
(1,1,2,2,3,3,4,5,6,7,8).                               (17)
```

Thus no one-hub two-jump routing exists.  The particular lower hub `A` in
`(14)` supplies the six far states missed by the chosen upper hub `B`.

## 5. A lawful two-chart atlas at one-pole scale

The failure in Section 2 belongs to the single combined functional `H`, not
to the full lawful response bank.  Number the 22 positive upset rows exactly
as in THM-3238 and write

```text
f_i=r_(N_i,U_i)
```

for the unscaled response in `(5)` of that theorem.  Every stitched multiplier
is positive, so this rescaling does not change a strict comparison.  Thus row
2 is the degree-14 upset with minimal generators

```text
(3,2,1^9), (2,2,2,1^8),
```

and row 10 is the degree-11 principal upset generated by `(2,1^9)`.  Define

```text
C_i={sigma != Q: some one-pole tau <=_Q sigma has f_i(tau)>f_i(sigma)}. (18)
```

The exact target-bank split for rows 2 and 10 is

```text
|C_2|=4215, |C_10|=4014, |C_2 intersect C_10|=3911,
|C_2\C_10|=304, |C_10\C_2|=103, uncovered=0.          (19)
```

Consequently every nonreset state has a strict reset-distance-monotone
one-pole ascent in at least one of these two lawful coordinates.  Such a move
decreases `d(sigma,Q)` by exactly one.  Reapplying `(19)` therefore routes every
state to `Q` in exactly `d(sigma,Q)` moves.

This pair is not an isolated coincidence.  The complete covering-pair census
is

```text
(2,7) (2,10) (2,13) (2,17) (2,19) (2,21)
(3,9) (3,11) (3,16) (3,22)
(7,14) (7,18) (7,22)
(10,11) (10,16) (10,22)
(11,13) (11,17) (11,21)
(12,13) (12,19)
(13,14) (13,18) (13,22)
(14,19)
(16,17) (16,21)
(17,22) (18,19) (19,22) (21,22).                    (20)
```

There are exactly 31 pairs.  No one-row atlas exists within these 22 THM-3238
rows: the best single row is row 22, which covers 4,227 of the 4,318 nonreset
states and misses 91.  This does not exhaust every possible lawful upset
response outside the stitched bank.

The selector in `(19)` is load-bearing.  Its choice of row may change after
every pole edit.  The certificate `(19)` does not furnish one positive lawful
dual, a convex combination of rows, a probability/Markov transport, or a
scalar height function.  It furnishes a finite state-dependent local gauge
atlas whose transition problem remains separate.

### Transition geometry and sharp switch depth

The abstract set-cover nerve of `{C_2,C_10}` has two vertices and one edge, so
it is contractible.  This is a statement about the finite cover, not a nerve
lemma for a topological or simplicial model.  Separately, the overlap is
connected in the physical one-pole graph.  The three induced-region statistics
`(vertices, components, edges, cycle rank)` are

```text
C_2\C_10: (304,29,  603,  328),
C_10\C_2: (103, 6,  179,   82),
C_2∩C_10: (3911,1,18554,14644).                         (21)
```

This simple nerve hides a genuine directional obstruction.  For a state in
`C_i`, let `D_i(sigma)` be the set of pole coordinates whose `Q`-directed
one-pole move strictly raises `f_i`.  Across the 3,911 overlap states, the
histogram of `|D_2 intersect D_10|` is

```text
common directions       0     1    2   3  4
overlap states        3453   294  113  45  6.            (22)
```

Thus 3,453 overlap states have no common local ascent direction.  More
globally, minimize the number of row-label switches along a complete
reset-distance-monotone one-pole ascent route.  Exact dynamic programming
along the `Q`-distance grading gives

```text
minimum switches        0     1  2
nonreset states        716  3600  2.                     (23)
```

There are 534 nonreset states admitting an all-row-2 route and 182 admitting
an all-row-10 route; those state sets are disjoint.  The two states whose
minimum is exactly two chart-label switches are

```text
(1,2,2,3,3,4,4,5,5,6,7,8),
(2,2,3,3,4,4,5,5,6,7,8).                               (24)
```

The number in `(23)` is a finite transition-depth invariant of this directed
atlas.  It is not a Čech class or monodromy representation: defining either
would require canonical transition maps on the overlap, while `(22)` shows
that even a common ascent direction is usually absent.

### A sharp constant-gauge clutch obstruction

The switching cannot be removed merely by replacing the two charts with one
fixed positive blend.  Put

```text
H_lambda=lambda f_2+f_10,  lambda>0.                    (25)
```

A state is trapped for `(25)` when every `Q`-directed one-pole edge `e`
satisfies `Delta_e H_lambda<=0`; equality belongs to the trap because ascent
is strict.  At

```text
A=(2,2,2,3,3,4,5,5,6,7,8),                            (26)
```

the only edge with positive `Delta f_2` is delete-2, with

```text
(Delta f_2,Delta f_10)
 =(647427527551915200,-82016379613632).
```

The other two directed edges have both differences negative.  Hence `A` is
trapped throughout

```text
0<lambda<=U=427168643821/3372018372666225.              (27)
```

At

```text
B=(1,3,3,4,5),                                         (28)
```

the insert-6, insert-7, and insert-8 edges all have negative `Delta f_2` and
positive `Delta f_10`.  Insert-6 gives the largest zero, so `B` is trapped
throughout

```text
lambda>=L=44548722230872990/295146673301558860624447.   (29)
```

Direct rational comparison gives `L<U`.  Thus `(27)` and `(29)` cover every
positive `lambda`, proving that no constant positive combination of rows 2
and 10 supplies a reset-distance-monotone one-pole ascent everywhere.

Two trap states are sharp.  No single state can remain trapped in both pure
limits: `(19)` says every state has an `f_2`-raising or `f_10`-raising edge,
which breaks the trap for sufficiently large or sufficiently small positive
`lambda`, respectively.  This is a clutch obstruction inside the two-row
cone only; it does not rule out every scalar combination of all 22 rows.

## 6. Exact proof and hostile boundary

The companion performs only integer and finite combinatorial operations.  It
reconstructs the complete capacity box `(3)`, checks the state digest `(7)`,
decompresses and verifies the exact rank permutation `(8)`, and then exhausts:

1. every one-pole neighbor for `(9)`;
2. every reset-distance-monotone upward pair for `(11)` and `(12)`;
3. the sharp state `(13)`;
4. every far-state incidence with the two hubs in `(14)`; and
5. every possible hub within radius 10 of `Q` for the minimality claim `(17)`;
6. the elder-rule union-find merge tree `(11a)--(11c)`; and
7. all 22 lawful row covers, all 231 unordered row pairs, and the distinguished
   split `(19)`; and
8. the induced graphs `(21)`, direction histogram `(22)`, and switch dynamic
   program `(23)--(24)`; and
9. the exact edge differences and complementary intervals `(25)--(29)`.

The independently extracted rank certificate is deliberately separate from
the expensive product-Gamma reconstruction in THM-3238.  The hostile audits
compared it against that full exact vector; matching only the certificate's
own digest would not have established provenance.  They also independently
reconstructed the complete 22-row exact bank before comparing every row-cover,
direction, and displayed edge-difference value.  Thus none of the row-atlas
certificates is promoted from a self-digest alone.

## 7. Scope and holotopy reading

This is a theorem about the fixed support `(1,3)`, bank `I2`, physical
submultiset metric `(4)`, folded reset-distance preorder `(5)`, and exact
THM-3238 functional.
It proves neither an analogous radius for the other maintained faces nor a
uniform statement for arbitrary radial coefficients.  It also gives no
probability transport, Markov deletion law, single scalar local gauge,
simplicial collapse, deformation retract, or proof of the Gaussian Moment
Conjecture.

The useful holotopy lesson is more precise.  A global exposed section can fail
to localize on the elementary deletion cover in two inequivalent ways.  One
may enlarge metric scale while keeping one height (`H` first has a unique sink
at radius 10), or retain elementary scale while enlarging the gauge (`f_2`
and `f_10` cover every one-pole obstruction).  The latter repair is strictly
local but nonlinear because its chart changes with the state, sharply twice
on `(24)`.  THM-3160's flat order-independent pole holotopy, the `H`-oriented
coarse routing, and the row-selector atlas are therefore three different
structures.  A genuine identification still requires a lawful transition or
positive common gauge; none is supplied here, and `(25)--(29)` rule out the
cheapest constant gauge inside the distinguished two-row cone.

QED.
