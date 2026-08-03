---
id: THM-3244
title: "Unique-reset Rips routing and deletion-graph non-Morse boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.  On
  the complete 4,319-state THM-3238 physical bank, the exact exposed-reset
  functional has 32 nonreset one-pole local maxima, so global exposure does
  not induce a local deletion flow.  Its one-pole superlevel merge tree has
  33 births, 32 finite bars, and as many as 19 simultaneous components; the
  second-highest state persists as a separate component for 2,096 rank steps.
  Nevertheless, rows 2 and 10 of the lawful THM-3238 response bank form a
  state-dependent one-pole atlas: at every nonreset state at least one row has
  a strict reset-monotone ascent.  Iteration reaches the reset in exactly the
  edit distance.  There are exactly 31 covering row pairs, while no single
  row covers all states.  The ordinary two-chart nerve is contractible, but
  3,453 of the 3,911 overlap states have no common ascent direction.  Minimal
  chart-switch depth has histogram 716/3,600/2 at depths 0/1/2, so two
  switches are sharply necessary.
  The reset-monotone escape radius is at most 10, with a unique sharp
  radius-10 trap; at radius 10 the reset is the unique directed sink.  More
  strongly, two explicit hubs route every state to the reset in at most two
  strict, reset-monotone radius-10 jumps, and no one-hub two-level atlas can
  do so.  This is a finite directed-graph theorem, not a simplicial collapse
  or a result for another support/bank.
source: root/multiscale-newton-flag/2026-08-03
audit: >
  The exact companion reconstructs the canonical 4,319 physical
  submultisets, pins the immutable THM-3238 script and output, verifies an
  independently extracted exact rank certificate for all THM-3238 integer
  coordinates, and exhausts every one-pole neighbor, every reset-monotone
  upward pair, every possible radius-10 hub, and the complete elder-rule
  union-find merge tree of the one-pole superlevel filtration.  A second
  independently extracted bit certificate verifies all 22 lawful one-pole
  row covers, the distinguished row-(2,10) atlas, the 31-pair census, and the
  no-single-row boundary.  A direction-bit certificate additionally checks
  the overlap graph, common-direction histogram, and exact switch-depth
  dynamic program.  Normal and optimized runs reproduce the stored output.
  Independent hostile theorem audit is pending.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3216-depth-nine-degree-fourteen-unique-reset-face-and-omega-cone-boundary
script: 04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py
output: 05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out
script_sha256: b05b35ce006f8adc52822ac325b2769f1ac979ed03b833e4f379c3730ef2e636
output_sha256: e79fbeb23a53940dd11435d70a36b708b2b8feee4d61bfa6d27e6a026748db5f
hash_basis: LF-normalized bytes
---

# THM-3244 -- unique-reset Rips routing and deletion-graph non-Morse boundary

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

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
`Q`-monotone from `sigma`, and write `tau <=_Q sigma`, when

```text
|n_j(tau)-n_j(Q)| <= |n_j(sigma)-n_j(Q)|                (5)
```

for every `j`.  For an integer `r>=1`, define the directed reset-Rips graph
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

## 2. One-pole exposure is not Morse

Exhausting every legal insertion and deletion in the nonempty physical bank
gives exactly

```text
32                                                           (9)
```

nonreset states whose `H`-value is strictly larger than the value at every
one-pole neighbor.  Hence the unique global maximum `Q` does not produce a
strict one-pole ascent, a deletion gradient, or an `H`-induced discrete-Morse
matching.

The reset-monotone obstruction is larger still.  For `sigma!=Q`, define its
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

Thus 87 nonreset states have no radius-one reset-monotone ascent, even though
only 32 of them are local maxima in the full one-pole graph.  This separates
failure of the chosen reset orthant from failure of every local direction.

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
`R_10^+(H,Q)`.  Radius 10 is sharp for this fixed metric, order, and monotone
orthant.

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
A is a strict reset-monotone radius-10 ascent from  36 states,
B is a strict reset-monotone radius-10 ascent from 127 states,
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
as in THM-3238 and write `f_i` for the unscaled response of row `i`.  Thus row
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

Consequently every nonreset state has a strict reset-monotone one-pole ascent
in at least one of these two lawful coordinates.  A one-pole monotone move
decreases `d(sigma,Q)` by exactly one.  Reapplying `(19)` therefore routes
every state to `Q` in exactly `d(sigma,Q)` moves.

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

There are exactly 31 pairs.  No one-row atlas exists: the best single row is
row 22, which covers 4,227 of the 4,318 nonreset states and misses 91.

The selector in `(19)` is load-bearing.  Its choice of row may change after
every pole edit.  Thus `(19)` is not one positive lawful dual, a convex
combination of rows, a probability/Markov transport, or a scalar height
function.  It is a finite local gauge atlas whose transition problem remains
separate.

### Transition geometry and sharp switch depth

The ordinary nerve of the cover `{C_2,C_10}` has two vertices and one edge,
so it is contractible.  In fact the overlap is connected in the physical
one-pole graph.  The three induced-region statistics `(vertices, components,
edges, cycle rank)` are

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
`Q`-monotone one-pole ascent route.  Exact dynamic programming along the
`Q`-distance grading gives

```text
minimum switches        0     1  2
nonreset states        716  3600  2.                     (23)
```

There are 534 all-row-2 routes and 182 all-row-10 routes.  The two states
sharply requiring two switches are

```text
(1,2,2,3,3,4,4,5,5,6,7,8),
(2,2,3,3,4,4,5,5,6,7,8).                               (24)
```

The number in `(23)` is a finite transition-depth invariant of this directed
atlas.  It is not a Čech class or monodromy representation: defining either
would require canonical transition maps on the overlap, while `(22)` shows
that even a common ascent direction is usually absent.

## 6. Exact proof and hostile boundary

The companion performs only integer and finite combinatorial operations.  It
reconstructs the complete capacity box `(3)`, checks the state digest `(7)`,
decompresses and verifies the exact rank permutation `(8)`, and then exhausts:

1. every one-pole neighbor for `(9)`;
2. every reset-monotone upward pair for `(11)` and `(12)`;
3. the sharp state `(13)`;
4. every far-state incidence with the two hubs in `(14)`; and
5. every possible hub within radius 10 of `Q` for the minimality claim `(17)`;
6. the elder-rule union-find merge tree `(11a)--(11c)`; and
7. all 22 lawful row covers, all 231 unordered row pairs, and the distinguished
   split `(19)`; and
8. the induced graphs `(21)`, direction histogram `(22)`, and switch dynamic
   program `(23)--(24)`.

The independently extracted rank certificate is deliberately separate from
the expensive product-Gamma reconstruction in THM-3238.  An independent
audit must compare it against that full exact vector; matching only the
certificate's own digest would not establish provenance.  The same boundary
applies to the row-cover and direction certificates: their finite arithmetic
is self-contained, but their lawful-response provenance requires comparison
with the exact THM-3238 rows.

## 7. Scope and holotopy reading

This is a theorem about the fixed support `(1,3)`, bank `I2`, physical
submultiset metric `(4)`, reset orthant `(5)`, and exact THM-3238 functional.
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
on `(24)`.  THM-3160's
flat order-independent pole holotopy, the `H`-oriented coarse routing, and the
row-selector atlas are therefore three different structures.  A genuine
identification still requires a lawful transition or positive common gauge;
none is supplied here.

QED.
