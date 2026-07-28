---
id: THM-2680
title: "Dilation-reversed two-edge fibre products, Fibonacci cospan, and depth-three clock obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On the canonical THM-2616 guard-safe raw carrier, D(x)={13x} identifies a
  current rail-owner clock with the following shallow clock and satisfies
  j(Dx)=h(x).  This reverses the printed clock arrows and makes the first
  genuine two-edge object E_0 intersect D^{-1}E_1, with both source labels
  retained.  Exact refined-grid computation finds two positive source-pair
  families whose s_1=6 and s_1=7 zero columns are inherited one-edge
  exceptions.  The complete safe/safe atlas has 17,160 positive labelled
  rows.  Safe/danger and danger/safe add 4,488 and 3,696, while all 27,640
  formal danger/danger chains are base-support zeros.  Every positive sector
  pair occupies the same ten clock triples.  Their supports are pairwise
  disjoint and total 25,344.  THM-2682 independently proves, by the central
  arrival clock identity, that every three-event D-chain on this carrier is
  empty; THM-2684 extends the clock-diagonal no-go to both endpoint arrivals
  and the full THM-2584 rail bank.  Different parent carriers/handoffs,
  configuration switching, source transport there, and LRC(14) remain open.
source: root-2026-07-28-dilation-clock-fibre-product
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2670-sharp-graph-clock-incidence-atlas-and-physical-gluing-boundary
related:
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2644-odd-torsor-purity-return-gate-and-nonlinear-fixed-branch-decoder
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2682-central-arrival-clock-trap-and-three-event-dilation-nilpotence
  - THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence
script: 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
output: 05-knowledge/results/lrc14_dilation_reversed_clock_fibre_product_probe.out
script_sha256: 56347893caa028a8d9f1b72c6e886dd38c617019f53e746e966bee65cb91a2ad
output_sha256: 1e7fb8babc74b52a4c18e05ff3604fae3282641df6bb4392a14db5bde266701a
secondary_script: 04-computation/lrc14_clock_handoff_d_pullback_scout.py
secondary_output: 05-knowledge/results/lrc14_clock_handoff_d_pullback_scout.out
secondary_script_sha256: 1870550958dbff99122c1c222fc6719de6b971341b50dba59bb889e1a7baab0d
secondary_output_sha256: 75e1733bee6d1bd25d857e0d989c0f0730a9d3676ccd158bacc51576fbecaffd
referee_engine: 04-computation/lrc14_thm2670_physical_reversed_two_edge_scout.py
referee_script: 04-computation/lrc14_clock_handoff_d_pullback_independent_referee.py
referee_output: 05-knowledge/results/lrc14_clock_handoff_d_pullback_independent_referee.out
referee_engine_sha256: a017c6e96f99f23f8067db28084f52dd6509f421b0c617e0643a8250cce0c474
referee_script_sha256: 2df48269f02e9dfc0f063be6fa499732d51c4f85244b3dba6f141ebc8f7026f7
referee_output_sha256: 7146a0c771c15e3e88d7b21cab8fcc3635c94588547a571bba276f4b0d987c2f
hash_basis: LF-normalized bytes
---

# THM-2680 -- dilation gives a physical two-edge fibre product

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2670 isolates the typing error in the old clock products: different
rail-owner charts are disjoint, so ordinary Boolean multiplication joins
labels but not witnesses.  There is nevertheless a canonical map between
scales.  Pulling the following event back by base-thirteen dilation produces
a genuine common-`x` two-edge object.  It is nonempty, sparse, and sharply
source-sensitive.

## 1. Dilation reverses the clock arrows

Let

```text
D(x)={13x}.                                                (1)
```

Write `O_b` for the THM-2616 rail-owner clock cell labelled `b` and `S_a`
for its shallow clock cell labelled `a`.  Exact half-open interval arithmetic
gives

```text
O_b = D^(-1)(S_b),                  b in F_7.              (2)
```

The seven owner cells are pairwise disjoint.  The companion checks all `42`
ordered disjoint pairs and all seven identities `(2)`, including the wrapped
clock-zero presentation.

The retained numerical states are consecutive base-thirteen digits:

```text
j(x)=floor(13^6 x) mod 13,
h(x)=floor(13 {13^6 x}).                                  (3)
```

Hence, away from the measure-zero endpoints fixed by the half-open convention,

```text
j(Dx)=h(x).                                                (4)
```

Thus a `D`-compatible candidate chronology traverses THM-2670's printed
`(owner,shallow)` arrow backwards.  If its chronological clock path is

```text
a -> b -> c,                                               (5)
```

the current stored atom has `(owner,shallow)=(b,a)` and the following atom
has `(owner,shallow)=(c,b)`.

## 2. The correctly typed fibre product

For a rail occurrence `alpha`, source `s`, binary labels
`epsilon,kappa`, and sharp-graph states `j,h`, let

```text
E^(s)_(b,a;alpha,j,h,epsilon,kappa) subset O_b             (6)
```

be the unintegrated guard-safe atom used in THM-2670.  It is the intersection
of the THM-2616 middle rail, present packet, sharp nonzero deep root,
half-tooth, predecessor digit, carry half, and delayed-word digit.  Every
factor is nonnegative.

The direct intersection of successive stored atoms lies in disjoint owner
cells and is empty.  The `D`-typed two-edge object is instead

```text
F=E^(s0)_(b,a;alpha,j,h,epsilon,kappa)
  intersect
  D^(-1) E^(s1)_(c,b;alpha',h,k,epsilon',kappa').          (7)
```

The intermediate equality in `(7)` is exactly `(4)`.  Both source labels are
retained independently; no source-transition law is assumed.

The exact computation represents `(7)` on the refined integer grid `13T`.
If

```text
y={13^6x},                                                 (8)
```

the two delayed words reduce to the exact interval intersection

```text
Q_A(y) intersect Q_B({13y}).                              (9)
```

The following atom's remaining step profile is pulled back branch by branch
under `D`; the current profile is rescaled to the same grid.  Multiplying the
two positive integer weights preserves support.  The resulting numerator is
therefore an exact nonemptiness certificate, but is **not** asserted to be a
canonical transition mass.  Eight independent one-edge controls verify that
refining both grids multiplies the inherited numerator by exactly `13`.

## 3. Two positive source-pair families

Use the chronological triple notation `(shallow,current owner,next owner)`.
For

```text
(a,b,c)=(3,1,0),                                           (10)
```

the number of positive rail-labelled atom pairs in `(7)` is

```text
20  for 132 of the 144 ordered source pairs,
 0  exactly for (s0,s1)=(s,6),  s=1,...,12.               (11)
```

An exact positive numerator certificate at `(s0,s1)=(1,1)` is

```text
126816337986097204341478787325120                         (12)
```

and current/following labels

```text
(alpha,j,h,epsilon,kappa)=(2,5,2,0,0),
(alpha',j',h',epsilon',kappa')=(0,2,6,1,1).               (13)
```

In particular `h=j'=2`, as required.

For

```text
(a,b,c)=(4,0,1),                                           (14)
```

the corresponding census is

```text
14  for 132 ordered source pairs,
 0  exactly for (s0,s1)=(s,7),  s=1,...,12.               (15)
```

The first `(1,1)` positive numerator certificate is

```text
90407221954729505049108909024000                          (16)
```

and labels

```text
(1,11,3,0,0),                 (3,3,5,1,0).                (17)
```

Again `h=j'=3`.  The exceptional columns in `(11)` and `(15)` are inherited
one-edge zeros: the following banks `(owner,shallow)=(0,1)` at `s1=6` and
`(1,0)` at `s1=7`, respectively, are already empty before the `D`-pullback.
Thus these columns are not a new source-drift obstruction.  Both source
labels remain independently retained, and no source-transport rule follows.

## 4. Sharp hostile and sparse clock atlas

The nearby triple

```text
(a,b,c)=(3,1,2)                                             (18)
```

has zero fibre-product atoms for all `144` ordered source pairs, but for an
inherited reason: its following one-edge cell `(owner,shallow)=(2,1)` is
universally empty.  This is a useful zero-edge control, not a hostile to
`D`-gluing itself.

The sharp hostile is instead

```text
(a,b,c)=(0,1,0).                                           (19)
```

Its `D`-fibre product is empty for all `144` source pairs.  At `(s0,s1)=(1,1)`,
however, both one-edge relations are positive and their formal Boolean
product has support `29`.  Thus the failure in `(19)` occurs only after the
actual witness pullback and cannot be explained by a missing edge or an
incompatible printed state label.

At the generic fixed source pair `(s0,s1)=(1,1)`, exhaust all

```text
7*6*6=252                                                  (20)
```

ordered clock triples with distinct adjacent clocks.  Exactly `146` have
both one-edge relations positive, and all `146` also have positive formal
Boolean product.  Only ten have a positive physical `D`-fibre product.  The
tuple records `(a,b,c,positive rail-labelled atom-pair occurrence count)`:

```text
(3,0,3,6),  (3,0,4,10), (3,0,5,10), (3,0,6,12),
(3,1,0,20), (4,0,1,14), (4,0,2,14), (4,0,3,14),
(4,0,4,12), (4,6,0,18).                                  (21)
```

Only shallow clocks `3,4` occur in `(21)`.  This concentration is explained
and promoted by THM-2682: the central arrival return preserves precisely
those two half-clock labels and kills every three-event continuation.

## 5. Complete cospan and the empty physical two-simplex

The independent complete atomwise audit extends `(21)` across both delayed
guard sectors.  The ordered sector-pair census is

```text
safe   -> safe:     17,160 positive labelled rows,
safe   -> danger:    4,488 positive labelled rows,
danger -> safe:      3,696 positive labelled rows,
danger -> danger:        0 positive labelled rows.       (22)
```

For the safe/safe leg, there are `32,736` base-positive labelled pairs;
`17,688` also have nonempty joint delayed support, and `528` of those have
zero exact prefix integral.  The final numerator content is

```text
811099927240490160.                                      (23)
```

Each of the ten triples in `(21)` has `132` active safe/safe source pairs.
The first five exclude precisely `s_1=6`, and the last five exclude precisely
`s_1=7`, with no restriction on `s_0`.  The respective physical state counts
per active pair are

```text
6,10,10,12,20,14,14,14,12,18.                           (24)
```

The mixed finite-universe accounting is

```text
safe -> danger:
  formal=67,729, base-positive=4,752,
  physical=4,488, skew-zero=264;

danger -> safe:
  formal=141,593, base-positive=5,016,
  word-empty=1,320, physical=3,696.                       (25)
```

Both mixed directions realize all `144` ordered source pairs somewhere.
Their per-clock positive-row counts in the order `(21)` are

```text
safe -> danger:  264,528,528,528,264,528,528,528,264,528,
danger -> safe:  660,660,660,660,396,132,132,132,132,132. (26)
```

The danger square has a stronger zero than THM-2670's formal `U^3=0`:
all `27,640` formal danger/danger state chains have

```text
supp(E_0) intersect D^(-1)supp(E_1)=empty                (27)
```

before either delayed word is consulted.

At the richest triple `(3,1,0)`, the four sector legs have physical state
supports

```text
safe->safe=20, safe->danger=2, danger->safe=3,
danger->danger=0.                                        (28)
```

The three positive supports are disjoint and give the guard-free support
`25`.  After retaining only sector type, `(28)` is exactly

```text
          next safe   next danger
safe          1            1
danger        1            0,                            (29)
```

the Fibonacci adjacency matrix.  This binary grammar does not iterate,
because every positive leg in `(22)` occupies **exactly** the same ten clock
triples `(21)`.  Let `C` denote this common set.  Then

```text
{a:(a,b,c) in C}={3,4},
{b:(a,b,c) in C}={0,1,6}.                                (30)
```

Consequently there are no `(a,b,c),(b,c,d)` in `C`.  This gives an independent
finite-atlas proof of THM-2682's structural conclusion: any positive three-edge
event

```text
E_0 intersect D^(-1)E_1 intersect D^(-2)E_2              (31)
```

would project to both adjacent positive two-edge events, a contradiction.
Thus the full physical guard-cospan handoff has a nonempty one-skeleton but
no two-simplex.  THM-2682 explains this more strongly by the central-arrival
clock trap, before the sector restrictions.  Neither proof needs a
source-transition hypothesis.

## 6. What this changes, and what remains

THM-2670's fixed-step sevenfold Boolean products vanish in the stored arrow
orientation.  Equations `(10)`--`(17)` prove that this cannot be promoted to
a blanket physical chronology no-go: after applying the correctly oriented
`D` handoff, genuine two-edge common-witness fibre products exist.  The
inherited source exceptions `(11)`, `(15)` and the sharp physical/formal gap
`(19)`--`(21)` simultaneously show why one positive overlap is far from a
transition.

Equation `(31)` is now proved empty, so iterating this particular quotient is
finished.  A replacement must change the object: a common endpoint fibre
product, a larger source/rail natural extension, a configuration-switching
atlas, or a phase-retaining correspondence not equal to `D`.  No Bockstein
unit, endpoint owner, Perron-sheet selector, global positive transition,
holonomy trivialization, scalar row exclusion, or LRC(14) conclusion follows.

## 7. Reproduction

Run

```bash
python3 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
python3 -O 04-computation/lrc14_dilation_reversed_clock_fibre_product_probe.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_dilation_reversed_clock_fibre_product_probe.out.
```

The computation uses exact integer interval arithmetic only.  Its finite
universe, positive controls, all-source hostile, refined-grid formula checks,
and current hashes are recorded above.

The independent cospan/danger companion is

```bash
python 04-computation/lrc14_clock_handoff_d_pullback_scout.py
python -O 04-computation/lrc14_clock_handoff_d_pullback_scout.py
```

Both runs byte-match
`05-knowledge/results/lrc14_clock_handoff_d_pullback_scout.out`.  It rebuilds
the `7,436`, `1,636`, and `8,360` positive atom occurrences, checks the exact
`13T` prefix formula, exhausts `(27)`, and resolves all four legs in `(28)`.
Its LF-normalized hash pair is the secondary pair in the frontmatter.  A
separate implementation independently reconstructed the global rows
`(22)`--`(26)` and the nonconcatenation `(30)`.  Its frozen commands are

```bash
python 04-computation/lrc14_clock_handoff_d_pullback_independent_referee.py
python -O 04-computation/lrc14_clock_handoff_d_pullback_independent_referee.py
```

and both byte-match
`05-knowledge/results/lrc14_clock_handoff_d_pullback_independent_referee.out`.
The referee imports the independently hashed physical two-edge engine named
in the frontmatter; all three LF hashes are recorded there.

QED.
