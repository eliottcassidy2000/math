---
id: THM-2684
title: "Three-tooth rail envelope, diagonal arrival law, and full dilation nilpotence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The complete
  positive THM-2584 route-two rail bank has the exact three-tooth support envelope
  B_0=[0,1/28), B_6=[13/28,15/28), B_12=[27/28,1), with no other arrival
  digit present.  Under D(x)={13x}, each tooth meets only its own image tooth,
  so every varying-arrival physical handoff is empty.  The three raw diagonal
  return cylinders are positive, but on each of them the consecutive shallow
  clocks c_7(Dx),c_7(D^2x) agree.  Hence a legal two-event rail-clock product
  exists but every clock-legal three-event rail product on the full inherited
  bank is empty.  This closes the endpoint-carrier escape left by THM-2682,
  not other parent carriers, handoffs, boundary semantics, edge grammars, or
  LRC(14).
source: root-2026-07-28-three-tooth-rail-envelope
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
  - THM-2682-central-arrival-clock-trap-and-three-event-dilation-nilpotence
related:
  - THM-2586-arrival-conditioned-ancestry-sheet-and-later-root-realization
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2670-sharp-graph-clock-incidence-atlas-and-physical-gluing-boundary
script: 04-computation/lrc14_alternate_arrival_physical_rail_handoff.py
output: 05-knowledge/results/lrc14_alternate_arrival_physical_rail_handoff.out
script_sha256: ac4d7f0a9ff67205505306c19c9792d4af89c1a797339f323ed282376cc3b39d
output_sha256: de5e5dd90c9bf59c901aeed289bc55e2a3c5eb2daa288230f960057bf15eb67f
hash_basis: LF-normalized bytes
---

# THM-2684 -- the inherited rails form a three-tooth diagonal dynamical system

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2682 closes three-event dilation chronology on the central arrival-six
carrier.  Its general clock-return calculation identifies endpoint arrivals
`0,12` as the only phase-evasion labels that are also present in the inherited
rail tensor.  That left a precise question: do those two endpoints furnish a
physical bypass, or are they artifacts of projecting rail and phase support
separately?

They are projection artifacts.  Once all positive rails are viewed as one
geometric object, their support is a three-tooth set with a diagonal Markov
law.  A second dilation contracts each surviving return into a forbidden
clock diagonal.  This gives one mechanism for the central and both endpoint
obstructions.

## 1. The full supported rail bank

Retain every positive route-two rail from THM-2584 as the weighted half-open
set

```text
R_(m,s,b,t)(x),                                          (1)
```

where `m` is the arrival digit, `s in {1,...,12}` is the source, `b` is the
owner clock, and `t` is the deep label.  Exact reconstruction gives

```text
(m,t)=(0,0):       81 positive rails,
(m,t)=(6,0),(6,12):81+81 positive rails,
(m,t)=(12,12):     81 positive rails,                    (2)
```

and no positive rail for any other `(m,t)`.  Thus there are `324` rails in
the complete inherited bank.

Forget the labels only after taking the union of their actual positive
weighted supports.  The result is unexpectedly simple:

```text
B_0  = union_(s,b,t) supp R_(0,s,b,t)  =[0,1/28),
B_6  = union_(s,b,t) supp R_(6,s,b,t)  =[13/28,15/28),
B_12 = union_(s,b,t) supp R_(12,s,b,t) =[27/28,1).       (3)
```

These are exact equalities of half-open sets, not bounding intervals and not
measure identities.  Every point of each right-hand side belongs to at least
one positive weighted rail piece; nothing outside the three displayed teeth
belongs to the bank.

## 2. Dilation has identity arrival adjacency

Put `D(x)={13x}`.  None of the three intervals in (3) wraps internally under
this map, and direct endpoint arithmetic gives

```text
D(B_0)  =[0,13/28),
D(B_6)  =[1/28,27/28),
D(B_12) =[15/28,1).                                     (4)
```

Comparing (3) and (4), positive-length intersection occurs exactly on the
diagonal:

```text
length(B_j intersect D(B_i))>0  iff i=j,
                             i,j in {0,6,12}.             (5)
```

The closures touch at the four displayed boundary values, but the half-open
sets themselves are disjoint; another endpoint convention could add only
null points.  Consequently every one of the six varying-arrival
arrows in THM-2682's projected phase atlas disappears at the actual rail
fibre-product layer, before clocks, delayed words, or units are multiplied.

The full cellwise calculation refines (5).  For arrivals `0,6,12`, the exact
numbers of positive raw rail-key pairs are respectively

```text
3,744, 6,696, 3,744,                                    (6)
```

while all six off-diagonal arrival pairs have zero.  After cutting by the
current shallow clock and requiring both stored clock edges to be nonconstant,
the corresponding legal two-event object counts are

```text
2,376, 4,224, 2,376.                                    (7)
```

They cover `132,144,132` source pairs respectively; the two endpoint banks
are identified by reflection.  In particular, the two-event relation is
genuinely nonempty.

## 3. Every third event lies on a forbidden clock diagonal

Because (5) applies to each adjacent pair, a three-event chain can only have
one of the arrival words

```text
000, 666, 12 12 12.                                     (8)
```

The raw three-rail unions are not empty.  Distributing intersection over the
finite rail unions in (3) gives exactly

```text
G_0  =B_0  intersect D^(-1)B_0  intersect D^(-2)B_0
     =[0,1/4732),
G_6  =B_6  intersect D^(-1)B_6  intersect D^(-2)B_6
     =[2365/4732,2367/4732),
G_12 =B_12 intersect D^(-1)B_12 intersect D^(-2)B_12
     =[4731/4732,1),                                    (9)
```

where `4732=28*13^2`.  Thus a proof by raw triple emptiness would be false.

Let

```text
c_7(u)=floor(7{u}+1/2) mod 7.                            (10)
```

On `G_0`, both `Dx` and `D^2x` lie in the clock-zero neighbourhood of zero;
on `G_12`, both lie in the clock-zero neighbourhood of one.  On `G_6`, the
central affine return fixes `1/2` and preserves its two sides, exactly as in
THM-2682.  Hence uniformly on all three positive raw cylinders,

```text
c_7(Dx)=c_7(D^2x).                                      (11)
```

Dilation covariance identifies these two quantities with

```text
shallow(E_0)=c_7(Dx),
owner(E_0)=shallow(E_1)=c_7(D^2x).                       (12)
```

But the stored incidence edge of `E_0` is required to be nonconstant.
Equations (11)--(12) force every raw triple in (9) to be clock-illegal, so no
raw triple belongs to the legal event product.  Combining this with (5) proves

```text
some clock-legal two-rail product is positive;
every clock-legal three-rail product is empty.            (13)
```

Thus support-product nilpotence length three holds at the rail-clock layer.
THM-2680 separately supplies a positive full two-event atom, while the empty
rail-clock gate here forces every downstream full three-event atom product to
be empty.  Every longer product contains a forbidden three-event prefix.

## 4. The abstract three-tooth mechanism

The endpoint arithmetic is not peculiar to `13/7`.  For any `q>=2`, set

```text
p=2q-1,                 delta=1/(4q),
B_-=[0,delta),
B_0=[1/2-delta,1/2+delta),
B_+=[1-delta,1).                                        (14)
```

For `D_p(x)={px}` one has

```text
D_p(B_-)=[0,1/2-delta),
D_p(B_0)=[delta,1-delta),
D_p(B_+)=[1/2+delta,1),                                 (15)
```

so the positive adjacency of these three teeth is the identity matrix.  The
three diagonal two-return cylinders are obtained from (14) by replacing
`delta` with `delta/p^2`.  On each cylinder the consecutive nearest-`q`
clocks of `D_px` and `D_p^2x` agree: both endpoint cylinders lie inside clock
zero, while the central cylinder either stays in one cell (`q` even) or
preserves the side of the boundary at `1/2` (`q` odd).

This abstract lemma explains the numbers in (3)--(11).  What is special to
the LRC calculation is the exact fact that its full positive rail envelope
is precisely the `q=7` instance of (14).

## 5. Scope and reproduction

The conclusion occurs at the nonnegative route-two rail layer after imposing
the inherited nonconstant clock-edge grammar.  Any later restriction by
present packets, sharp states, tooth/carry labels,
delayed words, Bockstein units, endpoint sections, or scalar charts is a
subset and cannot restore an empty three-event support.  Conversely, this
theorem does not rule out a different parent carrier, a handoff other than
`D`, a boundary correspondence, an edge grammar allowing repetition, or an
independent AP/Euler/rank-box proof.  No scalar row is excluded and no
`LRC(14)` conclusion follows.

Run

```bash
python3 04-computation/lrc14_alternate_arrival_physical_rail_handoff.py
python3 -O 04-computation/lrc14_alternate_arrival_physical_rail_handoff.py
```

Both modes must byte-match

```text
05-knowledge/results/lrc14_alternate_arrival_physical_rail_handoff.out.
```

The exact companion reconstructs all `324` weighted rail profiles from the THM-2584
tensor, proves the three envelope identities (3), checks the image/intersection
law (4)--(5), exhausts every labelled two-event object and all `27` arrival
words of length three, and retains positive raw-triple controls so that a
clock-diagonal zero cannot be mistaken for geometric emptiness.  Normal and
optimized executions byte-match the frozen transcript and hashes.  The
independent hostile audit rederived (3)--(11), checked the two grid
refinements and covariance indexing, and verified the `7` clock-covariance,
`324` direct/sequential `D^2`-pullback, `1,200` sparse/ordinary-intersection,
and `8,976` reflected-object controls.

QED.
