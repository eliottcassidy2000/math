# LRC14 Moon-Core Proof Skeleton

## What Changed

The current family classifier had five live strict-counterexample families.
One of them should now be treated as closed:

```text
L1 apex-multiple residual:
  S = R union 14Q, |Q|>=7, |R|<=6
  closed by THM-571 gamma descent.
```

The old raw 14-shift sieve failed when two half-step residuals covered both
parities.  THM-571 fixes exactly that by descending to the 7-phase: two
half-step residuals mean at least nine speeds are multiples of 7; scale by 7,
use the below-frontier LRC input, then pigeonhole seven lifts.

So a strict counterexample must now be in the **Moon core**:

```text
qdiv>14
top-balanced
|S cap 14Z| <= 6
zero strict-Haar open front
not wide-gK8 discharged
not unit-petal/GW-strip discharged
not K33 state-lift discharged
not fixed-margin packet reducible
```

## Conditional Proof

HYP-2964 records the conditional theorem:

```text
LRC14 follows from four remaining lemmas:

A. global gK8 concentration extremality for wide packets;
B. zero-open K33 packet -> H=7 TournamentStateLift;
C. fixed-margin packet exhaustiveness for the bounded Moon core;
D. covering boundary-moment positivity outside K33/state-lift.
```

The proof skeleton is:

```text
strict counterexample
  -> qdiv>14                         (THM-523)
  -> top-balanced                     (HYP-2906)
  -> |S cap 14Z|<=6                   (THM-571)
  -> zero regular-open Haar mass      (otherwise open witness)
  -> not wide                         (gK8 concentration)
  -> not K33                          (state lift)
  -> fixed-margin Moon-core packet
  -> audited safe family or moment+
  -> contradiction.
```

## New Object

I introduced a proof-rank:

```text
rho_moon(P) =
  (wide_depth,
   max(0, |S cap 14Z|-6),
   zero_open_flag,
   k33_unlifted_flag,
   source_unknown_flag,
   fixed_margin_component_size)
```

Known reductions kill coordinates of this rank.  The desired proof is that no
packet has minimal positive `rho_moon`.

## Tournament Analysis

Vertices are proof gates, not runners:

```text
Q q-witness
P one-large peel
A apex-majority gamma descent
H Haar-open migration
W wide gK8 concentration
K K33 state lift
F fixed-margin packet theorem
Y boundary-moment positivity
U unnamed source-kernel falsifier
R raw residue/tournament shadow
```

Tie path:

```text
Q > P > A > H > W > K > F > Y > U > R.
```

This is a proof-order tournament: it preserves strict-counterexample status and
source ownership, while raw residue/tournament shadows sit last.

## Handoff

The proof target is no longer broad LRC14.  It is:

```text
prove the Moon core has no fixed-margin source-spectrum component
except those already routed to open fronts, gK8 positivity, or K33 state lift.
```

That is the place to spend the next serious proof effort.
