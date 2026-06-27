---
id: HYP-2964
title: LRC14 moon-core proof skeleton after apex-majority elimination
status: PROOF-SKELETON / conditional closure route, not a proof
source: codex-2026-06-24-S154
related:
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2960
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2929
  - HYP-2911
  - HYP-2910
  - HYP-2909
  - HYP-2908
  - HYP-2906
  - HYP-2905
  - HYP-2812
  - HYP-2809
  - THM-523
  - THM-568
  - THM-570
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2964: LRC14 Moon-Core Proof Skeleton After Apex-Majority Elimination

This pass tries to prove LRC14 by compressing the current proof stack, not by
adding another label.  The main concrete gain is that HYP-2961's live family
`L1` should be reclassified as discharged once THM-571 is accepted:

```text
L1 apex-multiple residual:
  S = R union 14Q, |Q|>=7, |R|<=6
  -> closed by THM-571 gamma descent.
```

Therefore a strict counterexample has to live in a smaller **Moon core**:

```text
qdiv(S)>14,
|S cap 14Z| <= 6,
top-balanced after HYP-2906,
not discharged by Haar-open migration,
not discharged by unit-petal/GW strips,
not discharged by wide gK8 concentration,
not discharged by K33/state lift,
and not reducible by the fixed-margin labelled-packet theorem.
```

That is the new proof target.

## Theorem Stack

Assume `S` is a primitive 13-speed strict counterexample:

```text
M(S) < 1/14.
```

The existing stack forces the following reductions.

### Gate 1: q-witness

By THM-523:

```text
qdiv(S) > 14.
```

All `qdiv<=14` rows are non-strict.  AP/Goddyn-Wong remain equality atoms, not
strict counterexamples.

### Gate 2: one-large peel

By HYP-2906:

```text
s_13 <= 13*s_12.
```

Otherwise the largest speed peels using the LRC13 margin on the seed.

### Gate 3: apex-majority elimination

Let

```text
M14(S) = S cap 14Z.
```

If `|M14(S)|>=7`, then write

```text
S = R union 14Q,     |Q|>=7, |R|<=6.
```

THM-570 handles the fourteen-shift sieve unless at least two residual speeds
are half-step speeds (`gcd(r,14)=7`).  HYP-2911 records the exact guardrail:
two half-step speeds can cover both parities in a raw `14`-lift sieve.  THM-571
then closes the actual apex-majority branch by descending to the `7`-phase:
with at least two half-step residuals, at least nine speeds are multiples of
`7`; after scaling by `7`, LRC<=13 supplies a strict `7`-phase, and at most
four remaining non-`7` speeds each forbid at most one of seven lifts.

Thus every strict counterexample satisfies:

```text
|M14(S)| <= 6.
```

This removes HYP-2961's `L1` live family from the strict-counterexample list.

### Gate 4: open-front migration

If the regular-open safe set

```text
O(S) = {t : ||s*t|| > 1/14 for all s in S}
```

has positive Haar measure, then `S` is not a counterexample.  HYP-2955,
HYP-2960, HYP-2962, and HYP-2963 all support the same packet-migration
principle: outside AP/GW boundary atoms, local packets seen so far migrate to
positive Haar fronts or to named C27/K33 routes.

So the Moon core also satisfies:

```text
mu(O(S)) = 0.
```

### Gate 5: wide concentration

If `S` is in the wide/decorrelated branch and HYP-2812's gK8 concentration
extremality is promoted to a theorem, then:

```text
10*p0(S) <= L_yK8(S) <= 10*cap,
```

and the wide branch is safe.  The proof route is that gK8 charges only the
extreme miss counts `q0,q6`, while decorrelation moves mass toward the middle.

Therefore the remaining Moon core is bounded/comparable or is an exception to
gK8 concentration.  The latter is now the named `L2` theorem target, not a
generic counterexample.

### Gate 6: K33 state lift

If a zero-open packet carries the nonunit `D9`/K33 label, then the intended
endpoint is HYP-2908/THM-572:

```text
zero-open K33 packet -> TournamentStateLift with packet value H=7
                     -> contradiction by THM-572.
```

Thus a Moon-core counterexample must either be non-K33 or expose the missing
state-lift construction.

### Gate 7: fixed-margin labelled-packet theorem

HYP-2962/HYP-2963 propose the final categorical step.  Map the residual to a
fixed-margin labelled packet:

```text
P(S) =
  q-cover counts,
  exact M/Farey branch,
  Haar open-vs-boundary carrier,
  source-spectrum pullback,
  C27/unital owner labels,
  K33/state-lift labels,
  gK8/L_y moment image.
```

If every fixed-margin packet component is connected by packet-preserving swaps
to one of the audited safe families, then no Moon-core counterexample remains.
The only possible failure is a missing non-scalar sector, HYP-2956's `F7` or
HYP-2961's unnamed source-kernel `L5`.

## Conditional LRC14 Theorem

The current stack proves the following conditional theorem.

```text
If the following lemmas hold:

  A. gK8 concentration extremality:
     all wide packets satisfy L_yK8 <= 10*cap;

  B. K33 state-lift extraction:
     every zero-open nonunit K33 packet constructs the THM-572 H=7 lift;

  C. fixed-margin packet exhaustiveness:
     every bounded Moon-core packet component connects to q-witness,
     AP/GW equality, unit-petal, K33, or covering-moment positive families;

  D. covering boundary-moment positivity:
     every remaining qdiv>14 zero-open packet has positive boundary-moment
     / gK8 / L_y image unless it is a K33 state-lift packet;

then LRC14 holds.
```

Proof sketch:

```text
strict counterexample
  -> qdiv>14                         by THM-523
  -> top-balanced                     by HYP-2906
  -> |S cap 14Z|<=6                   by THM-571
  -> zero regular-open Haar mass      otherwise open witness
  -> not wide                         by Lemma A
  -> not K33                          by Lemma B
  -> fixed-margin Moon-core packet    by source-spectrum packetization
  -> audited safe family or moment+   by Lemmas C,D
  -> contradiction.
```

This is not yet a proof because A-D are still theorem targets.  But the proof
burden is now much smaller and better typed than "prove LRC14."

## New Creative Move: The Moon-Core Rank

Define a proof-rank for a packet:

```text
rho_moon(P) =
  (wide_depth,
   m14_excess = max(0, |S cap 14Z|-6),
   zero_open_flag,
   k33_unlifted_flag,
   source_unknown_flag,
   fixed_margin_component_size)
```

ordered lexicographically after discharging q-witness rows.  Known reductions
strictly lower or zero out one coordinate:

```text
THM-571              kills m14_excess;
Haar-open migration  kills zero_open_flag;
gK8 concentration    kills wide_depth;
K33 state lift       kills k33_unlifted_flag;
fixed-margin swaps   should shrink component_size or expose source_unknown.
```

The proposed proof strategy is to show that no packet can have minimal positive
`rho_moon`.  A minimal positive packet would be:

```text
bounded,
top-balanced,
qdiv>14,
|S cap 14Z|<=6,
zero-open,
non-K33 or unlifted K33,
not wide-positive,
and fixed-margin isolated.
```

HYP-2963's audit finds no such packet in the bounded adversarial bank.

## Tournament Analysis

Assumption challenge: a proof tournament should not use runners as vertices.
Here the vertices are proof gates and live theorem obligations.

Vertices:

```text
Q  q-witness
P  one-large peel / boundary-state induction
A  apex-majority gamma descent
H  Haar-open migration
W  wide gK8 concentration
K  K33 state lift
F  fixed-margin packet theorem
Y  boundary-moment positivity
U  unnamed source-kernel falsifier
R  raw residue/tournament shadow
```

Pairwise observable:

```text
Obs(X,Y) =
  (does X eliminate an infinite parameter?,
   does X preserve source-spectrum ownership?,
   does X have a named proof theorem?,
   does X reduce rho_moon before scalarization?)
```

Switch/gauge:

Orient `X -> Y` when `X` is forced earlier in the strict-counterexample proof
or eliminates an infinite family before `Y`; ties follow:

```text
Q > P > A > H > W > K > F > Y > U > R.
```

Fingerprint:

```text
score histogram: 9,8,7,6,5,4,3,2,1,0
directed 3-cycles: 0
SCCs: ten singleton SCCs
Hamiltonian paths: unique under the tie gauge
```

The tournament preserves strict-counterexample status, source ownership,
family parameter status, and named theorem obligations.  It destroys exact row
identity, which is acceptable only after the packet `P(S)` has been emitted.

## Handoff

The next highest-value proof task is no longer L1.  It is one of:

1. prove gK8 concentration extremality globally;
2. build the K33 zero-open state-lift extraction;
3. prove fixed-margin packet exhaustiveness for the bounded Moon core;
4. add the `Y(S)` moment image to HYP-2963's classifier and verify no
   source-spectrum unknown packets survive.

The most promising creative route is (3): imitate the fixed-margin swap-chain
proof architecture, but keep the Johnson-like non-scalar sectors as
source-spectrum/C27/K33/moment labels rather than scalar counts.
