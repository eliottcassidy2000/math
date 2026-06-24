---
id: HYP-2961
title: LRC14 counterexample family and sporadic classifier
status: CLASSIFICATION TARGET / proof grammar for all possible strict counterexamples, not a proof
source: codex-2026-06-24-S153
related:
  - HYP-2960
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2947
  - HYP-2929
  - HYP-2909
  - HYP-2908
  - HYP-2906
  - HYP-2905
  - HYP-2904
  - HYP-2903
  - HYP-2896
  - HYP-2895
  - HYP-2893
  - THM-523
  - THM-560
  - THM-563
  - THM-564
  - THM-566
  - THM-568
  - THM-572
  - OPEN-Q-108
---

# HYP-2961: LRC14 Counterexample Family And Sporadic Classifier

This is a rigorous classification target for **possible strict
counterexamples** to LRC14.  It is not a proof that the buckets are empty.

The point is to stop saying "a counterexample" as if it were one amorphous
object.  After the repo's q-witness, source-spectrum, packet-migration,
skeleton-gate, wide/decorrelation, and state-lift work, any strict
counterexample has to enter one of a small number of named families or a finite
sporadic bucket.

## Normal Form

An LRC14 row is a primitive 13-element set

```text
S = {s_1 < ... < s_13} subset Z_{>0},  gcd(S)=1.
```

Its lonely value is

```text
M(S) = sup_t min_{s in S} ||s*t||.
```

A **strict counterexample** means

```text
M(S) < 1/14.
```

Define the divisor gate

```text
qdiv(S) = min { d >= 1 : no element of S is divisible by d }.
```

THM-523's q-witness gives:

```text
M(S) >= 1/qdiv(S).
```

Therefore every strict counterexample satisfies:

```text
qdiv(S) > 14.
```

This is the first hard distinction.  `qdiv=14` rows can be tight at the
threshold, but they are not strict counterexamples.  AP and Goddyn-Wong are
equality atoms, not bad rows.

## Equality Atoms Versus Counterexample Candidates

The classifier keeps the tight/equality locus visible because every failed
counterexample route has tried to impersonate it.

Known equality atoms:

```text
E0  AP = {1,2,...,13}
E1  GW = {1,2,...,11,13,24}
```

They have:

```text
qdiv=14,
strict-open safe set empty,
closed denominator-14 boundary support,
AP/GW six-pair owner skeleton.
```

HYP-2960 and HYP-2955 support the local theorem target that the `qdiv=14`
zero-open boundary atom is exactly `E0/E1` after the skeleton gate.  But even
without proving global equality classification, strict counterexamples have
already left this branch, because `qdiv=14` gives `M>=1/14`.

## Classification Theorem Target

The desired theorem is:

> Every primitive strict LRC14 counterexample belongs to exactly one live
> residual class below.  All non-live classes are discharged by an existing
> witness, open-front, induction, packet-migration, wide/decorrelation, or
> state-lift theorem.

The classes are ordered by the first gate that does not discharge the row.

## Discharged Non-Counterexample Classes

### D0. Direct q-witness rows

```text
qdiv(S) <= 14.
```

These have `M(S) >= 1/14`.  They are not strict counterexamples.  The subcase
`qdiv=14` contains equality atoms and open-front rows, but no strict bad row.

### D1. Scale-separated one-large rows

Let `m=s_12` and `v=s_13`.  If

```text
v > 13*m,
```

then HYP-2906 peels the large speed using LRC13 and a one-interval comb
argument.  This discharges the one-huge-speed family.

### D2. Positive Haar-open migration rows

Rows with

```text
O(S) = {t : ||s*t|| > 1/14 for all s in S}
```

of positive Haar measure have `M(S)>1/14`.  HYP-2955 finds that AP-neighborhood
mutations overwhelmingly migrate here: one-swap through `add<=420`, two-swap
through `add<=60`, and three-swap through `add<=30` have no non-AP/GW
boundary-only rows.

### D3. Unit-petal and GW-strip rows

Rows carrying a unit C27 petal or a unit-petal/GW strip have positive open
fronts in the tested atlas:

```text
10->20, 13->26, P10+GW, and related unit-visible splices.
```

They are not strict counterexamples unless a future row keeps the same unit
labels while destroying all open fronts.  Such a row would leave D3 and enter
the live sporadic bucket below.

### D4. K33/state-lift labelled rows

Rows carrying the nonunit `D9` / K33 packet, such as

```text
12->36, P10+K33,
```

are positive-open in the local atlas.  Globally, if such a packet ever becomes
zero-open, the intended route is HYP-2908/THM-572: construct the forbidden
`H=7` tournament state lift.  Thus K33-labelled rows are a family, but not an
unnamed counterexample family; they have a named endpoint obligation.

### D5. Wide rows already covered by known concentration machinery

The wide branch has its own family grammar:

```text
binding single-far rows       -> THM-563 / signed single-far control
generalized doublets r=2      -> THM-564 / HYP-2807 / HYP-2808
genuine-wide gK8-positive     -> HYP-2810, gK8/L_y moment image
finite low-window checks      -> exact bounded audits
```

A strict counterexample in the wide branch must evade the current gK8/L_y
positive image and all finite-window checks.  That survivor is not a generic
wide row; it is live family `L2` below.

## Live Families

### L1. Apex-multiple residual family

This is the sharpened THM-568/HYP-2909 residual.  A row in this family has:

```text
qdiv(S) > 14,
S = R union Q,
Q subset 14*Z,
|Q| >= 7,
R is 14-free,
|R| <= 6,
top-balanced after HYP-2906,
prime-covering after THM-523.
```

Reason it is live: denominator-14 apex witnesses are blocked by `Q`, and the
simple `r<=6` comb-teeth positivity has run out.  HYP-2909 identifies the
remaining theorem as:

```text
S = R union Q, Q subset 14*Z, |Q|>=7
  => M(S) > 1/14.
```

This is an infinite family unless the second-moment/resonance defect supplies
an effective bound.

### L2. Genuine-wide zero-moment family

This is the wide/decorrelation survivor:

```text
span-wide normalized packet,
not binding single-far,
not a discharged generalized doublet,
gK8/L_y image not positive enough,
not in a verified finite window.
```

Reason it is live: all evidence points to gK8 closing this family, but the
classification should leave it named until the global `max_E L_yK8 <= 10cap`
or equivalent concentration theorem is promoted to a proof.

### L3. Bounded covering core family

After omit-prime, dilation, remove-large, and scale-cluster reductions fail,
the row is a primitive bounded covering atom:

```text
qdiv>14,
top-balanced,
prime-covering for 2,3,5,7,11,13,
14-covering,
not scale-separated,
not discharged by local Haar-open migration,
not a known unit/K33 packet,
not wide-gK8-positive.
```

This is the finite atom in HYP-2905.  It is a family before an effective
height/spread bound is proved; after such a bound it becomes the sporadic
search atlas.

### L4. K33 zero-open state-lift family

This is the nonunit C27/K33 family after all open fronts have disappeared:

```text
K33/D9 label present,
strict-open Haar mass zero,
state-lift packet not yet constructed,
qdiv>14 or exact-period covering owner survives.
```

Reason it is live: the family is named, but the proof still has to turn every
such zero-open packet into the HYP-2908/THM-572 forbidden state lift.

### L5. Unnamed source-kernel family

This is the emergency bucket created by HYP-2953:

```text
qdiv>14,
zero strict-open witness,
no positive boundary-moment/gK8 image,
no AP/GW boundary ownership,
no C27 unit-petal discharge,
no K33/state-lift label,
pullback-consistent source-spectrum packet.
```

Reason it is live: this is exactly what a genuinely new counterexample would
look like after all known quotients are retained.  It is the bucket the next
gauntlet should try to prove empty.

## Sporadic Buckets

A **sporadic** is not merely a small row.  It is a primitive row whose
parameters are bounded after the live families have been removed.

The sporadic buckets are:

```text
S0  equality sporadics:
    AP and GW; tight but not strict counterexamples.

S1  local positive-front sporadics:
    12->36, 10->20, 13->26, P10+GW, P10+K33, 8->16,
    12->26, and the finite AP-neighborhood leaders.  All are discharged
    by q-witness or positive Haar mass in current audits.

S2  bounded-core unknown sporadics:
    rows in L3 after an effective bounded-core theorem gives a finite atlas.

S3  exact-period shell-collapse sporadics:
    rows where THM-568 gives D=14h but the proof has not yet forced h=1
    or covering strictness.

S4  state-lift sporadics:
    K33-labelled rows whose packet is finite and zero-open but whose H=7
    lift has not yet been built.

S5  unnamed-source sporadics:
    finite rows in L5 after all infinite parameters have been bounded.
```

Current evidence:

```text
S0 = {AP, GW}.
S1 contains only discharged rows in tested banks.
S2/S3/S4/S5 are proof-obligation buckets; no actual strict counterexample is known.
```

## Exhaustive Decision Procedure Target

The classifier can be implemented as a deterministic decision tree:

```text
input primitive 13-row S

if qdiv(S) <= 14:
    return D0, with equality-subclassifier for AP/GW boundary atoms

if top_speed > 13 * second_top_speed:
    return D1

compute strict Haar front and boundary owner skeleton
if strict Haar mass > 0:
    return D2, refined by C27/K33/local labels

compute source-spectrum pullback labels
if unit C27 petal/GW strip:
    return D3
if K33/D9 label:
    return D4 unless zero-open, then L4

compute wide/decorrelation packet
if binding single-far or generalized-doublet or gK8-positive:
    return D5
if wide and gK8 not discharged:
    return L2

if many apex multiples:
    return L1

if bounded/top-balanced/covering:
    return L3 or bounded-core sporadic bucket S2

if all labels are absent but source-spectrum remains pullback-consistent:
    return L5/S5
```

The most important invariant is that every branch states what predicate it
preserves and what it destroys.  A branch is not allowed to classify by raw
residue, raw tournament class, or scalar `M` alone.

## Tournament Analysis

Assumption challenge: the tournament vertices do not need to be runners,
residues, or arcs.  For a counterexample classifier the natural vertices are
candidate families and proof exits.

Candidate vertex sets considered:

```text
runners,
gaps,
denominator clocks,
denominator-14 apex points,
cover arcs,
exact-period unit packets,
C27 transfers,
K33 incidence packets,
wide moment packets,
state-lift obligations,
families/sporadic buckets.
```

Chosen tournament vertices:

```text
Q   qdiv witness gate
P   scale-peel / boundary-state induction
H   Haar-open migration
B   AP/GW boundary skeleton gate
U   unit-petal/GW strip discharge
K   K33/state-lift label
W   wide gK8/decorrelation branch
A   apex-multiple residual L1
C   bounded covering core L3
N   unnamed source-kernel L5
R   raw residue/tournament shadow
```

Pairwise observable:

```text
Obs(X,Y) =
  (which bucket preserves the strict-counterexample predicate after quotient?,
   which bucket has a named proof exit?,
   which bucket leaves fewer infinite parameters?,
   which bucket destroys less source-spectrum ownership data?)
```

Switch/gauge:

Orient `X -> Y` when `X` is an earlier forced gate in every strict
counterexample decision tree, or when `X` preserves the source/kernel data
needed to route `Y`.  Ties follow:

```text
Q > P > H > B > U > K > W > A > C > N > R.
```

Fingerprint of this proof-order tournament:

```text
score histogram: 10,9,8,7,6,5,4,3,2,1,0
directed 3-cycles: 0
SCCs: eleven singleton SCCs
Hamiltonian paths: unique under the stated tie gauge
```

What it preserves:

```text
strict-counterexample status,
source-spectrum ownership,
open-front discharge,
state-lift address,
and whether an infinite family parameter remains.
```

What it destroys:

```text
exact row identity,
fine wall-crossing order,
and low-level runner labels after they are converted into packet owners.
```

## What Would Count As A Breakthrough

Any of the following would materially advance the LRC14 proof.

1.  Prove L1 empty by a second-moment theorem for `R union 14Q` with
    `|Q|>=7`.
2.  Prove L2 empty by a global gK8/L_y concentration theorem.
3.  Prove L3 finite and empty by a bounded covering core theorem, or reduce it
    to `{AP,GW}` equality plus positive slack.
4.  Prove L4 empty by constructing the HYP-2908/THM-572 state lift for every
    zero-open K33 packet.
5.  Extend HYP-2955 so the classifier emits a source-kernel label for every
    tested row, then verify that L5/S5 remains empty.

The classification target is deliberately falsifiable.  A real strict
counterexample must land in L1, L2, L3, L4, or L5.  If it does not, the
decision tree has dropped a predicate and must be repaired.
