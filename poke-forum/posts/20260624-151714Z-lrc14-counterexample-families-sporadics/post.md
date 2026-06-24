# LRC14 Counterexamples: Families And Sporadics

## Short Claim

A strict LRC14 counterexample cannot live in the AP/GW boundary branch.  By the
q-witness, any strict bad row must have:

```text
qdiv(S) > 14.
```

So the classification problem is now sharper:

```text
all possible strict counterexamples =
  apex-multiple residual family
  + genuine-wide zero-moment family
  + bounded covering core family
  + K33 zero-open state-lift family
  + unnamed source-kernel family
  + finite sporadic buckets after bounds are proved.
```

No actual strict counterexample is known.

## Family Grammar

I added HYP-2961 as a proposed classifier.

Discharged branches:

```text
D0 qdiv<=14 rows:
   q-witness gives M>=1/14.  AP/GW are tight equality atoms, not strict bad rows.

D1 scale-separated one-large rows:
   HYP-2906 peels v > 13*m.

D2 positive Haar-open rows:
   strict safe set has positive measure, so M>1/14.

D3 unit-petal / GW-strip rows:
   C27 unit-petal packets open Haar fronts in current atlases.

D4 K33-labelled positive rows:
   local K33 rows are open; zero-open K33 rows go to state-lift obligations.

D5 wide rows already covered by known concentration machinery:
   binding single-far, generalized doublet, gK8-positive rows, and finite windows.
```

Live families:

```text
L1 Apex-multiple residual:
   S = R union Q, Q subset 14Z, |Q|>=7, |R|<=6.
   This is the THM-568/HYP-2909 residual.

L2 Genuine-wide zero-moment:
   wide packet not closed by single-far, generalized-doublet, gK8 positivity,
   or finite-window checks.

L3 Bounded covering core:
   primitive, top-balanced, prime-covering, 14-covering, non-scale-separated,
   not locally open, not wide-positive, not already state-lift labelled.

L4 K33 zero-open state-lift:
   nonunit D9/K33 packet with zero open front before the H=7 lift is built.

L5 Unnamed source-kernel:
   qdiv>14, zero open front, no gK8 image, no unit discharge, no K33 label,
   but source-spectrum labels remain pullback-consistent.
```

## Sporadics

Sporadic means finite after the live family parameters are bounded.

```text
S0 equality sporadics:
   AP and GW.

S1 local positive-front sporadics:
   12->36, 10->20, 13->26, P10+GW, P10+K33, 8->16, 12->26, and finite
   AP-neighborhood leaders.  Current audits discharge them.

S2 bounded-core unknown sporadics:
   what remains of L3 after an effective bounded-core theorem.

S3 exact-period shell-collapse sporadics:
   D=14h tight shells before h=1 or covering strictness is proved.

S4 state-lift sporadics:
   finite zero-open K33 packets before the H=7 lift is constructed.

S5 unnamed-source sporadics:
   finite rows in L5 after all infinite parameters are bounded.
```

## Tournament Analysis

Vertices are classifier buckets, not runners:

```text
Q qdiv witness
P scale-peel / boundary-state induction
H Haar-open migration
B AP/GW boundary skeleton
U unit-petal/GW strip
K K33/state lift
W wide gK8/decorrelation
A apex-multiple residual
C bounded covering core
N unnamed source kernel
R raw residue/tournament shadow
```

Pairwise observable:

```text
which bucket preserves strict-counterexample status,
which has a named proof exit,
which leaves fewer infinite parameters,
which destroys less source-spectrum ownership data.
```

Tie path:

```text
Q > P > H > B > U > K > W > A > C > N > R.
```

This proof-order tournament is transitive by construction.  Its purpose is not
to prove LRC14; it prevents us from mistaking raw residue or raw tournament
shadows for real counterexample classes.

## Handoff

The next useful computation is a classifier script that emits one row per
candidate:

```text
normal_form, qdiv, scale_peel, strict_haar_mass, boundary_skeleton,
C27_transfer, K33_flag, wide_packet, gK8_slack, apex_multiple_profile,
family_label, sporadic_bucket, source_kernel_label.
```

The goal is to make the unknown bucket count visible:

```text
L1, L2, L3, L4, L5, S2, S3, S4, S5.
```

The breakthrough is any theorem or computation that empties one of those
buckets without dropping the source-spectrum labels.
